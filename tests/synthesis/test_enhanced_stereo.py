"""Synthetic stereo-state boundaries; no model, embedding, alignment or yield claim."""

import json
from types import SimpleNamespace
from pathlib import Path

import pytest
from rdkit import Chem

from chembounce_synthesis.assets import convert_catalog
from chembounce_synthesis.catalog import Catalog
from chembounce_synthesis.common import absolute_smiles, canonical, molecule, require_absolute_stereo
from chembounce_synthesis.conformers import ConformerSettings, _check_stereo, reference_from_sdf
from chembounce_synthesis.evidence import _component_smiles
from chembounce_synthesis.forward_support import strict_smiles
from chembounce_synthesis.identity import chemical_identity
from chembounce_synthesis.pocket import BoundPoseConstraints
from chembounce_synthesis.reactions import Rules
from chembounce_synthesis.route_audit import audit_route, normalize_product, stereo_state
from chembounce_synthesis.routes import SearchLimits, _proposal_graph, verify_route
from chembounce_synthesis.similarity_constraints import BoundSimilarity, pose_binding

PLAIN = "C[C@H](F)[C@@H](O)Cl"
ABS = PLAIN + " |a:1,3|"
UNCERTAIN = [PLAIN + " |o1:1,3|", PLAIN + " |&1:1,3|", PLAIN + " |a:1,&1:3|", PLAIN + " |o1:1,&1:3|"]
AMIDE = "[C:1](=[O:2])[O;H1:3].[N;H1,H2:4]>>[C:1](=[O:2])[N:4]"


@pytest.fixture(autouse=True)
def no_scientific_calls(monkeypatch):
    from rdkit.Chem import AllChem, rdMolAlign
    from chembounce_synthesis import conformers

    def denied(*args, **kwargs):
        raise AssertionError("Fixture must not embed or align")

    monkeypatch.setattr(AllChem, "EmbedMultipleConfs", denied)
    monkeypatch.setattr(AllChem, "EmbedMolecule", denied)
    monkeypatch.setattr(rdMolAlign, "AlignMol", denied)
    if conformers.rdShapeAlign is not None:
        monkeypatch.setattr(conformers.rdShapeAlign, "AlignMol", denied)


@pytest.mark.parametrize("smiles", UNCERTAIN)
@pytest.mark.parametrize("boundary", [molecule, canonical, chemical_identity, _component_smiles])
def test_text_ingress_rejects_before_stereo_information_loss(smiles, boundary):
    with pytest.raises(ValueError, match="nonabsolute enhanced stereochemistry"):
        boundary(smiles)


@pytest.mark.parametrize("smiles", UNCERTAIN)
def test_raw_molecule_boundaries_and_diagnostics(smiles):
    mol = Chem.MolFromSmiles(smiles)
    before = Chem.MolToCXSmiles(mol)
    for operation in (
        lambda: require_absolute_stereo(mol),
        lambda: absolute_smiles(mol),
        lambda: normalize_product(mol),
        lambda: _check_stereo(mol, ConformerSettings(reject_unspecified_stereo=False)),
        lambda: BoundSimilarity._check_state(mol),
        lambda: pose_binding(mol, 0),
    ):
        with pytest.raises(ValueError, match="nonabsolute enhanced stereochemistry"):
            operation()
    state = stereo_state(mol)
    assert not state["fully_specified"] and state["unsupported_enhanced_stereo"]
    assert all(group["atom_indices"] for group in state["unsupported_enhanced_stereo"])
    assert Chem.MolToCXSmiles(mol) == before


@pytest.mark.parametrize("smiles", [PLAIN, ABS, "N[C@@H]([2H])C(=O)O", "C/C=C/C", "[13CH3][C@H](F)Cl"])
def test_supported_absolute_stereo_has_unchanged_serialization(smiles):
    mol = Chem.MolFromSmiles(smiles)
    expected = Chem.MolToSmiles(mol, isomericSmiles=True)
    assert canonical(smiles) == absolute_smiles(mol, isomericSmiles=True) == expected
    assert _component_smiles(smiles) == expected
    assert stereo_state(mol)["fully_specified"]
    assert _check_stereo(mol, ConformerSettings()) == 0


@pytest.mark.parametrize("smiles", UNCERTAIN)
def test_core_input_parser_and_raw_molecule_library_reject(smiles):
    import utils
    import chemistry
    from chembounce import search_similar_scaffolds

    with pytest.raises(utils.SmilesInputError, match="nonabsolute enhanced"):
        utils.init_mol(smiles)
    with pytest.raises(ValueError, match="nonabsolute enhanced"):
        chemistry.prepare_replacement(Chem.MolFromSmiles(smiles), Chem.MolFromSmiles("C"))
    with pytest.raises(ValueError, match="nonabsolute enhanced"):
        search_similar_scaffolds(Chem.MolFromSmiles("CO"), [Chem.MolFromSmiles(smiles)])


def catalog_text(path, smiles):
    path.write_text(f"id\tsmiles\tsource\tsnapshot\ninput\t{smiles}\tsynthetic test\t2026-10-17\n")
    return path


@pytest.mark.parametrize("smiles", UNCERTAIN)
def test_catalog_tsv_cannot_promote_relative_or_mixture_state(tmp_path, smiles):
    with pytest.raises(ValueError, match="nonabsolute enhanced"):
        Catalog(catalog_text(tmp_path / "catalog.tsv", smiles))


@pytest.mark.parametrize("suffix", [".smi", ".smiles", ".txt", ".csv", ".tsv", ".sdf"])
@pytest.mark.parametrize("smiles", UNCERTAIN[:2])
def test_converter_retains_rejection_and_does_not_treat_cx_as_identifier(tmp_path, suffix, smiles):
    source = tmp_path / ("input" + suffix)
    if suffix == ".sdf":
        with Chem.SDWriter(str(source)) as writer:
            writer.SetForceV3000(True)
            for text, ident in ((smiles, "uncertain"), ("CO", "good"), (ABS, "absolute")):
                mol = Chem.MolFromSmiles(text)
                mol.SetProp("id", ident)
                writer.write(mol)
    elif suffix in (".csv", ".tsv"):
        import csv

        with source.open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t" if suffix == ".tsv" else ",")
            writer.writerows([["smiles", "id"], [smiles, "uncertain"], ["CO", "good"], [ABS, "absolute"]])
    else:
        source.write_text(f"{smiles} uncertain\nCO good\n{ABS} absolute\n")
    original = source.read_bytes()
    output = tmp_path / "converted.tsv"
    result = convert_catalog([source], output, "synthetic only", "2026-10-17")
    assert result["records"] == 2 and len(result["rejected_records"]) == 1
    assert "nonabsolute enhanced" in result["rejected_records"][0]["reason"]
    assert result["rejected_records"][0]["input_index"] == 0
    assert result["rejected_records"][0]["record"] == 0
    assert set(Catalog(output).blocks) == {"good", "absolute"}
    assert source.read_bytes() == original


@pytest.mark.parametrize("suffix", [".smi", ".sdf"])
def test_converter_preserves_explicit_abs_and_real_identifier(tmp_path, suffix):
    source = tmp_path / ("input" + suffix)
    if suffix == ".smi":
        source.write_text(ABS + " supplied-id\n")
    else:
        with Chem.SDWriter(str(source)) as writer:
            writer.SetForceV3000(True)
            mol = Chem.MolFromSmiles(ABS)
            mol.SetProp("id", "supplied-id")
            writer.write(mol)
    out = tmp_path / "output.tsv"
    result = convert_catalog([source], out, "synthetic only", "2026-10-17")
    assert result["records"] == 1 and not result["rejected_records"]
    assert Catalog(out).get("supplied-id").smiles == canonical(PLAIN)


@pytest.mark.parametrize("smiles", UNCERTAIN[:2])
def test_reference_sdf_and_pocket_reject_before_pose_processing(tmp_path, smiles):
    source = tmp_path / "reference.sdf"
    with Chem.SDWriter(str(source)) as writer:
        writer.SetForceV3000(True)
        writer.write(Chem.MolFromSmiles(smiles))
    with pytest.raises(ValueError, match="nonabsolute enhanced"):
        reference_from_sdf(source, canonical(PLAIN))
    # Only the raw supplier check should execute: no feature/coordinate/receptor work.
    hit_pose = SimpleNamespace(bound_pose=True, mol=Chem.MolFromSmiles(PLAIN), relative_energies={0: 0})
    pocket = {
        "receptor_pdb": "not-opened.pdb",
        "reference_ligand": {"chain": "A", "residue_number": 1, "residue_name": "LIG"},
    }
    with pytest.raises(ValueError, match="nonabsolute enhanced"):
        BoundPoseConstraints(pocket, hit_pose, source)


def reaction_assets(tmp_path, amino):
    path = tmp_path / "catalog.tsv"
    path.write_text(
        "id\tsmiles\tsource\tsnapshot\n"
        + "".join(
            f"{i}\t{s}\tsynthetic only\t2026-10-17\n"
            for i, s in {"a": "CC(=O)O", "g": amino, "n": "CN"}.items()
        )
    )
    rulefile = tmp_path / "rules.json"
    rulefile.write_text(
        json.dumps(
            {"schema_version": 1, "rules": [{"id": "amide", "smarts": AMIDE, "source": "synthetic only"}]}
        )
    )
    return Catalog(path), Rules(rulefile)


@pytest.mark.parametrize("group", ["o1:1", "&1:1"])
def test_previous_two_step_false_positive_is_stopped_at_input(tmp_path, group):
    # Before the fix this input lost its enhanced group, then both amide steps
    # and strict cached forward validation passed with a single absolute state.
    with pytest.raises(ValueError, match="nonabsolute enhanced"):
        reaction_assets(tmp_path, "N[C@@H](C)C(=O)O |" + group + "|")


def test_two_step_absolute_control_and_route_claim_ingress(tmp_path):
    catalog, rules = reaction_assets(tmp_path, "N[C@@H](C)C(=O)O |a:1|")
    first = rules.apply("amide", (catalog.get("a").smiles, catalog.get("g").smiles))[0][0]
    final = rules.apply("amide", (first, catalog.get("n").smiles))[0][0]
    route = {
        "schema_version": 1,
        "root": "r1",
        "nodes": [{"id": i, "building_block_id": i} for i in ("a", "g", "n")]
        + [
            {"id": "r0", "rule_id": "amide", "reactants": ["a", "g"], "product": first},
            {"id": "r1", "rule_id": "amide", "reactants": ["r0", "n"], "product": final},
        ],
    }
    audit = audit_route(route, catalog, rules)
    # Concrete origins are runtime-capability dependent; no false unsupported-runtime success.
    if all(
        a.HasProp("react_idx")
        for a in rules.get("amide")
        .reaction.RunReactants((molecule("CC(=O)O"), molecule("CN")))[0][0]
        .GetAtoms()
    ):
        assert audit["strict_eligible"]
    assert len(audit["steps"]) == 2
    assert (
        audit["steps"][1]["reactant_slots"][0]["parent_binding_sha256"]
        == audit["steps"][0]["step_binding_sha256"]
    )
    route["nodes"][-1]["product"] = UNCERTAIN[0]
    with pytest.raises(ValueError, match="nonabsolute enhanced"):
        verify_route(route, catalog, rules)
    proposal = {
        "reactions": [
            {"reaction_template": "amide", "product": UNCERTAIN[0], "reactants": ["CC(=O)O", "CN"]}
        ],
        "building_blocks": ["CC(=O)O", "CN"],
    }
    with pytest.raises(ValueError, match="nonabsolute enhanced"):
        _proposal_graph(proposal, rules, SearchLimits())


@pytest.mark.parametrize("role", ["hit", "candidate", "design"])
def test_similarity_preflight_rejects_without_embedding(role):
    from chembounce_synthesis.selective_scoring import preflight_design

    if role == "hit":
        operation = lambda: BoundSimilarity(UNCERTAIN[0])
    elif role == "candidate":
        operation = lambda: BoundSimilarity(PLAIN).preflight_candidate(UNCERTAIN[0])
    else:
        operation = lambda: preflight_design(UNCERTAIN[0])
    with pytest.raises(ValueError, match="nonabsolute enhanced"):
        operation()


def test_forward_input_still_disallows_cx_annotations():
    for smiles in UNCERTAIN + [ABS]:
        with pytest.raises(ValueError):
            strict_smiles(smiles)
    assert strict_smiles(PLAIN) == canonical(PLAIN)


ATROP_BASE = "CC1=CC=CC(I)=C1N1C(C)=CC=C1Br"


@pytest.mark.parametrize("suffix", ["wU:7.7", "wD:7.7", "$label$", "(0,0,;1,0,)", "a:1,3,wU:7.7"])
def test_nonabs_cx_is_rejected_before_a_runtime_can_ignore_it(suffix):
    text = ATROP_BASE + " |" + suffix + "|"
    for operation in (canonical, _component_smiles):
        with pytest.raises(ValueError, match="Unsupported CXSMILES"):
            operation(text)


@pytest.mark.parametrize(
    "text",
    [
        PLAIN + " label",
        ABS + " extra",
        ABS + " |a:1|",
        PLAIN + " |a:999|",
        PLAIN + " |a:1,1|",
        PLAIN + " |a:0|",
        PLAIN + " |a:1,3,999|",
        PLAIN + " |a:1| garbage |",
        "F[Pt@SP1](Cl)(Br)I |a:1|",
    ],
)
def test_exact_text_consumes_all_metadata_and_validates_abs_references(text):
    with pytest.raises(ValueError):
        canonical(text)


@pytest.mark.parametrize("version", ["v2000", "v3000"])
def test_atrop_sdf_rejected_even_when_parser_loses_axial_state(tmp_path, version):
    source = Path(__file__).parent / "fixtures" / ("atrop-up-" + version + ".sdf")
    raw = list(Chem.SDMolSupplier(str(source), removeHs=False))[0]
    # Literal synthetic fixture produced by RDKit2025.09.6 from the documented
    # axial example. Old2023 parses it without an atrop BondStereo but keeps
    # unconsumed wedge properties. This test must pass in both runtimes.
    with pytest.raises(ValueError, match="atropisomeric|wedging"):
        absolute_smiles(raw)
    assert not stereo_state(raw)["fully_specified"]
    with pytest.raises(ValueError, match="atropisomeric|wedging"):
        reference_from_sdf(source, ATROP_BASE)
    # A subsequent ordinary record survives, with the original rejected ordinal.
    combined = tmp_path / "mixed.sdf"
    good = tmp_path / "good.sdf"
    with Chem.SDWriter(str(good)) as writer:
        mol = Chem.MolFromSmiles("CO")
        mol.SetProp("id", "good")
        writer.write(mol)
    combined.write_bytes(source.read_bytes() + good.read_bytes())
    output = tmp_path / "catalog.tsv"
    result = convert_catalog([combined], output, "synthetic", "2026-10-17")
    assert result["records"] == 1 and result["rejected_records"][0]["record"] == 0
    assert (
        "atropisomeric" in result["rejected_records"][0]["reason"]
        or "wedging" in result["rejected_records"][0]["reason"]
    )
    assert set(Catalog(output).blocks) == {"good"}


@pytest.mark.parametrize("suffix", ["wU:7.7", "wD:7.7"])
def test_atrop_smi_and_second_annotation_cannot_hide_in_identifier(tmp_path, suffix):
    source = tmp_path / "input.smi"
    source.write_text(ATROP_BASE + " |" + suffix + "| bad\n" + ABS + " |" + suffix + "| hidden\nCO good\n")
    output = tmp_path / "catalog.tsv"
    result = convert_catalog([source], output, "synthetic", "2026-10-17")
    assert result["records"] == 1
    assert [r["record"] for r in result["rejected_records"]] == [0, 1]
    assert all("Unsupported CXSMILES" in r["reason"] for r in result["rejected_records"])


def test_raw_explicit_atrop_bond_guard_without_chemistry_generation():
    # Parsing is not an atrop capability on older RDKit: skip only this explicit
    # raw-enum fixture; the literal SDF and pre-parser string tests run everywhere.
    mol = Chem.MolFromSmiles(ATROP_BASE + " |wU:7.7|")
    if not any(str(b.GetStereo()).startswith("STEREOATROP") for b in mol.GetBonds()):
        pytest.skip("Runtime parser lacks explicit atrop bond states")
    for operation in (absolute_smiles, normalize_product, lambda m: _check_stereo(m, ConformerSettings())):
        with pytest.raises(ValueError, match="atropisomeric"):
            operation(mol)


def test_both_packaged_profiles_are_honestly_source_mismatched():
    from chembounce_synthesis import profiles
    from chembounce_synthesis.workflow import _config

    for name in (profiles.PROFILE, profiles.SELECTIVE_PROFILE):
        config = _config(
            {
                "schema_version": 1,
                "catalog": "unused.tsv",
                "rules": "unused.json",
                "provider": {"kind": "saved", "path": "unused.json"},
                "recommendation_profile": name,
                "synthesis_evidence": {
                    "profile": "computational-route-support-v1",
                    "predictor": {"kind": "reactiont5", "model_path": "unused"},
                },
            }
        )
        provenance = profiles.profile_provenance(config)
        assert provenance["status"] == "custom_unvalidated"
        assert {"chembounce_synthesis/common.py", "chembounce_synthesis/conformers.py"} <= provenance[
            "runtime_source_mismatches"
        ].keys()
        assert not provenance["material_overrides"]


def test_unterminated_smi_annotation_retains_record_rejection(tmp_path):
    source = tmp_path / "input.smi"
    source.write_text(PLAIN + " |a:1\nCO good\n")
    result = convert_catalog([source], tmp_path / "out.tsv", "synthetic", "2026-10-17")
    assert result["records"] == 1
    assert result["rejected_records"][0]["record"] == 0
    assert "Unsupported CXSMILES" in result["rejected_records"][0]["reason"]


@pytest.mark.parametrize("annotation", ["o1:1", "&1:1", "wU:1.0", "a:1"])
def test_reaction_template_cx_not_silently_ignored_by_old_parser(tmp_path, annotation):
    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [{"id": "r", "smarts": AMIDE + " |" + annotation + "|", "source": "synthetic"}],
            }
        )
    )
    with pytest.raises(ValueError, match="CX reaction SMARTS annotations are unsupported"):
        Rules(path)
