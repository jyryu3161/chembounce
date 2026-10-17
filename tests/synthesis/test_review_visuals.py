"""Tiny synthetic display fixtures; no inferred 3D coordinates or scientific claims."""

import copy
import hashlib
import json
import struct
import zlib

import pytest
from rdkit import Chem

from chembounce_synthesis import review_display as d, review_visuals as p, report
from chembounce_synthesis.common import molecule, canonical
from chembounce_synthesis.similarity_constraints import pose_binding
from .test_report import bundle, save


def manual_pose(smiles, cid, offset=0):
    mol = molecule(smiles)
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.SetId(cid)
    conf.Set3D(True)
    for index in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(index, (index * 0.7 + offset, index % 3 * 0.5, index % 2 * 0.4))
    mol.AddConformer(conf, assignId=False)
    return pose_binding(mol, cid)


def real_graph_bundle():
    """Producer-shaped metadata with manually assigned nonphysical test coordinates."""
    docs = bundle()
    trials = docs["refinement.json"]["candidates"][0]["trials"]
    row = docs["results.json"]["results"][0]
    hit = docs["manifest.json"]["hit"]
    for trial in trials:
        seed = trial["seed"]
        candidate = manual_pose(row["smiles"], seed, 0.2)
        reference = manual_pose(hit, seed + 100)
        trial["selective_pose_evidence"] = {
            name: copy.deepcopy(candidate)
            for name in ("hit_candidate", "design_candidate", "design_reference")
        }
        witness = trial["similarity_constraints"]["witness"]
        witness["reference_pose"] = reference
        witness["candidate_pose_sha256"] = candidate["sha256"]
    rep = trials[1]
    row["similarity_witness_sha256"] = d.binding(rep["similarity_constraints"]["witness"])
    docs["scores.json"]["scores"][0].update(copy.deepcopy(rep))
    return docs


def attach(docs):
    hit = docs["manifest.json"]["hit"]
    artifact = p.build_visuals(hit, docs["results.json"]["results"], docs["refinement.json"]["candidates"])
    docs[d.ARTIFACT] = artifact
    raw = json.dumps(artifact).encode()
    docs["manifest.json"]["review_visuals"] = {
        "version": d.VERSION,
        "status": "available",
        "path": d.ARTIFACT,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "publication_completion_within_budget": None,
        "timing_scope": "final_persistence_unverified",
        "prepublication_elapsed_seconds": 0.0,
    }
    return docs


def rebind(docs):
    docs["manifest.json"]["review_visuals"]["sha256"] = hashlib.sha256(
        json.dumps(docs[d.ARTIFACT]).encode()
    ).hexdigest()


@pytest.fixture
def visual_docs():
    return attach(real_graph_bundle())


@pytest.mark.parametrize(
    "smiles",
    [
        "CC(=O)Nc1ccccc1",
        "C[C@H](F)Cl",
        "C[C@H](F)Cl |a:1|",
        "[H][C@](N)(C)C(=O)O |a:1|",
        "[2H][C@](F)(Cl)Br",
        "F/C=C/Cl",
        "C[NH2+]C",
    ],
)
def test_actual_png_identity_and_static_encoding(smiles):
    exact = canonical(smiles)
    result = p.depiction(exact)
    assert result["status"] == "available" and result["smiles"] == exact
    assert result["graph"]["smiles"] == exact
    import base64

    assert d.checked_png(base64.b64decode(result["png_base64"])) == (480, 320)


def test_producer_reader_real_indexed_graph_and_separate_seed_pairs(visual_docs, tmp_path):
    before = copy.deepcopy(visual_docs)
    save(tmp_path, visual_docs)
    view = report.write_report(tmp_path)
    assert view["artifact_consistency"] == "consistent_recorded_bundle", view
    assert not view["visual_issues"]
    visual = view["candidates"][0]["visuals"]
    assert [s["seed"] for s in visual["seeds"]] == [42, 43, 44]
    assert all(s["pair"]["status"] == "available" for s in visual["seeds"])
    assert visual["representative_seed"] == 43
    assert visual_docs == before
    page = (tmp_path / "review.html").read_text()
    assert "data:image/png;base64," in page and "<script>" in page and "Content-Security-Policy" in page
    assert "Reference pose" in page and "new ResizeObserver" in page


@pytest.mark.parametrize(
    "mutation",
    [
        lambda a: a["candidates"].clear(),
        lambda a: a["candidates"][0]["seeds"].pop(),
        lambda a: a["candidates"][0]["seeds"].reverse(),
        lambda a: a["candidates"][0].update(representative_seed=42),
        lambda a: a["candidates"][0]["key"].__setitem__(0, "other-route"),
        lambda a: a["candidates"][0]["seeds"][0]["pair"].update(candidate_pose_sha256="f" * 64),
        lambda a: a["candidates"][0]["depictions"]["original_hit"].update(png_sha256="f" * 64),
        lambda a: a["candidates"][0]["depictions"]["original_hit"].update(width=2),
        lambda a: a["candidates"][0]["seeds"][0]["pair"]["reference_graph"]["atoms"][0].update(
            atomic_number=8
        ),
    ],
)
def test_rebound_corruption_is_visual_diagnostic_not_qualification_change(visual_docs, tmp_path, mutation):
    mutation(visual_docs[d.ARTIFACT])
    rebind(visual_docs)
    save(tmp_path, visual_docs)
    view = report.build_report(tmp_path)
    assert view["visual_issues"]
    assert view["candidates"][0]["visuals"]["status"] == "inconsistent"
    assert view["artifact_consistency"] == "consistent_recorded_bundle"
    assert "<script>" not in report.render_html(view)


def test_graph_rebound_projection_is_checked(visual_docs, tmp_path):
    graph = visual_docs[d.ARTIFACT]["candidates"][0]["seeds"][0]["pair"]["reference_graph"]
    graph["atoms"][0].update(atomic_number=8, element="O")
    graph["sha256"] = d.binding({k: v for k, v in graph.items() if k != "sha256"})
    rebind(visual_docs)
    save(tmp_path, visual_docs)
    assert "topology" in report.build_report(tmp_path)["visual_issues"][0]


def chunk(kind, data):
    return (
        struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)
    )


@pytest.mark.parametrize("kind", [b"acTL", b"fcTL", b"fdAT"])
def test_png_animation_rejected(kind):
    import base64

    raw = base64.b64decode(p.depiction("CCO")["png_base64"])
    with pytest.raises(ValueError, match="Animated"):
        d.checked_png(raw[:33] + chunk(kind, b"\0" * 8) + raw[33:])


@pytest.mark.parametrize(
    "change",
    [lambda b: b + b, lambda b: b + b"x", lambda b: b[:40], lambda b: b[:45] + bytes([b[45] ^ 1]) + b[46:]],
)
def test_png_partial_trailing_and_crc(change):
    import base64

    raw = base64.b64decode(p.depiction("CCO")["png_base64"])
    with pytest.raises(ValueError):
        d.checked_png(change(raw))


def test_zero_expensive_calls_and_inputs_preserved(monkeypatch):
    from chembounce_synthesis import conformers, selective_scoring
    from rdkit.Chem import AllChem, rdMolAlign, rdDistGeom

    docs = real_graph_bundle()
    before = copy.deepcopy(docs)

    def forbidden(*a, **k):
        pytest.fail("Scientific computation invoked by display")

    for module, name in [
        (conformers, "ensemble"),
        (selective_scoring, "score_candidate"),
        (AllChem, "EmbedMultipleConfs"),
        (AllChem, "EmbedMolecule"),
        (rdDistGeom, "EmbedMolecule"),
        (rdDistGeom, "EmbedMultipleConfs"),
        (rdMolAlign, "AlignMol"),
    ]:
        monkeypatch.setattr(module, name, forbidden)
    if getattr(conformers, "rdShapeAlign", None) is not None:
        monkeypatch.setattr(conformers.rdShapeAlign, "AlignMol", forbidden)
    artifact = p.build_visuals(
        docs["manifest.json"]["hit"], docs["results.json"]["results"], docs["refinement.json"]["candidates"]
    )
    assert artifact["candidates"][0]["seeds"][0]["pair"]["status"] == "available"
    assert docs == before


def test_expired_and_overrunning_native_operation(monkeypatch):
    docs = real_graph_bundle()
    args = (
        docs["manifest.json"]["hit"],
        docs["results.json"]["results"],
        docs["refinement.json"]["candidates"],
    )
    monkeypatch.setattr(p, "depiction", lambda _: pytest.fail("expired work started"))
    with pytest.raises(TimeoutError):
        p.build_visuals(*args, deadline_monotonic=1, clock=lambda: 1)
    now = [0.0]

    def slow(_):
        now[0] = 6.0
        return {"status": "available"}

    monkeypatch.setattr(p, "depiction", slow)
    with pytest.raises(TimeoutError):
        p.build_visuals(*args, clock=lambda: now[0])


def test_capability_failure_does_not_drop_seed_coverage(monkeypatch):
    docs = real_graph_bundle()
    monkeypatch.setattr(p, "depiction", lambda _: (_ for _ in ()).throw(RuntimeError("Cairo unavailable")))
    artifact = p.build_visuals(
        docs["manifest.json"]["hit"], docs["results.json"]["results"], docs["refinement.json"]["candidates"]
    )
    assert all(v["status"] == "unavailable" for v in artifact["candidates"][0]["depictions"].values())
    assert len(artifact["candidates"][0]["seeds"]) == 3


@pytest.mark.parametrize("failure", ["build", "artifact", "manifest"])
def test_publication_failure_preserves_every_core_byte(monkeypatch, tmp_path, failure):
    docs = real_graph_bundle()
    docs["manifest.json"]["review_visuals"] = {
        "version": d.VERSION,
        **p.unavailable("optional_display_not_published"),
    }
    save(tmp_path, docs)
    before = {name: (tmp_path / name).read_bytes() for name in docs}
    original = p._atomic
    if failure == "build":
        monkeypatch.setattr(p, "build_visuals", lambda *a, **k: (_ for _ in ()).throw(MemoryError()))
    else:

        def fail(path, *args):
            if path.name == (d.ARTIFACT if failure == "artifact" else "manifest.json"):
                raise OSError("controlled persistence failure")
            return original(path, *args)

        monkeypatch.setattr(p, "_atomic", fail)
    outcome = p.publish_review_visuals(
        tmp_path,
        docs["manifest.json"]["hit"],
        docs["results.json"]["results"],
        docs["refinement.json"]["candidates"],
    )
    assert outcome["status"] == "unavailable"
    assert before == {name: (tmp_path / name).read_bytes() for name in docs}
    assert report.build_report(tmp_path)["candidates"][0]["visuals"]["status"] == "unavailable"


def test_legacy_and_unavailable_representative_are_not_replaced(visual_docs, tmp_path):
    visual_docs[d.ARTIFACT]["candidates"][0]["seeds"][1]["pair"] = p.unavailable(
        "saved geometry not available"
    )
    rebind(visual_docs)
    save(tmp_path, visual_docs)
    view = report.build_report(tmp_path)
    assert not view["visual_issues"]
    page = report.render_html(view)
    assert "selected value='43'" in page
    assert "43 — unavailable" in page
    for name in (d.ARTIFACT,):
        (tmp_path / name).unlink()
    visual_docs.pop(d.ARTIFACT)
    visual_docs["manifest.json"].pop("review_visuals")
    save(tmp_path, visual_docs)
    view = report.build_report(tmp_path)
    assert view["candidates"][0]["visuals"]["status"] == "not_exported"
    assert "<script>" not in report.render_html(view)


def test_workflow_optional_export_uses_actual_selected_all_seed_evidence(
    tmp_path, assets, proposal, strict_predictor, caplog
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    from chembounce_synthesis import recommend

    archive = tmp_path / "proposal.json"
    archive.write_text(
        json.dumps({"schema_version": 1, "proposals": [{"design": "CC(=O)Nc1ccncc1", "response": proposal}]})
    )
    model = tmp_path / "synthetic-model-directory"
    model.mkdir()
    config = {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(archive)},
        "synthesis_evidence": {"predictor": {"kind": "reactiont5", "model_path": str(model)}},
        "search": {"building_block_top_k": 1, "beam_width": 1},
        "max_designs": 1,
        "max_design_attempts": 4,
        "max_scored_products": 1,
        "max_results": 1,
        "refinement": {"top_n": 1, "max_total_conformers": 24},
    }
    caplog.set_level("INFO", logger="chembounce_synthesis.workflow")
    out = tmp_path / "new-run"
    rows = recommend("CC(=O)Nc1ccccc1", ["CC(=O)Nc1ccncc1"], config, out)
    assert len(rows) == 1 and rows[0]["recommendation_status"] == "qualified"
    assert "post-return observation" in caplog.text and "elapsed_seconds" in caplog.text
    view = report.write_report(out)
    assert view["artifact_consistency"] == "consistent_recorded_bundle", view
    assert not view["visual_issues"], view["visual_issues"]
    visual = view["candidates"][0]["visuals"]
    assert visual["status"] == "available", visual
    assert [s["seed"] for s in visual["seeds"]] == [42, 43, 44]
    assert all(s["pair"]["status"] == "available" for s in visual["seeds"])
    saved = json.loads((out / "refinement.json").read_text())["candidates"][0]["trials"]
    for displayed, original in zip(visual["seeds"], saved):
        assert (
            displayed["coordinates"]["reference"]
            == original["similarity_constraints"]["witness"]["reference_pose"]["coordinates"]
        )
        assert (
            displayed["coordinates"]["candidate"]
            == original["selective_pose_evidence"]["hit_candidate"]["coordinates"]
        )
    # This controlled workflow uses genuine producer data; no model was loaded.


@pytest.mark.parametrize(
    "edit",
    [
        lambda g: g["atoms"][0].update(index=True),
        lambda g: g["atoms"][0].update(isotope=-1),
        lambda g: g["bonds"].append(copy.deepcopy(g["bonds"][0])),
        lambda g: g["bonds"][0].update(end=999),
    ],
)
def test_graph_schema_even_after_selfhash_rebound(visual_docs, tmp_path, edit):
    g = visual_docs[d.ARTIFACT]["candidates"][0]["seeds"][0]["pair"]["reference_graph"]
    edit(g)
    g["sha256"] = d.binding({k: v for k, v in g.items() if k != "sha256"})
    rebind(visual_docs)
    save(tmp_path, visual_docs)
    assert report.build_report(tmp_path)["visual_issues"]


def test_topology_excludes_stereo_not_claimed_stereo_authentication():
    mol = molecule("C[C@H](F)Cl")
    graph = p.graph_record(mol, canonical("C[C@H](F)Cl"))
    old = graph["atom_order_sha256"]
    graph["atoms"][1]["cip"] = "S" if graph["atoms"][1]["cip"] == "R" else "R"
    graph["sha256"] = d.binding({k: v for k, v in graph.items() if k != "sha256"})
    assert d.checked_graph(graph)["atom_order_sha256"] == old


def test_saved_pose_and_caller_graph_remain_unchanged():
    from chembounce_synthesis.similarity_constraints import restore_pose

    saved = manual_pose(canonical("[2H][C@](F)(Cl)Br"), 3)
    mol = restore_pose(saved)
    before = mol.ToBinary()
    p.graph_record(mol, saved["smiles"])
    assert mol.ToBinary() == before and pose_binding(mol, 3) == saved


def test_oversized_and_malformed_artifact_has_separate_diagnostics(tmp_path):
    docs = real_graph_bundle()
    docs["manifest.json"]["review_visuals"] = {
        "version": d.VERSION,
        "status": "available",
        "path": d.ARTIFACT,
        "sha256": "f" * 64,
    }
    save(tmp_path, docs)
    (tmp_path / d.ARTIFACT).write_text('{"schema_version":1,"x":NaN}')
    view = report.build_report(tmp_path)
    assert view["visual_issues"] and view["artifact_consistency"] == "consistent_recorded_bundle"
    assert "<script>" not in report.render_html(view)


def test_visual_symlink_and_snapshot_change_preserve_input(tmp_path, visual_docs, monkeypatch):
    save(tmp_path, visual_docs)
    original = (tmp_path / d.ARTIFACT).read_bytes()
    target = tmp_path.parent / (tmp_path.name + "-outside.json")
    target.write_bytes(original)
    (tmp_path / d.ARTIFACT).unlink()
    (tmp_path / d.ARTIFACT).symlink_to(target)
    with pytest.raises(ValueError, match="Input presence changed"):
        report.build_report(tmp_path)
    assert target.read_bytes() == original


def test_png_capability_and_limit_failures_do_not_mutate_input(monkeypatch):
    from rdkit.Chem.Draw import rdMolDraw2D

    docs = real_graph_bundle()
    before = copy.deepcopy(docs)
    monkeypatch.delattr(rdMolDraw2D, "MolDraw2DCairo")
    artifact = p.build_visuals(
        docs["manifest.json"]["hit"], docs["results.json"]["results"], docs["refinement.json"]["candidates"]
    )
    assert all(image["status"] == "unavailable" for image in artifact["candidates"][0]["depictions"].values())
    assert docs == before


def test_atomic_final_replace_overrun_is_honest_content_not_budget_success(tmp_path, monkeypatch):
    docs = real_graph_bundle()
    save(tmp_path, docs)
    now = [0.0]
    replace = p.os.replace

    def delayed(source, destination):
        replace(source, destination)
        if destination.name == "manifest.json":
            now[0] = 6.0

    monkeypatch.setattr(p.os, "replace", delayed)
    result = p.publish_review_visuals(
        tmp_path,
        docs["manifest.json"]["hit"],
        docs["results.json"]["results"],
        docs["refinement.json"]["candidates"],
        clock=lambda: now[0],
    )
    assert result == {"status": "available", "elapsed_seconds": 6.0, "deadline_overrun": True}
    descriptor = json.loads((tmp_path / "manifest.json").read_text())["review_visuals"]
    assert descriptor["publication_completion_within_budget"] is None
    assert descriptor["timing_scope"] == "final_persistence_unverified"
    assert not report.build_report(tmp_path)["visual_issues"]


def test_expiry_during_fsync_suppresses_commit_and_keeps_core(tmp_path, monkeypatch):
    docs = real_graph_bundle()
    save(tmp_path, docs)
    before = {name: (tmp_path / name).read_bytes() for name in docs}
    now = [0.0]
    calls = []

    def delayed(fd):
        now[0] = 6.0

    monkeypatch.setattr(p.os, "fsync", delayed)
    monkeypatch.setattr(p.os, "replace", lambda *a: calls.append(a))
    result = p.publish_review_visuals(
        tmp_path,
        docs["manifest.json"]["hit"],
        docs["results.json"]["results"],
        docs["refinement.json"]["candidates"],
        clock=lambda: now[0],
    )
    assert result["status"] == "unavailable" and result["deadline_overrun"]
    assert not calls and before == {name: (tmp_path / name).read_bytes() for name in docs}


def test_exact_limits_and_missing_entries_cannot_be_silently_extended(visual_docs, tmp_path):
    visual_docs[d.ARTIFACT]["limits"]["atoms"] += 1
    rebind(visual_docs)
    save(tmp_path, visual_docs)
    assert report.build_report(tmp_path)["visual_issues"]


def test_coordinates_beyond_bound_rejected_without_work(monkeypatch):
    docs = real_graph_bundle()
    trial = docs["refinement.json"]["candidates"][0]["trials"][0]
    pose = trial["selective_pose_evidence"]["hit_candidate"]
    pose["coordinates"][0][0] = 1e100
    pose["sha256"] = d.binding({k: v for k, v in pose.items() if k != "sha256"})
    trial["similarity_constraints"]["witness"]["candidate_pose_sha256"] = pose["sha256"]
    with pytest.raises(ValueError, match="coordinate bounds"):
        p.saved_pair(trial, docs["manifest.json"]["hit"], docs["results.json"]["results"][0]["smiles"])


@pytest.mark.parametrize("kind", [b"zTXt", b"iTXt", b"iCCP"])
def test_compressed_ancillary_payloads_are_not_passed_to_browser(kind):
    import base64

    data = base64.b64decode(p.depiction("CCO")["png_base64"])
    with pytest.raises(ValueError, match="Unsupported PNG"):
        d.checked_png(data[:33] + chunk(kind, b"opaque-compressed-metadata") + data[33:])


def test_complete_artifact_limit_fails_without_partial_advertisement(tmp_path, monkeypatch):
    docs = real_graph_bundle()
    docs["manifest.json"]["review_visuals"] = {
        "version": d.VERSION,
        **p.unavailable("optional_display_not_published"),
    }
    save(tmp_path, docs)
    before = {name: (tmp_path / name).read_bytes() for name in docs}
    monkeypatch.setitem(p.LIMITS, "artifact_bytes", 1000)
    result = p.publish_review_visuals(
        tmp_path,
        docs["manifest.json"]["hit"],
        docs["results.json"]["results"],
        docs["refinement.json"]["candidates"],
    )
    assert result["status"] == "unavailable" and "byte limit" in result["reason"]
    assert not (tmp_path / d.ARTIFACT).exists()
    assert before == {name: (tmp_path / name).read_bytes() for name in docs}


def test_view_render_escapes_caller_supplied_visual_labels(visual_docs, tmp_path):
    save(tmp_path, visual_docs)
    view = report.build_report(tmp_path)
    view["candidates"][0]["visuals"]["seeds"][0]["seed"] = "</option><img src=x onerror=evil()>"
    page = report.render_html(view)
    assert "<img src=x" not in page and "&lt;/option&gt;" in page


def test_deferred_optional_file_mutation_is_detected(visual_docs, tmp_path, monkeypatch):
    save(tmp_path, visual_docs)
    original = report._assemble

    def mutate(*args):
        result = original(*args)
        (tmp_path / d.ARTIFACT).write_text("{}")
        return result

    monkeypatch.setattr(report, "_assemble", mutate)
    with pytest.raises(ValueError, match="changed"):
        report.write_report(tmp_path)
    assert not (tmp_path / "review.html").exists()


def test_png_declared_dimensions_and_inflated_raster_bound():
    header = struct.pack(">IIBBBBB", 1025, 1, 8, 2, 0, 0, 0)
    blob = (
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", header)
        + chunk(b"IDAT", zlib.compress(b"\0" * 3076))
        + chunk(b"IEND", b"")
    )
    with pytest.raises(ValueError, match="dimensions"):
        d.checked_png(blob)
    header = struct.pack(">IIBBBBB", 1, 1, 8, 2, 0, 0, 0)
    blob = (
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", header)
        + chunk(b"IDAT", zlib.compress(b"\0" * 10000))
        + chunk(b"IEND", b"")
    )
    with pytest.raises(ValueError, match="decompressed size"):
        d.checked_png(blob)


def test_failed_retained_trial_coverage_is_explicit_not_omitted(tmp_path):
    docs = real_graph_bundle()
    record = docs["refinement.json"]["candidates"][0]
    record["status"] = "failed"
    record["trials"][1] = {"seed": 43, "status": "failed", "reason": "controlled absent pose"}
    attach(docs)
    save(tmp_path, docs)
    view = report.build_report(tmp_path)
    assert view["artifact_consistency"] == "incomplete_or_inconsistent"
    assert not view["visual_issues"]
    seeds = view["candidates"][0]["visuals"]["seeds"]
    assert [s["seed"] for s in seeds] == [42, 43, 44]
    assert seeds[1]["pair"]["status"] == "unavailable"
    assert "retained_trial_not_successful" in seeds[1]["pair"]["reason"]


@pytest.mark.parametrize("also_bad_metric", [False, True])
def test_corrupt_compressed_png_remains_independent_display_diagnostic(
    visual_docs, tmp_path, also_bad_metric
):
    import base64

    image = visual_docs[d.ARTIFACT]["candidates"][0]["depictions"]["original_hit"]
    raw = base64.b64decode(image["png_base64"])
    malformed = raw[:33] + chunk(b"IDAT", b"not-a-zlib-stream") + chunk(b"IEND", b"")
    image["png_base64"] = base64.b64encode(malformed).decode()
    image["png_sha256"] = hashlib.sha256(malformed).hexdigest()
    if also_bad_metric:
        visual_docs["results.json"]["results"][0]["hit_shape"] = 0.0
    rebind(visual_docs)
    save(tmp_path, visual_docs)
    view = report.write_report(tmp_path)
    assert view["visual_issues"] and view["candidates"][0]["visuals"]["status"] == "inconsistent"
    if also_bad_metric:
        assert (
            view["candidates"][0]["issues"] and view["artifact_consistency"] == "incomplete_or_inconsistent"
        )
    else:
        assert "compressed image data" in view["visual_issues"][0]
        assert view["artifact_consistency"] == "consistent_recorded_bundle"
