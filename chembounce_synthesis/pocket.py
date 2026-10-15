"""Optional geometric checks in an explicitly verified bound-reference frame.

Distances, van der Waals overlaps and user-specified anchors are geometric
observations. They are not docking energies, binding affinity or activity.
"""

import copy
import math
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

from .common import canonical, digest
from .conformers import FEATURES

PROTEIN_RESIDUES = set(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL SEC PYL".split()
)
WATER_RESIDUES = {"HOH", "WAT", "DOD"}


def _number(value, name, *, positive=False):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value < 0
        or (positive and value == 0)
    ):
        raise ValueError(f"pocket.{name} must be a finite {'positive' if positive else 'nonnegative'} number")


def _selector(value):
    required = {"chain", "residue_number", "residue_name"}
    if (
        not isinstance(value, dict)
        or not required <= set(value)
        or set(value) - required - {"insertion_code"}
    ):
        raise ValueError(
            "Residue selectors require chain, residue_number, residue_name and optional insertion_code"
        )
    if not isinstance(value["chain"], str) or len(value["chain"]) > 1:
        raise ValueError("PDB chain identifiers must be zero or one characters")
    if isinstance(value["residue_number"], bool) or not isinstance(value["residue_number"], int):
        raise ValueError("Residue number must be an integer")
    if not isinstance(value["residue_name"], str) or not 1 <= len(value["residue_name"]) <= 3:
        raise ValueError("PDB residue names require one to three characters")
    code = value.get("insertion_code", "")
    if not isinstance(code, str) or len(code) > 1:
        raise ValueError("PDB insertion code must be zero or one characters")
    return (value["chain"], value["residue_number"], code, value["residue_name"])


def validate_pocket_config(config):
    allowed = {
        "receptor_pdb",
        "reference_ligand",
        "model",
        "reference_altloc",
        "coordinate_tolerance",
        "include_waters",
        "include_hetero",
        "protein_chains",
        "include_residues",
        "exclude_residues",
        "vdw_overlap_tolerance",
        "max_clashes",
        "contact_distance",
        "max_reported_pairs",
        "anchors",
    }
    if not isinstance(config, dict) or set(config) - allowed:
        raise ValueError("pocket must be an object with supported settings")
    if not config.get("receptor_pdb") and not config.get("anchors"):
        raise ValueError("pocket requires receptor_pdb or explicit anchors")
    if "receptor_pdb" in config:
        if not isinstance(config["receptor_pdb"], str) or not config["receptor_pdb"].strip():
            raise ValueError("pocket.receptor_pdb must be a path")
        _selector(config.get("reference_ligand"))
    elif "reference_ligand" in config:
        raise ValueError("reference_ligand requires receptor_pdb")
    for key in ("include_waters", "include_hetero"):
        if key in config and not isinstance(config[key], bool):
            raise ValueError(f"pocket.{key} must be true or false")
    for key in ("model", "max_reported_pairs"):
        if key in config and (
            isinstance(config[key], bool) or not isinstance(config[key], int) or config[key] < 1
        ):
            raise ValueError(f"pocket.{key} must be a positive integer")
    if config.get("max_clashes") is not None and (
        isinstance(config["max_clashes"], bool)
        or not isinstance(config["max_clashes"], int)
        or config["max_clashes"] < 0
    ):
        raise ValueError("pocket.max_clashes must be a nonnegative integer")
    for key in ("coordinate_tolerance", "contact_distance", "vdw_overlap_tolerance"):
        if key in config:
            _number(config[key], key, positive=key != "vdw_overlap_tolerance")
    for key in ("include_residues", "exclude_residues"):
        if key in config:
            if not isinstance(config[key], list):
                raise ValueError(f"pocket.{key} must be a list")
            for selector in config[key]:
                _selector(selector)
    if "protein_chains" in config and (
        not isinstance(config["protein_chains"], list)
        or any(not isinstance(chain, str) or len(chain) > 1 for chain in config["protein_chains"])
    ):
        raise ValueError("pocket.protein_chains must be a list of PDB chain IDs")
    if config.get("reference_altloc") is not None and (
        not isinstance(config["reference_altloc"], str) or len(config["reference_altloc"]) != 1
    ):
        raise ValueError("reference_altloc must be one explicit alternate-location character")
    anchors = config.get("anchors", [])
    if not isinstance(anchors, list) or len(anchors) > 128:
        raise ValueError("pocket.anchors must be a list of at most 128 explicit anchors")
    names = set()
    for anchor in anchors:
        if not isinstance(anchor, dict) or not isinstance(anchor.get("id"), str) or not anchor["id"]:
            raise ValueError("Every pocket anchor requires a nonempty id")
        if anchor["id"] in names:
            raise ValueError("Pocket anchor IDs must be unique")
        names.add(anchor["id"])
        kind = anchor.get("kind")
        fields = {"id", "kind", "max_distance", "required"}
        if kind == "atom":
            fields |= {"reference_atom_index", "candidate_smarts", "candidate_atom_index"}
        elif kind == "feature":
            fields |= {"reference_feature_index"}
        else:
            raise ValueError("Pocket anchor kind must be atom or feature")
        if set(anchor) != fields:
            raise ValueError(f"Anchor {anchor['id']} requires exactly {sorted(fields)}")
        _number(anchor["max_distance"], "anchor.max_distance", positive=True)
        if not isinstance(anchor["required"], bool):
            raise ValueError("Anchor required must be true or false")
        index_names = (
            ("reference_atom_index", "candidate_atom_index")
            if kind == "atom"
            else ("reference_feature_index",)
        )
        for key in index_names:
            if isinstance(anchor[key], bool) or not isinstance(anchor[key], int) or anchor[key] < 0:
                raise ValueError(f"Anchor {key} must be a nonnegative integer")
        if kind == "atom":
            query = (
                Chem.MolFromSmarts(anchor["candidate_smarts"])
                if isinstance(anchor["candidate_smarts"], str)
                else None
            )
            if query is None or anchor["candidate_atom_index"] >= query.GetNumAtoms():
                raise ValueError("Anchor candidate SMARTS/atom index is invalid")
    return copy.deepcopy(config)


def _residue_dict(key):
    return dict(zip(("chain", "residue_number", "insertion_code", "residue_name"), key))


def _read_pdb(path, model, reference_selector, reference_altloc):
    residues, current_model, explicit_models = {}, 1, False
    table = Chem.GetPeriodicTable()
    for line in path.read_text().splitlines():
        record = line[:6].strip()
        if record == "MODEL":
            explicit_models = True
            current_model = int(line[10:14].strip())
            continue
        if record not in {"ATOM", "HETATM"} or current_model != model:
            continue
        try:
            residue = (line[21:22].strip(), int(line[22:26]), line[26:27].strip(), line[17:20].strip())
            element = line[76:78].strip()
            if not element:
                if record == "ATOM" and residue[3] in PROTEIN_RESIDUES:
                    element = line[12:16].strip().lstrip("0123456789")[0]
                else:
                    raise ValueError("Explicit element columns are required for nonprotein PDB atoms")
            number = table.GetAtomicNumber(element.capitalize())
            xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            occupancy = float(line[54:60].strip() or "1")
            if not all(math.isfinite(x) for x in xyz) or not math.isfinite(occupancy) or occupancy < 0:
                raise ValueError("Nonfinite coordinates or invalid occupancy")
            if number <= 1 or occupancy == 0:
                continue
            atom = {
                "serial": int(line[6:11]),
                "atom_name": line[12:16].strip(),
                "altloc": line[16:17].strip(),
                "record": record,
                "residue": residue,
                "element": element.capitalize(),
                "atomic_number": number,
                "position": xyz,
                "occupancy": occupancy,
            }
        except (ValueError, RuntimeError, IndexError) as exc:
            raise ValueError(f"Invalid PDB atom record: {line[:30]!r}: {exc}") from exc
        residues.setdefault(residue, []).append(atom)
    if not explicit_models and model != 1:
        raise ValueError("PDB without MODEL records only supports model=1")
    if not residues:
        raise ValueError("Selected PDB model contains no heavy atoms")
    atoms, altlocs = [], []
    for residue, entries in residues.items():
        alternatives = {atom["altloc"] for atom in entries} - {""}
        if residue == reference_selector and reference_altloc is not None:
            if reference_altloc not in alternatives:
                raise ValueError("Requested reference ligand alternate location is absent")
            selected = reference_altloc
        else:
            # Select a coherent residue alternate, not a hybrid of per-atom choices.
            selected = min(
                alternatives,
                key=lambda alt: (
                    -sum(a["occupancy"] for a in entries if a["altloc"] == alt),
                    alt != "A",
                    alt,
                ),
                default="",
            )
        if selected:
            altlocs.append({**_residue_dict(residue), "selected_altloc": selected})
        by_name = {}
        for atom in sorted(entries, key=lambda a: (-a["occupancy"], a["altloc"] != "", a["serial"])):
            if atom["altloc"] in {"", selected}:
                by_name.setdefault(atom["atom_name"], atom)
        atoms.extend(by_name.values())
    return atoms, altlocs


def _feature_records(mol, cid):
    return [
        {
            "feature_index": i,
            "family": feature.GetFamily(),
            "atom_indices": list(feature.GetAtomIds()),
            "position": list(feature.GetPos()),
        }
        for i, feature in enumerate(FEATURES.GetFeaturesForMol(mol, confId=cid))
    ]


class BoundPoseConstraints:
    def __init__(self, config, hit_pose, reference_sdf):
        self.config = validate_pocket_config(config)
        if reference_sdf is None or hit_pose is None or not hit_pose.bound_pose:
            raise ValueError("Pocket constraints require an explicit bound reference_sdf and its fixed pose")
        source = Path(reference_sdf).resolve()
        supplied = list(Chem.SDMolSupplier(str(source), removeHs=False))
        if len(supplied) != 1 or supplied[0] is None:
            raise ValueError("Pocket reference SDF requires exactly one valid molecule")
        self.reference = Chem.Mol(hit_pose.mol)
        self.reference_cid = next(iter(hit_pose.relative_energies))
        original = Chem.RemoveHs(supplied[0])
        if canonical(Chem.MolToSmiles(original)) != canonical(Chem.MolToSmiles(self.reference)):
            raise ValueError("Pocket reference file and scoring reference identities differ")
        if not np.allclose(
            original.GetConformer().GetPositions(),
            self.reference.GetConformer(self.reference_cid).GetPositions(),
            atol=1e-8,
            rtol=0,
        ):
            raise ValueError("Pocket reference scoring coordinates differ from the supplied SDF frame")
        self.reference_features = _feature_records(self.reference, self.reference_cid)
        self.anchors = copy.deepcopy(config.get("anchors", []))
        ref_positions = self.reference.GetConformer(self.reference_cid).GetPositions()
        for anchor in self.anchors:
            if anchor["kind"] == "atom":
                index = anchor["reference_atom_index"]
                if (
                    index >= self.reference.GetNumAtoms()
                    or self.reference.GetAtomWithIdx(index).GetAtomicNum() <= 1
                ):
                    raise ValueError("Reference atom anchors require a valid heavy-atom index")
                anchor["position"] = ref_positions[index].tolist()
            else:
                index = anchor["reference_feature_index"]
                if index >= len(self.reference_features):
                    raise ValueError("Reference feature anchor index is absent")
                feature = self.reference_features[index]
                anchor.update(position=feature["position"], family=feature["family"])
        self.receptor, self.tree = [], None
        self.provenance = {
            "reference_sdf": str(source),
            "reference_sdf_sha256": digest(source),
            "reference_features": self.reference_features,
            "settings": copy.deepcopy(config),
            "effective_settings": {
                "model": 1,
                "reference_altloc": None,
                "coordinate_tolerance": 0.25,
                "include_waters": False,
                "include_hetero": False,
                "protein_chains": None,
                "include_residues": [],
                "exclude_residues": [],
                "vdw_overlap_tolerance": 0.4,
                "max_clashes": None,
                "contact_distance": 4.0,
                "max_reported_pairs": 1000,
                "anchors": [],
                **copy.deepcopy(config),
            },
            "rdkit_version": rdBase.rdkitVersion,
            "coordinate_frame": "unchanged_supplied_bound_reference",
            "interpretation": "Heavy-atom geometry and explicit anchors; not binding affinity or activity",
        }
        if config.get("receptor_pdb"):
            self._load_receptor(config)

    def _load_receptor(self, config):
        path = Path(config["receptor_pdb"]).resolve()
        selector = _selector(config["reference_ligand"])
        atoms, altlocs = _read_pdb(path, config.get("model", 1), selector, config.get("reference_altloc"))
        bound = [atom for atom in atoms if atom["residue"] == selector]
        heavy = [atom.GetIdx() for atom in self.reference.GetAtoms() if atom.GetAtomicNum() > 1]
        if len(bound) != len(heavy):
            raise ValueError("Selected PDB ligand and reference SDF have different heavy-atom counts")
        positions = self.reference.GetConformer(self.reference_cid).GetPositions()[heavy]
        target = np.array([atom["position"] for atom in bound])
        costs = np.linalg.norm(positions[:, None, :] - target[None, :, :], axis=2)
        elements = np.array([self.reference.GetAtomWithIdx(index).GetAtomicNum() for index in heavy])
        mismatch = elements[:, None] != np.array([atom["atomic_number"] for atom in bound])[None, :]
        costs[mismatch] = 1e9
        rows, cols = linear_sum_assignment(costs)
        distance = float(costs[rows, cols].max())
        if distance > config.get("coordinate_tolerance", 0.25):
            raise ValueError(
                f"Reference SDF/PDB ligand coordinate frames do not match (maximum distance {distance:.3f} A)"
            )
        includes = {_selector(value) for value in config.get("include_residues", [])}
        excludes = {_selector(value) for value in config.get("exclude_residues", [])} | {selector}
        chains = config.get("protein_chains")
        counts = {
            "reference_or_excluded": 0,
            "water": 0,
            "other_hetero": 0,
            "protein": 0,
            "included_explicitly": 0,
        }
        for atom in atoms:
            residue = atom["residue"]
            if residue in excludes:
                counts["reference_or_excluded"] += 1
                continue
            if residue in includes:
                selected, category = True, "included_explicitly"
            elif residue[3] in WATER_RESIDUES:
                selected, category = config.get("include_waters", False), "water"
            elif atom["record"] == "ATOM" and residue[3] in PROTEIN_RESIDUES:
                selected, category = chains is None or residue[0] in chains, "protein"
            else:
                selected, category = (
                    config.get("include_hetero", False) and atom["record"] == "HETATM",
                    "other_hetero",
                )
            if selected:
                self.receptor.append(atom)
                counts[category] += 1
        if not self.receptor:
            raise ValueError("Explicit receptor filters selected no heavy atoms")
        self.tree = cKDTree([atom["position"] for atom in self.receptor])
        table = Chem.GetPeriodicTable()
        self.receptor_radii = np.array([table.GetRvdw(atom["atomic_number"]) for atom in self.receptor])
        self.provenance.update(
            receptor_pdb=str(path),
            receptor_pdb_sha256=digest(path),
            receptor_heavy_atoms=len(self.receptor),
            selection_counts=counts,
            selected_altlocs=altlocs,
            frame_max_atom_distance=distance,
            altloc_policy="coherent per-residue maximum summed occupancy; ties prefer A then lexical",
            receptor_filters={
                "model": config.get("model", 1),
                "include_waters": config.get("include_waters", False),
                "include_hetero": config.get("include_hetero", False),
                "protein_chains": chains,
                "include_residues": config.get("include_residues", []),
                "exclude_residues": config.get("exclude_residues", []),
                "reference_ligand_automatically_excluded": config["reference_ligand"],
            },
        )

    def validate_reference(self, hit_pose):
        if hit_pose is None or not hit_pose.bound_pose:
            raise ValueError("Pocket constraints cannot be used with a freely generated hit ensemble")
        if (
            canonical(Chem.MolToSmiles(hit_pose.mol)) != canonical(Chem.MolToSmiles(self.reference))
            or len(hit_pose.relative_energies) != 1
        ):
            raise ValueError("Pocket constraints require the same single bound-reference chemical state")
        cid = next(iter(hit_pose.relative_energies))
        if not np.allclose(
            hit_pose.mol.GetConformer(cid).GetPositions(),
            self.reference.GetConformer(self.reference_cid).GetPositions(),
            atol=1e-8,
            rtol=0,
        ):
            raise ValueError("Pocket constraints require the original fixed bound-reference coordinates")

    def _anchor_assessment(self, mol, cid):
        if not self.anchors:
            return []
        positions = mol.GetConformer(cid).GetPositions()
        features, possibilities, candidates = _feature_records(mol, cid), [], {}
        for anchor in self.anchors:
            options = {}
            if anchor["kind"] == "atom":
                query = Chem.MolFromSmarts(anchor["candidate_smarts"])
                matches = mol.GetSubstructMatches(query, useChirality=True, maxMatches=10001)
                if len(matches) > 10000:
                    raise ValueError("Pocket atom-anchor matching exceeded its explicit 10000-match bound")
                for match in matches:
                    index = match[anchor["candidate_atom_index"]]
                    if mol.GetAtomWithIdx(index).GetAtomicNum() > 1:
                        options[("atom", index)] = (positions[index], [index])
            else:
                for feature in features:
                    if feature["family"] == anchor["family"]:
                        options[("feature", feature["family"], tuple(feature["atom_indices"]))] = (
                            feature["position"],
                            feature["atom_indices"],
                        )
            candidates.update(options)
            possibilities.append(options)
        identifiers = list(candidates)
        costs = np.full((len(self.anchors), len(identifiers) + len(self.anchors)), 1e9)
        distances = []
        for row, (anchor, options) in enumerate(zip(self.anchors, possibilities)):
            values = {}
            for col, identifier in enumerate(identifiers):
                if identifier in options:
                    distance = float(np.linalg.norm(np.array(options[identifier][0]) - anchor["position"]))
                    values[identifier] = distance
                    if distance <= anchor["max_distance"]:
                        costs[row, col] = distance / anchor["max_distance"]
            costs[row, len(identifiers) + row] = 1e6 if anchor["required"] else 1000
            distances.append(values)
        rows, cols = linear_sum_assignment(costs)
        reports = []
        for row, col in zip(rows, cols):
            anchor = self.anchors[row]
            matched = col < len(identifiers) and costs[row, col] < 1e3
            identifier = identifiers[col] if matched else None
            reports.append(
                {
                    "id": anchor["id"],
                    "kind": anchor["kind"],
                    "required": anchor["required"],
                    "max_distance": anchor["max_distance"],
                    "reference_position": anchor["position"],
                    "matched": bool(matched),
                    "distance": distances[row][identifier] if matched else None,
                    "nearest_distance": min(distances[row].values(), default=None),
                    "candidate_atom_indices": candidates[identifier][1] if matched else [],
                    "passed": bool(matched or not anchor["required"]),
                }
            )
        return reports

    def evaluate(self, mol, cid):
        positions = mol.GetConformer(cid).GetPositions()
        if not np.isfinite(positions).all():
            raise ValueError("Pocket assessment received nonfinite candidate coordinates")
        clashes, contacts, clash_count, max_overlap = [], {}, 0, 0.0
        tolerance = self.config.get("vdw_overlap_tolerance", 0.4)
        contact_distance = self.config.get("contact_distance", 4.0)
        if self.tree is not None:
            table = Chem.GetPeriodicTable()
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() <= 1:
                    continue
                index, radius = atom.GetIdx(), table.GetRvdw(atom.GetAtomicNum())
                cutoff = max(contact_distance, radius + float(self.receptor_radii.max()))
                for other in self.tree.query_ball_point(positions[index], cutoff):
                    receptor = self.receptor[other]
                    distance = float(np.linalg.norm(positions[index] - receptor["position"]))
                    overlap = float(radius + self.receptor_radii[other] - distance)
                    if distance <= contact_distance:
                        residue = receptor["residue"]
                        contacts[residue] = min(contacts.get(residue, float("inf")), distance)
                    max_overlap = max(max_overlap, overlap)
                    if overlap > tolerance:
                        clash_count += 1
                        if len(clashes) < self.config.get("max_reported_pairs", 1000):
                            clashes.append(
                                {
                                    "candidate_atom_index": index,
                                    "receptor_atom_serial": receptor["serial"],
                                    "receptor_atom_name": receptor["atom_name"],
                                    "residue": _residue_dict(receptor["residue"]),
                                    "distance": distance,
                                    "vdw_overlap": overlap,
                                }
                            )
        anchors = self._anchor_assessment(mol, cid)
        reasons = [
            f"Required anchor {anchor['id']} not satisfied" for anchor in anchors if not anchor["passed"]
        ]
        if self.config.get("max_clashes") is not None and clash_count > self.config["max_clashes"]:
            reasons.append(
                f"Heavy-atom clash pairs {clash_count} exceed max_clashes={self.config['max_clashes']}"
            )
        return {
            "passed": not reasons,
            "reasons": reasons,
            "candidate_conformer_id": cid,
            "clashes": {
                "assessed": self.tree is not None,
                "pair_count": clash_count,
                "pairs": clashes,
                "reported_pairs_truncated": clash_count > len(clashes),
                "maximum_vdw_overlap": max_overlap,
                "vdw_overlap_tolerance": tolerance,
                "radii_source": "RDKit periodic-table van der Waals radii",
            },
            "contacts": {
                "distance_cutoff": contact_distance,
                "residues": [
                    {**_residue_dict(key), "minimum_distance": value}
                    for key, value in sorted(contacts.items())
                ],
            },
            "anchors": anchors,
            "coordinate_frame": self.provenance["coordinate_frame"],
        }


def constraints_from_config(config, hit_pose, reference_sdf=None):
    return BoundPoseConstraints(config, hit_pose, reference_sdf)
