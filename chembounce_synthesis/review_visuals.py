"""Optional new-run display producer. No embedding, alignment or scientific scoring.

The finite allowance is cooperative: native calls and persistence can overrun.
Core files are already complete; only an atomic optional manifest extension is
made. Known expiry prevents starting a commit; a commit returning late can
remain content-available, with final persistence timing explicitly unverified.
"""

import base64
import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time

from .review_display import ARTIFACT, LIMITS, VERSION, binding, checked_graph, checked_png, require


def unavailable(reason):
    return {"status": "unavailable", "reason": str(reason)[:1024] or "display_unavailable"}


def _deadline(deadline, clock):
    if clock() >= deadline:
        raise TimeoutError("display_deadline_exhausted")


def graph_record(mol, smiles):
    """Capture exact indexed graph; topology intentionally does not include stereo."""
    from rdkit import Chem
    from .common import canonical, absolute_smiles, require_absolute_stereo
    from .similarity_constraints import indexed_topology

    require_absolute_stereo(mol)
    require(canonical(absolute_smiles(mol, isomericSmiles=True)) == smiles, "Display exact isomeric identity")
    require(
        mol.GetNumAtoms() <= LIMITS["atoms"] and mol.GetNumBonds() <= LIMITS["bonds"], "Display graph limit"
    )
    # Stereo assignment is on a defensive copy; never change a caller or saved pose.
    mol = Chem.Mol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=False, force=True)
    value = {
        "smiles": smiles,
        "atom_order_sha256": indexed_topology(mol),
        "atoms": [
            {
                "index": a.GetIdx(),
                "atomic_number": a.GetAtomicNum(),
                "element": a.GetSymbol(),
                "isotope": a.GetIsotope(),
                "formal_charge": a.GetFormalCharge(),
                "chiral_tag": str(a.GetChiralTag()),
                "cip": a.GetProp("_CIPCode") if a.HasProp("_CIPCode") else None,
            }
            for a in mol.GetAtoms()
        ],
        "bonds": [
            {
                "begin": min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                "end": max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                "type": str(b.GetBondType()),
                "stereo": str(b.GetStereo()),
                "stereo_atoms": list(b.GetStereoAtoms()),
            }
            for b in mol.GetBonds()
        ],
    }
    value["sha256"] = binding(value)
    checked_graph(value)
    return value


def depiction(smiles):
    """Real static PNG from the exact structure; 2D work only, on a private copy."""
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D
    from .common import molecule, canonical, absolute_smiles

    mol = molecule(smiles)
    graph = graph_record(mol, smiles)
    draw = Chem.Mol(mol)
    draw.RemoveAllConformers()
    # PrepareMolForDrawing creates only depiction coordinates, no 3D inference.
    draw = rdMolDraw2D.PrepareMolForDrawing(draw, kekulize=True, addChiralHs=False)
    require(
        canonical(absolute_smiles(draw, isomericSmiles=True)) == smiles, "Depiction changed exact identity"
    )
    drawer = rdMolDraw2D.MolDraw2DCairo(480, 320)
    drawer.drawOptions().includeMetadata = False
    drawer.DrawMolecule(draw)
    drawer.FinishDrawing()
    data = bytes(drawer.GetDrawingText())
    width, height = checked_png(data)
    require(
        canonical(absolute_smiles(draw, isomericSmiles=True)) == smiles, "Depiction changed exact identity"
    )
    return {
        "status": "available",
        "smiles": smiles,
        "graph": graph,
        "png_base64": base64.b64encode(data).decode("ascii"),
        "png_sha256": hashlib.sha256(data).hexdigest(),
        "width": width,
        "height": height,
    }


def saved_pair(trial, hit, product):
    from .similarity_constraints import restore_pose

    require(trial.get("status") == "success", "retained_trial_not_successful")
    witness = trial["similarity_constraints"]["witness"]
    reference = witness["reference_pose"]
    candidate = trial["selective_pose_evidence"]["hit_candidate"]
    require(
        reference["smiles"] == hit
        and candidate["smiles"] == product
        and witness["candidate_pose_sha256"] == candidate["sha256"],
        "Saved pair identity/witness mismatch",
    )
    graphs = []
    for saved in (reference, candidate):
        require(len(saved["coordinates"]) <= LIMITS["atoms"], "Pose atom limit")
        graph = graph_record(restore_pose(saved), saved["smiles"])
        checked_graph(graph, saved)
        graphs.append(graph)
    return {
        "status": "available",
        "reference_graph": graphs[0],
        "candidate_graph": graphs[1],
        "reference_pose_sha256": reference["sha256"],
        "candidate_pose_sha256": candidate["sha256"],
        "witness_sha256": binding(witness),
    }


def build_visuals(hit, rows, refinements, *, deadline_monotonic=None, clock=time.monotonic):
    """Build complete bounded coverage or raise; input records are never changed."""
    start = clock()
    deadline = (
        min(start + LIMITS["allowance_seconds"], deadline_monotonic)
        if deadline_monotonic is not None
        else start + LIMITS["allowance_seconds"]
    )
    _deadline(deadline, clock)
    require(len(rows) <= LIMITS["candidates"], "Display candidate limit")
    from rdkit import rdBase

    entries = []
    for row in rows:
        _deadline(deadline, clock)
        records = [
            r for r in refinements if r.get("product") == row["smiles"] and r.get("design") == row["design"]
        ]
        require(len(records) <= 1, "Ambiguous display refinement")
        trials = records[0]["trials"] if records else []
        require(len(trials) <= LIMITS["seeds_per_candidate"], "Display retained seed limit")
        entry = {
            "key": [row["route_id"], row["smiles"], row["design"]],
            "row_sha256": binding(row),
            "representative_seed": row.get("pose_seed"),
            "depictions": {},
            "seeds": [],
        }
        for role, smiles in (("original_hit", hit), ("final_product", row["smiles"])):
            try:
                _deadline(deadline, clock)
                image = depiction(smiles)
                _deadline(deadline, clock)
                entry["depictions"][role] = image
            except Exception as exc:
                entry["depictions"][role] = unavailable(f"{type(exc).__name__}: {exc}")
        for trial in trials:
            try:
                _deadline(deadline, clock)
                pair = saved_pair(trial, hit, row["smiles"])
                _deadline(deadline, clock)
            except Exception as exc:
                pair = unavailable(f"{type(exc).__name__}: {exc}")
            entry["seeds"].append({"seed": trial["seed"], "trial_sha256": binding(trial), "pair": pair})
        entries.append(entry)
        require(
            len(json.dumps(entries, allow_nan=False).encode()) <= LIMITS["artifact_bytes"] - 4096,
            "Display artifact byte limit",
        )
    _deadline(deadline, clock)
    return {
        "schema_version": 1,
        "version": VERSION,
        "limits": dict(LIMITS),
        "hit": hit,
        "producer": {
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "rdkit_version": rdBase.rdkitVersion,
        },
        "generation_elapsed_seconds": clock() - start,
        "candidates": entries,
    }


def _atomic(path, data, deadline, clock):
    """Same-directory replace; check after fsync and immediately before replace."""
    _deadline(deadline, clock)
    fd, name = tempfile.mkstemp(prefix=".review-display-", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        _deadline(deadline, clock)
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def publish_review_visuals(out, hit, rows, refinements, *, deadline_monotonic=None, clock=time.monotonic):
    """Best-effort optional export after core writes. Returns separate display cost.

    Existing manifest has an unavailable default. If work fails or expires before commit,
    no available descriptor is installed. A complete unadvertised orphan is safe.
    Available means complete content, not verified budget success.
    An already-started filesystem replace cannot be retroactively interrupted;
    any overrun is returned honestly. No scientific result is removed or retried.
    """
    out = Path(out)
    start = clock()
    deadline = (
        min(start + LIMITS["allowance_seconds"], deadline_monotonic)
        if deadline_monotonic is not None
        else start + LIMITS["allowance_seconds"]
    )
    status = unavailable("optional_display_not_published")
    try:
        _deadline(deadline, clock)
        manifest_path = out / "manifest.json"
        original = manifest_path.read_bytes()
        manifest = json.loads(original)
        # Core output durability precedes optional work; no core byte changes.
        for filename in (
            "manifest.json",
            "results.json",
            "routes.json",
            "scores.json",
            "refinement.json",
            "qualification.json",
            "summary.json",
        ):
            _deadline(deadline, clock)
            with (out / filename).open("rb") as stream:
                os.fsync(stream.fileno())
        _deadline(deadline, clock)
        directory_fd = os.open(out, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
        artifact = build_visuals(hit, rows, refinements, deadline_monotonic=deadline, clock=clock)
        data = (json.dumps(artifact, sort_keys=True, allow_nan=False) + "\n").encode()
        require(len(data) <= LIMITS["artifact_bytes"], "Display artifact byte limit")
        _deadline(deadline, clock)
        require(not (out / ARTIFACT).exists(), "Display artifact already exists")
        _atomic(out / ARTIFACT, data, deadline, clock)
        _deadline(deadline, clock)
        require(manifest_path.read_bytes() == original, "Core manifest changed during display")
        updated = copy.deepcopy(manifest)
        updated["review_visuals"] = {
            "version": VERSION,
            "status": "available",
            "path": ARTIFACT,
            "sha256": hashlib.sha256(data).hexdigest(),
            "publication_completion_within_budget": None,
            "timing_scope": "final_persistence_unverified",
            "prepublication_elapsed_seconds": clock() - start,
        }
        _atomic(
            manifest_path,
            (json.dumps(updated, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode(),
            deadline,
            clock,
        )
        status = {"status": "available"}
    except Exception as exc:
        status = unavailable(f"{type(exc).__name__}: {exc}")
    return {**status, "elapsed_seconds": clock() - start, "deadline_overrun": clock() > deadline}
