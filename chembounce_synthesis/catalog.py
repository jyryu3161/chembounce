"""Versioned building-block snapshots; no live stock claims or pickle loading."""

import csv
import json
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path

from rdkit import Chem, DataStructs

from .common import FP, canonical, digest, molecule, positive


@dataclass(frozen=True)
class BuildingBlock:
    id: str
    smiles: str
    source: str
    snapshot: str
    availability: str = "catalogued"


class Catalog:
    def __init__(self, filename):
        self.path = Path(filename).resolve()
        self.sha256 = digest(self.path)
        self.import_provenance = None
        manifest_path = self.path.with_suffix(self.path.suffix + ".manifest.json")
        if manifest_path.is_file():
            manifest = json.loads(manifest_path.read_text())
            if manifest.get("output_sha256") != self.sha256:
                raise ValueError("Catalog checksum differs from its import manifest; regenerate the snapshot")
            self.import_provenance = {
                "path": str(manifest_path),
                "sha256": digest(manifest_path),
                "inputs": manifest.get("inputs", []),
                "records": manifest.get("records"),
                "rejected_record_count": len(manifest.get("rejected_records", [])),
                "id_policy": manifest.get("id_policy"),
                "availability_policy": manifest.get("availability_policy"),
            }
        self.blocks = {}
        self.by_smiles = {}
        self.mols, self.fps = {}, {}
        with self.path.open(newline="", encoding="utf-8-sig") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if not {"id", "smiles", "source", "snapshot"}.issubset(reader.fieldnames or []):
                raise ValueError("Catalog TSV requires id, smiles, source, snapshot columns")
            for line, row in enumerate(reader, 2):
                if any(
                    not isinstance(row.get(key), str) or not row[key].strip()
                    for key in ("id", "source", "snapshot")
                ):
                    raise ValueError(f"Empty catalog provenance at line {line}")
                ident = row["id"].strip()
                if ident in self.blocks:
                    raise ValueError(f"Duplicate building-block ID: {ident}")
                mol = molecule(row["smiles"])
                smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                status = (row.get("availability") or "catalogued").strip()
                if status not in {"catalogued", "in_stock", "unavailable"}:
                    raise ValueError(f"Unknown availability: {status}")
                block = BuildingBlock(ident, smiles, row["source"], row["snapshot"], status)
                self.blocks[ident] = block
                if status != "unavailable":
                    self.by_smiles.setdefault(smiles, []).append(ident)
                    if smiles not in self.mols:
                        # Converted snapshots already use canonical atom order.
                        # Reparse only noncanonical input; avoid parsing every
                        # molecule twice for a large, normalized catalog.
                        if row["smiles"] != smiles:
                            mol = molecule(smiles)
                        self.mols[smiles] = mol
                        self.fps[smiles] = FP.GetFingerprint(mol)
        if not self.mols:
            raise ValueError("Catalog contains no eligible building blocks")
        # A cache on the class method would keep entire discarded catalogs alive
        # through its self keys. Own the bounded cache within each instance.
        self.compatible = lru_cache(maxsize=512)(self.compatible)

    def get(self, ident):
        block = self.blocks.get(ident)
        if block is None or block.availability == "unavailable":
            raise ValueError(f"Building block is absent or unavailable: {ident}")
        return block

    def compatible(self, role_smarts):
        query = Chem.MolFromSmarts(role_smarts)
        if query is None:
            raise ValueError("Invalid reaction-role SMARTS")
        return tuple(s for s, m in self.mols.items() if m.HasSubstructMatch(query, useChirality=True))

    def nearest(self, smiles, role_smarts, k=8):
        """Reactant-role filtering then Morgan search with a diverse tail.

        This is a bounded shortlist, not a global 3D nearest-neighbour claim.
        Exact matches are always first; duplicate vendor SKUs do not occupy slots.
        """
        positive(k, "building_block_top_k")
        query = canonical(smiles)
        eligible = self.compatible(role_smarts)
        if not eligible:
            return []
        sims = DataStructs.BulkTanimotoSimilarity(
            FP.GetFingerprint(molecule(query)), [self.fps[s] for s in eligible]
        )
        ordered = sorted(zip(eligible, sims), key=lambda x: (x[0] != query, -x[1], x[0]))
        close_count = max(1, (3 * k) // 4)
        chosen = [s for s, _ in ordered[:close_count]]
        pool = [s for s, _ in ordered[close_count : max(256, 8 * k)]]
        while len(chosen) < k and pool:
            # Max-min diversity in the coarse shortlist, deterministic tie breaking.
            item = max(
                pool,
                key=lambda s: (
                    1 - max(DataStructs.BulkTanimotoSimilarity(self.fps[s], [self.fps[c] for c in chosen])),
                    s,
                ),
            )
            chosen.append(item)
            pool.remove(item)
        return [self.blocks[sorted(self.by_smiles[s])[0]] for s in chosen]

    def evidence(self, ident):
        return asdict(self.get(ident)) | {"catalog_sha256": self.sha256}
