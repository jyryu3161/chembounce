#!/usr/bin/env python3
"""Run original 1.0 attachment code; requires the legacy Python dependencies.

The RDBASE compatibility shim permits the old module to import on pip RDKit.
It does not modify the tagged source. No external scaffold database is required.
"""

import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
LEGACY_REF = "1.0-original"
BOOTSTRAP = r"""
import os, json, rdkit
os.environ['RDBASE'] = os.path.dirname(rdkit.__file__)
import chembounce as cb
import utils
from rdkit import Chem
smiles = 'C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N'
mol, canonical = utils.init_mol(smiles)
report = {'rdkit': rdkit.__version__, 'input': smiles, 'compatibility_shim': 'RDBASE set to pip rdkit location', 'fragments': []}
for frag in sorted(cb.sg.get_all_murcko_fragments(mol, False, 1), key=Chem.MolToSmiles):
    row = {'scaffold': Chem.MolToSmiles(frag)}
    try:
        row['remainders'] = cb.remove_substructures(mol, frag)
    except Exception as exc:
        row.update(error_type=type(exc).__name__, message=str(exc))
    report['fragments'].append(row)
props = utils._get_molecular_prop_extd_(mol)
report['properties'] = props
report['default_filter_failures'] = utils.determ_props(props, utils._default_thrs_(), True)[1]
print(json.dumps(report))
"""


def main():
    with tempfile.TemporaryDirectory() as directory:
        names = subprocess.check_output(
            ["git", "ls-tree", "-r", "--name-only", LEGACY_REF], cwd=ROOT, text=True
        ).splitlines()
        for name in names:
            if name.endswith(".py") and "/" not in name:
                content = subprocess.check_output(["git", "show", LEGACY_REF + ":" + name], cwd=ROOT)
                (Path(directory) / name).write_bytes(content)
        output = subprocess.check_output([sys.executable, "-c", BOOTSTRAP], cwd=directory, text=True)
    report = json.loads(output)
    report["tag_commit"] = subprocess.check_output(
        ["git", "rev-list", "-n", "1", LEGACY_REF], cwd=ROOT, text=True
    ).strip()
    path = Path(__file__).with_name("legacy_adp.json")
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
