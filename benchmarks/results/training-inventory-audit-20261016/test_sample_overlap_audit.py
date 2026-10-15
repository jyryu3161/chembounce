"""Small audit-contract regressions, not benchmark scientific evidence."""
import argparse
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("sample_overlap_audit", Path(__file__).with_name("sample_overlap_audit.py"))
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class AuditTests(unittest.TestCase):
    def test_exact_and_family_distinct_and_duplicates(self):
        rows = [{"source": "x", "role": "query", "smiles": s} for s in
                ("C[C@H](O)C(=O)O", "C[C@@H](O)C(=O)O", "C[C@H](O)C(=O)O", "not_smiles")]
        result = audit.compare(["C[C@H](O)C(=O)O", "", "C[C@H](O)C(=O)O"], rows, audit.make_normalizer())
        total = next(x for x in result["groups"] if x["source"] == x["role"] == "ALL")
        self.assertEqual(total["exact"]["distinct"], 2)
        self.assertEqual(total["family"]["distinct"], 1)
        self.assertEqual(total["exact"]["matching_rows"], 2)
        self.assertEqual(total["family"]["matching_rows"], 3)
        self.assertEqual(total["exact"]["unassessed_rows"], 1)
        self.assertEqual(result["sample"]["exact"]["duplicate_rows"], 1)
        self.assertEqual(result["sample"]["empty_lines"], 1)

    def test_family_failure_does_not_erase_exact(self):
        def normalize(smiles):
            return {"exact": smiles, "family": None, "exact_error": None, "family_error": "limit"}
        result = audit.compare(["C"], [{"source": "x", "role": "query", "smiles": "C"}], normalize)
        total = next(x for x in result["groups"] if x["source"] == x["role"] == "ALL")
        self.assertEqual(total["exact"]["matching_distinct"], 1)
        self.assertEqual(total["family"]["unassessed_rows"], 1)
        self.assertEqual(len(result["sample_normalization_failures"]), 1)

    def test_panel_and_controls_fail_closed_on_wrong_link(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            panel = {"compounds": [{"id": "a", "target_id": "T", "label": "active", "split": "test", "smiles": "CC"}],
                     "queries": [{"id": "q", "compound_id": "a"}]}
            (root / "panel.json").write_text(json.dumps(panel))
            (root / "controls.json").write_text(json.dumps([{"query_id": "q", "compound_id": "a", "target_id": "WRONG", "label": "active"}]))
            sources = [{"id": "x", "kind": "controls", "path": "controls.json", "panel": "panel.json"}]
            with self.assertRaisesRegex(ValueError, "target lookup"):
                audit.project_sources(sources, root)

    def test_snapshot_pins_source_bytes_and_preserves_failed_run_candidates(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            payload = json.dumps([{"method_id": "m", "query_id": "q", "target_id": "T", "seed": 42,
                "status": "timeout", "response": {"candidates": [{"smiles": "CC"}]}, "evaluation": {"candidates": []}}]).encode()
            (root / "runs.json").write_bytes(payload)
            result = audit.project_sources([{"id": "x", "kind": "runs", "path": "runs.json"}], root)
            self.assertEqual(result["sources"][0]["sha256"], hashlib.sha256(payload).hexdigest())
            self.assertEqual(result["rows"][0]["status"], "timeout")
            self.assertEqual(result["rows"][0]["role"], "returned_candidate")

    def test_archive_hash_and_regular_member_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            archive = root / "sample.tar.gz"
            with tarfile.open(archive, "w:gz") as stream:
                member = tarfile.TarInfo("data/train.smi")
                member.size = 3
                stream.addfile(member, io.BytesIO(b"CC\n"))
            args = argparse.Namespace(archive=archive, archive_sha256="bad", member="data/train.smi", source_url="https://example.invalid/source", output=root / "sample.smi")
            with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
                audit.extract(args)
            args.archive_sha256 = audit.digest(archive)
            audit.extract(args)
            self.assertEqual(Path(args.output).read_bytes(), b"CC\n")
            with self.assertRaisesRegex(ValueError, "overwrite"):
                audit.extract(args)


if __name__ == "__main__":
    unittest.main()
