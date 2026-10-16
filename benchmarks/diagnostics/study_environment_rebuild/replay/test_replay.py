"""Synthetic filesystem/report fixtures only; never run scientific helpers."""

import contextlib
import copy
import importlib.util
import io
import json
import math
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("study_replay_compare", HERE / "compare.py")
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)


class ReplayTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.assets = self.root / "fresh"
        self.reference = copy.deepcopy(compare.load_reference("research"))
        self.report = self.fixture("research", self.reference)

    def fixture(self, profile, reference):
        expected, _, pins, _ = reference
        source = self.assets / "sources/chembounce_current"
        venv = self.assets / "venvs" / profile
        venv.mkdir(parents=True, exist_ok=True)
        hashes = pins["sources"]["chembounce_current"]["expected_files_sha256"]
        for relative in hashes:
            file = source / relative
            file.parent.mkdir(parents=True, exist_ok=True)
            file.write_text("Synthetic source fixture: " + relative)
            hashes[relative] = compare.sha(file)
        report = copy.deepcopy(expected)
        for name in expected["module_origins"]:
            if name in compare.PROJECT_FILES:
                file = source / compare.PROJECT_FILES[name]
            elif name in compare.UPSTREAM_FILES:
                file = self.assets / "sources" / (profile + "_upstream") / compare.UPSTREAM_FILES[name]
                file.parent.mkdir(parents=True, exist_ok=True)
                file.write_text("Synthetic upstream fixture: " + name)
            else:
                file = venv / "lib/python3.10/site-packages" / name / "__init__.py"
                file.parent.mkdir(parents=True, exist_ok=True)
                file.write_text("Synthetic site fixture: " + name)
            digest = compare.sha(file)
            expected["module_origins"][name] = digest
            report["module_origins"][name] = {"path": str(file), "sha256": digest}
        return report

    def check(self, report=None, reference=None, profile="research"):
        return compare.compare_report(self.report if report is None else report, profile,
                                      self.assets, self.reference if reference is None else reference)

    def test_all_four_profiles_synthetic_origin_layout(self):
        for profile, count in (("research", 8), ("forward-cpu", 8), ("crem", 7), ("mol2mol", 9)):
            with self.subTest(profile=profile):
                reference = copy.deepcopy(compare.load_reference(profile))
                report = self.fixture(profile, reference)
                result = self.check(report, reference, profile)
                self.assertEqual(result["module_origin_checks"], count)
                self.assertEqual(result["conformer_scalar_checks"], 58)
                self.assertFalse(result["probe_execution_proven"])

    def test_required_fields_and_extra_field(self):
        for field in self.report:
            with self.subTest(field=field):
                value = copy.deepcopy(self.report)
                del value[field]
                with self.assertRaises(ValueError):
                    self.check(value)
        self.report["unrecognized"] = True
        with self.assertRaises(ValueError):
            self.check()

    def test_exact_typed_values(self):
        mutations = [
            lambda x: x["identities"][0].update(has_scaffold=1),
            lambda x: x["fingerprints"]["query_bit_similarities"].__setitem__(0, 1),
            lambda x: x["conformer"].update(seed=float(x["conformer"]["seed"])),
            lambda x: x["conformer"].update(minimize_status=False),
            lambda x: x.update(python="Different Python build"),
            lambda x: x.update(rdkit_version="2026.09.1"),
            lambda x: x.update(profile="mol2mol"),
            lambda x: x["reaction_products"].append("extra"),
            lambda x: x["fingerprints"]["count_sparse"][0].update(extra=1),
        ]
        for index, mutation in enumerate(mutations):
            with self.subTest(mutation=index):
                report = copy.deepcopy(self.report)
                mutation(report)
                with self.assertRaises(ValueError):
                    self.check(report)

    def test_exact_small_float_difference_rejected(self):
        self.report["fingerprints"]["query_count_similarities"][0] -= 1e-12
        with self.assertRaises(ValueError):
            self.check()

    def test_tolerance_inclusive_boundary_and_next_float(self):
        self.reference[0]["conformer"]["energy"] = 0.0
        for value in (compare.TOLERANCE, -compare.TOLERANCE):
            self.report["conformer"]["energy"] = value
            self.assertEqual(self.check()["maximum_absolute_conformer_difference"], compare.TOLERANCE)
        self.report["conformer"]["energy"] = math.nextafter(compare.TOLERANCE, math.inf)
        with self.assertRaises(ValueError):
            self.check()

    def test_nonfinite_anywhere_rejected(self):
        for number in (math.nan, math.inf, -math.inf):
            for location in ("energy", "position", "fingerprint"):
                with self.subTest(number=number, location=location):
                    report = copy.deepcopy(self.report)
                    if location == "energy":
                        report["conformer"]["energy"] = number
                    elif location == "position":
                        report["conformer"]["positions"][0][0] = number
                    else:
                        report["fingerprints"]["query_count_similarities"][0] = number
                    with self.assertRaises(ValueError):
                        self.check(report)

    def test_58_values_and_shapes_required(self):
        for mutation in (
            lambda x: x["positions"].pop(),
            lambda x: x["positions"][0].pop(),
            lambda x: x.update(energy=True),
            lambda x: x["positions"][0].__setitem__(0, "0.0"),
        ):
            report = copy.deepcopy(self.report)
            mutation(report["conformer"])
            with self.assertRaises(ValueError):
                self.check(report)

    def test_module_membership_and_fields(self):
        for mutation in (
            lambda x: x.pop("numpy"),
            lambda x: x.update(extra=x["numpy"]),
            lambda x: x["numpy"].pop("sha256"),
            lambda x: x["numpy"].update(extra=True),
            lambda x: x["numpy"].update(path="numpy/__init__.py"),
            lambda x: x["numpy"].update(sha256="0" * 64),
        ):
            report = copy.deepcopy(self.report)
            mutation(report["module_origins"])
            with self.assertRaises(ValueError):
                self.check(report)

    def test_old_source_fallback_rejected_even_with_equal_bytes(self):
        old = self.assets / "sources/chembounce_frozen/chembounce_synthesis/identity.py"
        old.parent.mkdir(parents=True)
        old.write_bytes(Path(self.report["module_origins"]["chembounce_synthesis.identity"]["path"]).read_bytes())
        self.report["module_origins"]["chembounce_synthesis.identity"]["path"] = str(old)
        with self.assertRaises(ValueError):
            self.check()

    def test_old_venv_fallback_rejected_even_with_equal_bytes(self):
        old = self.assets / "venvs/original-research/lib/python3.10/site-packages/numpy/__init__.py"
        old.parent.mkdir(parents=True)
        old.write_bytes(Path(self.report["module_origins"]["numpy"]["path"]).read_bytes())
        self.report["module_origins"]["numpy"]["path"] = str(old)
        with self.assertRaises(ValueError):
            self.check()

    def test_wrong_relative_origin_inside_declared_root(self):
        file = Path(self.report["module_origins"]["numpy"]["path"])
        wrong = file.with_name("wrong.py")
        wrong.write_bytes(file.read_bytes())
        self.report["module_origins"]["numpy"]["path"] = str(wrong)
        with self.assertRaises(ValueError):
            self.check()

    def test_missing_or_changed_module_bytes(self):
        file = Path(self.report["module_origins"]["numpy"]["path"])
        file.write_text("tampered")
        with self.assertRaises(ValueError):
            self.check()
        file.unlink()
        with self.assertRaises(OSError):
            self.check()

    def test_all_41_source_files_required(self):
        file = self.assets / "sources/chembounce_current/benchmarks/research/controls.py"
        self.assertTrue(file.is_file())
        file.write_text("tampered unimported scientific source")
        with self.assertRaises(ValueError):
            self.check()

    def test_source_file_symlink_escape(self):
        file = self.assets / "sources/chembounce_current/chemistry.py"
        outside = self.root / "old-chemistry.py"
        outside.write_bytes(file.read_bytes())
        file.unlink()
        file.symlink_to(outside)
        with self.assertRaises(ValueError):
            self.check()

    def test_json_duplicate_and_nonfinite_tokens(self):
        file = self.root / "bad.json"
        for data in ('{"a":1,"a":2}', '{"a":NaN}', '{"a":Infinity}', '{"a":1e999}'):
            file.write_text(data)
            with self.assertRaises(ValueError):
                compare.read_json(file)

    def copied_artifacts(self):
        target = self.root / "artifacts/replay"
        shutil.copytree(HERE, target)
        shutil.copyfile(HERE.parent / "pins.json", target.parent / "pins.json")
        return target

    def test_reference_artifact_tampering_and_missing(self):
        target = self.copied_artifacts()
        for relative in ("provenance.json", "probe_scientific.py", "audit_reinvent_installed.py",
                         "criteria.json", "expected/research.json", "../pins.json"):
            with self.subTest(file=relative):
                file = target / relative
                original = file.read_bytes()
                file.write_bytes(original + b" ")
                with self.assertRaises(ValueError):
                    compare.load_reference("research", target)
                file.write_bytes(original)
        (target / "expected/crem.json").unlink()
        with self.assertRaises(OSError):
            compare.load_reference("research", target)

    def invoke(self, output, report=None):
        file = self.root / "report.json"
        file.write_text(json.dumps(self.report if report is None else report))
        with patch.object(compare, "load_reference", return_value=self.reference), contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return compare.main(["--profile", "research", "--assets-root", str(self.assets),
                                 "--report", str(file), "--output", str(output)])

    def test_cli_pass_and_exclusive_output_preservation(self):
        output = self.root / "result.json"
        self.assertEqual(self.invoke(output), 0)
        original = output.read_bytes()
        self.assertEqual(json.loads(original)["status"], "passed_saved_report_checks")
        self.assertEqual(self.invoke(output), 2)
        self.assertEqual(output.read_bytes(), original)

    def test_cli_failure_is_recorded_and_not_overwritten(self):
        output = self.root / "failure.json"
        self.report.pop("identities")
        self.assertEqual(self.invoke(output), 1)
        original = output.read_bytes()
        self.assertEqual(json.loads(original)["status"], "failed")
        self.assertEqual(self.invoke(output), 2)
        self.assertEqual(output.read_bytes(), original)

    def test_cli_missing_input_records_failure(self):
        output = self.root / "missing.json"
        with contextlib.redirect_stdout(io.StringIO()):
            result = compare.main(["--profile", "research", "--assets-root", str(self.assets),
                                   "--report", str(self.root / "absent.json"), "--output", str(output)])
        self.assertEqual(result, 1)
        self.assertEqual(json.loads(output.read_text())["status"], "failed")

    def test_cli_dangling_symlink_output_preserved(self):
        output = self.root / "result-link.json"
        missing = self.root / "must-not-create.json"
        output.symlink_to(missing)
        self.assertEqual(self.invoke(output), 2)
        self.assertTrue(output.is_symlink())
        self.assertFalse(missing.exists())


if __name__ == "__main__":
    unittest.main()
