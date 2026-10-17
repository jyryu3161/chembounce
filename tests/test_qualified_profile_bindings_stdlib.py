"""Packaged artifact/binding checks independent of synthetic runtime fixtures.

Execute only the adapter's exact binding functions plus stdlib-only profiles.py.
This authenticates packaged bytes, not scientific runtime compatibility or efficacy.
Runnable with unittest without importing the scientific package.
"""

import ast
import hashlib
import importlib.util
from pathlib import Path
import sys
from types import ModuleType
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


class PackagedBenchmarkBindings(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location(
            "_benchmark_profiles_stdlib", ROOT / "chembounce_synthesis/profiles.py"
        )
        cls.profiles = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.profiles)
        tree = ast.parse((ROOT / "benchmarks/research/qualified_adapter.py").read_text())
        nodes = [
            node
            for node in tree.body
            if (
                isinstance(node, ast.FunctionDef)
                and node.name in {"benchmark_profile", "selective_benchmark_profile"}
            )
            or (
                isinstance(node, ast.Assign)
                and any(
                    isinstance(target, ast.Name) and target.id == "PROFILE_PROTOCOLS"
                    for target in node.targets
                )
            )
        ]
        cls.namespace = {}
        exec(
            compile(ast.Module(body=nodes, type_ignores=[]), "qualified_adapter.py:bindings", "exec"),
            cls.namespace,
        )

    def call(self, name, *args):
        package = ModuleType("chembounce_synthesis")
        package.profiles = self.profiles
        with patch.dict(
            sys.modules, {"chembounce_synthesis": package, "chembounce_synthesis.profiles": self.profiles}
        ):
            return self.namespace[name](*args)

    def test_packaged_v2_and_v5_have_exact_distinct_authorities(self):
        for version, name in (
            ("v2", self.profiles.SELECTIVE_PROFILE),
            ("v3", self.profiles.LANDMARK_SELECTIVE_PROFILE),
        ):
            with self.subTest(version=version):
                raw = (ROOT / "chembounce_synthesis/data" / (name + ".json")).read_bytes()
                digest = hashlib.sha256(raw).hexdigest()
                self.assertEqual(digest, self.profiles.PROFILE_SHA256[name])
                binding = {
                    "schema": "qualified-benchmark-protocol-" + version,
                    "recommendation_profile": name,
                    "profile_artifact_sha256": digest,
                }
                self.assertEqual(self.call("benchmark_profile", binding), name)
                self.assertTrue(self.call("selective_benchmark_profile", name))
                for field, value in (
                    ("schema", "qualified-benchmark-protocol-v4"),
                    ("profile_artifact_sha256", "0" * 64),
                    ("recommendation_profile", self.profiles.DIAGNOSTIC_SELECTIVE_PROFILE),
                    ("extra", True),
                ):
                    with self.subTest(field=field), self.assertRaises(ValueError):
                        self.call("benchmark_profile", {**binding, field: value})

    def test_no_default_upgrade_and_no_intermediate_profiles(self):
        with patch.object(self.profiles, "DEFAULT_RECOMMENDATION_PROFILE", "future-profile"):
            self.assertEqual(self.call("benchmark_profile"), self.profiles.PROFILE)
        self.assertFalse(self.call("selective_benchmark_profile", self.profiles.PROFILE))
        for name in (
            self.profiles.GUARDED_SELECTIVE_PROFILE,
            self.profiles.DIAGNOSTIC_SELECTIVE_PROFILE,
            "future-profile",
        ):
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.call("selective_benchmark_profile", name)


if __name__ == "__main__":
    unittest.main()
