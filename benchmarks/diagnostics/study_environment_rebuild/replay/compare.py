"""Validate saved fixed CPU probe reports; this does not execute the probe."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PROVENANCE_SHA256 = "a3a353110dff10433c1b329db1d041e2ccfe60bce560d949d63df7b1a1c51dc6"
PROFILES = ("research", "forward-cpu", "crem", "mol2mol")
TOLERANCE = 1e-10
PROJECT_FILES = {
    "chembounce_synthesis.identity": "chembounce_synthesis/identity.py",
    "chembounce_synthesis.common": "chembounce_synthesis/common.py",
    "chemistry": "chemistry.py",
}
UPSTREAM_FILES = {
    "crem": "crem/__init__.py",
    "crem.crem": "crem/crem.py",
    "reinvent": "reinvent/__init__.py",
    "reinvent.models.transformer.mol2mol.mol2mol": "reinvent/models/transformer/mol2mol/mol2mol.py",
    "reinvent.models.transformer.core.dataset.dataset": "reinvent/models/transformer/core/dataset/dataset.py",
}
SITE_MODULES = {"numpy", "rdkit", "torch", "scipy", "transformers"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "Duplicate JSON key: " + key)
        result[key] = value
    return result


def finite_tree(value):
    if type(value) is float:
        require(math.isfinite(value), "Nonfinite JSON numeric value")
    elif type(value) is dict:
        for item in value.values():
            finite_tree(item)
    elif type(value) is list:
        for item in value:
            finite_tree(item)


def read_json(path):
    raw = Path(path).read_bytes()
    value = json.loads(raw, object_pairs_hook=unique_object)
    finite_tree(value)  # Reject NaN/Infinity and finite-looking JSON such as 1e999.
    return value, hashlib.sha256(raw).hexdigest()


def typed_equal(actual, expected, location):
    require(type(actual) is type(expected), "Type mismatch: " + location)
    if type(expected) is dict:
        require(actual.keys() == expected.keys(), "Field membership mismatch: " + location)
        for key in expected:
            typed_equal(actual[key], expected[key], location + "." + key)
    elif type(expected) is list:
        require(len(actual) == len(expected), "Length mismatch: " + location)
        for index, (left, right) in enumerate(zip(actual, expected)):
            typed_equal(left, right, f"{location}[{index}]")
    else:
        require(actual == expected, "Value mismatch: " + location)


def load_reference(profile, directory=HERE):
    require(profile in PROFILES, "Unknown profile")
    directory = Path(directory)
    provenance, digest = read_json(directory / "provenance.json")
    require(digest == PROVENANCE_SHA256, "Changed replay provenance")
    files = {
        **provenance["helpers_sha256"],
        **provenance["expected_snapshots_sha256"],
        "criteria.json": provenance["criteria_sha256"],
    }
    for relative, expected in files.items():
        require(sha(directory / relative) == expected, "Changed replay artifact: " + relative)
    pins, digest = read_json(directory.parent / "pins.json")
    require(digest == provenance["source_pins_sha256"], "Changed source pins")
    criteria, _ = read_json(directory / "criteria.json")
    require(
        criteria["comparison"]["conformer_absolute_tolerance"] == TOLERANCE, "Unexpected comparison tolerance"
    )
    require(criteria["comparison"]["conformer_relative_tolerance"] == 0, "Unexpected relative tolerance")
    expected, _ = read_json(directory / "expected" / (profile + ".json"))
    return expected, criteria, pins, provenance


def rooted(path, root, label):
    resolved = Path(path).resolve(strict=True)
    require(resolved.is_relative_to(root), "Origin outside declared root: " + label)
    return resolved


def compare_report(report, profile, assets_root, reference):
    expected, criteria, pins, provenance = reference
    require(type(report) is dict, "Probe report must be an object")
    finite_tree(report)
    require(report.keys() == expected.keys(), "Probe report field membership mismatch")
    require(report["profile"] == profile, "Wrong report profile")
    for key in criteria["comparison"]["exact_fields"]:
        typed_equal(report[key], expected[key], key)
    conformer, baseline = report["conformer"], expected["conformer"]
    require(
        type(conformer) is dict and conformer.keys() == baseline.keys(), "Conformer field membership mismatch"
    )
    for key in criteria["comparison"]["conformer_exact_fields"]:
        typed_equal(conformer[key], baseline[key], "conformer." + key)
    positions = conformer["positions"]
    require(type(positions) is list and len(positions) == 19, "Expected 19 coordinate rows")
    require(
        all(type(row) is list and len(row) == 3 for row in positions), "Expected three coordinates per atom"
    )
    values = [conformer["energy"], *[x for row in positions for x in row]]
    targets = [baseline["energy"], *[x for row in baseline["positions"] for x in row]]
    require(
        len(values) == 58 and all(type(v) is float and math.isfinite(v) for v in values),
        "Expected 58 finite floating-point conformer values",
    )
    differences = [abs(a - b) for a, b in zip(values, targets)]
    require(max(differences) <= TOLERANCE, "Conformer absolute tolerance exceeded")

    assets = Path(assets_root).resolve(strict=True)
    require(assets.is_dir(), "Assets root must be a directory")
    source = rooted(assets / "sources/chembounce_current", assets, "study source")
    venv = rooted(assets / "venvs" / profile, assets, "profile environment")
    source_spec = pins["sources"]["chembounce_current"]
    require(source_spec["commit"] == provenance["execution_source_commit"], "Wrong source commit pin")
    require(len(source_spec["expected_files_sha256"]) == 41, "Expected all 41 study source pins")
    for relative, digest in source_spec["expected_files_sha256"].items():
        path = rooted(source / relative, source, relative)
        require(sha(path) == digest, "Changed study source: " + relative)
    origins = report["module_origins"]
    require(
        type(origins) is dict and origins.keys() == expected["module_origins"].keys(),
        "Module membership mismatch",
    )
    for name, expected_hash in expected["module_origins"].items():
        entry = origins[name]
        require(type(entry) is dict and set(entry) == {"path", "sha256"}, "Invalid origin: " + name)
        require(
            type(entry["path"]) is str and Path(entry["path"]).is_absolute(),
            "Origin must be an absolute file path: " + name,
        )
        if name in PROJECT_FILES:
            origin_root, relative = source, PROJECT_FILES[name]
        elif name in UPSTREAM_FILES:
            require(profile in ("crem", "mol2mol"), "Unexpected upstream module")
            origin_root = rooted(assets / "sources" / (profile + "_upstream"), assets, name)
            relative = UPSTREAM_FILES[name]
        else:
            require(name in SITE_MODULES, "Unexpected site module: " + name)
            origin_root, relative = venv, "lib/python3.10/site-packages/" + name + "/__init__.py"
        target = rooted(origin_root / relative, origin_root, name)
        actual = rooted(entry["path"], origin_root, name)
        require(actual == target, "Wrong module file: " + name)
        require(
            type(entry["sha256"]) is str and entry["sha256"] == expected_hash,
            "Reported module hash mismatch: " + name,
        )
        require(sha(actual) == expected_hash, "Current module bytes mismatch: " + name)
    return {
        "status": "passed_saved_report_checks",
        "profile": profile,
        "conformer_scalar_checks": 58,
        "maximum_absolute_conformer_difference": max(differences),
        "absolute_tolerance": TOLERANCE,
        "module_origin_checks": len(origins),
        "study_source_files_checked": 41,
        "assets_root": str(assets),
        "probe_execution_proven": False,
        "scope": "Saved fixed CPU report checked against its own profile snapshot and current declared files; no probe execution, independent OS or general numerical-equivalence claim.",
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=PROFILES, required=True)
    parser.add_argument("--assets-root", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    # Exclusive creation also rejects existing files and dangling symlinks.
    try:
        output = args.output.open("x")
    except OSError as exc:
        print("Cannot create new result: " + str(exc), file=sys.stderr)
        return 2
    with output:
        try:
            reference = load_reference(args.profile)
            report, report_hash = read_json(args.report)
            result = compare_report(report, args.profile, args.assets_root, reference)
            result.update(
                report_sha256=report_hash,
                comparator_sha256=sha(__file__),
                provenance_sha256=PROVENANCE_SHA256,
            )
            code = 0
        except Exception as exc:
            result = {
                "status": "failed",
                "profile": args.profile,
                "error": type(exc).__name__ + ": " + str(exc),
                "probe_execution_proven": False,
            }
            code = 1
        result["checked_utc"] = datetime.now(timezone.utc).isoformat()
        output.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"status": result["status"], "output": str(args.output)}))
    return code


if __name__ == "__main__":
    raise SystemExit(main())
