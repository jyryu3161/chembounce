"""Disabled telemetry structure for reviewed current source, with exact94b history.

Canonical AST retains every semantic field. Only empty type_params on the three
Python3.12 declaration nodes are omitted, matching their absence on3.10/3.11.
Nonempty type parameters and unknown fields remain compared. Historical94b
digests below are retained as attribution; the reviewed P1 v1 digest is kept
as history, while current P1 v2 digests bind the new retrieval source. No new
scientific call is normalized away.
"""

import ast
import hashlib
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BASE = "94b5c0149f2025de7d6e48e01244c8f1257339a5"


def dump(node):
    return ast.dump(node, include_attributes=False)


def expression(text):
    return dump(ast.parse(text, mode="eval").body)


def statement(text):
    return dump(ast.parse(text).body[0])


OBSERVATION_IFS = {
    statement(
        'if telemetry is not None and summary["deadline_reached"]:\n    telemetry_phase.interrupt("workflow_deadline_observed")'
    ),
    statement('if telemetry is not None:\n    telemetry.skip("common.refinement")'),
    statement(
        'if telemetry is not None:\n    for name in ("generation", "allocation"):\n        telemetry.unavailable(f"native.{name}", "legacy_native_path_not_instrumented")'
    ),
    statement(
        'if "telemetry" in options:\n    raise ValueError("Telemetry is process-local, not a generator_options callback")'
    ),
}
FORWARD = expression('{"telemetry": telemetry} if telemetry is not None else {}')


class DisabledTelemetry(ast.NodeTransformer):
    """Remove exact approved observation shapes; leave mutations detectable."""

    def visit_ImportFrom(self, node):
        if node.module in {"phase_telemetry", "chembounce_synthesis.phase_telemetry"}:
            assert [(n.name, n.asname) for n in node.names] == [
                ("phase_span", None),
                ("validate_telemetry", None),
            ]
            return None
        return node

    def visit_FunctionDef(self, node):
        node = self.generic_visit(node)
        if node.args.kwonlyargs and node.args.kwonlyargs[-1].arg == "telemetry":
            assert dump(node.args.kw_defaults[-1]) == dump(ast.Constant(value=None))
            node.args.kwonlyargs.pop()
            node.args.kw_defaults.pop()
        return node

    def visit_Expr(self, node):
        if dump(node) in {
            statement('validate_telemetry(telemetry, "native")'),
            statement('validate_telemetry(telemetry, "common")'),
        }:
            return None
        return self.generic_visit(node)

    def visit_With(self, node):
        if len(node.items) == 1:
            item = node.items[0]
            call = item.context_expr
            if (
                isinstance(call, ast.Call)
                and isinstance(call.func, ast.Name)
                and call.func.id == "phase_span"
            ):
                assert not call.keywords and len(call.args) == 2
                assert dump(call.args[0]) == dump(ast.Name(id="telemetry", ctx=ast.Load()))
                assert isinstance(call.args[1], ast.Constant) and isinstance(call.args[1].value, str)
                assert item.optional_vars is None or dump(item.optional_vars) == dump(
                    ast.Name(id="telemetry_phase", ctx=ast.Store())
                )
                return self.generic_visit(node).body
        return self.generic_visit(node)

    def visit_If(self, node):
        if dump(node) in OBSERVATION_IFS:
            return None
        return self.generic_visit(node)

    def visit_Call(self, node):
        node.keywords = [k for k in node.keywords if not (k.arg is None and dump(k.value) == FORWARD)]
        return self.generic_visit(node)


def canonical(node):
    if isinstance(node, ast.AST):
        fields = []
        for name, value in ast.iter_fields(node):
            if (
                name == "type_params"
                and isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
                and value == []
            ):
                continue
            fields.append([name, canonical(value)])
        return [type(node).__name__, fields]
    if isinstance(node, list):
        return [canonical(value) for value in node]
    return [type(node).__name__, repr(node)]


def structural_digest(source):
    tree = DisabledTelemetry().visit(ast.parse(source))
    return hashlib.sha256(
        json.dumps(canonical(tree), separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()


# Baseline byte SHA256 and canonical AST SHA256 generated from git show BASE:path
# during authoring. Tests need neither Git nor ancestor objects/shallow history.
BASELINES = {
    "chembounce.py": {
        "source_sha256": "daaa7c955be76a89e75fc1f4527b64e889255bdca836509d048570ade2188299",
        "ast_sha256": "0ecfe5751d5292211ab11f215ef67d622ee2b885fbb2fe99df2d094e4abf5bbd",
    },
    "chembounce_synthesis/workflow.py": {
        "source_sha256": "866cb7c9f741a953028f8847354300fe201aaaa34735855d9e77ab35afa74605",
        "ast_sha256": "a93dcc699c86f62d2523316c36397ec0cb4fe336bcc0b324362b34df9f9b4501",
    },
    "benchmarks/research/qualified_adapter.py": {
        "source_sha256": "acc6b411bb973800779527fc54fbb9180ad5b2669f68c231c56adb7f67f372d8",
        "ast_sha256": "1071adec654a21e3b16e3da989ce07e6dffeecb3e5ae86ead1ffdc94e9cd0314",
    },
}

# Approved feature source is a new semantic baseline for the two changed files.
# The qualified adapter retains its exact94b disabled-telemetry AST pin; its
# instrumented source bytes have always differed from that pre-telemetry pin.
ROUTE_FIRST_V1_REVIEWED = {
    "commit": "0454d8643bc03b061565439a4ad77649c114114d",
    "scope": "Reviewed opt-in P1 route-first generation; not exact94b semantic parity",
    "files": {
        "chembounce.py": {
            "source_sha256": "b1bc441565b248add217f2402c791538cbe31d193440acd781ab8569f7f2c14f",
            "ast_sha256": "c9984d40f8eeed5250b15007b9a8992dffcbea8e3d3b14c8a5445fbf470cb872",
        },
        "chembounce_synthesis/workflow.py": {
            "source_sha256": "0e8082aacede68db0898456d677065c0c05a40f64eaba933a37d040b44f41b2c",
            "ast_sha256": "ce54ab6f45ac9f63923397112387e96f7de65308796093fff1ec511dcbe88da6",
        },
    },
}

CURRENT_REVIEWED = {
    "commit": "826b2af90d75f5243851395250fe53e96224f356",
    "scope": "Reviewed opt-in P1 v2 catalog retrieval and bounded stereo-anchor postings; not exact94b parity",
    "files": {
        "chembounce.py": {
            "source_sha256": "04bf2033dc16a720074c63fb5ab16e80a2ef7168b9dfe5e49562b3684aa19dff",
            "ast_sha256": "f1d90fc41f97e3b565442a4c210deda3e77c918ed6c25d89803778b522f71159",
        },
        "chembounce_synthesis/workflow.py": {
            "source_sha256": "bed15cb634abc67e9c9bfa8ddf578b8e26076d8308e1d74c19eb3d68d929df1e",
            "ast_sha256": "1cd153a3e57c579160a025633a417cae78006e0423a2d2346a27594c606de30e",
        },
    },
}

MAPPED_NH_REVIEWED = {
    "base_commit": "c558d70b29729072d1fed3f1614a5e52fd53d5af",
    "scope": "One reviewed workflow manifest source-provenance tuple addition for route_first_nh.py; preceding P1 v2 source pin remains historical",
    "workflow_source_sha256": "ebc033185ba5bd8dedbb67712e2d0e9ac1a45391c3d9996e7d8133e468b74ed9",
    "workflow_ast_sha256": "78df6c387a5be9f9ef732e9118248adc7a492bb939ee170180352ea9a6573a98",
}

MARKED_SECOND_REVIEWED = {
    "base_commit": "e6cf009",
    "scope": "The route-first manifest adds only two marked-site helper source hashes; prior source pins remain historical",
    "workflow_source_sha256": "79580260de8c2132ecb60049c0716b8200ba14e8e428047cc1c94033d903ef81",
    "workflow_ast_sha256": "f2fb8737c77c71f41731c2d1aec90b12ee69f51e25601a952217743768ab5d11",
}


def without_marked_manifest_sources(source):
    """Remove only the two newly reviewed helper names from the manifest tuple."""
    new = expression(
        '{name: _SOURCE_HASHES[name] for name in ('
        '"route_first.py", "route_first_nh.py", "route_first_marker.py", '
        '"route_first_marker_second.py", "routes.py", "reactions.py", '
        '"catalog.py", "workflow.py")}'
    )

    class RemoveMarkedNames(ast.NodeTransformer):
        count = 0

        def visit_DictComp(self, node):
            if dump(node) == new:
                self.count += 1
                return ast.parse(
                    '{name: _SOURCE_HASHES[name] for name in ('
                    '"route_first.py", "route_first_nh.py", "routes.py", '
                    '"reactions.py", "catalog.py", "workflow.py")}',
                    mode="eval",
                ).body
            return self.generic_visit(node)

    transformer = RemoveMarkedNames()
    tree = transformer.visit(ast.parse(source))
    assert transformer.count == 1
    return ast.unparse(tree)


def structural_digest_without_mapped_nh_manifest_source(source):
    """Remove only the reviewed helper name from one exact manifest tuple."""
    new = expression(
        "{name: _SOURCE_HASHES[name] for name in ("
        '"route_first.py", "route_first_nh.py", "routes.py", "reactions.py", '
        '"catalog.py", "workflow.py")}'
    )

    class RemoveReviewedManifestName(ast.NodeTransformer):
        count = 0

        def visit_DictComp(self, node):
            if dump(node) == new:
                self.count += 1
                return ast.parse(
                    "{name: _SOURCE_HASHES[name] for name in ("
                    '"route_first.py", "routes.py", "reactions.py", "catalog.py", "workflow.py")}',
                    mode="eval",
                ).body
            return self.generic_visit(node)

    transformer = RemoveReviewedManifestName()
    tree = transformer.visit(ast.parse(source))
    assert transformer.count == 1
    tree = DisabledTelemetry().visit(tree)
    return hashlib.sha256(
        json.dumps(canonical(tree), separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()


@pytest.mark.parametrize("path,expected", BASELINES.items())
def test_disabled_current_reviewed_source_structure_preserved(path, expected):
    source = (ROOT / path).read_bytes()
    if path in CURRENT_REVIEWED["files"]:
        expected = CURRENT_REVIEWED["files"][path]
        if path == "chembounce_synthesis/workflow.py":
            assert hashlib.sha256(source).hexdigest() == MARKED_SECOND_REVIEWED["workflow_source_sha256"]
            assert structural_digest(source.decode()) == MARKED_SECOND_REVIEWED["workflow_ast_sha256"]
            prior = without_marked_manifest_sources(source.decode())
            assert structural_digest(prior) == MAPPED_NH_REVIEWED["workflow_ast_sha256"]
            assert (
                structural_digest_without_mapped_nh_manifest_source(prior) == expected["ast_sha256"]
            )
            return
        assert hashlib.sha256(source).hexdigest() == expected["source_sha256"]
    assert structural_digest(source.decode()) == expected["ast_sha256"]


def test_only_known_empty_version_fields_are_ignored():
    for source in ("def f():\n    return 1", "async def f():\n    return 1", "class C:\n    pass"):
        node = ast.parse(source).body[0]
        node._fields = tuple(f for f in node._fields if f != "type_params")
        original = canonical(node)
        node._fields += ("type_params",)
        node.type_params = []
        assert canonical(node) == original
        node.type_params = [ast.Name(id="T", ctx=ast.Load())]
        assert canonical(node) != original
        node.type_params = []
        node._fields += ("future_semantic_field",)
        node.future_semantic_field = []
        assert canonical(node) != original


@pytest.mark.parametrize(
    "changed",
    [
        "def f():\n    return 2",
        "def f():\n    if telemetry is not None:\n        scientific_call()\n    return 1",
        'def f():\n    validate_telemetry(telemetry, "native", extra=True)\n    return 1',
        'def f():\n    if telemetry is not None and summary["deadline_reached"]:\n        telemetry_phase.interrupt("workflow_deadline_observed")\n        scientific_call()\n    return 1',
    ],
)
def test_real_semantic_changes_are_not_normalized(changed):
    assert structural_digest(changed) != structural_digest("def f():\n    return 1")
