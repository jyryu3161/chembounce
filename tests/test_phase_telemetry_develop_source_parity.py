"""Disabled telemetry parity against exact94b, independent of Git history.

Canonical AST retains every semantic field. Only empty type_params on the three
Python3.12 declaration nodes are omitted, matching their absence on3.10/3.11.
Nonempty type parameters and unknown fields remain compared. Digests were made
from the exact recorded baseline bytes, not from transformed post-port source.
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
        'if telemetry is not None and summary["deadline_reached"]:\n    phase.interrupt("workflow_deadline_observed")'
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
                    ast.Name(id="phase", ctx=ast.Store())
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


@pytest.mark.parametrize("path,expected", BASELINES.items())
def test_disabled_source_structure_preserved(path, expected):
    assert structural_digest((ROOT / path).read_text()) == expected["ast_sha256"]


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
        'def f():\n    if telemetry is not None and summary["deadline_reached"]:\n        phase.interrupt("workflow_deadline_observed")\n        scientific_call()\n    return 1',
    ],
)
def test_real_semantic_changes_are_not_normalized(changed):
    assert structural_digest(changed) != structural_digest("def f():\n    return 1")
