"""Synthetic contract tests. No checkpoint load, neural inference or GPU access."""

import copy
import json
import os
from pathlib import Path
import random
import sys
import time

import psutil
import pytest

from benchmarks.external.mol2mol import observation, selection, support
from benchmarks.research.runner import bounded_method


class Vocabulary:
    symbols = {0: "PAD", 1: "^", 2: "$", 3: "C", 4: "N"}

    def __getitem__(self, token):
        return {symbol: index for index, symbol in self.symbols.items()}[token]

    def __contains__(self, token):
        return token in self.symbols.values()

    def decode(self, sequence):
        return [self.symbols[int(x)] for x in sequence]


class Tokenizer:
    def tokenize(self, value):
        return ["^", *value, "$"]

    def untokenize(self, tokens):
        result = []
        for token in tokens:
            if token == "$":
                break
            if token != "^":
                result.append(token)
        return "".join(result)


def test_observer_preserves_object_result_rng_and_restores_class_method():
    vocabulary = Vocabulary()
    result = ["same object"]
    vocabulary.decode = lambda seq: result
    original = vocabulary.decode
    state = random.getstate()
    with observation.observe_decode(vocabulary) as captured:
        assert vocabulary.decode([1, 3, 2]) is result
    assert captured == [[1, 3, 2]]
    assert vocabulary.decode is original
    assert random.getstate() == state
    del vocabulary.decode
    with observation.observe_decode(vocabulary):
        vocabulary.decode([1, 3, 2])
    assert "decode" not in vars(vocabulary)


def test_observer_preserves_exception_identity_and_restoration():
    vocabulary = Vocabulary()
    error = RuntimeError("official decode error")

    def failing(sequence):
        raise error

    vocabulary.decode = failing
    with pytest.raises(RuntimeError) as raised:
        with observation.observe_decode(vocabulary) as captured:
            vocabulary.decode([1])
    assert raised.value is error and vocabulary.decode is failing and captured == []


def bound(captured=None, outputs=None):
    return observation.bind_outputs(
        captured or [[1, 3, 2], [1, 3, 2], [1, 3, 2, 0], [1, 4, 4, 4]],
        [[1, 3, 2], [1, 3, 2]],
        ["C", "C"],
        outputs or ["C", "NNN"],
        [1.25, 4.5],
        eos=2,
        bos=1,
        max_length=4,
        tokenizer=Tokenizer(),
        vocabulary=Vocabulary(),
    )


def test_decode_binding_keeps_each_nll_and_eos_padding_or_truncation():
    rows = bound()
    assert rows[0] == {
        "batch_row": 0,
        "smiles": "C",
        "native_nll": 1.25,
        "token_ids": [1, 3, 2, 0],
        "terminated": True,
        "truncated": False,
        "eos_position": 2,
        "padded_length": 4,
    }
    assert rows[1]["native_nll"] == 4.5 and rows[1]["truncated"]
    assert rows[1]["eos_position"] is None


@pytest.mark.parametrize(
    "captured,outputs,reason",
    [
        ([[1, 4, 2], [1, 3, 2], [1, 3, 2, 0], [1, 4, 4, 4]], None, "input decode order"),
        ([[1, 3, 2]], None, "row count"),
        (None, ["N", "NNN"], "binding mismatch"),
        ([[1, 3, 2], [1, 3, 2], [1, 3, 2, 4], [1, 4, 4, 4]], None, "non-padding"),
        ([[1, 3, 2], [1, 3, 2], [1, 3, 2, 0], [1, 4, 4]], None, "length limit"),
    ],
)
def test_decode_contract_change_fails_closed(captured, outputs, reason):
    with pytest.raises(ValueError, match=reason):
        bound(captured, outputs)


def test_unsupported_input_tokens_are_rejected_before_official_encode_can_drop_them():
    assert observation.validate_tokens("CN", Tokenizer(), Vocabulary(), 4) == ["^", "C", "N", "$"]
    with pytest.raises(ValueError, match="Unsupported input tokens"):
        observation.validate_tokens("CO", Tokenizer(), Vocabulary(), 4)
    with pytest.raises(ValueError, match="positional capacity"):
        observation.validate_tokens("CCN", Tokenizer(), Vocabulary(), 4)


def row(smiles, nll, *, truncated=False):
    return {"smiles": smiles, "native_nll": nll, "terminated": not truncated, "truncated": truncated}


def test_common_filter_preserves_all_raw_outcomes_and_ranks_families_by_native_nll():
    raw = [
        row("broken", 0),
        row("c1ccncc1.CC", 0),
        row("CCCC", 0),
        row("c1ccccc1", 0.1),
        row("Cc1ccccc1", 0.2),
        row("c1ccncc1", 3),
        row("c1ccncc1", 2),
        row("C[C@H]1CCCN1", 1),
        row("C[C@@H]1CCCN1", 4),
        row("c1ccsc1", 0.01, truncated=True),
        row("c1ccoc1", "nan"),
        row("C1CCOC1", 2),
        row("C1CCSC1", 2),
    ]
    original = copy.deepcopy(raw)
    candidates, records = selection.select("c1ccccc1", raw)
    assert raw == original and len(records) == len(raw)
    assert [x["native_nll"] for x in candidates] == [1, 2, 2]
    assert [x["smiles"] for x in candidates[1:]] == ["C1CCOC1", "C1CCSC1"]
    assert records[0]["selection_status"] == "invalid_disconnected_or_identity_failure"
    assert records[1]["selection_status"] == "invalid_disconnected_or_identity_failure"
    assert records[2]["selection_status"] == "acyclic"
    assert records[3]["selection_status"] == "input_family"
    assert records[4]["selection_status"] == "unchanged_normalized_murcko"
    assert records[5]["selection_status"] == "duplicate_family"
    assert records[6]["selection_status"] == "output_cap"
    assert records[8]["selection_status"] == "duplicate_family"
    assert records[9]["selection_status"] == "unterminated_at_model_limit"
    assert records[10]["selection_status"] == "nonfinite_native_nll"
    assert all(c["route_evidence"]["level"] == "none" and "route" not in c for c in candidates)


def test_common_identity_normalizes_input_state_and_fails_closed_for_identity_error(monkeypatch):
    candidates, rows = selection.select("c1ccncc1", [row("c1cc[nH+]cc1", 1)])
    assert not candidates and rows[0]["selection_status"] == "input_family"
    import chembounce_synthesis.identity as identity

    original = identity.chemical_identity

    def fails_one(smiles):
        if smiles == "C1CCOC1":
            raise ValueError("tautomer limit")
        return original(smiles)

    monkeypatch.setattr(identity, "chemical_identity", fails_one)
    candidates, rows = selection.select("c1ccccc1", [row("C1CCOC1", 0), row("C1CCSC1", 1)])
    assert len(candidates) == 1 and candidates[0]["smiles"] == "C1CCSC1"
    assert "tautomer limit" in rows[0]["detail"]
    with pytest.raises(TimeoutError):
        selection.select("c1ccccc1", [row("C1CCSC1", 1)], deadline=time.monotonic() - 1)


def request():
    return {
        "schema_version": 1,
        "kind": "external",
        "method_id": "mol2mol_pretrained",
        "query": {"id": "q", "smiles": "c1ccccc1", "target_id": "public_id"},
        "seed": 42,
        "settings": dict(support.SETTINGS),
        "assets": {
            "mol2mol_bundle": "/unused",
            **support.asset_paths(
                {"model": "/unused/model", "adapter_root": "/unused/adapter", "common_root": "/unused/common"}
            ),
        },
        "budget": {"wall_seconds": 360, "max_memory_mb": 24000, "max_candidates": 3, "max_evaluations": 256},
        "output_dir": "/unused",
        "deadline_monotonic": time.monotonic() + 360,
        "finalization_margin_seconds": 2,
    }


@pytest.mark.parametrize(
    "place,key,value",
    [
        ("settings", "check_nll", True),
        ("settings", "target_smiles_path", "/do-not-read"),
        ("query", "labels", ["active"]),
        ("query", "reference_sdf", "/do-not-read"),
        ("settings", "temperature", 0.5),
        ("settings", "input_randomization", True),
        ("assets", "test_data", "/do-not-read"),
        ("budget", "max_candidates", 4),
    ],
)
def test_whitelist_rejects_extra_information_and_posthoc_setting_changes(place, key, value):
    valid = request()
    support.validate_request(valid)
    valid[place][key] = value
    with pytest.raises(ValueError):
        support.validate_request(valid)


def test_hash_tamper_rejected_and_import_origin_requires_pinned_tracked_file(tmp_path, monkeypatch):
    source = tmp_path / "official.py"
    source.write_text("# original\n")
    entries = {source.name: support.digest(source)}
    support.check_hashes(tmp_path, entries)
    source.write_text("# changed\n")
    with pytest.raises(ValueError, match="Pinned file changed"):
        support.check_hashes(tmp_path, entries)
    import types

    fake = types.ModuleType("reinvent.fixture")
    fake.__file__ = str(source)
    monkeypatch.setitem(sys.modules, fake.__name__, fake)
    bundle = {"upstream": str(tmp_path), "source_sha256": entries}
    with pytest.raises(ValueError, match="not a pinned tracked"):
        support.upstream_origins(bundle)
    bundle["source_sha256"][source.name] = support.digest(source)
    assert support.upstream_origins(bundle)[fake.__name__]["relative_path"] == source.name
    fake.__file__ = str(Path(__file__))
    with pytest.raises(ValueError, match="outside pinned checkout"):
        support.upstream_origins(bundle)


@pytest.mark.skipif(os.name != "posix", reason="POSIX owned group")
@pytest.mark.parametrize("mode", ["timeout", "nonzero", "success"])
def test_outer_runner_controls_adapter_sampler_grandchild_and_rejects_partial(tmp_path, mode):
    # The fixture exercises our real nested spawning helper under the real outer
    # runner. It does not load a checkpoint or substitute a scientific baseline.
    nested = tmp_path / "nested.py"
    nested.write_text("""import json,os,signal,subprocess,sys,time
from pathlib import Path
out=Path(sys.argv[1]); mode=sys.argv[2]
code="import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(30)"
p=subprocess.Popen([sys.executable,'-c',code])
(out/'group.json').write_text(json.dumps({'sampler':os.getpid(),'group':os.getpgrp(),'grandchild':p.pid}))
(out/'response.json').write_text(json.dumps({'candidates':[{'smiles':'C1CCOC1'}],'evaluations':1,'evaluation_unit':'fixture'}))
if mode=='timeout': time.sleep(30)
if mode=='nonzero': sys.exit(9)
""")
    outer = tmp_path / "outer.py"
    root = Path(__file__).resolve().parents[1]
    outer.write_text(f"""import json,os,sys
from pathlib import Path
sys.path.insert(0,{str(root)!r})
from benchmarks.external.mol2mol.adapter import run_child
r=json.loads(Path(sys.argv[1]).read_text()); out=Path(r['output_dir'])
(out/'adapter-group.json').write_text(json.dumps({{'pid':os.getpid(),'group':os.getpgrp()}}))
run_child([sys.executable,{str(nested)!r},str(out),{mode!r}],deadline=r['deadline_monotonic'],log_dir=out)
""")
    output = tmp_path / "attempt"
    req = {
        "output_dir": str(output),
        "budget": {"wall_seconds": 1.5, "max_candidates": 3, "max_evaluations": 256, "max_memory_mb": 24000},
    }
    method = {"kind": "external", "command": [sys.executable, str(outer), "{request}", "{response}"]}
    result = bounded_method(method, req, output)
    groups = json.loads((output / "group.json").read_text())
    outer_group = json.loads((output / "adapter-group.json").read_text())
    assert groups["group"] == outer_group["group"] == outer_group["pid"]
    if mode == "success":
        assert result["status"] == "success" and len(result["response"]["candidates"]) == 1
    else:
        assert result["status"] in {"timeout", "worker_error"} and result["response"] == {}
    deadline = time.monotonic() + 3
    while time.monotonic() < deadline:
        try:
            if psutil.Process(groups["grandchild"]).status() == psutil.STATUS_ZOMBIE:
                break
        except psutil.NoSuchProcess:
            break
        time.sleep(0.02)
    else:
        pytest.fail("Outer hard cleanup left sampler descendant alive")
    assert result["cleanup"]["leader_reaped"]


@pytest.fixture
def official_shaped_fixture(monkeypatch):
    """Exercise orchestration, not a neural substitute or checkpoint equivalence."""
    from types import SimpleNamespace as NS
    from benchmarks.external.mol2mol import worker

    clock = {"now": 0.0, "calls": 0, "seeds": [], "duration": 2.0}
    monkeypatch.setattr(worker.time, "monotonic", lambda: clock["now"])
    monkeypatch.setattr(worker, "upstream_origins", lambda bundle: {})

    class Tensor:
        def __init__(self, rows):
            self.rows = rows

        def tolist(self):
            return self.rows

        def to(self, device):
            return self

    class Model:
        def __init__(self):
            self.vocabulary = Vocabulary()
            self.tokenizer = Tokenizer()
            self.max_sequence_length = 4
            self.network = NS(src_embed=[None, NS(pe=NS(shape=(1, 5000, 8)))])

        @classmethod
        def load_from_file(cls, path, mode, device):
            assert mode == "inference"
            return cls()

        def set_temperature(self, value):
            assert value == 1.0

        def sample(self, src, mask, strategy):
            assert strategy == "multinomial"
            clock["calls"] += 1
            clock["now"] += clock["duration"]
            inputs = [self.tokenizer.untokenize(self.vocabulary.decode(row)) for row in src.rows]
            outputs = [self.tokenizer.untokenize(self.vocabulary.decode([1, 4, 2])) for _ in src.rows]
            return inputs, outputs, [float(i) for i in range(len(src.rows))]

    class Dataset:
        def __init__(self, strings, vocabulary, tokenizer):
            self.rows = [[vocabulary[t] for t in tokenizer.tokenize(s)] for s in strings]

        def __len__(self):
            return len(self.rows)

        def __getitem__(self, index):
            return self.rows[index]

        @staticmethod
        def collate_fn(rows):
            return Tensor(rows), Tensor([])

    def seed(value):
        clock["seeds"].append(value)

    monkeypatch.setattr(worker.random, "seed", seed)
    cuda = NS(
        is_available=lambda: True,
        device_count=lambda: 1,
        set_device=lambda d: None,
        manual_seed_all=seed,
        reset_peak_memory_stats=lambda d: None,
        get_device_name=lambda d: "fixture-only",
        max_memory_allocated=lambda d: 123,
        max_memory_reserved=lambda d: 456,
    )
    torch = NS(
        set_num_threads=lambda n: None,
        set_num_interop_threads=lambda n: None,
        manual_seed=seed,
        cuda=cuda,
        device=lambda value: value,
        __version__="fixture",
        version=NS(cuda="fixture"),
    )
    numpy = NS(random=NS(seed=seed))
    monkeypatch.setattr(worker, "official_api", lambda bundle: (torch, numpy, Model, Dataset, {}))
    return worker, clock


@pytest.mark.parametrize(
    "deadline,expected_batches,stop", [(1000, 8, "raw_cap"), (20, 2, "early_finish_reserve")]
)
def test_batch_cap_seed_once_and_fixed_anytime_reserve(
    official_shaped_fixture, tmp_path, deadline, expected_batches, stop
):
    worker, clock = official_shaped_fixture
    path = tmp_path / "raw.json"
    # Input is longer than generated output cap; encoder capacity is the true input limit.
    worker.sample(
        {"deadline_monotonic": deadline, "settings": dict(support.SETTINGS), "seed": 42, "smiles": "CCCCC"},
        {"model": "unused-fixture"},
        path,
    )
    result = support.read_json(path)
    assert clock["calls"] == len(result["batches"]) == expected_batches
    assert clock["seeds"] == [42, 42, 42, 42]  # Python, numpy, torch, CUDA; never per batch
    assert len(result["rows"]) == expected_batches * 32
    assert result["stop_reason"] == stop
    assert result["input_token_capacity"] == 5000 and result["max_sequence_length"] == 4
    assert [r["native_nll"] for r in result["rows"][:32]] == list(map(float, range(32)))


def test_unsupported_input_fails_before_any_sample_and_crossing_deadline_fails(
    official_shaped_fixture, tmp_path
):
    worker, clock = official_shaped_fixture
    job = {"deadline_monotonic": 100, "settings": dict(support.SETTINGS), "seed": 42, "smiles": "CO"}
    with pytest.raises(ValueError, match="Unsupported input"):
        worker.sample(job, {"model": "unused-fixture"}, tmp_path / "bad.json")
    assert clock["calls"] == 0
    job["smiles"] = "C"
    clock["duration"] = 101
    with pytest.raises(TimeoutError, match="no accepted partial"):
        worker.sample(job, {"model": "unused-fixture"}, tmp_path / "late.json")
    assert len(support.read_json(tmp_path / "late.json")["rows"]) == 32


def test_checkpoint_tamper_rejected_before_source_import_or_pickle(tmp_path, monkeypatch):
    model = tmp_path / "not-official.prior"
    model.write_bytes(b"not a pickle")
    monkeypatch.setattr(support, "check_hashes", lambda *args: pytest.fail("Reached source/import stage"))
    with pytest.raises(ValueError, match="refusing pickle load"):
        support.verify_bundle({"upstream_commit": support.UPSTREAM_COMMIT, "model": str(model)})


def test_nonfinite_deadline_and_unowned_direct_adapter_invocation_fail(tmp_path, monkeypatch):
    from benchmarks.external.mol2mol import adapter

    invalid = request()
    invalid["deadline_monotonic"] = float("nan")
    with pytest.raises(ValueError, match="Finite outer runner"):
        support.validate_request(invalid)
    monkeypatch.setattr(adapter.os, "getpgrp", lambda: 1)
    monkeypatch.setattr(adapter.os, "getpid", lambda: 2)
    with pytest.raises(RuntimeError, match="owned process group"):
        adapter.run(request(), tmp_path / "response.json")


def test_runner_checks_adapter_source_assets_before_external_execution(tmp_path):
    from benchmarks.research.runner import load_protocol

    fake = {
        "model": str(tmp_path / "model.prior"),
        "adapter_root": str(tmp_path / "adapter"),
        "common_root": str(tmp_path / "common"),
    }
    paths = support.asset_paths(fake)
    paths["mol2mol_bundle"] = str(tmp_path / "bundle.json")
    assets = {}
    for name, filename in paths.items():
        path = Path(filename)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture, never executed\n")
        assets[name] = {"artifact": str(path), "sha256": support.digest(path)}
    method = {
        "id": "mol2mol_pretrained",
        "kind": "external",
        "assets": assets,
        "settings": dict(support.SETTINGS),
        "command": [sys.executable, paths["adapter_source_adapter"], "{request}", "{response}"],
    }
    protocol = {
        "schema_version": 1,
        "seeds": [42, 43, 44],
        "budget": request()["budget"],
        "methods": [method],
    }
    path = tmp_path / "protocol.json"
    support.write_json(path, protocol)
    load_protocol(path, {})
    Path(paths["adapter_source_adapter"]).write_text("changed adapter could ignore self-checks\n")
    with pytest.raises(ValueError, match="[Hh]ash|SHA|sha|[Dd]igest"):
        load_protocol(path, {})


@pytest.mark.parametrize(
    "name",
    [
        "mol2mol-job.json",
        "mol2mol-raw.json",
        "mol2mol-selection.json",
        "sampler.stdout.log",
        "sampler.stderr.log",
        "response.json",
        "mol2mol-response-provenance.json",
    ],
)
def test_fresh_attempt_guard_preserves_every_owned_evidence_file(tmp_path, monkeypatch, name):
    from benchmarks.external.mol2mol import adapter

    existing = tmp_path / name
    existing.write_bytes(b"immutable original evidence\n")
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    r = request()
    r["output_dir"] = str(tmp_path)
    monkeypatch.setattr(adapter.os, "getpgrp", lambda: 123)
    monkeypatch.setattr(adapter.os, "getpid", lambda: 123)
    monkeypatch.setattr(
        adapter, "read_json", lambda p: pytest.fail("Read model bundle before freshness guard")
    )
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        adapter.run(r, tmp_path / "response.json")
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before


def test_fresh_attempt_guard_rejects_dangling_link_and_stale_atomic_output(tmp_path):
    from benchmarks.external.mol2mol.adapter import fresh_output

    link = tmp_path / "mol2mol-raw.json"
    link.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError):
        fresh_output(tmp_path, tmp_path / "response.json")
    link.unlink()
    stale = tmp_path / "mol2mol-raw.json.123.tmp"
    stale.write_text("incomplete old evidence")
    with pytest.raises(FileExistsError):
        fresh_output(tmp_path, tmp_path / "response.json")
    stale.unlink()
    assert fresh_output(tmp_path, tmp_path / "response.json") == tmp_path


@pytest.mark.parametrize("seed", [42.0, "42", True, 45])
def test_seed_requires_exact_declared_integer(seed):
    r = request()
    r["seed"] = seed
    with pytest.raises(ValueError, match="Expected proposal seed"):
        support.validate_request(r)


@pytest.mark.parametrize(
    "available,count,reason", [(False, 0, "no silent CPU fallback"), (True, 2, "exactly one authorized GPU")]
)
def test_declared_gpu_runtime_fails_without_model_sampling(
    official_shaped_fixture, tmp_path, available, count, reason, monkeypatch
):
    worker, clock = official_shaped_fixture
    api = worker.official_api({})
    api[0].cuda.is_available = lambda: available
    api[0].cuda.device_count = lambda: count
    monkeypatch.setattr(worker, "official_api", lambda bundle: api)
    with pytest.raises(RuntimeError, match=reason):
        worker.sample(
            {"deadline_monotonic": 100, "settings": dict(support.SETTINGS), "seed": 42, "smiles": "C"},
            {"model": "unused"},
            tmp_path / "raw.json",
        )
    assert clock["calls"] == 0 and not (tmp_path / "raw.json").exists()
