"""Synthetic administrative tests; no real assets/scientific imports or preparation."""

import copy
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace, ModuleType

import pytest

from benchmarks.research import polar_capacity_prepare as p
from benchmarks.research import polar_capacity_development as d
from benchmarks.research import polar_capacity_observer as observer
from benchmarks.research import polar_capacity_supervisor as s


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(value, (dict, list)):
        s.write(path, value)
    else:
        path.write_bytes(value if isinstance(value, bytes) else value.encode())
    return str(path.resolve()), d.digest(path)


def authority(tmp_path, monkeypatch):
    """Real pin/archive/git validation, controlled authority constants only.

    Temporary source trees and runtime records stand for the deployment authority.
    validate_plan/load_plan/hash/archive checks are never replaced.
    """
    pins = {}

    def add(path, value):
        name, sha = put(path, value)
        pins[name] = sha
        return name

    realroot = Path(p.__file__).resolve().parents[2]
    roles = {}
    for role in ("baseline", "capacity", "checker"):
        root = tmp_path / role
        package = root / "chembounce_synthesis"
        for name in (
            "common.py",
            "identity.py",
            "conformers.py",
            "selective_scoring.py",
            "similarity_constraints.py",
        ):
            add(package / name, "# controlled source, never imported\n")
        for name in ("qualified-development-v1.json", d.PROFILE + ".json"):
            add(package / "data" / name, (realroot / "chembounce_synthesis/data" / name).read_bytes())
        if role == "checker":
            for name in (
                "polar_capacity_prepare.py",
                "polar_capacity_development.py",
                "polar_capacity_observer.py",
                "polar_capacity_supervisor.py",
                "qualified_collector.py",
                "qualified_adapter.py",
            ):
                add(root / "benchmarks/research" / name, "# controlled source inventory\n")
        subprocess.run(["git", "init", "-q", str(root)], check=True)
        subprocess.run(["git", "-C", str(root), "add", "."], check=True)
        subprocess.run(
            [
                "git",
                "-C",
                str(root),
                "-c",
                "user.name=Synthetic",
                "-c",
                "user.email=synthetic@example.invalid",
                "commit",
                "-qm",
                role,
            ],
            check=True,
        )
        roles[role] = dict(
            path=str(root),
            commit=subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip(),
            tree=subprocess.check_output(
                ["git", "-C", str(root), "rev-parse", "HEAD^{tree}"], text=True
            ).strip(),
            files={name: sha for name, sha in pins.items() if Path(name).is_relative_to(root)},
        )
    monkeypatch.setattr(d, "ARMS", {r: roles[r]["commit"] for r in ("baseline", "capacity")})
    monkeypatch.setattr(
        p, "__file__", str(Path(roles["checker"]["path"]) / "benchmarks/research/polar_capacity_prepare.py")
    )
    monkeypatch.setattr(
        d,
        "__file__",
        str(Path(roles["checker"]["path"]) / "benchmarks/research/polar_capacity_development.py"),
    )
    catalog = add(tmp_path / "assets/catalog.tsv", "not parsed")
    sidecar = add(tmp_path / "assets/catalog.tsv.manifest.json", {})
    rules = add(tmp_path / "assets/rules.json", {})
    feature = add(tmp_path / "assets/BaseFeatures.fdef", "not perceived")
    model_file = add(tmp_path / "model/model.safetensors", b"tiny synthetic bytes, never loaded")
    model = str(Path(model_file).parent)
    slots = []
    for name in ("attempt-0000", "attempt-0009"):
        root = tmp_path / "historical" / name
        pool = ["literal-design-" + str(i) for i in range(12)]
        query = dict(id=name, smiles="unparsed-literal")
        paths = [
            add(root / "native/design-pool.json", dict(designs=pool)),
            add(root / "common/request.json", dict(query=query)),
            add(
                root / "common/synthesis/proposals.json",
                dict(
                    proposals=[
                        dict(design=pool[i], sampling={"seed": seed})
                        for i, seed in zip([1, 1, 2, 2, 3, 3], [42, 43] * 3)
                    ]
                ),
            ),
        ]
        slots.append(
            dict(
                historical_slot=name,
                bindings={x: pins[x] for x in paths},
                model_metadata=dict(path=model, model_files_sha256={"model.safetensors": pins[model_file]}),
            )
        )
    locked = dict(
        recommendation_profile=d.PROFILE,
        mode="hop",
        minimum_scores=dict(hit_tanimoto=0.35, hit_shape=0.7, hit_pharmacophore=0.5),
        similarity_constraints=dict(
            version="selective-similarity-constraints-v2", polar_recall=0.75, landmarks=[]
        ),
        ranking_policy="original-hit-v1",
        route_search_policy="original-hit-v1",
        route_search_cache_size=4096,
        max_design_attempts=12,
        max_designs=3,
        max_scored_products=12,
        max_results=3,
        search=dict(
            beam_width=6,
            building_block_top_k=3,
            max_combinations=256,
            max_steps=12,
            reaction_product_limit=64,
        ),
        expansion=dict(building_block_top_k=6, max_routes=12, min_tanimoto=0.0, rounds=1),
    )
    design_path = add(
        tmp_path / "design.json",
        dict(
            locked_configuration=locked,
            development_bindings=dict(slots=slots),
            assets=dict(
                catalog=dict(path=catalog, recorded_sha256=pins[catalog]),
                rules=dict(path=rules, recorded_sha256=pins[rules]),
                catalog_sidecar=dict(path=sidecar, sha256=pins[sidecar]),
            ),
        ),
    )
    gate = add(tmp_path / "gate.json", dict(scope="synthetic"))
    monkeypatch.setattr(d, "DESIGN_SHA", pins[design_path])
    monkeypatch.setattr(d, "GATE_SHA", pins[gate])
    runtime = add(tmp_path / "runtime/metadata.txt", "synthetic runtime inventory")
    python = str(Path(sys.executable).resolve())
    pins[python] = d.digest(python)
    monkeypatch.setattr(d, "required_runtime_files", lambda: {runtime, python})
    contexts = {
        r: dict(cwd=x["path"], PYTHONPATH=x["path"], executable=python, python="synthetic", packages={})
        for r, x in roles.items()
    }
    review = add(tmp_path / "review.json", dict(status="synthetic"))
    clearance = add(tmp_path / "clearance.json", dict(status="synthetic"))
    spec = dict(
        schema=p.SCHEMA,
        design=dict(path=design_path, sha256=pins[design_path]),
        gate=dict(path=gate, sha256=pins[gate]),
        sources={r: roles[r] for r in ("baseline", "capacity")},
        checker=roles["checker"],
        contexts=contexts,
        pins=pins,
        python=str(Path(sys.executable).absolute()),
        feature_definition=dict(path=feature, sha256=pins[feature]),
        generation_source=dict(path=runtime, sha256=pins[runtime]),
        output=str(tmp_path / "preparation"),
        run_output=str(tmp_path / "study"),
        workloads=dict(
            identities=[dict(pid=999999, start_ticks=1)], command_markers=["synthetic-do-not-launch"]
        ),
        reviews=[dict(path=review, sha256=pins[review])],
        clearance=dict(path=clearance, sha256=pins[clearance]),
    )
    path = tmp_path / "spec.json"
    s.write(path, spec)
    auth = tmp_path / "auth.json"
    s.write(
        auth,
        dict(
            schema=p.AUTH,
            spec_sha256=d.digest(path),
            source_review=spec["reviews"][0],
            literal_review=spec["reviews"][0],
            clearance=spec["clearance"],
            allow_asset_hashing=True,
            allow_import_discovery=True,
        ),
    )
    args = SimpleNamespace(
        spec=str(path),
        spec_sha256=d.digest(path),
        authorization=str(auth),
        authorization_sha256=d.digest(auth),
        output=spec["output"],
    )
    return spec, args


def test_real_spec_authority_and_archive(tmp_path, monkeypatch):
    spec, args = authority(tmp_path, monkeypatch)
    assert p.authorize(args) == spec
    queries = p.queries_from_design(spec)
    assert list(queries) == ["attempt-0000", "attempt-0009"]
    for q in queries.values():
        d.validate_config(q["config"], q["proposals"])
        assert len(d.validate_archive(q)) == 12
        assert q["config"]["minimum_scores"] == dict(hit_tanimoto=0.35, hit_shape=0.7, hit_pharmacophore=0.5)


@pytest.mark.parametrize(
    "mutation",
    ["extra", "profile", "source", "runtime", "output", "pin", "dirty", "auth", "asset", "archive"],
)
def test_spec_and_inputs_fail_closed(tmp_path, monkeypatch, mutation):
    spec, args = authority(tmp_path, monkeypatch)
    if mutation == "extra":
        spec["execution_authorized"] = True
    if mutation == "profile":
        spec["schema"] = "polar-capacity-development-launch-v1"
    if mutation == "source":
        spec["sources"]["baseline"]["commit"] = "bad"
    if mutation == "runtime":
        spec["contexts"]["baseline"]["cwd"] = "/wrong"
    if mutation == "output":
        spec["run_output"] = spec["output"]
    if mutation == "pin":
        spec["pins"].pop(spec["feature_definition"]["path"])
    if mutation == "dirty":
        (Path(spec["sources"]["baseline"]["path"]) / "unknown.py").write_text("pass")
    if mutation == "asset":
        Path(spec["feature_definition"]["path"]).write_text("changed")
    if mutation == "archive":
        design = d.read(spec["design"]["path"])
        path = next(
            x for x in design["development_bindings"]["slots"][0]["bindings"] if x.endswith("proposals.json")
        )
        Path(path).write_text("{}")
    if mutation == "auth":
        auth = d.read(args.authorization)
        auth["schema"] = "polar-capacity-development-launch-v1"
        s.write(args.authorization, auth)
        args.authorization_sha256 = d.digest(args.authorization)
        with pytest.raises(ValueError):
            p.authorize(args)
    else:
        with pytest.raises((ValueError, d.IntegrityError, FileExistsError)):
            p.validate_spec(spec)


def synthetic_discoveries(spec, queries):
    config = {key: value["config"] for key, value in queries.items()}
    source = str(Path(spec["checker"]["path"]) / "benchmarks/research/polar_capacity_observer.py")
    return {
        r: dict(
            configs=copy.deepcopy(config),
            observer=[dict(source_path=source, source_sha256=spec["pins"][source], implementation=[])],
        )
        for r in d.ARMS
    }


def test_actual_existing_prepare_with_full_synthetic_authority(tmp_path, monkeypatch):
    spec, _ = authority(tmp_path, monkeypatch)
    queries = p.queries_from_design(spec)
    plan = p.assemble_plan(spec, queries, synthetic_discoveries(spec, queries))
    assert d.validate_plan(plan) == plan  # no validator/hash/archive stubs
    path = tmp_path / "candidate.json"
    s.write(path, plan)
    for key, value in d.environment(spec["checker"]["path"]).items():
        monkeypatch.setenv(key, value)
    output = tmp_path / "existing-prepare"
    assert (
        d.main(["prepare", "--plan", str(path), "--plan-sha256", d.digest(path), "--output", str(output)])
        == 0
    )
    terminal = d.read(output / "terminal.json")
    assert terminal["status"] == "prepared_not_authorized"
    ledger = d.read(output / "ledger.json")
    assert ledger["rows"] == [
        dict(index=i, status="not_started", request=d.request_for(plan, d.digest(path), i)) for i in range(4)
    ]
    assert not Path(spec["run_output"]).exists()
    assert not list(tmp_path.rglob("*launch*authorization*"))
    # Rebound candidate mutation still reaches actual immutable-plan validation.
    plan["queries"]["attempt-0000"]["config"]["minimum_scores"]["hit_shape"] = 0.69
    s.write(path, plan)
    assert (
        d.main(
            [
                "prepare",
                "--plan",
                str(path),
                "--plan-sha256",
                d.digest(path),
                "--output",
                str(tmp_path / "rejected"),
            ]
        )
        == 1
    )
    assert not (tmp_path / "rejected/ledger.json").exists()


def test_arm_config_difference_rejected(tmp_path, monkeypatch):
    spec, _ = authority(tmp_path, monkeypatch)
    queries = p.queries_from_design(spec)
    discovery = synthetic_discoveries(spec, queries)
    discovery["capacity"]["configs"]["attempt-0000"]["max_results"] = 1
    with pytest.raises(ValueError, match="different scientific"):
        p.assemble_plan(spec, queries, discovery)


@pytest.mark.parametrize("which", ["python", "caught_python", "observed", "asset", "network"])
def test_guard_denies_before_sentinel_without_replacing_callables(tmp_path, which):
    # Fresh stdlib subprocess: no real scientific dependency import.
    output = tmp_path / "probe"
    output.mkdir()
    sentinel = tmp_path / "sentinel"
    asset = tmp_path / "forbidden.bin"
    asset.write_bytes(b"synthetic")
    code = """import sys,types,pathlib
from benchmarks.research.polar_capacity_prepare import ProbeGuard
from benchmarks.research.polar_capacity_observer import Observer
out, sentinel, asset, which = map(str,sys.argv[1:])
g=ProbeGuard(out,{},[asset]);g.install()
m=types.ModuleType('chembounce_synthesis.fake');m.__file__=__file__ if '__file__' in globals() else sys.argv[0]
exec("def score_candidate():\\n open("+repr(sentinel)+",'w').write('forbidden')",m.__dict__)
original=m.score_candidate
if which in ('python','caught_python'):
 try: m.score_candidate()
 except BaseException: pathlib.Path(sentinel).write_text('caught')
elif which=='asset': open(asset,'rb').read()
elif which=='network':
 import socket
 socket.getaddrinfo('invalid.example',80)
else:
 class Denying(Observer):
  def emit(self,row):
   if row.get('event')=='started':g.fatal('observed_scientific_call_denied')
   super().emit(row)
 m.__name__='synthetic_owner';sys.modules[m.__name__]=m
 m.__file__=str(pathlib.Path(Observer.__module__.replace('.','/')).with_suffix('.py').resolve())
 o=Denying(pathlib.Path(out)/'events',module_prefix='synthetic_owner');o.install(m,'score_candidate',kind='science');m.score_candidate()
assert m.score_candidate is original
"""
    result = subprocess.run(
        [sys.executable, "-c", code, str(output), str(sentinel), str(asset), which],
        cwd=Path(__file__).resolve().parents[1],
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        capture_output=True,
        text=True,
        timeout=10,
    )
    assert result.returncode == 79, result.stderr
    assert not sentinel.exists()
    events = [json.loads(x) for x in (output / "guard-events.jsonl").read_text().splitlines()]
    assert len(events) == 1


@pytest.mark.parametrize("mutation", [None, "nonidempotent", "custom_policy", "wrong_inventory"])
def test_probe_uses_existing_observer_and_restores_cache(tmp_path, monkeypatch, mutation):
    """Fake module imports exercise real observer/guard without scientific code."""
    spec, _ = authority(tmp_path, monkeypatch)
    queries = p.queries_from_design(spec)
    monkeypatch.setattr(
        p,
        "discover_runtime",
        lambda spec: {spec["generation_source"]["path"]: spec["generation_source"]["sha256"]},
    )
    monkeypatch.setattr(p, "imported_file_pins", lambda: {})
    from functools import lru_cache

    package = ModuleType("chembounce_synthesis")
    workflow = ModuleType("chembounce_synthesis.workflow")
    conformers = ModuleType("chembounce_synthesis.conformers")
    conformers.__file__ = str(Path(observer.__file__).resolve())

    @lru_cache(maxsize=4)
    def ensemble(value):
        raise AssertionError("No scientific calls allowed")

    conformers.ensemble = ensemble
    workflow._config = lambda x: x
    workflow._resolve_paths = lambda x, base: x
    workflow.effective_policy = lambda x: {"calibration": {"status": "engineering_policy_matched"}}
    package.workflow, package.conformers = workflow, conformers
    rdkit = ModuleType("rdkit")
    rdkit.RDConfig = SimpleNamespace(RDDataDir=str(Path(spec["feature_definition"]["path"]).parent))
    for name, obj in [
        ("chembounce_synthesis", package),
        (workflow.__name__, workflow),
        (conformers.__name__, conformers),
        ("rdkit", rdkit),
    ]:
        monkeypatch.setitem(sys.modules, name, obj)
    monkeypatch.setattr(
        d, "context_check", lambda *a: None
    )  # isolated probe interface; authority covered separately
    monkeypatch.setattr(d, "runtime", lambda: spec["contexts"]["baseline"])

    def install(obs):
        obs.install(conformers, "ensemble", kind="ensemble")
        obs.seal()

    monkeypatch.setattr(observer, "install_scientific", install)
    out = tmp_path / "probe"
    out.mkdir()
    if mutation:
        expected = None
        if mutation == "nonidempotent":

            def changing(config):
                config["extra"] = config.get("extra", 0) + 1
                return config

            workflow._config = changing
        elif mutation == "custom_policy":
            workflow.effective_policy = lambda x: {"calibration": {"status": "custom_unvalidated"}}
        else:
            expected = [{"wrong": "inventory"}]
        with pytest.raises(ValueError):
            p.probe(spec, "baseline", queries, out, expected)
        assert conformers.ensemble is ensemble and sys.getprofile() is None
        return
    result = p.probe(spec, "baseline", queries, out)
    assert result["scientific_calls"] == 0 and result["cache_before"] == result["cache_after"]
    assert conformers.ensemble is ensemble
    second = tmp_path / "verify"
    second.mkdir()
    verified = p.probe(spec, "baseline", queries, second, result["observer"])
    assert verified["observer"] == result["observer"]
    assert verified["configs"] == result["configs"]
    assert result["guard_events"] == []


@pytest.mark.parametrize("elapsed,write_delay,expected", [(299, 0, 0), (300, 0, 0), (301, 0, 1), (299, 2, 1)])
def test_one_parent_allocation_includes_final_write(tmp_path, monkeypatch, elapsed, write_delay, expected):
    spec, args = authority(tmp_path, monkeypatch)
    clock = [0.0]
    monkeypatch.setattr(d, "runtime", lambda: spec["contexts"]["checker"])
    monkeypatch.setattr(p.time, "monotonic", lambda: clock[0])
    monkeypatch.setattr(d, "cpu_guard", lambda: None)
    # Full structural authority retained; expensive pins delegated to fake child.
    real_write = s.write

    def writing(path, value):
        real_write(path, value)
        if Path(path).name == "preparation-receipt.json" and value["status"] != "running":
            clock[0] += write_delay

    monkeypatch.setattr(s, "write", writing)

    def supervise(argv, **kw):
        assert kw["allocation"] == 300 and kw["watchdog"] == 330 and kw["grace"] == 30
        kw["output"].mkdir()
        real_write(
            kw["output"] / "terminal.json",
            dict(status="candidate_requires_independent_review", study_launch_authorization_emitted=False),
        )
        clock[0] = elapsed
        return dict(status="completed", cleanup={"complete": True})

    monkeypatch.setattr(s, "supervise", supervise)
    assert p._run(args, 0.0) == expected
    receipt = d.read(Path(spec["output"]) / "preparation-receipt.json")
    assert receipt["status"] == (
        "allocation_exceeded" if expected else "candidate_requires_independent_review"
    )
    assert receipt["study_launch_authorization_emitted"] is False


def test_nested_phase_uses_remaining_deadline(tmp_path, monkeypatch):
    monkeypatch.setattr(p.time, "monotonic", lambda: 250)
    args = SimpleNamespace(spec="s", spec_sha256="h", authorization="a", authorization_sha256="h")
    spec = dict(workloads={})

    def supervise(argv, **kw):
        assert kw["allocation"] == kw["watchdog"] == 50 and kw["grace"] == 0
        kw["output"].mkdir()
        s.write(kw["output"] / "terminal.json", {"status": "synthetic"})
        return {"status": "completed"}

    monkeypatch.setattr(s, "supervise", supervise)
    assert p.phase(args, spec, "probe", tmp_path / "phase", str(tmp_path), 300) == {"status": "synthetic"}
    with pytest.raises(TimeoutError):
        p.phase(args, spec, "probe", tmp_path / "never", str(tmp_path), 250)


@pytest.mark.parametrize("offset", [0.0, 100.0, 10000.0])
def test_process_start_clock_conversion_is_conservative(monkeypatch, offset):
    monkeypatch.setattr(p.time, "monotonic", lambda: 500.0)
    monkeypatch.setattr(p.time, "clock_gettime", lambda kind: 500.01 + offset)
    monkeypatch.setattr(s, "identity", lambda pid: {"start_ticks": int((201 + offset) * 100)})
    monkeypatch.setattr(p.os, "sysconf", lambda name: 100)
    assert p.startup_monotonic() == pytest.approx(200.99)
    assert 500.0 - p.startup_monotonic() == pytest.approx(299.01)


def test_runtime_hashes_discovered_not_required_as_inputs(tmp_path, monkeypatch):
    spec, _ = authority(tmp_path, monkeypatch)
    new = tmp_path / "runtime/new-runtime.py"
    new.write_text("# prospective bytes")
    monkeypatch.setattr(d, "required_runtime_files", lambda: {str(new), spec["generation_source"]["path"]})
    monkeypatch.setattr(p, "generation_path", lambda: spec["generation_source"]["path"])
    monkeypatch.setattr(p.importlib.metadata, "distribution", lambda name: SimpleNamespace(files=[]))
    assert str(new) not in spec["pins"]
    discovered = p.discover_runtime(spec)
    assert discovered[str(new)] == d.digest(new)
    queries = p.queries_from_design(spec)
    arms = synthetic_discoveries(spec, queries)
    plan = p.assemble_plan(spec, queries, arms, discovered)
    assert plan["pins"][str(new)] == d.digest(new) and str(new) not in spec["pins"]
    spec["pins"][str(new)] = "wrong"
    with pytest.raises(ValueError):
        p.discover_runtime(spec)


@pytest.mark.parametrize("kind", ["unlisted", "alias_into_root", "symlink_out_of_root"])
def test_guard_denies_full_model_root_and_symlinks(tmp_path, kind):
    model = tmp_path / "model"
    model.mkdir()
    outside = tmp_path / "outside.bin"
    outside.write_text("tiny")
    member = model / "unlisted.bin"
    member.write_text("tiny")
    if kind == "alias_into_root":
        target = tmp_path / "alias.bin"
        target.symlink_to(member)
    elif kind == "symlink_out_of_root":
        target = model / "alias.bin"
        target.symlink_to(outside)
    else:
        target = member
    output = tmp_path / "probe"
    output.mkdir()
    code = "from benchmarks.research.polar_capacity_prepare import ProbeGuard;import sys;g=ProbeGuard(sys.argv[1],{},[],[sys.argv[2]]);g.install();open(sys.argv[3],'rb').read()"
    result = subprocess.run(
        [sys.executable, "-c", code, str(output), str(model), str(target)],
        cwd=Path(__file__).resolve().parents[1],
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        capture_output=True,
        timeout=10,
    )
    assert result.returncode == 79, result.stderr
    event = json.loads((output / "guard-events.jsonl").read_text())
    assert event["kind"] == "asset_read_denied"


def test_known_blocked_import_probe_is_logged_not_permitted(tmp_path):
    # Synthetic source/callsite authority; no socket operation is invoked.
    source = tmp_path / "connection.py"
    lines = ["# fixture"] * 137
    lines[124] = "def _has_ipv6(host):"
    lines[125] = "    try:"
    lines[126] = "        guard.audit('socket.bind', (sock, (host, 0)))"
    lines[127] = "    except RuntimeError: return False"
    lines[136] = "HAS_IPV6 = _has_ipv6('::1')"
    source.write_text("\n".join(lines) + "\n")
    output = tmp_path / "probe"
    output.mkdir()
    code = """import sys,socket,pathlib,types
from benchmarks.research import polar_capacity_prepare as p
path=pathlib.Path(sys.argv[2]);p.IPV6_SOURCE_SHA=p.d.digest(path)
guard=p.ProbeGuard(sys.argv[1],{str(path):p.IPV6_SOURCE_SHA},[])
sock=types.SimpleNamespace(family=socket.AF_INET6,type=socket.SOCK_STREAM,proto=0)
ns=dict(guard=guard,sock=sock);exec(compile(path.read_text(),str(path),'exec'),ns)
assert ns['HAS_IPV6'] is False
assert len(guard.events)==1 and guard.events[0]['expected_blocked_ipv6_probe'] is True
guard.stream.close()
"""
    result = subprocess.run(
        [sys.executable, "-c", code, str(output), str(source)],
        cwd=Path(__file__).resolve().parents[1],
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        capture_output=True,
        text=True,
        timeout=10,
    )
    assert result.returncode == 0, result.stderr
    event = json.loads((output / "guard-events.jsonl").read_text())
    assert event["kind"] == "network_denied" and event["expected_blocked_ipv6_probe"]


def test_child_dispatch_requires_actual_owned_parent(tmp_path):
    args = SimpleNamespace(
        command="coordinate",
        output=str(tmp_path / "coordinator"),
        deadline=100,
        input=None,
        expected=None,
        expected_sha256=None,
        role=None,
    )
    spec = dict(output=str(tmp_path), workloads={})
    Path(args.output).mkdir()
    s.write(
        Path(args.output) / "receipt.json",
        dict(status="running", supervisor_identity=dict(pid=999999, start_ticks=1)),
    )
    with pytest.raises(ValueError, match="live owned"):
        p.owned_entry(args, spec)


def test_owned_child_cannot_gain_supervisor_startup_time(tmp_path, monkeypatch):
    output = tmp_path / "coordinator"
    output.mkdir()
    parent = dict(pid=10, start_ticks=20)
    child = dict(pid=11, start_ticks=21)
    args = SimpleNamespace(
        command="coordinate",
        output=str(output),
        deadline=310.0,
        input=None,
        expected=None,
        expected_sha256=None,
        role=None,
        spec_sha256="bound",
    )
    spec = dict(output=str(tmp_path), workloads={})
    monkeypatch.setattr(p.os, "getppid", lambda: 10)
    monkeypatch.setattr(p.os, "getpid", lambda: 11)
    monkeypatch.setattr(s, "identity", lambda pid: parent if pid == 10 else child)
    monkeypatch.setattr(p.time, "monotonic", lambda: 50.0)
    monkeypatch.setattr(d, "workload_check", lambda blocked: [])
    s.write(
        output / "receipt.json",
        dict(
            status="running",
            supervisor_identity=parent,
            identity=child,
            started_monotonic=10.0,
            allocation_seconds=300.0,
            argv=[sys.executable, str(Path(p.__file__).resolve()), *sys.argv[1:]],
        ),
    )
    s.write(
        tmp_path / "preparation-receipt.json",
        dict(status="running", spec_sha256="bound", started_monotonic=0.0),
    )
    assert p.owned_entry(args, spec) == 300.0


@pytest.mark.parametrize("event", sorted(p.NETWORK_EVENTS - {"socket.bind"}))
def test_all_declared_network_audit_events_are_durably_denied(tmp_path, event):
    output = tmp_path / "probe"
    output.mkdir()
    code = "from benchmarks.research.polar_capacity_prepare import ProbeGuard;import sys;g=ProbeGuard(sys.argv[1],{},[]);g.audit(sys.argv[2],('synthetic-no-operation',))"
    result = subprocess.run(
        [sys.executable, "-c", code, str(output), event],
        cwd=Path(__file__).resolve().parents[1],
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        capture_output=True,
        timeout=10,
    )
    assert result.returncode == 79, result.stderr
    row = json.loads((output / "guard-events.jsonl").read_text())
    assert row["kind"] == "network_denied" and row["event"] == event
    assert row["expected_blocked_ipv6_probe"] is False
