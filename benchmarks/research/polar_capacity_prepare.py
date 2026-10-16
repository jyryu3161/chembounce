"""Bounded administrative discovery, never scientific launch authorization.

Only a separately reviewed spec/authorization permits actual preparation.
Imports here are administrative. Scientific imports occur only in owned probes.
"""

import argparse
from contextlib import contextmanager
import copy
import json
import importlib.metadata
import os
from pathlib import Path
import signal
import socket
import sys
import time

if __package__:
    from . import polar_capacity_development as d
    from . import polar_capacity_supervisor as s
else:
    import polar_capacity_development as d
    import polar_capacity_supervisor as s

SCHEMA = "polar-capacity-preparation-spec-v1"
AUTH = "polar-capacity-preparation-authorization-v1"
SUCCESS_SECONDS, WATCHDOG_SECONDS, CEILING_SECONDS = 300, 330, 360
IPV6_SOURCE_SHA = "2633bbdb69731e5ccb5cf4e4afd65605d86c7979cc5633126f50c92d5ad74a74"
NETWORK_EVENTS = {
    "socket.connect",
    "socket.bind",
    "socket.sendto",
    "socket.sendmsg",
    "socket.getaddrinfo",
    "socket.gethostbyname",
    "socket.gethostbyaddr",
    "socket.getnameinfo",
}


def exact(value, keys):
    if not isinstance(value, dict) or set(value) != set(keys.split()):
        raise ValueError("Unknown/missing preparation fields")


def absolute(value):
    if not isinstance(value, str) or not Path(value).is_absolute() or str(Path(value).resolve()) != value:
        raise ValueError("Preparation paths must be resolved absolute paths")
    return Path(value)


def validate_spec(spec, *, files=True):
    exact(
        spec,
        "schema design gate sources checker contexts pins python feature_definition generation_source output run_output workloads reviews clearance",
    )
    if spec["schema"] != SCHEMA or set(spec["sources"]) != set(d.ARMS):
        raise ValueError("Wrong preparation schema/source arms")
    if spec["design"]["sha256"] != d.DESIGN_SHA or spec["gate"]["sha256"] != d.GATE_SHA:
        raise ValueError("Wrong reviewed design gate")
    if set(spec["contexts"]) != {"baseline", "capacity", "checker"}:
        raise ValueError("Actual three runtime contexts required")
    paths = [absolute(spec[k]) for k in ("output", "run_output")] + [Path(spec["python"])]
    if not paths[2].is_absolute():
        raise ValueError("Absolute interpreter invocation required")
    if paths[2].absolute() != Path(sys.executable).absolute():
        raise ValueError("Actual interpreter invocation differs from preparation spec")
    if paths[0] == paths[1] or paths[0] in paths[1].parents or paths[1] in paths[0].parents:
        raise ValueError("Preparation and study outputs must be disjoint")
    if paths[1].exists():
        raise FileExistsError("Study output already exists")
    if (
        not spec["reviews"]
        or not spec["workloads"].get("identities")
        or not spec["workloads"].get("command_markers")
    ):
        raise ValueError("Review and current known-workload coverage required")
    required = {str(paths[2].resolve())}
    for record in [
        spec["design"],
        spec["gate"],
        spec["feature_definition"],
        spec["generation_source"],
        spec["clearance"],
        *spec["reviews"],
    ]:
        exact(record, "path sha256")
        absolute(record["path"])
        if spec["pins"].get(record["path"]) != record["sha256"]:
            raise ValueError("Missing preparation authority pin")
        required.add(record["path"])
    for role, source in [*spec["sources"].items(), ("checker", spec["checker"])]:
        exact(source, "path commit tree files")
        root = absolute(source["path"])
        if role in d.ARMS and source["commit"] != d.ARMS[role]:
            raise ValueError("Immutable comparison source changed")
        for output in paths[:2]:
            if output == root or root in output.parents:
                raise ValueError("Output must not modify a source checkout")
        expected = spec["contexts"][role]
        if (
            expected.get("cwd") != str(root)
            or expected.get("PYTHONPATH") != str(root)
            or expected.get("executable") != str(paths[2].resolve())
        ):
            raise ValueError("Runtime context differs from actual launch paths")
        if any(spec["pins"].get(p) != sha for p, sha in source["files"].items()):
            raise ValueError("Source pins incomplete")
        if files:
            d.tracked_inventory(source)
    checker_script = str(Path(spec["checker"]["path"]) / "benchmarks/research/polar_capacity_prepare.py")
    if checker_script != str(Path(__file__).resolve()) or checker_script not in spec["checker"]["files"]:
        raise ValueError("Preparation executable is outside checker authority")
    if not required <= spec["pins"].keys():
        raise ValueError("Missing executable/authority pin")
    if files:
        d.verify_pins(spec["pins"])
    return spec


def authorize(args, *, files=True):
    if d.digest(args.spec) != args.spec_sha256 or d.digest(args.authorization) != args.authorization_sha256:
        raise ValueError("Preparation spec/authorization digest differs")
    spec, auth = d.read(args.spec), d.read(args.authorization)
    exact(
        auth,
        "schema spec_sha256 source_review literal_review clearance allow_asset_hashing allow_import_discovery external_owner",
    )
    if auth["schema"] != AUTH or auth["spec_sha256"] != args.spec_sha256:
        raise ValueError("Not preparation authority")
    if auth["allow_asset_hashing"] is not True or auth["allow_import_discovery"] is not True:
        raise ValueError("Administrative capabilities not authorized")
    if (
        auth["source_review"] not in spec["reviews"]
        or auth["literal_review"] not in spec["reviews"]
        or auth["clearance"] != spec["clearance"]
    ):
        raise ValueError("Preparation approval evidence differs")
    owner = auth["external_owner"]
    exact(owner, "source review success_seconds watchdog_seconds cleanup_seconds scope")
    if (
        owner["success_seconds"] != 300
        or owner["watchdog_seconds"] != 330
        or owner["cleanup_seconds"] != 30
        or owner["scope"] != "owned-parent-and-descendants-through-cleanup-and-final-writes"
        or owner["review"] not in spec["reviews"]
    ):
        raise ValueError("Reviewed external whole-process enforcement is required")
    for record in (owner["source"], owner["review"]):
        exact(record, "path sha256")
        absolute(record["path"])
        if spec["pins"].get(record["path"]) != record["sha256"]:
            raise ValueError("External owner source/review pin missing or different")
    return validate_spec(spec, files=files)


def assemble_config(design, archive, model):
    """Build v2 policy explicitly; never import historical v1 scalar settings."""
    locked = design["locked_configuration"]
    config = {
        k: copy.deepcopy(locked[k])
        for k in (
            "recommendation_profile",
            "mode",
            "minimum_scores",
            "similarity_constraints",
            "ranking_policy",
            "route_search_policy",
            "route_search_cache_size",
            "max_design_attempts",
            "max_designs",
            "max_scored_products",
            "max_results",
            "search",
            "expansion",
        )
    }
    config.update(
        schema_version=1,
        catalog=design["assets"]["catalog"]["path"],
        rules=design["assets"]["rules"]["path"],
        provider={"kind": "saved", "path": archive},
        allocation_policy="balanced",
        property_thresholds={},
    )
    config["conformers"] = dict(
        scoring_objective="original-hit-v1",
        num_conformers=8,
        seed=42,
        energy_window=5.0,
        strain_penalty=0.05,
        max_iterations=500,
        pharmacophore_sigma=1.0,
        reject_unspecified_stereo=True,
    )
    config["refinement"] = dict(
        top_n=6,
        num_conformers=8,
        seeds=[42, 43, 44],
        threshold_policy="all_seeds",
        ranking_statistic="minimum",
        max_total_conformers=144,
        max_score_range=None,
        max_metric_range=None,
        adaptive_max_conformers=None,
        conformer_growth_factor=2,
        borderline_margin=0.0,
    )
    config["synthesis_evidence"] = dict(
        profile="computational-route-support-v1",
        policy={"require_forward_agreement": True},
        max_predictor_calls=24,
        max_input_tokens=12288,
        max_output_tokens=18432,
        predictor=dict(
            kind="reactiont5",
            model_path=model,
            model_revision="933114058cb2604dc1bf536dbebdfcefbe83d4fc",
            device="cpu",
            top_k=3,
            max_new_tokens=256,
            max_input_tokens=512,
        ),
    )
    d.validate_config(config, archive)
    return config


def queries_from_design(spec):
    design = d.read(spec["design"]["path"])
    for key, record in design["development_bindings"].items():
        if key != "slots" and spec["pins"].get(record["path"]) != record["sha256"]:
            raise ValueError("Missing historical collection/plan/diagnosis authority")
    queries = {}
    for slot in design["development_bindings"]["slots"]:
        if any(spec["pins"].get(p) != sha for p, sha in slot["bindings"].items()):
            raise ValueError("Full reviewed historical binding set required")

        def original(suffix):
            matches = [p for p in slot["bindings"] if p.endswith(suffix)]
            if len(matches) != 1:
                raise ValueError("Ambiguous historical binding")
            return matches[0]

        archive = original("/common/synthesis/proposals.json")
        request = original("/common/request.json")
        row = dict(
            request=request,
            design_pool=original("/native/design-pool.json"),
            proposals=archive,
            query=d.read(request)["query"],
            config=assemble_config(design, archive, slot["model_metadata"]["path"]),
        )
        d.validate_archive(row)
        queries[slot["historical_slot"]] = row
    if set(queries) != {"attempt-0000", "attempt-0009"}:
        raise ValueError("Wrong historical slice")
    return queries


class ProbeGuard:
    """Fixed-program defense in depth; not a native-code security sandbox.

    Original callable objects are never replaced by guards. Profiling blocks
    Python entrypoints; observed wrappers deny every started call before its
    original invocation. Arbitrary native extension calls may evade profiling.
    """

    def __init__(self, output, pins, denied_paths, denied_roots=()):
        self.stream = (Path(output) / "guard-events.jsonl").open("x")
        self.pins = pins
        denied_paths = tuple(denied_paths)
        self.denied = {str(Path(os.path.abspath(p))) for p in denied_paths} | {
            str(Path(p).resolve()) for p in denied_paths
        }
        self.roots = {Path(os.path.abspath(p)) for p in denied_roots} | {
            Path(p).resolve() for p in denied_roots
        }
        self.events = []
        self.phase = "imports"
        self.active = True
        self.previous = sys.getprofile()
        if self.previous is not None:
            raise ValueError("Unexpected preexisting profiler")

    def record(self, kind, **fields):
        row = dict(kind=kind, phase=self.phase, **fields)
        self.events.append(row)
        self.stream.write(json.dumps(row, sort_keys=True) + "\n")
        self.stream.flush()
        os.fsync(self.stream.fileno())

    def fatal(self, kind, **fields):
        self.record(kind, **fields)
        # A swallowed Python exception would disable setprofile and permit
        # subsequent work. Terminate this disposable probe before the call.
        os._exit(79)

    def audit(self, event, args):
        if not self.active:
            return
        if event == "open" and isinstance(args[0], (str, bytes, os.PathLike)):
            lexical = Path(os.path.abspath(os.fsdecode(args[0])))
            resolved = lexical.resolve()
            if any(str(path) in self.denied for path in (lexical, resolved)) or any(
                path == root or root in path.parents for path in (lexical, resolved) for root in self.roots
            ):
                self.fatal("asset_read_denied", path=str(resolved), lexical=str(lexical))
        if event in NETWORK_EVENTS:
            frame = sys._getframe(1)
            parent = frame.f_back
            sock = args[0] if event == "socket.bind" else None
            expected = bool(
                event == "socket.bind"
                and self.phase == "imports"
                and len(args) == 2
                and args[1] == ("::1", 0)
                and sock.family == socket.AF_INET6
                and sock.type == socket.SOCK_STREAM
                and sock.proto == 0
                and frame.f_code.co_name == "_has_ipv6"
                and frame.f_lineno == 127
                and parent
                and parent.f_code.co_name == "<module>"
                and parent.f_lineno == 137
                and not any(row["kind"] == "network_denied" for row in self.events)
                and self.pins.get(str(Path(frame.f_code.co_filename).resolve())) == IPV6_SOURCE_SHA
            )
            self.record(
                "network_denied",
                event=event,
                expected_blocked_ipv6_probe=expected,
                caller=[frame.f_code.co_filename, frame.f_code.co_name, frame.f_lineno],
                address=repr(args[1:]),
            )
            if not expected:
                os._exit(79)
            raise RuntimeError("IPv6 capability bind denied before operation")

    def profile(self, frame, event, arg):
        if not self.active:
            return
        module, name = frame.f_globals.get("__name__", ""), frame.f_code.co_name
        if event == "call":
            scientific = module.startswith("chembounce_synthesis") and (
                name
                in {
                    "run",
                    "score_candidate",
                    "refine_candidate",
                    "prepare_reference",
                    "ensemble",
                    "_overlay",
                    "predict_raw",
                    "predict",
                    "propose",
                    "apply",
                    "search",
                }
                or (
                    name == "__init__"
                    and type(frame.f_locals.get("self")).__name__
                    in {
                        "Catalog",
                        "Rules",
                        "ReactionT5Forward",
                        "PinnedReactionT5Forward",
                        "TokenReplay",
                        "SynLlamaProvider",
                    }
                )
            )
            loader = module.startswith(("transformers", "tokenizers", "safetensors")) and name in {
                "from_pretrained",
                "_from_pretrained",
                "generate",
                "load_file",
                "load_model",
                "__call__",
            }
            constructor = name == "__init__" and (
                module.startswith(("transformers.tokenization", "tokenizers"))
                or module == "torch.nn.modules.module"
            )
            if (
                scientific
                or loader
                or constructor
                or (module.startswith("torch") and name in {"load", "load_state_dict"})
            ):
                self.fatal("scientific_call_denied", module=module, name=name)
        elif event == "c_call":
            owner, native = getattr(arg, "__module__", "") or "", getattr(arg, "__name__", "")
            if native in {
                "EmbedMolecule",
                "EmbedMultipleConfs",
                "AlignMol",
                "GetBestAlignmentTransform",
                "GetAlignmentTransform",
                "safe_open",
            } or owner.startswith("tokenizers"):
                self.fatal("native_call_denied", module=owner, name=native)

    def install(self):
        sys.addaudithook(self.audit)
        sys.setprofile(self.profile)

    def close(self):
        if sys.getprofile() != self.profile:
            self.fatal("profiler_changed")
        sys.setprofile(self.previous)
        self.active = False
        self.stream.close()


def generation_path():
    return str(
        Path(
            importlib.metadata.distribution("transformers").locate_file("transformers/generation/utils.py")
        ).resolve()
    )


def discover_runtime(spec):
    """Prospective hashes; new hashes remain candidates for later review."""
    files = d.required_runtime_files()
    # Bind the transitive import capability probe and its distribution metadata.
    dist = importlib.metadata.distribution("urllib3")
    for relative in dist.files or ():
        path = Path(dist.locate_file(relative)).resolve()
        if path.is_file() and path.suffix != ".pyc":
            files.add(str(path))
    generation = spec["generation_source"]
    if (
        generation["path"] != generation_path()
        or generation["path"] not in files
        or d.digest(generation["path"]) != generation["sha256"]
    ):
        raise ValueError("Declared generation implementation differs")
    result = {path: d.digest(path) for path in sorted(files)}
    for path, sha in result.items():
        if path in spec["pins"] and spec["pins"][path] != sha:
            raise ValueError("Runtime bytes differ from existing authority")
    return result


def startup_monotonic():
    # /proc start ticks and CLOCK_BOOTTIME share an epoch. Subtract boottime
    # process age from the earlier monotonic reading; this is conservative
    # across read latency and any suspend since process startup.
    before = time.monotonic()
    boot = time.clock_gettime(time.CLOCK_BOOTTIME)
    start = s.identity(os.getpid())["start_ticks"] / os.sysconf("SC_CLK_TCK")
    if start > boot:
        raise ValueError("Process startup clock is inconsistent")
    return before - (boot - start)


def imported_file_pins():
    """Include observed transitive import files, beyond declared distributions."""
    paths = set()
    for module in list(sys.modules.values()):
        path = getattr(module, "__file__", None)
        if path:
            resolved = Path(path).resolve()
            if not resolved.is_file():
                raise ValueError("Imported file cannot be pinned: " + str(resolved))
            paths.add(str(resolved))
    return {path: d.digest(path) for path in sorted(paths)}


def probe(spec, role, queries, output, expected=None):
    source = {**spec["sources"][role], "runtime": spec["contexts"][role]}
    d.context_check(source, {**spec["checker"], "runtime": spec["contexts"]["checker"]})
    # Hashing installed runtime files is an authorized administrative discovery
    # stage, not model loading. Production assets are never hashed in probes.
    assets = {q["config"][k] for q in queries.values() for k in ("catalog", "rules")}
    assets |= {q["config"]["catalog"] + ".manifest.json" for q in queries.values()}
    model_roots = {q["config"]["synthesis_evidence"]["predictor"]["model_path"] for q in queries.values()}
    guard = ProbeGuard(output, dict(spec["pins"]), assets, model_roots)
    guard.phase = "runtime_hashing"
    guard.install()
    try:
        runtime_pins = discover_runtime(spec)
        guard.pins.update(runtime_pins)
        guard.phase = "imports"
        from chembounce_synthesis import workflow, conformers
        from rdkit import RDConfig

        if __package__:
            from .polar_capacity_observer import Observer, install_scientific
        else:
            from polar_capacity_observer import Observer, install_scientific

        class NoWorkObserver(Observer):
            def seal(self):
                if self.expected is None:
                    # Discovery has no prior inventory. This local seal only
                    # records installed instrumentation; it is not approval.
                    # A separate fresh verifier must reconstruct and match it,
                    # followed by independent literal review before execution.
                    self.expected = copy.deepcopy(self.inventory)
                super().seal()

            def emit(self, row):
                if row.get("event") == "started":
                    guard.fatal("observed_scientific_call_denied", observed_kind=row.get("kind"))
                super().emit(row)

        feature = str((Path(RDConfig.RDDataDir) / "BaseFeatures.fdef").resolve())
        if (
            feature != spec["feature_definition"]["path"]
            or d.digest(feature) != spec["feature_definition"]["sha256"]
        ):
            raise ValueError("Actual feature definition differs")
        ensemble = conformers.ensemble
        before = ensemble.cache_info()._asdict()
        configs = {}
        for key, query in queries.items():
            original = copy.deepcopy(query["config"])
            config = workflow._resolve_paths(workflow._config(copy.deepcopy(original)), Path.cwd())
            again = workflow._resolve_paths(workflow._config(copy.deepcopy(config)), Path.cwd())
            if config != again or original != query["config"]:
                raise ValueError("Non-idempotent/mutating config resolution")
            if workflow.effective_policy(config)["calibration"]["status"] != "engineering_policy_matched":
                raise ValueError("Engineering policy is not matched")
            d.validate_config(config, query["proposals"])
            configs[key] = config
        observer = NoWorkObserver(Path(output) / "observer-events.jsonl", expected_inventory=expected)
        try:
            install_scientific(observer)
            originals = list(observer.restores)
            actual = copy.deepcopy(observer.inventory)
        finally:
            observer.close()
        if any(getattr(owner, name) is not original for owner, name, original, wrapped in originals):
            raise ValueError("Observer did not restore original aliases")
        if (
            conformers.ensemble is not ensemble
            or ensemble.cache_info()._asdict() != before
            or observer.events
        ):
            raise ValueError("Probe changed cache identity/state or performed work")
        guard.phase = "inventory"
        for path, sha in imported_file_pins().items():
            if path in runtime_pins and runtime_pins[path] != sha:
                raise ValueError("Runtime import changed during probe")
            if path in spec["pins"] and spec["pins"][path] != sha:
                raise ValueError("Imported bytes differ from existing authority")
            runtime_pins[path] = sha
        d.context_check(source, {**spec["checker"], "runtime": spec["contexts"]["checker"]})
        return dict(
            status="discovered_not_approved",
            runtime=d.runtime(),
            runtime_files=sorted(runtime_pins),
            runtime_pins=runtime_pins,
            configs=configs,
            observer=actual,
            cache_before=before,
            cache_after=ensemble.cache_info()._asdict(),
            restored_aliases=[owner.__name__ + "." + name for owner, name, _, _ in originals],
            scientific_calls=0,
            guard_events=guard.events,
            guard_limit="Fixed reviewed import/config/introspection program; profiling is not a native-code sandbox",
        )
    finally:
        guard.close()


def source_record(spec, role):
    original = spec["checker"] if role == "checker" else spec["sources"][role]
    return {
        **original,
        "runtime": spec["contexts"][role],
        "synthesis_source": {
            Path(p).name: sha
            for p, sha in original["files"].items()
            if Path(p).parent == Path(original["path"]) / "chembounce_synthesis" and p.endswith(".py")
        },
    }


def assemble_plan(spec, queries, discoveries, runtime_pins=None):
    if discoveries["baseline"]["configs"] != discoveries["capacity"]["configs"]:
        raise ValueError("Arms resolved different scientific configs")
    queries = copy.deepcopy(queries)
    for key in queries:
        queries[key]["config"] = discoveries["baseline"]["configs"][key]
    pins = dict(spec["pins"])
    for mapping in [runtime_pins or {}, *(value.get("runtime_pins", {}) for value in discoveries.values())]:
        for path, sha in mapping.items():
            if path in pins and pins[path] != sha:
                raise ValueError("Discovered runtime conflicts with preexisting authority")
            pins[path] = sha
    return dict(
        schema="polar-capacity-development-plan-v2",
        design=spec["design"],
        gate=spec["gate"],
        pins=pins,
        sources={r: source_record(spec, r) for r in d.ARMS},
        checker=source_record(spec, "checker"),
        queries=queries,
        observer={r: discoveries[r]["observer"] for r in d.ARMS},
        python=spec["python"],
        output=spec["run_output"],
        workloads=spec["workloads"],
        administrative_runtime=spec["contexts"]["checker"],
        feature_definition=spec["feature_definition"],
    )


def owned_entry(args, spec):
    base = Path(spec["output"]) / "coordinator"
    allowed = {
        "coordinate": {base},
        "handoff": {base / "handoff"},
        "probe": {base / (str(args.role) + suffix) for suffix in ("-discover", "-verify")},
    }
    if Path(args.output) not in allowed[args.command]:
        raise ValueError("Unexpected owned preparation output")
    if args.command in {"probe", "handoff"}:
        expected_input = base / ("queries.json" if args.command == "probe" else "candidate-plan.json")
        if args.input != str(expected_input):
            raise ValueError("Unexpected owned preparation input")
    if args.command == "probe":
        verifying = Path(args.output).name.endswith("-verify")
        expected = str(base / (str(args.role) + "-expected.json")) if verifying else None
        if args.role not in d.ARMS or args.expected != expected or bool(args.expected_sha256) != verifying:
            raise ValueError("Unexpected discovery/verification authority")
    receipt = d.read(Path(args.output) / "receipt.json")
    if (
        receipt.get("status") != "running"
        or not s.same(s.identity(os.getppid()), receipt.get("supervisor_identity"))
        or not s.same(s.identity(os.getpid()), receipt.get("identity"))
    ):
        raise ValueError("No live owned preparation parent")
    if args.deadline != receipt["started_monotonic"] + receipt["allocation_seconds"] or receipt["argv"][
        1:
    ] != [str(Path(__file__).resolve()), *sys.argv[1:]]:
        raise ValueError("Preparation child command/deadline differs")
    preparation = d.read(Path(spec["output"]) / "preparation-receipt.json")
    if preparation.get("status") != "running" or preparation.get("spec_sha256") != args.spec_sha256:
        raise ValueError("No active whole-preparation authority")
    # Supervisor startup may shift its own local deadline by a few ticks;
    # never transfer that overhead into the preparation's successful lease.
    deadline = min(args.deadline, preparation["started_monotonic"] + SUCCESS_SECONDS)
    if time.monotonic() >= deadline or d.workload_check(spec["workloads"]):
        raise ValueError("Preparation lease expired or known workload live")
    return deadline


def child_args(args, command, output, *extra):
    return [
        sys.executable,
        str(Path(__file__).resolve()),
        command,
        "--spec",
        args.spec,
        "--spec-sha256",
        args.spec_sha256,
        "--authorization",
        args.authorization,
        "--authorization-sha256",
        args.authorization_sha256,
        "--output",
        str(output),
        "--deadline",
        "{deadline}",
        *extra,
    ]


def phase(args, spec, command, directory, source, deadline, *extra):
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise TimeoutError("Preparation allocation exhausted")
    result = s.supervise(
        child_args(args, command, directory, *extra),
        cwd=source,
        env=d.environment(source),
        output=directory,
        allocation=remaining,
        watchdog=remaining,
        grace=0,
        blocked=spec["workloads"],
    )
    if result["status"] != "completed" or time.monotonic() > deadline:
        raise RuntimeError("Administrative child failed or exceeded shared allocation")
    return d.read(directory / "terminal.json")


def coordinator(args, spec):
    output = Path(args.output)
    runtime_pins = discover_runtime(spec)
    queries = queries_from_design(spec)
    qpath = output / "queries.json"
    s.write(qpath, queries)
    discoveries = {}
    for role in d.ARMS:
        common = ("--role", role, "--input", str(qpath), "--input-sha256", d.digest(qpath))
        first = phase(
            args,
            spec,
            "probe",
            output / (role + "-discover"),
            spec["sources"][role]["path"],
            args.deadline,
            *common,
        )
        expected = output / (role + "-expected.json")
        s.write(expected, first["observer"])
        second = phase(
            args,
            spec,
            "probe",
            output / (role + "-verify"),
            spec["sources"][role]["path"],
            args.deadline,
            *common,
            "--expected",
            str(expected),
            "--expected-sha256",
            d.digest(expected),
        )
        for key in (
            "configs",
            "observer",
            "runtime",
            "runtime_files",
            "runtime_pins",
            "cache_before",
            "cache_after",
            "scientific_calls",
        ):
            if first[key] != second[key]:
                raise ValueError("Independent probe discovery/verification differs: " + key)
        discoveries[role] = first
    plan = assemble_plan(spec, queries, discoveries, runtime_pins)
    d.validate_plan(plan)
    path = output / "candidate-plan.json"
    s.write(path, plan)
    phase(
        args,
        spec,
        "handoff",
        output / "handoff",
        spec["checker"]["path"],
        args.deadline,
        "--input",
        str(path),
        "--input-sha256",
        d.digest(path),
    )
    d.verify_pins(plan["pins"])
    if Path(spec["run_output"]).exists():
        raise ValueError("Preparation touched study output")
    return dict(
        status="candidate_requires_independent_review",
        candidate_plan=str(path),
        candidate_sha256=d.digest(path),
        study_launch_authorization_emitted=False,
        prepared_ledger=str(output / "handoff/existing-prepare/ledger.json"),
        pins_before_after_matched=True,
        deadline=args.deadline,
    )


@contextmanager
def administrative_alarm(deadline):
    state = {"expired": False}

    def expired(number, frame):
        state["expired"] = True
        raise TimeoutError("Administrative preparation360-second ceiling exhausted")

    old = signal.signal(signal.SIGALRM, expired)
    timer = signal.setitimer(signal.ITIMER_REAL, max(0.001, deadline - time.monotonic()))
    if timer != (0.0, 0.0):
        signal.setitimer(signal.ITIMER_REAL, *timer)
        signal.signal(signal.SIGALRM, old)
        raise ValueError("Preexisting alarm conflicts with preparation")
    try:
        yield state
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


def run(args):
    # Linux process start includes imports and parent startup; tick flooring is
    # conservative, unlike beginning the allocation after validation/hashing.
    started = startup_monotonic()
    with administrative_alarm(started + CEILING_SECONDS) as expiry:
        return _run(args, started, expiry)


def _run(args, started, expiry=None):
    expiry = expiry if expiry is not None else {"expired": False}

    def expired():
        # Sticky signal state survives a supervising layer catching TimeoutError.
        # No new parent I/O/work after expiry; the reviewed external owner must
        # own ceiling enforcement, including parent persistence and cleanup.
        return expiry["expired"] or time.monotonic() >= started + CEILING_SECONDS

    if expired():
        return 1
    d.cpu_guard()
    if expired():
        return 1
    spec = authorize(args, files=False)
    if expired():
        return 1
    if d.runtime() != spec["contexts"]["checker"]:
        raise ValueError("Parent runtime/cwd differs from declared preparation context")
    if expired():
        return 1
    output = absolute(args.output)
    if str(output) != spec["output"]:
        raise ValueError("Preparation output differs")
    if expired():
        return 1
    output.mkdir(exist_ok=False)
    receipt = {
        "status": "running",
        "started_monotonic": started,
        "success_allocation_seconds": 300,
        "watchdog_seconds": 330,
        "cleanup_seconds": 30,
        "spec_sha256": args.spec_sha256,
        "study_launch_authorization_emitted": False,
    }

    def persist():
        if expired():
            raise TimeoutError("No persistence after administrative ceiling expiry")
        try:
            s.write(output / "preparation-receipt.json", receipt)
        except BaseException:
            receipt["status"] = "persistence_interrupted"
            if expired():
                raise
            s.write(output / "preparation-receipt.json", receipt)
            raise

    try:
        persist()
        if expired():
            return 1
        remaining = started + SUCCESS_SECONDS - time.monotonic()
        if remaining <= 0:
            raise TimeoutError("No preparation allocation remains")
        child = s.supervise(
            child_args(args, "coordinate", output / "coordinator"),
            cwd=spec["checker"]["path"],
            env=d.environment(spec["checker"]["path"]),
            output=output / "coordinator",
            allocation=remaining,
            watchdog=started + WATCHDOG_SECONDS - time.monotonic(),
            grace=30,
            blocked=spec["workloads"],
        )
        if expired():
            return 1
        receipt["child"] = child
        receipt["status"] = (
            "candidate_requires_independent_review" if child["status"] == "completed" else "failed"
        )
        if child["status"] == "completed":
            result = d.read(output / "coordinator/terminal.json")
            if (
                result.get("status") != "candidate_requires_independent_review"
                or result.get("study_launch_authorization_emitted") is not False
            ):
                raise ValueError("Coordinator returned invalid authority")
            receipt["candidate"] = result
    except BaseException as exc:
        receipt["status"] = "failed"
        receipt["error"] = dict(type=type(exc).__name__, message=str(exc))
    finally:
        if expired():
            # Keep the previous durable running/partial evidence. Neither an
            # alarm nor failed cleanup grants time for another correction write.
            return 1
        receipt["parent_elapsed_before_final_write"] = time.monotonic() - started
        if receipt["parent_elapsed_before_final_write"] > SUCCESS_SECONDS:
            receipt["status"] = "allocation_exceeded"
        try:
            persist()
        except BaseException:
            if expired():
                return 1
            raise
        if expired():
            return 1
        if time.monotonic() - started > SUCCESS_SECONDS:
            receipt["status"] = "allocation_exceeded"
            receipt["parent_elapsed_after_write"] = time.monotonic() - started
            try:
                persist()
            except BaseException:
                if expired():
                    return 1
                raise
            if expired():
                return 1
    return 0 if receipt["status"] == "candidate_requires_independent_review" else 1


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["prepare", "coordinate", "probe", "handoff"])
    for name in ("spec", "spec-sha256", "authorization", "authorization-sha256", "output"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--deadline", type=float)
    parser.add_argument("--role", choices=list(d.ARMS))
    for name in ("input", "input-sha256", "expected", "expected-sha256"):
        parser.add_argument("--" + name)
    args = parser.parse_args(argv)
    if args.command == "prepare":
        return run(args)
    # Authenticate ownership before expensive hashing/imports or output mutation.
    spec = authorize(args, files=False)
    args.deadline = owned_entry(args, spec)
    output = Path(args.output)
    if (output / "terminal.json").exists():
        raise FileExistsError("Preparation child terminal already exists")
    result = {"status": "failed"}
    try:
        # Probe installs its logging guard before scientific imports. Other
        # administrative paths use the existing network-denial/environment check.
        if args.command == "probe":
            if (
                os.environ.get("CUDA_VISIBLE_DEVICES") != ""
                or any(os.environ.get(k) != "1" for k in d.THREADS)
                or any(os.environ.get(k) != "1" for k in ("HF_HUB_OFFLINE", "TRANSFORMERS_OFFLINE"))
            ):
                raise ValueError("Probe requires CPU1/offline environment")
        else:
            d.cpu_guard()
        # The coordinator/handoff explicitly hash authorized asset bytes.
        # Import probes only check source/runtime authority and install their
        # asset-opening guard; they do not rehash production assets.
        if args.command != "probe":
            authorize(args)
        else:
            d.tracked_inventory(spec["sources"][args.role])
            d.tracked_inventory(spec["checker"])
        if args.command == "coordinate":
            result = coordinator(args, spec)
        else:
            if not args.input or d.digest(args.input) != args.input_sha256:
                raise ValueError("Bound child input differs")
            if args.command == "probe":
                expected = None
                if args.expected:
                    if d.digest(args.expected) != args.expected_sha256:
                        raise ValueError("Expected inventory digest differs")
                    expected = d.read(args.expected)
                result = probe(spec, args.role, d.read(args.input), output, expected)
            else:
                target = output / "existing-prepare"
                code = d.main(
                    [
                        "prepare",
                        "--plan",
                        args.input,
                        "--plan-sha256",
                        args.input_sha256,
                        "--output",
                        str(target),
                    ]
                )
                if code != 0 or d.read(target / "terminal.json").get("status") != "prepared_not_authorized":
                    raise ValueError("Existing preparation failed")
                ledger = d.read(target / "ledger.json")
                plan = d.read(args.input)
                expected_rows = [
                    dict(index=i, status="not_started", request=d.request_for(plan, args.input_sha256, i))
                    for i in range(4)
                ]
                if ledger != dict(plan_sha256=args.input_sha256, rows=expected_rows):
                    raise ValueError("Existing prepare ledger differs")
                result = dict(status="prepared_not_authorized", plan_sha256=args.input_sha256)
        if time.monotonic() > args.deadline:
            raise TimeoutError("Preparation child exceeded shared deadline")
    except BaseException as exc:
        result = dict(status="failed", error=dict(type=type(exc).__name__, message=str(exc)))
    s.write(output / "terminal.json", result)
    return 1 if result["status"] == "failed" else 0


if __name__ == "__main__":
    raise SystemExit(main())
