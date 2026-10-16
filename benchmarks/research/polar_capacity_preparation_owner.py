"""Fixed external owner for administrative preparation; never a study launcher.

Only the reviewed, hash-bound preparation spec/authorization is accepted. The
owner clock includes its process startup and all finalization; it is not the
child's lifetime. No scientific imports or asset constructors occur here.
"""

import argparse
from contextlib import contextmanager
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

if __package__:
    from . import polar_capacity_prepare as p
else:
    import polar_capacity_prepare as p

d, s = p.d, p.s
LITERAL_SCHEMA = "polar-capacity-preparation-literal-review-v1"
REVIEW_STATUS = "approved_authority_payload"
OWNER_SCHEMA = "polar-capacity-preparation-owner-v1"
ENV_KEYS = (
    *d.THREADS,
    "CUDA_VISIBLE_DEVICES",
    "HF_HUB_OFFLINE",
    "TRANSFORMERS_OFFLINE",
    "PYTHONDONTWRITEBYTECODE",
    "PYTHONPATH",
)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def binding(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def review_payload(spec_core, auth_core, spec_path, authorization_path):
    """Prospective payload from envelope-free records: no fabricated hash needed."""
    invocation = {"spec": str(p.absolute(spec_path)), "authorization": str(p.absolute(authorization_path))}
    literal_path = str(Path(invocation["authorization"]).with_name("literal-review.json"))
    if (
        invocation["spec"] == invocation["authorization"]
        or "spec_sha256" in auth_core
        or "literal_review" in auth_core
        or literal_path in spec_core["pins"]
        or any(row.get("path") == literal_path for row in spec_core["reviews"])
    ):
        raise ValueError("Prospective payload must exclude exactly the literal-review envelope")
    result = {
        "spec": copy.deepcopy(spec_core),
        "authorization": copy.deepcopy(auth_core),
        "invocation": invocation,
    }
    canonical(result)  # reject nonfinite content before independent review
    return result


def authority_payload(spec, auth, spec_path, authorization_path):
    """Remove only the designated literal-review envelope, never other values.

    The independently reviewed payload precedes final spec/auth byte hashes.
    Literal review content is validated separately against this exact digest.
    """
    spec, auth = copy.deepcopy(spec), copy.deepcopy(auth)
    invocation = {"spec": str(p.absolute(spec_path)), "authorization": str(p.absolute(authorization_path))}
    if invocation["spec"] == invocation["authorization"]:
        raise ValueError("Spec and authorization paths must differ")
    review = auth["literal_review"]
    p.exact(review, "path sha256")
    path = str(p.absolute(review["path"]))
    if sum(item.get("path") == path for item in spec["reviews"]) != 1 or review not in spec["reviews"]:
        raise ValueError("Exactly one designated literal review is required")
    protected = {
        spec[key]["path"]
        for key in ("design", "gate", "clearance", "feature_definition", "generation_source")
    }
    protected.update(invocation.values())
    protected.update(
        (
            auth["source_review"]["path"],
            auth["external_owner"]["source"]["path"],
            auth["external_owner"]["review"]["path"],
        )
    )
    protected.add(str(Path(spec["python"]).resolve()))
    for source in [*spec["sources"].values(), spec["checker"]]:
        protected.update(source["files"])
    if path in protected or spec["pins"].get(path) != review["sha256"]:
        raise ValueError("Literal review collides with protected authority or lacks its exact pin")
    # Arbitrary scientific/asset pins cannot be reused as the review envelope:
    # the payload review must be a dedicated file beside the declared auth.
    expected = str(Path(invocation["authorization"]).with_name("literal-review.json"))
    if path != expected:
        raise ValueError("Literal review must use the dedicated invocation-adjacent path")
    spec["reviews"].remove(review)
    del spec["pins"][path]
    del auth["spec_sha256"]
    del auth["literal_review"]
    return review_payload(spec, auth, invocation["spec"], invocation["authorization"])


def paths(spec):
    preparation = p.absolute(spec["output"])
    owner = preparation.with_name(preparation.name + ".owner")
    return owner, owner / "adapter"


def environment_record(source):
    expected = d.environment(source)
    result = {key: os.environ.get(key) for key in ENV_KEYS}
    if result != {key: expected[key] for key in ENV_KEYS}:
        raise ValueError("Owner/adapter requires exact CPU1, CUDA-hidden, offline environment")
    return result


def source_authority(spec, auth):
    root = Path(spec["checker"]["path"])
    if (
        p.s is not s
        or p.d is not d
        or d.supervise is not s.supervise
        or s.supervise.__globals__ is not s.__dict__
    ):
        raise ValueError("Preparation and owner must share the exact supervisor module")
    own = str(Path(__file__).resolve())
    record = auth["external_owner"]["source"]
    if (
        record != {"path": own, "sha256": d.digest(own)}
        or spec["checker"]["files"].get(own) != record["sha256"]
    ):
        raise ValueError("Executing owner is outside reviewed source authority")
    for module, name in (
        (p, "polar_capacity_prepare.py"),
        (d, "polar_capacity_development.py"),
        (s, "polar_capacity_supervisor.py"),
    ):
        actual = str(Path(module.__file__).resolve())
        if actual != str(root / "benchmarks/research" / name) or spec["pins"].get(actual) != d.digest(actual):
            raise ValueError("Imported preparation/tooling source differs")
    d.tracked_inventory(spec["checker"])


def authorize(args):
    spec = p.authorize(args, files=False)
    auth = d.read(args.authorization)
    source_authority(spec, auth)
    if d.runtime() != spec["contexts"]["checker"]:
        raise ValueError("Actual owner/adapter runtime context differs")
    environment_record(spec["checker"]["path"])
    if args.output != spec["output"]:
        raise ValueError("Preparation output differs from authority")
    payload = authority_payload(spec, auth, args.spec, args.authorization)
    literal = auth["literal_review"]
    if d.digest(literal["path"]) != literal["sha256"]:
        raise ValueError("Literal review bytes differ")
    review = d.read(literal["path"])
    p.exact(review, "schema status payload_sha256")
    if review != {"schema": LITERAL_SCHEMA, "status": REVIEW_STATUS, "payload_sha256": binding(payload)}:
        raise ValueError("Final authority does not reconstruct reviewed payload")
    for record in (auth["source_review"], auth["external_owner"]["review"], auth["clearance"]):
        if d.digest(record["path"]) != record["sha256"]:
            raise ValueError("Source review/current clearance bytes differ")
    owner, _ = paths(spec)
    for other in (
        Path(spec["run_output"]),
        Path(spec["output"]),
        *[Path(x["path"]) for x in [*spec["sources"].values(), spec["checker"]]],
    ):
        if owner == other or owner in other.parents or other in owner.parents:
            raise ValueError("Owner output overlaps source/preparation/study paths")
    if Path(spec["run_output"]).exists():
        raise FileExistsError("Study output must remain absent")
    return spec, binding(payload)


class OwnerClock:
    """First supervisor t0 is owner startup; later values are actual monotonic.

    Supervisor setup is deliberately charged to the one owner allocation. This
    scope is neither child wall time nor an invented earlier child launch.
    """

    def __init__(self, anchor, clock=None):
        self.anchor, self.clock, self.first = anchor, clock or time.monotonic, True

    def __call__(self):
        if self.first:
            self.first = False
            return self.anchor
        return self.clock()


def adapter_argv(args, anchor):
    return [
        sys.executable,
        str(Path(__file__).resolve()),
        "adapter",
        "--spec",
        args.spec,
        "--spec-sha256",
        args.spec_sha256,
        "--authorization",
        args.authorization,
        "--authorization-sha256",
        args.authorization_sha256,
        "--output",
        args.output,
        "--anchor",
        repr(anchor),
    ]


def owner_argv(args):
    command = adapter_argv(args, 0.0)
    command[2] = "run"
    return command[:-2]


def process_context(pid):
    root = Path("/proc") / str(pid)
    return {
        "argv": [os.fsdecode(item) for item in (root / "cmdline").read_bytes().split(b"\0") if item],
        "cwd": str((root / "cwd").resolve()),
    }


def active_identity(pid):
    result = s.identity(pid)
    if result is None or result["state"] == "Z":
        raise ValueError("Live owned identity is unavailable")
    return result


def adapter_authority(args, spec, payload_sha):
    owner, directory = paths(spec)
    owner_receipt = d.read(owner / "owner-receipt.json")
    receipt = d.read(directory / "receipt.json")
    anchor = args.anchor
    if not isinstance(anchor, float) or not math.isfinite(anchor) or anchor < 0:
        raise ValueError("Invalid owner anchor")
    parent, current = active_identity(os.getppid()), active_identity(os.getpid())
    expected_argv = adapter_argv(args, anchor)
    parent_context, child_context = process_context(parent["pid"]), process_context(current["pid"])
    actual_argv = [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]]
    if (
        owner_receipt.get("status") != "running"
        or owner_receipt.get("schema") != OWNER_SCHEMA
        or owner_receipt.get("spec_sha256") != args.spec_sha256
        or owner_receipt.get("authorization_sha256") != args.authorization_sha256
        or owner_receipt.get("payload_sha256") != payload_sha
        or owner_receipt.get("started_monotonic") != anchor
        or owner_receipt.get("adapter_argv") != expected_argv
        or owner_receipt.get("owner_argv") != owner_argv(args)
        or parent_context != {"argv": owner_argv(args), "cwd": spec["checker"]["path"]}
        or child_context != {"argv": expected_argv, "cwd": spec["checker"]["path"]}
        or owner_receipt.get("environment") != environment_record(spec["checker"]["path"])
        or owner_receipt.get("environment_sha256") != binding(dict(os.environ))
        or not s.same(parent, owner_receipt.get("identity"))
        or not s.same(parent, receipt.get("supervisor_identity"))
        or not s.same(current, receipt.get("identity"))
        or current["ppid"] != parent["pid"]
        or receipt.get("status") != "running"
        or receipt.get("started_monotonic") != anchor
        or receipt.get("allocation_seconds") != 300
        or receipt.get("watchdog_seconds") != 330
        or receipt.get("kill_grace_seconds") != 30
        or receipt.get("argv") != expected_argv
        or actual_argv != expected_argv
        or receipt.get("cwd") != spec["checker"]["path"]
        or str(Path.cwd().resolve()) != spec["checker"]["path"]
        or Path(sys.executable).absolute() != Path(spec["python"]).absolute()
    ):
        raise ValueError("Adapter lacks exact live owner/command/context authority")
    if time.monotonic() >= anchor + 300 or s.workload_check(spec["workloads"]):
        raise ValueError("Owner allocation exhausted or known workload live")
    if Path(spec["output"]).exists() or (directory / "terminal.json").exists():
        raise FileExistsError("Preparation/adapter output reuse")
    return anchor


def complete(receipt, allocation):
    def finite(value):
        return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)

    resources = receipt.get("resources", {})
    identity = receipt.get("identity")
    return (
        receipt.get("status") == "completed"
        and type(receipt.get("returncode")) is int
        and receipt["returncode"] == 0
        and receipt.get("allocation_seconds") == allocation
        and finite(allocation)
        and 0 < allocation <= 300
        and receipt.get("cleanup", {}).get("complete") is True
        and receipt["cleanup"].get("remaining") == []
        and finite(receipt.get("started_monotonic"))
        and finite(receipt.get("wall_seconds"))
        and 0 <= receipt["wall_seconds"] <= allocation
        and isinstance(identity, dict)
        and identity.get("pid", 0) > 0
        and identity.get("start_ticks", -1) >= 0
        and isinstance(receipt.get("terminal_sha256"), str)
        and len(receipt["terminal_sha256"]) == 64
        and resources.get("samples", 0) > 0
        and resources.get("gaps") == []
        and finite(resources.get("sampled_tree_peak_rss_bytes"))
        and resources["sampled_tree_peak_rss_bytes"] >= 0
    )


def validate_prepared(args, spec, anchor, *, check_files=True):
    """Bind actual existing prepare artifacts, not a success flag or fake plan."""
    base = Path(spec["output"])
    coordinator = base / "coordinator"
    handoff = coordinator / "handoff"
    existing = handoff / "existing-prepare"
    files = {}

    def read(path):
        path = Path(path)
        files[str(path)] = d.digest(path)
        return d.read(path)

    receipt = read(base / "preparation-receipt.json")
    terminal = read(coordinator / "terminal.json")
    coordinator_receipt = read(coordinator / "receipt.json")
    if (
        receipt.get("status") != "candidate_requires_independent_review"
        or receipt.get("started_monotonic") != anchor
        or receipt.get("spec_sha256") != args.spec_sha256
        or receipt.get("success_allocation_seconds") != 300
        or receipt.get("watchdog_seconds") != 330
        or receipt.get("cleanup_seconds") != 30
        or receipt.get("study_launch_authorization_emitted") is not False
        or receipt.get("child") != coordinator_receipt
        or receipt.get("candidate") != terminal
        or not complete(coordinator_receipt, coordinator_receipt.get("allocation_seconds", -1))
        or coordinator_receipt["started_monotonic"] < anchor
        or coordinator_receipt["started_monotonic"] + coordinator_receipt["wall_seconds"] > anchor + 300
        or coordinator_receipt.get("terminal_sha256") != files[str(coordinator / "terminal.json")]
        or not 0 <= receipt.get("parent_elapsed_before_final_write", float("inf")) <= 300
    ):
        raise ValueError("Preparation receipt/coordinator evidence is incomplete")
    plan_path = coordinator / "candidate-plan.json"
    plan_sha = d.digest(plan_path)
    files[str(plan_path)] = plan_sha
    if terminal != dict(
        status="candidate_requires_independent_review",
        candidate_plan=str(plan_path),
        candidate_sha256=plan_sha,
        study_launch_authorization_emitted=False,
        prepared_ledger=str(existing / "ledger.json"),
        pins_before_after_matched=True,
        deadline=anchor + 300,
    ):
        raise ValueError("Coordinator candidate/ledger/deadline binding differs")
    plan = (
        d.load_plan(plan_path, plan_sha)
        if check_files
        else d.validate_plan(d.read(plan_path), check_files=False)
    )
    for key, expected in {
        "design": spec["design"],
        "gate": spec["gate"],
        "python": spec["python"],
        "output": spec["run_output"],
        "workloads": spec["workloads"],
        "administrative_runtime": spec["contexts"]["checker"],
        "feature_definition": spec["feature_definition"],
        "sources": {role: p.source_record(spec, role) for role in d.ARMS},
        "checker": p.source_record(spec, "checker"),
    }.items():
        if plan[key] != expected:
            raise ValueError("Prepared plan differs from input authority: " + key)
    if any(plan["pins"].get(path) != sha for path, sha in spec["pins"].items()):
        raise ValueError("Prepared plan changed protected input pins")
    existing_terminal = read(existing / "terminal.json")
    p.exact(existing_terminal, "status plan_sha256 process_lifetime_ru_maxrss_kib")
    usage = existing_terminal.pop("process_lifetime_ru_maxrss_kib")
    if not isinstance(usage, (int, float)) or not math.isfinite(usage) or usage < 0:
        raise ValueError("Existing prepare lifetime RSS metadata is invalid")
    handoff_terminal = read(handoff / "terminal.json")
    handoff_receipt = read(handoff / "receipt.json")
    if (
        handoff_terminal != {"status": "prepared_not_authorized", "plan_sha256": plan_sha}
        or not complete(handoff_receipt, handoff_receipt.get("allocation_seconds", -1))
        or handoff_receipt["started_monotonic"] < anchor
        or handoff_receipt["started_monotonic"] + handoff_receipt["wall_seconds"] > anchor + 300
        or handoff_receipt.get("terminal_sha256") != files[str(handoff / "terminal.json")]
        or existing_terminal != handoff_terminal
        or read(existing / "ledger.json")
        != dict(
            plan_sha256=plan_sha,
            rows=[
                dict(index=i, status="not_started", request=d.request_for(plan, plan_sha, i))
                for i in range(4)
            ],
        )
        or Path(spec["run_output"]).exists()
    ):
        raise ValueError("Existing prepare handoff/four-slot ledger differs")
    return dict(
        candidate_plan=str(plan_path),
        candidate_sha256=plan_sha,
        files_sha256=files,
        planned_slots=4,
        study_launch_authorization_emitted=False,
    )


class CeilingExpired(TimeoutError):
    pass


@contextmanager
def bounded_persistence(anchor, expiry):
    """Owner-local guard around the existing writer, including supervise finally.

    Does not change bytes/semantics of an allowed write. Restore the original
    callable even on interruption. Kernel calls already in progress may block;
    unsafe cleanup or missing final authority is never success.
    """
    original = s.write

    def expired():
        return expiry["expired"] or time.monotonic() >= anchor + 360

    def write(*args, **kwargs):
        if expired():
            raise CeilingExpired("No new persistence after owner ceiling expiry")
        return original(*args, **kwargs)

    s.write = write
    try:
        yield expired
    finally:
        s.write = original


def adapter(args):
    if args.anchor is None or not math.isfinite(args.anchor) or time.monotonic() >= args.anchor + 300:
        raise ValueError("No live adapter allocation")
    spec, payload = authorize(args)
    anchor = adapter_authority(args, spec, payload)
    _, directory = paths(spec)
    with p.administrative_alarm(anchor + 360) as expiry, bounded_persistence(anchor, expiry) as expired:
        result = {"status": "failed", "study_launch_authorization_emitted": False}
        try:
            code = p._run(args, anchor, expiry)
            if expired():
                return 1
            if code != 0 or time.monotonic() > anchor + 300:
                raise ValueError("Preparation did not finish inside owner allocation")
            result = dict(
                status="prepared_candidate_validated", evidence=validate_prepared(args, spec, anchor)
            )
            if time.monotonic() > anchor + 300:
                raise TimeoutError("Adapter validation exceeded shared allocation")
        except BaseException as exc:
            result = dict(status="failed", error=dict(type=type(exc).__name__, message=str(exc)))
        if expired():
            return 1
        s.write(directory / "terminal.json", result)
        return (
            0
            if result["status"] == "prepared_candidate_validated" and time.monotonic() <= anchor + 300
            else 1
        )


def _run(args, anchor, expiry):
    with bounded_persistence(anchor, expiry) as expired:
        if expired() or time.monotonic() >= anchor + 300:
            return 1
        spec, payload = authorize(args)
        if expired():
            return 1
        d.cpu_guard()
        if process_context(os.getpid()) != {"argv": owner_argv(args), "cwd": spec["checker"]["path"]}:
            raise ValueError("Owner must use the exact reviewed absolute-script invocation")
        owner, directory = paths(spec)
        if owner.exists() or Path(spec["output"]).exists():
            raise FileExistsError("Owner/preparation output reuse")
        if time.monotonic() >= anchor + 300 or s.workload_check(spec["workloads"]):
            raise ValueError("No current quiet preparation allocation")
        owner.mkdir()
        receipt = dict(
            schema=OWNER_SCHEMA,
            status="running",
            identity=active_identity(os.getpid()),
            started_monotonic=anchor,
            elapsed_scope="owner process startup through all preparation and finalization; not child lifetime",
            spec_sha256=args.spec_sha256,
            authorization_sha256=args.authorization_sha256,
            payload_sha256=payload,
            adapter_argv=adapter_argv(args, anchor),
            owner_argv=owner_argv(args),
            environment=environment_record(spec["checker"]["path"]),
            environment_sha256=binding(d.environment(spec["checker"]["path"])),
            success_seconds=300,
            watchdog_seconds=330,
            cleanup_seconds=30,
            parent_before=s.parent_memory(),
            study_launch_authorization_emitted=False,
        )
        first_check = True

        def check_workloads(blocked):
            nonlocal first_check
            if first_check:
                first_check = False
                if time.monotonic() >= anchor + 300:
                    raise TimeoutError("Supervisor setup exhausted owner allocation before launch")
            return s.workload_check(blocked)

        try:
            s.write(owner / "owner-receipt.json", receipt)
            if expired():
                return 1
            child = s.supervise(
                receipt["adapter_argv"],
                cwd=spec["checker"]["path"],
                env=d.environment(spec["checker"]["path"]),
                output=directory,
                allocation=300,
                watchdog=330,
                grace=30,
                blocked=spec["workloads"],
                clock=OwnerClock(anchor),
                check_workload=check_workloads,
            )
            if expired():
                return 1
            receipt["child"] = child
            if not complete(child, 300) or time.monotonic() > anchor + 300:
                raise ValueError("Adapter return/cleanup/owner allocation failed")
            terminal = d.read(directory / "terminal.json")
            evidence = validate_prepared(args, spec, anchor, check_files=False)
            if terminal != dict(status="prepared_candidate_validated", evidence=evidence) or child.get(
                "terminal_sha256"
            ) != d.digest(directory / "terminal.json"):
                raise ValueError("Adapter attestation differs from retained actual evidence")
            receipt.update(status="prepared_not_authorized", evidence=evidence)
        except BaseException as exc:
            receipt.update(status="failed", error=dict(type=type(exc).__name__, message=str(exc)))
        if expired():
            return 1
        receipt["parent_after"] = s.parent_memory()
        receipt["elapsed_before_final_write"] = time.monotonic() - anchor
        if receipt["elapsed_before_final_write"] > 300:
            receipt["status"] = "allocation_exceeded"
        try:
            s.write(owner / "owner-receipt.json", receipt)
        except BaseException:
            if expired():
                return 1
            receipt["status"] = "persistence_interrupted"
            s.write(owner / "owner-receipt.json", receipt)
            raise
        if expired():
            return 1
        if time.monotonic() > anchor + 300:
            receipt["status"] = "allocation_exceeded"
            s.write(owner / "owner-receipt.json", receipt)
        return 0 if receipt["status"] == "prepared_not_authorized" else 1


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["run", "adapter"])
    for name in ("spec", "spec-sha256", "authorization", "authorization-sha256", "output"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--anchor", type=float)
    args = parser.parse_args(argv)
    if args.command == "adapter":
        return adapter(args)
    if args.anchor is not None:
        raise ValueError("Owner derives its own process-start anchor")
    anchor = p.startup_monotonic()
    with p.administrative_alarm(anchor + 360) as expiry:
        return _run(args, anchor, expiry)


if __name__ == "__main__":
    raise SystemExit(main())
