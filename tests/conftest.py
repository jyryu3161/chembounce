import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

ADP = "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N"


# Deliberately requested only by synthetic benchmark orchestration/collector tests.
# These temporary artifacts exercise real profile validation, not production match.
import pytest


@pytest.fixture
def synthetic_historical_profile_bindings(tmp_path, monkeypatch):
    import copy
    import hashlib
    import json

    from chembounce_synthesis import profiles

    original_resource = profiles._resource
    original_digests = dict(profiles.PROFILE_SHA256)
    paths, artifacts, packaged = {}, {}, {}
    with monkeypatch.context() as scoped:
        for name in (profiles.PROFILE, profiles.SELECTIVE_PROFILE):
            raw = original_resource(name).read_bytes()
            packaged[name] = hashlib.sha256(raw).hexdigest()
            artifact = copy.deepcopy(profiles.load_profile(name))
            current = profiles.numerical_runtime(name)
            artifact["provenance"]["runtime"] = {
                key: value for key, value in current.items() if key != "source_sha256"
            }
            artifact["provenance"]["source_sha256"] = current["source_sha256"]
            artifact["provenance"]["synthetic_fixture"] = {
                "purpose": "SYNTHETIC orchestration/collector coverage; NOT calibration or production binding",
                "derived_from_packaged_sha256": packaged[name],
            }
            artifact["limitations"] = [
                "SYNTHETIC TEST ARTIFACT; no empirical fit, source release or scientific performance claim",
                *artifact.get("limitations", []),
            ]
            if name == profiles.PROFILE:
                artifact["fit_sha256"] = "a" * 64
                artifact["provenance"]["procedure_sha256"] = "b" * 64
                artifact["review"] = {key: "c" * 64 for key in artifact["review"]}
                artifact["calibration_spec"]["version"] = "SYNTHETIC TEST ONLY; NOT A FIT"
            path = tmp_path / (name + "-SYNTHETIC.json")
            path.write_text(json.dumps(artifact, sort_keys=True))
            paths[name], artifacts[name] = path, artifact
            scoped.setitem(profiles.PROFILE_SHA256, name, hashlib.sha256(path.read_bytes()).hexdigest())
        scoped.setattr(profiles, "_resource", lambda name: paths[name])
        yield {"paths": paths, "artifacts": artifacts, "packaged_sha256": packaged}
    assert profiles._resource is original_resource
    assert profiles.PROFILE_SHA256 == original_digests
    for name, digest in packaged.items():
        assert hashlib.sha256(original_resource(name).read_bytes()).hexdigest() == digest
