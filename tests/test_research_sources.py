import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import shutil
import threading

import pytest

from benchmarks.research.sources import acquire_snapshot, archive_snapshot, main, verify_snapshot


def manifest_fixture(tmp_path, *, as_list=False, url="https://example.org/raw.json"):
    root = tmp_path / "original"
    root.mkdir()
    source = root / "raw.json"
    source.write_text('{"source":"fixed fixture"}\n')
    row = {
        "name": "test",
        "path": str(source),
        "url": url,
        "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "bytes": source.stat().st_size,
        "license": "fixture-only",
    }
    data = [row] if as_list else {"sampling": "retained provenance", "sources": [row]}
    manifest = root / "manifest.json"
    manifest.write_text(json.dumps(data))
    return manifest, source, data


@pytest.mark.parametrize("as_list", [False, True])
def test_packaged_sources_are_portable_and_keep_original_provenance(tmp_path, as_list):
    manifest, source, original = manifest_fixture(tmp_path, as_list=as_list)
    output = tmp_path / "snapshot"
    result = acquire_snapshot([manifest], output)
    portable = output / result["source_manifests"][0]["portable"]
    rewritten = json.loads(portable.read_text())
    row = rewritten[0] if as_list else rewritten["sources"][0]
    assert (portable.parent / row["path"]).read_bytes() == source.read_bytes()
    assert row["license"] == "fixture-only"
    if not as_list:
        assert rewritten["sampling"] == original["sampling"]
    retained = output / result["source_manifests"][0]["original"]
    assert retained.read_bytes() == manifest.read_bytes()
    moved = tmp_path / "relocated"
    shutil.move(output, moved)
    shutil.rmtree(manifest.parent)
    assert verify_snapshot(moved)["status"] == "verified"
    assert main(["verify", str(moved)]) == 0


def test_snapshot_detects_tampering_and_refuses_overwrite(tmp_path):
    manifest, _, _ = manifest_fixture(tmp_path)
    output = tmp_path / "snapshot"
    result = acquire_snapshot([manifest], output)
    with pytest.raises(ValueError, match="empty"):
        acquire_snapshot([manifest], output)
    raw = next(row for row in result["artifacts"] if row["role"] == "raw_source")
    (output / raw["artifact"]).write_text("changed")
    with pytest.raises(ValueError, match="hash mismatch"):
        verify_snapshot(output)


def test_package_checks_retained_hash_before_writing_snapshot(tmp_path):
    manifest, source, _ = manifest_fixture(tmp_path)
    source.write_text("changed")
    output = tmp_path / "snapshot"
    with pytest.raises(ValueError, match="Retained source hash mismatch"):
        acquire_snapshot([manifest], output)
    assert not output.exists()


def test_archive_is_deterministic_and_contains_portable_snapshot(tmp_path):
    manifest, _, _ = manifest_fixture(tmp_path)
    output = tmp_path / "snapshot"
    acquire_snapshot([manifest], output)
    first = archive_snapshot(output, tmp_path / "one.zip")
    second = archive_snapshot(output, tmp_path / "two.zip")
    assert first["sha256"] == second["sha256"]
    with pytest.raises(ValueError, match="outside"):
        archive_snapshot(output, output / "nested.zip")


@pytest.fixture
def server():
    class Handler(BaseHTTPRequestHandler):
        body = b'{"source":"fixed fixture"}\n'

        def do_GET(self):
            self.send_response(200)
            self.end_headers()
            self.wfile.write(self.body)

        def log_message(self, *args):
            pass

    service = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=service.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{service.server_port}/raw.json", Handler
    finally:
        service.shutdown()
        service.server_close()
        thread.join()


def test_fetch_verifies_bytes_and_detects_mutable_url_drift(tmp_path, server):
    url, handler = server
    manifest, source, _ = manifest_fixture(tmp_path, url=url)
    source.unlink()
    output = tmp_path / "snapshot"
    assert acquire_snapshot([manifest], output, mode="fetch")["status"] == "complete"
    handler.body = b"upstream content changed"
    failed = tmp_path / "changed"
    with pytest.raises(ValueError, match="Downloaded source hash mismatch"):
        acquire_snapshot([manifest], failed, mode="fetch")
    assert json.loads((failed / "snapshot.json").read_text())["status"] == "incomplete"
    assert not list(failed.rglob("*.part"))
    with pytest.raises(ValueError, match="incomplete"):
        verify_snapshot(failed)


def test_snapshot_verification_rejects_path_escape(tmp_path):
    manifest, _, _ = manifest_fixture(tmp_path)
    output = tmp_path / "snapshot"
    acquire_snapshot([manifest], output)
    path = output / "snapshot.json"
    data = json.loads(path.read_text())
    data["artifacts"][0]["artifact"] = "../original/raw.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="escapes"):
        verify_snapshot(output)


def test_fetch_size_limit_cannot_leave_a_complete_or_partial_artifact(tmp_path, server):
    url, _ = server
    manifest, _, _ = manifest_fixture(tmp_path, url=url)
    output = tmp_path / "oversized"
    with pytest.raises(ValueError, match="exceeds max_bytes"):
        acquire_snapshot([manifest], output, mode="fetch", max_bytes=5)
    assert json.loads((output / "snapshot.json").read_text())["status"] == "incomplete"
    assert not list(output.rglob("*.part"))
