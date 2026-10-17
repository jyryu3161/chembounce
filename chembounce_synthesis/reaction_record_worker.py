"""Internal Linux bounded native-parser worker; not a public data format/API."""

import base64
import json
import os
from pathlib import Path
import resource
import sys

# Invoked by absolute file path; bind imports to the installed package/worktree.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from chembounce_synthesis import reaction_records as records


def main():
    def deny_network(event, args):
        if event in {"socket.connect", "socket.getaddrinfo", "socket.bind"}:
            raise RuntimeError("Offline parser refuses Python network operations")
    sys.addaudithook(deny_network)
    try:
        # Establish the optional runtime before bounding additional decoder allocation.
        records._dependencies()
        from chembounce_synthesis import common  # noqa: F401
        page_size = os.sysconf("SC_PAGE_SIZE")
        current_virtual = int(Path("/proc/self/statm").read_text().split()[0]) * page_size
        cap = current_virtual + records.LIMITS["parser_extra_address_space_bytes"]
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        if hard != resource.RLIM_INFINITY:
            cap = min(cap, hard)
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
        resource.setrlimit(resource.RLIMIT_CPU, (25, 26))
        resource.setrlimit(resource.RLIMIT_FSIZE, (records.MAX_FILE, records.MAX_FILE))
        raw = sys.stdin.buffer.read(records.MAX_FILE + 1)
        records.require(len(raw) <= records.MAX_FILE, "request_bound", "Input cap")
        request = records.json_read(raw)
        snapshot = request["snapshot"]
        dataset = base64.b64decode(request["dataset"], validate=True)
        descriptor = base64.b64decode(request["descriptor"], validate=True)
        rows, raw_rows, metadata = records._replay_local(snapshot, dataset, descriptor)
        value = {"records": rows, "raw_records": [base64.b64encode(r).decode() for r in raw_rows],
                 "metadata": metadata, "parser_observation": {"baseline_virtual_bytes": current_virtual,
                 "address_space_cap_bytes": cap, "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
                 "cpu_limit_seconds": 25, "output_file_limit_bytes": records.MAX_FILE}}
        output = records.encoded(value)
        records.require(len(output) <= records.MAX_FILE, "worker_output_bound", "Output cap")
    except MemoryError:
        output = b'{"error":{"code":"parser_memory_limit","message":"Additional address-space cap reached"}}'
    except Exception as exc:
        output = json.dumps({"error": {"code": getattr(exc, "code", "parser_error"),
                                       "message": str(exc)[:2000]}}).encode()
    sys.stdout.buffer.write(output)


if __name__ == "__main__":
    main()
