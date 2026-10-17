"""Saved-snapshot report accounting and HTML safety; no chemistry or models."""

import copy
import hashlib
from html.parser import HTMLParser
import json

import pytest

from benchmarks.research import qualified_report as report
from test_qualified_analysis import collection, plan, qualified, row


class Tables(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.rows = []
        self.tags = []
        self.cell = None
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))
        if tag == "tr":
            self.current = []
        if tag in ("th", "td"):
            self.cell = ""

    def handle_data(self, data):
        if self.cell is not None:
            self.cell += data

    def handle_endtag(self, tag):
        if tag in ("th", "td"):
            self.current.append(self.cell)
            self.cell = None
        if tag == "tr":
            self.rows.append(self.current)


def render(p, c):
    raw = json.dumps(c).encode()
    return report.render_report(
        json.dumps(p).encode(), raw, collection_sha256=hashlib.sha256(raw).hexdigest()
    )


def test_partial_hierarchy_status_costs_and_overlapping_failures():
    p = plan()  # Two targets with unequal query counts, two seeds, two methods.
    r = row(p, outputs=3)
    r["qualified"] = [qualified("same", 1), qualified("same", 2), qualified("other", 3)]
    failure = row(p, target="B", query="q3", outputs=0, status="error")
    failure["costs"]["wall_seconds"] = 7
    failure["costs"]["peak_rss_mb"] = 120
    failure["failures"] = [
        {"stage": "route", "reason": "budget"},
        {"stage": "route", "reason": "budget"},
        {"stage": "route", "reason": "unsupported"},
    ]
    no_candidates = row(p, query="q2", outputs=0)
    c = collection(p, [r, no_candidates, failure])
    before = copy.deepcopy((p, c))
    html = render(p, c)
    rows = Tables(html).rows
    # A:q1 averages 1/2; A:q2=0; B:q3=0. Equal targets gives 1/8, not pooled 1/6.
    assert ["a", "6", "1/8 (12.50%)", "1/12 (8.33%)", "0 (0.00%)", "1"] in rows
    assert ["b", "6", "0 (0.00%)", "0 (0.00%)", "0 (0.00%)", "0"] in rows
    assert "3 recorded / 12 planned; 9 missing" in html
    assert "Partial saved snapshot" in html
    assert ["a", "1", "1", "1", "0", "0", "0", "3", "0", "3"] in rows
    assert ["route", "3", "1"] in rows
    assert ["route", "budget", "2"] in rows
    assert ["a", "wall_seconds", "sum", "27", "3 / 6", "incomplete"] in rows
    assert ["a", "peak_rss_mb", "maximum", "120", "3 / 6", "incomplete"] in rows
    assert ["a", "gpu_seconds", "sum", "unknown", "0 / 6", "incomplete"] in rows
    assert ["a", "peak_gpu_mb", "maximum", "unknown", "0 / 6", "incomplete"] in rows
    assert ["a", "metrics", "hit_shape", "3", "0.7", "0.7", "0.7"] in rows
    assert ["b", "metrics", "hit_shape", "0", "unavailable", "unavailable", "unavailable"] in rows
    assert ["a", "A", "q2", "42", "recorded", "no_candidates", "0", "0", "0", "none recorded"] in rows
    assert [
        "a",
        "A",
        "q1",
        "43",
        "missing",
        "not_started",
        "0",
        "0",
        "0",
        "collection: missing_planned_attempt",
    ] in rows
    assert "not a sequential funnel" in html
    assert "not campaign duration" in html
    assert "Includes duplicate-family selected outputs" in html
    assert (p, c) == before


def test_complete_records_do_not_imply_campaign_completion_and_explicit_not_started():
    p = plan()
    records = [
        row(p, method=m, target=q["target_id"], query=q["query_id"], seed=s, outputs=0, status="not_started")
        for m in p["methods"]
        for q in p["queries"]
        for s in p["seeds"]
    ]
    c = collection(p, records)
    c["administrative_setup"] = {"status": "not_started", "wall_seconds": "UNVALIDATED_SENTINEL"}
    c["campaign_wall_seconds_through_checkpoint"] = "UNVALIDATED_SENTINEL"
    html = render(p, c)
    assert "All planned slots have records" in html
    assert "12 recorded / 12 planned; 0 missing" in html
    assert "Record completeness does not assert operational campaign completion" in html
    assert ["a", "0", "0", "0", "0", "0", "0", "6", "6", "0"] in Tables(html).rows
    assert "UNVALIDATED_SENTINEL" not in html


def test_escape_artifact_strings_and_no_external_assets():
    attack = '<script src="https://evil.invalid/x">&\'</script>'
    p = plan()
    p["run_id"] = attack
    p["methods"][0] = attack
    p["queries"][0]["query_id"] = attack
    r = row(p, method=attack, query=attack, outputs=0, status="error")
    r["failures"] = [{"stage": attack, "reason": attack}]
    r["costs"]["native_units"] = {attack: 3}
    html = render(p, collection(p, [r]))
    parsed = Tables(html)
    assert attack not in html
    assert "&lt;script" in html
    assert any(attack in cells for cells in parsed.rows)
    assert not any(tag in {"script", "img", "link", "iframe", "a"} for tag, _ in parsed.tags)
    assert not any("src" in attrs or "href" in attrs for _, attrs in parsed.tags)
    assert "default-src 'none'" in html


@pytest.mark.parametrize(
    "mutation, message",
    [
        (lambda p, c: c["attempts"].append(copy.deepcopy(c["attempts"][0])), "Duplicate planned attempt"),
        (lambda p, c: c.update(plan_sha256="f" * 64), "Stale plan binding"),
        (lambda p, c: c["attempts"][0].update(query_sha256="f" * 64), "Stale query binding"),
        (lambda p, c: c.update(collector_source_sha256="f" * 64), "Collector source mismatch"),
    ],
)
def test_reject_invalid_binding_before_render(mutation, message):
    p = plan()
    c = collection(p, [row(p)])
    mutation(p, c)
    with pytest.raises(ValueError, match=message):
        render(p, c)


def test_cli_digest_provenance_and_exclusive_output(tmp_path):
    p = plan()
    plan_path, collection_path, output = [tmp_path / n for n in ("plan.json", "collection.json", "view.html")]
    plan_path.write_text(json.dumps(p))
    collection_path.write_text(json.dumps(collection(p, [])))
    digest = hashlib.sha256(collection_path.read_bytes()).hexdigest()
    args = [
        "--plan",
        str(plan_path),
        "--collection",
        str(collection_path),
        "--output",
        str(output),
        "--collection-sha256",
    ]
    with pytest.raises(ValueError, match="Collection byte digest mismatch"):
        report.main(args + ["0" * 64])
    assert not output.exists()
    with pytest.raises(SystemExit):
        report.main(args[:-1])
    assert not output.exists()
    report.main(args + [digest])
    original = output.read_bytes()
    assert digest.encode() in original
    assert hashlib.sha256(plan_path.read_bytes()).hexdigest().encode() in original
    with pytest.raises(FileExistsError):
        report.main(args + [digest])
    assert output.read_bytes() == original
    collection_path.write_text("{}")
    with pytest.raises(ValueError, match="Collection byte digest mismatch"):
        report.main(args + [digest])
    assert output.read_bytes() == original
