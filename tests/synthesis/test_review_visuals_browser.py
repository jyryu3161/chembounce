"""Optional existing-browser checks; never download browsers or fetch resources."""

import hashlib
import json
import math
import os
from pathlib import Path

import pytest

from chembounce_synthesis import report
from .test_review_visuals import attach, real_graph_bundle, rebind
from .test_report import save

playwright = pytest.importorskip("playwright.sync_api", reason="Optional locally installed Playwright")


@pytest.fixture
def browser():
    with playwright.sync_playwright() as runtime:
        executable = os.environ.get("CHEMBOUNCE_CHROMIUM_EXECUTABLE", runtime.chromium.executable_path)
        if not Path(executable).is_file():
            pytest.skip("No existing Chromium; this test never installs one")
        instance = runtime.chromium.launch(executable_path=executable, headless=True, args=["--no-sandbox"])
        yield instance
        instance.close()


def open_local(browser, path, width=1100, height=900):
    context = browser.new_context(viewport={"width": width, "height": height}, has_touch=width < 600)
    requests, errors = [], []
    context.route("http://**/*", lambda route: (requests.append(route.request.url), route.abort()))
    context.route("https://**/*", lambda route: (requests.append(route.request.url), route.abort()))
    page = context.new_page()
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.goto(path.as_uri())
    return context, page, requests, errors


@pytest.mark.parametrize("width", [1100, 390])
def test_offline_shared_camera_seed_and_accessible_controls(browser, tmp_path, width):
    docs = attach(real_graph_bundle())
    save(tmp_path, docs)
    report.write_report(tmp_path)
    context, page, requests, errors = open_local(browser, tmp_path / "review.html", width)
    root = page.locator(".pose-viewer")
    assert root.get_attribute("data-display-status") == "available"
    before = root.get_attribute("data-payload")
    original = json.loads(root.get_attribute("data-projected"))
    assert page.locator(".pose-viewer select").input_value() == "43"
    # Both molecules share camera and normalization, so vector lengths across
    # molecules obey the exact same zoom ratio, rather than independent centering.
    page.locator("[data-zoom='1.2']").click()
    zoomed = json.loads(root.get_attribute("data-projected"))
    for role in ("reference", "candidate"):
        for axis in (0, 1):
            delta0 = original[role][1][axis] - original[role][0][axis]
            delta1 = zoomed[role][1][axis] - zoomed[role][0][axis]
            assert delta1 == pytest.approx(delta0 * 1.2)
    dx0 = original["candidate"][0][0] - original["reference"][0][0]
    dx1 = zoomed["candidate"][0][0] - zoomed["reference"][0][0]
    assert dx0 != 0 and dx1 == pytest.approx(dx0 * 1.2)
    canvas = page.locator("canvas")
    canvas.focus()
    page.keyboard.press("ArrowRight")
    rotated = json.loads(root.get_attribute("data-projected"))
    assert rotated != zoomed
    vector = lambda positions: [positions["candidate"][0][i] - positions["reference"][0][i] for i in range(3)]
    assert vector(rotated) != pytest.approx(vector(zoomed))
    scale = min(canvas.evaluate("c=>c.width"), 360) * 0.64 * 1.2
    distance = lambda positions: math.sqrt(
        sum((value / scale if i < 2 else value) ** 2 for i, value in enumerate(vector(positions)))
    )
    assert distance(rotated) == pytest.approx(distance(zoomed))
    assert root.get_attribute("data-payload") == before
    page.locator("[data-action=reset]").click()
    assert json.loads(root.get_attribute("data-projected")) == original
    if width < 600:
        canvas.scroll_into_view_if_needed()
        box = canvas.bounding_box()
        x = box["x"] + 50
        y = box["y"] + 50
        session = context.new_cdp_session(page)
        session.send("Input.dispatchTouchEvent", {"type": "touchStart", "touchPoints": [{"x": x, "y": y}]})
        session.send(
            "Input.dispatchTouchEvent", {"type": "touchMove", "touchPoints": [{"x": x + 40, "y": y + 20}]}
        )
        session.send("Input.dispatchTouchEvent", {"type": "touchEnd", "touchPoints": []})
        assert json.loads(root.get_attribute("data-projected")) != original
        page.locator("[data-action=reset]").click()
    page.locator(".pose-viewer select").select_option("42")
    assert "Saved seed 42" in page.locator(".pose-note").inner_text()
    page.locator(".pose-viewer input").first.uncheck()
    assert root.get_attribute("data-payload") == before
    page.locator(".pose-viewer input").first.check()
    page.evaluate(
        "document.querySelector('canvas').dispatchEvent(new WheelEvent('wheel',{deltaY:-1e300,cancelable:true}))"
    )
    for _ in range(100):
        page.keyboard.press("ArrowDown")
    assert "null" not in root.get_attribute("data-projected")
    assert page.evaluate("document.documentElement.scrollWidth<=window.innerWidth")
    assert not requests and not errors
    assert page.locator("img").evaluate_all("nodes=>nodes.every(i=>i.complete&&i.naturalWidth===480)")
    # Fixed CSP disallows an unapproved inline script even on a local page.
    page.evaluate(
        "() => { const s=document.createElement('script');s.textContent='window.badInjected=true';document.body.appendChild(s); }"
    )
    assert page.evaluate("window.badInjected===undefined")
    context.close()


def test_unavailable_representative_requires_explicit_choice(browser, tmp_path):
    docs = attach(real_graph_bundle())
    docs["review_visuals.json"]["candidates"][0]["seeds"][1]["pair"] = {
        "status": "unavailable",
        "reason": "fixture missing saved pose",
    }
    rebind(docs)
    save(tmp_path, docs)
    report.write_report(tmp_path)
    context, page, requests, errors = open_local(browser, tmp_path / "review.html")
    root = page.locator(".pose-viewer")
    assert root.get_attribute("data-display-status") == "unavailable"
    assert page.locator("select").input_value() == "43"
    assert "Seed 43: unavailable" in page.locator(".pose-note").inner_text()
    page.locator("select").select_option("44")
    assert root.get_attribute("data-display-status") == "available"
    page.locator("select").select_option("43")
    assert root.get_attribute("data-display-status") == "unavailable"
    assert root.get_attribute("data-projected") is None
    assert not requests and not errors
    context.close()


@pytest.mark.parametrize("coordinate", [0.0, 1e-300, 1e6])
def test_degenerate_and_extreme_bounded_camera(browser, tmp_path, coordinate):
    docs = real_graph_bundle()
    for trial in docs["refinement.json"]["candidates"][0]["trials"]:
        witness = trial["similarity_constraints"]["witness"]
        poses = [witness["reference_pose"], *trial["selective_pose_evidence"].values()]
        for pose in poses:
            pose["coordinates"] = [[coordinate] * 3 for _ in pose["coordinates"]]
            pose["sha256"] = report._binding({k: v for k, v in pose.items() if k != "sha256"})
        witness["candidate_pose_sha256"] = trial["selective_pose_evidence"]["hit_candidate"]["sha256"]
    rep = docs["refinement.json"]["candidates"][0]["trials"][1]
    docs["results.json"]["results"][0]["similarity_witness_sha256"] = report._binding(
        rep["similarity_constraints"]["witness"]
    )
    docs["scores.json"]["scores"][0].update(rep)
    attach(docs)
    save(tmp_path, docs)
    report.write_report(tmp_path)
    context, page, requests, errors = open_local(browser, tmp_path / "review.html")
    root = page.locator(".pose-viewer")
    assert root.get_attribute("data-display-status") == "available"
    before = hashlib.sha256(root.get_attribute("data-payload").encode()).hexdigest()
    for _ in range(20):
        page.locator("[data-zoom='1.2']").click()
    page.locator("canvas").focus()
    page.keyboard.press("ArrowLeft")
    assert "null" not in root.get_attribute("data-projected")
    assert hashlib.sha256(root.get_attribute("data-payload").encode()).hexdigest() == before
    assert not requests and not errors
    context.close()
