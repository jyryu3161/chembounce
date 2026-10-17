/* Fixed offline saved-pose viewer. One shared frame; never align or mutate poses. */
"use strict";
(() => {
  function freeze(v) { if (v && typeof v === "object") { Object.values(v).forEach(freeze); Object.freeze(v); } return v; }
  document.querySelectorAll(".pose-viewer").forEach(root => {
    const data = freeze(JSON.parse(root.dataset.payload));
    const canvas = root.querySelector("canvas"), ctx = canvas.getContext("2d");
    const select = root.querySelector("select"), note = root.querySelector(".pose-note");
    const toggles = root.querySelectorAll("input[type=checkbox]");
    let yaw = 0, pitch = 0, zoom = 1, panX = 0, panY = 0, drag = null;
    const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
    const angle = x => ((x % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
    function reset() { yaw = pitch = panX = panY = 0; zoom = 1; draw(); }
    function current() { return data.seeds.find(s => String(s.seed) === select.value); }
    function draw() {
      const width = Math.max(240, Math.min(1000, canvas.clientWidth || 600)), height = 360;
      canvas.width = width; canvas.height = height;
      ctx.fillStyle = "#f7f9fc"; ctx.fillRect(0, 0, width, height);
      const seed = current();
      if (!seed || seed.pair.status !== "available") {
        note.textContent = `Seed ${select.value}: unavailable — ${seed ? seed.pair.reason : "representative seed not retained"}`;
        root.dataset.displayStatus = "unavailable";
        delete root.dataset.projected;
        return;
      }
      root.dataset.displayStatus = "available";
      note.textContent = `Saved seed ${seed.seed}. Shape ${seed.metrics.shape}; pharmacophore ${seed.metrics.pharmacophore}; polar recall ${seed.metrics.polar_recall === null ? "not applicable" : seed.metrics.polar_recall}. Reference pose ${seed.pair.reference_pose_sha256}; candidate pose ${seed.pair.candidate_pose_sha256}.`;
      const all = seed.coordinates.reference.concat(seed.coordinates.candidate);
      const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
      all.forEach(p => p.forEach((v, i) => { lo[i] = Math.min(lo[i], v); hi[i] = Math.max(hi[i], v); }));
      const origin = lo.map((v, i) => v / 2 + hi[i] / 2);
      const extent = Math.max(1e-6, ...hi.map((v, i) => v - lo[i]));
      const scale = Math.min(width, height) * .64 * zoom;
      function project(p) {
        const [x, y, z] = p.map((v, i) => (v - origin[i]) / extent);
        const a = x * Math.cos(yaw) + z * Math.sin(yaw), b = -x * Math.sin(yaw) + z * Math.cos(yaw);
        return [width / 2 + a * scale + panX, height / 2 - (y * Math.cos(pitch) - b * Math.sin(pitch)) * scale + panY, y * Math.sin(pitch) + b * Math.cos(pitch)];
      }
      const projected = {};
      ["reference", "candidate"].forEach((role, index) => {
        const graph = seed.pair[role + "_graph"], points = seed.coordinates[role].map(project);
        projected[role] = points;
        if (!toggles[index].checked) return;
        const color = index ? "#b43b68" : "#1664a3";
        ctx.strokeStyle = color; ctx.lineWidth = index ? 2 : 3; ctx.globalAlpha = .8;
        graph.bonds.forEach(b => {
          const a = points[b.begin], c = points[b.end];
          ctx.setLineDash(b.type === "AROMATIC" ? [4, 3] : []);
          ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(c[0], c[1]); ctx.stroke();
        });
        ctx.setLineDash([]); ctx.globalAlpha = 1;
        graph.atoms.forEach((a, i) => {
          const p = points[i]; ctx.fillStyle = color;
          ctx.beginPath(); ctx.arc(p[0], p[1], 4, 0, 2 * Math.PI); ctx.fill();
          if (a.element !== "C" || a.isotope || a.formal_charge) {
            const charge = a.formal_charge ? (a.formal_charge > 0 ? "+" : "") + a.formal_charge : "";
            ctx.font = "12px sans-serif"; ctx.fillText((a.isotope || "") + a.element + charge, p[0] + 6, p[1] - 5);
          }
        });
      });
      // Observable display-only coordinates for accessibility/testing; saved arrays stay frozen.
      root.dataset.projected = JSON.stringify(projected);
    }
    select.addEventListener("change", reset);
    toggles.forEach(t => t.addEventListener("change", draw));
    root.querySelector("[data-action=reset]").addEventListener("click", reset);
    root.querySelectorAll("[data-zoom]").forEach(b => b.addEventListener("click", () => { zoom = clamp(zoom * Number(b.dataset.zoom), .25, 8); draw(); }));
    canvas.addEventListener("pointerdown", e => { drag = {x:e.clientX, y:e.clientY}; canvas.setPointerCapture(e.pointerId); });
    canvas.addEventListener("pointerup", () => { drag = null; });
    canvas.addEventListener("pointercancel", () => { drag = null; });
    canvas.addEventListener("pointermove", e => {
      if (!drag) return;
      const dx = clamp(e.clientX - drag.x, -1000, 1000), dy = clamp(e.clientY - drag.y, -1000, 1000);
      if (e.shiftKey) { panX = clamp(panX + dx, -2000, 2000); panY = clamp(panY + dy, -720, 720); }
      else { yaw = angle(yaw + dx * .01); pitch = angle(pitch + dy * .01); }
      drag = {x:e.clientX, y:e.clientY}; draw();
    });
    canvas.addEventListener("wheel", e => { e.preventDefault(); zoom = clamp(zoom * Math.exp(clamp(-e.deltaY * .001, -1, 1)), .25, 8); draw(); }, {passive:false});
    canvas.addEventListener("keydown", e => {
      if (e.key === "ArrowLeft") yaw = angle(yaw - .1);
      else if (e.key === "ArrowRight") yaw = angle(yaw + .1);
      else if (e.key === "ArrowUp") pitch = angle(pitch - .1);
      else if (e.key === "ArrowDown") pitch = angle(pitch + .1);
      else if (e.key === "+" || e.key === "=") zoom = clamp(zoom * 1.2, .25, 8);
      else if (e.key === "-") zoom = clamp(zoom / 1.2, .25, 8);
      else if (e.key === "0") { reset(); e.preventDefault(); return; }
      else return;
      e.preventDefault(); draw();
    });
    new ResizeObserver(draw).observe(canvas);
    draw();
  });
})();
