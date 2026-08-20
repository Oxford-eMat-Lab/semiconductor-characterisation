/*
 * make_figures.js
 * ----------------
 * Generates CGV/figures.pptx - one slide per diagram used in
 * cgv_analysis.ipynb - using native, editable PowerPoint shapes
 * (rectangles, lines, text boxes), so the figures can be tinkered with
 * directly in PowerPoint rather than being flat images. Mirrors the
 * pattern (and helper functions) used in EQE/TLM/KPSPV's decks.
 *
 * Run:  node make_figures.js
 * Then export each slide to figures/<name>.jpg (see export_figures.sh).
 *
 * NOTE on subscripts: there is no Unicode subscript character for x, s,
 * ox, it, m or p, so every subscript here is a real PowerPoint subscript
 * run. Write them in the source as "V~FB~" and the label() helper splits
 * them.
 */
const pptxgen = require("pptxgenjs");

const pres = new pptxgen();
pres.layout = "LAYOUT_16x9"; // 10" x 5.625"

// ---- palette (kept close to the other decks) ---------------------------
const C = {
  metal: "8C8C8C",   // gate metal / contacts
  diel: "F2C14E",    // oxide
  semi: "DCE9F7",    // neutral semiconductor
  scr: "BBD3EE",     // space-charge (depletion) region
  inv: "F6C9CE",     // inversion layer
  edge: "3A4A5C",    // outlines and wiring
  text: "1F2937",
  muted: "6B7280",
  pos: "C0392B",     // positive charge
  neg: "1F6FB2",     // negative charge / electrons
  hf: "1F6FB2",       // HF curve
  lf: "0F5132",       // LF curve
  dd: "C0392B",       // deep-depletion curve
  panel: "F4F6F8",
  meas: "D95F02",     // the measured quantity
  trap: "7A5195",     // interface trap
};
const FONT = "Arial";

// ---- small helpers (shared shape across the project's decks) -----------
function runs(text, opts) {
  const base = {
    fontSize: opts.fontSize, bold: opts.bold, italic: opts.italic,
    color: opts.color, fontFace: FONT, breakLine: false,
  };
  const out = [];
  text.split("\n").forEach((lineText, li) => {
    const parts = lineText.split(/~([^~]+)~/);
    parts.forEach((part, i) => {
      if (part === "") return;
      const o = Object.assign({}, base);
      if (i % 2) o.subscript = true;
      out.push({ text: part, options: o });
    });
    if (li < text.split("\n").length - 1 && out.length) {
      out[out.length - 1].options.breakLine = true;
    }
  });
  return out;
}

function label(slide, text, x, y, w, opts = {}) {
  const o = {
    fontSize: opts.fontSize || 11,
    bold: opts.bold || false,
    italic: opts.italic || false,
    color: opts.color || C.text,
  };
  slide.addText(runs(text, o), {
    x, y, w, h: opts.h || 0.28,
    align: opts.align || "left",
    fontFace: FONT, margin: 0, valign: "middle",
  });
}
function title(slide, text) {
  slide.addText(text, {
    x: 0.4, y: 0.18, w: 9.2, h: 0.5,
    fontSize: 22, bold: true, color: C.text, fontFace: FONT, margin: 0,
  });
}
function caption(slide, text) {
  label(slide, text, 0.4, 5.02, 9.2, {
    fontSize: 11, italic: true, color: C.muted, h: 0.45,
  });
}
function box(slide, x, y, w, h, fill, opts = {}) {
  return slide.addShape(opts.shape || pres.ShapeType.rect, {
    x, y, w, h,
    fill: fill === null ? { type: "none" } : { color: fill },
    line: { color: opts.line || C.edge, width: opts.lineWidth || 1 },
  });
}
function line(slide, x, y, dx, dy, color, width = 1, opts = {}) {
  const ln = { color, width };
  if (opts.dash) ln.dashType = opts.dash;
  if (opts.arrow || opts.both) ln.endArrowType = "triangle";
  if (opts.both) ln.beginArrowType = "triangle";
  slide.addShape(pres.ShapeType.line, { x, y, w: dx, h: dy, line: ln });
}
function poly(slide, pts, color, width = 1.5, dash) {
  for (let i = 0; i < pts.length - 1; i++) {
    line(slide, pts[i][0], pts[i][1],
      pts[i + 1][0] - pts[i][0], pts[i + 1][1] - pts[i][1],
      color, width, { dash });
  }
}
function textBox(slide, text, x, y, w, h, fill, opts = {}) {
  box(slide, x, y, w, h, fill, { line: opts.line });
  slide.addText(runs(text, {
    fontSize: opts.fontSize || 10.5, bold: opts.bold || false,
    italic: false, color: opts.color || C.text,
  }), { x, y, w, h, align: "center", valign: "middle", fontFace: FONT,
        margin: 0.02 });
}
function curve(slide, pts, color, width = 2.2, dash) {
  poly(slide, pts, color, width, dash);
}

/* =====================================================================
 * SLIDE 1 - fig_cgv_mos_structure
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "The MOS capacitor: what is actually under the probes");

  const x0 = 1.55, w = 4.15;
  // gate metal
  box(s, x0, 1.15, w, 0.42, C.metal);
  label(s, "gate metal (or poly-Si)", x0 + 0.15, 1.15, w - 0.3,
    { h: 0.42, fontSize: 10.5, color: "FFFFFF", align: "center" });
  // oxide
  box(s, x0, 1.57, w, 0.30, C.diel);
  label(s, "oxide, t~ox~", x0 + w + 0.12, 1.57, 1.15,
    { h: 0.30, fontSize: 10, color: C.text });
  // semiconductor with depletion region drawn in
  box(s, x0, 1.87, w, 1.9, C.semi);
  box(s, x0, 1.87, w, 0.75, C.scr);
  label(s, "depletion region\n(fixed dopant charge)", x0 + 0.15, 1.90, w - 1.5,
    { h: 0.55, fontSize: 9.5, color: C.text });
  label(s, "W", x0 + w - 0.45, 1.90, 0.3, { h: 0.30, fontSize: 11,
    italic: true, bold: true });
  label(s, "quasi-neutral bulk (N~A~ or N~D~)", x0 + 0.15, 3.15, w - 0.3,
    { h: 0.30, fontSize: 9.5, color: C.muted });
  // back contact
  box(s, x0, 3.77, w, 0.34, C.metal);
  label(s, "back ohmic contact", x0 + 0.15, 3.77, w - 0.3,
    { h: 0.34, fontSize: 10, color: "FFFFFF", align: "center" });

  // bias source: a clean wire loop from the gate, down the left margin,
  // through the V_G source box, and back to the substrate contact.
  const busX = x0 - 0.85, vgY = 2.45, vgH = 0.44, vgW = 0.62;
  line(s, x0, 1.36, busX - x0, 0, C.edge, 1.4);        // gate lead
  line(s, busX, 1.36, 0, vgY - 1.36, C.edge, 1.4);     // bus, top half
  box(s, busX - vgW / 2, vgY, vgW, vgH, "FFFFFF", { line: C.edge });
  label(s, "V~G~", busX - vgW / 2, vgY, vgW, { h: vgH, fontSize: 12,
    bold: true, align: "center", italic: true });
  line(s, busX, vgY + vgH, 0, 3.94 - (vgY + vgH), C.edge, 1.4); // bus, bottom half
  line(s, busX, 3.94, x0 - busX, 0, C.edge, 1.4);      // contact lead

  // right-hand panel: charge picture
  const rx = 6.35, rw = 3.15;
  box(s, rx, 1.05, rw, 3.35, "FFFFFF", { line: "C7CDD4" });
  label(s, "charge balance (Eq. 3, 7)", rx + 0.15, 1.14, rw - 0.3,
    { h: 0.24, fontSize: 11, bold: true });
  label(s, "Gate charge Q~G~ = -Q~sc~", rx + 0.15, 1.46, rw - 0.3,
    { h: 0.26, fontSize: 10.5 });
  label(s, "Q~sc~ = space-charge density\nin the semiconductor (Eq. 3)",
    rx + 0.15, 1.75, rw - 0.3, { h: 0.5, fontSize: 9.5, color: C.muted });
  label(s, "For p-type in depletion:", rx + 0.15, 2.30, rw - 0.3,
    { h: 0.22, fontSize: 10 });
  label(s, "Q~sc~ = -q N~A~ W  <  0", rx + 0.15, 2.55, rw - 0.3,
    { h: 0.26, fontSize: 10.5, italic: true });
  label(s, "so the gate carries positive charge,\nbalanced by exposed " +
    "negative acceptor\nions in the depletion region.",
    rx + 0.15, 2.86, rw - 0.3, { h: 0.75, fontSize: 9.5, color: C.muted });
  label(s, "+ + +  gate", rx + 0.15, 3.66, rw - 0.3, { h: 0.24,
    fontSize: 9.5, color: C.pos });
  label(s, "-  -  -  ionised acceptors", rx + 0.15, 3.90, rw - 0.3,
    { h: 0.24, fontSize: 9.5, color: C.neg });

  caption(s, "Three terminals in, one number out per bias point: an LCR " +
    "meter reads C~m~ and G~m~ at the gate, and everything in this " +
    "notebook is inference from that single admittance.");
}

/* =====================================================================
 * SLIDE 2 - fig_cgv_equivalent_circuit
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "What the meter reads vs. what you want");

  const y = 2.2, x1 = 0.9, x2 = 3.0, x3 = 5.3, x4 = 7.6, xEnd = 8.9;
  label(s, "meter", x1 - 0.5, y - 0.6, 1.2, { h: 0.24, fontSize: 10.5,
    bold: true, color: C.muted });
  // terminal
  line(s, x1 - 0.35, y, 0.35, 0, C.edge, 1.6);
  // Rs
  box(s, x1, y - 0.22, 0.9, 0.44, "FFFFFF", { line: C.edge });
  label(s, "R~s~", x1, y - 0.22, 0.9, { h: 0.44, fontSize: 12, bold: true,
    align: "center", italic: true });
  line(s, x1 + 0.9, y, x2 - (x1 + 0.9), 0, C.edge, 1.6);
  // Cox
  box(s, x2, y - 0.22, 0.9, 0.44, C.diel);
  label(s, "C~ox~", x2, y - 0.22, 0.9, { h: 0.44, fontSize: 12, bold: true,
    align: "center", italic: true });
  line(s, x2 + 0.9, y, x3 - (x2 + 0.9), 0, C.edge, 1.6);
  // Cs parallel Gp block
  box(s, x3, y - 0.85, 1.5, 1.7, "FFFFFF", { line: "C7CDD4", shape: "roundRect" });
  label(s, "semiconductor\nbranch", x3, y - 0.80, 1.5, { h: 0.4,
    fontSize: 9, color: C.muted, align: "center" });
  box(s, x3 + 0.15, y - 0.28, 0.55, 0.4, C.semi);
  label(s, "C~s~", x3 + 0.15, y - 0.28, 0.55, { h: 0.4, fontSize: 11,
    bold: true, align: "center", italic: true });
  box(s, x3 + 0.8, y - 0.28, 0.55, 0.4, "FFF6E5", { line: C.trap });
  label(s, "G~p~", x3 + 0.8, y - 0.28, 0.55, { h: 0.4, fontSize: 11,
    bold: true, align: "center", italic: true, color: C.trap });
  label(s, "‖", x3 + 0.68, y - 0.28, 0.16, { h: 0.4, fontSize: 13,
    align: "center" });
  line(s, x3 + 1.5, y, x4 - (x3 + 1.5), 0, C.edge, 1.6);
  line(s, x4, y, xEnd - x4, 0, C.edge, 1.6, { arrow: false });
  line(s, xEnd, y, 0, -1.2, C.edge, 1.4, { dash: "dash" });
  label(s, "ground /\nback contact", xEnd - 0.3, y - 1.55, 1.3,
    { h: 0.35, fontSize: 9, color: C.muted });

  label(s, "measured: C~m~, G~m~  (Eq. 16)", x1, y + 0.55, 3.0,
    { h: 0.26, fontSize: 10.5, color: C.meas, bold: true });
  label(s, "wanted: G~p~/ω  →  D~it~  (Eq. 17)", x3, y + 0.55, 2.7,
    { h: 0.26, fontSize: 10.5, color: C.trap, bold: true });

  box(s, 0.6, 3.55, 8.8, 1.15, C.panel, { line: "C7CDD4" });
  label(s, "Two corrections stand between C~m~/G~m~ and physics:", 0.8, 3.65,
    8.4, { h: 0.24, fontSize: 11, bold: true });
  label(s, "1.  R~s~ distorts C~m~ and G~m~ even in accumulation - " +
    "extract it from the accumulation plateau (Eq. 19) and undo it (Eq. 20).",
    0.8, 3.92, 8.4, { h: 0.30, fontSize: 10 });
  label(s, "2.  C~ox~ still sits in series with the semiconductor branch - " +
    "the admittance transform (Eq. 16) removes it to leave the true " +
    "G~p~/ω.", 0.8, 4.22, 8.4, { h: 0.40, fontSize: 10 });

  caption(s, "Skip either step, or use an inaccurate C~ox~, and the error " +
    "propagates straight into every D~it~ you extract.");
}

/* =====================================================================
 * SLIDE 3 - fig_cgv_curve_shapes
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "Why HF, LF and deep-depletion C-V curves differ");

  const ox = 1.0, oy = 4.35, w = 6.3, h = 3.05;
  // axes
  line(s, ox, oy, w, 0, C.edge, 1.6, { arrow: true });
  line(s, ox, oy, 0, -h, C.edge, 1.6, { arrow: true });
  label(s, "V~G~", ox + w - 0.15, oy + 0.08, 0.5, { h: 0.22,
    fontSize: 10.5, italic: true, color: C.muted });
  label(s, "C", ox - 0.35, oy - h - 0.30, 0.3, { h: 0.24, fontSize: 12,
    italic: true, bold: true });

  const yCox = oy - h + 0.25, yCmin = oy - h + 1.55, yCdd = oy - h + 2.35;
  const xAcc = ox + 0.35, xFB = ox + 1.7, xDep = ox + 3.1, xInv = ox + 4.7,
    xEnd = ox + w - 0.25;

  // Cox reference line
  line(s, ox, yCox, w, 0, "C7CDD4", 1, { dash: "dash" });
  label(s, "C~ox~", ox - 0.55, yCox - 0.12, 0.5, { h: 0.24, fontSize: 9.5,
    color: C.muted, align: "right" });

  // HF curve: flat accumulation, drop through depletion, flat low at inversion
  curve(s, [[xAcc, yCox], [xFB, yCox], [xDep, yCmin + 0.5], [xInv, yCmin],
    [xEnd, yCmin]], C.hf, 2.4);
  // LF curve: same until inversion, then rises back to Cox
  curve(s, [[xInv - 0.25, yCmin + 0.12], [xInv + 0.25, yCmin - 0.35],
    [xEnd, yCox + 0.05]], C.lf, 2.4);
  // deep-depletion: keeps falling past Cmin
  curve(s, [[xInv - 0.25, yCmin + 0.12], [xEnd, yCdd]], C.dd, 2.2,
    "dash");

  // VFB and VT markers
  line(s, xFB, oy, 0, -h, C.muted, 1, { dash: "dash" });
  label(s, "V~FB~", xFB - 0.3, oy + 0.10, 0.6, { h: 0.22, fontSize: 10,
    italic: true, color: C.muted, align: "center" });
  line(s, xInv, oy, 0, -h, C.muted, 1, { dash: "dash" });
  label(s, "V~T~", xInv - 0.3, oy + 0.10, 0.6, { h: 0.22, fontSize: 10,
    italic: true, color: C.muted, align: "center" });

  label(s, "accumulation", xAcc, oy - h + 0.35, 1.3, { h: 0.2,
    fontSize: 9, color: C.muted });
  label(s, "depletion", xFB + 0.25, oy - h + 1.0, 1.2, { h: 0.2,
    fontSize: 9, color: C.muted });
  label(s, "inversion", xInv + 0.15, oy - h + 0.35, 1.2, { h: 0.2,
    fontSize: 9, color: C.muted });

  // legend
  const lx = ox + w + 0.35;
  [["LF (equilibrium)", C.lf, false], ["HF (1 MHz)", C.hf, false],
    ["deep depletion\n(swept too fast)", C.dd, true]].forEach(
    ([txt, col, dash], i) => {
      const ly = 1.3 + i * 0.62;
      line(s, lx, ly, 0.35, 0, col, 2.4, dash ? { dash: "dash" } : {});
      label(s, txt, lx + 0.45, ly - 0.18, 1.85, { h: 0.4, fontSize: 9.5,
        color: col });
    });

  caption(s, "All three agree in accumulation and through depletion; they " +
    "split only where minority carriers matter - whether they can be " +
    "supplied fast enough sets which curve you get (Sec. 6).");
}

/* =====================================================================
 * SLIDE 4 - fig_cgv_doping_profile
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "The C-V curve's slope is a doping profile");

  // left: cross-section with growing depletion width at 3 bias points
  const x0 = 0.65, w = 4.0;
  box(s, x0, 1.15, w, 0.35, C.metal);
  box(s, x0, 1.50, w, 0.18, C.diel);
  box(s, x0, 1.68, w, 1.75, C.semi);
  [[1.68, 0.45], [1.68, 0.95], [1.68, 1.55]].forEach(
    ([y0, dw], i) => {
      box(s, x0, y0, w, dw, C.scr, { line: i === 2 ? C.edge : "9FB4C8" });
    });
  label(s, "W(V~1~)", x0 + w + 0.1, 2.02, 1.1, { h: 0.2, fontSize: 9,
    color: C.muted });
  label(s, "W(V~2~)", x0 + w + 0.1, 2.52, 1.1, { h: 0.2, fontSize: 9,
    color: C.muted });
  label(s, "W(V~3~)", x0 + w + 0.1, 3.12, 1.1, { h: 0.2, fontSize: 9,
    color: C.muted, bold: true });

  // right: 1/C^2 vs V schematic
  const rx = 5.6, ry = 4.35, rw = 3.6, rh = 3.0;
  line(s, rx, ry, rw, 0, C.edge, 1.6, { arrow: true });
  line(s, rx, ry, 0, -rh, C.edge, 1.6, { arrow: true });
  label(s, "V", rx + rw - 0.2, ry + 0.08, 0.3, { h: 0.2, fontSize: 10.5,
    italic: true });
  label(s, "1/C²", rx - 0.5, ry - rh - 0.05, 0.6, { h: 0.24,
    fontSize: 10.5, italic: true });
  curve(s, [[rx + 0.3, ry - 0.3], [rx + rw - 0.3, ry - rh + 0.3]], C.meas,
    2.4);
  label(s, "steeper slope\n→ lower N", rx + 1.6, ry - 1.7, 1.8,
    { h: 0.5, fontSize: 9.5, color: C.muted });
  label(s, "N(W) = 2 / [qε~s~A² |d(1/C²)/dV|]  (Eq. 12/13)",
    rx - 0.1, ry + 0.35, rw + 0.6, { h: 0.24, fontSize: 9.5, italic: true });

  box(s, 0.5, 3.60, 4.35, 1.28, C.panel, { line: "C7CDD4" });
  label(s, "why |·|, not the bare derivative", 0.65, 3.67, 4.05, { h: 0.22,
    fontSize: 10, bold: true });
  label(s, "The edge of the depletion region is where N(W) is sampled " +
    "(Eq. 12/13). The textbook sign assumes V rising means more reverse " +
    "bias; this notebook's V~G~ deepens depletion directly, the opposite " +
    "sense - N(W) can never be negative, so take the absolute value " +
    "(Sec. 9).", 0.65, 3.91, 4.05, { h: 0.95, fontSize: 9.1,
    color: C.muted });

  caption(s, "Sweep V, extract W(V) from C, then N(W) from the local " +
    "slope of 1/C² - a depth profile built one derivative at a " +
    "time, well-resolved only where the slope is well-resolved.");
}

/* =====================================================================
 * SLIDE 5 - fig_cgv_conductance_method
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "The conductance peak: how a trap dissipates energy");

  // left: band diagram sketch with a trap level exchanging charge
  const bx = 0.7, by = 1.2, bw = 3.6, bh = 2.6;
  box(s, bx, by, bw, bh, "FFFFFF", { line: "C7CDD4" });
  line(s, bx + 0.4, by + 0.4, bw - 0.8, 0, C.text, 2.0);
  label(s, "E~c~", bx + 0.05, by + 0.28, 0.5, { h: 0.24, fontSize: 10,
    italic: true });
  line(s, bx + 0.4, by + bh - 0.5, bw - 0.8, 0, C.text, 2.0);
  label(s, "E~v~", bx + 0.05, by + bh - 0.62, 0.5, { h: 0.24, fontSize: 10,
    italic: true });
  line(s, bx + 0.4, by + bh / 2, bw - 0.8, 0, C.muted, 1.2, { dash: "dash" });
  label(s, "E~F~", bx + 0.05, by + bh / 2 - 0.12, 0.5, { h: 0.24,
    fontSize: 10, italic: true, color: C.muted });
  box(s, bx + bw / 2 - 0.08, by + bh / 2 - 0.55, 0.16, 0.16, C.trap);
  label(s, "E~it~ - trap level\nat the interface", bx + bw / 2 + 0.2,
    by + bh / 2 - 0.62, 1.8, { h: 0.42, fontSize: 9.5, color: C.trap });
  line(s, bx + bw / 2, by + bh - 0.5, 0, -(bh / 2 - 0.55 + 0.5 - bh / 2),
    C.trap, 1.2, { arrow: true, dash: "dash" });
  label(s, "capture / emission,\nrate set by τ~it~ (Eq. 18)",
    bx + 0.3, by + bh - 0.9, bw - 0.7, { h: 0.4, fontSize: 9,
    color: C.muted, align: "center" });

  // right: Gp/omega vs omega peak
  const rx = 5.1, ry = 4.35, rw = 4.1, rh = 3.1;
  line(s, rx, ry, rw, 0, C.edge, 1.6, { arrow: true });
  line(s, rx, ry, 0, -rh, C.edge, 1.6, { arrow: true });
  label(s, "log ω", rx + rw - 0.5, ry + 0.08, 0.6, { h: 0.2,
    fontSize: 10, italic: true });
  label(s, "G~p~/ω", rx - 0.55, ry - rh - 0.05, 0.6, { h: 0.24,
    fontSize: 10, italic: true });
  // single-level model: tall, narrow peak. The real, broadened peak is
  // both LOWER (larger y-on-page value) and WIDER - do not swap these.
  const peakX = rx + rw * 0.42, peakYsingle = ry - rh * 0.85,
    peakYbroad = ry - rh * 0.52;
  curve(s, [[rx + 0.2, ry - 0.15], [peakX - 0.7, ry - rh * 0.35],
    [peakX, peakYsingle], [peakX + 0.7, ry - rh * 0.35],
    [rx + rw - 0.3, ry - 0.15]], C.hf, 2.2, "dash");
  curve(s, [[rx + 0.2, ry - 0.10], [peakX - 1.3, ry - rh * 0.28],
    [peakX, peakYbroad], [peakX + 1.3, ry - rh * 0.28],
    [rx + rw - 0.3, ry - 0.10]], C.trap, 2.4);
  line(s, peakX, ry, 0, -rh + 0.1, C.muted, 1, { dash: "dash" });
  label(s, "ωτ~it~ = 1", peakX - 0.5, ry + 0.10, 1.0, { h: 0.22,
    fontSize: 9.5, color: C.muted, align: "center" });

  // legend sits below the axis, clear of both curves
  const lx = rx + 0.3, ly = ry + 0.32;
  line(s, lx, ly + 0.10, 0.35, 0, C.hf, 2.2, { dash: "dash" });
  label(s, "single level (Eq. 17)", lx + 0.45, ly, 2.0, { h: 0.22,
    fontSize: 9.5, color: C.hf });
  line(s, lx + 2.1, ly + 0.10, 0.35, 0, C.trap, 2.4);
  label(s, "real (broadened)", lx + 2.55, ly, 1.6, { h: 0.22,
    fontSize: 9.5, color: C.trap });

  caption(s, "A trap only shows up in G~p~ while it is actively " +
    "exchanging charge with a band - too slow or too fast relative to " +
    "ω and it looks capacitive instead, which is why the peak sits " +
    "at ωτ~it~ = 1 and not at higher or lower frequency.");
}

/* =====================================================================
 * SLIDE 6 - how this deck works (not exported)
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "How this deck works");
  const mono = (text, x, y, w, opts = {}) => {
    s.addText(text, {
      x, y, w, h: 0.24,
      fontSize: opts.fontSize || 10.5, bold: opts.bold || false,
      color: opts.color || C.text, fontFace: "Consolas", margin: 0,
      valign: "middle",
    });
  };

  box(s, 0.45, 0.95, 4.45, 2.05, C.panel, { line: "C7CDD4" });
  label(s, "1  ·  Edit, then export", 0.62, 1.06, 4.1, { h: 0.24,
    fontSize: 12, bold: true });
  label(s, "Edit any slide in this deck, then export every slide as JPEG\n" +
    "at about 150 dpi  (File → Export → JPEG → All slides).",
    0.62, 1.34, 4.15, { h: 0.55, fontSize: 10.5 });
  label(s, "Simpler: run this in the CGV folder - it renders every " +
    "slide\nand copies the JPEGs to both places for you:",
    0.62, 1.96, 4.15, { h: 0.45, fontSize: 10, italic: true,
      color: C.muted });
  mono("./export_figures.sh", 0.62, 2.46, 4.1, { bold: true,
    color: "0F5132" });

  box(s, 5.10, 0.95, 4.45, 2.05, C.panel, { line: "C7CDD4" });
  label(s, "2  ·  Save each JPEG in BOTH folders", 5.27, 1.06, 4.1,
    { h: 0.24, fontSize: 12, bold: true });
  label(s, "Same file name in each. Miss the second one and the website\n" +
    "keeps showing the old picture.", 5.27, 1.34, 4.15, { h: 0.45,
      fontSize: 10.5 });
  mono("CGV/figures/<name>.jpg", 5.27, 1.88, 4.1, { bold: true });
  label(s, "used by the Jupyter notebook", 5.27, 2.12, 4.1, { h: 0.2,
    fontSize: 9.5, color: C.muted });
  mono("docs/assets/<name>.jpg", 5.27, 2.40, 4.1, { bold: true });
  label(s, "used by the documentation website", 5.27, 2.62, 4.1, { h: 0.2,
    fontSize: 9.5, color: C.muted });

  label(s, "Slide → file name", 0.45, 3.20, 4.0, { h: 0.24,
    fontSize: 12, bold: true });
  [["Slide 1", "fig_cgv_mos_structure.jpg"],
   ["Slide 2", "fig_cgv_equivalent_circuit.jpg"],
   ["Slide 3", "fig_cgv_curve_shapes.jpg"],
   ["Slide 4", "fig_cgv_doping_profile.jpg"],
   ["Slide 5", "fig_cgv_conductance_method.jpg"]].forEach(([sl, nm], i) => {
    const y = 3.54 + i * 0.30;
    if (i % 2 === 0) box(s, 0.45, y, 9.10, 0.30, "F4F6F8",
      { line: "F4F6F8" });
    label(s, sl, 0.62, y, 1.2, { h: 0.30, fontSize: 10.5, color: C.muted });
    mono(nm, 1.85, y + 0.03, 4.5);
  });

  label(s, "This .pptx is deliberately not tracked by git - it is your " +
    "local working copy. The exported .jpg files are what the repository " +
    "and the website use, so those are the ones to commit.",
    0.45, 5.08, 9.10, { h: 0.45, fontSize: 10, italic: true,
      color: C.muted });
}

pres.writeFile({ fileName: "figures.pptx" }).then(() => {
  console.log("wrote figures.pptx");
});
