/*
 * make_figures.js
 * ----------------
 * Generates KPSPV/figures.pptx - one slide per diagram used in
 * kpspv_analysis.ipynb - using native, editable PowerPoint shapes
 * (rectangles, lines, text boxes), so the figures can be tinkered with
 * directly in PowerPoint rather than being flat images.
 *
 * Run:  node make_figures.js
 * Then export each slide to figures/<name>.jpg (see export_figures.sh).
 *
 * NOTE on subscripts: there is no Unicode subscript character for b, c, s,
 * F, p or t, so every subscript here is a real PowerPoint subscript run.
 * Write them in the source as "V~b~" and the label() helper splits them.
 */
const pptxgen = require("pptxgenjs");

const pres = new pptxgen();
pres.layout = "LAYOUT_16x9"; // 10" x 5.625"

// ---- palette -----------------------------------------------------------
const C = {
  probe: "595959",   // metal probe tip
  metal: "8C8C8C",   // stage / rear metal
  diel: "F2C14E",    // dielectric film
  semi: "DCE9F7",    // neutral semiconductor
  scr: "BBD3EE",     // space-charge region
  edge: "3A4A5C",    // outlines and wiring
  text: "1F2937",
  muted: "6B7280",
  light: "E8A33D",   // illumination / mechanical drive
  pos: "C0392B",     // positive charge
  neg: "1F6FB2",     // negative charge
  band: "1F2937",    // band edges
  fermi: "0F5132",   // Fermi level, backing voltage
  vac: "7570B3",     // vacuum level
  panel: "F4F6F8",
  spv: "D95F02",     // the measured quantity
  states: "7A5195",  // interface states
};
const FONT = "Arial";

// ---- small helpers -----------------------------------------------------
// "V~b~" -> [{text:"V"}, {text:"b", options:{subscript:true}}]
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
function charges(slide, x, y, n, dy, sign, color) {
  for (let k = 0; k < n; k++) {
    label(slide, sign, x, y + k * dy, 0.18,
      { h: 0.18, fontSize: 13, bold: true, color, align: "center" });
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

/* =====================================================================
 * SLIDE 1 - fig_kpspv_bands_cpd
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "Fermi-level alignment and the contact potential difference");

  const PW = 2.90, GAP = 0.32, PHI_P = 1.00, PHI_S = 0.50;

  const panel = (x0, tag, sub, vacOffset, connected, showCharge, nulled) => {
    box(s, x0, 1.00, PW, 3.30, "FFFFFF", { line: "C7CDD4" });
    label(s, tag, x0 + 0.12, PW - 0.24 > 0 ? 1.04 : 1.04, PW - 0.24,
      { h: 0.24, fontSize: 12, bold: true });
    label(s, sub, x0 + 0.12, 1.30, PW - 0.24,
      { h: 0.42, fontSize: 9.5, color: C.muted });
    const xp = x0 + 0.46, xs = x0 + 1.72, wcol = 0.62;
    const yvP = 1.92, yvS = yvP + vacOffset, ybody = 3.72;
    [[xp, yvP, PHI_P, "Φ~p~"], [xs, yvS, PHI_S, "Φ~s~"]].forEach(
      ([x, yv, phi, sym]) => {
        line(s, x - 0.08, yv, wcol + 0.16, 0, C.vac, 1.6, { dash: "dash" });
        line(s, x - 0.08, yv + phi, wcol + 0.16, 0, C.fermi, 2.4);
        line(s, x + wcol / 2, yv, 0, phi, C.text, 1.2, { both: true });
        label(s, sym, x + wcol / 2 + 0.07, yv + phi / 2 - 0.10, 0.45,
          { h: 0.2, fontSize: 11, italic: true });
      });
    label(s, "E~vac~", xp - 0.54, yvP - 0.13, 0.52,
      { h: 0.2, fontSize: 9.5, color: C.vac, italic: true, align: "right" });
    label(s, "E~F~", xp - 0.54, yvP + PHI_P - 0.02, 0.52,
      { h: 0.2, fontSize: 9.5, color: C.fermi, italic: true, align: "right" });
    box(s, xp, ybody, wcol, 0.28, C.probe);
    box(s, xs, ybody, wcol, 0.28, C.semi);
    label(s, "probe", xp - 0.09, ybody + 0.30, wcol + 0.18,
      { h: 0.2, fontSize: 9.5, align: "center", color: C.muted });
    label(s, "sample", xs - 0.09, ybody + 0.30, wcol + 0.18,
      { h: 0.2, fontSize: 9.5, align: "center", color: C.muted });
    if (connected) {
      const ymid = ybody + 0.14;
      line(s, xp, ymid, -0.26, 0, C.edge, 1.4);
      line(s, xp - 0.26, ymid, 0, -0.55, C.edge, 1.4);
      line(s, xs + wcol, ymid, 0.26, 0, C.edge, 1.4);
      line(s, xs + wcol + 0.26, ymid, 0, -0.55, C.edge, 1.4);
      if (nulled) {
        const yb = ymid - 0.55;
        line(s, xp - 0.26, yb, 0.62, 0, C.edge, 1.4);
        line(s, xs + wcol + 0.26, yb, -0.62, 0, C.edge, 1.4);
        [[0.0, 0.20], [0.07, 0.10]].forEach(([dy, w]) => {
          line(s, (xp + xs + wcol) / 2 - w / 2, yb + dy, w, 0, C.fermi, 2.2);
        });
        label(s, "V~b~", (xp + xs + wcol) / 2 - 0.45, yb - 0.26, 0.9,
          { h: 0.2, fontSize: 11, bold: true, color: C.fermi,
            align: "center" });
      } else {
        line(s, xp - 0.26, ymid - 0.55, xs + wcol + 0.52 - xp, 0, C.edge,
          1.4);
      }
    }
    if (showCharge) {
      charges(s, xp + wcol + 0.04, ybody - 0.02, 3, 0.11, "−", C.neg);
      charges(s, xs - 0.22, ybody - 0.02, 3, 0.11, "+", C.pos);
    }
    if (vacOffset) {
      const xa = xs - 0.22;
      line(s, xp + wcol + 0.08, yvP, xa - xp - wcol - 0.08, 0, C.spv, 0.9,
        { dash: "sysDot" });
      line(s, xa, yvP, 0, vacOffset, C.spv, 2.0, { both: true });
      label(s, "qV~CPD~", xa - 0.80, yvP + vacOffset / 2 - 0.10, 0.72,
        { h: 0.22, fontSize: 10, bold: true, color: C.spv, align: "right" });
    }
  };

  panel(0.22, "A  ·  isolated",
    "Vacuum levels line up.\nThe Fermi levels do not.", 0.0, false, false,
    false);
  panel(0.22 + PW + GAP, "B  ·  connected",
    "Electrons move until the Fermi\nlevels match. Now E~vac~ does not.",
    PHI_P - PHI_S, true, true, false);
  panel(0.22 + 2 * (PW + GAP), "C  ·  nulled",
    "A backing voltage removes the\nfield, and is what is read.", 0.0, true,
    false, true);

  label(s, "CPD  =  (Φ~s~ − Φ~p~) / q   =   − V~b~ (null)",
    0.22, 4.42, 9.56,
    { h: 0.34, fontSize: 15, bold: true, align: "center", color: C.spv });
  caption(s, "Neither work function changes on connection — only the " +
    "offset between the two vacuum levels, and that offset is the CPD.");
}

/* =====================================================================
 * SLIDE 2 - fig_kpspv_instrument
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "The vibrating capacitor, and where the measurement is made");

  // ---- left: the Kelvin circuit ----
  box(s, 0.30, 0.92, 4.35, 3.85, "FFFFFF", { line: "C7CDD4" });
  label(s, "i(t) = (V~b~ + CPD) · dC/dt", 0.42, 1.08, 4.1,
    { h: 0.32, fontSize: 13.5, bold: true, align: "center", color: C.spv });
  label(s, "a DC potential, read as an AC current", 0.42, 1.40, 4.1,
    { h: 0.24, fontSize: 10, italic: true, color: C.muted,
      align: "center" });
  const xt = 1.55, wt = 1.55, yt = 2.30;
  box(s, xt, yt, wt, 0.24, C.probe);
  label(s, "probe tip", xt + wt + 0.14, yt - 0.02, 1.2,
    { h: 0.24, fontSize: 10.5, bold: true });
  line(s, xt + wt / 2, 1.86, 0, 0.36, C.light, 2.2, { both: true });
  label(s, "d~1~ sin ωt", xt + wt / 2 + 0.14, 1.88, 1.1,
    { h: 0.22, fontSize: 10, italic: true, color: C.light });
  line(s, xt - 0.20, yt + 0.24, 0, 0.52, C.muted, 1.0, { both: true });
  label(s, "d~0~", xt - 0.60, yt + 0.42, 0.36,
    { h: 0.22, fontSize: 11, italic: true, color: C.muted,
      align: "right" });
  label(s, "C(t)", xt + wt / 2 + 0.34, yt + 0.34, 0.6,
    { h: 0.22, fontSize: 11, italic: true, color: C.muted });
  box(s, 1.30, 3.06, 2.05, 0.28, C.semi);
  box(s, 1.30, 3.34, 2.05, 0.14, C.metal);
  label(s, "sample", 3.48, 3.06, 1.1, { h: 0.24, fontSize: 10.5,
    bold: true });
  label(s, "stage", 3.48, 3.30, 1.1, { h: 0.22, fontSize: 9.5,
    color: C.muted });
  line(s, xt, yt + 0.12, -0.93, 0, C.edge, 1.4);
  line(s, 0.62, yt + 0.12, 0, 0.62, C.edge, 1.4);
  [[0.0, 0.30], [0.09, 0.15]].forEach(([dy, w]) => {
    line(s, 0.62 - w / 2, yt + 0.74 + dy, w, 0, C.fermi, 2.4);
  });
  label(s, "V~b~", 0.14, yt + 0.72, 0.40,
    { h: 0.24, fontSize: 11.5, bold: true, color: C.fermi,
      align: "right" });
  line(s, 0.62, yt + 0.83, 0, 1.30, C.edge, 1.4);
  line(s, 0.62, 4.13, 1.28, 0, C.edge, 1.4);
  box(s, 1.90, 3.90, 0.62, 0.46, "FFFFFF",
    { shape: pres.ShapeType.triangle });
  label(s, "preamp", 1.72, 4.40, 1.0,
    { h: 0.22, fontSize: 9.5, color: C.muted, align: "center" });
  line(s, 2.52, 4.13, 0.83, 0, C.edge, 1.4);
  line(s, 3.35, 4.13, 0, -0.65, C.edge, 1.4);

  // ---- right: the off-null line ----
  box(s, 4.95, 0.92, 4.75, 3.85, "FFFFFF", { line: "C7CDD4" });
  label(s, "V~ptp~ = k (V~b~ + CPD)", 5.05, 1.08, 4.55,
    { h: 0.32, fontSize: 13.5, bold: true, align: "center", color: C.spv });
  label(s, "measure away from the null, then extrapolate to it", 5.05, 1.40,
    4.55, { h: 0.24, fontSize: 10, italic: true, color: C.muted,
      align: "center" });
  const OX = 5.42, OW = 3.90, YAX = 3.20, YTOP = 1.95, YBOT = 4.45;
  const x0 = OX + 0.14, x1 = OX + OW - 0.14, y0 = 4.28, y1 = 2.12;
  const xnull = x0 + ((y0 - YAX) / (y0 - y1)) * (x1 - x0);
  [OX + 0.02, OX + OW - 0.86].forEach((xg) => {   // drawn first: behind
    box(s, xg, YTOP + 0.02, 0.84, YBOT - YTOP - 0.04, "E9F4EC",
      { line: "CBE4D3" });
  });
  box(s, OX, YAX - 0.16, OW, 0.32, "FBE4E4", { line: "FBE4E4" });
  line(s, OX, YAX, OW, 0, C.edge, 1.2);
  line(s, xnull, YTOP, 0, YBOT - YTOP, C.muted, 1.1, { dash: "dash" });
  poly(s, [[x0, y0], [x1, y1]], C.neg, 2.6);
  [0.03, 0.09, 0.15, 0.79, 0.85, 0.91].forEach((t) => {
    box(s, x0 + t * (x1 - x0) - 0.045, y0 + t * (y1 - y0) - 0.045, 0.09,
      0.09, C.pos, { line: C.pos, shape: pres.ShapeType.ellipse });
  });
  label(s, "off-null\nwindow", OX + 0.04, YTOP + 0.10, 0.80,
    { h: 0.40, fontSize: 9, align: "center", bold: true, color: "0F5132" });
  label(s, "off-null\nwindow", OX + OW - 0.84, YBOT - 0.52, 0.80,
    { h: 0.40, fontSize: 9, align: "center", bold: true, color: "0F5132" });
  label(s, "noise floor", OX + 1.05, YAX - 0.01, 0.9,
    { h: 0.2, fontSize: 9, color: C.pos });
  label(s, "null:  V~b~ = −CPD", xnull - 0.82, YBOT + 0.06, 1.64,
    { h: 0.24, fontSize: 10.5, bold: true, color: C.muted,
      align: "center" });
  label(s, "V~b~", OX + OW - 0.02, YAX + 0.08, 0.34,
    { h: 0.22, fontSize: 11, italic: true });
  label(s, "V~ptp~", OX - 0.44, YTOP - 0.02, 0.42,
    { h: 0.22, fontSize: 11, italic: true, align: "right" });
  caption(s, "The gradient k drifts with spacing, area and gain; the zero " +
    "crossing does not. So the crossing is the measurement — and it is " +
    "also where the signal-to-noise ratio is worst.");
}

/* =====================================================================
 * SLIDE 3 - fig_kpspv_surface_bands
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "A charged dielectric bends the bands of the semiconductor");

  const Y_EC = 2.34, Y_EI = 2.78, Y_EV = 3.92, Y_EF = 2.48;

  // Page y grows downward, so smaller y is higher energy. bend > 0
  // therefore draws the bands DOWN in energy at the surface - E_c falling
  // towards E_F, i.e. accumulation on n-type. bend < 0 draws them up in
  // energy, towards inversion.
  const bands = (x0, sign, tag, sub, bend) => {
    box(s, x0, 1.00, 4.62, 3.55, "FFFFFF", { line: "C7CDD4" });
    label(s, tag, x0 + 0.14, 1.04, 4.34, { h: 0.24, fontSize: 12,
      bold: true });
    label(s, sub, x0 + 0.14, 1.28, 4.34, { h: 0.24, fontSize: 9.5,
      color: C.muted });
    const xd = x0 + 0.34, xi = x0 + 0.86, xb = x0 + 1.86, xr = x0 + 4.16;
    box(s, xd, 1.60, 0.52, 2.62, C.diel);
    box(s, xi, 1.60, 1.00, 2.62, C.scr, { line: C.scr });
    label(s, "dielectric", xd - 0.16, 4.28, 0.84,
      { h: 0.2, fontSize: 9, color: C.muted, align: "center" });
    label(s, "space charge", xi, 4.28, 1.00,
      { h: 0.2, fontSize: 9, color: C.muted, align: "center" });
    label(s, "neutral silicon (n-type)", xb - 0.14, 4.28, 1.74,
      { h: 0.2, fontSize: 9, color: C.muted, align: "center" });
    [[Y_EC, C.band, "E~c~", 2.0, undefined],
     [Y_EI, C.muted, "E~i~", 1.0, "dash"],
     [Y_EV, C.band, "E~v~", 2.0, undefined]].forEach(
      ([yb2, col, name, lw, dash]) => {
        poly(s, [[xr, yb2], [xb, yb2], [xi, yb2 + bend]], col, lw, dash);
        label(s, name, xr + 0.06, yb2 - 0.10, 0.34,
          { h: 0.2, fontSize: 10, italic: true, color: col });
      });
    line(s, xi, Y_EF, xr - xi, 0, C.fermi, 2.0, { dash: "sysDash" });
    label(s, "E~F~", xr + 0.06, Y_EF - 0.10, 0.34,
      { h: 0.2, fontSize: 10, italic: true, color: C.fermi });
    const pos = sign === "+";
    charges(s, xd + 0.24, 1.74, 4, 0.19, sign, pos ? C.pos : C.neg);
    label(s, "Q~f~", xd - 0.02, 1.42, 0.56,
      { h: 0.2, fontSize: 11, bold: true, color: pos ? C.pos : C.neg,
        align: "center" });
    charges(s, xi + 0.62, 1.74, 4, 0.19, pos ? "−" : "+",
      pos ? C.neg : C.pos);
    label(s, "Q~sc~", xi + 0.42, 1.42, 0.56,
      { h: 0.2, fontSize: 10, bold: true, color: pos ? C.neg : C.pos,
        align: "center" });
    const top = Y_EC + bend + 0.20, bot = Y_EV + bend - 0.20;
    for (let k = 0; k < 5; k++) {
      line(s, xi - 0.07, top + (k * (bot - top)) / 4, 0.14, 0, C.states,
        1.8);
    }
    label(s, "D~it~", xi - 0.18, (top + bot) / 2 - 0.10, 0.34,
      { h: 0.2, fontSize: 10, bold: true, color: C.states,
        align: "right" });
    line(s, xi, Y_EC, xb - xi, 0, C.spv, 0.9, { dash: "sysDot" });
    line(s, xi + 0.76, Y_EC, 0, bend, C.spv, 2.2, { both: true });
    label(s, "qφ~s~", xi + 0.84, Y_EC + bend * 0.5 - 0.10, 0.5,
      { h: 0.2, fontSize: 11, bold: true, color: C.spv });
  };

  bands(0.22, "+", "Positive film charge  →  accumulation",
    "electrons pile up at the surface;  φ~s~ > 0", 0.32);
  bands(5.16, "−", "Negative film charge  →  inversion",
    "holes pile up at the surface;  φ~s~ < 0", -0.60);
  caption(s, "Charge in the film is mirrored by charge in the silicon and " +
    "in the interface states. Q~f~ + Q~sc~ + Q~it~ = 0 is what fixes " +
    "φ~s~.");
}

/* =====================================================================
 * SLIDE 4 - fig_kpspv_spv_bands
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "Light flattens the bands — and the change in CPD is the SPV");

  const DY_EC = 2.40, DY_EI = 3.18, DY_EV = 3.96, DY_EF = 2.72;

  const spvPanel = (x0, tag, sub, bend, split) => {
    box(s, x0, 1.00, 4.62, 3.50, "FFFFFF", { line: "C7CDD4" });
    label(s, tag, x0 + 0.14, 1.04, 4.34, { h: 0.24, fontSize: 12,
      bold: true });
    label(s, sub, x0 + 0.14, 1.28, 4.34, { h: 0.24, fontSize: 9.5,
      color: C.muted });
    const xd = x0 + 0.34, xi = x0 + 0.82, xb = x0 + 1.92, xr = x0 + 4.16;
    box(s, xd, 1.62, 0.48, 2.58, C.diel);
    box(s, xi, 1.62, 1.10, 2.58, C.scr, { line: C.scr });
    charges(s, xd + 0.20, 1.76, 4, 0.19, "−", C.neg);
    [[DY_EC, C.band, "E~c~", 2.0, undefined],
     [DY_EI, C.muted, "E~i~", 1.0, "dash"],
     [DY_EV, C.band, "E~v~", 2.0, undefined]].forEach(
      ([yb2, col, name, lw, dash]) => {
        poly(s, [[xr, yb2], [xb, yb2], [xi, yb2 + bend]], col, lw, dash);
        label(s, name, xr + 0.06, yb2 - 0.10, 0.34,
          { h: 0.2, fontSize: 10, italic: true, color: col });
      });
    if (split) {
      // on n-type at moderate injection E_Fn barely moves; E_Fp does
      line(s, xi, DY_EF - 0.06, xr - xi, 0, C.neg, 1.8,
        { dash: "sysDash" });
      line(s, xi, DY_EF + split, xr - xi, 0, C.pos, 1.8,
        { dash: "sysDash" });
      label(s, "E~Fn~", xr + 0.06, DY_EF - 0.16, 0.42,
        { h: 0.2, fontSize: 10, italic: true, color: C.neg });
      label(s, "E~Fp~", xr + 0.06, DY_EF + split - 0.10, 0.42,
        { h: 0.2, fontSize: 10, italic: true, color: C.pos });
      [0.42, 0.90, 1.38].forEach((dx) => {
        line(s, xb + dx, DY_EV - 0.06, 0, DY_EC - DY_EV + 0.12, C.light,
          2.0, { arrow: true });
      });
      label(s, "hν", xb + 0.28, DY_EC - 0.30, 0.4,
        { h: 0.2, fontSize: 10, bold: true, italic: true, color: C.light });
    } else {
      line(s, xi, DY_EF, xr - xi, 0, C.fermi, 2.0, { dash: "sysDash" });
      label(s, "E~F~", xr + 0.06, DY_EF - 0.10, 0.34,
        { h: 0.2, fontSize: 10, italic: true, color: C.fermi });
    }
    line(s, xi, DY_EC, xb - xi, 0, C.spv, 0.9, { dash: "sysDot" });
    line(s, xi + 0.80, DY_EC, 0, bend, C.spv, 2.2, { both: true });
    label(s, "qφ~s~", xi + 0.88, DY_EC + bend * 0.5 - 0.10, 0.6,
      { h: 0.2, fontSize: 11, bold: true, color: C.spv });
    label(s, "dielectric", xd - 0.16, 4.24, 0.80,
      { h: 0.2, fontSize: 9, color: C.muted, align: "center" });
    label(s, "n-type silicon", xb - 0.20, 4.24, 1.50,
      { h: 0.2, fontSize: 9, color: C.muted, align: "center" });
  };

  spvPanel(0.22, "Dark", "the bands are bent by the film charge", -0.62,
    null);
  spvPanel(5.16, "Illuminated",
    "injected carriers screen the charge; the bands flatten", -0.14, 0.78);
  line(s, 4.90, 2.62, 0.30, 0, C.light, 2.6, { arrow: true });
  label(s, "SPV  =  CPD~dark~ − CPD~light~  =  φ~s~(light) − φ~s~(dark)",
    0.22, 4.58, 9.56,
    { h: 0.32, fontSize: 15, bold: true, align: "center", color: C.spv });
  caption(s, "The probe work function, the calibration and the dielectric " +
    "term are identical in the two readings, so they cancel. Only the " +
    "band bending survives.");
}

/* =====================================================================
 * SLIDE 5 - fig_kpspv_setup
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "A macro-scale Kelvin probe / SPV bench");

  box(s, 0.55, 1.05, 5.85, 3.55, "FBFCFD", { line: C.muted,
    lineWidth: 1.6 });
  label(s, "earthed enclosure  (electrical screening + darkness)", 0.62,
    4.28, 5.7, { h: 0.24, fontSize: 10, italic: true, color: C.muted,
      align: "center" });

  box(s, 2.55, 1.35, 1.75, 0.52, "E8EDF2");
  label(s, "vibrating head", 2.55, 1.42, 1.75, { h: 0.36, fontSize: 10,
    bold: true, align: "center" });
  line(s, 3.42, 1.87, 0, 0.44, C.edge, 1.6);
  box(s, 2.90, 2.31, 1.05, 0.20, C.probe);
  line(s, 4.12, 2.05, 0, 0.34, C.light, 2.2, { both: true });
  label(s, "d~0~ ≈ 0.1–1 mm", 4.36, 2.42, 1.5, { h: 0.22, fontSize: 9.5,
    italic: true, color: C.muted });
  box(s, 1.55, 2.75, 2.70, 0.30, C.semi);
  box(s, 1.55, 3.05, 3.75, 0.22, C.metal);
  label(s, "sample", 1.55, 2.78, 2.70, { h: 0.24, fontSize: 10.5,
    bold: true, align: "center" });
  box(s, 4.42, 2.75, 0.88, 0.30, "F5D97B");
  label(s, "ref.", 4.42, 2.78, 0.88, { h: 0.24, fontSize: 9.5, bold: true,
    align: "center" });
  label(s, "reference sample (Au / HOPG),\nmeasured before and after",
    3.55, 3.35, 2.6, { h: 0.40, fontSize: 9.5, italic: true,
      color: C.muted, align: "center" });
  label(s, "earthed x–y stage", 1.55, 3.05, 3.75, { h: 0.22, fontSize: 9,
    color: "FFFFFF", align: "center" });

  line(s, 0.95, 1.95, 1.55, 0.72, C.light, 3.0, { arrow: true });
  box(s, 0.62, 1.55, 0.75, 0.42, "FDF0DC", { line: C.light });
  label(s, "fibre", 0.62, 1.58, 0.75, { h: 0.36, fontSize: 9.5,
    bold: true, align: "center", color: C.light });

  box(s, 6.75, 1.05, 2.85, 1.35, C.panel, { line: "C7CDD4" });
  label(s, "control electronics", 6.88, 1.14, 2.6, { h: 0.24,
    fontSize: 11.5, bold: true });
  label(s, "• vibration drive, ω\n• backing voltage V~b~\n" +
    "• preamp + off-null fit\n• spacing regulation from k",
    6.88, 1.40, 2.6, { h: 0.92, fontSize: 10 });
  box(s, 6.75, 2.60, 2.85, 1.00, C.panel, { line: "C7CDD4" });
  label(s, "light source", 6.88, 2.69, 2.6, { h: 0.24, fontSize: 11.5,
    bold: true });
  label(s, "• halogen or LED, shuttered\n• intensity swept for saturation\n" +
    "• wavelength swept for L~n~", 6.88, 2.95, 2.6, { h: 0.62,
      fontSize: 10 });
  line(s, 6.72, 1.72, -0.30, 0, C.edge, 1.4, { arrow: true });
  line(s, 6.72, 3.05, -0.30, 0, C.light, 1.4, { arrow: true });
  box(s, 6.75, 3.80, 2.85, 0.80, "FFF6E5", { line: C.light });
  label(s, "outputs", 6.88, 3.86, 2.6, { h: 0.22, fontSize: 10.5,
    bold: true });
  label(s, "CPD~dark~ ,  CPD~light~  →  SPV", 6.88, 4.10, 2.6,
    { h: 0.30, fontSize: 11.5, bold: true, color: C.spv });
  caption(s, "Nothing touches the sample. The probe averages over its own " +
    "diameter, typically 2 mm, so what it reports is an area average.");
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
  label(s, "Simpler: run this in the KPSPV folder — it renders every " +
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
  mono("KPSPV/figures/<name>.jpg", 5.27, 1.88, 4.1, { bold: true });
  label(s, "used by the Jupyter notebook", 5.27, 2.12, 4.1, { h: 0.2,
    fontSize: 9.5, color: C.muted });
  mono("docs/assets/<name>.jpg", 5.27, 2.40, 4.1, { bold: true });
  label(s, "used by the documentation website", 5.27, 2.62, 4.1, { h: 0.2,
    fontSize: 9.5, color: C.muted });

  label(s, "Slide → file name", 0.45, 3.20, 4.0, { h: 0.24, fontSize: 12,
    bold: true });
  [["Slide 1", "fig_kpspv_bands_cpd.jpg"],
   ["Slide 2", "fig_kpspv_instrument.jpg"],
   ["Slide 3", "fig_kpspv_surface_bands.jpg"],
   ["Slide 4", "fig_kpspv_spv_bands.jpg"],
   ["Slide 5", "fig_kpspv_setup.jpg"]].forEach(([sl, nm], i) => {
    const y = 3.54 + i * 0.30;
    if (i % 2 === 0) box(s, 0.45, y, 9.10, 0.30, "F4F6F8",
      { line: "F4F6F8" });
    label(s, sl, 0.62, y, 1.2, { h: 0.30, fontSize: 10.5, color: C.muted });
    mono(nm, 1.85, y + 0.03, 4.5);
  });

  label(s, "This .pptx is deliberately not tracked by git — it is your " +
    "local working copy. The exported .jpg files are what the repository " +
    "and the website use, so those are the ones to commit.",
    0.45, 5.08, 9.10, { h: 0.45, fontSize: 10, italic: true,
      color: C.muted });
}

pres.writeFile({ fileName: "figures.pptx" }).then(() => {
  console.log("wrote figures.pptx");
});
