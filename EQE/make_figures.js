/*
 * make_figures.js
 * ----------------
 * Generates EQE/figures.pptx - one slide per diagram used in
 * eqe_analysis.ipynb - using native, editable PowerPoint shapes
 * (rectangles, lines, text boxes), so the figures can be tinkered with
 * directly in PowerPoint rather than being flat images.
 *
 * Run:  node make_figures.js
 * Then export each slide to figures/<name>.jpg (see export_figures.sh).
 */
const pptxgen = require("pptxgenjs");

const pres = new pptxgen();
pres.layout = "LAYOUT_16x9"; // 10" x 5.625"

// ---- palette -----------------------------------------------------------
const C = {
  arc: "5B9BD5", // anti-reflection coating
  emitter: "A8C4E0", // n+ emitter
  base: "DCE9F7", // p-type base
  bsf: "D3BFA6", // p+ back surface field
  diel: "F2C14E", // dielectric rear passivation
  metal: "595959", // metallisation
  edge: "3A4A5C", // outlines
  text: "1F2937",
  muted: "6B7280",
  light: "E8A33D", // incident light
  refl: "D95F02", // reflection loss
  para: "7570B3", // parasitic absorption
  trans: "1B9E77", // transmitted / collected
  loss: "C0392B", // recombination loss
  panel: "F4F6F8",
};
const FONT = "Arial";

// ---- small helpers -----------------------------------------------------
function title(slide, text) {
  slide.addText(text, {
    x: 0.4, y: 0.22, w: 9.2, h: 0.5,
    fontSize: 22, bold: true, color: C.text, fontFace: FONT, margin: 0,
  });
}
function label(slide, text, x, y, w, opts = {}) {
  slide.addText(text, {
    x, y, w, h: opts.h || 0.28,
    fontSize: opts.fontSize || 11,
    bold: opts.bold || false,
    italic: opts.italic || false,
    color: opts.color || C.text,
    align: opts.align || "left",
    fontFace: FONT, margin: 0, valign: "middle",
  });
}
function box(slide, x, y, w, h, fill, opts = {}) {
  slide.addShape(pres.ShapeType.rect, {
    x, y, w, h,
    fill: { color: fill },
    line: { color: opts.line || C.edge, width: opts.lineWidth || 1 },
  });
}
function arrow(slide, x, y, w, h, color, width = 2) {
  slide.addShape(pres.ShapeType.line, {
    x, y, w, h,
    line: { color, width, endArrowType: "triangle" },
  });
}
function plainLine(slide, x, y, w, h, color, width = 1, dash) {
  const line = { color, width };
  if (dash) line.dashType = dash;
  slide.addShape(pres.ShapeType.line, { x, y, w, h, line });
}
function textBox(slide, text, x, y, w, h, fill, opts = {}) {
  slide.addShape(pres.ShapeType.rect, {
    x, y, w, h,
    fill: { color: fill },
    line: { color: opts.line || C.edge, width: 1 },
  });
  slide.addText(text, {
    x, y, w, h,
    fontSize: opts.fontSize || 10.5,
    bold: opts.bold || false,
    color: opts.color || C.text,
    align: "center", valign: "middle", fontFace: FONT, margin: 0.02,
  });
}

/* =====================================================================
 * SLIDE 1 - fig_cell_structure
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "Crystalline silicon solar cell — cross-section");

  const X = 2.0, W = 5.0;
  const yARC = 1.32, hARC = 0.16;
  const yEm = yARC + hARC, hEm = 0.30;
  const yBase = yEm + hEm, hBase = 1.95;
  const yBSF = yBase + hBase, hBSF = 0.22;
  const yRear = yBSF + hBSF, hRear = 0.32;

  // incident light (arrows fall between the contact fingers)
  [2.2, 3.5, 5.2, 6.6].forEach((ax) => {
    arrow(s, ax, 0.86, 0, 0.40, C.light, 2.25);
  });
  label(s, "incident light", 0.45, 0.90, 1.6, { fontSize: 11, color: C.light, bold: true });

  // layer stack
  box(s, X, yARC, W, hARC, C.arc);
  box(s, X, yEm, W, hEm, C.emitter);
  box(s, X, yBase, W, hBase, C.base);
  box(s, X, yBSF, W, hBSF, C.bsf);
  box(s, X, yRear, W, hRear, C.metal);

  // front contact fingers
  [0.45, 2.1, 3.75].forEach((dx) => {
    box(s, X + dx, yARC - 0.26, 0.55, 0.26, C.metal);
  });

  // right-hand labels
  const L = X + W + 0.18;
  label(s, "front contact fingers (Ag)", L, yARC - 0.15, 3.0, { fontSize: 10.5 });
  label(s, "anti-reflection coating", L, yARC + 0.06, 3.0, { fontSize: 10.5 });
  label(s, "n⁺ emitter", L, yEm + 0.10, 3.0, { fontSize: 10.5, bold: true });
  label(s, "p-type base  (absorber)", L, yBase + hBase / 2 - 0.14, 3.0, { fontSize: 10.5, bold: true });
  label(s, "p⁺ back surface field", L, yBSF + 0.00, 3.0, { fontSize: 10.5 });
  label(s, "rear metallisation (Al)", L, yRear + 0.05, 3.0, { fontSize: 10.5 });

  // depth axis
  arrow(s, X - 0.42, yEm, 0, hEm + hBase, C.text, 1.5);
  label(s, "z", X - 0.68, yEm + 0.8, 0.3, { fontSize: 13, italic: true, bold: true });
  label(s, "0", X - 0.68, yEm - 0.06, 0.3, { fontSize: 10, color: C.muted });

  // thickness bracket
  plainLine(s, X - 0.20, yBase, 0.10, 0, C.muted, 1);
  plainLine(s, X - 0.20, yBSF, 0.10, 0, C.muted, 1);
  plainLine(s, X - 0.15, yBase, 0, hBase, C.muted, 1);
  label(s, "W", X - 0.42, yBase + hBase / 2 - 0.12, 0.3, { fontSize: 12, italic: true, color: C.muted });

  label(s,
    "Wavelength selects depth: blue light is absorbed within the emitter and coating, " +
    "red light reaches the rear.",
    0.4, 4.80, 9.2, { fontSize: 11, italic: true, color: C.muted });
}

/* =====================================================================
 * SLIDE 2 - fig_optical_losses  (EQE vs IQE)
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "Where the incident photons go — EQE vs IQE");

  // front surface line
  const ySurf = 1.95;
  plainLine(s, 1.0, ySurf, 6.2, 0, C.edge, 1.5);
  label(s, "front surface", 7.30, ySurf - 0.16, 1.4, { fontSize: 10, color: C.muted });

  // incident
  arrow(s, 2.1, 1.05, 0.55, 0.80, C.light, 3);
  label(s, "incident photons  Φ₀", 0.42, 1.14, 1.7, { fontSize: 11.5, bold: true, color: C.light });

  // reflected
  arrow(s, 2.68, 1.88, 0.62, -0.78, C.refl, 3);
  label(s, "reflected  R" + "ₑₓₜ", 3.35, 1.02, 1.6, { fontSize: 11.5, bold: true, color: C.refl });

  // parasitic absorption (front layers)
  box(s, 1.0, ySurf, 6.2, 0.34, "E5E1F2", { line: C.para });
  label(s, "coating + emitter:  parasitic absorption  A" + "ₑₓₜ",
    1.12, ySurf + 0.03, 5.0, { fontSize: 10.5, color: C.para, bold: true });

  // transmitted into absorber
  arrow(s, 2.35, ySurf + 0.36, 0, 0.42, C.trans, 3);
  label(s, "T" + "ₑₓₜ" + " = 1 − R" + "ₑₓₜ" + " − A" + "ₑₓₜ",
    2.55, ySurf + 0.40, 2.2, { fontSize: 11.5, bold: true, color: C.trans });

  // absorber
  box(s, 1.0, 2.80, 6.2, 1.05, C.base);
  label(s, "absorber:  photons absorbed → electron-hole pairs generated",
    1.15, 2.90, 5.6, { fontSize: 11, bold: true });

  // split collected / recombined
  arrow(s, 2.4, 3.60, -0.55, 0.42, C.trans, 2.5);
  arrow(s, 5.6, 3.60, 0.55, 0.42, C.loss, 2.5);
  textBox(s, "collected\n→ photocurrent", 1.05, 4.05, 1.85, 0.62, "DFF3EC", { color: "0F5132", bold: true, fontSize: 10.5 });
  textBox(s, "recombined\n→ lost", 5.6, 4.05, 1.55, 0.62, "FBE4E2", { color: "8A2A22", bold: true, fontSize: 10.5 });

  // definitions panel
  box(s, 7.35, 2.80, 2.25, 1.87, C.panel, { line: "C7CDD4" });
  label(s, "EQE  =  collected / incident", 7.48, 2.92, 2.05, { fontSize: 10.5, bold: true });
  label(s, "IQE  =  collected / entered", 7.48, 3.20, 2.05, { fontSize: 10.5, bold: true });
  plainLine(s, 7.48, 3.50, 1.98, 0, "C7CDD4", 1);
  label(s, "EQE = T" + "ₑₓₜ" + " × IQE        (10)", 7.48, 3.60, 2.05, { fontSize: 10, color: C.muted });
  label(s, "IQE = EQE / T" + "ₑₓₜ" + "         (11)", 7.48, 3.90, 2.05, { fontSize: 10, color: C.muted });
  label(s, "≈ EQE / (1 − R" + "ₑₓₜ" + ") only", 7.60, 4.18, 2.05, { fontSize: 9.5, color: C.muted, italic: true });
  label(s, "where A" + "ₑₓₜ" + " ≈ 0", 7.60, 4.40, 2.05, { fontSize: 9.5, color: C.muted, italic: true });
  
  label(s,
    "EQE counts every incident photon, so it is charged for reflection. " +
    "IQE divides that out and reports only how well the cell converts the light it actually receives.",
    0.4, 4.98, 9.2, { fontSize: 11, italic: true, color: C.muted, h: 0.5 });
}

/* =====================================================================
 * SLIDE 3 - fig_collection_efficiency
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "Generation depth and carrier collection");

  const X = 1.25, Y = 1.35, W = 7.0, H = 1.85;

  box(s, X, Y, W, H, C.base);
  box(s, X, Y, 0.16, H, C.emitter); // front / junction
  box(s, X + W - 0.16, Y, 0.16, H, C.metal); // rear

  label(s, "front / junction", X - 0.05, Y - 0.32, 1.7, { fontSize: 10, color: C.muted });
  label(s, "rear surface", X + W - 1.35, Y - 0.32, 1.5, { fontSize: 10, color: C.muted, align: "right" });

  // penetration arrows for three wavelengths
  const lambdas = [
    { nm: "400 nm", frac: 0.10, dy: 0.42, col: "3B5BA5" },
    { nm: "700 nm", frac: 0.45, dy: 0.92, col: "2E7D32" },
    { nm: "1000 nm", frac: 0.93, dy: 1.42, col: "C0392B" },
  ];
  lambdas.forEach((L) => {
    arrow(s, X + 0.16, Y + L.dy, (W - 0.32) * L.frac, 0, L.col, 2.5);
    label(s, L.nm, X + 0.24, Y + L.dy - 0.36, 1.0, { fontSize: 9.5, bold: true, color: L.col });
  });

  // depth axis
  arrow(s, X, Y + H + 0.16, W, 0, C.text, 1.25);
  label(s, "depth  z", X + W / 2 - 0.5, Y + H + 0.22, 1.2, { fontSize: 11, italic: true });
  label(s, "0", X - 0.04, Y + H + 0.22, 0.3, { fontSize: 10, color: C.muted });
  label(s, "W", X + W - 0.22, Y + H + 0.22, 0.3, { fontSize: 10, color: C.muted });

  // collection efficiency wedge (schematic: high at front, decaying to rear)
  const yEta = Y + H + 0.72;
  s.addShape(pres.ShapeType.rtTriangle, {
    x: X, y: yEta, w: W, h: 0.78,
    fill: { color: "DFF3EC" }, line: { color: "1B9E77", width: 1 },
  });
  s.addText(
    [
      { text: "collection efficiency  \u03B7" },
      { text: "c", options: { subscript: true } },
    ],
    { x: X + 0.12, y: yEta - 0.30, w: 2.6, h: 0.28, fontSize: 11, bold: true,
      color: "0F5132", fontFace: FONT, margin: 0, valign: "middle" });
  label(s, "high", X + 0.12, yEta + 0.42, 0.8, { fontSize: 10, color: "0F5132" });
  label(s, "low", X + W - 0.75, yEta + 0.42, 0.6, { fontSize: 10, color: "0F5132", align: "right" });

  label(s,
    "Blue light generates carriers where collection is efficient but front losses are largest; " +
    "red light generates them deep, where the diffusion length and rear surface decide the outcome.",
    0.4, 5.05, 9.2, { fontSize: 11, italic: true, color: C.muted });
}

/* =====================================================================
 * SLIDE 4 - fig_albsf_vs_perc
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "Al-BSF vs PERC — the rear side sets the near-IR response");

  function cell(x, variant) {
    const W = 3.7;
    const yARC = 1.35, hARC = 0.14;
    const yEm = yARC + hARC, hEm = 0.24;
    const yBase = yEm + hEm, hBase = 1.35;
    const yRearLayer = yBase + hBase, hRearLayer = 0.22;
    const yMetal = yRearLayer + hRearLayer, hMetal = 0.28;

    box(s, x, yARC, W, hARC, C.arc);
    box(s, x, yEm, W, hEm, C.emitter);
    box(s, x, yBase, W, hBase, C.base);

    if (variant === "albsf") {
      box(s, x, yRearLayer, W, hRearLayer, C.bsf);
      label(s, "p⁺ Al-BSF", x + W + 0.06, yRearLayer - 0.02, 0.95, { fontSize: 9 });
    } else {
      // dielectric with local contact openings
      box(s, x, yRearLayer, W, hRearLayer, C.diel);
      [0.75, 1.85, 2.95].forEach((dx) => {
        box(s, x + dx, yRearLayer, 0.28, hRearLayer, C.metal, { lineWidth: 0.75 });
      });
      label(s, "dielectric +\nlocal contacts", x + W + 0.05, yRearLayer - 0.10, 0.88, { fontSize: 9, h: 0.5 });
    }
    box(s, x, yMetal, W, hMetal, C.metal);

    return { yBase, hBase, yRearLayer, yMetal, hMetal, W };
  }

  const xL = 0.55, xR = 5.35;
  const g = cell(xL, "albsf");
  cell(xR, "perc");

  label(s, "Al-BSF", xL, 0.92, 3.7, { fontSize: 14, bold: true, align: "center" });
  label(s, "PERC", xR, 0.92, 3.7, { fontSize: 14, bold: true, align: "center" });

  // near-IR light behaviour: absorbed/lost at rear vs reflected back
  // Al-BSF: the ray reaches the rear and is absorbed / lost there
  arrow(s, xL + 0.35, g.yBase + 0.12, 2.55, g.yRearLayer - 0.03 - (g.yBase + 0.12),
        C.loss, 2.25);
  label(s, "lost at rear", xL + 1.15, g.yRearLayer - 0.30, 1.6,
        { fontSize: 9.5, color: C.loss, bold: true });

  // PERC: the ray reaches the rear and is specularly reflected back up,
  // drawn as a true V with its vertex on the rear reflector
  const vX = xR + 2.10, vY = g.yRearLayer - 0.03;
  arrow(s, xR + 0.35, g.yBase + 0.12, vX - (xR + 0.35), vY - (g.yBase + 0.12),
        "1B9E77", 2.25);
  s.addShape(pres.ShapeType.line, {
    x: vX, y: g.yBase + 0.22, w: 1.25, h: vY - (g.yBase + 0.22),
    flipV: true,
    line: { color: "1B9E77", width: 2.25, endArrowType: "triangle" },
  });
  label(s, "reflected back", xR + 0.95, vY - 0.30, 1.9,
        { fontSize: 9.5, color: "1B9E77", bold: true });

  // captions
  textBox(s,
    "high rear recombination (large S)\nweak rear reflection",
    xL, 4.05, 3.7, 0.62, "FBE4E2", { color: "8A2A22", fontSize: 10.5 });
  textBox(s,
    "passivated rear (small S)\nrear reflector returns near-IR light",
    xR, 4.05, 3.7, 0.62, "DFF3EC", { color: "0F5132", fontSize: 10.5 });

  label(s,
    "Both cells share the same front optics, so their EQE and IQE curves differ only above ~900 nm — " +
    "the signature of a rear-side difference.",
    0.4, 4.95, 9.2, { fontSize: 11, italic: true, color: C.muted, h: 0.45 });
}

/* =====================================================================
 * SLIDE 5 - fig_measurement_setup
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "Differential spectral responsivity (DSR) measurement");

  const yRow = 0.95, bh = 0.55, bw = 1.30;

  // --- optical chain, left to right ---
  textBox(s, "Xe / halogen\nlamp", 0.45, yRow, bw, bh, "FDF2DC", { fontSize: 9.5 });
  arrow(s, 1.78, yRow + bh / 2, 0.25, 0, C.text, 2);
  textBox(s, "mono-\nchromator", 2.06, yRow, bw, bh, "FDF2DC", { fontSize: 9.5 });
  arrow(s, 3.39, yRow + bh / 2, 0.25, 0, C.text, 2);
  textBox(s, "chopper", 3.67, yRow, bw, bh, "FDF2DC", { fontSize: 9.5 });
  arrow(s, 5.00, yRow + bh / 2, 0.25, 0, C.text, 2);
  textBox(s, "beam\nsplitter", 5.28, yRow, bw, bh, "FDF2DC", { fontSize: 9.5 });
  arrow(s, 6.61, yRow + bh / 2, 0.25, 0, C.muted, 1.75);
  textBox(s, "monitor\nphotodiode", 6.89, yRow, 1.55, bh, C.panel,
    { fontSize: 9.5, color: C.muted });

  // --- monochromatic beam down onto the cell ---
  arrow(s, 5.93, yRow + bh, 0, 1.05, C.light, 2.5);
  label(s, "monochromatic,\nchopped", 6.05, yRow + bh + 0.16, 1.7,
    { fontSize: 9.5, color: C.light, h: 0.42 });

  // --- bias lamp ---
  textBox(s, "bias lamp\n(white light)", 0.60, 1.95, 1.45, bh, "FFF6D9", { fontSize: 9.5 });
  arrow(s, 2.08, 2.42, 1.55, 0.28, C.light, 2.25);

  // --- the cell ---
  const yCell = 2.62;
  box(s, 3.30, yCell, 3.30, 0.28, C.base);
  box(s, 3.30, yCell + 0.28, 3.30, 0.18, C.metal);
  label(s, "solar cell on temperature-controlled chuck", 3.30, yCell + 0.52, 3.30,
    { fontSize: 9.5, align: "center", color: C.muted });

  // --- signal chain ---
  arrow(s, 6.63, yCell + 0.14, 0.40, 0, C.text, 2);
  textBox(s, "trans-impedance\namplifier", 7.05, 2.45, 1.60, bh, "E9EEF5", { fontSize: 9.5 });
  arrow(s, 7.85, 3.00, 0, 0.33, C.text, 2);
  textBox(s, "lock-in amplifier", 7.05, 3.35, 1.60, 0.48, "E9EEF5", { fontSize: 9.5 });

  // signal out of the lock-in, arrowhead on the left end
  s.addShape(pres.ShapeType.line, {
    x: 6.42, y: 3.59, w: 0.60, h: 0,
    line: { color: C.text, width: 2, beginArrowType: "triangle" },
  });
  s.addText(
    [
      { text: "\u0394j" },
      { text: "sc", options: { subscript: true } },
      { text: "  (AC signal)" },
    ],
    { x: 4.75, y: 3.45, w: 1.60, h: 0.28, fontSize: 10.5, bold: true,
      align: "right", color: C.text, fontFace: FONT, margin: 0, valign: "middle" });

  // --- reference cell note ---
  textBox(s, "reference cell of known SR  →  sets the absolute scale",
    0.45, 3.95, 4.30, 0.45, C.panel, { fontSize: 10, color: C.muted });

  // --- key equation, with proper subscripts ---
  box(s, 0.45, 4.58, 4.30, 0.55, "EAF3EF", { line: "1B9E77" });
  s.addText(
    [
      { text: "s\u0303 (\u03BB, E" },
      { text: "bias", options: { subscript: true } },
      { text: ")  =  \u0394j" },
      { text: "sc", options: { subscript: true } },
      { text: " / \u0394E" },
      { text: "\u03BB", options: { subscript: true } },
      { text: "        (15)" },
    ],
    { x: 0.60, y: 4.58, w: 4.05, h: 0.55, fontSize: 12.5, bold: true,
      color: "0F5132", fontFace: FONT, margin: 0, valign: "middle" });

  label(s,
    "The bias lamp sets a realistic injection level; the lock-in recovers only the small chopped signal.",
    5.00, 4.62, 4.60, { fontSize: 10.5, italic: true, color: C.muted, h: 0.48 });
}

/* =====================================================================
 * SLIDE 6 - how to update these figures
 * This slide documents the workflow; it is NOT exported as a figure.
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "How to update these figures");

  const mono = (txt, x, y, w, opts = {}) => s.addText(txt, {
    x, y, w, h: opts.h || 0.26, fontSize: opts.fontSize || 10.5,
    bold: opts.bold || false, color: opts.color || C.text,
    fontFace: "Courier New", margin: 0, valign: "middle",
  });

  // --- step 1 -----------------------------------------------------------
  box(s, 0.45, 0.95, 4.45, 2.05, C.panel, { line: "C7CDD4" });
  label(s, "1  ·  Edit, then export", 0.62, 1.06, 4.1, { fontSize: 12, bold: true });
  label(s, "Edit any slide in this deck, then export every slide as JPEG\nat about 150 dpi  (File \u2192 Export \u2192 JPEG \u2192 All slides).",
    0.62, 1.36, 4.15, { fontSize: 10.5, h: 0.55 });
  label(s, "Simpler: run this in the EQE folder \u2014 it renders every slide\nand copies the JPEGs to both places for you:",
    0.62, 1.96, 4.15, { fontSize: 10, italic: true, color: C.muted, h: 0.45 });
  mono("./export_figures.sh", 0.62, 2.46, 4.1, { bold: true, color: "0F5132" });

  // --- step 2 -----------------------------------------------------------
  box(s, 5.10, 0.95, 4.45, 2.05, C.panel, { line: "C7CDD4" });
  label(s, "2  ·  Save each JPEG in BOTH folders", 5.27, 1.06, 4.1, { fontSize: 12, bold: true });
  label(s, "Same file name in each. Miss the second one and the website\nkeeps showing the old picture.",
    5.27, 1.36, 4.15, { fontSize: 10.5, h: 0.45 });
  mono("EQE/figures/<name>.jpg", 5.27, 1.88, 4.1, { bold: true });
  label(s, "used by the Jupyter notebook", 5.27, 2.12, 4.1, { fontSize: 9.5, color: C.muted });
  mono("docs/assets/<name>.jpg", 5.27, 2.40, 4.1, { bold: true });
  label(s, "used by the documentation website", 5.27, 2.62, 4.1, { fontSize: 9.5, color: C.muted });

  // --- slide -> file name table ----------------------------------------
  label(s, "Slide \u2192 file name", 0.45, 3.20, 4.0, { fontSize: 12, bold: true });
  const row = (i, slide, name) => {
    const y = 3.54 + i * 0.30;
    if (i % 2 === 0) box(s, 0.45, y, 9.10, 0.30, "F4F6F8", { line: "F4F6F8" });
    label(s, slide, 0.62, y, 1.2, { fontSize: 10.5, color: C.muted });
    mono(name, 1.85, y, 4.5, { fontSize: 10.5 });
  };
  row(0, "Slide 1", "fig_cell_structure.jpg");
  row(1, "Slide 2", "fig_optical_losses.jpg");
  row(2, "Slide 3", "fig_collection_efficiency.jpg");
  row(3, "Slide 4", "fig_albsf_vs_perc.jpg");
  row(4, "Slide 5", "fig_measurement_setup.jpg");

  label(s,
    "This .pptx is deliberately not tracked by git \u2014 it is your local working copy. " +
    "The exported .jpg files are what the repository and the website use, so those are the ones to commit.",
    0.45, 5.08, 9.10, { fontSize: 10, italic: true, color: C.muted, h: 0.45 });
}


pres.writeFile({ fileName: "figures.pptx" }).then(() => {
  console.log("wrote figures.pptx");
});
