/*
 * make_figures.js
 * ----------------
 * Generates TLM/figures.pptx - one slide per diagram used in
 * tlm_analysis.ipynb - using native, editable PowerPoint shapes
 * (rectangles, lines, text boxes), so the figures can be tinkered with
 * directly in PowerPoint rather than being flat images.
 *
 * Run:  node make_figures.js
 * Then export each slide to figures/<name>.jpg (see export_figures.sh).
 *
 * Note on notation: Unicode has no subscript "c", "s", "T" or "λ", so
 * subscripted symbols are built with pptxgenjs rich-text arrays
 * ([{text:"R"},{text:"C",options:{subscript:true}}]), never with Unicode
 * subscript characters - those silently render as the wrong letter.
 */
const pptxgen = require("pptxgenjs");

const pres = new pptxgen();
pres.layout = "LAYOUT_16x9"; // 10" x 5.625"

// ---- palette -----------------------------------------------------------
const C = {
  metal: "595959",      // metal contact
  metalLite: "8C8C8C",
  semi: "DCE9F7",       // doped semiconductor sheet
  semiDeep: "A8C4E0",   // the sheet seen edge-on
  sub: "F0EDE6",        // insulating substrate / isolating layer
  edge: "3A4A5C",       // outlines
  text: "1F2937",
  muted: "6B7280",
  current: "1B9E77",    // current flow
  crowd: "D95F02",      // current crowding / transfer length
  fit: "7570B3",        // fitted line
  data: "C0392B",       // data points
  panel: "F4F6F8",
  grid: "D5DBE1",
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
function arrow(slide, x, y, w, h, color, width = 2, opts = {}) {
  slide.addShape(pres.ShapeType.line, {
    x, y, w, h,
    flipV: opts.flipV || false,
    line: {
      color, width,
      endArrowType: "triangle",
      beginArrowType: opts.double ? "triangle" : "none",
      dashType: opts.dash,
    },
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
// dimension line with ticks at both ends and a caption above it
function dimension(slide, x1, x2, y, caption, opts = {}) {
  const color = opts.color || C.muted;
  plainLine(slide, x1, y - 0.09, 0, 0.18, color, 1);
  plainLine(slide, x2, y - 0.09, 0, 0.18, color, 1);
  arrow(slide, x1, y, x2 - x1, 0, color, 1, { double: true });
  slide.addText(caption, {
    x: x1, y: y - 0.36, w: x2 - x1, h: 0.24,
    fontSize: opts.fontSize || 10, color, align: "center",
    fontFace: FONT, margin: 0, valign: "middle",
  });
}
// subscripted symbol as a rich-text run array
function sym(base, subscript, rest) {
  const runs = [{ text: base, options: { italic: true } },
                { text: subscript, options: { subscript: true, italic: true } }];
  if (rest) runs.push({ text: rest });
  return runs;
}

/* =====================================================================
 * SLIDE 1 - fig_tlm_resistance_chain
 * Where the measured resistance comes from: the series chain of Eq. (1).
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "What a two-probe measurement actually measures");

  const yTop = 1.80, hSemi = 0.80;
  const xL = 1.5, xR = 8.5;
  const wPad = 1.3;

  // substrate + semiconductor sheet
  box(s, xL, yTop + hSemi, xR - xL, 0.45, C.sub, { line: C.grid });
  label(s, "insulating substrate", xL + 0.15, yTop + hSemi + 0.09, 2.4,
        { fontSize: 10, color: C.muted, italic: true });
  box(s, xL, yTop, xR - xL, hSemi, C.semi);
  label(s, "doped semiconductor layer", xL + 2.2, yTop + hSemi - 0.32, 2.8,
        { fontSize: 11, color: C.text, align: "center" });

  // two metal contacts sitting on top
  box(s, xL, yTop - 0.28, wPad, 0.28, C.metal, { line: C.edge });
  box(s, xR - wPad, yTop - 0.28, wPad, 0.28, C.metal, { line: C.edge });
  label(s, "metal contact", xL - 0.05, yTop - 0.66, 1.4,
        { fontSize: 10, color: C.text, align: "center" });
  label(s, "metal contact", xR - wPad - 0.05, yTop - 0.66, 1.4,
        { fontSize: 10, color: C.text, align: "center" });

  // current path through the chain: down through one contact, along the
  // sheet, and back up through the other
  const yPath = yTop + 0.32;
  arrow(s, xL + 0.62, yTop - 0.28, 0, 0.60, C.current, 2.25);
  plainLine(s, xL + 0.62, yPath, xR - xL - 1.24, 0, C.current, 2.25);
  arrow(s, xR - 0.62, yTop - 0.28, 0, 0.60, C.current, 2.25, { flipV: true });
  label(s, "current path", 0.20, yPath - 0.14, 1.15,
        { fontSize: 10, color: C.current, bold: true, align: "right" });

  // the series chain of resistances
  const yChain = 3.35, hChain = 0.5;
  const seg = [
    { w: 0.9, fill: C.metalLite, t: sym("R", "m") },
    { w: 1.15, fill: C.crowd, t: sym("R", "C") },
    { w: 2.9, fill: C.semiDeep, t: sym("R", "semi") },
    { w: 1.15, fill: C.crowd, t: sym("R", "C") },
    { w: 0.9, fill: C.metalLite, t: sym("R", "m") },
  ];
  let x = xL;
  seg.forEach((g) => {
    box(s, x, yChain, g.w, hChain, g.fill);
    s.addText(g.t, {
      x, y: yChain, w: g.w, h: hChain, fontSize: 13, bold: true,
      color: C.text, align: "center", valign: "middle", fontFace: FONT,
      margin: 0,
    });
    x += g.w;
  });

  dimension(s, xL, x, yChain + 0.92, "");
  s.addText([{ text: "measured  " }, ...sym("R", "T"),
             { text: "  =  2" }, ...sym("R", "m"), { text: "  +  2" },
             ...sym("R", "C"), { text: "  +  " }, ...sym("R", "semi")], {
    x: xL, y: yChain + 0.98, w: x - xL, h: 0.34, fontSize: 14,
    color: C.text, align: "center", valign: "middle", fontFace: FONT,
    margin: 0,
  });

  label(s,
    "The metal term is negligible; the two contact terms are what TLM is built to isolate.",
    xL, yChain + 1.42, x - xL,
    { fontSize: 11, color: C.muted, align: "center", italic: true });
}

/* =====================================================================
 * SLIDE 2 - fig_tlm_structure
 * The test pattern: one strip, many contacts, several spacings.
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "The TLM test structure (top view)");

  const yStrip = 1.60, hStrip = 1.45;
  const xL = 1.35, xR = 9.30;

  // caption above the strip
  s.addText([{ text: "isolated doped layer, sheet resistance " },
             { text: "R", options: { italic: true } },
             { text: "S", options: { subscript: true, italic: true } }], {
    x: xL, y: yStrip - 0.40, w: 4.2, h: 0.28, fontSize: 11, color: C.muted,
    italic: true, fontFace: FONT, margin: 0, valign: "middle",
  });

  // the isolated doped strip
  box(s, xL, yStrip, xR - xL, hStrip, C.semi);

  // contacts with progressively larger gaps
  const wPad = 0.55;
  const gaps = [0.25, 0.45, 0.72, 1.05, 1.50];
  const xs = [];
  let x = xL + 0.25;
  for (let i = 0; i <= gaps.length; i++) {
    xs.push(x);
    box(s, x, yStrip + 0.11, wPad, hStrip - 0.22, C.metal, { line: C.edge });
    if (i < gaps.length) x += wPad + gaps[i];
  }

  // contact width W, marked outside the strip on the left
  const xW = xL - 0.20;
  plainLine(s, xW - 0.07, yStrip + 0.11, 0.14, 0, C.muted, 1);
  plainLine(s, xW - 0.07, yStrip + hStrip - 0.11, 0.14, 0, C.muted, 1);
  arrow(s, xW, yStrip + 0.11, 0, hStrip - 0.22, C.muted, 1, { double: true });
  s.addText([{ text: "W", options: { italic: true } }], {
    x: 0.30, y: yStrip + hStrip / 2 - 0.26, w: 0.72, h: 0.28, fontSize: 14,
    color: C.muted, align: "right", fontFace: FONT, margin: 0,
    valign: "middle",
  });
  label(s, "contact width", 0.05, yStrip + hStrip / 2 + 0.02, 0.97,
        { fontSize: 8.5, color: C.muted, align: "right" });

  // contact length d, marked on the last contact
  const xd = xs[xs.length - 1];
  dimension(s, xd, xd + wPad, yStrip - 0.24, "", { fontSize: 9 });
  s.addText([{ text: "d", options: { italic: true } }], {
    x: xd - 0.25, y: yStrip - 0.64, w: wPad + 0.5, h: 0.26, fontSize: 13,
    color: C.muted, align: "center", fontFace: FONT, margin: 0,
    valign: "middle",
  });

  // spacings L1..L5 below
  const yDim = yStrip + hStrip + 0.60;
  for (let i = 0; i < gaps.length; i++) {
    const a = xs[i] + wPad, b = xs[i + 1];
    dimension(s, a, b, yDim, "", { fontSize: 9 });
    s.addText([{ text: "L", options: { italic: true } },
               { text: String(i + 1), options: { subscript: true } }], {
      x: (a + b) / 2 - 0.35, y: yDim + 0.10, w: 0.70, h: 0.26,
      fontSize: 11, color: C.muted, align: "center", fontFace: FONT,
      margin: 0, valign: "middle",
    });
  }

  label(s,
    "Each neighbouring pair gives one resistance at one spacing. "
    + "Same sheet, same contacts — only L changes.",
    0.6, yDim + 0.70, 8.8,
    { fontSize: 12, color: C.text, align: "center" });
}

/* =====================================================================
 * SLIDE 3 - fig_tlm_crowding
 * Why the contact area is not the drawn area.
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "Current crowding and the transfer length");

  // ---- left panel: cross-section under one contact ----
  const yTop = 2.15, hSemi = 0.85;
  const xL = 0.55, xR = 5.05;
  const xPadL = 2.45;

  label(s, "leading edge of the contact", xL, yTop - 0.98, 1.85,
        { fontSize: 9.5, color: C.crowd, align: "center" });
  plainLine(s, xL + 0.92, yTop - 0.70, xPadL - (xL + 0.92), 0.32, C.crowd, 1,
            "dash");

  box(s, xPadL, yTop - 0.28, xR - xPadL - 0.10, 0.28, C.metal, { line: C.edge });
  label(s, "metal contact", xPadL + 0.35, yTop - 0.62, 1.9,
        { fontSize: 10, color: C.text, align: "center" });

  box(s, xL, yTop, xR - xL, hSemi, C.semi);
  label(s, "doped layer", xL + 0.10, yTop + hSemi - 0.30, 1.4,
        { fontSize: 10, color: C.muted, italic: true });

  // lateral current arriving from the left
  plainLine(s, xL + 0.18, yTop + 0.44, xPadL - xL - 0.55, 0, C.current, 2.25);
  arrow(s, xPadL - 0.45, yTop + 0.44, 0.42, 0, C.current, 2.25);

  // vertical injection arrows, shortening away from the leading edge
  [0.08, 0.30, 0.56, 0.90, 1.32, 1.82].forEach((dx) => {
    const hArr = 0.42 * Math.exp(-dx / 0.55) + 0.05;
    arrow(s, xPadL + dx, yTop, 0, hArr, C.crowd, 2.0, { flipV: true });
  });

  // the transfer length marked under the contact
  dimension(s, xPadL, xPadL + 0.55, yTop + hSemi + 0.38, "");
  s.addText([{ text: "L", options: { italic: true } },
             { text: "T", options: { subscript: true, italic: true } }], {
    x: xPadL - 0.25, y: yTop + hSemi + 0.46, w: 1.05, h: 0.28,
    fontSize: 13, color: C.crowd, align: "center", bold: true,
    fontFace: FONT, margin: 0, valign: "middle",
  });
  label(s, "most of the current has crossed within one transfer length",
        xL - 0.10, yTop + hSemi + 0.86, (xR - xL) + 0.2,
        { fontSize: 10, color: C.muted, align: "center", italic: true });

  // ---- right panel: the exponential decay ----
  const px = 6.35, py = 1.80, pw = 3.05, ph = 1.75;
  s.addText([{ text: "current still flowing in the semiconductor, " },
             { text: "I", options: { italic: true } },
             { text: "(x)", options: { italic: true } }], {
    x: px - 0.6, y: py - 0.38, w: pw + 1.2, h: 0.28, fontSize: 10,
    color: C.muted, align: "center", fontFace: FONT, margin: 0,
    valign: "middle",
  });

  box(s, px, py, pw, ph, C.panel, { line: C.grid });
  plainLine(s, px, py + ph, pw, 0, C.edge, 1);
  plainLine(s, px, py, 0, ph, C.edge, 1);

  // exponential decay drawn as a polyline of short segments
  const N = 26, LTpx = pw * 0.22;
  for (let i = 0; i < N; i++) {
    const x1 = (i / N) * pw, x2 = ((i + 1) / N) * pw;
    const y1 = ph * Math.exp(-x1 / LTpx), y2 = ph * Math.exp(-x2 / LTpx);
    plainLine(s, px + x1, py + ph - y1, x2 - x1, y1 - y2, C.crowd, 2.25);
  }
  plainLine(s, px + LTpx, py + ph - ph * Math.exp(-1), 0, ph * Math.exp(-1),
            C.muted, 1, "dash");
  s.addText([{ text: "x = L", options: { italic: true } },
             { text: "T", options: { subscript: true, italic: true } }], {
    x: px + LTpx - 0.48, y: py + ph + 0.06, w: 0.96, h: 0.24, fontSize: 10,
    color: C.muted, align: "center", fontFace: FONT, margin: 0,
    valign: "middle",
  });
  s.addText([{ text: "I", options: { italic: true } },
             { text: "(x) = ", options: { italic: true } },
             { text: "I", options: { italic: true } },
             { text: "0", options: { subscript: true, italic: true } },
             { text: " exp( –x / " },
             { text: "L", options: { italic: true } },
             { text: "T", options: { subscript: true, italic: true } },
             { text: " )" }], {
    x: px + 0.85, y: py + 0.30, w: 2.15, h: 0.3, fontSize: 13, color: C.crowd,
    fontFace: FONT, margin: 0, valign: "middle",
  });
  label(s, "distance into the contact, x", px, py + ph + 0.36, pw,
        { fontSize: 10, color: C.muted, align: "center" });

  // the two competing transport paths
  s.addText([{ text: "L", options: { italic: true, bold: true } },
             { text: "T", options: { subscript: true, italic: true, bold: true } },
             { text: " = √( ρ", options: { bold: true } },
             { text: "C", options: { subscript: true, bold: true } },
             { text: " / R", options: { bold: true, italic: true } },
             { text: "S", options: { subscript: true, bold: true, italic: true } },
             { text: " )", options: { bold: true } },
             { text: "   —  lateral transport in the sheet competing with vertical transport across the interface" }], {
    x: 0.55, y: 4.62, w: 8.9, h: 0.4, fontSize: 12, color: C.text,
    align: "center", fontFace: FONT, margin: 0, valign: "middle",
  });
}

/* =====================================================================
 * SLIDE 4 - fig_tlm_plot
 * The extraction plot and what each feature means.
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "Reading four parameters off one straight line");

  const ox = 2.50, oy = 4.05;          // origin of the plot (L = 0, R = 0)
  const pw = 2.70, ph = 2.20;
  const b = 0.55;                       // y-intercept, inches
  const m = (ph - b) / pw;              // slope, inches per inch
  const xInt = -b / m;                  // x-intercept, inches (negative)

  // axes, with the L axis extended left of the origin for the x-intercept
  plainLine(s, ox + xInt - 0.30, oy, pw - xInt + 0.30, 0, C.edge, 1.5);
  plainLine(s, ox, oy + 0.30, 0, -(ph + 0.55), C.edge, 1.5);
  s.addText([{ text: "contact spacing, " }, { text: "L", options: { italic: true } }], {
    x: ox, y: oy + 0.42, w: pw, h: 0.28, fontSize: 12, color: C.text,
    align: "center", fontFace: FONT, margin: 0, valign: "middle",
  });
  s.addText([{ text: "total\nresistance, " }, { text: "R", options: { italic: true } },
             { text: "T", options: { subscript: true, italic: true } }], {
    x: 0.12, y: oy - ph - 0.05, w: 1.30, h: 0.60, fontSize: 12,
    color: C.text, align: "right", fontFace: FONT, margin: 0,
    valign: "middle",
  });

  // the fitted line, from the x-intercept to the right edge of the plot
  plainLine(s, ox + xInt, oy, pw - xInt, -ph, C.fit, 2.5);

  // data points along it
  [0.50, 1.05, 1.60, 2.15, 2.65].forEach((dx) => {
    const dy = b + m * dx;
    s.addShape(pres.ShapeType.ellipse, {
      x: ox + dx - 0.065, y: oy - dy - 0.065, w: 0.13, h: 0.13,
      fill: { color: C.data }, line: { color: C.data, width: 1 },
    });
  });

  // y-intercept = 2 Rc
  plainLine(s, ox + xInt, oy - b, pw - xInt, 0, C.muted, 1, "dash");
  s.addShape(pres.ShapeType.ellipse, {
    x: ox - 0.06, y: oy - b - 0.06, w: 0.12, h: 0.12,
    fill: { color: C.fit }, line: { color: C.fit, width: 1 },
  });
  s.addText([{ text: "y-intercept = 2" }, { text: "R", options: { italic: true } },
             { text: "C", options: { subscript: true, italic: true } }], {
    x: ox + 0.14, y: oy - b - 0.36, w: 1.9, h: 0.28, fontSize: 12,
    color: C.muted, fontFace: FONT, margin: 0, valign: "middle",
  });

  // x-intercept = -2 LT
  s.addShape(pres.ShapeType.ellipse, {
    x: ox + xInt - 0.06, y: oy - 0.06, w: 0.12, h: 0.12,
    fill: { color: C.crowd }, line: { color: C.crowd, width: 1 },
  });
  s.addText([{ text: "x-intercept" }], {
    x: ox + xInt - 0.85, y: oy + 0.10, w: 1.7, h: 0.24, fontSize: 11,
    color: C.crowd, align: "center", fontFace: FONT, margin: 0,
    valign: "middle",
  });
  s.addText([{ text: "= –2" }, { text: "L", options: { italic: true } },
             { text: "T", options: { subscript: true, italic: true } }], {
    x: ox + xInt - 0.85, y: oy + 0.34, w: 1.7, h: 0.24, fontSize: 11,
    color: C.crowd, align: "center", fontFace: FONT, margin: 0,
    valign: "middle",
  });

  // slope triangle
  const sx = ox + 1.20, sy = oy - (b + m * 1.20);
  plainLine(s, sx, sy, 0.95, 0, C.muted, 1, "dash");
  plainLine(s, sx + 0.95, sy, 0, -m * 0.95, C.muted, 1, "dash");
  s.addText([{ text: "slope = " }, { text: "R", options: { italic: true } },
             { text: "S", options: { subscript: true, italic: true } },
             { text: " / " }, { text: "W", options: { italic: true } }], {
    x: ox + 0.10, y: oy - ph - 0.32, w: 2.2, h: 0.28, fontSize: 12,
    color: C.muted, fontFace: FONT, margin: 0, valign: "middle",
  });

  // the equation, and what falls out of it
  s.addText([{ text: "R", options: { italic: true } },
             { text: "T", options: { subscript: true, italic: true } },
             { text: " = (" }, { text: "R", options: { italic: true } },
             { text: "S", options: { subscript: true, italic: true } },
             { text: " / " }, { text: "W", options: { italic: true } },
             { text: ") ( " }, { text: "L", options: { italic: true } },
             { text: " + 2" }, { text: "L", options: { italic: true } },
             { text: "T", options: { subscript: true, italic: true } },
             { text: " )" }], {
    x: 5.85, y: 1.20, w: 3.75, h: 0.42, fontSize: 17, bold: true, color: C.text,
    align: "center", fontFace: FONT, margin: 0, valign: "middle",
  });

  const rows = [
    [[{ text: "slope × " }, { text: "W", options: { italic: true } }], "sheet resistance"],
    [[{ text: "intercept ÷ 2" }], "contact resistance"],
    [[{ text: "intercept ÷ 2·slope" }], "transfer length"],
    [[{ text: "R", options: { italic: true } },
      { text: "S", options: { subscript: true, italic: true } },
      { text: " × " }, { text: "L", options: { italic: true } },
      { text: "T", options: { subscript: true, italic: true } },
      { text: "2", options: { superscript: true } }], "specific contact resistivity"],
  ];
  rows.forEach((r, i) => {
    const y = 1.90 + i * 0.55;
    s.addShape(pres.ShapeType.rect, {
      x: 5.85, y, w: 1.75, h: 0.45,
      fill: { color: C.panel }, line: { color: C.edge, width: 1 },
    });
    s.addText(r[0], {
      x: 5.85, y, w: 1.75, h: 0.45, fontSize: 10.5, color: C.text,
      align: "center", valign: "middle", fontFace: FONT, margin: 0.02,
    });
    label(s, "→", 7.62, y + 0.09, 0.24, { fontSize: 11, color: C.muted });
    label(s, r[1], 7.88, y + 0.09, 1.75, { fontSize: 10.5, color: C.muted });
  });
  label(s,
    "One fit, four numbers — provided the line really is straight.",
    5.85, 4.18, 3.75, { fontSize: 10, color: C.muted, italic: true,
                        align: "center" });
}

/* =====================================================================
 * SLIDE 5 - fig_tlm_measurement
 * How one point on that line is measured.
 * ===================================================================== */
{
  const s = pres.addSlide();
  title(s, "Measuring one point: a four-wire sweep on one contact pair");

  const yStrip = 2.95, hStrip = 0.80;
  const xL = 1.85, xR = 8.55, wPad = 1.20;
  const xC1 = xL + 0.40;                 // left contact, spans xC1..xC1+wPad
  const xC2 = xR - wPad - 0.40;          // right contact

  // sample stack
  box(s, xL, yStrip + hStrip, xR - xL, 0.32, C.sub, { line: C.grid });
  box(s, xL, yStrip, xR - xL, hStrip, C.semi);
  box(s, xC1, yStrip - 0.24, wPad, 0.24, C.metal, { line: C.edge });
  box(s, xC2, yStrip - 0.24, wPad, 0.24, C.metal, { line: C.edge });

  // the instrument
  s.addShape(pres.ShapeType.rect, {
    x: 3.75, y: 0.90, w: 2.50, h: 0.60,
    fill: { color: C.panel }, line: { color: C.edge, width: 1 },
  });
  s.addText([{ text: "source-measure unit\n", options: { bold: true, fontSize: 11 } },
             { text: "sweeps V, records I", options: { fontSize: 9.5, color: C.muted } }], {
    x: 3.75, y: 0.90, w: 2.50, h: 0.60, color: C.text, align: "center",
    valign: "middle", fontFace: FONT, margin: 0.02,
  });

  // force leads (solid) - carry the current, land on the outer edges
  const yForce = 1.95;
  plainLine(s, 4.20, 1.50, 0, yForce - 1.50, C.text, 1.5);
  plainLine(s, xC1 + 0.20, yForce, 4.20 - (xC1 + 0.20), 0, C.text, 1.5);
  plainLine(s, xC1 + 0.20, yForce, 0, yStrip - 0.24 - yForce, C.text, 1.5);
  plainLine(s, 5.80, 1.50, 0, yForce - 1.50, C.text, 1.5);
  plainLine(s, 5.80, yForce, xC2 + wPad - 0.20 - 5.80, 0, C.text, 1.5);
  plainLine(s, xC2 + wPad - 0.20, yForce, 0, yStrip - 0.24 - yForce, C.text, 1.5);
  label(s, "force +", xC1 - 1.05, yForce - 0.26, 1.15,
        { fontSize: 9.5, color: C.text, align: "right" });
  label(s, "force –", xC2 + wPad - 0.10, yForce - 0.26, 1.15,
        { fontSize: 9.5, color: C.text, align: "left" });

  // sense leads (dashed) - measure the voltage, land on the inner edges
  const ySense = 2.48;
  plainLine(s, 4.55, 1.50, 0, ySense - 1.50, C.muted, 1.2, "dash");
  plainLine(s, xC1 + wPad - 0.20, ySense, 4.55 - (xC1 + wPad - 0.20), 0,
            C.muted, 1.2, "dash");
  plainLine(s, xC1 + wPad - 0.20, ySense, 0, yStrip - 0.24 - ySense,
            C.muted, 1.2, "dash");
  plainLine(s, 5.45, 1.50, 0, ySense - 1.50, C.muted, 1.2, "dash");
  plainLine(s, 5.45, ySense, xC2 + 0.20 - 5.45, 0, C.muted, 1.2, "dash");
  plainLine(s, xC2 + 0.20, ySense, 0, yStrip - 0.24 - ySense, C.muted, 1.2, "dash");
  label(s, "sense +", xC1 - 1.05, ySense - 0.26, 1.15,
        { fontSize: 9.5, color: C.muted, align: "right" });
  label(s, "sense –", xC2 + wPad - 0.10, ySense - 0.26, 1.15,
        { fontSize: 9.5, color: C.muted, align: "left" });

  // contact spacing dimension, between the facing contact edges
  dimension(s, xC1 + wPad, xC2, yStrip + hStrip + 0.42, "");
  s.addText([{ text: "contact spacing " }, { text: "L", options: { italic: true } }], {
    x: xC1, y: yStrip + hStrip + 0.50, w: (xC2 + wPad) - xC1, h: 0.28,
    fontSize: 11, color: C.muted, align: "center", fontFace: FONT,
    margin: 0, valign: "middle",
  });

  s.addText([{ text: "Four-wire probing keeps the lead and probe resistance out of the measured " },
             { text: "R", options: { italic: true } },
             { text: "T", options: { subscript: true, italic: true } },
             { text: " — otherwise it lands in the intercept and inflates ρ" },
             { text: "C", options: { subscript: true } },
             { text: "." }], {
    x: 0.60, y: 4.68, w: 8.8, h: 0.32, fontSize: 11, color: C.text,
    align: "center", fontFace: FONT, margin: 0, valign: "middle",
  });
  s.addText([{ text: "Repeat over every spacing, then fit " },
             { text: "R", options: { italic: true } },
             { text: "T", options: { subscript: true, italic: true } },
             { text: " versus " }, { text: "L", options: { italic: true } },
             { text: "." }], {
    x: 0.60, y: 5.08, w: 8.8, h: 0.30, fontSize: 11, color: C.muted,
    italic: true, align: "center", fontFace: FONT, margin: 0,
    valign: "middle",
  });
}

pres.writeFile({ fileName: "figures.pptx" }).then(() => {
  console.log("wrote figures.pptx (5 slides)");
});
