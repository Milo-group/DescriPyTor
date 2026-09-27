/* The explorer's renderer: one molecule, or a stack of superposed structures,
 * freely rotated, with the paper's descriptor constructions drawn on it.
 *
 *   Sterimol axes  any number, each with its own layers: the axis and L, B1 out to its
 *                  plane, the B1 plane and its parallel through the axis, loc B5, B5,
 *                  theta built from v and its shadow, and an end-on inset (Sterimol.sterimol)
 *   dipole frame   origin (atom or centroid), u v w, mu and its projections
 *                  (Sterimol.dipoleFrame, the port of calc_dipole_gaussian)
 *   charges        chosen sites with their values, or the whole molecule as a map
 *   overlays       more structures, superposed beforehand (superpose), in one depth order
 *
 * Every number comes from the gated code; this file only draws, and superposes.
 * The view is always a proper rotation, so a drawing never shows the mirror image.
 */
(function (root) {
  'use strict';
  const W = 640, H = 540;
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const dif = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const scl = (a, k) => [a[0] * k, a[1] * k, a[2] * k];
  const len = (a) => Math.sqrt(dot(a, a));
  const unit = (a) => { const n = len(a) || 1; return scl(a, 1 / n); };
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const f1 = (x) => (Math.round(x * 10) / 10).toFixed(1);
  const esc = (s) => String(s).replace(/&(?![#\w]+;)/g, '&amp;').replace(/</g, '&lt;').replace(/"/g, '&quot;');
  const sub = (k) => `<tspan dy="2.6" font-size="0.75em">${k}</tspan><tspan dy="-2.6">`;   // closed by the caller's </tspan>

  // ---------- elements ----------
  const COV = { H: .31, C: .76, N: .71, O: .66, F: .57, S: 1.05, Cl: 1.02, Br: 1.20, I: 1.39, P: 1.07, Cu: 1.32,
    Pd: 1.39, Ni: 1.24, Fe: 1.32, Co: 1.26, Rh: 1.42, Ir: 1.41, Pt: 1.36, Ru: 1.46, B: .84, Si: 1.11 };
  const BALL = { H: '#F3F5F8', C: '#A9B1BC', N: '#4A67C9', O: '#CF4A40', F: '#79AE63', S: '#DDB035', Cl: '#5E9F59',
    Br: '#A5552A', I: '#8C4BA8', P: '#DF8A2E', Cu: '#C8804A', Pd: '#7A8FA6', Ni: '#6FA38F', B: '#E7A987', Si: '#C9B99A',
    Fe: '#C4643C', Co: '#5E86C4', Mn: '#9C63B8', Cr: '#5FA5B5', Zn: '#7D8FA8', Ru: '#4E9E8F', Rh: '#9C6BA8', Ag: '#B9C0C8',
    Ir: '#5C7FA8', Pt: '#9AA8B8', Au: '#D3A73A', Mo: '#5E96A8', W: '#6E8CA8', Ti: '#9AA5AE', V: '#A0785C', Mg: '#84C46A',
    Ca: '#7FB36A', Li: '#C57BC0', Na: '#B07AD0', K: '#9A6FD0', Sn: '#8E9AA8', Pb: '#7B8390', Hg: '#9FA6B4', Cd: '#B79A7A' };
  const VDW = { H: 1.1, C: 1.7, N: 1.55, O: 1.52, F: 1.47, S: 1.8, Cl: 1.75, Br: 1.85, I: 1.98, P: 1.8, B: 1.92, Si: 2.1,
    Cu: 1.96, Ni: 1.97, Pd: 2.1, Fe: 2.04, Co: 2.0, Mn: 2.05, Cr: 2.06, Zn: 2.01, Ru: 2.13, Rh: 2.1, Ag: 2.11, Ir: 2.13,
    Pt: 2.09, Au: 2.14, Mo: 2.17, W: 2.18, Ti: 2.11, V: 2.07, Mg: 1.73, Ca: 2.31, Li: 1.82, Na: 2.27, K: 2.75,
    Sn: 2.17, Pb: 2.02, Hg: 2.23, Cd: 2.18 };
  const METAL = new Set(['Li', 'Na', 'K', 'Mg', 'Ca', 'Al', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Zr', 'Mo', 'Ru', 'Rh', 'Pd',
    'Ag', 'Cd', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Sn', 'Pb']);
  const covalent = (e) => COV[e] || .77;
  const ballColor = (e) => BALL[e] || '#9AA3AF';
  function shade(hex, k) {                 // k < 0 darkens, k > 0 lightens
    const n = parseInt(hex.slice(1), 16), c = [n >> 16, (n >> 8) & 255, n & 255];
    const t = k < 0 ? 0 : 255, a = Math.abs(k);
    return '#' + c.map((v) => Math.round(v + (t - v) * a).toString(16).padStart(2, '0')).join('');
  }

  // ---------- labels: the first candidate spot that overlaps least ----------
  const FS = { lab: 12, th: 13, vl: 11.5, key: 10.5, mark: 9, dpn: 14, alab: 8 };
  class Placer {
    // shift: {key: [dx, dy]} — where the reader has dragged a label, as an offset from where
    // the placer would have put it, so a nudge survives the figure being turned or re-measured
    constructor(w = W, h = H, k = 1, shift = {}, leaders = true) {
      this.w = w; this.h = h; this.k = k; this.circles = []; this.boxes = []; this.shift = shift; this.seen = {}; this.leaders = leaders;
    }
    key(s, cls) {                              // the same label across redraws, values ignored
      const kind = (cls + '|' + s)
        .replace(/-?\d+\.\d+/g, '')                                 // the value changes; the label is the same label
        .replace(/-?\d+(?=\s*(?:&#8491;|&#176;|°))/g, '')
        .replace(/<[^>]*>/g, '').replace(/&#?[0-9a-zA-Z]+;/g, '').replace(/[^\w|]+/g, '');
      this.seen[kind] = (this.seen[kind] || 0) + 1;
      return kind + '|' + this.seen[kind];
    }
    size(cls) { return (FS[cls.split(' ')[0]] || 12) * this.k; }
    box(x, y, s, cls, anchor) {
      const n = s.replace(/&#?[0-9a-zA-Z]+;/g, 'x').replace(/<[^>]+>/g, '').length;
      const fs = this.size(cls), w = 0.56 * fs * n;
      const x0 = x - (anchor === 'end' ? w : anchor === 'middle' ? w / 2 : 0);
      return [x0 - 4, y - fs, x0 + w + 4, y + fs / 2];
    }
    cost([x0, y0, x1, y1]) {
      let c = (x0 >= 4 && y0 >= 4 && x1 <= this.w - 4 && y1 <= this.h - 4) ? 0 : 1e4;
      for (const [cx, cy, r] of this.circles) {
        const d = Math.hypot(cx - Math.min(Math.max(cx, x0), x1), cy - Math.min(Math.max(cy, y0), y1));
        c += Math.max(0, r - d) * 10;
      }
      for (const [a0, b0, a1, b1] of this.boxes) c += Math.max(0, Math.min(x1, a1) - Math.max(x0, a0)) * Math.max(0, Math.min(y1, b1) - Math.max(y0, b0));
      return c;
    }
    put(cands, s, cls, attrs = '', anchor = null) {
      let best = null;
      // with no leader to follow, a label has to sit on the thing it names: pull it back towards it
      const pull = this.leaders === false && anchor ? 0.9 : 0;
      const score = (x, y, a) => {
        const b0 = this.box(x, y, s, cls, a);
        let c = this.cost(b0);
        if (pull) {
          const cx = Math.min(Math.max(anchor[0], b0[0]), b0[2]), cy = Math.min(Math.max(anchor[1], b0[1]), b0[3]);
          c += pull * Math.hypot(anchor[0] - cx, anchor[1] - cy);
        }
        if (!best || c < best.c) best = { c, x, y, a };
      };
      cands.forEach(([x, y, a]) => score(x, y, a));
      if (anchor && best.c > 80) {                    // every offered spot is buried: sweep rings further out
        for (const rad of [30, 44, 60, 80]) for (let k = 0; k < 12; k++) {
          const t = k * Math.PI / 6;
          score(anchor[0] + rad * Math.cos(t), anchor[1] + rad * Math.sin(t) + 4,
            Math.cos(t) < -0.3 ? 'end' : Math.cos(t) > 0.3 ? 'start' : 'middle');
        }
      }
      const key = this.key(s, cls), off = this.shift[key];
      if (off) { best = { c: best.c, x: best.x + off[0], y: best.y + off[1], a: best.a }; }
      const box = this.box(best.x, best.y, s, cls, best.a);
      this.boxes.push(box);
      const text = `<text class="${cls}" data-lab="${key}" x="${f1(best.x)}" y="${f1(best.y)}" text-anchor="${best.a}" font-size="${f1(this.size(cls))}"${attrs}>${s}</text>`;
      if (!anchor || this.leaders === false) return text;      // leaders off: the label stands on its own
      const cx = Math.min(Math.max(anchor[0], box[0]), box[2]), cy = Math.min(Math.max(anchor[1], box[1]), box[3]);
      if (!off && Math.hypot(anchor[0] - cx, anchor[1] - cy) < 18) return text;   // it already sits by what it names
      if (off && Math.hypot(anchor[0] - cx, anchor[1] - cy) < 12) return text;
      const col = (attrs.match(/fill="([^"]+)"/) || [])[1] || '#6B7482';
      return `<line class="leader" x1="${f1(cx)}" y1="${f1(cy)}" x2="${f1(anchor[0])}" y2="${f1(anchor[1])}" stroke="${col}"/>${text}`;
    }
    along(p, q, r = 3.5, step = 7) {
      const n = Math.max(2, Math.floor(Math.hypot(q[0] - p[0], q[1] - p[1]) / step));
      for (let i = 0; i < n; i++) { const t = i / (n - 1); this.circles.push([p[0] + (q[0] - p[0]) * t, p[1] + (q[1] - p[1]) * t, r]); }
    }
  }

  // ---------- colours and looks ----------
  const PAL = {                            // the paper's conventions (Fig. 2, case studies 1-3)
    axes: ['#1596A6', '#E07A1F', '#C9352B', '#7B57C9', '#3E8E41', '#B8338A'],
    sites: ['#C9352B', '#E58A1F', '#2E9E5B', '#7A3E9D', '#1596A6'],
    structs: ['#2F6FB5', '#D9822B', '#3E9E6A', '#B8455C', '#7A5CC4', '#8C7A3A', '#2A9FA8', '#C25BA7'],
    qpos: '#C9352B', qneg: '#2F5BD3',
  };
  // what the constructions are drawn in: the Descriptors colours, the theta figure's inks, or black
  const INKS = {
    colour: { label: 'Descriptor colours', axis: null, b1: '#6C4AB6', plane: '#3450A8', loc: '#C99A00', locText: '#8A6A00', b5: '#C9352B',
      theta: '#3450A8', right: '#434C59', u: '#2F5BD3', v: '#D1453B', w: '#2E9E5B', mu: '#7A3E9D', origin: '#E58A1F', text: '#1F2937', key: '#374151', eq: ' ' },
    paper: { label: 'θ construction (Fig. 2)', axis: '#151A21', b1: '#434C59', plane: '#3450A8', loc: '#151A21', locText: '#434C59', b5: '#93600F',
      theta: '#3450A8', right: '#434C59', u: '#2F5BD3', v: '#D1453B', w: '#2E9E5B', mu: '#7A3E9D', origin: '#E58A1F', text: '#434C59', key: '#434C59', eq: ' = ', draft: true },
    mono: { label: 'Black and white', axis: '#111111', b1: '#111111', plane: '#5C5C5C', loc: '#111111', locText: '#111111', b5: '#444444',
      theta: '#111111', right: '#444444', u: '#111111', v: '#555555', w: '#8A8A8A', mu: '#000000', origin: '#111111', text: '#111111', key: '#111111', eq: ' = ', draft: true, mono: true },
  };
  // atoms and bonds, after CYLview's family of looks (and the theta figure's own)
  const PRESETS = {
    construction: { label: 'Construction (θ figure)', atom: 'small', shade: 'paper', bonds: true, bondW: 0.3, bondGrey: '#6B7482', olc: '#8F99A7', bondOl: false },
    glossy: { label: 'Glossy', atom: 'ball', shade: 'glossy', bonds: true },
    matte: { label: 'Matte', atom: 'ball', shade: 'matte', bonds: true },
    flat: { label: 'Flat line art', atom: 'ball', shade: 'flat', bonds: true, outline: 1.9 },
    toon: { label: 'Toon', atom: 'ball', shade: 'toon', bonds: true, outline: 2.4 },
    tube: { label: 'Tubes', atom: 'bond', shade: 'glossy', bonds: true, bondW: 1.35, tube: true },
    spacefill: { label: 'Space-filling', atom: 'vdw', shade: 'glossy', bonds: false },
    wire: { label: 'Wireframe', atom: 'none', shade: 'flat', bonds: true, bondW: 0.32, wire: true },
  };
  const PALETTES = {
    soft: { label: 'Soft (θ figure)' },
    cylview: { label: 'CYLview-like', C: '#6E757E', H: '#F4F4F4', N: '#2F55D4', O: '#E0301E', S: '#E2BF2B', F: '#8FD06A', Cl: '#3DB34A',
      Br: '#9A3A2A', I: '#7B2FA0', P: '#F08C1E', Cu: '#C57A36', Ni: '#58A08C', Pd: '#5B7FA8', B: '#F2A58A', Si: '#C8B48C' },
    jmol: { label: 'Jmol', C: '#909090', H: '#FFFFFF', N: '#3050F8', O: '#FF0D0D', F: '#90E050', S: '#FFFF30', Cl: '#1FF01F', Br: '#A62929',
      I: '#940094', P: '#FF8000', Cu: '#C88033', Ni: '#50D050', Pd: '#006985', B: '#FFB5B5', Si: '#F0C8A0' },
    paper: { label: 'Paper (CYLview default)', C: '#BFC3C8', H: '#F7F7F7', N: '#4F5BEA', O: '#E3241B', S: '#E8C630', F: '#9BE07A',
      Cl: '#36C24A', Br: '#A0463A', I: '#8E3FB0', P: '#F2951E', Cu: '#C98A4B', Co: '#EAA3BE', Ni: '#6FB9A2', Pd: '#7C98B8', B: '#F4AE8C', Si: '#D0BD94' },
    print: { label: 'Greyscale', C: '#8A8A8A', H: '#F5F5F5', N: '#3F3F3F', O: '#626262', S: '#BDBDBD', F: '#D2D2D2', Cl: '#A8A8A8',
      Br: '#7A7A7A', I: '#555555', P: '#9A9A9A', B: '#E0E0E0', Si: '#B0B0B0', fallback: '#4A4A4A' },
  };
  const STYLE = { axisMark: 'arrow', look: 'descriptors', preset: 'glossy', palette: 'soft', ink: 'colour', hydrogens: 'all', bondColor: 'split', atomLabels: 'off', values: true, featureLabels: true, vdw: 'off', vdwOpacity: 0.32, vdwScale: 1, surfKind: 'cpk', surfIso: 0.003, surfStep: 0.4, leaders: true, planeNote: true, locMark: 'tick',
    dipoleArrow: 'gaussian', muPerD: 0, qRamp: 'bwr', qMax: 0, hiddenLines: true,
    atomScale: 1, bondScale: 1, outline: 1, fog: 0.3, perspective: 0, labelScale: 1, shadow: true, metalDash: true, bg: 'white',
    overlay: 'tint', overlayOpacity: 0.6 };
  // a look sets these together; everything else in a style is left as the viewer had it
  const LOOK_KEYS = ['preset', 'palette', 'ink', 'bondColor', 'shadow', 'fog', 'outline', 'atomScale', 'bondScale', 'perspective'];
  const LOOKS = {
    paper: { label: 'Paper (CYLview default)', preset: 'glossy', palette: 'paper', ink: 'colour', bondColor: 'black', shadow: false, fog: 0,
      outline: 0.6, atomScale: 0.72, bondScale: 0.62 },
    construction: { label: 'θ construction (Fig. 2)', preset: 'construction', palette: 'soft', ink: 'paper', bondColor: 'grey', shadow: false, fog: 0 },
    descriptors: { label: 'Descriptors (glossy)', preset: 'glossy', palette: 'soft', ink: 'colour', shadow: true, fog: 0.3 },
    cylview: { label: 'CYLview classic', preset: 'glossy', palette: 'cylview', ink: 'colour', shadow: true, fog: 0.45, outline: 1.2, perspective: 0.35 },
    matte: { label: 'Matte, deep', preset: 'matte', palette: 'cylview', ink: 'colour', shadow: true, fog: 0.7, perspective: 0.6 },
    lineart: { label: 'Flat line art', preset: 'flat', palette: 'soft', ink: 'colour', bondColor: 'grey', shadow: false, fog: 0 },
    toon: { label: 'Toon', preset: 'toon', palette: 'cylview', ink: 'colour', shadow: false, fog: 0 },
    tubes: { label: 'Tubes', preset: 'tube', palette: 'jmol', ink: 'colour', shadow: true, fog: 0.3 },
    spacefill: { label: 'Space-filling', preset: 'spacefill', palette: 'cylview', ink: 'colour', shadow: true, fog: 0.2 },
    print: { label: 'Print (greyscale)', preset: 'flat', palette: 'print', ink: 'mono', shadow: false, fog: 0 },
  };
  const lookStyle = (k) => Object.fromEntries(LOOK_KEYS.map((key) => [key, (LOOKS[k] || {})[key] ?? STYLE[key]]));
  const QNAME = { nbo: 'NBO', cm5: 'CM5', hirshfeld: 'Hirsh' };

  function spheres(colors, kind) {
    return [...colors].map((c) => {
      const id = `sg${kind}${c.slice(1)}`;
      if (kind === 'paper') return `<radialGradient id="${id}" cx="36%" cy="32%" r="72%"><stop offset="0" stop-color="${shade(c, 0.78)}"/><stop offset="0.38" stop-color="${c}"/><stop offset="1" stop-color="${shade(c, -0.42)}"/></radialGradient>`;
      if (kind === 'glossy') return `<radialGradient id="${id}" cx="34%" cy="30%" r="70%"><stop offset="0" stop-color="#FFFFFF"/><stop offset="0.13" stop-color="${shade(c, 0.6)}"/><stop offset="0.5" stop-color="${c}"/><stop offset="1" stop-color="${shade(c, -0.55)}"/></radialGradient>`;
      return `<radialGradient id="${id}" cx="40%" cy="36%" r="76%"><stop offset="0" stop-color="${shade(c, 0.35)}"/><stop offset="0.62" stop-color="${c}"/><stop offset="1" stop-color="${shade(c, -0.3)}"/></radialGradient>`;
    }).join('');
  }
  // bonds to draw: the package's own, plus non-metal pairs at covalent length it leaves unbonded (drawn dashed).
  // Metal bonds follow the package's rule (2.8 A, six per metal) and are never second-guessed here.
  function displayBonds(el, X, bonds) {
    const key = (a, b) => Math.min(a, b) * 100000 + Math.max(a, b);
    const has = new Set(bonds.map(([a, b]) => key(a, b))), extra = [];
    for (let i = 0; i < el.length; i++) {
      for (let j = i + 1; j < el.length; j++) {
        if ((el[i] === 'H' && el[j] === 'H') || METAL.has(el[i]) || METAL.has(el[j])) continue;
        const d = len(dif(X[i], X[j]));
        if (!has.has(key(i + 1, j + 1)) && d < 1.15 * (covalent(el[i]) + covalent(el[j]))) extra.push([i + 1, j + 1, d]);
      }
    }
    return { all: bonds.map(([a, b]) => [a, b, 'pkg']).concat(extra.map(([a, b]) => [a, b, 'extra'])), extra };
  }
  // the charge map: blue-white-red by default, or a purple-orange pair that survives red-green colour blindness
  const RAMPS = { bwr: { neg: PAL.qneg, pos: PAL.qpos, label: 'blue &#8211; white &#8211; red' },
    puor: { neg: '#542788', pos: '#B35806', label: 'purple &#8211; white &#8211; orange (colour-blind safe)' } };
  function divergent(q, qmax, ramp) {
    const t = Math.max(-1, Math.min(1, q / (qmax || 1))), r = RAMPS[ramp] || RAMPS.bwr;
    return shade(t < 0 ? r.neg : r.pos, 0.9 - 0.9 * Math.abs(t));
  }

  // ---------- one Sterimol axis, everything it contributes, in molecule coordinates ----------
  const toWorld = (r, f) => [0, 1, 2].map((k) => r.origin[k] + f[0] * r.basis[0][k] + f[1] * r.basis[1][k] + f[2] * r.basis[2][k]);
  // over = {n, B1}: the construction at one of the scan's other trials, for replaying it
  let wedgeScale = 0.58;                  // radius of the theta sector as a share of |v|; style.thetaBold widens it
  function axisGeometry(r, over) {
    const F = r.frame, idx = r.atoms, R = r.radii, B1 = over ? over.B1 : r.B1;
    // orient the normal towards the supporting plane, so B1 out along it lands on the plane
    let n = (over ? over.n : r.normal).slice();
    const supOf = (nn) => Math.max(...idx.map((i, j) => F[i][0] * nn[0] + F[i][2] * nn[1] + R[j]));
    if (Math.abs(supOf(n) - B1) > Math.abs(supOf([-n[0], -n[1]]) - B1)) n = [-n[0], -n[1]];
    const sup = idx.map((i, j) => F[i][0] * n[0] + F[i][2] * n[1] + R[j]);
    const smax = Math.max(...sup), bc = F[idx[sup.indexOf(smax)]][1];
    const v = F[r.b5], pm = Math.hypot(v[0], v[2]) || 1, q = [v[0] / pm, v[2] / pm], t = [-n[1], n[0]];
    // frame points by (across the axis in the plane, along the axis, out along the normal)
    const at = (a, y, o) => [a * t[0] + o * n[0], y, a * t[1] + o * n[1]];
    const w = (f) => toWorld(r, f), base = w([0, 0, 0]), dir = (f) => dif(w(f), base);
    const vN = v[0] * n[0] + v[2] * n[1], vT = v[0] * t[0] + v[2] * t[1];
    const vw = w(v), Sw = w(at(vT, v[1], 0));              // v, and its shadow in the plane through the axis
    const u1 = unit(dif(Sw, base)), vh = unit(dif(vw, base));
    const wv = unit(dif(vh, scl(u1, dot(u1, vh))));
    const th = r.B1_B5_angle * Math.PI / 180, rr0 = wedgeScale * len(dif(vw, base));
    const ring = Array.from({ length: 32 }, (_, k) => { const a = th * k / 31; return add(base, add(scl(u1, rr0 * Math.cos(a)), scl(wv, rr0 * Math.sin(a)))); });
    const Nw = dir([n[0], 0, n[1]]), Aw = dir([0, 1, 0]), Qw = dir([q[0], 0, q[1]]);
    const dn = scl(Nw, Math.sign(vN) || 1), k0 = 0.3, k2 = 0.22;
    const Yw = w([0, v[1], 0]), E5 = add(Yw, scl(Qw, r.B5));
    const lo = Math.min(0, vT) - 0.55, hi = Math.max(0, vT) + 0.55, ph = 0.45 + 0.08 * r.L;
    // phi: read where B5 is read, across the axis, from the B1 contact round to B5
    const uN = unit(Nw), rp = 0.45 * r.B5, qn = dif(Qw, scl(uN, dot(uN, Qw)));
    const wq = len(qn) > 1e-6 ? unit(qn) : unit(dir([t[0], 0, t[1]]));       // B5 on the normal: any plane will do
    const phA = Math.acos(Math.max(-1, Math.min(1, dot(uN, Qw))));
    return {
      phiAt: Yw, phiN: add(Yw, scl(uN, rp)), phiQ: add(Yw, scl(Qw, rp)),
      phiRing: Array.from({ length: 32 }, (_, k) => {
        const a = phA * k / 31;
        return add(Yw, add(scl(uN, rp * Math.cos(a)), scl(wq, rp * Math.sin(a))));
      }),
      n, bc, touch: idx.filter((_, j) => sup[j] >= smax - 0.03),
      lAtoms: idx.filter((i, j) => F[i][1] + R[j] >= r.L - 0.01),      // L is a max over the fragment: equivalent atoms tie
      base, tip: w([0, r.L, 0]),
      b1From: w([0, bc, 0]), b1To: w(at(0, bc, B1)), plane: [w(at(-0.55, bc, B1)), w(at(0.55, bc, B1))], B1,
      loc: Yw, b5To: E5, v: vw, S: Sw, u1, rr0, ring,
      right: [add(Sw, scl(u1, -k0)), add(add(Sw, scl(u1, -k0)), scl(dn, k0)), add(Sw, scl(dn, k0))],
      right5: [add(Yw, scl(Qw, k2)), add(add(Yw, scl(Qw, k2)), scl(Aw, k2)), add(Yw, scl(Aw, k2))],
      ref: [w(at(lo, -0.3, 0)), w(at(hi, -0.3, 0)), w(at(hi, r.L + 0.3, 0)), w(at(lo, r.L + 0.3, 0))],
      patch: [w(at(-ph, bc - ph, B1)), w(at(ph, bc - ph, B1)), w(at(ph, bc + ph, B1)), w(at(-ph, bc + ph, B1))],
      axisDir: Aw, normal: Nw,
    };
  }
  // the view the theta figure is drawn in: the axis up, turned so v stands in front of the plane,
  // then tipped toward the viewer. Rows: screen x, screen y, toward the viewer; a proper rotation.
  function axisView(ax, psiDeg = 38, epsDeg = 17) {
    const g = ax.g || axisGeometry(ax.res), Y = g.axisDir, N = g.normal, T = cross(Y, N), v = dif(g.v, g.base);
    let psi = psiDeg * Math.PI / 180;
    if (-Math.sin(psi) * dot(T, v) + Math.cos(psi) * dot(N, v) < 0) psi += Math.PI;
    const cp = Math.cos(psi), sp = Math.sin(psi), e = epsDeg * Math.PI / 180, ce = Math.cos(e), se = Math.sin(e);
    const xr = add(scl(T, cp), scl(N, sp)), zr = add(scl(T, -sp), scl(N, cp));
    return [xr, add(scl(Y, ce), scl(zr, -se)), add(scl(Y, se), scl(zr, ce))];
  }

  // ---------- superposition (Horn 1987, unit quaternions): never a mirror ----------
  function eig4top(A) {                    // Jacobi; the eigenvector of the largest eigenvalue
    const a = A.map((r) => r.slice()), V = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
    for (let sweep = 0; sweep < 64; sweep++) {
      let off = 0;
      for (let p = 0; p < 4; p++) for (let q = p + 1; q < 4; q++) off += a[p][q] * a[p][q];
      if (off < 1e-24) break;
      for (let p = 0; p < 4; p++) for (let q = p + 1; q < 4; q++) {
        if (a[p][q] === 0) continue;
        const th = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        const t = (th >= 0 ? 1 : -1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < 4; k++) { const x = a[k][p], y = a[k][q]; a[k][p] = c * x - s * y; a[k][q] = s * x + c * y; }
        for (let k = 0; k < 4; k++) { const x = a[p][k], y = a[q][k]; a[p][k] = c * x - s * y; a[q][k] = s * x + c * y; }
        for (let k = 0; k < 4; k++) { const x = V[k][p], y = V[k][q]; V[k][p] = c * x - s * y; V[k][q] = s * x + c * y; }
      }
    }
    let b = 0;
    for (let k = 1; k < 4; k++) if (a[k][k] > a[b][b]) b = k;
    return V.map((row) => row[b]);
  }
  // the rigid motion carrying the mobile points B onto A (least squares); apply() moves any point of B's structure
  function superpose(A, B) {
    const m = A.length, cen = (Z) => [0, 1, 2].map((k) => Z.reduce((s, p) => s + p[k], 0) / m), ca = cen(A), cb = cen(B);
    const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < m; i++) { const b = dif(B[i], cb), a = dif(A[i], ca); for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) M[r][c] += b[r] * a[c]; }
    const [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]] = M;
    const [w, x, y, z] = eig4top([[xx + yy + zz, yz - zy, zx - xz, xy - yx], [yz - zy, xx - yy - zz, xy + yx, zx + xz],
      [zx - xz, xy + yx, -xx + yy - zz, yz + zy], [xy - yx, zx + xz, yz + zy, -xx - yy + zz]]);
    const R = [[w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y)],
      [2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x)],
      [2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z]];
    const apply = (p) => { const d = dif(p, cb); return add(ca, R.map((row) => dot(row, d))); };
    let sq = 0;
    for (let i = 0; i < m; i++) { const d = dif(apply(B[i]), A[i]); sq += dot(d, d); }
    return { R, apply, rmsd: Math.sqrt(sq / m) };
  }

  // ---------- the end-on inset: looking along the axis from the base ----------
  function inset(ax, x0, y0, side, ink, k, V = true) {
    const g = ax.g, r = ax.res, n = g.n, B1 = r.B1, B5 = r.B5, F = r.frame, idx = r.atoms, R = r.radii, a5 = idx.indexOf(r.b5);
    const sc = (side / 2 - 8) / (B5 * 1.06), cx = x0 + side / 2, cy = y0 + side / 2;
    const pt = (p) => [cx + p[0] * sc, cy - p[1] * sc];
    const L2 = (p, q, col, wd, extra = '') => `<line x1="${f1(p[0])}" y1="${f1(p[1])}" x2="${f1(q[0])}" y2="${f1(q[1])}" stroke="${col}" stroke-width="${wd}" stroke-linecap="round"${extra}/>`;
    const o = [`<rect class="ins-box" x="${f1(x0)}" y="${f1(y0)}" width="${side}" height="${side}" rx="3"/>`,
      `<circle cx="${f1(cx)}" cy="${f1(cy)}" r="${f1(B5 * sc)}" fill="none" stroke="${ink.b5}" stroke-width=".9" stroke-dasharray="3 4" opacity=".75"/>`];
    idx.forEach((i, j) => {
      const [x, y] = pt([F[i][0], F[i][2]]);
      o.push(`<circle cx="${f1(x)}" cy="${f1(y)}" r="${f1(R[j] * sc)}" fill="#D9DFE8" fill-opacity="${j === a5 ? '.5' : '.35'}" stroke="${j === a5 ? ink.b5 : '#8F99A7'}" stroke-width="${j === a5 ? 1.2 : 0.5}"/>`);
    });
    const t = [-n[1], n[0]];
    o.push(L2(pt([n[0] * B1 - t[0] * B5, n[1] * B1 - t[1] * B5]), pt([n[0] * B1 + t[0] * B5, n[1] * B1 + t[1] * B5]), ink.b1, 1.4),
      `<circle cx="${f1(cx)}" cy="${f1(cy)}" r="1.8" fill="${ink.text}"/>`);
    const p5 = [F[r.b5][0], F[r.b5][2]], m5 = Math.hypot(p5[0], p5[1]) || 1, u5 = [p5[0] / m5, p5[1] / m5];
    o.push(L2([cx, cy], pt([n[0] * B1, n[1] * B1]), ink.b1, 1.3, ' opacity=".85"'), L2([cx, cy], pt([u5[0] * B5, u5[1] * B5]), ink.b5, 1.3, ' opacity=".85"'));
    const tk = [-u5[1] * 0.09 * B5, u5[0] * 0.09 * B5];
    o.push(L2(pt([u5[0] * B5 - tk[0], u5[1] * B5 - tk[1]]), pt([u5[0] * B5 + tk[0], u5[1] * B5 + tk[1]]), ink.b5, 1.3));
    const fs = f1(FS.mark * k);
    if ((ax.show || {}).phi && r.rho > 1e-6) {   // phi: the azimuth of B5 from the B1 normal, in this plane
      const d1 = [n[0], -n[1]], d2 = [u5[0], -u5[1]], rr = Math.min(0.3 * B5 * sc, side * 0.16);
      const P = (d, m) => [cx + d[0] * m, cy + d[1] * m], sweep = n[1] * u5[0] - n[0] * u5[1] > 0 ? 1 : 0;
      o.push(`<path class="sx-phi" d="M${f1(P(d1, rr)[0])},${f1(P(d1, rr)[1])} A${f1(rr)} ${f1(rr)} 0 0 ${sweep} ${f1(P(d2, rr)[0])},${f1(P(d2, rr)[1])}" fill="none" stroke="${ink.theta}"/>`,
        // the value goes with the box's other values, not on the arc: there is no room inside it
        `<text class="mark" x="${f1(x0 + 6)}" y="${f1(y0 + side - 27)}" font-size="${fs}" fill="${ink.theta}">&#966;${V ? ' ' + r.phi.toFixed(1) + '&#176;' : ''}</text>`);
    }
    o.push(`<text class="mark" x="${f1(cx)}" y="${f1(y0 - 5)}" text-anchor="middle" font-size="${fs}" fill="${ink.text}">end-on, ${esc(ax.label || '')}</text>`,
      `<text class="mark" x="${f1(x0 + 6)}" y="${f1(y0 + side - 16)}" font-size="${fs}" fill="${ink.b1}">B${sub(1)}${V ? ' ' + B1.toFixed(2) : ''}</tspan></text>`,
      `<text class="mark" x="${f1(x0 + 6)}" y="${f1(y0 + side - 5)}" font-size="${fs}" fill="${ink.b5}">B${sub(5)}${V ? ' ' + B5.toFixed(2) + ' &#8491;' : ''}</tspan></text>`);
    return o;
  }

  // B1 against the rotation the scan walks: the minimum the number comes from, and how flat it is there
  function profilePlot(ax, x0, y0, w, h, ink, k, V) {   // and where the replay currently is
    const r = ax.res, prof = r.profile || [];
    if (!prof.length) return [];
    const B = prof.map((p) => p.B1), lo = Math.min(...B), hi = Math.max(...B), span = (hi - lo) || 1;
    const px = (i) => x0 + 8 + (w - 16) * i / (prof.length - 1);
    const py = (b) => y0 + h - 16 - (h - 26) * (b - lo) / span;
    const flat = prof.filter((p) => p.B1 <= lo + 0.05).length;            // degrees within 0.05 A of the minimum
    const o = [`<rect class="ins-box" x="${f1(x0)}" y="${f1(y0)}" width="${f1(w)}" height="${f1(h)}" rx="3"/>`];
    o.push(`<path class="pf-band" d="M${f1(x0 + 8)},${f1(py(lo + 0.05))} L${f1(x0 + w - 8)},${f1(py(lo + 0.05))} L${f1(x0 + w - 8)},${f1(py(lo))} L${f1(x0 + 8)},${f1(py(lo))} Z" fill="${ink.b1}"/>`);
    o.push(`<path class="pf-line" d="M${prof.map((p, i) => f1(px(i)) + ',' + f1(py(p.B1))).join(' L')}" stroke="${ink.b1}"/>`);
    if (ax.preview) {                                  // the trial being replayed
      const j = prof.findIndex((p) => p.deg === ax.preview.deg);
      if (j >= 0) o.push(`<line class="pf-now" x1="${f1(px(j))}" y1="${f1(y0 + 6)}" x2="${f1(px(j))}" y2="${f1(y0 + h - 14)}" stroke="${ink.b5}"/>`,
        `<circle class="pf-dot" cx="${f1(px(j))}" cy="${f1(py(prof[j].B1))}" r="2.6" fill="${ink.b5}"/>`);
    }
    const i0 = prof.reduce((b, p, i) => (p.B1 < prof[b].B1 ? i : b), 0);
    o.push(`<circle class="pf-dot" cx="${f1(px(i0))}" cy="${f1(py(r.B1))}" r="3" fill="${ink.b1}"/>`,
      `<line class="pf-drop" x1="${f1(px(i0))}" y1="${f1(py(r.B1))}" x2="${f1(px(i0))}" y2="${f1(y0 + h - 14)}" stroke="${ink.b1}"/>`);
    const fs = f1(FS.mark * k);
    o.push(`<line class="pf-drop" x1="${f1(x0 + 8)}" y1="${f1(py(r.B1))}" x2="${f1(x0 + w - 8)}" y2="${f1(py(r.B1))}" stroke="${ink.b1}"/>`);
    o.push(`<text class="mark" x="${f1(x0 + w / 2)}" y="${f1(y0 - 5)}" text-anchor="middle" font-size="${fs}" fill="${ink.text}">B${sub(1)} over the scan, ${esc(ax.label || '')}</tspan></text>`,
      `<text class="mark" x="${f1(x0 + w / 2)}" y="${f1(y0 + h - 4)}" text-anchor="middle" font-size="${fs}" fill="${ink.text}">${V ? `min ${r.B1.toFixed(2)} &#8491; at ${r.scanDeg}&#176;, flat ${flat}&#176;` : `flat ${flat}&#176; of 90&#176;`}</text>`);
    return o;
  }

  /* mol: {el, X, bonds, displayBonds?, types?, charges?}
   * L:   {axes: [{res, color, label, tag, show: {axis, L, b1, plane, loc, b5 (implies loc), theta, drop, inset}}],
   *       dipole: {res, show: {frame, mu, u, v, w}, originSet}, sites: [{i, type, color}], map: {type} | null,
   *       fade: Set, legend: bool, style: {...STYLE}, primaryName: '',
   *       overlays: [{id, name, el, X (already superposed), bonds: [[a, b, kind]], types?, color, rmsd?}]}
   * view: {R: rows = screen x, screen y, toward the viewer; zoom}   opt: {interactive, pick: {a, b}} */
  function scene(mol, L, view, opt = {}) {
    const sty = Object.assign({}, STYLE, L.style || {});
    const pre = PRESETS[sty.preset] || PRESETS.glossy, pal = Object.assign({}, PALETTES[sty.palette] || PALETTES.soft, sty.atomColors || {});   // style.atomColors: per-element overrides
    const ink0 = Object.assign({}, INKS[sty.ink] || INKS.colour, sty.inkColors || {});
    const ink = ink0;                       // axes may shadow this with their own (ax.ink)
    const K = sty.labelScale || 1, eq = ink.eq, V = sty.values !== false;   // V: labels carry their values, or only the symbols
    const el = mol.el, X = mol.X, n = el.length;
    const axes = L.axes || [], dp = L.dipole && L.dipole.res, fade = L.fade || new Set();
    const dp2 = L.dipole2 && L.dipole2.res;
    const on = (ax, f) => (f === 'axis' ? (ax.show || {}).axis !== false : !!(ax.show || {})[f]);
    const c = [0, 1, 2].map((k) => X.reduce((s, p) => s + p[k], 0) / n);

    // ---- everything the fit has to hold ----
    const ends = [];
    axes.forEach((ax) => {
      wedgeScale = sty.thetaBold ? 0.95 : 0.58;
      const g = ax.g = axisGeometry(ax.res);
      if (on(ax, 'axis') || on(ax, 'L')) ends.push(g.tip);
      if (on(ax, 'b1')) ends.push(g.b1To);
      if (on(ax, 'b5')) ends.push(g.b5To);
      if (on(ax, 'plane')) ends.push(...g.ref, ...g.patch);
    });
    let muScale = 0;
    [dp, L.dipole2 && L.dipole2.res].forEach((dpx) => {      // both frames share one arrow scale
      if (!dpx) return;
      if (!muScale) {
        muScale = sty.muPerD > 0 ? sty.muPerD : Math.min(0.6, 2.6 / (len(dpx.mu) || 1));
        if (sty.dipoleArrow === 'chemistry') muScale = -muScale;      // the textbook arrow points + -> -, Gaussian's vector - -> +
      }
      ends.push(add(dpx.origin, scl(dpx.mu, muScale)));
      dpx.basis.forEach((b) => ends.push(add(dpx.origin, scl(b, 1.7))));
      const DL0 = dpx === dp ? L.dipole : L.dipole2;
      if (DL0 && DL0.atOrigin) ends.push([0, 0, 0], scl(dpx.mu, muScale), [1.55, 0, 0], [0, 1.55, 0], [0, 0, 1.55]);   // the raw arrow and axes have to fit too
    });
    const structs = [{ el, X, bonds: mol.displayBonds || mol.bonds.map(([a, b]) => [a, b, 'pkg']), types: mol.types, primary: true }]
      .concat((L.overlays || []).filter((s) => s && s.X && s.X.length).map((s) => Object.assign({ primary: false }, s)));
    // fadeHide: atoms no axis measures are dropped rather than ghosted, and leave the framing too --
    // on a stacked structure too, each judged by the axis that is measured on it
    const drop = L.fadeHide ? fade : new Set();
    const dropOf = new Map();
    if (L.fadeHide) structs.forEach((t) => {
      if (t.primary || !t.id) return;
      const keep = new Set();
      axes.forEach((ax) => {
        if (ax.onStruct !== t.id) return;
        ax.res.atoms.forEach((i) => keep.add(i));
        keep.add(ax.res.atoms.length ? (L.axes && L.axes[ax.j] ? L.axes[ax.j].a - 1 : -1) : -1);
      });
      if (keep.size) dropOf.set(t.id, new Set(t.el.map((_, i) => i).filter((i) => !keep.has(i))));
    });
    const allX = [].concat(...structs.map((t, k) => {
      const d = k === 0 ? drop : dropOf.get(t.id);
      return d && d.size ? t.X.filter((_, i) => !d.has(i)) : t.X;
    }));

    // ---- projection: the view's rotation, optional perspective, fitted to the page ----
    const Rv = view.R, kp = 0.35 * (sty.perspective || 0);
    const P0 = (p) => { const d = dif(p, c); return [dot(Rv[0], d), dot(Rv[1], d), dot(Rv[2], d)]; };
    const zr = Math.max(1, ...allX.map((p) => Math.abs(P0(p)[2])));
    const proj = (p) => { const q = P0(p), f = 1 / (1 - kp * Math.max(-1.5, Math.min(0.9, q[2] / zr))); return [q[0] * f, q[1] * f, q[2], f]; };
    const insets = axes.filter((ax) => on(ax, 'inset') || on(ax, 'phi')), profiles = axes.filter((ax) => on(ax, 'profile'));
    const IS = 132, PW = 168, PH = 88, insW = insets.length || profiles.length ? IS - 50 : 0;   // they sit in the right margin, partly
    const fit = allX.concat(ends).map(proj), xs = fit.map((q) => q[0]), ys = fit.map((q) => q[1]);
    const cx = (Math.max(...xs) + Math.min(...xs)) / 2, cy = (Math.max(...ys) + Math.min(...ys)) / 2;
    const fitS = Math.min((W - 150 - insW) / (Math.max(...xs) - Math.min(...xs) + 1.2),
      (H - 110) / (Math.max(...ys) - Math.min(...ys) + 1.2));
    const s = (view.zoom || 1) * (view.perA > 0 ? view.perA : fitS);   // perA: hold one scale across panels instead of fitting each
    L.fitScale = fitS;                                                  // what fitting this molecule would have given
    const X0 = (W - insW) / 2, Y0 = H / 2 - 12;
    const pr = (p) => { const q = proj(p); return [X0 + (q[0] - cx) * s, Y0 - (q[1] - cy) * s, q[2], q[3]]; };
    structs.forEach((t) => { t.P = t.X.map(pr); });
    const P = structs[0].P;

    // ---- colours, hydrogens, sizes ----
    const qmap = L.map && mol.charges && mol.charges[L.map.type];
    const qmax = qmap ? (sty.qMax > 0 ? sty.qMax : Math.max(...qmap.map(Math.abs))) : 0;
    const elCol = (e) => pal[e] || pal.fallback || ballColor(e);
    const zAll = [].concat(...structs.map((t) => t.P.map((q) => q[2]))), zmin = Math.min(...zAll), zmax = Math.max(...zAll);
    const fogged = (col, z) => { const k = Math.round(sty.fog * (zmax - z) / ((zmax - zmin) || 1) * 5) / 5 * 0.72; return k > 0 ? shade(col, k) : col; };
    const shades = new Set(), markers = new Set();
    const fillRef = (col) => { if (pre.shade === 'flat' || pre.shade === 'toon') return col; shades.add(col); return `url(#sg${pre.shade}${col.slice(1)})`; };
    // never hide an atom a drawn number stands on: B5's atom, L's (all tied), B1's tangent atoms, the dipole's frame atoms, charge sites
    const keepH = new Set([...(L.sites || []).map((q) => q.i), ...axes.flatMap((ax) => [ax.res.b5, ...ax.g.lAtoms, ...ax.g.touch]),
      ...((L.dipole && L.dipole.keep) || []).map((k) => k - 1)]);
    structs.forEach((t) => {                // hydrogens: all, polar only, or none -- never one a layer is drawn on
      t.hidden = new Set();
      if (sty.hydrogens === 'all') return;
      const nb = t.el.map(() => []);
      t.bonds.forEach(([a, b]) => { nb[a - 1].push(t.el[b - 1]); nb[b - 1].push(t.el[a - 1]); });
      t.el.forEach((e, i) => { if (e === 'H' && (sty.hydrogens === 'none' || nb[i].every((x) => x === 'C')) && !(t.primary && keepH.has(i))) t.hidden.add(i); });
    });
    drop.forEach((i) => structs[0].hidden.add(i));
    structs.forEach((t) => { const d = dropOf.get(t.id); if (d) d.forEach((i) => t.hidden.add(i)); });
    const hidden = structs[0].hidden;
    const CPKT = (root.Sterimol && root.Sterimol.CPK) || {};
    const rAng = (e, ty) => (pre.atom === 'vdw' ? 0.92 * (VDW[e] || 1.8) : pre.atom === 'bond' ? 0.17 * sty.bondScale : pre.atom === 'none' ? 0
      : pre.atom === 'small' ? 0.13 * (CPKT[ty] || CPKT[e] || 1.7) + 0.04 : 0.2 + 0.3 * covalent(e));
    const rbOf = (t, i) => s * sty.atomScale * rAng(t.el[i], t.types && t.types[i]) * t.P[i][3];
    const rb = (i) => rbOf(structs[0], i);
    const bw0 = 0.26 * s * sty.bondScale * (pre.bondW || 1);
    const ol = 0.8 * sty.outline * (pre.outline || 1);
    const olc = pre.olc || (pre.shade === 'flat' || pre.shade === 'toon' ? '#16191E' : '#3B424C');
    const grey = pre.bondGrey || '#A7AFBA';

    // ---- the structures, far to near, all in one depth order ----
    const items = [], hits = [], lab = new Placer(W, H, K, L.labShift || {}, sty.leaders !== false), pick = opt.pick || {};
    structs.forEach((t) => {
      const op = t.primary ? 1 : sty.overlayOpacity, opa = op < 1 ? ` opacity="${f1(op * 100) / 100}"` : '';
      const tint = !t.primary && sty.overlay === 'tint';
      const baseCol = (i) => (t.primary ? (qmap ? divergent(qmap[i], qmax, sty.qRamp) : elCol(t.el[i])) : tint ? t.color : elCol(t.el[i]));
      const colOf = (i) => fogged(baseCol(i), t.P[i][2]);
      const ghost = (i) => t.primary && fade.has(i);
      if (pre.bonds) t.bonds.forEach(([a, b, kind]) => {
        const i = a - 1, j = b - 1;
        if (t.hidden.has(i) || t.hidden.has(j)) return;
        const p = t.P[i], q = t.P[j], m = [(p[0] + q[0]) / 2, (p[1] + q[1]) / 2], bw = bw0 * (p[3] + q[3]) / 2;
        const dash = kind === 'extra' || (sty.metalDash && (METAL.has(t.el[i]) || METAL.has(t.el[j])))
          ? ` stroke-dasharray="${f1(Math.max(2, bw * 0.9))} ${f1(Math.max(2, bw * 0.75))}"` : '';
        const one = (sty.bondColor === 'grey' || sty.bondColor === 'black') && !tint, mono = sty.bondColor === 'black' ? '#1D2025' : grey;
        const ci = one ? mono : shade(colOf(i), pre.wire ? 0 : 0.15), cj = one ? mono : shade(colOf(j), pre.wire ? 0 : 0.15);
        const seg = (x0, y0, x1, y1, stl) => `<line x1="${f1(x0)}" y1="${f1(y0)}" x2="${f1(x1)}" y2="${f1(y1)}" style="${stl}"${dash}/>`;
        const w = (wd, col) => `stroke:${col};stroke-width:${f1(wd)}px;stroke-linecap:round`;
        let g = `<g class="bond${ghost(i) || ghost(j) ? ' ghost' : ''}${kind === 'extra' ? ' xb' : ''}"${opa}>`;
        if (!pre.wire && pre.bondOl !== false && ol > 0) g += seg(p[0], p[1], q[0], q[1], w(bw + 2 * ol, olc));
        g += one ? seg(p[0], p[1], q[0], q[1], w(bw, ci)) : seg(p[0], p[1], m[0], m[1], w(bw, ci)) + seg(m[0], m[1], q[0], q[1], w(bw, cj));
        if (pre.tube) g += seg(p[0], p[1], q[0], q[1], `stroke:#FFFFFF;stroke-opacity:.42;stroke-width:${f1(bw * 0.28)}px;stroke-linecap:round`);
        items.push([(p[2] + q[2]) / 2 - 1e-3, g + '</g>']);
        if (t.primary && !ghost(i) && !ghost(j)) lab.along(p, q, bw / 2 + 1.5, 8);
      });
      t.el.forEach((e, i) => {
        if (t.hidden.has(i)) return;
        const [x, y, z] = t.P[i], r = rbOf(t, i), gh = ghost(i);
        let g = '';
        if (r > 0) {
          g += `<circle class="ball${gh ? ' ghost' : ''}" cx="${f1(x)}" cy="${f1(y)}" r="${f1(r)}" fill="${fillRef(colOf(i))}" style="stroke:${olc};stroke-width:${f1(ol * (pre.atom === 'vdw' ? 0.7 : 1))}px"${opa}/>`;
          if (pre.shade === 'toon') {
            const hx = x - 0.34 * r, hy = y - 0.38 * r;
            g += `<ellipse${gh ? ' class="ghost"' : ''} cx="${f1(hx)}" cy="${f1(hy)}" rx="${f1(0.26 * r)}" ry="${f1(0.17 * r)}" transform="rotate(-35 ${f1(hx)} ${f1(hy)})" fill="#FFFFFF" opacity="${f1(0.85 * op * 100) / 100}"/>`;
          }
        }
        if (t.primary && (i === pick.a || i === pick.b)) g += `<circle class="pickring ${i === pick.a ? 'pa' : 'pb'}" cx="${f1(x)}" cy="${f1(y)}" r="${f1(r + 4.5)}"/>`;
        if (t.primary && pick.hi && pick.hi.has(i)) g += `<circle class="pickring hi" cx="${f1(x)}" cy="${f1(y)}" r="${f1(r + 6.5)}" stroke="${pick.hiColor || '#3450A8'}"/>`;
        if (opt.interactive) hits.push([z, `<circle class="hit" fill="transparent" data-i="${i}"${t.primary ? '' : ` data-s="${esc(t.id)}"`} cx="${f1(x)}" cy="${f1(y)}" r="${f1(Math.max(r + 2, 7))}"><title>${t.primary ? '' : esc(t.name) + ' · '}${e}${i + 1}</title></circle>`]);
        items.push([z, g]);
        if (!gh) lab.circles.push([x, y, Math.max(r, 3) + 1.5]);
      });
    });

    // ---- the van der Waals surface, as its own layer ----
    // Opaque circles inside one translucent group: overlapping spheres then read as a single
    // surface instead of darkening where they meet. CPK radii, the same ones Sterimol measures with.
    const vdwLayer = (() => {
      if (!sty.vdw || sty.vdw === 'off') return '';
      const t = structs[0], op = Math.max(0.05, Math.min(1, sty.vdwOpacity == null ? 0.32 : sty.vdwOpacity));
      const sc = sty.vdwScale || 1;
      const skin = (i) => fogged(qmap ? divergent(qmap[i], qmax, sty.qRamp) : elCol(t.el[i]), 0);
      if (sty.surfKind === 'density' && typeof Promol !== 'undefined') {
        // the promolecular isosurface: meshed once per molecule and isovalue, then projected
        const skip = new Set([...t.hidden, ...fade]);
        const m = Promol.cached(t.el, t.X, { iso: sty.surfIso || Promol.ISO0, step: sty.surfStep || 0.4, skip });
        if (!m.quads.length) return '';
        const Q = m.verts.map(pr);
        const nv = m.norms.map((n) => [dot(Rv[0], n), -dot(Rv[1], n), dot(Rv[2], n)]);   // screen axes: x right, y down, z at the viewer
        const LX = -0.42, LY = -0.57, LZ = 0.71;
        const faces = [];
        m.quads.forEach((q) => {
          let nx = 0, ny = 0, nz = 0;
          for (const k of q) { nx += nv[k][0]; ny += nv[k][1]; nz += nv[k][2]; }
          if (nz <= 0) return;                            // the far side is hidden by the near side
          const nn = Math.hypot(nx, ny, nz) || 1;
          const lit = Math.max(-1, Math.min(1, (nx * LX + ny * LY + nz * LZ) / nn));
          const z = (Q[q[0]][2] + Q[q[1]][2] + Q[q[2]][2] + Q[q[3]][2]) / 4;
          // a density has no element identity: one neutral skin, unless a charge map is on
          const base = qmap ? divergent(qmap[m.near[q[0]]], qmax, sty.qRamp) : '#DFE3E8';
          const col = shade(fogged(base, z), 0.3 * lit - 0.14);
          faces.push([z, `<path d="M${q.map((k) => f1(Q[k][0]) + ',' + f1(Q[k][1])).join(' L')} Z" fill="${col}" stroke="${col}"/>`]);
        });
        return `<g class="vdw density" opacity="${f1(op * 100) / 100}">`
          + faces.sort((a, b) => a[0] - b[0]).map(([, d]) => d).join('') + '</g>';
      }
      const balls = t.el.map((e, i) => [t.P[i], i, e])
        .filter(([, i]) => !t.hidden.has(i) && !fade.has(i))
        .sort((p, q) => p[0][2] - q[0][2])
        .map(([P, i, e]) => {
          const r = s * sc * (VDW[e] || 1.8) * P[3];
          const col = fogged(qmap ? divergent(qmap[i], qmax, sty.qRamp) : elCol(e), P[2]);
          return `<circle cx="${f1(P[0])}" cy="${f1(P[1])}" r="${f1(r)}" fill="${fillRef(col)}"/>`;
        }).join('');
      return `<g class="vdw" opacity="${f1(op * 100) / 100}">${balls}</g>`;
    })();

    // ---- pass 1: every construction line (and the obstacles it leaves for labels) ----
    const under = [], over = [], jobs = [], legend = [];
    // what a construction line can disappear behind: the measured structure's drawn atoms
    const occ = [];
    if (sty.hiddenLines !== false) structs[0].P.forEach((q, i) => { if (!hidden.has(i) && !fade.has(i)) occ.push([q[0], q[1], rbOf(structs[0], i), q[2]]); });
    const behind = (p) => occ.some(([x, y, rr, z]) => {          // is this projected point behind a sphere's front surface?
      const d = Math.hypot(p[0] - x, p[1] - y);
      if (d >= rr) return false;
      const rw = rr / s, dw = d / s;
      return p[2] + 0.02 < z + Math.sqrt(Math.max(0, rw * rw - dw * dw));
    });
    function runsOf(a, b) {                                      // split a 3D segment into visible and hidden runs
      const N = 32, pts = [];
      for (let k = 0; k <= N; k++) { const t = k / N; pts.push(pr([a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t])); }
      const hid = pts.map(behind), out = [];
      let k0 = 0;
      for (let k = 1; k <= N; k++) if (hid[k] !== hid[k0] || k === N) { out.push([pts[k0], pts[k], hid[k0]]); k0 = k; }
      return out;
    }
    const mk = (col) => { markers.add(col); return `url(#mk${col.slice(1)})`; };
    const seg1 = (p, q, cls, col, extra = '') => `<line class="${cls}" x1="${f1(p[0])}" y1="${f1(p[1])}" x2="${f1(q[0])}" y2="${f1(q[1])}" stroke="${col}"${extra}/>`;
    const ln = (a, b, cls, col, extra = '') => {
      if (sty.hiddenLines === false || !occ.length) return seg1(pr(a), pr(b), cls, col, extra);
      const rs = runsOf(a, b);                                   // faint where it passes behind an atom
      return rs.map(([p, q, h], k) => seg1(p, q, cls + (h ? ' behind' : ''), col, k === rs.length - 1 ? extra : '')).join('');
    };
    const path = (pts, cls, attrs, close = true) => `<path class="${cls}" d="M${pts.map((p) => { const q = pr(p); return f1(q[0]) + ',' + f1(q[1]); }).join(' L')}${close ? ' Z' : ''}"${attrs}/>`;
    const fill = (col) => ` fill="${col}"`;

    let planeSaid = false;
    axes.forEach((ax) => {
      // ax.ink recolours one axis's whole construction, so two structures can be told apart
      const ink = ax.ink ? Object.assign({}, ink0, ax.ink) : ink0;
      const g = ax.g, r = ax.res, col = ink.axis || ax.color, pb = pr(g.base), ptip = pr(g.tip);
      if (ax.preview) {                   // one of the scan's other rotations, replayed
        const gp = axisGeometry(r, ax.preview);
        under.push(path(gp.patch, 'sx-patch trial', ` fill="${ink.b1}" stroke="${ink.b1}"`));
        over.push(ln(gp.b1From, gp.b1To, 'sx-b1 trial', ink.b1), ln(gp.plane[0], gp.plane[1], 'sx-trace trial', ink.b1));
        const pt = pr(gp.b1To);
        jobs.push([0, () => lab.put([[pt[0] + 10, pt[1] - 6, 'start'], [pt[0] - 10, pt[1] - 6, 'end'], [pt[0], pt[1] + 16, 'middle']],
          `B${sub(1)} ${ax.preview.B1.toFixed(2)} &#8491; at ${ax.preview.deg}&#176;</tspan>`, 'lab', fill(ink.b1), pt)]);
      }
      if (on(ax, 'plane')) {
        under.push(path(g.ref, 'sx-ref', ` fill="${ink.plane}" stroke="${ink.plane}"`), path(g.patch, 'sx-patch', ` fill="${ink.b1}" stroke="${ink.b1}"`));
        g.touch.forEach((i) => { if (!hidden.has(i)) over.push(`<circle class="sx-touch" cx="${f1(P[i][0])}" cy="${f1(P[i][1])}" r="${f1(rb(i) + 3)}" stroke="${ink.b1}"/>`); });
        const pc = g.patch.map(pr), rc = g.ref.map(pr), edge = rc.map((p, i) => [(p[0] + rc[(i + 1) % 4][0]) / 2, (p[1] + rc[(i + 1) % 4][1]) / 2]);
        jobs.push([6, () => lab.put(pc.flatMap(([x, y]) => [[x, y - 5, 'middle'], [x, y + 13, 'middle'], [x - 6, y + 3, 'end'], [x + 6, y + 3, 'start']]),
          `B${sub(1)} plane</tspan>`, 'mark', fill(ink.b1))]);
        if (!planeSaid && sty.planeNote !== false) {   // one annotation however many axes show their plane
          planeSaid = true;
          jobs.push([7, () => lab.put(edge.flatMap(([x, y]) => [[x, y + 15, 'middle'], [x, y - 6, 'middle']]), `&#8741; B${sub(1)}, through the axis</tspan>`, 'mark', fill(ink.plane), edge[0])]);
        }
      }
      // how a measured vector is marked. The graduations sit in molecule space (every 1 A of
      // the real length) but are drawn in screen space, so they stay square to the line and
      // legible whatever the view. Used by the axis, B1, B5 and v alike.
      const marked = (P1, P2, cls, colour, opt) => {
        const o = opt || {}, style2 = o.mark || sty.axisMark || 'arrow', out = [];
        const q1 = pr(P1), q2 = pr(P2);
        const dx = q2[0] - q1[0], dy = q2[1] - q1[1], dn = Math.hypot(dx, dy) || 1;
        const u2 = [dx / dn, dy / dn], n2 = [-u2[1], u2[0]];
        const L3 = len(dif(P2, P1)) || 1, dir = scl(dif(P2, P1), 1 / L3);   // the real length, in A
        const at = (a1) => pr(add(P1, scl(dir, a1)));
        const seg = (z1, z2, c, extra) => `<line class="${c}" x1="${f1(z1[0])}" y1="${f1(z1[1])}" x2="${f1(z2[0])}" y2="${f1(z2[1])}" stroke="${colour}"${extra || ''}/>`;
        const head = o.head === false ? '' : ` marker-end="${mk(colour)}"`;
        if (style2 === 'tape') {
          const w = o.w || 3.4, band = (z1, z2, on2) => {
            const q = [[z1[0] + n2[0] * w, z1[1] + n2[1] * w], [z2[0] + n2[0] * w, z2[1] + n2[1] * w],
              [z2[0] - n2[0] * w, z2[1] - n2[1] * w], [z1[0] - n2[0] * w, z1[1] - n2[1] * w]];
            return `<path class="sx-tape${on2 ? ' on' : ''}" d="M${q.map((z) => f1(z[0]) + ',' + f1(z[1])).join(' L')} Z" fill="${on2 ? colour : 'none'}" stroke="${colour}"/>`;
          };
          let prev = q1, i = 0;
          for (let a1 = 1; a1 <= L3; a1++) { const z = at(a1); out.push(band(prev, z, i++ % 2 === 0)); prev = z; }
          out.push(band(prev, q2, i % 2 === 0));
        } else if (style2 === 'ruler') {
          out.push(seg(q1, q2, cls, head));
          for (let a1 = 1; a1 <= L3 - 0.05; a1++) {
            const z = at(a1), h = a1 % 5 === 0 ? 7.5 : 4.5;
            out.push(seg([z[0] - n2[0] * h, z[1] - n2[1] * h], [z[0] + n2[0] * h, z[1] + n2[1] * h], 'sx-grad'));
            if (V && o.numbers !== false && a1 % 2 === 0) out.push(
              `<text class="mark" x="${f1(z[0] + n2[0] * 13)}" y="${f1(z[1] + n2[1] * 13 + 3)}" text-anchor="middle" font-size="${f1(FS.mark * K * 0.82)}" fill="${colour}">${a1}</text>`);
          }
        } else if (style2 === 'dim') {
          const e = o.w ? o.w * 2 : 7;
          out.push(seg([q1[0] - n2[0] * e, q1[1] - n2[1] * e], [q1[0] + n2[0] * e, q1[1] + n2[1] * e], 'sx-dimend'),
            seg([q2[0] - n2[0] * e, q2[1] - n2[1] * e], [q2[0] + n2[0] * e, q2[1] + n2[1] * e], 'sx-dimend'),
            seg(q1, q2, cls, ` marker-start="${mk(colour)}"${head}`));
        } else if (style2 === 'pin') {
          const w = o.w || 3.6;
          out.push(`<path class="sx-pin" d="M${f1(q1[0] + n2[0] * w)},${f1(q1[1] + n2[1] * w)} L${f1(q2[0])},${f1(q2[1])} L${f1(q1[0] - n2[0] * w)},${f1(q1[1] - n2[1] * w)} Z" fill="${colour}"/>`,
            `<circle class="sx-pindot" cx="${f1(q1[0])}" cy="${f1(q1[1])}" r="${f1(w * 0.9)}" fill="${colour}"/>`);
        } else {
          out.push(ln(P1, P2, cls, colour, head));
        }
        return out;
      };
      if (on(ax, 'axis')) { over.push(...marked(g.base, g.tip, 'sx-axis', col)); lab.along(pb, ptip); }
      if (on(ax, 'L')) {
        jobs.push([5, () => lab.put([[ptip[0], ptip[1] - 12, 'middle'], [ptip[0] + 10, ptip[1] + 4, 'start'], [ptip[0] - 10, ptip[1] + 4, 'end'], [ptip[0], ptip[1] + 18, 'middle']],
          V ? `L${eq}${r.L.toFixed(2)} &#8491;` : 'L', 'lab', fill(ink.text), ptip)]);
      }
      if (on(ax, 'phi')) {                  // the azimuth of B5 from the B1 contact, in the plane across the axis
        over.push(path([g.phiAt].concat(g.phiRing), 'sx-wedge', fill(ink.theta)), path(g.phiRing, 'sx-arc', ` stroke="${ink.theta}"`, false),
          ln(g.phiAt, g.phiN, 'sx-shadow', ink.b1), ln(g.phiAt, g.phiQ, 'sx-shadow', ink.b5));
        jobs.push([1, () => {
          const pm = pr(g.phiRing[16]), pa = pr(g.phiAt), dx = pm[0] - pa[0], dy = pm[1] - pa[1], dn = Math.hypot(dx, dy) || 1;
          const cands = [12, 22].flatMap((d) => [[pm[0] + dx / dn * d, pm[1] + dy / dn * d + 4, dx < -0.3 * dn ? 'end' : dx > 0.3 * dn ? 'start' : 'middle']]);
          return lab.put(cands.concat([[pm[0], pm[1] - 10, 'middle'], [pm[0], pm[1] + 16, 'middle']]),
            V ? `&#966;${eq}${r.phi.toFixed(1)}&#176;` : '&#966;', 'th', fill(ink.theta), pm);
        }]);
      }
      if (on(ax, 'theta')) {                // the angle at the base between v and its shadow in the plane through the axis
        const pS = pr(g.S), pv = pr(g.v), full = ax.show.theta !== 'wedge' && ax.show.drop !== false;
        over.push(path([g.base].concat(g.ring), 'sx-wedge', fill(ink.theta)), path(g.ring, 'sx-arc', ` stroke="${ink.theta}"`, false),
          ln(g.base, g.S, 'sx-shadow', ink.theta));
        if (full) over.push(ln(g.S, g.v, 'sx-perp', ink.theta), path(g.right, 'sx-right', ` stroke="${ink.right}"`, false));
        over.push(...marked(g.base, g.v, ink.draft ? 'sx-v' : 'sx-v dash', ink.b5, { w: 2.8 }),
          `<circle class="sx-ring" cx="${f1(pv[0])}" cy="${f1(pv[1])}" r="${f1(rb(r.b5) + 3)}" stroke="${ink.b5}"/>`,
          `<circle cx="${f1(pb[0])}" cy="${f1(pb[1])}" r="2.4" fill="${ink.theta}"/>`);
        lab.along(pb, pS); lab.along(pS, pv); lab.along(pb, pv); lab.circles.push([pS[0], pS[1], 9]);
        jobs.push([0, () => {                // theta gets the best spot: across the shadow from v, by the arc
          const cands = [];
          [0.8, 1.05, 0.55, 1.3].forEach((f) => {
            const pa = pr(add(g.base, scl(g.u1, f * g.rr0))), ux = pa[0] - pb[0], uy = pa[1] - pb[1], un = Math.hypot(ux, uy) || 1;
            let px = -uy / un, py = ux / un;
            if (px * (pv[0] - pb[0]) + py * (pv[1] - pb[1]) > 0) { px = -px; py = -py; }
            [14, 24].forEach((d) => cands.push([pa[0] + d * px, pa[1] + d * py + 4, px < -0.3 ? 'end' : px > 0.3 ? 'start' : 'middle']));
          });
          return lab.put(cands, V ? `&#952;${eq}${r.B1_B5_angle.toFixed(1)}&#176;` : '&#952;', 'th', fill(ink.theta), pr(g.ring[Math.floor(g.ring.length / 2)]));
        }]);
        if (full) jobs.push([4, () => {
          const pm = pr(add(g.base, scl(dif(g.v, g.base), 0.62))), nx = -(pm[1] - pb[1]), ny = pm[0] - pb[0], nn = Math.hypot(nx, ny) || 1;
          return lab.put([1, -1, 1.7, -1.7].map((sd) => [pm[0] + sd * 12 * nx / nn, pm[1] + sd * 12 * ny / nn + 4, 'middle']), 'v', 'vl', fill(ink.b5));
        }]);
      }
      if (on(ax, 'b1')) {
        over.push(...marked(g.b1From, g.b1To, 'sx-b1', ink.b1, { w: 2.8, numbers: false }));
        if (!on(ax, 'plane')) { over.push(ln(g.plane[0], g.plane[1], 'sx-trace', ink.b1)); lab.along(pr(g.plane[0]), pr(g.plane[1]), 3.5); }
        const p1 = pr(g.b1From), p2 = pr(g.b1To), pm = [(p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2];
        lab.along(p1, p2);
        jobs.push([2, () => lab.put([[p2[0] + 8, p2[1] - 6, 'start'], [p2[0] - 8, p2[1] - 6, 'end'], [p2[0], p2[1] + 16, 'middle'], [p2[0], p2[1] - 12, 'middle'],
          [pm[0], pm[1] - 9, 'middle'], [pm[0], pm[1] + 15, 'middle'], [p2[0] + 14, p2[1] + 16, 'start'], [p2[0] - 14, p2[1] + 16, 'end'],
          [p2[0] + 14, p2[1] - 16, 'start'], [p2[0] - 14, p2[1] - 16, 'end'], [p2[0], p2[1] + 30, 'middle'], [p2[0], p2[1] - 26, 'middle']],
          `B${sub(1)}${V ? eq + r.B1.toFixed(2) + ' &#8491;' : ''}</tspan>`, 'lab', fill(ink.b1), p2)]);
      }
      if (on(ax, 'loc') || on(ax, 'b5')) {  // loc B5: the level along the axis where the widest point sits (B5 always shows it)
        const p = pr(g.loc), pb = pr(g.base), pc = pr(c);
        if (on(ax, 'theta')) over.push(ln(g.S, g.loc, 'sx-locline', ink.draft ? ink.loc : ink.locText));
        const seg = (a, b, cls, col) => `<line class="${cls}" x1="${f1(a[0])}" y1="${f1(a[1])}" x2="${f1(b[0])}" y2="${f1(b[1])}" stroke="${col}"/>`;
        const ux = p[0] - pb[0], uy = p[1] - pb[1], un = Math.hypot(ux, uy) || 1;
        const ax1 = [ux / un, uy / un];
        let nx = -ax1[1], ny = ax1[0];                       // across the axis, on the side away from the molecule
        if (nx * (pc[0] - p[0]) + ny * (pc[1] - p[1]) > 0) { nx = -nx; ny = -ny; }
        const mark = sty.locMark || 'tick';
        let M = [p[0] + nx * 12, p[1] + ny * 12];
        if (mark === 'run') {                                // the run of axis up to it, ending in a collar
          over.unshift(ln(g.base, g.loc, 'sx-runhalo', '#FFFFFF'), ln(g.base, g.loc, 'sx-run', ink.loc));
          const e1 = unit(cross(g.axisDir, Math.abs(g.axisDir[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0])), e2 = cross(g.axisDir, e1);
          const ring = Array.from({ length: 33 }, (_, k) => { const t = 2 * Math.PI * k / 32; return add(g.loc, add(scl(e1, 0.5 * Math.cos(t)), scl(e2, 0.5 * Math.sin(t)))); });
          over.push(path(ring, 'sx-collarhalo', ' stroke="#FFFFFF"', false), path(ring, 'sx-collar', ` stroke="${ink.loc}"`, false));
          lab.along(pb, p, 4);
          M = [(pb[0] + p[0]) / 2, (pb[1] + p[1]) / 2 + 0];
        } else if (mark === 'caret') {                       // a pointer at the axis from the label's side
          const apex = [p[0] + nx * 4, p[1] + ny * 4], w = 5.5, h = 13;
          const b1 = [apex[0] + nx * h - ax1[0] * w, apex[1] + ny * h - ax1[1] * w];
          const b2 = [apex[0] + nx * h + ax1[0] * w, apex[1] + ny * h + ax1[1] * w];
          over.push(`<path class="sx-caret" d="M${f1(apex[0])},${f1(apex[1])} L${f1(b1[0])},${f1(b1[1])} L${f1(b2[0])},${f1(b2[1])} Z" fill="${ink.loc}"/>`);
          M = [p[0] + nx * 18, p[1] + ny * 18];
        } else if (mark === 'cross') {                       // the X of the published panels
          const k = 6;
          over.push(`<path class="sx-cross" d="M${f1(p[0] - k)},${f1(p[1] - k)} L${f1(p[0] + k)},${f1(p[1] + k)} M${f1(p[0] - k)},${f1(p[1] + k)} L${f1(p[0] + k)},${f1(p[1] - k)}" stroke="${ink.loc}"/>`);
        } else {                                             // a tick across the axis
          const t0 = [p[0] - nx * 7.5, p[1] - ny * 7.5], t1 = [p[0] + nx * 7.5, p[1] + ny * 7.5];
          over.push(seg(t0, t1, 'sx-tickhalo', '#FFFFFF'), seg(t0, t1, 'sx-loctick', ink.loc));
        }
        lab.circles.push([p[0], p[1], 9]);
        jobs.push([1, () => lab.put([1, -1, 1.7, -1.7].map((sd) => [M[0] + nx * sd * 4, M[1] + ny * sd * 4 + 4, 'middle'])
          .concat([[14, -8], [-14, -8], [14, 16], [-14, 16], [26, -20], [-26, -20], [26, 28], [-26, 28]].map(([dx, dy]) => [p[0] + dx, p[1] + dy, dx > 0 ? 'start' : 'end'])),
          `loc B${sub(ax.tag ? '5,' + ax.tag : 5)}${V ? eq + r.loc_B5.toFixed(2) + ' &#8491;' : ''}</tspan>`, 'lab', fill(ink.locText), p)]);
      }
      if (on(ax, 'b5')) {                   // B5: from loc B5 on the axis, across it, to the far side of the atom
        const p = pr(g.loc), e = pr(g.b5To), ex = e[0] - p[0], ey = e[1] - p[1], en = Math.hypot(ex, ey) || 1;
        if (ink.draft) {
          over.push(ln(g.loc, g.b5To, 'sx-b5 dim', ink.b5),
            `<line class="sx-b5 dim" x1="${f1(e[0] - 5 * ey / en)}" y1="${f1(e[1] + 5 * ex / en)}" x2="${f1(e[0] + 5 * ey / en)}" y2="${f1(e[1] - 5 * ex / en)}" stroke="${ink.b5}"/>`,
            path(g.right5, 'sx-right', ` stroke="${ink.right}"`, false));
        } else over.push(...marked(g.loc, g.b5To, 'sx-b5', ink.b5, { w: 2.8, numbers: false }));
        lab.along(p, e);
        jobs.push([3, () => lab.put([[e[0] + 10 * ex / en, e[1] + 10 * ey / en + 4, ex > 0 ? 'start' : 'end'], [e[0], e[1] - 10, 'middle'], [e[0], e[1] + 16, 'middle'],
          [(p[0] + e[0]) / 2, (p[1] + e[1]) / 2 - 8, 'middle'], [(p[0] + e[0]) / 2, (p[1] + e[1]) / 2 + 15, 'middle'],
          [e[0] - 26, e[1] - 14, 'middle'], [e[0] + 26, e[1] - 14, 'middle'], [e[0] - 26, e[1] + 20, 'middle'], [e[0] + 26, e[1] + 20, 'middle']],
          `B${sub(5)}${V ? eq + r.B5.toFixed(2) + ' &#8491;' : ''}</tspan>`, 'lab', fill(ink.b5), e)]);
      }
      legend.push({ col, text: `${esc(ax.label || '')}  B${sub(1)} ${r.B1.toFixed(2)}</tspan> · B${sub(5)} ${r.B5.toFixed(2)}</tspan> · L ${r.L.toFixed(2)} · loc B${sub(5)} ${r.loc_B5.toFixed(2)}</tspan> · &#952; ${r.B1_B5_angle.toFixed(1)}&#176;${on(ax, 'phi') ? ` &#183; &#966; ${r.phi.toFixed(1)}&#176;` : ''}` });
    });
    // ---------- the same parameter, read from a position: cone angle and %Vbur ----------
    // Both are computed here (like the promolecular surface) so a session only carries the
    // position, not the numbers. Vbur.* is the gated port of metal_complex / morfeus.
    const M0 = structs[0];
    (L.cones || []).forEach((cn) => {
      if (typeof Vbur === 'undefined') return;
      const r = Vbur.coneAngle(M0.el, M0.X, cn.apex, { atoms: cn.atoms });
      if (!r) return;
      const col = cn.color || ink.theta, A = r.apex, u = r.axis, al = r.half * Math.PI / 180;
      const reach = Math.max(...r.atoms.map((i) => Math.hypot(M0.X[i][0] - A[0], M0.X[i][1] - A[1], M0.X[i][2] - A[2])));
      const s0 = (cn.reach && cn.reach <= 2 ? cn.reach * reach : cn.reach) || reach;   // slant length drawn, or a fraction of the reach
      const e1 = unit(Math.abs(u[0]) < 0.9 ? [0, -u[2], u[1]] : [-u[2], 0, u[0]]);
      const e2 = [u[1] * e1[2] - u[2] * e1[1], u[2] * e1[0] - u[0] * e1[2], u[0] * e1[1] - u[1] * e1[0]];
      const C = add(A, scl(u, s0 * Math.cos(al))), rho = s0 * Math.sin(al);
      const rim = Array.from({ length: 96 }, (_, k) => {
        const t = 2 * Math.PI * k / 96;
        return add(C, add(scl(e1, rho * Math.cos(t)), scl(e2, rho * Math.sin(t))));
      });
      const pA = pr(A), pRim = rim.map(pr);
      // the two silhouette points: furthest either side of the projected axis
      const pc = pr(C), ax2 = [pc[0] - pA[0], pc[1] - pA[1]], an = Math.hypot(ax2[0], ax2[1]) || 1;
      const side = pRim.map((q) => ((q[0] - pA[0]) * ax2[1] - (q[1] - pA[1]) * ax2[0]) / an);
      const kL = side.indexOf(Math.min(...side)), kR = side.indexOf(Math.max(...side));
      // both arcs join the two silhouette points; the cone's outline is the one that bulges
      // away from the apex on screen (the other one collapses onto it when seen side-on)
      const walk = (from, to) => { const out = []; for (let k = from; ; k = (k + 1) % pRim.length) { out.push(pRim[k]); if (k === to) break; } return out; };
      const far = (arr) => { const q = arr[Math.floor(arr.length / 2)]; return Math.hypot(q[0] - pA[0], q[1] - pA[1]); };
      const a1 = walk(kL, kR), a2 = walk(kR, kL);
      const arc = far(a1) >= far(a2) ? a1 : a2;
      over.push(`<path class="cone-face" d="M${f1(pA[0])},${f1(pA[1])} ${arc.map((q) => 'L' + f1(q[0]) + ',' + f1(q[1])).join(' ')} Z" fill="${col}"/>`,
        `<path class="cone-rim" d="M${pRim.map((q) => f1(q[0]) + ',' + f1(q[1])).join(' L')} Z" fill="none" stroke="${col}"/>`,
        ln(A, add(A, scl(u, s0)), 'cone-axis', col),
        ln(A, rim[kL], 'cone-edge', col), ln(A, rim[kR], 'cone-edge', col),
        `<circle cx="${f1(pA[0])}" cy="${f1(pA[1])}" r="3" fill="${col}"/>`);
      if (sty.featureLabels !== false) jobs.push([2, () => {
        const q = pr(rim[kL]);
        return lab.put([[q[0] - 10, q[1] - 6, 'end'], [q[0] + 10, q[1] - 6, 'start'], [q[0], q[1] - 14, 'middle'], [q[0], q[1] + 16, 'middle']],
          V ? `&#952;${eq}${r.angle.toFixed(1)}&#176;` : '&#952;', 'th', fill(col), q);
      }]);
      legend.push({ col, text: `cone at ${esc(M0.el[cn.apex - 1])}${cn.apex}  &#952; ${r.angle.toFixed(1)}&#176;` });
    });
    (L.vburs || []).forEach((bv) => {
      if (typeof Vbur === 'undefined') return;
      const c = M0.X[bv.atom - 1], R = bv.radius || 3.5;
      const r = Vbur.buriedVolume(M0.el, M0.X, c, R);
      const col = bv.color || ink.b5, pc = pr(c), rad = R * s * (pc[3] || 1);
      over.push(`<circle class="bv-ball" cx="${f1(pc[0])}" cy="${f1(pc[1])}" r="${f1(rad)}" fill="${col}" stroke="${col}"/>`);
      if (bv.dots !== false) {                            // the occupied grid points, outer shell
        const gp = Vbur.grid().pts, step = bv.step || 6;
        const pts = [];
        for (let j = 11 * 4000; j < 12 * 4000; j += step) {
          if (!r.buried[j]) continue;
          const q = pr([c[0] + R * gp[3 * j], c[1] + R * gp[3 * j + 1], c[2] + R * gp[3 * j + 2]]);
          if (q[2] > 0) pts.push(q);                      // front hemisphere only
        }
        over.push(`<g class="bv-dots" fill="${col}">${pts.map((q) => `<circle cx="${f1(q[0])}" cy="${f1(q[1])}" r="1.5"/>`).join('')}</g>`);
      }
      if (sty.featureLabels !== false) jobs.push([2, () => lab.put(
        [[pc[0], pc[1] - rad - 8, 'middle'], [pc[0] + rad + 8, pc[1], 'start'], [pc[0] - rad - 8, pc[1], 'end'], [pc[0], pc[1] + rad + 14, 'middle']],
        V ? `%V<tspan baseline-shift="sub" font-size="72%">bur</tspan>${eq}${r.percent.toFixed(1)}` : `%V<tspan baseline-shift="sub" font-size="72%">bur</tspan>`,
        'lab', fill(col), pc)]);
      legend.push({ col, text: `%V<tspan baseline-shift="sub" font-size="72%">bur</tspan> at ${esc(M0.el[bv.atom - 1])}${bv.atom}, r ${R} &#8491;  ${r.percent.toFixed(1)}%` });
    });
    const sideBoxes = insets.map((ax) => ({ ax, kind: 'inset', h: IS })).concat(profiles.map((ax) => ({ ax, kind: 'profile', h: PH })));
    let sideY = 34;
    sideBoxes.forEach((b) => { b.y = sideY; sideY += b.h + 30; lab.boxes.push([W - Math.max(IS, PW) - 16, b.y - 16, W, b.y + b.h + 4]); });

    // one dipole frame, or two when a stacked structure carries its own
    let muSaid = false;                     // the arrow scale is the same for both: say it once
    [[L.dipole, dp], [L.dipole2, dp2]].forEach(([DL, dpx]) => {
      if (!dpx) return;
      const ink = DL.colors ? Object.assign({}, ink0, DL.colors) : ink0;
      const sh = DL.show || {}, O = dpx.origin, po = pr(O);
      const axs = [['v', dpx.basis[0], dpx.dipole_x], ['u', dpx.basis[1], dpx.dipole_y], ['w', dpx.basis[2], dpx.dipole_z]];
      if ((DL.originSet || []).length > 1) DL.originSet.forEach((a) => over.push(ln(O, X[a - 1], 'dp-spoke', ink.origin)));
      if (sh.frame !== false) axs.forEach(([k, b]) => {
        const tip = add(O, scl(b, 1.7)), p = pr(add(O, scl(b, 1.95)));
        over.push(ln(O, tip, 'dp-ax', ink[k], ` marker-end="${mk(ink[k])}"`),
          `<text class="dpn" x="${f1(p[0])}" y="${f1(p[1] + 5)}" text-anchor="middle" font-size="${f1(FS.dpn * K)}" fill="${ink[k]}">${k}&#770;</text>`);
        lab.along(po, pr(tip)); lab.circles.push([p[0], p[1], 8]);
      });
      const muTip = add(O, scl(dpx.mu, muScale));
      if (sh.mu !== false) lab.along(po, pr(muTip), 5);      // the mu arrow keeps the projection labels off it
      axs.forEach(([k, b, comp]) => {
        if (!sh[k]) return;
        const foot = add(O, scl(b, comp * muScale));
        if (sh.frame === false) {
          const far = add(O, scl(b, Math.sign(comp || 1) * Math.max(1.7, Math.abs(comp * muScale) + 0.5)));
          over.push(ln(O, far, 'dp-guide', ink[k]));
          lab.along(po, pr(far));
          const pn = pr(add(far, scl(b, Math.sign(comp || 1) * 0.3))), pf = pr(far);
          jobs.push([8, () => lab.put([[pn[0], pn[1] + 5, 'middle'], [pf[0] + 12, pf[1] + 4, 'start'], [pf[0] - 12, pf[1] + 4, 'end'],
            [pf[0], pf[1] - 12, 'middle'], [pf[0], pf[1] + 20, 'middle']], `${k}&#770;`, 'dpn', fill(ink[k]))]);
        }
        over.push(ln(O, foot, 'dp-comp', ink[k], ` marker-end="${mk(ink[k])}"`), ln(muTip, foot, 'dp-drop', ink.mu));
        const p = pr(foot), mid = pr(add(O, scl(b, comp * muScale * 0.55)));
        jobs.push([9, () => lab.put([[mid[0], mid[1] - 11, 'middle'], [mid[0], mid[1] + 20, 'middle'],
          [p[0] + 9, p[1] + 4, 'start'], [p[0] - 9, p[1] + 4, 'end'], [p[0], p[1] + 17, 'middle'], [p[0], p[1] - 10, 'middle']],
          `&#956;<tspan dy="2.6" font-size="0.75em">${k}</tspan>${V ? `<tspan dy="-2.6"> ${comp.toFixed(2)} D</tspan>` : ''}`, 'lab mu', fill(ink[k]))]);
      });
      if (sh.mu !== false) {
        over.push(ln(O, muTip, 'dp-mu', ink.mu, ` marker-end="${mk(ink.mu)}"`));
        if (sty.dipoleArrow === 'chemistry') {                        // the crossed tail that marks the positive end
          const po2 = pr(O), pt2 = pr(muTip), ax2 = [pt2[0] - po2[0], pt2[1] - po2[1]];
          const an = Math.hypot(ax2[0], ax2[1]) || 1, cx2 = po2[0] + ax2[0] / an * 7, cy2 = po2[1] + ax2[1] / an * 7;
          over.push(`<line class="dp-mu" x1="${f1(cx2 - -ax2[1] / an * 5)}" y1="${f1(cy2 - ax2[0] / an * 5)}" x2="${f1(cx2 + -ax2[1] / an * 5)}" y2="${f1(cy2 + ax2[0] / an * 5)}" stroke="${ink.mu}"/>`);
        }
        lab.along(po, pr(muTip));
        const p = pr(muTip);
        jobs.push([10, () => lab.put([[p[0] + 9, p[1] - 4, 'start'], [p[0] - 9, p[1] - 4, 'end'], [p[0], p[1] - 12, 'middle'], [p[0], p[1] + 18, 'middle']],
          V ? `&#956; ${dpx.total.toFixed(2)} D` : '&#956;', 'lab mu', fill(ink.mu), p)]);
      }
      over.push(`<circle class="dp-origin" cx="${f1(po[0])}" cy="${f1(po[1])}" r="4.2" fill="${ink.origin}"/>`);
      if (DL.rawAxes) {                      // the file's own axes: what the feather's dip_x/y/z are measured on
        const O0 = DL.atOrigin ? [0, 0, 0] : O, po0 = pr(O0);   // with atOrigin they sit where the file puts them
        [[[1, 0, 0], 'x'], [[0, 1, 0], 'y'], [[0, 0, 1], 'z']].forEach(([e, k], n) => {
          const tip = add(O0, scl(e, 1.55)), pt = pr(tip);
          over.push(ln(O0, tip, 'dp-raw', ink.text, ` marker-end="${mk(ink.text)}"`));
          lab.along(po0, pt);
          jobs.push([8, () => lab.put([[pt[0] + 10, pt[1] + 4, 'start'], [pt[0] - 10, pt[1] + 4, 'end'],
            [pt[0], pt[1] - 9, 'middle'], [pt[0], pt[1] + 16, 'middle']],
            `${k}&#770;`, 'dpn', fill(ink.text), pt)]);
        });
        legend.push({ col: ink.text, text: `on the file&#8217;s axes &#956;<tspan dy="2.6" font-size="0.75em">x</tspan><tspan dy="-2.6"> ${dpx.mu[0].toFixed(2)}</tspan>`
          + ` &#183; &#956;<tspan dy="2.6" font-size="0.75em">y</tspan><tspan dy="-2.6"> ${dpx.mu[1].toFixed(2)}</tspan>`
          + ` &#183; &#956;<tspan dy="2.6" font-size="0.75em">z</tspan><tspan dy="-2.6"> ${dpx.mu[2].toFixed(2)}</tspan>` });
      }
      if (DL.atOrigin) {                     // the same vector where the file puts it: no frame, no origin of ours
        const Z = [0, 0, 0], pz = pr(Z), zt = add(Z, scl(dpx.mu, muScale)), pzt = pr(zt);
        over.push(ln(Z, zt, 'dp-mu', ink.mu, ` marker-end="${mk(ink.mu)}"`),                 // the vector itself, in mu ink
          `<circle class="dp-origin raw" cx="${f1(pz[0])}" cy="${f1(pz[1])}" r="3.4" fill="${ink.text}"/>`);
        lab.along(pz, pzt);
        jobs.push([10, () => lab.put([[pzt[0] + 9, pzt[1] - 4, 'start'], [pzt[0] - 9, pzt[1] - 4, 'end'],
          [pzt[0], pzt[1] - 12, 'middle'], [pzt[0], pzt[1] + 18, 'middle']],
          `&#956; at the file origin${V ? `, ${dpx.total.toFixed(2)} D` : ''}`, 'lab', fill(ink.mu), pzt)]);
      }
      const perD = Math.abs(muScale);
      if (!muSaid) legend.push({ muSaid: (muSaid = true), col: ink.mu, text: `&#956; arrow ${perD.toFixed(2)} &#8491;/D, ${sty.dipoleArrow === 'chemistry' ? '+ &#8594; &#8722; (chemistry)' : '&#8722; &#8594; + (Gaussian vector)'}` });
      legend.push({ col: ink.mu, text: `&#956;<tspan dy="2.6" font-size="0.75em">u</tspan><tspan dy="-2.6"> ${dpx.dipole_y.toFixed(2)}</tspan> · &#956;<tspan dy="2.6" font-size="0.75em">v</tspan><tspan dy="-2.6"> ${dpx.dipole_x.toFixed(2)}</tspan> · &#956;<tspan dy="2.6" font-size="0.75em">w</tspan><tspan dy="-2.6"> ${dpx.dipole_z.toFixed(2)}</tspan> · |&#956;| ${dpx.total.toFixed(2)} D` });
    });

    (L.sites || []).forEach((q0) => {        // charge sites; the rings are obstacles for every label
      const q = mol.charges && mol.charges[q0.type] && mol.charges[q0.type][q0.i];
      if (q == null || !P[q0.i]) return;
      const [x, y] = P[q0.i], r = rb(q0.i), colr = ink.mono ? ink.text : q0.color, d = r + 10;
      over.push(`<circle class="q-ring" cx="${f1(x)}" cy="${f1(y)}" r="${f1(r + 4)}" stroke="${colr}"/>`);
      lab.circles.push([x, y, r + 8]);
      jobs.push([11, () => lab.put([[x + d, y - d * 0.4, 'start'], [x - d, y - d * 0.4, 'end'], [x, y - d - 4, 'middle'], [x, y + d + 12, 'middle'],
        [x + d, y + d, 'start'], [x - d, y + d, 'end']],
      `q<tspan dy="2.6" font-size="0.72em">${QNAME[q0.type] || q0.type}</tspan>${V ? `<tspan dy="-2.6"> ${q >= 0 ? '+' : '&#8722;'}${Math.abs(q).toFixed(3)}</tspan>` : ''}`, 'lab q', fill(colr), [x, y])]);
    });

    // ---- pass 2: labels, most important first, now that every line is known ----
    // featureLabels off: the constructions are drawn, nothing is written over the molecule
    const dropLab = (sty.dropLabels || []).map((s) => new RegExp(s));   // style.dropLabels: patterns of labels not to write
    const labels = sty.featureLabels === false ? []
      : jobs.map((j, k) => [j[0], k, j[1]]).sort((a, b) => a[0] - b[0] || a[1] - b[1]).map(([, , f]) => f()).filter((t) => !dropLab.some((re) => re.test(String(t))));

    const alabs = [];
    if (sty.atomLabels && sty.atomLabels !== 'off') {
      el.forEach((e, i) => {
        if (fade.has(i) || hidden.has(i) || (sty.atomLabels === 'heavy' && e === 'H')) return;
        alabs.push(`<text class="alab" x="${f1(P[i][0])}" y="${f1(P[i][1] + 3 * K)}" text-anchor="middle" font-size="${f1(FS.alab * K)}" fill="#1F2937">${sty.atomLabels === 'num' ? i + 1 : e + (i + 1)}</text>`);
      });
    }

    // keys sit just under the drawing, not pinned to the bottom edge
    const yDrawn = Math.max(...structs.flatMap((t) => t.P.map((p, i) => (t.hidden.has(i) ? -Infinity : p[1] + rbOf(t, i)))), ...ends.map((q) => pr(q)[1]));
    const keys = [];
    if (qmap) {
      const x0 = 22, y0 = Math.min(H - 30, yDrawn + 26), w = 120;
      const ramp = RAMPS[sty.qRamp] || RAMPS.bwr;
      keys.push(`<defs><linearGradient id="qbar"><stop offset="0" stop-color="${ramp.neg}"/><stop offset=".5" stop-color="#FFFFFF"/><stop offset="1" stop-color="${ramp.pos}"/></linearGradient></defs>`,
        `<rect x="${x0}" y="${y0}" width="${w}" height="8" rx="2" fill="url(#qbar)" stroke="#B9C0CA" stroke-width=".6"/>`,
        `<text class="key" x="${x0}" y="${y0 - 5}" font-size="${f1(FS.key * K)}" fill="${ink.key}">q<tspan dy="2.6" font-size="0.72em">${QNAME[L.map.type] || L.map.type}</tspan></text>`,
        `<text class="key" x="${x0}" y="${y0 + 20}" font-size="${f1(FS.key * K)}" fill="${ink.key}">&#8722;${qmax.toFixed(2)} e</text>`,
        `<text class="key" x="${x0 + w}" y="${y0 + 20}" text-anchor="end" font-size="${f1(FS.key * K)}" fill="${ink.key}">+${qmax.toFixed(2)} e${sty.qMax > 0 ? ' (pinned)' : ''}</text>`);
    }
    if (structs.length > 1) {                // which structure is which
      if (L.primaryName) legend.push({ col: '#6B7482', text: `${esc(L.primaryName)} (measured)` });
      structs.slice(1).forEach((t) => legend.push({ col: sty.overlay === 'tint' ? t.color : '#B9C0CA', text: esc(t.name) + (t.rmsd != null ? ` · RMSD ${t.rmsd.toFixed(3)} &#8491;` : '') }));
    }
    if (L.legend !== false && legend.length) {
      const dy = 15 * K;
      legend.forEach((e, k) => {
        const y = Math.min(H - 12, yDrawn + 26 + (legend.length - 1) * dy) - (legend.length - 1 - k) * dy;
        keys.push(`<circle cx="${W - 14}" cy="${f1(y - 3.5 * K)}" r="${f1(3.6 * K)}" fill="${e.col}"/>`,
          `<text class="key" x="${W - 22}" y="${f1(y)}" text-anchor="end" font-size="${f1(FS.key * K)}" fill="${ink.key}">${e.text}</text>`);
      });
    }

    const defs = spheres(shades, pre.shade) +
      [...markers].map((col) => `<marker id="mk${col.slice(1)}" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="5.2" markerHeight="5.2" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="${col}"/></marker>`).join('') +
      (sty.shadow ? '<filter id="softshadow" x="-20%" y="-20%" width="140%" height="140%"><feDropShadow dx="1.4" dy="2.4" stdDeviation="2.3" flood-color="#1B2230" flood-opacity=".22"/></filter>' : '');
    return [`<svg class="fig scene${sty.thetaBold ? ' theta-bold' : ''}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(opt.title || 'Molecule with descriptor constructions')}">`,
      `<defs>${defs}</defs>`, ...under,
      sty.vdw === 'back' ? vdwLayer : '',
      sty.shadow ? '<g filter="url(#softshadow)">' : '<g>', ...items.sort((p, q) => p[0] - q[0]).map(([, g]) => g), '</g>',
      sty.vdw === 'front' ? vdwLayer : '',
      ...alabs, ...over, ...labels, ...keys,
      ...sideBoxes.flatMap((b) => (b.kind === 'inset' ? inset(b.ax, W - IS - 10, b.y, IS, ink, K, V)
        : profilePlot(b.ax, W - PW - 10, b.y, PW, PH, ink, K, V))),
      ...hits.sort((p, q) => p[0] - q[0]).map(([, h]) => h), '</svg>'].filter(Boolean).join('\n');
  }

  const SCENE_CSS = `
.scene text{font-family:"IBM Plex Sans",Helvetica,Arial,sans-serif;paint-order:stroke;stroke:#FFFFFF;stroke-width:3.2px;stroke-linejoin:round}
.scene text.th,.scene text.dpn{font-family:"STIX Two Text",Georgia,"Times New Roman",serif;font-style:italic}
.scene text.dpn{font-weight:600;stroke-width:3px}
.scene text.vl{font-style:italic}
.scene text.mark,.scene text.alab{font-family:"IBM Plex Mono",Menlo,Consolas,monospace}
.scene text.alab{stroke-width:2.4px}
.scene text.mu,.scene text.q{font-weight:600}
.scene .ghost{opacity:.18}
.scene .behind{opacity:.3}
.scene .sx-axis{fill:none;stroke-width:2.2;stroke-dasharray:7 4;stroke-linecap:round}
.scene .sx-b1{stroke-width:2;stroke-linecap:round}
.scene .sx-trace{stroke-width:3;stroke-linecap:round;opacity:.55}
.scene .sx-b5{stroke-width:1.9;stroke-dasharray:5 3.5;stroke-linecap:round}
.scene .sx-b5.dim{stroke-dasharray:none;stroke-width:1.5}
.scene .sx-loctick{stroke-width:3.4;stroke-linecap:round}
.scene .sx-tickhalo{stroke-width:7;stroke-linecap:round;opacity:.9}
.scene .sx-caret{stroke:#FFFFFF;stroke-width:1.1}
.scene .sx-cross{fill:none;stroke-width:3.2;stroke-linecap:round}
.scene .sx-run{stroke-width:3.8;stroke-linecap:round}
.scene .sx-runhalo{stroke-width:7.5;stroke-linecap:round;opacity:.85}
.scene .sx-collar{fill:none;stroke-width:2.6}
.scene .sx-collarhalo{fill:none;stroke-width:5.5;opacity:.85}
.scene .sx-locline{stroke-width:1;stroke-dasharray:4 4;opacity:.75}
.scene .sx-ref{fill-opacity:.06;stroke-opacity:.35;stroke-width:1}
.scene .sx-patch{fill-opacity:.08;stroke-opacity:.7;stroke-width:1}
.scene .trial{opacity:.55;stroke-dasharray:5 3}
.scene .sx-patch.trial{fill-opacity:.12;stroke-dasharray:4 3}
.scene .sx-wedge{fill-opacity:.2;stroke:none}
.scene .sx-arc{fill:none;stroke-width:2.2;stroke-dasharray:2 4.5;stroke-linecap:round}
.scene .sx-phi{stroke-width:1.3;stroke-linecap:round}
.scene .sx-grad{stroke-width:1.5;stroke-linecap:round}
.scene .sx-dimend{stroke-width:1.8;stroke-linecap:round}
.scene .sx-tape{stroke-width:1;fill-opacity:0}
.scene .sx-tape.on{fill-opacity:.55}
.scene .sx-pin{fill-opacity:.85}
.scene .cone-face{fill-opacity:.13}
.scene .cone-rim{stroke-width:1.4;stroke-dasharray:4 3;opacity:.85}
.scene .cone-axis{stroke-width:1.6;stroke-dasharray:6 4}
.scene .cone-edge{stroke-width:1.8}
.scene .bv-ball{fill-opacity:.07;stroke-width:1.4;stroke-dasharray:5 4}
.scene .bv-dots{fill-opacity:.55}
.scene.theta-bold .sx-wedge{fill-opacity:.55}
.scene.theta-bold .sx-arc{stroke-width:4.2;stroke-dasharray:none}
.scene.theta-bold .sx-shadow{stroke-width:2.8}
.scene.theta-bold .sx-v,.scene.theta-bold .sx-v.dash{stroke-width:3}
.scene .sx-shadow{stroke-width:1.8}
.scene .sx-perp{stroke-width:1.8;stroke-dasharray:1.5 4;stroke-linecap:round}
.scene .sx-right{fill:none;stroke-width:1.4}
.scene .sx-v{stroke-width:1.8}
.scene .sx-v.dash{stroke-width:2;stroke-dasharray:5 3.5;stroke-linecap:round}
.scene .sx-ring{fill:none;stroke-width:1.3}
.scene .sx-touch{fill:none;stroke-width:1.1}
.scene .dp-ax{stroke-width:2.4;stroke-linecap:round}
.scene .dp-mu{stroke-width:3.6;stroke-linecap:round}
.scene .dp-raw{stroke-width:2.4;stroke-dasharray:6 4;opacity:.8;stroke-linecap:round}
.scene .dp-origin.raw{opacity:.8}
.scene .dp-comp{stroke-width:3.4;stroke-linecap:round}
.scene .dp-guide{stroke-width:1.4;stroke-dasharray:4 3;opacity:.9}
.scene .dp-drop{stroke-width:1.1;stroke-dasharray:2 3;opacity:.8}
.scene .dp-spoke{stroke-width:.9;stroke-dasharray:2 3;opacity:.8}
.scene .dp-origin{stroke:#FFFFFF;stroke-width:1.6}
.scene .q-ring{fill:none;stroke-width:2.6}
.scene .ins-box{fill:#FFFFFF;fill-opacity:.92;stroke:#D3D9E1;stroke-width:1}
.scene .vdw circle{stroke:none}
.scene .vdw.density path{stroke-width:.6;stroke-linejoin:round}
.scene .pf-line{fill:none;stroke-width:1.6;stroke-linejoin:round}
.scene .pf-band{fill-opacity:.13;stroke:none}
.scene .pf-drop{stroke-width:.9;stroke-dasharray:3 2;opacity:.7}
.scene .pf-now{stroke-width:1.4;opacity:.85}
.scene .leader{stroke-width:.9;opacity:.6}
.scene text[data-lab]{cursor:grab}
.scene text[data-lab]:active{cursor:grabbing}
.scene .pickring{fill:none;stroke-width:2.2}
.scene .pickring.pa{stroke:#151A21} .scene .pickring.pb{stroke:#3450A8}
.scene .pickring.hi{stroke-width:2.6;stroke-dasharray:4 3}`;

  // a figure file that stands on its own: no hit targets, its own styles, a white page unless transparent.
  // mm gives it a real physical size (so it lands in a manuscript at the intended column width);
  // state embeds the session that drew it, which the page can read back.
  function standalone(svg, { bg = 'white', scale = 2, mm = 0, state = null } = {}) {
    svg = svg.replace(/<circle class="hit"[\s\S]*?<\/circle>/g, '').replace(/<circle class="pickring[^>]*?(?:\/>|><\/circle>)/g, '');   // written either way, as a string or re-serialized by the page
    const m = svg.match(/viewBox="0 0 (\d+) (\d+)"/), w = m ? +m[1] : W, h = m ? +m[2] : H;
    const size = mm > 0 ? `width="${f1(mm)}mm" height="${f1(mm * h / w)}mm"` : `width="${w * scale}" height="${h * scale}"`;
    const meta = state ? `<metadata id="thetascene">${String(JSON.stringify(state)).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')}</metadata>` : '';
    return svg.replace(/<svg([^>]*)>/, `<svg$1 ${size}>${meta}<style>${SCENE_CSS}</style>` +
      (bg === 'transparent' ? '' : '<rect width="100%" height="100%" fill="#FFFFFF"/>'));
  }
  // what the drawing measures in points once it is placed at a given width, for the type-size check
  function typeCheck(svg, mm) {
    const pt = (u) => u * (mm / W) * 2.8346457;
    const sizes = [...svg.matchAll(/font-size="([\d.]+)"/g)].map((x) => +x[1]);
    const strokes = [...svg.matchAll(/stroke-width[:="]+([\d.]+)/g)].map((x) => +x[1]).filter((x) => x > 0.05);
    return { minTextPt: sizes.length ? pt(Math.min(...sizes)) : null, minLinePt: strokes.length ? pt(Math.min(...strokes)) : null,
      widthMm: mm, heightMm: mm * H / W };
  }

  const api = { scene, axisGeometry, axisView, displayBonds, superpose, standalone, typeCheck, Placer, shade, f1, esc, ballColor, covalent,
    PAL, INKS, PRESETS, PALETTES, LOOKS, LOOK_KEYS, lookStyle, STYLE, SCENE_CSS, W, H };
  root.ThetaScene = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);
