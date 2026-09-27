// Every fixture drawn in every look, with every layer on, checked for the things that break silently.
//   node smoke.js
// Values are gated elsewhere (gate.js); this checks the drawing: no NaN or undefined in the SVG,
// balanced tags, the layers actually present, and the loc B5 marks and label modes all rendering.
const fs = require('fs');
const path = require('path');
const S = require('./sterimol.js');
global.Promol = require('./promol.js');
const SC = require('./scene.js');

const fx = (f) => path.join(__dirname, 'fixtures', f);
const MOLS = JSON.parse(fs.readFileSync(path.join(__dirname, 'theta_explorer.html'), 'utf8')
  .split('const MOLS = ')[1].split('const REFS = ')[0].trim().replace(/;$/, ''));

function mol(entry) {
  const m = S.parseXYZ(entry.xyz);
  m.bonds = S.connectivity(m.el, m.X, entry.thr); m.types = S.atomTypes(m.el, m.bonds);
  m.displayBonds = SC.displayBonds(m.el, m.X, m.bonds).all;
  m.dipole = entry.dipole || null; m.charges = entry.charges || null;
  return m;
}
const selList = (t) => String(t || '').split(/[\s,;]+/).filter(Boolean).flatMap((p) => {
  const m = p.match(/^(\d+)(?:-(\d+))?$/); if (!m) return [];
  const a = +m[1], b = m[2] ? +m[2] : a; return Array.from({ length: Math.abs(b - a) + 1 }, (_, k) => Math.min(a, b) + k);
});
const selArg = (t) => { const l = selList(t); return l.length > 1 ? l : l.length ? l[0] : null; };

let checks = 0, bad = 0;
const golden = {};                                     // preset drawings, kept as regression goldens
const fail = (where, why) => { bad++; console.log(`  ${where}: ${why}`); };
function check(where, svg, want) {
  checks++;
  const m = svg.match(/NaN|undefined|Infinity/);
  if (m) return fail(where, `${m[0]} in the drawing near "${svg.slice(Math.max(0, m.index - 40), m.index + 40).replace(/\s+/g, ' ')}"`);
  const open = (svg.match(/<(circle|line|path|text|g|ellipse|rect|marker|radialGradient)\b/g) || []).length;
  const close = (svg.match(/<\/(circle|line|path|text|g|ellipse|rect|marker|radialGradient)>|\/>/g) || []).length;
  if (close < open) return fail(where, `${open} elements opened, ${close} closed`);
  for (const [cls, n] of Object.entries(want || {})) {
    const got = (svg.match(new RegExp(`class="[^"]*\\b${cls}\\b`, 'g')) || []).length;
    if (got < n) return fail(where, `expected at least ${n} .${cls}, found ${got}`);
  }
}

const ALL = { axis: true, L: true, b1: true, plane: true, loc: true, b5: true, theta: true, phi: true, drop: true, inset: true };
const pca = (X) => {                                   // a stable view, not the point of the test
  const n = X.length, c = [0, 1, 2].map((k) => X.reduce((s, p) => s + p[k], 0) / n);
  const d = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const C = [0, 1, 2].map((i) => [0, 1, 2].map((j) => X.reduce((s, p) => s + (p[i] - c[i]) * (p[j] - c[j]), 0)));
  const pw = (M, v) => { for (let k = 0; k < 60; k++) { const w = M.map((r) => d(r, v)); const q = Math.hypot(...w) || 1; v = w.map((x) => x / q); } return v; };
  const v1 = pw(C, [1, 0.3, 0.2]), l1 = d(v1, C.map((r) => d(r, v1)));
  let v2 = pw(C.map((r, i) => r.map((x, j) => x - l1 * v1[i] * v1[j])), [0.2, 1, 0.3]);
  v2 = v2.map((x, i) => x - d(v2, v1) * v1[i]); const q = Math.hypot(...v2) || 1; v2 = v2.map((x) => x / q);
  return [v1, v2, [v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]]];
};

for (const entry of MOLS) {
  const m = mol(entry), R = pca(m.X);
  const a2 = entry.a || 1, b2 = entry.b || 2, a = a2, b = b2;
  const axes = [{ res: S.sterimol(m, a, b), color: '#1596A6', show: ALL, label: 'test', tag: '1' }];
  // every look, with every Sterimol layer on
  for (const look of Object.keys(SC.LOOKS)) {
    const style = Object.assign({}, SC.STYLE, SC.lookStyle(look));
    check(`${entry.name} / ${look}`, SC.scene(m, { axes, legend: true, style }, { R, zoom: 1 }),
      { ball: 1, 'sx-axis': 1, 'sx-b1': 1, 'sx-wedge': 1, 'ins-box': 1 });
  }
  // the phi arc must subtend the phi it is labelled with: chord = 2 r sin(phi/2)
  {
    const svg = SC.scene(m, { axes, style: Object.assign({}, SC.STYLE, { values: true }) }, { R, zoom: 1 });
    const a = svg.match(/class="sx-phi" d="M([-\d.]+),([-\d.]+) A([\d.]+) [\d.]+ 0 0 \d ([-\d.]+),([-\d.]+)"/);
    checks++;
    const phi = S.sterimol(m, a2, b2).phi;
    if (!a) { if (phi > 0.05) fail(entry.name, 'phi is ' + phi.toFixed(1) + ' but no arc was drawn'); }
    else {
      const rr = +a[3], chord = Math.hypot(+a[4] - +a[1], +a[5] - +a[2]);
      const want = 2 * rr * Math.sin(phi * Math.PI / 360);
      if (Math.abs(chord - want) > 0.4) fail(entry.name, `phi arc chord ${chord.toFixed(2)} px, ${phi.toFixed(1)} deg wants ${want.toFixed(2)}`);
    }
  }
  // the five ways to mark the axis
  for (const axisMark of ['arrow', 'ruler', 'tape', 'dim', 'pin']) {
    const style = Object.assign({}, SC.STYLE, { axisMark });
    const want = { arrow: 'sx-axis', ruler: 'sx-grad', tape: 'sx-tape', dim: 'sx-dimend', pin: 'sx-pin' }[axisMark];
    check(`${entry.name} / axis ${axisMark}`, SC.scene(m, { axes, style }, { R, zoom: 1 }), { [want]: 1 });
  }
  // the four loc B5 marks, and labels with and without their values
  for (const locMark of ['tick', 'caret', 'run', 'cross']) {
    for (const values of [true, false]) {
      const style = Object.assign({}, SC.STYLE, { locMark, values });
      check(`${entry.name} / loc ${locMark}${values ? '' : ' (symbols)'}`, SC.scene(m, { axes, style }, { R, zoom: 1 }),
        { [locMark === 'run' ? 'sx-run' : locMark === 'caret' ? 'sx-caret' : locMark === 'cross' ? 'sx-cross' : 'sx-loctick']: 1 });
    }
  }
  // the promolecular surface: a closed mesh, at the radius the density says
  {
    const svg = SC.scene(m, { axes, style: Object.assign({}, SC.STYLE, { vdw: 'front', surfKind: 'density' }) }, { R, zoom: 1 });
    checks++;
    const g = svg.match(/<g class="vdw density"[^>]*>([\s\S]*?)<\/g>/);
    if (!g) fail(entry.name, 'no density surface');
    else if ((g[1].match(/<path/g) || []).length < 200) fail(entry.name, `density surface has only ${(g[1].match(/<path/g) || []).length} faces`);
  }

  // the van der Waals layer: one group, one circle per drawn atom, and nothing at all when off
  for (const where of ['front', 'back']) {
    const svg = SC.scene(m, { axes, style: Object.assign({}, SC.STYLE, { vdw: where }) }, { R, zoom: 1 });
    checks++;
    const g = svg.match(/<g class="vdw" opacity="[\d.]+">([\s\S]*?)<\/g>/);
    if (!g) { fail(entry.name, `no vdW layer for vdw:${where}`); continue; }
    const n = (g[1].match(/<circle/g) || []).length;
    if (n !== m.el.length) fail(entry.name, `vdW layer has ${n} spheres for ${m.el.length} atoms`);
    const at = svg.indexOf('<g class="vdw"'), mol = svg.indexOf('class="ball');
    if (where === 'front' ? at < mol : at > mol) fail(entry.name, `the vdW layer is on the wrong side for vdw:${where}`);
  }

  // feature labels off: the constructions stay, nothing is written over the molecule
  {
    const bare = SC.scene(m, { axes, style: Object.assign({}, SC.STYLE, { featureLabels: false }) }, { R, zoom: 1 });
    checks++;
    if (/data-lab=/.test(bare)) fail(entry.name, 'a label survived featureLabels:false');
    else if (!/class="sx-b1/.test(bare)) fail(entry.name, 'featureLabels:false dropped the construction too');
  }

  // a trial rotation of the scan, replayed, and a held scale instead of a fit
  const prof = axes[0].res.profile;
  const trial = prof[Math.floor(prof.length / 3)];
  check(`${entry.name} / scan trial at ${trial.deg} deg`,
    SC.scene(m, { axes: [Object.assign({}, axes[0], { show: Object.assign({}, ALL, { profile: true }), preview: trial })] }, { R, zoom: 1 }),
    { trial: 3, 'pf-now': 1 });
  const held = { axes, style: SC.STYLE };
  const fitted = SC.scene(m, held, { R, zoom: 1 });
  checks++;
  if (!(held.fitScale > 0)) fail(entry.name, `no fitted scale came back (${held.fitScale})`);
  else if (SC.scene(m, held, { R, zoom: 1, perA: held.fitScale }) !== fitted)
    fail(entry.name, 'holding the scale at the fitted value drew something different');

  // hydrogens, perspective, labels, an overlay of the same structure
  const ov = [{ id: 's2', name: 'overlay', el: m.el, X: m.X.map((p) => [p[0] + 0.4, p[1], p[2]]), bonds: m.displayBonds, types: m.types, color: '#D9822B', rmsd: 0.4 }];
  check(`${entry.name} / stacked, polar H, perspective, labels`,
    SC.scene(m, { axes, overlays: ov, primaryName: entry.name, legend: true, style: Object.assign({}, SC.STYLE, { hydrogens: 'polar', perspective: 1, atomLabels: 'sym', labelScale: 1.4 }) }, { R, zoom: 1 }), { ball: 2 });
  // every preset of this fixture, as the page applies it
  for (const pr of entry.presets || []) {
    const L = pr.layers;
    const A = L.axes.map((x) => ({ res: S.sterimol(m, x.a, x.b), color: x.color, show: x.show, tag: x.tag, label: `${x.a}-${x.b}` }));
    const fade = new Set(selList(L.fade).map((i) => i - 1));
    if (L.fadeOff && A.length) { const keep = new Set(); A.forEach((x, i) => { x.res.atoms.forEach((k) => keep.add(k)); keep.add(L.axes[i].a - 1); }); m.el.forEach((_, i) => { if (!keep.has(i)) fade.add(i); }); }
    const dip = L.dip.on && m.dipole ? { res: S.dipoleFrame(m.X, m.dipole, selArg(L.dip.origin), selArg(L.dip.y), selArg(L.dip.plane)),
      originSet: Array.isArray(selArg(L.dip.origin)) ? selArg(L.dip.origin) : [], show: L.dip.show, colors: L.dip.colors || null } : null;
    const style = Object.assign({}, SC.STYLE, pr.style && pr.style.look ? SC.lookStyle(pr.style.look) : {}, pr.style || {},
      dip && dip.colors ? { inkColors: dip.colors } : {});
    const svg = SC.scene(m, { axes: A, dipole: dip, sites: L.sites.map((q) => ({ i: q.i - 1, type: q.type, color: q.color })),
      map: L.map ? { type: L.map } : null, fade, legend: L.legend, style }, { R: pr.view === 'axis' && A.length ? SC.axisView(A[Math.max(0, L.active)]) : R, zoom: 1 });
    check(`${entry.name} / preset "${pr.name.split(':')[0]}"`, svg);
    golden[`${entry.name}_${pr.name.split(':')[0].trim().replace(/[^\w]+/g, '_')}`] = svg;
  }
}

// the preset drawings must not drift: compare against stored goldens (node smoke.js --update to re-bless)
const GOLD = path.join(__dirname, 'goldens');
const update = process.argv.includes('--update');
if (update) fs.mkdirSync(GOLD, { recursive: true });
for (const [name, svg] of Object.entries(golden)) {
  const f = path.join(GOLD, name + '.svg');
  checks++;
  if (update) { fs.writeFileSync(f, svg); continue; }
  if (!fs.existsSync(f)) { fail(name, 'no golden yet — run node smoke.js --update'); continue; }
  const was = fs.readFileSync(f, 'utf8');
  if (was !== svg) {
    const i = [...was].findIndex((c, k) => c !== svg[k]);
    fail(name, `drawing changed at character ${i}: golden "${was.slice(i, i + 60).replace(/\s+/g, ' ')}" vs now "${svg.slice(i, i + 60).replace(/\s+/g, ' ')}"`);
  }
}

// the substructure layer's two pure pieces (RDKit itself only runs in a browser)
const SM = require('./smarts.js');
[['c1ccccc1', '[#6]1[#6][#6][#6][#6][#6]1'], ['[NX3]C=O', '[NX3]C~O'], ['ClBr', 'ClBr'], ['[#6X3]~[#8X1]', '[#6X3]~[#8X1]'],
  ['C1=NCCO1', 'C1~NCCO1']].forEach(([q, want]) => {
  checks++;
  const got = SM.anyBonds(q);
  if (got !== want) fail('anyBonds', `${q} became ${got}, expected ${want}`);
});
{
  const mb = SM.molblock(['C', 'H', 'Cu'], [[0, 0, 0], [1.1, 0, 0], [0, 2.2, 0]], [[1, 2], [1, 3]]).split(String.fromCharCode(10));
  checks++;
  if (!/^  3  1/.test(mb[3])) fail('molblock', `counts line is "${mb[3]}"`);
  const val = (line) => +line.slice(48, 51);                 // the V2000 valence field
  if (val(mb[4]) !== 1 || val(mb[5]) !== 1 || val(mb[6]) !== 15)
    fail('molblock', `valences ${val(mb[4])}/${val(mb[5])}/${val(mb[6])}, expected 1/1/15 (the metal bond is dropped)`);
}

// superposition still recovers a rigid motion exactly
const m0 = mol(MOLS[0]);
const q = [0.3, -0.5, 0.7, 0.2], n = Math.hypot(...q), [w, x, y, z] = q.map((v) => v / n);
const Rm = [[w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y)],
  [2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x)],
  [2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z]];
const moved = m0.X.map((p) => Rm.map((r, k) => r[0] * p[0] + r[1] * p[1] + r[2] * p[2] + [3, -7, 1.5][k]));
const f = SC.superpose(m0.X, moved);
checks++;
if (f.rmsd > 1e-9) fail('superpose', `rmsd ${f.rmsd} on a rigid motion`);

console.log(`${checks} drawings checked, ${bad} bad`);
process.exit(bad ? 1 : 0);
