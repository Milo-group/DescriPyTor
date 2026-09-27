// Every paper preset of the built page, rendered to SVG with the page's own code.
//   python build.py && node render_presets.js <out_dir>     (then PNGs via headless Chrome)
// The molecules and presets are read out of theta_explorer.html itself, so the files
// are exactly what the page shows when a preset is applied (in the default look).
const fs = require('fs');
const path = require('path');
const S = require('./sterimol.js');
const SC = require('./scene.js');

const out = process.argv[2] || path.join(__dirname, '_presets');
fs.mkdirSync(out, { recursive: true });
const html = fs.readFileSync(path.join(__dirname, 'theta_explorer.html'), 'utf8');
const i0 = html.indexOf('const MOLS = ') + 'const MOLS = '.length;
const MOLS = JSON.parse(html.slice(i0, html.indexOf('const REFS = ', i0)).trim().replace(/;$/, ''));

const selList = (t) => String(t || '').split(/[\s,;]+/).filter(Boolean).flatMap((p) => {
  const m = p.match(/^(\d+)(?:-(\d+))?$/); if (!m) return [];
  const a = +m[1], b = m[2] ? +m[2] : a; return Array.from({ length: Math.abs(b - a) + 1 }, (_, k) => Math.min(a, b) + k);
});
const selArg = (t) => { const l = selList(t); return l.length > 1 ? l : l.length ? l[0] : null; };
function pca(X) {                          // same as the page's principalAxes
  const n = X.length, c = [0, 1, 2].map((k) => X.reduce((s, p) => s + p[k], 0) / n);
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const C = [0, 1, 2].map((i) => [0, 1, 2].map((j) => X.reduce((s, p) => s + (p[i] - c[i]) * (p[j] - c[j]), 0)));
  const pw = (M, v) => { for (let k = 0; k < 80; k++) { const w = M.map((r) => dot(r, v)); const q = Math.hypot(...w) || 1; v = w.map((x) => x / q); } return v; };
  const v1 = pw(C, [1, 0.3, 0.2]), l1 = dot(v1, C.map((r) => dot(r, v1)));
  let v2 = pw(C.map((r, i) => r.map((x, j) => x - l1 * v1[i] * v1[j])), [0.2, 1, 0.3]);
  v2 = v2.map((x, i) => x - dot(v2, v1) * v1[i]); const q = Math.hypot(...v2) || 1; v2 = v2.map((x) => x / q);
  return [v1, v2, [v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]]];
}

const written = [];
for (const m of MOLS) {
  const mol = S.parseXYZ(m.xyz);
  mol.bonds = S.connectivity(mol.el, mol.X); mol.types = S.atomTypes(mol.el, mol.bonds);
  mol.displayBonds = SC.displayBonds(mol.el, mol.X, mol.bonds).all;
  mol.dipole = m.dipole || null; mol.charges = m.charges || null;
  const lbl = (k) => mol.el[k - 1] + k;
  (m.presets || []).forEach((pr, k) => {
    const L = pr.layers;
    const axes = L.axes.map((ax) => ({ res: S.sterimol(mol, ax.a, ax.b), color: ax.color, show: ax.show, label: `${lbl(ax.a)}→${lbl(ax.b)}`, a: ax.a }));
    const fade = new Set(selList(L.fade).map((i) => i - 1));
    if (L.fadeOff && axes.length) {        // as the page's layerInputs
      const keep = new Set(); axes.forEach((x) => { x.res.atoms.forEach((i) => keep.add(i)); keep.add(x.a - 1); });
      mol.el.forEach((_, i) => { if (!keep.has(i)) fade.add(i); });
    }
    const inputs = {
      axes,
      dipole: L.dip.on && mol.dipole ? { res: S.dipoleFrame(mol.X, mol.dipole, selArg(L.dip.origin), selArg(L.dip.y), selArg(L.dip.plane)),
        originSet: Array.isArray(selArg(L.dip.origin)) ? selArg(L.dip.origin) : [], show: L.dip.show } : null,
      sites: L.sites.map((q) => ({ i: q.i - 1, type: q.type, color: q.color })),
      map: L.map ? { type: L.map } : null, fade, legend: L.legend,
    };
    const R = pr.view === 'axis' && axes.length ? SC.axisView(axes[Math.max(0, L.active)]) : pca(mol.X);
    const svg = SC.standalone(SC.scene(mol, inputs, { R, zoom: 1 }));
    const file = `${m.name}_${k + 1}_${pr.name.split(':')[0].replace(/[^\w]+/g, '_')}`.replace(/_+$/, '');
    fs.writeFileSync(path.join(out, file + '.svg'), svg);
    fs.writeFileSync(path.join(out, file + '.html'), `<!doctype html><html><body style="margin:0;background:#fff">${svg}</body></html>`);
    written.push(file);
  });
}
console.log(written.join('\n'));
