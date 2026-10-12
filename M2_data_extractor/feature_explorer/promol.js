/* The promolecular electron density, and its isosurface.
 *
 * rho(r) = sum over atoms of c_e * exp(-alpha_e * r), in bohr: a sum of spherical atomic
 * densities with no wavefunction behind it -- the approximation NCIPLOT uses.
 *
 * The three fitted elements come from Gaussian SCF cubes (Getting_started_with_examples/
 * cube_example): c and alpha were fitted on Ad_1_a in the 3e-4..3e-2 au shell that sets the
 * surface, then tested on Bn_1_a and Cy_1_a, which they had never seen. Against those two the
 * rho = 0.003 isosurface lands within 0.04 A of the SCF one (Jaccard 0.97) and its Sterimol
 * L/B1/B5 within 0.1 A. Every other element falls back to a tail of alpha = 2.0 bohr^-1 scaled
 * so that a lone atom's 0.003 surface sits at its CPK radius -- for the fitted three that rule
 * would have given 1.59/1.08/1.50 A against their CPK 1.70/1.10/1.52, so it is the same picture.
 *
 * The isosurface is meshed with naive surface nets: one vertex per cell that the surface cuts,
 * quads across the grid edges it crosses. Blobby by construction, which suits a density.
 */
(function (root) {
  const BOHR = 0.529177249;
  const ISO0 = 0.003;                        // the isovalue cube_sterimol reads the surface at
  const FIT = { H: [0.376, 2.359], C: [1.169, 1.988], O: [4.120, 2.546] };
  const VDW = { H: 1.1, C: 1.7, N: 1.55, O: 1.52, F: 1.47, S: 1.8, Cl: 1.75, Br: 1.85, I: 1.98, P: 1.8, B: 1.92, Si: 2.1,
    Cu: 1.96, Ni: 1.97, Pd: 2.1, Fe: 2.04, Co: 2.0, Mn: 2.05, Cr: 2.06, Zn: 2.01, Ru: 2.13, Rh: 2.1, Ag: 2.11, Ir: 2.13,
    Pt: 2.09, Au: 2.14, Mo: 2.17, W: 2.18, Ti: 2.11, V: 2.07, Mg: 1.73, Ca: 2.31, Li: 1.82, Na: 2.27, K: 2.75,
    Sn: 2.17, Pb: 2.02, Hg: 2.23, Cd: 2.18 };

  // [c, alpha] for an element: fitted where we have a cube, else a tail that puts the lone-atom
  // 0.003 surface on the CPK radius
  const atomDens = (e) => FIT[e] || [ISO0 * Math.exp(2.0 * (VDW[e] || 1.8) / BOHR), 2.0];
  const atomRadius = (e, iso) => { const [c, a] = atomDens(e); return Math.log(c / iso) / a * BOHR; };

  function field(el, X, iso) {               // rho - iso at a point, and its gradient, in Angstrom
    const P = el.map((e) => atomDens(e));
    return {
      at(x, y, z) {
        let s = 0;
        for (let i = 0; i < el.length; i++) {
          const dx = x - X[i][0], dy = y - X[i][1], dz = z - X[i][2];
          const d = Math.sqrt(dx * dx + dy * dy + dz * dz) / BOHR;
          s += P[i][0] * Math.exp(-P[i][1] * d);
        }
        return s - iso;
      },
      grad(x, y, z) {
        let gx = 0, gy = 0, gz = 0;
        for (let i = 0; i < el.length; i++) {
          const dx = x - X[i][0], dy = y - X[i][1], dz = z - X[i][2];
          const d = Math.sqrt(dx * dx + dy * dy + dz * dz) / BOHR || 1e-9;
          const k = -P[i][1] * P[i][0] * Math.exp(-P[i][1] * d) / (d * BOHR * BOHR);
          gx += k * dx; gy += k * dy; gz += k * dz;
        }
        return [gx, gy, gz];
      },
    };
  }

  const CORNERS = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
  const EDGES = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]];

  // -> {verts: [[x,y,z]], norms: [[x,y,z]], near: [atom index], quads: [[a,b,c,d]]}
  function mesh(el, X, { iso = ISO0, step = 0.4, pad = 2.6, skip = null } = {}) {
    const keep = el.map((_, i) => !(skip && skip.has(i)));
    const use = el.filter((_, i) => keep[i]), Xu = X.filter((_, i) => keep[i]);
    const idx = el.map((_, i) => i).filter((i) => keep[i]);
    if (!use.length) return { verts: [], norms: [], near: [], quads: [] };
    const lo = [0, 1, 2].map((k) => Math.min(...Xu.map((p) => p[k])) - pad);
    const hi = [0, 1, 2].map((k) => Math.max(...Xu.map((p) => p[k])) + pad);
    const n = [0, 1, 2].map((k) => Math.max(2, Math.ceil((hi[k] - lo[k]) / step) + 1));
    const f = field(use, Xu, iso);
    const at = (i, j, k) => [lo[0] + i * step, lo[1] + j * step, lo[2] + k * step];
    const V = new Float32Array(n[0] * n[1] * n[2]);
    const vi = (i, j, k) => (i * n[1] + j) * n[2] + k;
    for (let i = 0; i < n[0]; i++) for (let j = 0; j < n[1]; j++) for (let k = 0; k < n[2]; k++) {
      const p = at(i, j, k); V[vi(i, j, k)] = f.at(p[0], p[1], p[2]);
    }
    const cell = new Int32Array((n[0] - 1) * (n[1] - 1) * (n[2] - 1)).fill(-1);
    const ci = (i, j, k) => (i * (n[1] - 1) + j) * (n[2] - 1) + k;
    const verts = [], norms = [], near = [];
    for (let i = 0; i < n[0] - 1; i++) for (let j = 0; j < n[1] - 1; j++) for (let k = 0; k < n[2] - 1; k++) {
      const s = CORNERS.map(([a, b, c]) => V[vi(i + a, j + b, k + c)]);
      let neg = 0;
      for (const v of s) if (v < 0) neg++;
      if (neg === 0 || neg === 8) continue;
      let px = 0, py = 0, pz = 0, m = 0;     // the cell's vertex: mean of the crossings on its edges
      for (const [a, b] of EDGES) {
        if ((s[a] < 0) === (s[b] < 0)) continue;
        const t = s[a] / (s[a] - s[b]);
        px += CORNERS[a][0] + t * (CORNERS[b][0] - CORNERS[a][0]);
        py += CORNERS[a][1] + t * (CORNERS[b][1] - CORNERS[a][1]);
        pz += CORNERS[a][2] + t * (CORNERS[b][2] - CORNERS[a][2]);
        m++;
      }
      const p = [lo[0] + (i + px / m) * step, lo[1] + (j + py / m) * step, lo[2] + (k + pz / m) * step];
      const g = f.grad(p[0], p[1], p[2]), gn = Math.hypot(g[0], g[1], g[2]) || 1;
      let best = 0, bd = Infinity;
      for (let q = 0; q < Xu.length; q++) {
        const d = (p[0] - Xu[q][0]) ** 2 + (p[1] - Xu[q][1]) ** 2 + (p[2] - Xu[q][2]) ** 2;
        if (d < bd) { bd = d; best = q; }
      }
      cell[ci(i, j, k)] = verts.length;
      verts.push(p); norms.push([-g[0] / gn, -g[1] / gn, -g[2] / gn]); near.push(idx[best]);
    }
    // a quad for every grid edge the surface crosses: the four cells that share it
    const quads = [];
    const face = (a, b, c, d, flip) => { if (a >= 0 && b >= 0 && c >= 0 && d >= 0) quads.push(flip ? [d, c, b, a] : [a, b, c, d]); };
    for (let i = 1; i < n[0] - 1; i++) for (let j = 1; j < n[1] - 1; j++) for (let k = 1; k < n[2] - 1; k++) {
      const v0 = V[vi(i, j, k)] < 0;
      if (V[vi(i + 1, j, k)] < 0 !== v0) face(cell[ci(i, j - 1, k - 1)], cell[ci(i, j, k - 1)], cell[ci(i, j, k)], cell[ci(i, j - 1, k)], v0);
      if (V[vi(i, j + 1, k)] < 0 !== v0) face(cell[ci(i - 1, j, k - 1)], cell[ci(i, j, k - 1)], cell[ci(i, j, k)], cell[ci(i - 1, j, k)], !v0);
      if (V[vi(i, j, k + 1)] < 0 !== v0) face(cell[ci(i - 1, j - 1, k)], cell[ci(i, j - 1, k)], cell[ci(i, j, k)], cell[ci(i - 1, j, k)], v0);
    }
    return { verts, norms, near, quads, iso, step };
  }

  // one mesh per (molecule, isovalue, step, hidden atoms); rebuilt only when one of those changes
  const CACHE = new Map();
  function cached(el, X, opt) {
    const key = [el.join(''), X.length, X[0] && X[0].join(','), opt.iso, opt.step, opt.skip ? [...opt.skip].join(',') : ''].join('|');
    if (!CACHE.has(key)) {
      if (CACHE.size > 8) CACHE.delete(CACHE.keys().next().value);
      CACHE.set(key, mesh(el, X, opt));
    }
    return CACHE.get(key);
  }

  const api = { mesh, cached, field, atomDens, atomRadius, ISO0, FIT };
  root.Promol = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);
