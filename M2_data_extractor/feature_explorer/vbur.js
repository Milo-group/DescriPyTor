/* Percent buried volume and the exact cone angle, ported from the package.
 *
 * vbur: metal_complex.buried_volume -- the same Fibonacci 4000 x 12-shell grid,
 * Bondi radii scaled by 1.17, the atom sitting on the centre skipped, shells
 * weighted by r^2. Gated against the Python in gate_pos.js.
 *
 * cone: the exact cone angle (Bilbrey 2013) as the smallest cone from the apex
 * containing every atom's own tangent cone -- found by minimising
 * max_i (angle(u, m_i) + asin(r_i / d_i)) over the axis u, which is what the
 * tangent-set construction solves in closed form. Gated against morfeus.
 */
(function (root) {
  'use strict';
  // Cavallo/Bondi radii, as metal_complex.VDW_BURIAL has them
  const VDW = { H: 1.20, C: 1.70, N: 1.55, O: 1.52, F: 1.47, P: 1.80, S: 1.80,
    Cl: 1.75, Br: 1.85, I: 1.98, Ni: 1.63, Cu: 1.40 };
  // CRC radii, which is what morfeus.ConeAngle uses by default (radii_type="crc")
  const CONE_VDW = { H: 1.10, B: 1.92, C: 1.70, N: 1.55, O: 1.52, F: 1.47, Si: 2.10, P: 1.80,
    S: 1.80, Cl: 1.75, Br: 1.85, I: 1.98, Ni: 1.97, Cu: 1.96, Pd: 2.10, Au: 2.14 };

  let GRID = null;                       // 48000 unit points, built once
  function grid() {
    if (GRID) return GRID;
    const n = 4000, dirs = new Float64Array(n * 3), gold = Math.PI * (1 + Math.sqrt(5));
    for (let i = 0; i < n; i++) {
      const phi = Math.acos(1 - 2 * (i + 0.5) / n), th = gold * (i + 0.5), s = Math.sin(phi);
      dirs[3 * i] = Math.cos(th) * s; dirs[3 * i + 1] = Math.sin(th) * s; dirs[3 * i + 2] = Math.cos(phi);
    }
    const shells = Array.from({ length: 12 }, (_, k) => 0.15 + (1.0 - 0.15) * k / 11);
    const pts = new Float64Array(12 * n * 3), w = new Float64Array(12 * n);
    let p = 0;
    shells.forEach((r, k) => {
      for (let i = 0; i < n; i++) {
        pts[p] = r * dirs[3 * i]; pts[p + 1] = r * dirs[3 * i + 1]; pts[p + 2] = r * dirs[3 * i + 2];
        w[k * n + i] = r * r; p += 3;
      }
    });
    GRID = { pts, w, wsum: w.reduce((a, b) => a + b, 0), n: 12 * n };
    return GRID;
  }

  /** centre: [x,y,z]; returns {percent, buried: Uint8Array over the grid} */
  function buriedVolume(el, X, centre, radius, opt) {
    radius = radius || 3.5;
    const scale = (opt && opt.scale) || 1.17, includeH = !(opt && opt.includeH === false);
    const g = grid(), keep = [], r2 = [];
    let skip = -1, best = 1e-6;
    X.forEach((p, i) => {
      const d = Math.hypot(p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]);
      if (d < best) { best = d; skip = i; }
    });
    X.forEach((p, i) => {
      if (i === skip || (!includeH && el[i] === 'H')) return;
      const r = (VDW[el[i]] === undefined ? 1.7 : VDW[el[i]]) * scale;
      if (Math.hypot(p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]) > radius + r) return;
      keep.push(p); r2.push(r * r);
    });
    const buried = new Uint8Array(g.n);
    if (!keep.length) return { percent: 0, buried, radius, centre, skip };
    let acc = 0;
    for (let j = 0; j < g.n; j++) {
      const x = centre[0] + radius * g.pts[3 * j], y = centre[1] + radius * g.pts[3 * j + 1],
        z = centre[2] + radius * g.pts[3 * j + 2];
      for (let a = 0; a < keep.length; a++) {
        const dx = x - keep[a][0], dy = y - keep[a][1], dz = z - keep[a][2];
        if (dx * dx + dy * dy + dz * dz < r2[a]) { buried[j] = 1; acc += g.w[j]; break; }
      }
    }
    return { percent: 100 * acc / g.wsum, buried, radius, centre, skip };
  }

  /** Exact cone angle from an apex atom (1-based), over the atoms reachable from it. */
  function coneAngle(el, X, apex1, opt) {
    const apex = X[apex1 - 1], only = opt && opt.atoms;   // 1-based subset, else everything else
    const m = [], beta = [], idx = [], inside = [];
    X.forEach((p, i) => {
      if (i === apex1 - 1) return;
      if (only && only.indexOf(i + 1) < 0) return;
      const v = [p[0] - apex[0], p[1] - apex[1], p[2] - apex[2]];
      const d = Math.hypot(v[0], v[1], v[2]);
      const r = CONE_VDW[el[i]] === undefined ? 1.7 : CONE_VDW[el[i]];
      // the apex sits inside this atom's sphere, so it subtends no cone. morfeus refuses the
      // whole calculation here; we report them instead, and a caller decides.
      if (d <= r) { inside.push(i + 1); return; }
      m.push([v[0] / d, v[1] / d, v[2] / d]);
      beta.push(Math.asin(r / d));
      idx.push(i);
    });
    if (!m.length) return null;
    const worst = (u) => {                                 // half-angle needed for axis u
      let a = -1, k = -1;
      for (let i = 0; i < m.length; i++) {
        const c = Math.max(-1, Math.min(1, u[0] * m[i][0] + u[1] * m[i][1] + u[2] * m[i][2]));
        const t = Math.acos(c) + beta[i];
        if (t > a) { a = t; k = i; }
      }
      return [a, k];
    };
    const norm = (v) => { const n = Math.hypot(v[0], v[1], v[2]) || 1; return [v[0] / n, v[1] / n, v[2] / n]; };
    // The objective is a max of smooth terms, so it has kinks: a coordinate search stalls on
    // them (0.014 deg short of morfeus). Sweep the sphere, then Nelder-Mead in (theta, phi)
    // from the best starts, which walks along a kink.
    const dir = (p) => [Math.sin(p[0]) * Math.cos(p[1]), Math.sin(p[0]) * Math.sin(p[1]), Math.cos(p[0])];
    const f = (p) => worst(dir(p))[0];
    const N = 4000, gold = Math.PI * (1 + Math.sqrt(5)), cand = [];
    for (let i = 0; i < N; i++) {
      const th = Math.acos(1 - 2 * (i + 0.5) / N), ph = (gold * (i + 0.5)) % (2 * Math.PI);
      cand.push([f([th, ph]), [th, ph]]);
    }
    cand.sort((p, q) => p[0] - q[0]);
    const nelder = (p0) => {                              // 2-D Nelder-Mead, standard coefficients
      let S = [p0, [p0[0] + 0.08, p0[1]], [p0[0], p0[1] + 0.08]].map((p) => [f(p), p]);
      for (let it = 0; it < 600; it++) {
        S.sort((a, b) => a[0] - b[0]);
        const [lo, mid, hi] = S;
        if (Math.abs(hi[0] - lo[0]) < 1e-14) break;
        const c = [(lo[1][0] + mid[1][0]) / 2, (lo[1][1] + mid[1][1]) / 2];
        const move = (t) => [c[0] + t * (c[0] - hi[1][0]), c[1] + t * (c[1] - hi[1][1])];
        const r = move(1), fr = f(r);
        if (fr < lo[0]) { const e = move(2), fe = f(e); S[2] = fe < fr ? [fe, e] : [fr, r]; }
        else if (fr < mid[0]) S[2] = [fr, r];
        else {
          const q = move(-0.5), fq = f(q);
          if (fq < hi[0]) S[2] = [fq, q];
          else S = S.map(([, p]) => { const z = [(p[0] + lo[1][0]) / 2, (p[1] + lo[1][1]) / 2]; return [f(z), z]; });
        }
      }
      S.sort((a, b) => a[0] - b[0]);
      return S[0];
    };
    let a = cand[0][0], u = dir(cand[0][1]);
    for (let i = 0; i < 12 && i < cand.length; i++) {
      const [fv, p] = nelder(cand[i][1]);
      if (fv < a) { a = fv; u = dir(p); }
    }
    const [, k] = worst(u);
    return { angle: 2 * a * 180 / Math.PI, half: a * 180 / Math.PI, axis: u, apex, inside,
      tangent: idx.filter((_, i) => Math.abs(Math.acos(Math.max(-1, Math.min(1,
        u[0] * m[i][0] + u[1] * m[i][1] + u[2] * m[i][2]))) + beta[i] - a) < 1e-6),
      worstAtom: idx[k], radii: idx.map((i) => CONE_VDW[el[i]] === undefined ? 1.7 : CONE_VDW[el[i]]), atoms: idx };
  }

  const api = { buriedVolume, coneAngle, grid, VDW, CONE_VDW };
  if (typeof module === 'object' && module.exports) module.exports = api; else root.Vbur = api;
}(typeof globalThis !== 'undefined' ? globalThis : this));
