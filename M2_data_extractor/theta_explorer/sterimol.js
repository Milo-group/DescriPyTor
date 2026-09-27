/* Sterimol B1 / B5 / L / loc_B5 / theta, ported line for line from DescriPyTor.
 *
 * Source: M2_data_extractor/extractor_utils/sterimol_utils.py (get_sterimol_df and
 * everything it calls) and utils/help_functions.py (extract_connectivity, nob_atype,
 * CPK_RADII). Gated against the live Python on 118 axes by gate.js.
 *
 * Two details matter for bit-for-bit agreement and are reproduced deliberately:
 *   - the frame's third atom is the first element of a Python set, so PySet below
 *     emulates CPython 3.10's set layout for small ints;
 *   - transformed coordinates are rounded to 4 decimals (numpy round-half-even)
 *     before anything is measured, and the outputs are rounded again.
 */
(function (root) {
  'use strict';

  const CPK = {
    C: 1.50, C3: 1.60, 'C6/N6': 1.70, H: 1.00, N: 1.50, N4: 1.45, O: 1.35, O2: 1.35, P: 1.40,
    S: 1.70, S1: 1.00, F: 1.35, Cl: 1.80, S4: 1.40, Br: 1.95, I: 2.15, X: 1.92,
    Li: 1.82, Be: 1.53, Na: 2.27, Mg: 1.73, K: 2.75, Ca: 2.31, Rb: 3.03, Sr: 2.49, Cs: 3.43, Ba: 2.68,
    Al: 1.84, Ga: 1.87, Ge: 2.11, In: 1.93, Sn: 2.17, Tl: 1.96, Pb: 2.02, Bi: 2.07,
    Sc: 2.18, Ti: 2.11, V: 2.07, Cr: 2.06, Mn: 2.05, Fe: 2.04, Co: 2.00, Ni: 1.97, Cu: 1.96, Zn: 2.01,
    Y: 2.32, Zr: 2.23, Nb: 2.18, Mo: 2.17, Tc: 2.16, Ru: 2.13, Rh: 2.10, Pd: 2.10, Ag: 2.11, Cd: 2.18,
    Hf: 2.23, Ta: 2.22, W: 2.18, Re: 2.16, Os: 2.16, Ir: 2.13, Pt: 2.09, Au: 2.14, Hg: 2.23,
    La: 2.43, Ce: 2.42, Pr: 2.40, Nd: 2.39, Pm: 2.38, Sm: 2.36, Eu: 2.35, Gd: 2.34, Tb: 2.33, Dy: 2.31,
    Ho: 2.30, Er: 2.29, Tm: 2.27, Yb: 2.26, Lu: 2.24, Ac: 2.47, Th: 2.45, Pa: 2.43, U: 2.41, Np: 2.39,
    Pu: 2.37, Am: 2.35, Cm: 2.35,
  };
  const METALS = new Set(('Li Be Na Mg K Ca Rb Sr Cs Ba Fr Ra Al Ga Ge In Sn Sb Tl Pb Bi Sc Ti V Cr Mn Fe ' +
    'Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os Ir Pt Au Hg La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er ' +
    'Tm Yb Lu Ac Th Pa U Np Pu Am Cm').split(' '));
  const NOF = new Set(['N', 'O', 'F']);
  const HAL = new Set(['Cl', 'Br', 'F', 'I']);

  // ---------- small numeric helpers ----------
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const norm = (a) => Math.sqrt(dot(a, a));
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  // numpy.round(x, d): scale, rint (half to even), unscale
  function npRound(x, d) {
    const s = Math.pow(10, d), y = x * s, f = Math.floor(y), r = y - f;
    const k = r > 0.5 ? f + 1 : r < 0.5 ? f : (f % 2 === 0 ? f : f + 1);
    return k / s;
  }
  function solve3(A, b) {              // Cramer; A is 3x3, nonsingular here
    const det = (m) => m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
      m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    const D = det(A);
    return [0, 1, 2].map((c) => det(A.map((row, i) => row.map((v, j) => (j === c ? b[i] : v)))) / D);
  }

  // ---------- CPython 3.10 set layout, for small non-negative ints ----------
  class PySet {
    constructor() { this.mask = 7; this.table = new Array(8).fill(null); this.fill = 0; this.used = 0; }
    static of(items) { const s = new PySet(); for (const k of items) s.add(k); return s; }
    *[Symbol.iterator]() { for (const k of this.table) if (k !== null) yield k; }
    static insertClean(table, mask, key) {
      let perturb = key, i = key & mask;
      for (;;) {
        if (table[i] === null) { table[i] = key; return; }
        if (i + 9 <= mask) for (let j = 1; j <= 9; j++) if (table[i + j] === null) { table[i + j] = key; return; }
        perturb = Math.floor(perturb / 32);
        i = (i * 5 + 1 + perturb) & mask;
      }
    }
    resize(minused) {
      let size = 8;
      while (size <= minused) size <<= 1;
      const old = this.table;
      this.table = new Array(size).fill(null);
      this.mask = size - 1;
      for (const k of old) if (k !== null) PySet.insertClean(this.table, this.mask, k);
      this.fill = this.used;
    }
    add(key) {
      let perturb = key, i = key & this.mask;
      for (;;) {
        let probes = (i + 9 <= this.mask) ? 9 : 0, j = i;
        do {
          if (this.table[j] === null) {
            this.table[j] = key; this.fill++; this.used++;
            if (this.fill * 5 >= this.mask * 3) this.resize(this.used > 50000 ? this.used * 2 : this.used * 4);
            return;
          }
          if (this.table[j] === key) return;
          j++;
        } while (probes--);
        perturb = Math.floor(perturb / 32);
        i = (i * 5 + 1 + perturb) & this.mask;
      }
    }
    merge(other) {                     // set_merge
      if (other.used === 0) return;
      if ((this.fill + other.used) * 5 >= this.mask * 3) this.resize((this.used + other.used) * 2);
      if (this.fill === 0 && this.mask === other.mask && other.fill === other.used) {
        this.table = other.table.slice(); this.fill = other.fill; this.used = other.used; return;
      }
      if (this.fill === 0) {
        this.fill = other.used; this.used = other.used;
        for (const k of other) PySet.insertClean(this.table, this.mask, k);
        return;
      }
      for (const k of other) this.add(k);
    }
    union(other) { const r = new PySet(); r.merge(this); r.merge(other); return r; }
  }

  // ---------- input ----------
  function parseXYZ(text) {
    const lines = text.replace(/\r/g, '').split('\n');
    let start = 0, comment = '';
    if (/^\s*\d+\s*$/.test(lines[0] || '')) { start = 2; comment = (lines[1] || '').trim(); }
    const el = [], X = [];
    for (let i = start; i < lines.length; i++) {
      const t = lines[i].trim().split(/\s+/);
      if (t.length < 4 || !/^[A-Za-z]{1,2}$/.test(t[0])) continue;
      const xyz = t.slice(1, 4).map(Number);
      if (xyz.some((v) => !isFinite(v))) continue;
      el.push(t[0][0].toUpperCase() + t[0].slice(1).toLowerCase());
      X.push(xyz);
    }
    if (!el.length) throw new Error('No atoms found. Expected lines of: symbol x y z');
    return { el, X, comment };
  }

  // ---------- extract_connectivity (help_functions.py), literal ----------
  function connectivity(el, X, thr = 1.82, metalThr = 2.8, maxCoord = 6) {
    const n = el.length, kept = [], seen = new Set();
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        const d = norm(sub(X[i], X[j])), a1 = el[i], a2 = el[j];
        let rm = i === j;
        if ((a1 === 'H' && !NOF.has(a2)) || (a1 === 'H' && a2 === 'H') || ((a1 === 'H' || a2 === 'H') && d >= 1.5)) rm = true;
        const metal = METALS.has(a1) || METALS.has(a2);
        if (!metal) { if (d >= thr || d === 0) rm = true; } else if (d > metalThr || d === 0) rm = true;
        if ((HAL.has(a1) || HAL.has(a2)) && thr <= d && d < 2.6 && !metal) rm = false;
        if (rm) continue;
        const lo = Math.min(i, j), hi = Math.max(i, j), key = lo * n + hi;
        if (seen.has(key)) continue;       // drop_duplicates keeps the first surviving row
        seen.add(key);
        kept.push({ lo, hi, d, a1, a2 });
      }
    }
    // halogens keep their shortest bond. Keyed as the package keys them: by the
    // min-index column when the row's FIRST symbol is a halogen, max-index for the second.
    const byHal = new Map();
    kept.forEach((r, idx) => {
      if (HAL.has(r.a1)) { if (!byHal.has(r.lo)) byHal.set(r.lo, []); byHal.get(r.lo).push([idx, r.d]); }
      if (HAL.has(r.a2)) { if (!byHal.has(r.hi)) byHal.set(r.hi, []); byHal.get(r.hi).push([idx, r.d]); }
    });
    const drop = new Set();
    for (const list of byHal.values()) {
      if (list.length > 1) { list.sort((p, q) => p[1] - q[1]); list.slice(1).forEach(([idx]) => drop.add(idx)); }
    }
    const rows = kept.filter((_, idx) => !drop.has(idx));
    const isM = (r) => METALS.has(r.a1) || METALS.has(r.a2);
    const nonMetal = rows.filter((r) => !isM(r));
    const groups = new Map();
    rows.filter(isM).forEach((r) => {
      const m = METALS.has(r.a1) ? r.lo : r.hi;
      if (!groups.has(m)) groups.set(m, []);
      groups.get(m).push(r);
    });
    const metalKept = [];
    [...groups.keys()].sort((a, b) => a - b).forEach((m) => {
      const g = groups.get(m).map((r, k) => [r, k]);
      g.sort((p, q) => p[0].d - q[0].d || p[1] - q[1]);
      g.slice(0, maxCoord).forEach(([r]) => metalKept.push(r));
    });
    return nonMetal.concat(metalKept).map((r) => [r.lo + 1, r.hi + 1]);   // 1-based, like the package
  }

  // ---------- nob_atype ----------
  function atomTypes(el, bonds) {
    const nob = new Array(el.length).fill(0);
    bonds.forEach(([a, b]) => { nob[a - 1]++; nob[b - 1]++; });
    return el.map((s, i) => {
      const k = nob[i];
      if (['H', 'F', 'P', 'Cl', 'Br', 'I'].includes(s)) return s;
      if (s === 'O') return k < 1.5 ? 'O2' : 'O';
      if (s === 'S') return k < 2.5 ? 'S' : k < 5.5 ? 'S4' : 'S1';
      if (s === 'N') return k < 2.5 ? 'C6/N6' : 'N';
      if (s === 'C') return k < 2.5 ? 'C3' : k < 3.5 ? 'C6/N6' : 'C';
      if (METALS.has(s)) return s;
      return 'X';
    });
  }

  // ---------- direction_atoms_for_sterimol: the frame's third atom ----------
  function thirdAtom(bonds, origin, direction) {         // 1-based in, 1-based out
    const col = (c, v) => bonds.filter((b) => b[c] === v);
    const nbrDir = PySet.of(col(0, direction).map((b) => b[1])).union(PySet.of(col(1, direction).map((b) => b[0])));
    const nbrOrg = PySet.of(col(0, origin).map((b) => b[1])).union(PySet.of(col(1, origin).map((b) => b[0])));
    for (const a of nbrDir) if (a !== origin) return a;
    for (const a of nbrOrg) if (a !== direction) return a;
    throw new Error('Cannot fix the frame: neither atom has another neighbour to set the third direction.');
  }

  // ---------- get_molecule_connections(mode='all') + get_specific_bonded_atoms_df + filter ----------
  // `block`: 1-based atoms the walk may not pass, as metal_complex.sterimol's block= does. A
  // chelate needs it -- without it the walk leaves the substituent, goes round the backbone
  // into the other arm and picks up the ancillary, so the fragment is the whole complex.
  function fragment(bonds, n, origin, direction, block) {  // 1-based in, sorted 0-based out
    const stop = new Set(block || []);
    const adj = Array.from({ length: n + 1 }, () => []);
    bonds.forEach(([a, b]) => { adj[a].push(b); adj[b].push(a); });
    if (!adj[origin].includes(direction)) { adj[origin].push(direction); adj[direction].push(origin); }
    const reach = new Set([origin, direction]), stack = [direction];
    while (stack.length) {
      const u = stack.pop();
      for (const w of adj[u]) if (w !== origin && !stop.has(w) && !reach.has(w)) { reach.add(w); stack.push(w); }
    }
    const atoms = new Set();                             // only atoms on a real bond inside the set
    bonds.forEach(([a, b]) => { if (reach.has(a) && reach.has(b)) { atoms.add(a - 1); atoms.add(b - 1); } });
    stop.forEach((i) => atoms.delete(i - 1));
    if (!atoms.size) throw new Error('The two atoms share no bonded fragment.');
    return [...atoms].sort((a, b) => a - b);
  }

  // ---------- the frame: calc_new_base_atoms + calc_basis_vector + transform_row ----------
  function frame(X, o, d, t) {                           // 0-based
    const O = X[o], vy = sub(X[d], O), y = vy.map((v) => v / norm(vy));
    const vc = sub(X[t], O), cp = vc.map((v) => v / norm(vc.map((u) => u + 1e-8)));
    const ang = Math.acos(dot(cp, y) / (norm(cp) * norm(y)));
    const x = solve3([cp, y, cross(cp, y)], [Math.cos(ang - Math.PI / 2), 0, 0]);
    const B = [x, y, cross(x, y)];
    const F = X.map((p) => { const r = sub(p, O); return B.map((row) => npRound(dot(row, r), 4)); });
    F.basis = B; F.origin = O;                           // to map drawings back onto the molecule
    return F;
  }

  // ---------- get_b1s_list + b1s_for_loop_function ----------
  function scan(P, R) {
    const T = Array.from({ length: 100 }, (_, k) => (k === 99 ? 2 * Math.PI : k * (2 * Math.PI / 99)));
    const cloud = [];
    P.forEach(([x, z], i) => T.forEach((t) => cloud.push([x + R[i] * Math.cos(t), z + R[i] * Math.sin(t)])));
    const at = (deg) => {
      const th = deg * Math.PI / 180, c = Math.cos(th), s = Math.sin(th);
      let mx = -Infinity, nx = Infinity, my = -Infinity, ny = Infinity;
      for (const [x, z] of cloud) {
        const u = x * c - z * s, v = x * s + z * c;
        if (u > mx) mx = u; if (u < nx) nx = u; if (v > my) my = v; if (v < ny) ny = v;
      }
      const ext = [Math.abs(mx), Math.abs(nx), Math.abs(my), Math.abs(ny)];
      let k = 0;
      for (let q = 1; q < 4; q++) if (ext[q] < ext[k]) k = q;
      return { deg, B1: ext[k], n: k < 2 ? [c, -s] : [s, c] };
    };
    const argmin = (list) => list.reduce((b, r) => (r.B1 < b.B1 ? r : b));
    const sweep = Array.from({ length: 90 }, (_, i) => at(18 + i));      // 18..107 deg: the four half-extents make the period 90
    const coarse = argmin(sweep);
    const best = argmin([coarse.deg - 1, coarse.deg, coarse.deg + 1].map(at));
    best.profile = sweep.map((x) => ({ deg: x.deg, B1: x.B1, n: x.n }));   // every trial the scan made: its width and its direction
    return best;
  }

  // ---------- get_sterimol_df: the whole thing ----------
  function sterimol(mol, a, b, opt) {                    // a, b 1-based, like the package
    const { el, X } = mol, n = el.length;
    if (!(a >= 1 && a <= n && b >= 1 && b <= n) || a === b) throw new Error('Pick two different atoms, 1 to ' + n + '.');
    const bonds = mol.bonds || (mol.bonds = connectivity(el, X));
    const types = mol.types || (mol.types = atomTypes(el, bonds));
    const t = thirdAtom(bonds, a, b);
    const F = frame(X, a - 1, b - 1, t - 1);
    // `block` stops the walk at an atom; `drop` takes atoms out of the measurement even when
    // the walk reaches them (the package's drop_atoms). Either keeps a chelate or a ring from
    // looping the fragment back around the whole molecule.
    const drop = new Set((opt && opt.drop) || []);
    const idx = fragment(bonds, n, a, b, opt && opt.block).filter((i) => !drop.has(i + 1));
    if (!idx.length) throw new Error('Nothing left to measure once those atoms are dropped.');
    const R = idx.map((i) => CPK[types[i]]);
    const P = idx.map((i) => [F[i][0], F[i][2]]);
    const mag = P.map(([x, z]) => Math.sqrt(x * x + z * z));
    let k5 = 0;
    idx.forEach((_, j) => { if (R[j] + mag[j] > R[k5] + mag[k5]) k5 = j; });
    const L = Math.max(...idx.map((i, j) => F[i][1] + R[j]));
    const best = scan(P, R);
    const v = F[idx[k5]], n3 = [best.n[0], 0, best.n[1]], vn = norm(v);
    const theta = vn > 1e-10 ? Math.asin(Math.min(1, Math.max(0, Math.abs(dot(v, n3)) / vn))) * 180 / Math.PI : 0;
    const pr = Math.hypot(v[0], v[2]);
    // phi (the scan's own B1_B5_angle, which calc_sterimol overwrites with theta): the azimuth of
    // B5 from the B1 contact, so orient the normal at the contact the way the Python's sign does
    const supOf = (q) => Math.max(...idx.map((i, j) => F[i][0] * q[0] + F[i][2] * q[1] + R[j]));
    const sgn = Math.abs(supOf(best.n) - best.B1) <= Math.abs(supOf([-best.n[0], -best.n[1]]) - best.B1) ? 1 : -1;
    const cphi = pr > 0 ? Math.max(-1, Math.min(1, sgn * (v[0] * best.n[0] + v[2] * best.n[1]) / pr)) : 1;
    return {
      B1: npRound(best.B1, 4), B5: npRound(R[k5] + mag[k5], 4), L: npRound(L, 4),
      loc_B5: npRound(v[1], 4), B1_B5_angle: npRound(theta, 4),
      phi: Math.acos(cphi) * 180 / Math.PI, rho: vn > 0 ? pr / vn : 0,
      // what a drawing needs, in the Sterimol frame (y = axis, origin at atom a)
      frame: F, basis: F.basis, origin: F.origin, atoms: idx, radii: R, b5: idx[k5], normal: best.n, third: t, bonds, types, scanDeg: best.deg,
      rawB1: best.B1, profile: best.profile,
    };
  }

  // ---------- calc_dipole_gaussian (dipole_utils.py): the dipole in a declared frame ----------
  // Selections are 1-based atoms or arrays of them (centroid), as in the package:
  // origin -> centroid, Y toward the axis selection, X from the plane selection made
  // perpendicular to Y, Z = X x Y. Paper names: u = Y, v = X, w = Z.
  function dipoleFrame(X, dipole, origin, ySel, planeSel) {
    const cen = (sel) => {
      const ids = (Array.isArray(sel) ? sel : [sel]).map((i) => i - 1);
      return [0, 1, 2].map((k) => ids.reduce((s, i) => s + X[i][k], 0) / ids.length);
    };
    const unit = (v) => { const n = norm(v); return n < 1e-12 ? [0, 0, 0] : v.map((c) => c / n); };
    const O = cen(origin);
    const y = unit(sub(cen(ySel), O));
    const cp = unit(sub(cen(planeSel), O));
    let x = sub(cp, y.map((c) => c * dot(cp, y)));
    if (norm(x) < 1e-12) x = cross(Math.abs(y[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0], y);
    x = unit(x);
    const z = unit(cross(x, y));
    if (dot(x, cross(y, z)) < 0) x = x.map((c) => -c);
    const mu = dipole.slice(0, 3);
    const [cx, cy, cz] = [x, y, z].map((row) => dot(row, mu));
    return { origin: O, basis: [x, y, z], mu, dipole_x: cx, dipole_y: cy, dipole_z: cz,
      total: dipole.length > 3 ? dipole[3] : norm(mu) };
  }

  const api = { parseXYZ, connectivity, atomTypes, thirdAtom, fragment, frame, scan, sterimol, dipoleFrame, PySet, CPK, npRound };
  root.Sterimol = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);
