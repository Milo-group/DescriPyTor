/* Substructure matching, through RDKit-JS.
 *
 * Why this exists: atom numbers are per-file. The same oxazoline nitrogen is atom 12 in one
 * ligand and atom 31 in the next, so anything keyed to numbers (an axis, an alignment, a
 * dipole origin) cannot be carried across a series. A SMARTS match is keyed to the structure
 * instead, and RDKit returns its atoms in query order — which is exactly an atom mapping
 * between two different molecules.
 *
 * Bond orders: a geometry file has none. The molblock below is written with single bonds, so
 * queries are rewritten to match bonds as 'any' and aromatic atoms by element (c -> [#6]).
 * Degree queries carry what bond orders would have: a carbonyl carbon is [#6X3], its oxygen
 * [#8X1]. Turn the rewrite off to send the query through untouched.
 *
 * RDKit is loaded the first time a query is run, from ./rdkit/ beside the page if it is there,
 * otherwise from the CDN. Nothing else in the page depends on it.
 */
const SMARTS = (() => {
  const CDN = 'https://unpkg.com/@rdkit/rdkit@2026.3.6/dist/';
  const METAL = new Set(['Li', 'Na', 'K', 'Rb', 'Cs', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
    'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Al', 'Ga', 'In', 'Sn', 'Tl', 'Pb', 'Bi', 'La', 'Ce', 'Sm', 'Eu', 'Yb', 'U']);
  const AROM = { c: '[#6]', n: '[#7]', o: '[#8]', s: '[#16]', p: '[#15]', b: '[#5]' };

  let RD = null, pending = null;

  const script = (src) => new Promise((ok, no) => {
    const s = document.createElement('script');
    s.src = src; s.onload = () => ok(src); s.onerror = () => no(new Error('could not load ' + src));
    document.head.appendChild(s);
  });

  function ready() {                         // one load, however many queries
    if (RD) return Promise.resolve(RD);
    if (!pending) {
      // a page opened from a file:// URL cannot fetch the .wasm beside it, so the CDN goes first there
      const bases = location.protocol === 'file:' ? [CDN, './rdkit/'] : ['./rdkit/', CDN];
      pending = (async () => {
        let last = null;
        for (const base of bases) {
          try {
            if (typeof initRDKitModule === 'undefined') await script(base + 'RDKit_minimal.js');
            RD = await initRDKitModule({ locateFile: (f) => base + f });
            return RD;
          } catch (e) { last = e; }
        }
        pending = null;
        throw new Error('RDKit could not be loaded (' + (last && last.message) + '). Put RDKit_minimal.js and .wasm in a folder "rdkit" beside this page to work offline.');
      })();
    }
    return pending;
  }

  // the structure as RDKit sees it: every atom, every bond except the ones to a metal (so the
  // organic part still perceives rings and valences), single bonds throughout
  function molblock(el, X, bonds) {
    const keep = bonds.filter(([a, b]) => !METAL.has(el[a - 1]) && !METAL.has(el[b - 1]));
    const n3 = (v) => String(v).padStart(3);
    const c10 = (v) => v.toFixed(4).padStart(10);
    const deg = el.map(() => 0);
    keep.forEach(([a, b]) => { deg[a - 1]++; deg[b - 1]++; });
    const out = ['', '  theta-explorer', '', `${n3(el.length)}${n3(keep.length)}  0  0  0  0  0  0  0  0999 V2000`];
    // the valence field is the atom's own bond count (0 is written 15), so RDKit adds no implicit
    // hydrogens: every H is already in the file, and X3 then means three real neighbours
    el.forEach((s, i) => out.push(`${c10(X[i][0])}${c10(X[i][1])}${c10(X[i][2])} ${(s + '   ').slice(0, 3)} 0`
      + `${n3(0)}${n3(0)}${n3(0)}${n3(0)}${n3(deg[i] || 15)}${n3(0)}${n3(0)}${n3(0)}${n3(0)}${n3(0)}${n3(0)}`));
    keep.forEach(([a, b]) => out.push(`${n3(a)}${n3(b)}  1  0  0  0  0`));
    out.push('M  END');
    return out.join('\n');
  }

  // bond orders out, aromatic atoms down to elements: what a query has to look like against geometry
  function anyBonds(q) {
    let out = '', depth = 0, prev = '';
    for (const ch of q) {
      if (ch === '[') depth++;
      else if (ch === ']') depth--;
      if (!depth && ch !== '[' && ch !== ']') {
        if ('=#:-'.includes(ch)) { out += '~'; prev = ch; continue; }
        if (AROM[ch] && !/[A-Z]/.test(prev)) { out += AROM[ch]; prev = ch; continue; }   // not the l of Cl or the r of Br
      }
      out += ch; prev = ch;
    }
    return out;
  }

  // -> {matches: [[0-based atom, ...], ...]}; throws with something a chemist can act on
  async function match(mol, query, ignoreBondOrders = true) {
    const rd = await ready();
    const q = (query || '').trim();
    if (!q) return { matches: [], query: q };
    const used = ignoreBondOrders ? anyBonds(q) : q;
    const qm = rd.get_qmol(used);
    if (!qm || !qm.is_valid()) { if (qm) qm.delete(); throw new Error(`"${used}" is not a SMARTS pattern RDKit understands`); }
    const m = rd.get_mol(molblock(mol.el, mol.X, mol.bonds));
    if (!m || !m.is_valid()) { qm.delete(); if (m) m.delete(); throw new Error('RDKit could not read this structure'); }
    let hits = [];
    try {
      const raw = m.get_substruct_matches(qm);          // "{}" when nothing matched
      const parsed = raw ? JSON.parse(raw) : [];
      hits = Array.isArray(parsed) ? parsed : [];
    } finally { qm.delete(); m.delete(); }
    // one match per set of atoms: the three methyls of a tert-butyl are one match, not six
    const seen = new Set(), out = [];
    hits.map((h) => h.atoms).filter(Array.isArray).forEach((a) => {
      const key = a.slice().sort((p, q) => p - q).join(',');
      if (!seen.has(key)) { seen.add(key); out.push(a); }
    });
    return { matches: out, query: used };
  }

  const api = { match, molblock, anyBonds, ready, loaded: () => !!RD };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  return api;
})();
