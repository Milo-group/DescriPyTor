// The browser port must equal the Python: every value on every fixture axis.
//   python make_refs.py && node gate.js
const fs = require('fs');
const path = require('path');
const S = require('./sterimol.js');

const refs = JSON.parse(fs.readFileSync(path.join(__dirname, 'refs.json'), 'utf8'));
const KEYS = ['B1', 'B5', 'L', 'loc_B5', 'B1_B5_angle'];
let bad = 0;
for (const c of refs) {
  const mol = S.parseXYZ(fs.readFileSync(path.join(__dirname, 'fixtures', c.file), 'utf8'));
  const r = S.sterimol(mol, c.a, c.b);
  for (const k of KEYS) {
    if (r[k] !== c.ref[k]) { bad++; console.log(`${c.file} ${c.a}-${c.b} ${k}: JS ${r[k]} vs Python ${c.ref[k]}`); }
  }
  // phi is the one value the Python takes off the sampled circles, not the atom centre, so it
  // carries up to ~0.5 deg of that sampling; anything larger is a real disagreement
  if (c.phi !== undefined && Math.abs(r.phi - c.phi) > 0.6) {
    bad++; console.log(`${c.file} ${c.a}-${c.b} phi: JS ${r.phi.toFixed(3)} vs Python ${c.phi}`);
  }
}
console.log(`${refs.length} axes x ${KEYS.length + 1} values: ${bad} differ`);
process.exit(bad ? 1 : 0);
