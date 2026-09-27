// The browser port of %Vbur and the cone angle must equal the Python / morfeus.
//   <mol_ml python> make_pos_refs.py && node gate_pos.js
const fs = require('fs');
const path = require('path');
const S = require('./sterimol.js');
const V = require('./vbur.js');

const refs = JSON.parse(fs.readFileSync(path.join(__dirname, 'pos_refs.json'), 'utf8'));
let bad = 0, worstV = 0, worstC = 0;
for (const c of refs) {
  const m = S.parseXYZ(fs.readFileSync(path.join(__dirname, 'fixtures', c.file), 'utf8'));
  if (c.kind === 'vbur') {
    const got = V.buriedVolume(m.el, m.X, m.X[c.atom - 1], c.radius).percent;
    const d = Math.abs(got - c.ref);
    worstV = Math.max(worstV, d);
    if (d > 1e-6) { bad++; console.log(`${c.file} vbur atom ${c.atom} r=${c.radius}: JS ${got.toFixed(6)} vs Python ${c.ref}`); }
  } else {
    const r = V.coneAngle(m.el, m.X, c.atom);
    const d = Math.abs(r.angle - c.ref);
    worstC = Math.max(worstC, d);
    const tan = r.tangent.map((i) => i + 1).sort((a, b) => a - b).join(',');
    if (d > 0.01) { bad++; console.log(`${c.file} cone at ${c.atom}: JS ${r.angle.toFixed(4)} vs morfeus ${c.ref}`); }
    else if (c.tangent.length && tan !== c.tangent.join(',')) {
      bad++; console.log(`${c.file} cone at ${c.atom}: tangent atoms JS [${tan}] vs morfeus [${c.tangent}]`);
    }
  }
}
console.log(`${refs.length} values: ${bad} differ  (worst %Vbur ${worstV.toExponential(1)}, worst cone ${worstC.toFixed(4)} deg)`);
process.exit(bad ? 1 : 0);
