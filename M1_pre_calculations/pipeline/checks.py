"""Structure checks every stage runs on its output, with the package's own bond rule.

A non-metal pair is bonded below BOND_SCALE (1.15) x the sum of Pyykko covalent radii, as in
``utils.help_functions.extract_connectivity``. At build time the reference is written per
molecule: ``refs/<id>.bonds`` (the non-metal bond graph) and ``refs/<id>.donors``. The same rule,
in awk (``CHECK_AWK``), runs on the cluster after each stage, with no python there:

- FAIL: a non-metal bond broken or formed (fragmentation, H transfer, ring opening);
- FAIL: two atoms closer than 0.6 x their radii sum, or a coordinate that is not a number;
- FAIL: a donor or ancillary more than 1.25 x (r_M + r_X) from its metal (it came off);
- FAIL: a chelate bite angle outside 65-105 deg;
- WARN: a new metal contact (other than the donors) inside 1.25 x (r_M + r_X), e.g. an arene
  or agostic contact: chemistry to look at, not an error.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from utils.help_functions import BOND_SCALE, GeneralConstants, _METAL_ELEMENTS

RADII = dict(GeneralConstants.PYYKKO_RADII.value)
METALS = set(_METAL_ELEMENTS)
CLASH = 0.6          # x radii sum: closer than this is two atoms on top of each other
METAL_BOND = 1.25    # x radii sum: a metal-ligand bond / contact
BITE = (65.0, 105.0)


def nonmetal_bonds(symbols, X) -> list[tuple[int, int]]:
    """1-based (i, j), i < j, of every bonded non-metal pair."""
    X = np.asarray(X, float)
    out = []
    for i in range(len(symbols)):
        if symbols[i] in METALS:
            continue
        for j in range(i + 1, len(symbols)):
            if symbols[j] in METALS:
                continue
            lim = BOND_SCALE * (RADII.get(symbols[i], 0.77) + RADII.get(symbols[j], 0.77))
            if np.linalg.norm(X[i] - X[j]) < lim:
                out.append((i + 1, j + 1))
    return out


def write_refs(workdir: Path, mid: str, symbols, X, donors_1based=(), ancillary_1based=()):
    ref = Path(workdir) / "refs"
    ref.mkdir(parents=True, exist_ok=True)
    (ref / f"{mid}.bonds").write_text("".join(f"{i} {j}\n" for i, j in nonmetal_bonds(symbols, X)),
                                      encoding="utf-8", newline="\n")
    (ref / f"{mid}.donors").write_text(" ".join(map(str, donors_1based)) + "\n", encoding="utf-8", newline="\n")
    (ref / f"{mid}.anc").write_text(" ".join(map(str, ancillary_1based)) + "\n", encoding="utf-8", newline="\n")


def check_awk() -> str:
    table = " ".join(f"r[\"{k}\"]={v};" for k, v in sorted(RADII.items()))
    metals = " ".join(f"m[\"{k}\"]=1;" for k in sorted(METALS))
    return CHECK_AWK.replace("@@RADII@@", table).replace("@@METALS@@", metals) \
        .replace("@@SCALE@@", str(BOND_SCALE)).replace("@@CLASH@@", str(CLASH)) \
        .replace("@@METALBOND@@", str(METAL_BOND)).replace("@@BITELO@@", str(BITE[0])).replace("@@BITEHI@@", str(BITE[1]))


# awk -v ref=refs/<id>.bonds -v donors="<i j>" -v anc="<k l>" -f check_structure.awk <xyz>
# prints "FAIL <reason>" / "WARN <reason>" lines; nothing when the structure is fine.
CHECK_AWK = r"""
BEGIN {
  @@RADII@@
  @@METALS@@
  scale = @@SCALE@@; clash = @@CLASH@@; mb = @@METALBOND@@
  while ((getline l < ref) > 0) { split(l, p, " "); refb[p[1] " " p[2]] = 1; nref++ }
  nd = split(donors, dn, " ")
  na = split(anc, an, " ")
}
NR == 1 { n = $1; next }
NR == 2 { next }
NR <= n + 2 {
  k = NR - 2; s[k] = $1; x[k] = $2; y[k] = $3; z[k] = $4
  if ($2 !~ /^-?[0-9.]+([eE][-+]?[0-9]+)?$/ || $3 !~ /^-?[0-9.]+([eE][-+]?[0-9]+)?$/ || $4 !~ /^-?[0-9.]+([eE][-+]?[0-9]+)?$/)
    print "FAIL coordinate is not a number on atom " k
}
function rad(e) { return (e in r) ? r[e] : 0.77 }
function dist(a, b) { return sqrt((x[a]-x[b])^2 + (y[a]-y[b])^2 + (z[a]-z[b])^2) }
END {
  if (k != n) { print "FAIL xyz has " k " atoms, header says " n; exit }
  for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) {
    d = dist(i, j); rs = rad(s[i]) + rad(s[j])
    if (d < clash * rs) print "FAIL atoms " i " and " j " are " sprintf("%.2f", d) " A apart"
    if (!(s[i] in m) && !(s[j] in m)) {
      bonded = (d < scale * rs); key = i " " j
      if (bonded && !(key in refb)) print "FAIL new bond " s[i] i "-" s[j] j " (" sprintf("%.2f", d) " A)"
      if (!bonded && (key in refb)) print "FAIL broken bond " s[i] i "-" s[j] j " (" sprintf("%.2f", d) " A)"
    }
  }
  metal = 0; for (i = 1; i <= n; i++) if (s[i] in m) { metal = i; break }
  if (metal && nd > 0) {
    for (q = 1; q <= nd; q++) { isd[dn[q]] = 1
      d = dist(metal, dn[q])
      if (d > mb * (rad(s[metal]) + rad(s[dn[q]]))) print "FAIL donor " s[dn[q]] dn[q] " is " sprintf("%.2f", d) " A from the metal"
    }
    for (q = 1; q <= na; q++) { isa[an[q]] = 1
      d = dist(metal, an[q])
      if (d > mb * (rad(s[metal]) + rad(s[an[q]]))) print "FAIL ancillary " s[an[q]] an[q] " is " sprintf("%.2f", d) " A from the metal"
    }
    if (nd == 2) {
      ax = x[dn[1]]-x[metal]; ay = y[dn[1]]-y[metal]; az = z[dn[1]]-z[metal]
      bx = x[dn[2]]-x[metal]; by = y[dn[2]]-y[metal]; bz = z[dn[2]]-z[metal]
      c = (ax*bx + ay*by + az*bz) / (sqrt(ax*ax+ay*ay+az*az) * sqrt(bx*bx+by*by+bz*bz))
      if (c > 1) c = 1; if (c < -1) c = -1
      bite = atan2(sqrt(1 - c*c), c) * 180 / 3.14159265358979
      if (bite < @@BITELO@@ || bite > @@BITEHI@@) print "FAIL bite angle " sprintf("%.1f", bite) " deg"
    }
    for (i = 1; i <= n; i++) if (i != metal && !(i in isd) && !(i in isa) && !(s[i] in m) && s[i] != "H") {
      d = dist(metal, i)
      if (d < mb * (rad(s[metal]) + rad(s[i]))) print "WARN metal contact " s[i] i " at " sprintf("%.2f", d) " A"
    }
  }
}
"""
