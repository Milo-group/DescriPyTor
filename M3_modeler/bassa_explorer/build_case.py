"""Data for the BASSA model explorer.

Two steps per prepared case, because the sampler and the structure readers live in different environments
(inside the package the Python that runs the server has to have bassa_reg; structures are read by BASSA_EXPLORER_READER):
    <bassa_python>  build_case.py models <case>    rerun the corrected sampler on the training compounds, write <case>_models.json
    <reader_python> build_case.py page <case>      add structures, atom roles and the two reference lines, write <case>.js
Inside the DescriPyTor package: the page (index.html) is in this folder, the cases are in settings.CASES_DIR, and the molecule
view is the package's feature explorer in its figure mode. `models`, `page` and the paper parts of `attach` rebuild the prepared
cases of the DescriPyTor paper and need its folder (BASSA_EXPLORER_PAPER).
Cases: cs2, cs3_cp. Nothing here is a result: it shows the same run as the summary PDF (engine.fit, seed 7, 20,000 iterations).
"""
import json, re, sys
from pathlib import Path
import numpy as np, pandas as pd

from settings import S, CASES_DIR as C                                         # C: where the case files and manifest.js are

H = Path(__file__).resolve().parent; BASSA = Path(S["engine"]); V5 = Path(S["paper"]); TOP = 400
PKG = Path(S["descripytor"]); EXPLORER = PKG / "M2_data_extractor" / "feature_explorer"
FEATHERS = Path(S["cs2_feathers"]); CUCL = V5 / "zenodo_archive/geometries/case_study_3_corminboeuf/single_structure_xtb/CuCl"
CASES = {"cs2": dict(title="Cobalt C–H annulation · Ackermann (Dana et al. 2025)", unit="kcal/mol", response="ΔΔG‡ (more negative = more selective)"),
         "cs3_cp": dict(title="Copper cyclopropanation · Corminboeuf and Waser (Schoepfer et al. 2024)", unit="kcal/mol", response="ΔΔG‡ (higher = more selective)"),
         "wu_doyle": dict(title="Nickel Suzuki coupling of acetals · Wu and Doyle (2017)", unit="%", response="yield")}
WD = V5 / "notes/wu_doyle_cs"; SEEDS = (7, 11, 13)                              # the first seed is the run that is shown; the others only give the spread of the frequencies
LABELS = {"dipole_metal_y": "dipole along metal→donor", "dipole_metal_x": "dipole in-plane ⟂ metal→donor",
          "dipole_metal_z": "dipole out-of-plane ⟂ metal→donor", "dipole_total": "dipole magnitude", "vmin": "Vmin", "vbur_percent": "%Vbur", "cone_angle": "cone angle"}
LOOKS = {"wu_doyle": dict(hydrogens="none", labelScale=1.2, muPerD=5.0, vdwOpacity=0.08)}
SHOW = {"B1": "b1", "B5": "b5", "loc_B5": "loc", "theta": "theta", "phi": "theta", "B1_B5_angle": "theta", "L": "L",
        "loc_B1": "b1", "B1_tangent_atoms": "b1", "B1_B5_azimuth": "theta", "B1_B5_acute": "theta"}          # the last three: columns the extractor writes beside the classic five
ROLE = {"M": "M", "D": "D", "Cs": "C*", "R": "R", "Rp": "R′"}; PAIRS = {"MCs": "M Cs", "DCs": "D Cs", "CsR": "Cs R", "MD": "M D", "MR": "M R", "DR": "D R", "fromM": "M D", "sub": "Cs R", "RRp": "R Rp"}
A = r"[A-Za-z]*(\d+)"                                                                   # an atom in a column name: element letter and number
PI_MARKS = [(rf"pi_share\[{A}\|{A}=[A-Za-z]*\d+\]", "pi_share", "π on atom ", 1, " %"), (rf"pistar_occ\[{A}={A}\]", "pistar_occ", "π* occ ", 3, " e"),
            (rf"E2_pi_to_aryl\[{A}={A}\|{A}\]", "E2_out", "E2 π→aryl ", 1, " kcal/mol"), (rf"E2_aryl_to_pistar\[{A}={A}\|{A}\]", "E2_in", "E2 aryl→π* ", 1, " kcal/mol"),
            (rf"stretch\[{A}[=-]{A}\]", "stretch", "ν ", 0, " cm⁻¹"), (rf"ring_hi\[{A}\]", "ring_hi", "ring ν ", 0, " cm⁻¹"), (rf"ring_lo\[{A}\]", "ring_lo", "ring ν′ ", 0, " cm⁻¹"),
            (rf"E2_conjugation\[{A}-{A}\]", "E2_conj", "E2 conj. ", 1, " kcal/mol"), (rf"E2_O_to_pistar\[{A}\|{A}={A}\]", "E2_O", "E2 O→π* ", 1, " kcal/mol"), (rf"polarizability\[{A}->{A}\]", "alpha", "α∥ ", 0, " bohr³")]
PI_NAME = {"pi_share": "share of the π bond on atom %d (bond to %d)", "pistar_occ": "π* occupancy, %d=%d (NBO)", "E2_out": "donation π(%d=%d) → aryl on %d (NBO E2)", "E2_in": "donation aryl → π*(%d=%d), aryl on %d (NBO E2)",
           "E2_conj": "conjugation across %d–%d (NBO E2, both directions)", "stretch": "stretching frequency of the %d–%d bond",
           "ring_hi": "ring stretch, higher mode (ring at atom %d)", "ring_lo": "ring stretch, lower mode (ring at atom %d)", "E2_O": "donation lone pairs of atom %d → π*(%d=%d) (NBO E2)", "alpha": "polarizability along %d→%d"}
PI_SAY = {"pi_share": "The share of the π bonding orbital that sits on atom %d, in percent (its partner is atom %d): how the π bond is polarized.",
          "pistar_occ": "The number of electrons in the π antibonding orbital between atoms %d and %d. Zero for an isolated double bond; it fills as neighbours donate into it.",
          "E2_out": "The NBO second-order energy of donation from the π bond %d=%d into the π* orbitals of the aryl group attached at atom %d, summed. Larger means the double bond gives more to the ring.",
          "E2_in": "The NBO second-order energy of donation from the aryl group's π orbitals and lone pairs into the π* of %d=%d (aryl attached at atom %d), summed.",
          "stretch": "The frequency of the vibration in which atoms %d and %d move against each other along their bond. A stiffer, more localized bond vibrates at a higher frequency.",
          "ring_hi": "The higher of the two ring C–C stretching modes of the ring that holds atom %d (for a benzene ring, the pair near 1650 cm⁻¹). It does not depend on how the ring atoms are numbered.",
          "ring_lo": "The lower of the two ring C–C stretching modes of the ring that holds atom %d. It does not depend on how the ring atoms are numbered.",
          "E2_conj": "The sum of the two donation energies across the %d–%d bond: how strongly the double bond and the aryl ring are conjugated.",
          "E2_O": "The NBO second-order energy of donation from the lone pairs of atom %d into the π* of %d=%d.", "alpha": "The polarizability of the molecule along the direction from atom %d to atom %d."}


# the vibration columns as the extractor names them (Frequency_Stretch_1_2, Frequency_Bending_5-6, cross_9_0, ...)
PI_MARKS += [(r"Frequency_Stretch_(\d+)_(\d+)", "x_stretch", "ν ", 0, " cm⁻¹"), (r"[A-Za-z]+_Stretch_(\d+)_(\d+)", "x_stretch_amp", "stretch ", 3, ""),
             (r"Frequency_Bending_(\d+)-(\d+)", "x_bend", "δ ", 0, " cm⁻¹"), (r"Cross_mag_Bending_(\d+)-(\d+)", "x_bend_mag", "bend ", 3, ""),
             (r"cross_angle_(\d+)(?:_\d+)?", "x_cross_angle", "cross ∠ ", 1, "°"), (r"para_angle_(\d+)(?:_\d+)?", "x_para_angle", "para ∠ ", 1, "°"),
             (r"cross_(\d+)(?:_\d+)?", "x_cross", "ring ν cross ", 0, " cm⁻¹"), (r"para_(\d+)(?:_\d+)?", "x_para", "ring ν para ", 0, " cm⁻¹")]
PI_NAME.update(x_stretch="stretching frequency of the %d–%d bond", x_stretch_amp="stretch of the %d–%d bond, second value", x_bend="bending frequency, atoms %d and %d",
               x_bend_mag="bending mode of atoms %d and %d, size of the motion", x_cross="ring vibration, cross mode (ring at atom %d)", x_para="ring vibration, para mode (ring at atom %d)",
               x_cross_angle="ring cross mode, angle to the ring plane (ring at atom %d)", x_para_angle="ring para mode, angle to the ring axis (ring at atom %d)")
PI_SAY.update(x_stretch="The frequency of the vibration in which atoms %d and %d move against each other along their bond.",
              x_stretch_amp="The second value the extractor writes for the stretching mode of the %d–%d bond (its intensity or amplitude).",
              x_bend="The frequency of the bending mode in which atoms %d and %d swing about the centre they share.",
              x_bend_mag="How large the bending motion of atoms %d and %d is in the mode the extractor chose.",
              x_cross="The frequency of the ring mode the extractor calls cross, for the ring that contains atom %d.", x_para="The frequency of the ring mode the extractor calls para, for the ring that contains atom %d.",
              x_cross_angle="The angle the extractor reports for the cross ring mode of the ring that contains atom %d.", x_para_angle="The angle the extractor reports for the para ring mode of the ring that contains atom %d.")
LABELS.update(iso="isotropic polarizability", aniso="polarizability anisotropy", energy="energy")


def drawing(c):
    """what the explorer should draw for a column name; None when the name is not one it can place.
    Platform names carry atom numbers. Metal-complex names carry roles (metal, donor, stereocentre, substituent) that each compound resolves itself."""
    m = re.fullmatch(r"ster_(\w+)\[[A-Za-z]*(\d+)->[A-Za-z]*(\d+)\]", c) or re.fullmatch(r"(B1|B5|L|loc_B5|loc_B1|B1_tangent_atoms|B1_B5_angle|B1_B5_azimuth|B1_B5_acute)_(\d+)-(\d+)", c)
    if m and m.group(1) in SHOW: return dict(kind="sterimol", a=int(m.group(2)), b=int(m.group(3)), flag=SHOW[m.group(1)])
    m = re.fullmatch(r"(%s)_(theta|phi|B1|B5|L|loc_B5)_(sym|asym)" % "|".join(PAIRS), c)
    if m: a, b = PAIRS[m.group(1)].split(); return dict(kind="sterimol", ra=a, rb=b, flag=SHOW[m.group(2)], arms=m.group(3))
    if c in SHOW and c != "phi": return dict(kind="sterimol", ra="M", rb="D", flag=SHOW[c], arms="one")
    m = re.fullmatch(r"dipole_metal_([xyz])", c)
    if m: return dict(kind="dipole", comp=dict(y="u", x="v", z="w")[m.group(1)], ro="M", ry="D", rp="Pl")
    if c == "dipole_total": return dict(kind="dipole_total", ro="M", ry="D", rp="Pl")
    m = re.fullmatch(r"q_(cm5|nbo|hirshfeld)\[[A-Za-z]*(\d+)\]", c) or re.fullmatch(r"(cm5|nbo|hirshfeld)_atom_(\d+)", c)
    if m and int(m.group(2)) >= 1: return dict(kind="charge", i=int(m.group(2)), type=m.group(1))      # atoms are numbered from 1
    m = re.fullmatch(r"mu_([xyz])\[[A-Za-z]*(\d+)->[A-Za-z]*(\d+)\|[A-Za-z]*(\d+)\]", c)
    if m: return dict(kind="dipole", comp=dict(y="u", x="v", z="w")[m.group(1)], origin=int(m.group(2)), y=int(m.group(3)), plane=int(m.group(4)))   # y runs origin -> y atom (the explorer's u)
    m = re.fullmatch(r"dipole_([xyz])_(\d+)-(\d+)-(\d+)", c)                 # the extractor's own names: dipole_x_1-2-19 is origin 1, y-axis atom 2, plane atom 19
    if m: return dict(kind="dipole", comp=dict(y="u", x="v", z="w")[m.group(1)], origin=int(m.group(2)), y=int(m.group(3)), plane=int(m.group(4)))
    m = re.fullmatch(r"total_dipole_(\d+)-(\d+)-(\d+)", c)
    g = re.fullmatch(r"(dipole_[xyz]|total_dipole)_\{([\d,]+)\}-(\d+)-(\d+)", c)        # the origin as the centroid of several atoms: dipole_x_{23,24,25,26,27,28}-23-1
    if g:
        o = ", ".join(g.group(2).split(","))
        if g.group(1) == "total_dipole": return dict(kind="dipole_total", origin=o, y=int(g.group(3)), plane=int(g.group(4)))
        return dict(kind="dipole", comp=dict(y="u", x="v", z="w")[g.group(1)[-1]], origin=o, y=int(g.group(3)), plane=int(g.group(4)))
    if m: return dict(kind="dipole_total", origin=int(m.group(1)), y=int(m.group(2)), plane=int(m.group(3)))
    m = re.fullmatch(r"(?:\d+_)?(bond_length|angle|dihedral)_\[?([\d, -]+)\]?", c)
    if m: return dict(kind=m.group(1), atoms=[int(v) for v in re.findall(r"\d+", m.group(2))])
    m = re.fullmatch(r"pi_occ\[[A-Za-z]*(\d+)=[A-Za-z]*(\d+)\]", c)
    if m: return dict(kind="mark", what="pi_occ", atoms=[int(m.group(1)), int(m.group(2))], prefix="π occ ", dec=3, unit=" e")
    for pat, what, prefix, dec, unit in PI_MARKS:
        m = re.fullmatch(pat, c)
        if m: return dict(kind="mark", what=what, atoms=[int(v) for v in m.groups()], prefix=prefix, dec=dec, unit=unit)
    m = re.fullmatch(r"atomic_dipole\[[A-Za-z]*(\d+)\|[A-Za-z]*(\d+)->[A-Za-z]*(\d+)\]", c)
    if m: return dict(kind="mark", what="atomic_dipole", at=int(m.group(1)), atoms=[int(m.group(2)), int(m.group(3))], prefix="μ(%s) " % m.group(1), dec=3, unit=" a.u.")
    m = re.fullmatch(r"(cm5|nbo|hirshfeld)_diff_(\d+)-(\d+)", c)
    if m: return dict(kind="charge_diff", type=m.group(1), i=int(m.group(2)), j=int(m.group(3)))
    m = re.fullmatch(r"mu_([xyz])\[ring->[A-Za-z]*(\d+)\]", c)          # shortcut: the ring is taken as the six atoms numbered from the named one (true for CS1's aryl ring, 23 to 28); a table numbered otherwise needs the ring's atoms in the name
    if m: n = int(m.group(2)); return dict(kind="dipole", comp=dict(y="u", x="v", z="w")[m.group(1)], origin="%d-%d" % (n, n + 5), y=n, plane=1)
    return {"mu_total": dict(kind="dipole_total"), "q_metal": dict(kind="atoms", roles=["M"]), "bite": dict(kind="atoms", roles=["D", "M"]), "MD_mean": dict(kind="atoms", roles=["M", "D"]), "MD_asym": dict(kind="atoms", roles=["M", "D"]),
            "a_R_Cstereo_N_sym": dict(kind="atoms", roles=["R", "Cs", "D"], angle=True), "a_R_Cstereo_N_asym": dict(kind="atoms", roles=["R", "Cs", "D"], angle=True),
            "q_donor_sym": dict(kind="atoms", roles=["D"]), "q_donor_asym": dict(kind="atoms", roles=["D"]), "wbo_MD_sym": dict(kind="atoms", roles=["M", "D"]), "wbo_MD_asym": dict(kind="atoms", roles=["M", "D"]),
            "wbo_metal_tot": dict(kind="atoms", roles=["M"]), "q_anc_sum": dict(kind="atoms", roles=["A"])}.get(c)


ARMS = {"sym": " · mean of the two arms", "asym": " · difference between the arms"}
GRAPH = {"kappa1": "Kier shape κ1", "kappa2": "Kier shape κ2", "zagreb2": "Zagreb M2", "wiener": "Wiener", "estrada": "Estrada", "global_eff": "global efficiency"}
PLAIN = {"lumo": "LUMO energy", "homo": "HOMO energy", "iso": "polarizability, isotropic", "aniso": "polarizability, anisotropy", "mu_bisector": "dipole along the donor bisector (M, D, D frame)",
         "mu_desym": "dipole across the two donors (M, D, D frame)", "mu_outofplane": "dipole out of the M, D, D plane", "q_spread": "spread of the atomic charges", "q_absmax": "largest atomic charge (absolute)",
         "q_metal": "charge on the metal", "bite": "bite angle D–M–D", "MD_mean": "metal–donor distance · mean of the two arms", "MD_asym": "metal–donor distance · difference between the arms",
         "wbo_metal_tot": "total bond order at the metal (Wiberg)", "q_anc_sum": "summed charge on the co-ligands of the metal"}
ARMED = {"a_R_Cstereo_N": "angle R–C*–D", "q_donor": "charge on the donor atoms", "wbo_MD": "metal–donor bond order (Wiberg)"}


def label(c, w):
    """a readable name for the table; `describe` gives the sentence for the hover note"""
    if c in LABELS: return LABELS[c]
    q = lambda s: s.replace("loc_B5", "loc B5").replace("loc_B1", "loc B1").replace("B1_tangent_atoms", "atoms touching the B1 plane,").replace("B1_B5_angle", "B1–B5 angle").replace("B1_B5_azimuth", "B1–B5 azimuth").replace("B1_B5_acute", "B1–B5 acute angle").replace("theta", "θ").replace("phi", "φ")
    if c in PLAIN: return PLAIN[c]
    m = re.fullmatch(r"(%s)_(sym|asym)" % "|".join(ARMED), c)
    if m: return ARMED[m.group(1)] + ARMS[m.group(2)]
    m = re.fullmatch(r"vbur_(noH_)?([\d.]+)", c)
    if m: return "%%Vbur at the metal, %s Å sphere%s" % (m.group(2), ", no H" if m.group(1) else "")
    m = re.fullmatch(r"tf_(\w+)", c)
    if m: g = m.group(1); k = re.fullmatch(r"chi(\d)(v?)", g); return "ligand graph: " + ("%sconnectivity χ%s" % ("valence " if k.group(2) else "", k.group(1)) if k else GRAPH.get(g, g))
    if not w: return c
    if w["kind"] == "sterimol" and "ra" in w: return "%s on %s→%s%s" % (q(re.sub(r"^(%s)_|_(sym|asym)$" % "|".join(PAIRS), "", c)), ROLE[w["ra"]], ROLE[w["rb"]], ARMS.get(w["arms"], ""))
    if w["kind"] == "sterimol": return "%s on %d→%d" % (q(re.sub(r"^ster_|\[.*|_\d+-\d+$", "", c)), w["a"], w["b"])
    if w["kind"] == "charge": return "%s charge, atom %d" % (QN[w["type"]], w["i"])
    if w["kind"] == "mark" and w["what"] in PI_NAME: return PI_NAME[w["what"]] % tuple(w["atoms"])
    if w["kind"] == "mark": return "π-bond occupancy, %d=%d (NBO)" % tuple(w["atoms"]) if w["what"] == "pi_occ" else "atomic dipole of atom %d, along %d→%d (Hirshfeld)" % (w["at"], *w["atoms"])
    if w["kind"] == "charge_diff": return "%s charge, atom %d minus atom %d" % (QN[w["type"]], w["i"], w["j"])
    if w["kind"] == "dipole": return "dipole %s %s→%d" % (dict(u="along", v="in-plane ⟂", w="out-of-plane ⟂")[w["comp"]], ("ring %s centre" % w["origin"].replace("-", "–") if "," not in w["origin"] else "centre of atoms %s" % w["origin"]) if isinstance(w["origin"], str) else w["origin"], w["y"])
    if w["kind"] == "dipole_total": return "dipole magnitude"
    return w["kind"].replace("_", " ") + " " + "–".join(map(str, w["atoms"]))


QN = dict(cm5="CM5", nbo="NBO", hirshfeld="Hirshfeld")
WORDS = "M: the metal. D: a donor atom bound to it. C*: the stereocentre next to that donor. R: the first atom of the substituent on C*. R′: the next heavy atom along that substituent. The chelate has two arms; a column is their mean or the difference between them."
STER = {"B1": "smallest width (Sterimol B1)", "B5": "largest width (Sterimol B5)", "L": "length (Sterimol L)", "loc_B5": "distance along the axis at which the widest point sits (loc B5)", "theta": "angle θ between the B1 direction and the widest point",
        "phi": "angle φ between the B1 and B5 directions", "B1_B5_angle": "angle between the B1 and B5 directions",
        "loc_B1": "distance along the axis at which the smallest width sits (loc B1)", "B1_tangent_atoms": "number of atoms that touch the B1 plane", "B1_B5_azimuth": "azimuth of the B5 direction, measured round the axis from the B1 direction",
        "B1_B5_acute": "acute angle between the B1 and B5 directions"}


def describe(c, w):
    """one or two sentences for the hover note of a descriptor: what is measured, on which atoms; the column name last"""
    t = label(c, w); col = " Column: %s." % c
    if w and w["kind"] == "sterimol":
        k = next((v for n, v in sorted(STER.items(), key=lambda x: -len(x[0])) if n in c), "Sterimol value")
        if "ra" in w: return "The %s, measured along the axis from %s to %s, over what lies beyond the second atom%s. %s%s" % (k, ROLE[w["ra"]], ROLE[w["rb"]], "" if w["arms"] == "one" else ARMS[w["arms"]].replace(" · ", "; the "), WORDS if w["arms"] != "one" else "", col)
        return "The %s, measured along the axis from atom %d to atom %d, over everything atom %d reaches without passing through atom %d.%s" % (k, w["a"], w["b"], w["b"], w["a"], col)
    if w and w["kind"] == "dipole":
        how = dict(u="along", v="perpendicular to, in the plane of,", w="perpendicular to, out of the plane of,")[w["comp"]]
        if "ro" in w: return "The component of the molecular dipole %s the metal to donor axis. %s%s" % (how, WORDS, col)
        start = ("the centre of ring %s" if "," not in w["origin"] else "the centre of atoms %s") % w["origin"] if isinstance(w["origin"], str) else "atom %d" % w["origin"]
        return "The component of the molecular dipole %s the axis from %s to atom %d; atom %d fixes the plane.%s" % (how, start, w["y"], w["plane"], col)
    if w and w["kind"] == "mark" and w["what"] in PI_SAY: return PI_SAY[w["what"]] % tuple(w["atoms"]) + col
    if w and w["kind"] == "mark" and w["what"] == "pi_occ": return "The number of electrons in the π bonding orbital between atoms %d and %d (the NBO labelled BD(2)). A full, localized π bond holds 2; a lower value means more of it has leaked into the neighbouring π system, so a more conjugated double bond.%s" % (*w["atoms"], col)
    if w and w["kind"] == "mark": return "The Hirshfeld atomic dipole of atom %d (how lopsided the electron density is around that atom), taken along the direction from atom %d to atom %d. In atomic units.%s" % (w["at"], *w["atoms"], col)
    if c.startswith("tf_"): return "%s: an index of how the free ligand's atoms are connected, divided by its number of heavy atoms. It carries size and branching, not three-dimensional shape; nothing is drawn for it.%s" % (t, col)
    if c.startswith("vbur"): return "%s: the share of a sphere around the metal that the ligand fills.%s" % (t, col)
    if w and w.get("roles"): return "%s. %s%s" % (t, WORDS, col)
    return t + "." + col


def scheme(case):
    """the reaction scheme of a prepared case as a data address, so it travels inside the case file"""
    import base64
    f = BASSA / "summary_figs" / ("reaction_%s.png" % {"cs3_cc_sym": "cs3_cc"}.get(case, case))
    return "data:image/png;base64," + base64.b64encode(f.read_bytes()).decode() if f.exists() else None


def more_roles(xyz, roles):
    """R′ (the next heavy atom along each arm's substituent) and the metal's co-ligands, added to the roles of a chelate; 1-based, None for an arm without one"""
    import tempfile
    sys.path.insert(0, str(PKG)); from M2_data_extractor import metal_complex as mc
    with tempfile.TemporaryDirectory() as t: f = Path(t) / "m.xyz"; f.write_text(xyz); c = mc.MetalComplex.from_xyz(str(f))
    rp = []
    for D, Cs, R in zip(roles["D"], roles["Cs"], roles["R"]):
        # shortcut: the package also leaves out atoms of the donor's ring; here only the donor and the metal are left out, which differs only when R sits on that ring
        heavy = [k for k in c.adj[R - 1] if k != Cs - 1 and c.symbols[k] != "H" and k not in (D - 1, c.metal)]; rp.append(max(heavy, key=lambda k: len(c.adj[k])) + 1 if heavy else None)
    return dict(roles, Rp=rp, A=[k + 1 for k in c.adj[c.metal] if k + 1 not in roles["D"]])



def folds(X, y, E, iters=20000):
    """the sampler rerun without each compound in turn: its prediction for that compound, and the descriptors' frequencies in that fold"""
    pred, freq = np.zeros(len(y)), []
    for i in range(len(y)):
        k = np.arange(len(y)) != i; f = E.fit(X[k], y[k], xt=X[~k], full=X, iters=iters); pred[i] = float(np.mean(f["pred"])); freq.append(f["freq"].reindex(X.columns).fillna(0.0).to_numpy())
    yv = y.to_numpy(float); F = np.array(freq); top = (np.argsort(-F, axis=1)[:, :3][:, :, None] == np.arange(X.shape[1])).any(1).sum(0)
    return dict(nested_q2_here=float(1 - ((yv - pred) ** 2).sum() / ((yv - yv.mean()) ** 2).sum()), fold_top=top.tolist(), n_folds=len(y),
                fold_freq=np.round(F, 3).tolist())                             # row = the compound left out, column = descriptor: who each frequency rests on


def compute(X, y, iters=20000, case="uploaded", unit="", response="response", title="", seeds=SEEDS[:1], with_folds=False):
    """one corrected-sampler run on a table: inclusion frequencies, the most visited models, what appears with what.
    More seeds give the spread of each frequency between chains; with_folds adds the leave-one-out reruns."""
    sys.path.insert(0, str(BASSA)); import engine as E
    X0 = X; f = E.fit(X, y, iters=iters, seed=seeds[0]); s = f["s"].to_numpy(); names = list(f["s"].index); X = X[names]; n = len(y)
    assert np.allclose(s.mean(1), f["freq"].to_numpy())
    Z = ((X - X.mean()) / X.std(ddof=0)).to_numpy(); yv = y.to_numpy(float); sst = ((yv - yv.mean()) ** 2).sum()
    keys, counts = np.unique(s.T, axis=0, return_counts=True); order = np.argsort(-counts); out = []
    for k in order[:TOP]:
        idx = np.flatnonzero(keys[k]); A = np.c_[np.ones(n), Z[:, idx]]; b = np.linalg.lstsq(A, yv, rcond=None)[0]; r = yv - A @ b
        h = np.einsum("ij,jk,ik->i", A, np.linalg.pinv(A.T @ A), A)                              # leverages: leave-one-out residual = r / (1 - h)
        out.append(dict(cols=idx.tolist(), share=float(counts[k] / s.shape[1]), r2=float(1 - (r ** 2).sum() / sst), q2=float(1 - ((r / np.clip(1 - h, 1e-6, None)) ** 2).sum() / sst)))   # a compound the model fits exactly (leverage 1) has no leave-one-out residual: the clip keeps Q2 finite and very low
    joint = (s.astype(float) @ s.T.astype(float)) / s.shape[1]; assert np.allclose(np.diag(joint), s.mean(1))
    lig = [str(i) for i in X.index]; draw = [drawing(c) for c in names]
    more = np.array([f["freq"].to_numpy()] + [E.fit(X0, y, iters=iters, seed=sd)["freq"].reindex(names).fillna(0.0).to_numpy() for sd in seeds[1:]])
    extra = dict(freq_lo=more.min(0).tolist(), freq_hi=more.max(0).tolist(), chains=len(seeds)); extra.update(folds(X, y, E, iters) if with_folds else {})
    return dict(**extra, case=case, title=title or case, n=n, names=names, freq=s.mean(1).tolist(), joint=np.round(joint, 4).tolist(), corr=np.round(np.nan_to_num(np.corrcoef(Z.T)), 3).tolist(), models=out, n_distinct=int(len(counts)),
                iterations=int(s.shape[1]), shown_share=float(counts[order[:TOP]].sum() / s.shape[1]), y=dict(zip(lig, yv.tolist())), X={c: [float("%.6g" % v) for v in X[c]] for c in names}, ligands=lig, unit=unit, response=response,
                draw=draw, label=[label(c, w) for c, w in zip(names, draw)], desc=[describe(c, w) for c, w in zip(names, draw)])


def models(case):
    sys.path.insert(0, str(BASSA)); import bayes_checks as BC
    d = compute(*BC.data(case), case=case, **CASES[case], seeds=SEEDS, with_folds=True); json.dump(d, open(C / f"{case}_models.json", "w"))
    print(case, len(d["names"]), "descriptors,", d["n_distinct"], "distinct models; the", TOP, "most visited cover %.0f %% of the run" % (100 * d["shown_share"]))


def structures(case, ligands):
    """per compound: xyz, the view, and what the drawing needs (dipole and charges, or the atom roles of a metal complex; all 1-based)"""
    sys.path.insert(0, str(BASSA)); sys.path.insert(0, str(EXPLORER)); sys.path.insert(0, str(PKG)); import summary_axes as SA
    if case == "cs2":
        import data as TD
        for lig in ligands:
            p = TD.payload(str(next(FEATHERS.glob(f"FL_lig_{lig[1:]}_str_*.feather"))), lig); yield lig, dict(xyz=p["xyz"], dipole=p.get("dipole"), charges=p.get("charges"), view=SA.view(p["xyz"], 6, 7))
    elif case == "wu_doyle":
        t = pd.read_csv(WD / "tables/wu_doyle_sterimol_dipole.csv").set_index("ligand"); role = pd.read_csv(WD / "tables/wu_doyle_table2.csv").set_index("ligand").role
        for lig in ligands:
            r = t.loc[lig]; txt = next((WD / "geometries/numbered").glob(r["id"] + "_*.xyz")).read_text().replace(chr(13), ""); v = np.array(SA.view(txt, 1, 2)); c, sn = np.cos(np.radians(-40)), np.sin(np.radians(-40))
            yield lig, dict(xyz=txt, dipole=[float(r.dipole_metal_x), float(r.dipole_metal_y), float(r.dipole_metal_z), float(r.dipole_total)], charges=None,     # as summary_axes.wu_doyle: the file is written in the Ni->P frame
                            roles=dict(M=[1], D=[2], Pl=[int(r.plane_carbon)]), view=[list(v[0]), list(c * v[1] + sn * v[2]), list(-sn * v[1] + c * v[2])])
    else:
        klass = pd.read_csv(BASSA / "f50_cs3_cp.csv", index_col=0).ligand_class
        for lig in ligands: yield lig, dict(chelate(CUCL / f"{lig}_CuCl_xtbopt.xyz", SA), klass=str(klass[lig]))


def chelate(f, SA):
    """a metal complex with two donors: its structure, the atoms that play each role (1-based, one per arm) and a view along C*->R"""
    sys.path.insert(0, str(PKG)); from M2_data_extractor import metal_complex as mc
    c = mc.MetalComplex.from_xyz(str(f)); arms = [mc.stereocentre(c.symbols, c.adj, a, b, metal=c.metal, coords=c.coords) for a, b in ((c.donor_1, c.donor_2), (c.donor_2, c.donor_1))]
    if None in arms: raise RuntimeError(f"{Path(f).name}: no stereocentre on an arm")
    txt = Path(f).read_text(); roles = more_roles(txt, dict(M=[c.metal + 1] * 2, D=[c.donor_1 + 1, c.donor_2 + 1], Cs=[a[0] + 1 for a in arms], R=[a[1] + 1 for a in arms]))
    return dict(xyz=txt, dipole=None, charges=None, roles=roles, view=SA.view(txt, roles["Cs"][0], roles["R"][0]))


def torsion(X, a, b, c, d):
    b0, b1, b2 = X[a] - X[b], X[c] - X[b], X[d] - X[c]; b1 = b1 / np.linalg.norm(b1); v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


def canon_cs1(sym, X, mu):
    """The rule of notes/analysis/bassa/build_wide.py (canon_cs1), copied because importing that script reruns it and
    rewrites its tables. Keep the two in step. Atom 13 is the benzyl ortho carbon nearer the O8-C9 bond, 14 the other,
    the meta carbons and the hydrogens follow; a mirror-image conformer (aryl twist negative) is reflected.
    Returns the order (new position -> old atom), the coordinates and the dipole."""
    perm = list(range(len(sym)))
    if abs(torsion(X, 7, 8, 11, 12)) > abs(torsion(X, 7, 8, 11, 13)):
        hyd = lambda c: [i for i, e in enumerate(sym) if e == "H" and np.linalg.norm(X[i] - X[c]) < 1.2]
        for a, b in ((12, 13), (14, 16)):
            for i, j in [(a, b)] + list(zip(hyd(a), hyd(b))): perm[i], perm[j] = perm[j], perm[i]
    X = X[perm].copy(); mu = np.array(mu, float)
    orth = [j for j, e in enumerate(sym) if e != "H" and j not in (0, 22) and np.linalg.norm(X[j] - X[22]) < (1.95 if e == "S" else 1.65)]
    if min((torsion(X, 2, 0, 22, j) for j in orth), key=abs) < 0: X[:, 2] *= -1; mu[2] *= -1
    return perm, X, mu


def renumber(case, mols):
    """Structures in the numbering the descriptor table uses. CS1's table was built after `build_wide.canon_cs1`: the benzyl
    ring's ortho and meta carbons in one order for every substrate, and three mirror-image conformers reflected. The raw
    feathers are not in that order, so without this the page marks the wrong ring atom on five substrates."""
    if case != "cs1": return mols
    for n, m in mols.items():
        rows = m["xyz"].strip().splitlines(); sym = [l.split()[0] for l in rows[2:]]; X = np.array([[float(v) for v in l.split()[1:4]] for l in rows[2:]])
        perm, X2, mu = canon_cs1(sym, X, m["dipole"][:3]); assert [sym[i] for i in perm] == sym, n
        assert abs(torsion(X2, 7, 8, 11, 12)) <= abs(torsion(X2, 7, 8, 11, 13)), n                 # atom 13 is now the ortho carbon nearer the O8-C9 bond
        m.update(xyz="\n".join(rows[:2] + ["%s %.8f %.8f %.8f" % (e, *x) for e, x in zip(sym, X2)]) + "\n", dipole=[*map(float, mu), m["dipole"][3]], charges={t: [v[i] for i in perm] for t, v in m["charges"].items()})
    return mols


def lines(case):
    """the reference numbers of a case from the report's own runs: the 1,000-shuffle chance line when it exists (deviation 27), else the 100-shuffle one"""
    out = {}; rep = BASSA / "sets/fixed_core" / f"bayes_{case}.json"; big = BASSA / "sets/fixed_core" / f"null1000_{case}.json"
    if rep.exists(): b = json.load(open(rep)); out.update(chance=b["null"]["null_top_q95"], p_top=b["null"]["p_top"], nested_q2=b["calibration"]["q2"], n_shuffles=b["null"]["n_perm"])
    if big.exists() and "all" in json.load(open(big)): a = json.load(open(big))["all"]; out.update(chance=a["q95"], chance_lo=a["q95_lo"], chance_hi=a["q95_hi"], p_top=a["p"], p_lo=a["p_lo"], p_hi=a["p_hi"], n_shuffles=a["n"])
    return out


def external(case, names):
    """the compounds the case study held out: the same columns, their response. Never fitted here."""
    if case == "cs2":
        t = pd.read_csv(BASSA / "f50_cs2.csv", index_col=0); ext = ["L3", "L12", "L15", "L16", "L22"]; X, y = t.loc[ext, names], t.loc[ext, "y"]
    elif case == "wu_doyle":
        a = pd.read_csv(WD / "tables/wu_doyle_sterimol_dipole.csv").set_index("ligand"); b = pd.read_csv(WD / "tables/wu_doyle_table2.csv").set_index("ligand"); ext = list(a.index[b.loc[a.index, "role"] == "heldout"])
        X = pd.concat([a, b.loc[a.index, ["vmin", "vbur_percent", "cone_angle"]]], axis=1).loc[ext, names]; y = b.loc[ext, "yield_percent"].astype(float)
    else: return None
    return dict(names=ext, X={c: [float("%.6g" % v) for v in X[c]] for c in names}, y=y.tolist())


def case_file(c):
    """where a case is read from: the cases folder, else the package (the bundled example)"""
    return C / f"{c}.js" if (C / f"{c}.js").exists() or not (H / f"{c}.js").exists() else H / f"{c}.js"


def manifest(add=None):
    """manifest.js lists the case files the page loads; prepared cases first, then the ones kept from the page"""
    C.mkdir(parents=True, exist_ok=True)
    f = C / "manifest.js"; have = json.loads(f.read_text(encoding="utf-8").split("=", 1)[1].rstrip(";\n")) if f.exists() else ["example"]
    if add and add not in have: have.append(add)
    have = [c for c in have if (C / f"{c}.js").exists() or (H / f"{c}.js").exists()]; f.write_text("window.CASE_FILES = " + json.dumps(have) + ";\n", encoding="utf-8"); return have


def standalone(cases):
    """the page, the explorer and the given cases in one HTML file that needs no server (and so cannot run new tables)"""
    page_ = (H / "index.html").read_text(encoding="utf-8"); ex = (EXPLORER / "feature_explorer.html").read_text(encoding="utf-8")
    a, b = '<script src="manifest.js"></script>', 'src="feature_explorer.html?embed=figure"'; assert page_.count(a) == 1 and page_.count(b) == 2
    data = "<script>window.STANDALONE = true; window.CASES = " + json.dumps(cases).replace("</", "<\\/") + ";</script>"
    return page_.replace(a, data).replace(b, 'srcdoc="' + ex.replace("&", "&amp;").replace('"', "&quot;") + '"')


def page(case):
    sys.path.insert(0, str(BASSA)); import summary_axes as SA
    d = json.load(open(C / f"{case}_models.json")); d["draw"] = [drawing(c) for c in d["names"]]; d["label"] = [label(c, w) for c, w in zip(d["names"], d["draw"])]; d.update(CASES[case])
    print("no drawing for:", [c for c, w in zip(d["names"], d["draw"]) if w is None]); d["ext"] = external(case, d["names"])
    d["mols"] = dict(structures(case, d["ligands"] + (d["ext"]["names"] if d["ext"] else [])))
    d.update(lines(case))                                                       # the two lines that stay on screen
    d["style"] = dict(SA.LOOK, values=True, atomLabels="all", shadow=False, **LOOKS.get(case, {}))
    (C / f"{case}.js").write_text("(window.CASES = window.CASES || {})[%s] = %s;\n" % (json.dumps(case), json.dumps(d)), encoding="utf-8")
    print("wrote", case + ".js", len(d["mols"]), "structures; cases:", manifest(case))


if __name__ == "__main__":
    for c, w in (("dipole_x_1-2-19", "dipole"), ("total_dipole_1-2-19", "dipole_total"), ("loc_B1_1-4", "sterimol"), ("B1_B5_azimuth_1-4", "sterimol"), ("Frequency_Bending_5-6", "mark"), ("Frequency_Stretch_1_2", "mark"), ("cross_9_0", "mark"), ("para_angle_9_0", "mark"), ("aniso", None),
                 ("B1_6-4", "sterimol"), ("ster_loc_B5[C6->C7]", "sterimol"), ("q_cm5[C6]", "charge"), ("nbo_atom_15", "charge"), ("mu_y[C4->N5|C6]", "dipole"), ("0_angle_[5, 6, 7]", "angle"), ("pi_occ[C1=C3]", "mark"), ("stretch[C1=C3]", "mark"), ("stretch[C3-O8]", "mark"), ("ring_hi[C12]", "mark"), ("E2_O_to_pistar[O8|C1=C3]", "mark"), ("pi_share[C3|C1=C3]", "mark"), ("polarizability[C1->C3]", "mark"), ("atomic_dipole[C1|C1->C3]", "mark"), ("aniso", None),
                 ("CsR_B1_sym", "sterimol"), ("fromM_L_asym", "sterimol"), ("MD_phi_asym", "sterimol"), ("MD_asym", "atoms"), ("q_metal", "atoms"), ("tf_chi1", None), ("L", "sterimol"), ("dipole_metal_y", "dipole"), ("cone_angle", None)):
        assert (drawing(c) or {}).get("kind") == w, c
    assert drawing("CsR_B1_sym")["ra"] == "Cs" and label("sub_theta_asym", drawing("sub_theta_asym")) == "θ on C*→R · difference between the arms" and label("RRp_loc_B5_sym", drawing("RRp_loc_B5_sym")) == "loc B5 on R→R′ · mean of the two arms" and "tf_" not in label("tf_chi3v", None)
    sys.path.insert(0, str(BASSA))
    if sys.argv[1] == "relabel":                                                # build_case.py relabel: names, hover notes, drawings, R′ and the reaction scheme again for every kept case, from its column names; no run is repeated
        for case in json.loads((C / "manifest.js").read_text(encoding="utf-8").split("=", 1)[1].rstrip(";\n ")):
            f = case_file(case); t = f.read_text(encoding="utf-8")
            if f.parent != C: continue                                          # the bundled example is not rewritten
            d = json.loads(t[t.index("] = ") + 4:-2]); d["draw"] = [drawing(c) for c in d["names"]]
            d["label"] = [label(c, w) for c, w in zip(d["names"], d["draw"])]; d["desc"] = [describe(c, w) for c, w in zip(d["names"], d["draw"])]; sc = scheme(case)
            if sc: d["scheme"] = sc
            for m in (d.get("mols") or {}).values():
                if m.get("roles") and "Cs" in m["roles"]: m["roles"] = more_roles(m["xyz"], {k: m["roles"][k] for k in ("M", "D", "Cs", "R")})
            f.write_text("(window.CASES = window.CASES || {})[%s] = %s;\n" % (json.dumps(case), json.dumps(d)), encoding="utf-8"); print(case, "relabelled;", sum(w is None for w in d["draw"]), "of", len(d["names"]), "columns not drawn;", "scheme" if sc else "no scheme")
        sys.exit()
    if sys.argv[1] == "attach":                                                 # build_case.py attach <kept case> <folder> ["title"] ["response label"] ["unit"]: structures (and names) for a case kept without them
        import read_structures as RS
        f = C / f"{sys.argv[2]}.js"; t = f.read_text(encoding="utf-8"); d = json.loads(t[t.index("] = ") + 4:-2]); got = RS.main(sys.argv[3], d["ligands"] + (d["ext"]["names"] if d.get("ext") else []))
        got_lines = lines(sys.argv[2]); d.update(got_lines); print("from the report's runs:", sorted(got_lines))
        try:                                                                    # a folder of chelate complexes: atom roles and a view, so the metal-referenced names are drawn on both arms
            import summary_axes as SA
            files = [x for x in Path(sys.argv[3]).iterdir() if x.suffix.lower() == ".xyz"]; got["found"] = {n: chelate(RS.match(n, files), SA) for n in got["found"]}; print("read as metal complexes with two donors")
        except Exception as e: print("not read as metal complexes (%s); plain structures kept" % str(e)[:80])
        tab = BASSA / f"f50_{sys.argv[2]}.csv"
        if tab.exists() and "ligand_class" in pd.read_csv(tab, nrows=1).columns:
            kl = pd.read_csv(tab, index_col=0).ligand_class
            for n in got["found"]: got["found"][n]["klass"] = str(kl[n])
        if "--raw" not in sys.argv: got["found"] = renumber(sys.argv[2], got["found"])
        for c, w in zip(d["names"], [drawing(c) for c in d["names"]]):         # where a column is a charge on a numbered atom, the structure shown must carry that value on that atom
            if w and w["kind"] == "charge" and all(m.get("charges") and w["type"] in m["charges"] for m in got["found"].values()):
                off = [n for k, n in enumerate(d["ligands"]) if n in got["found"] and abs(got["found"][n]["charges"][w["type"]][w["i"] - 1] - d["X"][c][k]) > 2e-4]
                if off: print("WARNING %s: the table and the structure disagree for %s" % (c, off))
        d["mols"] = got["found"]; d.update({k: v for k, v in zip(("title", "response", "unit"), [a for a in sys.argv[4:] if not a.startswith("--")])}); d["draw"] = [drawing(c) for c in d["names"]]; d["label"] = [label(c, w) for c, w in zip(d["names"], d["draw"])]
        f.write_text("(window.CASES = window.CASES || {})[%s] = %s;\n" % (json.dumps(sys.argv[2]), json.dumps(d)), encoding="utf-8"); print(len(got["found"]), "structures attached; none for", got["missing"])
    elif sys.argv[1] == "export":                                                 # build_case.py export cs2 cs3_cp -> explorer_standalone.html
        got = {}
        for c in sys.argv[2:] or manifest():
            t = case_file(c).read_text(encoding="utf-8"); got[c] = json.loads(t[t.index("] = ") + 4:-2])
        (C / "explorer_standalone.html").write_text(standalone(got), encoding="utf-8"); print("wrote explorer_standalone.html with", list(got))
    else: dict(models=models, page=page)[sys.argv[1]](sys.argv[2] if len(sys.argv) > 2 else "cs2")
