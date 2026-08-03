"""
Generate the animated SVG figures embedded in README.md.

Every figure is drawn from geometry this script actually computes, so the
picture and the formula agree with the implementation:

  frame.svg      dipole_utils._build_basis      (Gram-Schmidt frame)
  sterimol.svg   sterimol_utils.scan_b1_over_angles / get_extended_df_for_sterimol
  vibration.svg  vibrations_utils.calc_vibration_dot_product
  crossval.svg   modeling.LinearRegressionModel  (out-of-fold Q2)

Animation is SMIL, which renders inside <img> on GitHub. JavaScript never
runs in a Markdown README, so none is used here.

    python docs/animations/build_animations.py
"""

import numpy as np
from pathlib import Path

OUT = Path(__file__).resolve().parent

# Mid-tone palette: legible on both the light and the dark GitHub theme.
INK, FAINT, HAIR = "#7d8797", "#a7b0bd", "#c3cad4"
TEAL, BLUE, ORANGE, VIOLET, GREEN = "#12a5a5", "#4b8df8", "#e87d35", "#9b7bf0", "#2fa46a"
ATOM = "#9aa4b4"
MONO = "ui-monospace,SFMono-Regular,Menlo,Consolas,monospace"
SANS = "ui-sans-serif,system-ui,-apple-system,Segoe UI,Roboto,sans-serif"

CYCLE = 18.0  # seconds, shared by every figure


def svg(w, h, body, title, desc):
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" '
        f'height="{h}" font-family="{SANS}" role="img" aria-label="{desc}">\n'
        f"<title>{title}</title>\n<desc>{desc}</desc>\n"
        f'<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" '
        f'markerHeight="6" orient="auto-start-reverse">'
        f'<path d="M0 0 L10 5 L0 10 z" fill="context-stroke"/></marker></defs>\n'
        f"{body}\n</svg>\n"
    )


def fade(t0, t1, hold=0.35):
    """SMIL opacity keyframes: invisible, fade in at t0, fade out after t1."""
    k = [0.0, t0 - hold, t0, t1, t1 + hold, CYCLE]
    k = [max(0.0, min(CYCLE, x)) / CYCLE for x in k]
    return (
        f'<animate attributeName="opacity" values="0;0;1;1;0;0" '
        f'keyTimes="{";".join(f"{x:.4f}" for x in k)}" dur="{CYCLE}s" repeatCount="indefinite"/>'
    )


def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def txt(x, y, s, size=13, fill=INK, anchor="start", mono=False, weight="normal"):
    fam = f' font-family="{MONO}"' if mono else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}"{fam} fill="{fill}" '
        f'text-anchor="{anchor}" font-weight="{weight}">{esc(s)}</text>'
    )


def arrow(x1, y1, x2, y2, color, width=2.0, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" '
        f'stroke-width="{width}" marker-end="url(#ah)" stroke-linecap="round"{d}/>'
    )


def trace(d, color, width, dur, begin=0.0, length=1000):
    """A path that draws itself once per cycle."""
    return (
        f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}" '
        f'stroke-linecap="round" pathLength="{length}" stroke-dasharray="{length}" '
        f'stroke-dashoffset="{length}">'
        f'<animate attributeName="stroke-dashoffset" values="{length};0" dur="{dur}s" '
        f'begin="{begin}s" repeatCount="indefinite"/></path>'
    )


# --------------------------------------------------------------------------- #
# 1. Common reference frame  (dipole_utils._build_basis)
# --------------------------------------------------------------------------- #
def build_frame():
    W, H = 780, 400
    cx, cy, S = 235, 210, 46  # screen centre and px per Angstrom

    def P(v):  # model -> screen (y flips)
        return cx + v[0] * S, cy - v[1] * S

    # A benzene ring with a substituent, in the molecule's own arbitrary pose.
    ring = np.array([[np.cos(a), np.sin(a)] for a in np.deg2rad(np.arange(0, 360, 60) + 18)]) * 1.39
    sub = ring[0] * 2.35
    pose = np.deg2rad(-34)
    R = np.array([[np.cos(pose), -np.sin(pose)], [np.sin(pose), np.cos(pose)]])
    off = np.array([0.55, 0.30])
    ring, sub = ring @ R.T + off, sub @ R + off

    origin_set = [2, 3, 4]
    o = ring[origin_set].mean(axis=0)         # o = centroid(origin set)
    y_hat = ring[0] - o
    y_hat /= np.linalg.norm(y_hat)            # y = normalize(r_y - o)
    c = sub - o
    c /= np.linalg.norm(c)                    # c = normalize(r_plane - o)
    proj = float(c @ y_hat) * y_hat           # (c . y) y
    x_raw = c - proj                          # Gram-Schmidt residual
    x_hat = x_raw / np.linalg.norm(x_raw)

    b = []
    b.append(f'<rect width="{W}" height="{H}" fill="none"/>')
    b.append(txt(28, 40, "A common reference frame", 17, INK, weight="600"))
    b.append(txt(28, 62, "Descriptors only compare across a series if every molecule sits in the same frame.", 12.5, FAINT))

    # --- molecule (static) ---
    mol = ['<g>']
    for i in range(6):
        x1, y1 = P(ring[i]); x2, y2 = P(ring[(i + 1) % 6])
        mol.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{ATOM}" stroke-width="2.6"/>')
    x1, y1 = P(ring[0]); x2, y2 = P(sub)
    mol.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{ATOM}" stroke-width="2.6"/>')
    for i, p in enumerate(ring):
        px, py = P(p)
        col = TEAL if i in origin_set else (VIOLET if i == 0 else ATOM)
        mol.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="6.5" fill="{col}"/>')
    px, py = P(sub)
    mol.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="7.5" fill="{ORANGE}"/>')
    mol.append("</g>")
    b += mol

    ox, oy = P(o)
    yx, yy = P(o + y_hat * 1.9)
    cxs, cys = P(o + c * 2.1)
    prx, pry = P(o + proj * 2.1)
    xx, xy = P(o + x_hat * 1.5)
    zx, zy = P(o)

    # --- stage 1: the three selections ---
    b.append(f'<g opacity="0">{fade(0.6, 4.2)}')
    b.append(txt(28, 336, "1.  pick three groups of atoms", 13, TEAL, weight="600"))
    b.append(txt(28, 356, "origin set (teal) · y-atom (violet) · plane atom (orange)", 12, FAINT))
    for i in origin_set:
        px, py = P(ring[i])
        b.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="13" fill="none" stroke="{TEAL}" stroke-width="1.6" opacity=".8"/>')
    px, py = P(ring[0]); b.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="13" fill="none" stroke="{VIOLET}" stroke-width="1.6" opacity=".8"/>')
    px, py = P(sub); b.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="14" fill="none" stroke="{ORANGE}" stroke-width="1.6" opacity=".8"/>')
    b.append("</g>")

    # --- stage 2: origin ---
    b.append(f'<g opacity="0">{fade(4.4, 8.0)}')
    b.append(f'<circle cx="{ox:.1f}" cy="{oy:.1f}" r="5" fill="none" stroke="{INK}" stroke-width="2"/>')
    b.append(f'<circle cx="{ox:.1f}" cy="{oy:.1f}" r="11" fill="none" stroke="{INK}" stroke-width="1" opacity=".5"/>')
    b.append(txt(ox + 15, oy + 5, "o", 13, INK, mono=True))
    b.append(txt(28, 336, "2.  origin = centroid of the origin set", 13, TEAL, weight="600"))
    b.append(txt(28, 358, "o = (r₂ + r₃ + r₄) / 3        rᵢ → rᵢ − o", 12.5, FAINT, mono=True))
    b.append("</g>")

    # --- stage 3: y axis ---
    b.append(f'<g opacity="0">{fade(8.2, 11.4)}')
    b.append(arrow(ox, oy, yx, yy, VIOLET, 2.4))
    b.append(txt(yx + 8, yy - 4, "ŷ", 14, VIOLET, mono=True, weight="600"))
    b.append(txt(28, 336, "3.  the y axis points at the y-atom", 13, VIOLET, weight="600"))
    b.append(txt(28, 358, "ŷ = (r_y − o) / ‖r_y − o‖", 12.5, FAINT, mono=True))
    b.append("</g>")

    # --- stage 4: Gram-Schmidt ---
    b.append(f'<g opacity="0">{fade(11.6, 15.4)}')
    b.append(arrow(ox, oy, yx, yy, VIOLET, 2.0))
    b.append(arrow(ox, oy, cxs, cys, ORANGE, 2.2))
    b.append(txt(cxs + 8, cys + 2, "c", 13, ORANGE, mono=True, weight="600"))
    b.append(f'<line x1="{ox:.1f}" y1="{oy:.1f}" x2="{prx:.1f}" y2="{pry:.1f}" stroke="{VIOLET}" stroke-width="6" opacity=".22" stroke-linecap="round"/>')
    b.append(f'<line x1="{prx:.1f}" y1="{pry:.1f}" x2="{cxs:.1f}" y2="{cys:.1f}" stroke="{HAIR}" stroke-width="1.4" stroke-dasharray="4 4"/>')
    b.append(arrow(ox, oy, xx, xy, TEAL, 2.6))
    b.append(txt(xx - 4, xy + 20, "x̂", 14, TEAL, mono=True, weight="600"))
    b.append(txt(28, 336, "4.  strip the y-component off c — what is left is x", 13, TEAL, weight="600"))
    b.append(txt(28, 358, "x̂ = normalize( c − (c·ŷ) ŷ )        ẑ = x̂ × ŷ", 12.5, FAINT, mono=True))
    b.append("</g>")

    # --- stage 5: frame + transform ---
    b.append(f'<g opacity="0">{fade(15.6, 17.6)}')
    b.append(arrow(ox, oy, yx, yy, VIOLET, 2.4))
    b.append(arrow(ox, oy, xx, xy, TEAL, 2.4))
    b.append(f'<circle cx="{zx:.1f}" cy="{zy:.1f}" r="5.5" fill="none" stroke="{GREEN}" stroke-width="2.2"/>')
    b.append(f'<circle cx="{zx:.1f}" cy="{zy:.1f}" r="1.8" fill="{GREEN}"/>')
    b.append(txt(zx - 20, zy - 12, "ẑ", 13, GREEN, mono=True, weight="600"))
    b.append(txt(28, 336, "5.  every atom re-expressed in that frame", 13, GREEN, weight="600"))
    b.append(txt(28, 358, "r′ = B (r − o),   B = [ x̂ ; ŷ ; ẑ ]", 12.5, FAINT, mono=True))
    b.append("</g>")

    # --- right-hand panel: the whole recipe, always visible ---
    px0 = 470
    b.append(f'<line x1="{px0-28}" y1="96" x2="{px0-28}" y2="330" stroke="{HAIR}" stroke-width="1"/>')
    rows = [
        ("o", "centroid of the origin set", INK),
        ("ŷ", "(r_y − o) / ‖r_y − o‖", VIOLET),
        ("c", "(r_plane − o) / ‖r_plane − o‖", ORANGE),
        ("x̂", "normalize( c − (c·ŷ) ŷ )", TEAL),
        ("ẑ", "x̂ × ŷ", GREEN),
        ("r′", "B (r − o)", INK),
    ]
    b.append(txt(px0, 118, "the frame, in full", 12.5, FAINT, weight="600"))
    for i, (sym, rhs, col) in enumerate(rows):
        y = 150 + i * 30
        b.append(txt(px0, y, sym, 15, col, mono=True, weight="600"))
        b.append(txt(px0 + 30, y, "=", 13, HAIR, mono=True))
        b.append(txt(px0 + 48, y, rhs, 12.5, FAINT, mono=True))
        b.append(f'<g opacity="0">{fade(0.6 + i * 2.55, 2.9 + i * 2.55, 0.3)}'
                 f'<rect x="{px0-12}" y="{y-16}" width="292" height="24" rx="5" fill="{col}" opacity=".10"/></g>')
    return svg(W, H, "\n".join(b), "Building a common reference frame",
               "Animation: the origin is the centroid of the chosen atoms, y points at the y-atom, "
               "and x is what remains of the plane atom after its y-component is removed.")


# --------------------------------------------------------------------------- #
# 2. Sterimol B1 / B5  (sterimol_utils)
# --------------------------------------------------------------------------- #
def build_sterimol():
    W, H = 780, 466
    cx, cy, S = 210, 232, 30

    # Cross-section perpendicular to L: heavy atoms with their vdW radii (Angstrom).
    atoms = np.array([[0.00, 0.00, 1.70], [1.62, 0.58, 1.70],
                      [-1.28, 1.04, 1.70], [0.22, -1.70, 1.52]])
    pts, rad = atoms[:, :2], atoms[:, 2]

    th = np.linspace(0, 2 * np.pi, 361)
    u = np.stack([np.cos(th), np.sin(th)], axis=1)
    # Half-width of the radius-expanded cloud in direction u (the supporting line).
    halfwidth = (pts @ u.T + rad[:, None]).max(axis=0)
    B1, B5 = halfwidth.min(), float((np.linalg.norm(pts, axis=1) + rad).max())
    i1 = int(halfwidth.argmin())

    b = [txt(28, 40, "Sterimol B1 and B5", 17, INK, weight="600"),
         txt(28, 62, "Look down the L axis. Sweep a supporting line around the substituent.", 12.5, FAINT)]

    # axes
    b.append(f'<line x1="{cx-150}" y1="{cy}" x2="{cx+150}" y2="{cy}" stroke="{HAIR}" stroke-width="1"/>')
    b.append(f'<line x1="{cx}" y1="{cy-150}" x2="{cx}" y2="{cy+150}" stroke="{HAIR}" stroke-width="1"/>')

    # vdW disks
    for (x, y), r in zip(pts, rad):
        b.append(f'<circle cx="{cx+x*S:.1f}" cy="{cy-y*S:.1f}" r="{r*S:.1f}" fill="{ATOM}" opacity=".20"/>')
        b.append(f'<circle cx="{cx+x*S:.1f}" cy="{cy-y*S:.1f}" r="{r*S:.1f}" fill="none" stroke="{ATOM}" stroke-width="1" opacity=".55"/>')
    for (x, y) in pts:
        b.append(f'<circle cx="{cx+x*S:.1f}" cy="{cy-y*S:.1f}" r="4.5" fill="{ATOM}"/>')
    b.append(f'<circle cx="{cx}" cy="{cy}" r="3" fill="{INK}"/>')
    b.append(txt(cx - 16, cy + 18, "L", 12, INK, mono=True))

    # B5: the single furthest reach, independent of the sweep
    j5 = int((np.linalg.norm(pts, axis=1) + rad).argmax())
    d5 = pts[j5] / np.linalg.norm(pts[j5])
    b.append(arrow(cx, cy, cx + d5[0] * B5 * S, cy - d5[1] * B5 * S, ORANGE, 2.2))
    b.append(txt(cx + d5[0] * B5 * S - 26, cy - d5[1] * B5 * S - 10, f"B5 {B5:.2f} Å", 12, ORANGE, mono=True, weight="600"))

    # rotating supporting line + radius, length driven by the computed half-width
    step = halfwidth[::5]
    vals = ";".join(f"{v*S:.2f}" for v in step) + f";{halfwidth[0]*S:.2f}"
    b.append(f'<g transform="translate({cx},{cy})">')
    b.append('<g><animateTransform attributeName="transform" type="rotate" from="0" to="-360" '
             f'dur="{CYCLE}s" repeatCount="indefinite"/>')
    b.append(f'<line x1="0" y1="0" x2="60" y2="0" stroke="{BLUE}" stroke-width="2" marker-end="url(#ah)">'
             f'<animate attributeName="x2" values="{vals}" dur="{CYCLE}s" repeatCount="indefinite"/></line>')
    b.append(f'<line x1="60" y1="-92" x2="60" y2="92" stroke="{BLUE}" stroke-width="2.2" opacity=".85">'
             f'<animate attributeName="x1" values="{vals}" dur="{CYCLE}s" repeatCount="indefinite"/>'
             f'<animate attributeName="x2" values="{vals}" dur="{CYCLE}s" repeatCount="indefinite"/></line>')
    b.append("</g></g>")

    # B1 marker, revealed as the sweep passes the minimum
    t1 = CYCLE * i1 / 360.0
    d1 = u[i1]
    b.append(f'<g opacity="0">{fade(t1 + 0.15, CYCLE - 0.4, 0.25)}')
    b.append(arrow(cx, cy, cx + d1[0] * B1 * S, cy - d1[1] * B1 * S, TEAL, 2.6))
    b.append(txt(cx + d1[0] * B1 * S + 10, cy - d1[1] * B1 * S + 16, f"B1 {B1:.2f} Å", 12, TEAL, mono=True, weight="600"))
    b.append("</g>")

    # --- right: half-width against sweep angle ---
    gx0, gx1, gy0, gy1 = 430, 748, 118, 330
    lo, hi = halfwidth.min() - 0.35, halfwidth.max() + 0.25
    sx = lambda t: gx0 + (t / 360.0) * (gx1 - gx0)
    sy = lambda v: gy1 - (v - lo) / (hi - lo) * (gy1 - gy0)

    b.append(txt(gx0, 96, "half-width  b(θ)  as the line sweeps", 12.5, FAINT, weight="600"))
    b.append(f'<line x1="{gx0}" y1="{gy1}" x2="{gx1}" y2="{gy1}" stroke="{HAIR}" stroke-width="1"/>')
    b.append(f'<line x1="{gx0}" y1="{gy0}" x2="{gx0}" y2="{gy1}" stroke="{HAIR}" stroke-width="1"/>')
    for t in (0, 90, 180, 270, 360):
        b.append(f'<line x1="{sx(t):.1f}" y1="{gy1}" x2="{sx(t):.1f}" y2="{gy1+4}" stroke="{HAIR}" stroke-width="1"/>')
        b.append(txt(sx(t), gy1 + 18, f"{t}°", 10.5, FAINT, anchor="middle", mono=True))
    b.append(f'<line x1="{gx0}" y1="{sy(B1):.1f}" x2="{gx1}" y2="{sy(B1):.1f}" stroke="{TEAL}" stroke-width="1.2" stroke-dasharray="5 4" opacity=".8"/>')
    b.append(txt(gx1, sy(B1) - 8, f"B1 = min b(θ) = {B1:.2f} Å", 11.5, TEAL, anchor="end", mono=True))

    d = "M " + " L ".join(f"{sx(t):.1f} {sy(v):.1f}" for t, v in zip(np.degrees(th), halfwidth))
    b.append(trace(d, BLUE, 2.2, CYCLE))
    b.append(f'<circle r="4.5" fill="{BLUE}"><animateMotion dur="{CYCLE}s" repeatCount="indefinite" '
             f'path="{d}"/></circle>')
    b.append(f'<g opacity="0">{fade(t1 + 0.15, CYCLE - 0.4, 0.25)}'
             f'<circle cx="{sx(i1):.1f}" cy="{sy(B1):.1f}" r="5.5" fill="none" stroke="{TEAL}" stroke-width="2"/></g>')

    b.append(f'<line x1="28" y1="408" x2="752" y2="408" stroke="{HAIR}" stroke-width="1"/>')
    b.append(txt(28, 430, f"L  = maxᵢ (yᵢ + rᵢ)                    B5 = maxᵢ (‖(xᵢ,zᵢ)‖ + rᵢ) = {B5:.2f} Å", 12, FAINT, mono=True))
    b.append(txt(28, 450, f"B1 = min over θ of the smallest half-width = {B1:.2f} Å", 12, FAINT, mono=True))
    return svg(W, H, "\n".join(b), "How Sterimol B1 and B5 are found",
               "Animation: a supporting line sweeps around the substituent cross-section; the "
               "minimum of its distance from the axis is B1, and the furthest atom reach is B5.")


# --------------------------------------------------------------------------- #
# 3. Picking a vibrational mode  (vibrations_utils)
# --------------------------------------------------------------------------- #
def build_vibration():
    W, H = 780, 424
    rng = np.random.default_rng(11)
    freqs = np.array([412, 655, 830, 1010, 1188, 1345, 1502, 1638, 1745, 1980, 2240, 3055])
    align = np.array([.08, .17, .11, .32, .21, .44, .29, .93, .58, .19, .12, .71])
    lo, hi = 1400, 3500                       # the frequency window
    inwin = (freqs > lo) & (freqs < hi)
    pick = int(np.where(inwin, align, -1).argmax())

    b = [txt(28, 40, "Choosing the stretch that belongs to your bond", 17, INK, weight="600"),
         txt(28, 62, "Each normal mode is scored by how much of its motion lies along the bond.", 12.5, FAINT)]

    # --- left: the bond and the mode's displacement vectors ---
    ax, ay, bx, by = 92, 250, 268, 196
    ux, uy = bx - ax, by - ay
    n = np.hypot(ux, uy); ux, uy = ux / n, uy / n
    b.append(f'<g><animateTransform attributeName="transform" type="translate" '
             f'values="0 0;{ux*7:.2f} {uy*7:.2f};0 0;{-ux*7:.2f} {-uy*7:.2f};0 0" '
             f'dur="0.9s" repeatCount="indefinite"/>')
    b.append(f'<line x1="{ax}" y1="{ay}" x2="{bx}" y2="{by}" stroke="{ATOM}" stroke-width="3"/>')
    b.append(f'<circle cx="{ax}" cy="{ay}" r="13" fill="{ATOM}"/><circle cx="{bx}" cy="{by}" r="10" fill="{ATOM}"/>')
    b.append("</g>")
    b.append(arrow(ax, ay, ax - ux * 52, ay - uy * 52, ORANGE, 2.2))
    b.append(arrow(bx, by, bx + ux * 52, by + uy * 52, ORANGE, 2.2))
    b.append(txt(ax - ux * 52 - 26, ay - uy * 52 - 8, "dᵃ", 12.5, ORANGE, mono=True, weight="600"))
    b.append(txt(bx + ux * 52 + 8, by + uy * 52 + 4, "dᵇ", 12.5, ORANGE, mono=True, weight="600"))
    b.append(arrow(ax, ay, ax + ux * 92, ay + uy * 92, TEAL, 1.8, dash="5 4"))
    b.append(txt(ax + ux * 100, ay + uy * 100 + 16, "û", 13, TEAL, mono=True, weight="600"))
    b.append(txt(28, 322, "aₖ = |dᵃₖ · û| + |dᵇₖ · û|", 13, FAINT, mono=True))
    b.append(txt(28, 344, "pick argmaxₖ aₖ  subject to  1400 < νₖ < 3500 cm⁻¹", 12, FAINT, mono=True))

    # --- right: score per mode ---
    gx0, gx1, gy0, gy1 = 360, 748, 110, 286
    sx = lambda f: gx0 + (f - 300) / (3200 - 300) * (gx1 - gx0)
    sy = lambda v: gy1 - v * (gy1 - gy0)

    win_x1 = min(sx(hi), gx1)  # the window runs past the plotted range; clip it
    b.append(f'<rect x="{sx(lo):.1f}" y="{gy0}" width="{win_x1-sx(lo):.1f}" height="{gy1-gy0}" fill="{TEAL}" opacity=".07"/>')
    b.append(f'<line x1="{sx(lo):.1f}" y1="{gy0}" x2="{sx(lo):.1f}" y2="{gy1}" stroke="{TEAL}" stroke-width="1.3" stroke-dasharray="5 4"/>')
    b.append(txt(sx(lo) + 6, gy0 - 8, "frequency window", 11, TEAL, mono=True))
    b.append(f'<line x1="{gx0}" y1="{gy1}" x2="{gx1}" y2="{gy1}" stroke="{HAIR}" stroke-width="1"/>')
    b.append(txt(gx0, gy1 + 34, "ν  (cm⁻¹)", 11, FAINT, mono=True))
    b.append(txt(gx0 - 8, gy0 + 4, "aₖ", 11.5, FAINT, anchor="end", mono=True))

    per = CYCLE / len(freqs)
    for i, (f, a) in enumerate(zip(freqs, align)):
        x, y = sx(f), sy(a)
        col = TEAL if i == pick else ATOM
        b.append(f'<rect x="{x-6:.1f}" y="{y:.1f}" width="12" height="{gy1-y:.1f}" rx="2.5" fill="{col}" opacity=".45"/>')
        b.append(f'<g opacity="0">{fade(i*per, i*per + per*0.85, 0.12)}'
                 f'<rect x="{x-6:.1f}" y="{y:.1f}" width="12" height="{gy1-y:.1f}" rx="2.5" fill="{BLUE}"/>'
                 f'<text x="{x:.1f}" y="{y-9:.1f}" font-size="10.5" font-family="{MONO}" fill="{BLUE}" '
                 f'text-anchor="middle">{a:.2f}</text></g>')
        if i in (0, len(freqs) - 1):
            b.append(txt(x, gy1 + 16, str(f), 10, FAINT, anchor="middle", mono=True))

    xw, yw = sx(freqs[pick]), sy(align[pick])
    b.append(f'<g opacity="0">{fade(pick*per + per, CYCLE - 0.3, 0.25)}')
    b.append(f'<circle cx="{xw:.1f}" cy="{yw:.1f}" r="10" fill="none" stroke="{TEAL}" stroke-width="2"/>')
    b.append(txt(xw, yw - 22, f"{freqs[pick]} cm⁻¹", 11.5, TEAL, anchor="middle", mono=True, weight="600"))
    b.append(txt(gx0, 322, f"selected: ν = {freqs[pick]} cm⁻¹,  a = {align[pick]:.2f}", 12.5, TEAL, mono=True, weight="600"))
    b.append("</g>")
    b.append(f'<line x1="28" y1="366" x2="752" y2="366" stroke="{HAIR}" stroke-width="1"/>')
    b.append(txt(28, 390, f"The {freqs[-1]} cm⁻¹ mode scores {align[-1]:.2f} — higher than most — but it is a C–H", 11.5, FAINT))
    b.append(txt(28, 408, "stretch far outside the window, so the threshold is what keeps it out.", 11.5, FAINT))
    return svg(W, H, "\n".join(b), "Selecting the stretching mode for a bond",
               "Animation: each normal mode is scored by the projection of its displacement "
               "vectors onto the bond axis; the highest scorer inside the frequency window is kept.")


# --------------------------------------------------------------------------- #
# 4. Out-of-fold cross-validation  (M3_modeler)
# --------------------------------------------------------------------------- #
def build_crossval():
    W, H = 780, 470
    rng = np.random.default_rng(4)
    n, k = 20, 5
    y = rng.normal(0, 1, n)
    yhat = y * 0.83 + rng.normal(0, 0.42, n)
    q2 = 1 - ((y - yhat) ** 2).sum() / ((y - y.mean()) ** 2).sum()

    cw, gap, x0, y0, rh = 26, 4, 60, 130, 30
    b = [txt(28, 40, "Why the score is honest", 17, INK, weight="600"),
         txt(28, 62, "Each sample is predicted only by a model that never saw it.", 12.5, FAINT)]
    b.append(txt(x0, y0 - 16, "5 folds · 20 samples", 11.5, FAINT, mono=True))

    per = CYCLE / (k + 1.6)
    for f in range(k):
        yy = y0 + f * rh
        b.append(txt(x0 - 12, yy + 15, f"{f+1}", 11, FAINT, anchor="end", mono=True))
        for i in range(n):
            xx = x0 + i * (cw + gap)
            val = (i % k) == f
            col = ORANGE if val else ATOM
            op = ".85" if val else ".22"
            b.append(f'<rect x="{xx}" y="{yy}" width="{cw}" height="20" rx="3.5" fill="{col}" opacity="{op}"/>')
        b.append(f'<g opacity="0">{fade(f*per, CYCLE - 0.4, 0.2)}'
                 f'<rect x="{x0-4}" y="{yy-3}" width="{n*(cw+gap)+2}" height="26" rx="5" fill="none" '
                 f'stroke="{ORANGE}" stroke-width="1.4" opacity=".55"/></g>')

    # the pooled out-of-fold vector
    oy = y0 + k * rh + 26
    b.append(txt(x0 - 12, oy + 15, "ŷ", 12, TEAL, anchor="end", mono=True, weight="600"))
    b.append(txt(x0, oy - 10, "out-of-fold predictions", 11.5, TEAL, mono=True))
    for i in range(n):
        xx = x0 + i * (cw + gap)
        b.append(f'<rect x="{xx}" y="{oy}" width="{cw}" height="20" rx="3.5" fill="{HAIR}" opacity=".35"/>')
        b.append(f'<g opacity="0">{fade((i % k)*per + 0.3, CYCLE - 0.4, 0.15)}'
                 f'<rect x="{xx}" y="{oy}" width="{cw}" height="20" rx="3.5" fill="{TEAL}" opacity=".8"/></g>')

    b.append(f'<g opacity="0">{fade(k*per + 0.4, CYCLE - 0.3, 0.3)}')
    b.append(txt(x0, oy + 64, "Q² = 1 −  Σ(yᵢ − ŷᵢᵒᵒᶠ)² / Σ(yᵢ − ȳ)²", 15, INK, mono=True, weight="600"))
    b.append(txt(x0 + 400, oy + 64, f"=  {q2:.3f}", 15, TEAL, mono=True, weight="600"))
    b.append("</g>")
    b.append(txt(x0, oy + 104, "Scaling is fitted inside each fold, never on the full set — that is what", 11.5, FAINT))
    b.append(txt(x0, oy + 122, "keeps Q² from flattering the model. Samples passed to --leave-out sit", 11.5, FAINT))
    b.append(txt(x0, oy + 140, "outside this loop entirely.", 11.5, FAINT))
    return svg(W, H, "\n".join(b), "Out-of-fold cross-validation and Q-squared",
               "Animation: five folds each hold out a different fifth of the samples; the pooled "
               "held-out predictions are what Q-squared is computed from.")


if __name__ == "__main__":
    for name, fn in [("frame", build_frame), ("sterimol", build_sterimol),
                     ("vibration", build_vibration), ("crossval", build_crossval)]:
        p = OUT / f"{name}.svg"
        p.write_text(fn(), encoding="utf-8")
        print(f"{p.name}  {p.stat().st_size/1024:.1f} KB")
