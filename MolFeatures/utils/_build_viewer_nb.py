# -*- coding: utf-8 -*-
"""Builder for sterimol_dipole_viewer.ipynb — keeps the notebook JSON sane."""
import json, os

CELL1 = r'''# ── Cell 1: Imports + Sterimol engine (ported from steriplot) ─────────────────
import sys, os, re, base64, urllib.request, warnings, traceback as _tb
import numpy as np
import pandas as pd
import py3Dmol
import ipywidgets as widgets
from IPython.display import display, clear_output, publish_display_data
import matplotlib
matplotlib.use('module://matplotlib_inline.backend_inline')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ── path setup so project modules import cleanly ──────────────────────────────
_here         = os.path.abspath('')          # .../MolFeatures/utils
_mol_features = os.path.dirname(_here)        # .../MolFeatures
_m2           = os.path.join(_mol_features, 'M2_data_extractor')
for _p in [_mol_features, _m2]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

# Only the dipole extractor is pulled from the project — the Sterimol geometry
# engine below is the self-contained steriplot one (gives local/global/sub).
from feather_extractor import process_gaussian_dipole_text

warnings.filterwarnings('ignore')

# ── radii / colours ───────────────────────────────────────────────────────────
VDW = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'P': 1.80,
    'S': 1.80, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98, 'B': 1.92, 'Si': 2.10,
    'Co': 2.00, 'Ni': 2.00, 'Fe': 2.04, 'Cu': 1.96, 'Zn': 2.01, 'Pd': 2.10,
}
COVALENT_R = {
    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'P': 1.07,
    'S': 1.05, 'Cl': 1.02, 'Br': 1.20, 'I': 1.39, 'B': 0.84, 'Si': 1.11,
    'Co': 1.50, 'Ni': 1.24, 'Fe': 1.52, 'Cu': 1.32, 'Zn': 1.22, 'Pd': 1.39,
}
ELEM_COLOR = {
    'H': '#d0d0d0', 'C': '#333333', 'N': '#4455cc', 'O': '#cc3333',
    'F': '#22aaaa', 'Cl': '#22bb33', 'S': '#ccaa22', 'Br': '#994422',
    'P': '#ff8833', 'I': '#7733aa', 'B': '#ffb5b5', 'Si': '#f0c8a0',
}

# ── XYZ reading ───────────────────────────────────────────────────────────────
def read_xyz(path):
    atoms = []
    with open(path) as f:
        for line in f.readlines()[2:]:
            p = line.split()
            if len(p) >= 4:
                atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
    return atoms

# ── 3-atom molecular frame ────────────────────────────────────────────────────
# origin -> (0,0,0); axis atom -> +Z (the L axis); plane atom fixes the X/Y
# azimuth so x/y components are reproducible.  B1/B5/L magnitudes are unchanged.
def build_frame(coords, oi, ai, pi=None):
    coords = np.asarray(coords, float)
    origin = coords[oi].copy()
    P = coords - origin
    z = P[ai].astype(float)
    nz = np.linalg.norm(z)
    if nz < 1e-9:
        raise ValueError('Direction atom coincides with origin atom.')
    z = z / nz
    if pi is not None and pi not in (oi, ai):
        ref = P[pi].astype(float)
    else:
        ref = np.array([1.0, 0.0, 0.0]) if abs(z[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    x = ref - np.dot(ref, z) * z
    if np.linalg.norm(x) < 1e-9:
        ref = np.array([0.0, 1.0, 0.0]) if abs(z[1]) < 0.9 else np.array([1.0, 0.0, 0.0])
        x = ref - np.dot(ref, z) * z
    x = x / np.linalg.norm(x)
    y = np.cross(z, x)
    R = np.vstack([x, y, z])          # rows = new axes;  new_components = R @ old_vec
    return R, origin

def apply_frame(symbols, coords, R, origin):
    P = np.asarray(coords, float) - origin
    new = P @ R.T
    return [(s, float(c[0]), float(c[1]), float(c[2])) for s, c in zip(symbols, new)]

# ── substructure (BFS on geometric bond graph) ────────────────────────────────
def build_bond_graph(atoms, tol=1.3):
    n = len(atoms)
    g = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            si, xi, yi, zi = atoms[i]
            sj, xj, yj, zj = atoms[j]
            ri = COVALENT_R.get(si, 0.77)
            rj = COVALENT_R.get(sj, 0.77)
            d = np.sqrt((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)
            if d < tol * (ri + rj):
                g[i].append(j); g[j].append(i)
    return g

def find_substructure_atoms(atoms, oi, ai):
    g = build_bond_graph(atoms)
    visited = {oi}
    queue = [ai]
    while queue:
        node = queue.pop(0)
        if node in visited:
            continue
        visited.add(node)
        queue.extend(g[node])
    return visited

# ── per-slice radial (support-function) scan ──────────────────────────────────
def section_disks(atoms, z_slice, allowed=None):
    disks = []
    for idx, (sym, x, y, z) in enumerate(atoms):
        if allowed is not None and idx not in allowed:
            continue
        R = VDW.get(sym)
        if R is None:
            continue
        dz = z_slice - z
        if abs(dz) > R:
            continue
        disks.append({'atom_idx': idx, 'sym': sym, 'x': x, 'y': y,
                      'r': np.sqrt(max(R*R - dz*dz, 0.0))})
    return disks

def support_fn(disks, theta):
    ux, uy = np.cos(theta), np.sin(theta)
    best_val, best_disk = -np.inf, None
    for d in disks:
        val = d['x']*ux + d['y']*uy + d['r']
        if val > best_val:
            best_val, best_disk = val, d
    return best_val, best_disk

def polar_width_scan(atoms, z_slice, n_theta=360, allowed=None):
    disks = section_disks(atoms, z_slice, allowed)
    if not disks:
        return None
    thetas = np.linspace(0, 2*np.pi, n_theta, endpoint=False)
    radii = np.zeros(n_theta)
    sup = []
    for i, th in enumerate(thetas):
        radii[i], disk = support_fn(disks, th)
        sup.append(disk)
    k = int(np.argmax(radii))
    return {'z': z_slice, 'thetas': thetas, 'radii': radii, 'all_atoms': sup,
            'max_radius': float(radii[k]), 'argmax_radius': k}

def scan_sterimol(atoms, n_z=240, n_theta=360, frac=0.30, allowed=None):
    relevant = [(sym, z) for i, (sym, _, _, z) in enumerate(atoms)
                if allowed is None or i in allowed]
    L = max(z + VDW.get(sym, 1.7) for sym, z in relevant)
    z_vals = np.linspace(0.0, L, n_z)
    scans = [s for s in (polar_width_scan(atoms, z, n_theta, allowed) for z in z_vals)
             if s is not None]
    global_B5 = max(s['max_radius'] for s in scans)
    min_valid = frac * global_B5
    valid = [s for s in scans if s['max_radius'] >= min_valid]
    if not valid:
        raise ValueError('No valid B1 slices for this selection.')
    b5_scan = max(scans, key=lambda s: s['max_radius'])
    b1_scan = min(valid, key=lambda s: s['max_radius'])
    k5, k1 = b5_scan['argmax_radius'], b1_scan['argmax_radius']
    return {
        'B1': b1_scan['max_radius'], 'B5': b5_scan['max_radius'], 'L': float(L),
        'z_B1': float(b1_scan['z']), 'z_B5': float(b5_scan['z']),
        'theta_B1': float(b1_scan['thetas'][k1]), 'theta_B5': float(b5_scan['thetas'][k5]),
        'plus_atom_B1': b1_scan['all_atoms'][k1], 'plus_atom_B5': b5_scan['all_atoms'][k5],
        'scan_B1': b1_scan, 'scan_B5': b5_scan, 'all_scans': scans,
        'min_valid_radius': float(min_valid),
    }

def compute_global_envelope(all_scans, n_theta=360):
    thetas = np.linspace(0, 2*np.pi, n_theta, endpoint=False)
    g_radii = np.zeros(n_theta)
    for s in all_scans:
        g_radii = np.maximum(g_radii, s['radii'])
    return thetas, g_radii

def global_best_atom(scans, k_theta):
    best_val, best_atom = -np.inf, None
    for sc in scans:
        if sc['radii'][k_theta] > best_val:
            best_val, best_atom = sc['radii'][k_theta], sc['all_atoms'][k_theta]
    return best_atom

def make_global(res, g_th, g_rd):
    k5, k1 = int(np.argmax(g_rd)), int(np.argmin(g_rd))
    return {
        'B1': float(g_rd[k1]), 'B5': float(g_rd[k5]),
        'theta_B1': float(g_th[k1]), 'theta_B5': float(g_th[k5]),
        'plus_atom_B1': global_best_atom(res['all_scans'], k1),
        'plus_atom_B5': global_best_atom(res['all_scans'], k5),
    }

def angle_between_axes(t1, t2):
    d = abs((t2 - t1 + np.pi) % (2*np.pi) - np.pi)
    return np.degrees(min(d, np.pi - d))

print('Cell 1 OK — imports + Sterimol engine loaded.')
'''

CELL2 = r'''# ── Cell 2: Build UI ──────────────────────────────────────────────────────────
_state = {}

def _box(title, *children, color='#4455aa', bg='#e8ebf7'):
    hdr = widgets.HTML(f"<div style='font-weight:600;font-size:13px;color:#fff;"
                       f"background:{color};padding:5px 11px;border-radius:4px 4px 0 0'>{title}</div>")
    return widgets.VBox([hdr, *children],
                        layout=widgets.Layout(border=f'1px solid {color}',
                                              border_radius='4px', padding='0 0 7px 0',
                                              margin='5px 0'))

def _pad(*children):
    return widgets.VBox(list(children), layout=widgets.Layout(padding='7px 10px 0 10px'))

def _itext(desc, val=1, w='165px', dw='92px'):
    return widgets.BoundedIntText(value=val, min=0, max=9999, description=desc,
                                  layout=widgets.Layout(width=w),
                                  style={'description_width': dw})

def _ftext(desc, val=0.0, w='150px', dw='28px'):
    return widgets.FloatText(value=val, description=desc,
                             layout=widgets.Layout(width=w),
                             style={'description_width': dw})

def _chk(desc, val=True, w='180px'):
    return widgets.Checkbox(value=val, description=desc, indent=False,
                            layout=widgets.Layout(width=w))

_ORIG_OPTS = ['atom index', 'atom centroid', 'xyz point']
def _origin_row(label, mode_val='atom index', val='0'):
    m = widgets.Dropdown(options=_ORIG_OPTS, value=mode_val, description=label,
                         layout=widgets.Layout(width='255px'),
                         style={'description_width': '88px'})
    v = widgets.Text(value=val, placeholder='0  /  1,2  /  x y z',
                     layout=widgets.Layout(width='205px'),
                     description='value:', style={'description_width': '42px'})
    return m, v

def _parse_origin(mode, val_str, coords_arr, one_based=True):
    toks = [t.strip() for t in val_str.replace(',', ' ').split() if t.strip()]
    if not toks:
        return np.zeros(3)
    if mode == 'xyz point':
        if len(toks) != 3:
            raise ValueError(f'xyz point needs 3 numbers, got {val_str!r}')
        return np.array([float(t) for t in toks])
    off = 1 if one_based else 0
    idxs = [int(t) - off for t in toks]
    if mode == 'atom centroid':
        return coords_arr[idxs].mean(axis=0)
    return coords_arr[idxs[0]]

# ════════ Molecule ════════
w_file   = widgets.Text(value='', placeholder='path/to/molecule.xyz',
                        layout=widgets.Layout(width='430px'),
                        description='XYZ:', style={'description_width': '40px'})
btn_load = widgets.Button(description='Load', button_style='info',
                          icon='upload', layout=widgets.Layout(width='95px'))
out_load = widgets.Output()
sec_mol = _box('1 · Molecule', _pad(widgets.HBox([w_file, btn_load]), out_load))

# ════════ Sterimol basis ════════
w_oi = _itext('Origin atom:',    1, dw='90px')
w_ai = _itext('Direction atom:', 2, w='180px', dw='105px')
w_pi = _itext('Plane atom (opt):', 0, w='185px', dw='112px')
w_nz   = _itext('z-slices:',  240, w='150px', dw='62px')
w_nth  = _itext('θ-steps:',   360, w='150px', dw='62px')
w_frac = _ftext('B1 frac:', 0.30, w='150px', dw='52px')
w_so_mode, w_so_val = _origin_row('Arrow origin:', 'atom index', '1')
sec_ster = _box('2 · Sterimol basis  (1-based indices)',
                _pad(widgets.HBox([w_oi, w_ai, w_pi]),
                     widgets.HBox([w_nz, w_nth, w_frac]),
                     widgets.HTML("<small style='color:#666'>Plane atom = 0 → auto. "
                                  "Fixes the X/Y azimuth so dipole x/y are reproducible.</small>"),
                     widgets.HBox([w_so_mode, w_so_val])),
                color='#2266aa', bg='#e4eef7')

# ════════ Dipole ════════
w_dsrc = widgets.ToggleButtons(options=['Gaussian log', 'Manual'], value='Gaussian log',
                               style={'button_width': '130px'})
w_log  = widgets.Text(value='', placeholder='path/to/gaussian.log',
                      layout=widgets.Layout(width='470px'),
                      description='Log:', style={'description_width': '40px'})
w_logbox = widgets.VBox([w_log])
w_dx, w_dy, w_dz = _ftext('dx:'), _ftext('dy:'), _ftext('dz:')
w_dt   = _ftext('total:', w='160px', dw='42px')
w_dauto= _chk('auto total', True, '120px')
w_manbox = widgets.VBox([widgets.HBox([w_dx, w_dy, w_dz, w_dt, w_dauto])],
                        layout=widgets.Layout(display='none'))
w_do_mode, w_do_val = _origin_row('Arrow origin:', 'atom index', '1')

def _toggle_dsrc(ch):
    show_log = (ch['new'] == 'Gaussian log')
    w_logbox.layout.display = '' if show_log else 'none'
    w_manbox.layout.display = 'none' if show_log else ''
w_dsrc.observe(_toggle_dsrc, names='value')

sec_dip = _box('3 · Dipole',
               _pad(w_dsrc, w_logbox, w_manbox,
                    widgets.HTML("<small style='color:#666'>Vector is rotated into the "
                                 "Sterimol frame (the 'after basis change' value below).</small>"),
                    widgets.HBox([w_do_mode, w_do_val])),
               color='#8844aa', bg='#f0e8f7')

# ════════ Run ════════
btn_calc = widgets.Button(description='Calculate & Visualize', button_style='success',
                          icon='cogs', layout=widgets.Layout(width='250px', height='38px'))
out_sum = widgets.Output()

# ════════ View options (live) ════════
def _cap(text):
    return widgets.HTML(f"<span style='font-size:11px;font-weight:600;color:#338855'>{text}</span>",
                        layout=widgets.Layout(width='78px'))

def _slider(desc, val, lo, hi, step, dw='70px', w='220px'):
    return widgets.FloatSlider(value=val, min=lo, max=hi, step=step, description=desc,
                               continuous_update=False, readout_format='.2f',
                               layout=widgets.Layout(width=w),
                               style={'description_width': dw})

# -- definition / structure / background
w_mode   = widgets.Dropdown(options=[('Local (per-slice)', 'local'),
                                     ('Global (projected)', 'global')],
                            value='local', description='B1/B5 def:',
                            layout=widgets.Layout(width='235px'),
                            style={'description_width': '72px'})
w_struct = widgets.Dropdown(options=[('Full molecule', 'full'),
                                     ('Substructure', 'sub')],
                            value='full', description='Atoms:',
                            layout=widgets.Layout(width='205px'),
                            style={'description_width': '52px'})
w_bg     = widgets.Dropdown(options=[('White', 'white'), ('Light', '#f0f0f0'),
                                     ('Dark', '#202020'), ('Black', 'black')],
                            value='white', description='Bg:',
                            layout=widgets.Layout(width='150px'),
                            style={'description_width': '28px'})

# -- sterimol glyph toggles
w_v_ster     = _chk('B1/B5/L arrows', True, '145px')
w_v_ring     = _chk('Slice rings', True, '115px')
w_v_arc      = _chk('B1–B5 angle', True, '125px')
w_v_ster_lbl = _chk('Sterimol labels', True, '150px')

# -- atom / molecule
w_v_lbl    = _chk('Atom labels', True, '115px')
w_v_h      = _chk('Show H', True, '85px')
w_lblsize  = widgets.IntSlider(value=11, min=7, max=20, step=1, description='Label sz:',
                               continuous_update=False, layout=widgets.Layout(width='215px'),
                               style={'description_width': '60px'})
w_atomsc   = _slider('Atom size:', 0.21, 0.08, 0.45, 0.01, dw='62px', w='220px')

# -- dipole: 4 independent arrows
w_v_mux    = _chk('μx', True, '62px')
w_v_muy    = _chk('μy', True, '62px')
w_v_muz    = _chk('μz', True, '62px')
w_v_mutot  = _chk('μ total', True, '85px')
w_v_dip_lbl= _chk('μ labels', True, '95px')
w_v_dipfr  = widgets.Dropdown(options=[('Transformed (frame)', 'after'),
                                       ('Raw (Gaussian)', 'before')],
                              value='after', description='μ frame:',
                              layout=widgets.Layout(width='225px'),
                              style={'description_width': '58px'})
w_dipscale = _slider('μ scale:', 1.0, 0.2, 4.0, 0.1, dw='58px', w='230px')

view_bar = _box('4 · View  (updates instantly)',
    _pad(widgets.HBox([w_mode, w_struct, w_bg]),
         widgets.HBox([_cap('Sterimol'), w_v_ster, w_v_ring, w_v_arc, w_v_ster_lbl]),
         widgets.HBox([_cap('Atoms'), w_v_lbl, w_v_h, w_lblsize, w_atomsc]),
         widgets.HBox([_cap('Dipole'), w_v_mux, w_v_muy, w_v_muz, w_v_mutot, w_v_dip_lbl]),
         widgets.HBox([_cap(''), w_v_dipfr, w_dipscale])),
    color='#338855', bg='#e6f3ec')

out_3d  = widgets.Output()
out_2d  = widgets.Output()

ui = widgets.VBox([sec_mol, sec_ster, sec_dip,
                   widgets.HBox([btn_calc], layout=widgets.Layout(padding='4px 0')),
                   out_sum, view_bar,
                   widgets.HBox([out_3d, out_2d])])
display(ui)
print('Cell 2 OK — UI ready. Run Cell 3 to wire callbacks.')
'''

CELL3 = r'''# ── Cell 3: py3Dmol drawing + callbacks ───────────────────────────────────────

# ════════ 3Dmol.js bootstrap (steriplot pattern) ════════
_JS = {'src': None, 'uri': None}
def _get_js():
    if _JS['src'] is None:
        for url in ['https://cdn.jsdelivr.net/npm/3dmol/build/3Dmol-min.js',
                    'https://unpkg.com/3dmol/build/3Dmol-min.js',
                    'https://3dmol.csb.pitt.edu/build/3Dmol-min.js']:
            try:
                with urllib.request.urlopen(url, timeout=10) as r:
                    _JS['src'] = r.read().decode('utf-8', errors='replace')
                print(f'3Dmol.js loaded from {url}'); break
            except Exception:
                pass
        if _JS['src'] is None:
            print('WARNING: could not download 3Dmol.js; 3-D view may not render.')
    return _JS['src']

def _get_uri():
    js = _get_js()
    if js is None:
        return 'https://cdn.jsdelivr.net/npm/3dmol/build/3Dmol-min.js'
    if _JS['uri'] is None:
        _JS['uri'] = 'data:text/javascript;base64,' + base64.b64encode(js.encode()).decode('ascii')
    return _JS['uri']

def _new_view(width=560, height=520):
    uri = _get_uri()
    try:
        return py3Dmol.view(width=width, height=height, js=uri)
    except TypeError:
        return py3Dmol.view(width=width, height=height)

def _patch(html):
    uri = _get_uri()
    reset = ('<script>if(typeof window.$3Dmol==="undefined")'
             '{window.$3Dmolpromise=undefined;}</script>')
    html = re.sub(r'loadScriptAsync\((["\x27])(?:https?:)?//[^"\x27]*3[Dd]mol[^"\x27]*\1\)',
                  lambda m: 'loadScriptAsync(' + repr(uri) + ')', html)
    html = re.sub(r'(<script[^>]+src=)(["\x27])(?:https?:)?//[^"\x27]*3[Dd]mol[^"\x27]*\2',
                  lambda m: m.group(1) + repr(uri), html)
    return reset + html

def _show3d(html):
    html = _patch(html)
    publish_display_data({'text/html': html, 'application/3dmoljs_load.v0': html}, metadata={})

def _xyz_text(atoms):
    return f"{len(atoms)}\n\n" + ''.join(f"{s} {x:.6f} {y:.6f} {z:.6f}\n"
                                          for s, x, y, z in atoms)

# ════════ primitives (LABEL BOX FIX) ════════
_VIEWCFG = {'lbl_sz': 11}      # updated each render from the Label-size slider

def _pt(x, y, z): return {'x': float(x), 'y': float(y), 'z': float(z)}

def _cyl(v, a, b, color, r=0.025):
    v.addCylinder({'start': _pt(*a), 'end': _pt(*b), 'radius': r, 'color': color})

def _sph(v, p, color, r=0.12):
    v.addSphere({'center': _pt(*p), 'radius': r, 'color': color})

def _lbl(v, text, pos, color, sz=None, bg='white'):
    # Fix: explicit white bg + opacity + no border → no black boxes.
    if sz is None:
        sz = _VIEWCFG['lbl_sz']
    v.addLabel(str(text), {
        'position': _pt(*pos), 'backgroundColor': bg, 'backgroundOpacity': 0.7,
        'fontColor': color, 'fontSize': int(sz), 'borderThickness': 0.0,
        'inFront': True, 'showBackground': True,
    })

def _arrow(v, tail, tip, color, r=0.06, rr=2.3, mid=0.9):
    v.addArrow({'start': _pt(*tail), 'end': _pt(*tip), 'color': color,
                'radius': r, 'radiusRatio': rr, 'mid': mid})

# ════════ sterimol glyphs (Z = L axis; offset by arrow origin O) ════════
def _draw_L(v, L, O, show_lbl=True):
    _arrow(v, (O[0], O[1], O[2]), (O[0], O[1], O[2] + L), '#111111', r=0.05)
    if show_lbl:
        _lbl(v, f'L={L:.2f}', (O[0] + 0.15, O[1] + 0.15, O[2] + L), '#111111')

def _draw_ring(v, thetas, radii, z, color, O, r=0.02, n_spokes=18, spoke='#cfe0f5'):
    pts = [(O[0] + ri*np.cos(t), O[1] + ri*np.sin(t), O[2] + z)
           for t, ri in zip(thetas, radii)]
    n = len(pts)
    for i in range(n):
        _cyl(v, pts[i], pts[(i+1) % n], color, r=r)
    if n_spokes:
        step = max(1, n // n_spokes)
        for i in range(0, n, step):
            _cyl(v, (O[0], O[1], O[2] + z), pts[i], spoke, r=0.006)

def _draw_width(v, z, theta, b_val, color, O, label=None):
    ux, uy = np.cos(theta), np.sin(theta)
    end = (O[0] + b_val*ux, O[1] + b_val*uy, O[2] + z)
    _cyl(v, (O[0], O[1], O[2] + z), end, color, r=0.055)
    _sph(v, end, color, r=0.12)
    _sph(v, (O[0], O[1], O[2] + z), color, r=0.08)
    if label:
        _lbl(v, label, (end[0] + 0.1, end[1] + 0.1, end[2]), color)

def _draw_tick(v, z, color, O):
    t = 0.32
    _cyl(v, (O[0]-t, O[1], O[2]+z), (O[0]+t, O[1], O[2]+z), color, r=0.03)
    _cyl(v, (O[0], O[1]-t, O[2]+z), (O[0], O[1]+t, O[2]+z), color, r=0.03)

def _draw_arc(v, theta1, theta2, radius, color, O, z=0.0, n=40, show_lbl=True):
    d = (theta2 - theta1 + np.pi) % (2*np.pi) - np.pi
    if abs(d) > np.pi/2:
        d -= np.sign(d) * np.pi
    ts = np.linspace(theta1, theta1 + d, n)
    pts = [(O[0] + radius*np.cos(t), O[1] + radius*np.sin(t), O[2] + z) for t in ts]
    for i in range(len(pts)-1):
        _cyl(v, pts[i], pts[i+1], color, r=0.026)
    if show_lbl:
        mt = theta1 + d/2
        _lbl(v, f'{angle_between_axes(theta1, theta2):.1f}°',
             (O[0] + 1.35*radius*np.cos(mt), O[1] + 1.35*radius*np.sin(mt), O[2] + z), color)

def _hl_atom(v, atoms, disk, color):
    if disk is None:
        return
    a = atoms[disk['atom_idx']]
    _sph(v, (a[1], a[2], a[3]), color, r=VDW.get(a[0], 1.7) * 0.30)

# ════════ dipole arrows  (4 independent: μx, μy, μz, μ total) ════════
_MU_X, _MU_Y, _MU_Z, _MU_T = '#e23b3b', '#1faa4a', '#2f6fe2', '#8a2be2'

def _draw_dipole(v, vec, origin, span, scale_mult, show, show_lbl):
    dx, dy, dz = float(vec[0]), float(vec[1]), float(vec[2])
    mag = float(np.linalg.norm([dx, dy, dz]))
    if mag < 1e-9:
        return
    O = origin
    # total arrow length is a fixed fraction of the molecule, tuned by the slider;
    # components share the same scale so they stay proportional & comparable.
    base = 0.30 * span * float(scale_mult)
    sc = base / mag

    comps = [('x', (dx*sc, 0.0, 0.0), _MU_X, dx),
             ('y', (0.0, dy*sc, 0.0), _MU_Y, dy),
             ('z', (0.0, 0.0, dz*sc), _MU_Z, dz)]
    for nm, (cx, cy, cz), col, val in comps:
        if not show.get(nm):
            continue
        if abs(val) < 1e-6:
            continue
        tip = (O[0]+cx, O[1]+cy, O[2]+cz)
        _arrow(v, O, tip, col, r=0.055, rr=2.3)
        _sph(v, tip, col, r=0.11)
        if show_lbl:
            _lbl(v, f'μ{nm}={val:.2f}', (tip[0]+0.12, tip[1]+0.12, tip[2]+0.12), col, sz=10)

    if show.get('tot'):
        tip = (O[0]+dx*sc, O[1]+dy*sc, O[2]+dz*sc)
        _arrow(v, O, tip, _MU_T, r=0.09, rr=2.7)
        _sph(v, tip, _MU_T, r=0.14)
        if show_lbl:
            _lbl(v, f'μ={mag:.2f} D', (tip[0]+0.14, tip[1]+0.14, tip[2]+0.14), _MU_T, sz=11)

    _sph(v, O, '#333333', r=0.10)   # origin marker

# ════════ 3-D render (reads _state + live view widgets) ════════
def _render_3d():
    d = _state
    if 'atoms' not in d:
        return '<p style="color:#888">No results yet — press Calculate &amp; Visualize.</p>'
    atoms   = d['atoms']
    sub     = d['sub_indices']
    struct  = w_struct.value
    mode    = w_mode.value
    res     = d['res_sub'] if struct == 'sub' else d['res_full']
    glb     = d['glb_sub'] if struct == 'sub' else d['glb_full']
    span    = d['span']

    _VIEWCFG['lbl_sz'] = int(w_lblsize.value)
    asc = float(w_atomsc.value)
    slbl = w_v_ster_lbl.value

    v = _new_view()
    v.setBackgroundColor(w_bg.value)
    v.addModel(_xyz_text(atoms), 'xyz')
    v.setStyle({'stick': {'radius': 0.11}, 'sphere': {'scale': asc}})
    if struct == 'sub':
        v.setStyle({}, {'stick': {'radius': 0.07, 'color': 'gray'},
                        'sphere': {'scale': asc*0.66, 'color': 'gray'}})
        v.addStyle({'serial': [i+1 for i in sub]},
                   {'stick': {'radius': 0.14}, 'sphere': {'scale': asc*1.25}})
    if not w_v_h.value:
        v.setStyle({'elem': 'H'}, {'stick': {'radius': 0.04, 'color': '#cfcfcf'},
                                   'sphere': {'scale': asc*0.38, 'color': '#cfcfcf'}})

    if w_v_lbl.value:
        lcol = '#eeeeee' if w_bg.value in ('black', '#202020') else '#222222'
        lbg  = '#222222' if w_bg.value in ('black', '#202020') else 'white'
        for i, (sym, x, y, z) in enumerate(atoms):
            if sym == 'H' and not w_v_h.value:
                continue
            if struct == 'sub' and i not in sub:
                continue
            _lbl(v, f'{sym}{i+1}', (x+0.08, y+0.08, z+0.08), lcol, bg=lbg)

    Os = d['ster_origin']
    if w_v_ster.value:
        _draw_L(v, res['L'], Os, show_lbl=slbl)
        if mode == 'global':
            _draw_width(v, 0.0, glb['theta_B1'], glb['B1'], '#4169e1', Os,
                        f"B1={glb['B1']:.2f}" if slbl else None)
            _draw_width(v, 0.0, glb['theta_B5'], glb['B5'], '#dc143c', Os,
                        f"B5={glb['B5']:.2f}" if slbl else None)
            _hl_atom(v, atoms, glb['plus_atom_B1'], '#4169e1')
            _hl_atom(v, atoms, glb['plus_atom_B5'], '#dc143c')
            if w_v_arc.value:
                _draw_arc(v, glb['theta_B1'], glb['theta_B5'],
                          0.45*max(glb['B1'], glb['B5']), '#444', Os, show_lbl=slbl)
        else:
            if w_v_ring.value:
                _draw_ring(v, res['scan_B1']['thetas'], res['scan_B1']['radii'],
                           res['z_B1'], '#4169e1', Os, spoke='#bcd4f0')
                _draw_ring(v, res['scan_B5']['thetas'], res['scan_B5']['radii'],
                           res['z_B5'], '#dc143c', Os, spoke='#f5b8a8')
            _draw_tick(v, res['z_B1'], '#4169e1', Os)
            _draw_tick(v, res['z_B5'], '#dc143c', Os)
            _draw_width(v, res['z_B1'], res['theta_B1'], res['B1'], '#4169e1', Os,
                        f"B1={res['B1']:.2f}" if slbl else None)
            _draw_width(v, res['z_B5'], res['theta_B5'], res['B5'], '#dc143c', Os,
                        f"B5={res['B5']:.2f}" if slbl else None)
            _hl_atom(v, atoms, res['plus_atom_B1'], '#4169e1')
            _hl_atom(v, atoms, res['plus_atom_B5'], '#dc143c')
            if w_v_arc.value:
                _draw_arc(v, res['theta_B1'], res['theta_B5'],
                          0.45*max(res['B1'], res['B5']), '#444', Os, show_lbl=slbl)

    dip_show = {'x': w_v_mux.value, 'y': w_v_muy.value,
                'z': w_v_muz.value, 'tot': w_v_mutot.value}
    if any(dip_show.values()) and d.get('dip_after') is not None:
        vec = d['dip_after'] if w_v_dipfr.value == 'after' else d['dip_before']
        _draw_dipole(v, vec, d['dip_origin'], span,
                     scale_mult=w_dipscale.value, show=dip_show,
                     show_lbl=w_v_dip_lbl.value)

    v.zoomTo()
    return v._make_html()

def _update_view(change=None):
    if 'atoms' not in _state:
        return
    out_3d.clear_output()
    with out_3d:
        _show3d(_render_3d())
    out_2d.clear_output()
    with out_2d:
        _render_2d()

for _w in [w_mode, w_struct, w_bg,
           w_v_ster, w_v_ring, w_v_arc, w_v_ster_lbl,
           w_v_lbl, w_v_h, w_lblsize, w_atomsc,
           w_v_mux, w_v_muy, w_v_muz, w_v_mutot, w_v_dip_lbl, w_v_dipfr, w_dipscale]:
    _w.observe(_update_view, names='value')

# ════════ 2-D end-on projection (steriplot port) ════════
def _caliper(ax, theta, b_val, color, label, atom=None):
    ux, uy = np.cos(theta), np.sin(theta)
    nx, ny = -uy, ux
    pe = np.array([b_val*ux, b_val*uy])
    t = np.linspace(0, 2*np.pi, 240)
    ax.plot(b_val*np.cos(t), b_val*np.sin(t), '--', color=color, lw=0.8, alpha=0.4, zorder=1)
    jaw = b_val*0.12
    ax.plot([pe[0]-jaw*nx, pe[0]+jaw*nx], [pe[1]-jaw*ny, pe[1]+jaw*ny], '-', color=color, lw=1.6, zorder=6)
    ax.annotate('', xy=tuple(pe), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color=color, lw=1.7, mutation_scale=12))
    ax.plot(0, 0, 'k+', ms=9, mew=1.5, zorder=10)
    astr = f"  ({atom['sym']}{atom['atom_idx']})" if atom else ''
    return mlines.Line2D([], [], color=color, lw=2, marker='>', markersize=7,
                         markerfacecolor=color, label=f'{label} = {b_val:.2f} Å{astr}')

def _render_2d():
    d = _state
    if 'atoms' not in d:
        return
    atoms = d['atoms']; sub = d['sub_indices']
    struct = w_struct.value; mode = w_mode.value
    res = d['res_sub'] if struct == 'sub' else d['res_full']
    glb = d['glb_sub'] if struct == 'sub' else d['glb_full']
    g_th, g_rd = (d['env_sub'] if struct == 'sub' else d['env_full'])
    r_env = float(np.max(g_rd))

    fig, ax = plt.subplots(figsize=(5.0, 5.0))
    fig.patch.set_facecolor('#f8f8f8'); ax.set_facecolor('#f8f8f8')
    for i, (sym, x, y, _z) in enumerate(atoms):
        r = VDW.get(sym)
        if r is None:
            continue
        in_sub = i in sub
        if struct == 'sub' and not in_sub:
            fc, ec, fa, ea = '#e4e4e4', '#c0c0c0', 0.22, 0.32
        else:
            fc = ELEM_COLOR.get(sym, '#888888'); ec = fc; fa, ea = 0.20, 0.50
        ax.add_patch(plt.Circle((x, y), r, fc=fc, ec='none', alpha=fa, zorder=2))
        ax.add_patch(plt.Circle((x, y), r, fill=False, ec=ec, lw=0.55, alpha=ea, zorder=2))
        if sym != 'H' and (struct == 'full' or in_sub):
            ax.text(x, y, f'{sym}{i+1}', fontsize=6, ha='center', va='center',
                    color='#111', alpha=0.7, zorder=10)
    gx, gy = g_rd*np.cos(g_th), g_rd*np.sin(g_th)
    ax.plot(np.append(gx, gx[0]), np.append(gy, gy[0]), '--', color='#aaa', lw=1.0, alpha=0.7, zorder=3)
    handles = [mlines.Line2D([], [], color='#888', lw=1, ls='--', label='global envelope')]
    if mode == 'global':
        h1 = _caliper(ax, glb['theta_B1'], glb['B1'], '#4169e1', 'B1', glb['plus_atom_B1'])
        h5 = _caliper(ax, glb['theta_B5'], glb['B5'], '#dc143c', 'B5', glb['plus_atom_B5'])
        t1, t5, ttl = glb['theta_B1'], glb['theta_B5'], 'global projection'
    else:
        for scan, col, txt in [(res['scan_B1'], '#4169e1', f"B1 slice z={res['z_B1']:.2f}"),
                               (res['scan_B5'], '#dc143c', f"B5 slice z={res['z_B5']:.2f}")]:
            ex, ey = scan['radii']*np.cos(scan['thetas']), scan['radii']*np.sin(scan['thetas'])
            ax.fill(ex, ey, alpha=0.09, color=col, zorder=3)
            ax.plot(np.append(ex, ex[0]), np.append(ey, ey[0]), color=col, lw=1.3, zorder=4)
            handles.append(mlines.Line2D([], [], color=col, lw=1.5, label=txt))
        h1 = _caliper(ax, res['theta_B1'], res['B1'], '#4169e1', 'B1', res['plus_atom_B1'])
        h5 = _caliper(ax, res['theta_B5'], res['B5'], '#dc143c', 'B5', res['plus_atom_B5'])
        t1, t5, ttl = res['theta_B1'], res['theta_B5'], 'local slices'
    handles += [h1, h5]
    ang = angle_between_axes(t1, t5)
    _d = (t5 - t1 + np.pi) % (2*np.pi) - np.pi
    if abs(_d) > np.pi/2:
        _d -= np.sign(_d)*np.pi
    ar = 0.28*r_env
    at = np.linspace(t1, t1+_d, 80)
    ax.plot(ar*np.cos(at), ar*np.sin(at), '-', color='#666', lw=1.5, zorder=5)
    handles.append(mlines.Line2D([], [], color='#666', lw=1.5, label=f'angle = {ang:.1f}°'))
    ax.axhline(0, color='k', lw=0.3, alpha=0.18); ax.axvline(0, color='k', lw=0.3, alpha=0.18)
    ax.legend(handles=handles, fontsize=7, loc='upper right', framealpha=0.92, edgecolor='#ccc')
    _s = '  ·  sub' if struct == 'sub' else ''
    ax.set_title(f'End-on view (along L)\n{ttl}{_s}', fontsize=9, pad=5)
    ax.set_xlabel('x (Å)', fontsize=9); ax.set_ylabel('y (Å)', fontsize=9)
    lim = r_env*1.32
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
    plt.tight_layout(pad=0.6); plt.show()

# ════════ Load ════════
def on_load(b):
    out_load.clear_output()
    with out_load:
        try:
            path = w_file.value.strip()
            if not os.path.isfile(path):
                print(f'File not found:\n  {path}'); return
            raw = read_xyz(path)
            _state['raw_atoms'] = raw
            n = len(raw)
            for w in [w_oi, w_ai, w_pi]:
                w.max = n
            print(f'Loaded {n} atoms.')
            df = pd.DataFrame([(i+1, s, x, y, z) for i, (s, x, y, z) in enumerate(raw)],
                              columns=['idx', 'atom', 'x', 'y', 'z']).set_index('idx')
            display(df.style.set_caption('Atoms (1-based)')
                    .format({'x': '{:.4f}', 'y': '{:.4f}', 'z': '{:.4f}'}))
        except Exception as e:
            print(f'Load error: {e}'); _tb.print_exc()
btn_load.on_click(on_load)

# ════════ Calculate ════════
def on_calc(b):
    out_sum.clear_output(); out_3d.clear_output(); out_2d.clear_output()
    raw = _state.get('raw_atoms')
    if raw is None:
        with out_sum: print('Load a molecule first.'); return
    try:
        oi, ai = w_oi.value - 1, w_ai.value - 1          # 1-based → 0-based
        pi = (w_pi.value - 1) if w_pi.value and w_pi.value > 0 else None
        symbols = [a[0] for a in raw]
        coords  = np.array([[a[1], a[2], a[3]] for a in raw], float)

        R, origin = build_frame(coords, oi, ai, pi)
        atoms = apply_frame(symbols, coords, R, origin)
        new_coords = np.array([[a[1], a[2], a[3]] for a in atoms], float)
        span = float(np.ptp(new_coords, axis=0).max())

        sub_indices = find_substructure_atoms(atoms, oi, ai)
        nz, nth, frac = int(w_nz.value), int(w_nth.value), float(w_frac.value)

        res_full = scan_sterimol(atoms, nz, nth, frac, None)
        res_sub  = scan_sterimol(atoms, nz, nth, frac, sub_indices)
        env_full = compute_global_envelope(res_full['all_scans'], nth)
        env_sub  = compute_global_envelope(res_sub['all_scans'], nth)
        glb_full = make_global(res_full, *env_full)
        glb_sub  = make_global(res_sub,  *env_sub)
    except Exception as e:
        with out_sum: print(f'Sterimol error: {e}'); _tb.print_exc()
        return

    # ── dipole ──
    dip_before = dip_after = None
    dip_raw_df = None
    try:
        if w_dsrc.value == 'Gaussian log':
            lp = w_log.value.strip()
            if not os.path.isfile(lp):
                raise FileNotFoundError(f'Log not found: {lp}')
            with open(lp, errors='replace') as fh:
                dip_raw_df = process_gaussian_dipole_text(fh.read())
            dvec = dip_raw_df[['dip_x', 'dip_y', 'dip_z']].to_numpy(float)[0]
        else:
            dvec = np.array([w_dx.value, w_dy.value, w_dz.value], float)
        dip_before = dvec.copy()              # Gaussian / input frame
        dip_after  = R @ dvec                 # rotated into Sterimol frame
    except Exception as e:
        with out_sum: print(f'Dipole error: {e}'); _tb.print_exc()

    # ── arrow origins (display frame) ──
    try:
        ster_origin = _parse_origin(w_so_mode.value, w_so_val.value, new_coords)
        dip_origin  = _parse_origin(w_do_mode.value, w_do_val.value, new_coords)
    except Exception as e:
        with out_sum: print(f'Origin error: {e}'); return

    _state.update(atoms=atoms, sub_indices=sub_indices, span=span, R=R,
                  res_full=res_full, res_sub=res_sub,
                  glb_full=glb_full, glb_sub=glb_sub,
                  env_full=env_full, env_sub=env_sub,
                  dip_before=dip_before, dip_after=dip_after,
                  ster_origin=ster_origin, dip_origin=dip_origin)

    # ── summary ──
    with out_sum:
        def _row(tag, r, g):
            return (f"  {tag:<3} local={r:6.3f}   global={g:6.3f}")
        print('Sterimol (Å)        full molecule          substructure')
        print('                  local   global        local   global')
        print(f"  B1            {res_full['B1']:7.3f} {glb_full['B1']:7.3f}      "
              f"{res_sub['B1']:7.3f} {glb_sub['B1']:7.3f}")
        print(f"  B5            {res_full['B5']:7.3f} {glb_full['B5']:7.3f}      "
              f"{res_sub['B5']:7.3f} {glb_sub['B5']:7.3f}")
        print(f"  L             {res_full['L']:7.3f}              {res_sub['L']:7.3f}")
        print(f"  ∠(B1,B5)      {angle_between_axes(res_full['theta_B1'], res_full['theta_B5']):6.1f}°"
              f"             {angle_between_axes(res_sub['theta_B1'], res_sub['theta_B5']):6.1f}°")
        sub_lbl = ', '.join(f'{atoms[i][0]}{i+1}' for i in sorted(sub_indices))
        print(f"\n  Substructure ({len(sub_indices)}/{len(atoms)}): {sub_lbl}")
        if dip_after is not None:
            print('\nDipole (Debye)          x        y        z      total')
            b = dip_before; a = dip_after
            tot = float(np.linalg.norm(a))
            print(f"  before basis     {b[0]:8.3f} {b[1]:8.3f} {b[2]:8.3f} {np.linalg.norm(b):8.3f}")
            print(f"  after  basis     {a[0]:8.3f} {a[1]:8.3f} {a[2]:8.3f} {tot:8.3f}")

    with out_3d:
        _show3d(_render_3d())
    with out_2d:
        _render_2d()

btn_calc.on_click(on_calc)
_get_js()      # pre-fetch so first render is instant
print('Cell 3 OK — callbacks wired. Use the panel above.')
'''

MD = (
    "# Sterimol + Dipole Interactive Viewer\n\n"
    "Interactive **py3Dmol** view of **Sterimol** (B1/B5/L) and the molecular **dipole**, "
    "with the full steriplot feature set:\n\n"
    "- **Local** (per-slice) vs **Global** (projected envelope) B1/B5\n"
    "- **Full molecule** vs **Substructure** (BFS from the direction atom)\n"
    "- Dipole shown **before** and **after** the basis change (all of x/y/z/total)\n"
    "- **Interactive origin** for both the Sterimol arrows and the dipole arrows\n\n"
    "### Run order\n"
    "1. **Cell 1** – imports + Sterimol engine  \n"
    "2. **Cell 2** – build the UI  \n"
    "3. **Cell 3** – wire callbacks  \n\n"
    "Then: load an `.xyz` → set basis atoms / dipole source → **Calculate & Visualize**. "
    "The **View** panel re-renders instantly without recomputing.\n\n"
    "**Frame:** origin atom → (0,0,0), direction atom → **+Z = L axis**, "
    "optional plane atom fixes the X/Y azimuth (so dipole x/y are reproducible)."
)

def src(text):
    lines = text.split('\n')
    return [l + '\n' for l in lines[:-1]] + [lines[-1]]

nb = {
    'cells': [
        {'cell_type': 'markdown', 'metadata': {}, 'source': src(MD)},
        {'cell_type': 'code', 'execution_count': None, 'metadata': {}, 'outputs': [], 'source': src(CELL1)},
        {'cell_type': 'code', 'execution_count': None, 'metadata': {}, 'outputs': [], 'source': src(CELL2)},
        {'cell_type': 'code', 'execution_count': None, 'metadata': {}, 'outputs': [], 'source': src(CELL3)},
    ],
    'metadata': {
        'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
        'language_info': {'name': 'python', 'version': '3.9.0'},
    },
    'nbformat': 4, 'nbformat_minor': 4,
}

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sterimol_dipole_viewer.ipynb')
with open(out, 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
print('wrote', out)

# syntax-check each code cell
import ast
for i, c in enumerate(nb['cells']):
    if c['cell_type'] == 'code':
        code = ''.join(c['source'])
        ast.parse(code)
        print(f'cell {i}: syntax OK ({len(code)} chars)')
