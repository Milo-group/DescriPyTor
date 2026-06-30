import plotly.graph_objs as go
from plotly.offline import init_notebook_mode, iplot
import sys
import pandas as pd
from typing import *
from enum import Enum
import os
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
# Add the parent directory to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from utils  import help_functions as hf
except Exception as e :
    from .utils import help_functions as hf
# Now you can import from the parent directory



class GeneralConstants(Enum):
    """
    Holds constants for calculations and conversions
    1. covalent radii from Alvarez (2008) DOI: 10.1039/b801115j
    2. atomic numbers
    2. atomic weights
    """
    COVALENT_RADII= {
            'H': 0.31, 'He': 0.28, 'Li': 1.28,
            'Be': 0.96, 'B': 0.84, 'C': 0.76, 
            'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
            'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 
            'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
            'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 
            'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 
            'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 
            'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 
            'Br': 1.20, 'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95,
            'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
            'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39,
            'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39,
            'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40,
            'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04,
            'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98,
            'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92,
            'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
            'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62,
            'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36,
            'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46,
            'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50, 
            'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06,
            'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87,
            'Am': 1.80, 'Cm': 1.69
    }
    BONDI_RADII={
        'H': 1.10, 'C': 1.70, 'F': 1.47,
        'S': 1.80, 'B': 1.92, 'I': 1.98, 
        'N': 1.55, 'O': 1.52, 'Co': 2.00, 
        'Br': 1.83, 'Si': 2.10,'Ni': 2.00,
        'P': 1.80, 'Cl': 1.75, 
    }

def flatten_list(nested_list_arg: List[list]) -> List:
    """
    Flatten a nested list.
    turn [[1,2],[3,4]] to [1,2,3,4]
    """
    flat_list=[item for sublist in nested_list_arg for item in sublist]
    return flat_list

import re
from itertools import combinations

def plot_interactions(xyz_df, color, dipole_df=None, origin=None, sterimol_params=None):
    """
    Creates a 3D Plotly figure of the molecule (atoms + bonds) with optional dipole and Sterimol arrows.
    - Locks aspect to 'data' (no fake warping of small z).
    - Uses one coordinate array for all traces.
    - Correct visibility masks for toggles.
    """
    # ---------- helpers ----------
    def _parse_element(label: str) -> str:
        m = re.match(r"[A-Za-z]+", str(label).strip())
        if not m:
            s = str(label).strip()
        else:
            s = m.group(0)
        return s[0].upper() + s[1:].lower() if s else s

    def _resolve_origin(origin_arg, coords_arr):
        """
        origin_arg can be:
          - None: centroid of coords_arr
          - array-like shape (3,): xyz position
          - int: atom index (0-based)
          - list/array of ints: centroid of those atoms
        """
        if origin_arg is None:
            return coords_arr.mean(axis=0)
        ori = np.asarray(origin_arg)
        if ori.ndim == 1 and ori.shape[0] == 3 and np.isfinite(ori).all():
            return ori.astype(float)
        if np.issubdtype(ori.dtype, np.integer) and ori.ndim == 0:
            return coords_arr[int(ori)]
        if np.issubdtype(ori.dtype, np.integer) and ori.ndim == 1:
            return coords_arr[ori].mean(axis=0)
        raise ValueError("Invalid 'origin' argument. Use None, (3,), int, or list of ints.")

    def _planarity_metrics(coords_arr):
        """Optional diagnostic—unused in layout, handy to print if needed."""
        c = coords_arr.mean(0)
        X = coords_arr - c
        _, _, Vt = np.linalg.svd(X, full_matrices=False)
        n = Vt[-1]
        d = X @ n
        return dict(rms=float(np.sqrt((d**2).mean())), max_abs=float(np.max(np.abs(d))))

    # ---------- constants / lookups ----------
    try:
        atomic_radii_raw = GeneralConstants.COVALENT_RADII.value  # user env
        try:
            atomic_radii = dict(atomic_radii_raw)
        except Exception:
            atomic_radii = atomic_radii_raw
    except Exception:
        # minimal fallback
        atomic_radii = {"H": 0.31, "B": 0.85, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
                        "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
                        "Pd": 1.39, "Co": 1.26, "Ni": 1.24, "Fe": 1.24, "Cu": 1.32, "Zn": 1.22,
                        "Ag": 1.45, "Au": 1.36}

    cpk_colors = dict(
        C='black', F='green', H='white', N='blue', O='red', P='orange',
        S='yellow', Cl='green', Br='brown', I='purple',
        Ni='blue', Fe='red', Cu='orange', Zn='yellow', Ag='grey',
        Au='gold', Si='grey', B='pink', Pd='green', Co='pink'
    )

    # ---------- coordinates / atoms ----------
    coords = xyz_df[['x', 'y', 'z']].to_numpy(dtype=float)
    atoms_raw = xyz_df['atom'].astype(str).tolist()
    elements = [_parse_element(a) for a in atoms_raw]

    n = coords.shape[0]
    assert n == len(elements), f"Row count mismatch: coords={n}, atoms={len(elements)}"
    if not np.isfinite(coords).all():
        bad = np.argwhere(~np.isfinite(coords)).ravel().tolist()
        raise ValueError(f"Non-finite coordinates at rows: {bad}")

    # radii (safe defaults)
    default_r = 0.77
    radii = np.array([atomic_radii.get(el, default_r) for el in elements], dtype=float)

    # ---------- bonds (robust O(n^2)) ----------
    def get_bonds(thresh_scale=1.30, min_dist=0.10):
        bonds_local = {}
        for i, j in combinations(range(n), 2):
            dij = np.linalg.norm(coords[i] - coords[j])
            if dij <= min_dist:
                continue
            cutoff = (radii[i] + radii[j]) * thresh_scale
            if dij < cutoff:
                bonds_local[(i, j)] = round(dij, 2)
        return bonds_local

    bonds = get_bonds()

    # ---------- base traces ----------
    atom_colors = [cpk_colors.get(el, 'gray') for el in elements]
    hovertext = [f"{el} ({x:.3f},{y:.3f},{z:.3f})"
                 for el, (x, y, z) in zip(elements, coords)]
    atom_scatter = go.Scatter3d(
        x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
        mode='markers',
        marker=dict(color=atom_colors, size=5, line=dict(color='lightgray', width=2)),
        text=hovertext, name='atoms', hoverinfo='text'
    )

    bx, by, bz = [], [], []
    for (i, j) in bonds:
        xi, yi, zi = coords[i]
        xj, yj, zj = coords[j]
        bx += [xi, xj, None]
        by += [yi, yj, None]
        bz += [zi, zj, None]
    bond_trace = go.Scatter3d(
        x=bx, y=by, z=bz, mode='lines',
        line=dict(color=color, width=3), hoverinfo='none',
        name='bonds'
    )

    # ---------- annotation sets ----------
    annotations_idx = [
        dict(text=str(i + 1), x=coords[i, 0], y=coords[i, 1], z=coords[i, 2],
             showarrow=False, yshift=15, font=dict(color="blue"))
        for i in range(n)
    ]
    annotations_len = [
        dict(text=str(dist),
             x=(coords[i, 0] + coords[j, 0]) / 2,
             y=(coords[i, 1] + coords[j, 1]) / 2,
             z=(coords[i, 2] + coords[j, 2]) / 2,
             showarrow=False, yshift=10)
        for (i, j), dist in bonds.items()
    ]

    def add_traces():
        traces = []
        coords = xyz_df[['x','y','z']].to_numpy(float)
        center = coords.mean(axis=0)
        span = np.ptp(coords, axis=0)
        pad = np.median(span)*0.5   # how far beyond molecule you want to expand

        # create invisible boundary points
        x_dummy = [center[0] - pad, center[0] + pad]
        y_dummy = [center[1] - pad, center[1] + pad]
        z_dummy = [center[2] - pad, center[2] + pad]

        invisible_trace = go.Scatter3d(
            x=x_dummy,
            y=y_dummy,
            z=z_dummy,
            mode="markers",
            marker=dict(size=0.1, opacity=0),  # completely invisible
            hoverinfo="none",
            showlegend=False,
            name="frame_expander",
        )
        traces.append(invisible_trace)
            
        return traces
    
    # ---------- dipole arrows ----------
    def dipole_traces(dip_df, origin_arg, coords_for_scale, show_components=True):
        traces = []
        if dip_df is None or len(dip_df) == 0:
            return traces
        row = dip_df.iloc[0]
        vec = np.array([
            row.get("dipole_x", row.iloc[0]),
            row.get("dipole_y", row.iloc[1]),
            row.get("dipole_z", row.iloc[2]),
        ], dtype=float)
        if not np.all(np.isfinite(vec)):
            return traces

        tail = _resolve_origin(origin_arg, coords_for_scale)

        # auto-scale w.r.t. molecular span
        span = float(np.ptp(coords_for_scale, axis=0).max())
        dip_mag = float(np.linalg.norm(vec))
        scale = ((span / dip_mag) * 0.3 if dip_mag > 0 else 1.0) * 2.0

        comps = []
        if show_components:
            comps.extend([
                (np.array([vec[0] * scale, 0.0, 0.0]), "red",   "Dipole X"),
                (np.array([0.0, vec[1] * scale, 0.0]), "green", "Dipole Y"),
                (np.array([0.0, 0.0, vec[2] * scale]), "blue",  "Dipole Z"),
            ])
        comps.append((vec * scale, "purple", "Total Dipole"))

        for v, col, label in comps:
            L = float(np.linalg.norm(v))
            if L < 1e-10:
                continue
            end = tail + v
            shaft_end = end - 0.15 * (end - tail) / L

            traces.append(go.Scatter3d(
                x=[tail[0], shaft_end[0]],
                y=[tail[1], shaft_end[1]],
                z=[tail[2], shaft_end[2]],
                mode="lines",
                line=dict(color=col, width=4),
                name=label,
                showlegend=True,
            ))
            traces.append(go.Cone(
                x=[end[0]], y=[end[1]], z=[end[2]],
                u=[v[0]], v=[v[1]], w=[v[2]],
                anchor="tip",
                sizemode="scaled",
                sizeref=0.2,
                showscale=False,
                colorscale=[[0, col], [1, col]],
                name=label,
                showlegend=False,
            ))

        # mark origin
        traces.append(go.Scatter3d(
            x=[tail[0]], y=[tail[1]], z=[tail[2]],
            mode="markers",
            marker=dict(size=5, color="black", symbol="circle"),
            name="Origin (tail)"
        ))
        return traces

    # ---------- Sterimol arrows ----------
    def sterimol_traces(params, origin_arg):
        """
        Expects 'B1_coords','B5_coords','L_coords' as 2D vectors (x,y) in current global axes.
        Anchors at resolved 3D origin; lifts vectors with z=0 in global frame.
        """
        if params is None:
            return []
        required = ['B1_coords', 'B5_coords', 'L_coords', 'B1_value', 'B5_value', 'L_value']
        if not all(k in params for k in required):
            return []

        try:
            tail = _resolve_origin(origin_arg, coords)
        except Exception:
            tail = coords.mean(axis=0)

        vecs = [
            (np.asarray(params['B1_coords'], float), 'forestgreen', params['B1_value'], 'B1'),
            (np.asarray(params['B5_coords'], float), 'firebrick',   params['B5_value'], 'B5'),
            (np.asarray(params['L_coords'],  float), 'steelblue',   params['L_value'], 'L'),
        ]

        traces = []
        for vec2d, col, mag, name in vecs:
            vec3d = np.array([vec2d[0], vec2d[1], 0.0], float)
            L = float(np.linalg.norm(vec3d))
            if L < 1e-8:
                continue
            end = tail + vec3d
            shaft_end = end - 0.1 * (vec3d / L)

            traces.append(go.Scatter3d(
                x=[tail[0], shaft_end[0]],
                y=[tail[1], shaft_end[1]],
                z=[tail[2], shaft_end[2]],
                mode='lines',
                line=dict(color=col, width=4),
                name=f"{name}-shaft",
                showlegend=True
            ))
            traces.append(go.Cone(
                x=[end[0]], y=[end[1]], z=[end[2]],
                u=[vec3d[0]], v=[vec3d[1]], w=[vec3d[2]],
                anchor='tip',
                sizemode='absolute',
                sizeref=L * 0.07,   # head ~7% of length
                showscale=False,
                colorscale=[[0, col], [1, col]],
                name=name
            ))
        return traces

    # ---------- assemble data ----------
    data = [bond_trace, atom_scatter]

    dip_trs = dipole_traces(dipole_df, origin, coords_for_scale=coords) if dipole_df is not None else []
    for tr in dip_trs:
        tr.visible = False
    data.extend(dip_trs)

    add_trace=add_traces()
    data.extend(add_trace)
    for tr in add_trace:
        tr.visible = True
    ster_trs = sterimol_traces(sterimol_params, origin) if sterimol_params is not None else []
    for tr in ster_trs:
        tr.visible = False
    data.extend(ster_trs)

    # ---------- visibility masks ----------
    n_base = 2
    n_dip = len(dip_trs)
    n_ster = len(ster_trs)
    n_total = n_base + n_dip + n_ster

    vis_base_only = [True]*n_base + [False]*n_dip + [False]*n_ster
    vis_dip_only  = [True]*n_base + [True ]*n_dip + [False]*n_ster
    vis_ster_only = [True]*n_base + [False]*n_dip + [True ]*n_ster
    vis_both      = [True]*n_base + [True ]*n_dip + [True ]*n_ster
    # make add_trace vis
    # vis_add_only  = [True]*n_base + [False]*n_dip + [False]*n_ster + [True ]*n_add
    # ---------- buttons ----------
    buttons = [
        dict(label='Atom indices', method='relayout', args=[{'scene.annotations': annotations_idx}]),
        dict(label='Bond lengths', method='relayout', args=[{'scene.annotations': annotations_len}]),
        dict(label='Both',         method='relayout', args=[{'scene.annotations': annotations_idx + annotations_len}]),
        dict(label='Hide annotations', method='relayout', args=[{'scene.annotations': []}]),
    ]

    if n_dip > 0 and n_ster > 0:
        buttons.extend([
            dict(label='Show dipole',   method='update', args=[{'visible': vis_dip_only},  {}]),
            dict(label='Show Sterimol', method='update', args=[{'visible': vis_ster_only}, {}]),
            dict(label='Show both',     method='update', args=[{'visible': vis_both},      {}]),
            dict(label='Hide arrows',   method='update', args=[{'visible': vis_base_only}, {}]),
        ])
    elif n_dip > 0:
        buttons.extend([
            dict(label='Show dipole', method='update', args=[{'visible': vis_dip_only},  {}]),
            dict(label='Hide dipole', method='update', args=[{'visible': vis_base_only}, {}]),
        ])
    elif n_ster > 0:
        buttons.extend([
            dict(label='Show Sterimol', method='update', args=[{'visible': vis_ster_only}, {}]),
            dict(label='Hide Sterimol', method='update', args=[{'visible': vis_base_only}, {}]),
        ])

    updatemenus = [dict(buttons=buttons, direction='down', xanchor='left', yanchor='top')]

    return data, annotations_idx, updatemenus









def choose_conformers_input():
    string=input('Enter the conformers numbers: ')
    conformer_numbers=string.split(' ')
    conformer_numbers=[int(i) for i in conformer_numbers]
    return conformer_numbers


def unite_buttons(buttons_list, ref_index=0):
    buttons_keys=buttons_list[0].keys()
    united_button=dict.fromkeys(buttons_keys)
    for key in buttons_keys:
        if key=='args':
            all_annotations=[buttons[key][0]['scene.annotations'] for buttons in buttons_list]
            united_annotations=list(zip(*all_annotations))
            united_button[key]=[{'scene.annotations': united_annotations}]
        else:
            united_button[key]=buttons_list[ref_index][key]
    return united_button

def unite_updatemenus(updatemenus_list, ref_index=0):
    menus_keys=updatemenus_list[ref_index][0].keys()
    united_updatemenus_list=dict.fromkeys(menus_keys)
    for key in menus_keys:
        if key=='buttons':
            buttons_list=[updatemenus[0].get(key) for updatemenus in updatemenus_list]
            buttons_num=len(buttons_list[0])   
            segregated_buttons=[]
            for i in range(buttons_num):
                type_buttons=[buttons[i] for buttons in buttons_list]
                segregated_buttons.append(type_buttons)
            buttons=[unite_buttons(buttons) for buttons in segregated_buttons]
            united_updatemenus_list[key]=buttons
        else:
            united_updatemenus_list[key]=updatemenus_list[ref_index][0][key]
    return [united_updatemenus_list]


def compare_molecules(coordinates_df_list: List[pd.DataFrame],conformer_numbers:List[int]=None):
    if conformer_numbers is None:
        conformer_numbers=choose_conformers_input()
    # Create a subplot with 3D scatter plot
    colors_list=['red','purple','blue','green','yellow','orange','brown','black','pink','cyan','magenta']
    new_coodinates_df_list=[coordinates_df_list[i] for i in conformer_numbers]
    new_coodinates_df_list=renumbering_df_list(new_coodinates_df_list)
    # coordinates_df_list=renumber_xyz_by_mcs(coordinates_df_list)  ##needs fixing , renumbering not working.
    xyz_df=(coordinates_df_list[conformer_numbers[0]])
    data_main, annotations_id_main, updatemenus = plot_interactions(xyz_df,'grey')
    updatemenus_list=[updatemenus]
    # Iterate through the conformer numbers and create traces for each conformer
    for  conformer_number,color in zip((conformer_numbers[1:]),colors_list):
        xyz_df = coordinates_df_list[conformer_number]
        data, annotations_id, updatemenus_main = plot_interactions(xyz_df,color)
        data_main += data
        annotations_id_main += annotations_id
        updatemenus_list.append(updatemenus_main)
    # Set axis parameters
    updatemenus = unite_updatemenus(updatemenus_list)

    coords = xyz_df[['x','y','z']].to_numpy(float)
    span = np.ptp(coords, axis=0)            # total extent along each axis
    center = coords.mean(axis=0)
    pad = max(span) * 0.8                    # extra space around molecule

    x_range = [center[0] - span[0]/2 - pad, center[0] + span[0]/2 + pad]
    y_range = [center[1] - span[1]/2 - pad, center[1] + span[1]/2 + pad]
    z_range = [center[2] - span[2]/2 - pad, center[2] + span[2]/2 + pad]

    axis_params = dict(
        showgrid=True,
        showbackground=False,
        showticklabels=True,
        zeroline=True,
        titlefont=dict(color='black', size=12),
    )

    # --- Layout with fixed large frame ---
    layout = dict(
        scene=dict(
            aspectmode='data',
            xaxis=dict(**axis_params, title="X Axis", range=x_range),
            yaxis=dict(**axis_params, title="Y Axis", range=y_range),
            zaxis=dict(**axis_params, title="Z Axis", range=z_range),
            annotations=annotations_id_main,
        ),
        scene_camera=dict(
            eye=dict(x=1.8, y=1.8, z=1.8),
            center=dict(x=0, y=0, z=0),
            projection=dict(type='orthographic')
        ),
        margin=dict(r=0, l=0, b=0, t=30),
        showlegend=True,
        updatemenus=updatemenus,
        uirevision="keep-camera"
    )

    fig = go.Figure(data=data_main, layout=layout)
    fig.show()
    return fig

import dash
from dash import html, dcc, Output, Input, State
import plotly.graph_objects as go
import pandas as pd

def show_single_molecule(molecule_name,xyz_df=None,dipole_df=None, origin=None,sterimol_params=None,color='black'):
    if xyz_df is None:
        xyz_df=hf.get_df_from_file(hf.choose_filename()[0])
    # Create a subplot with 3D scatter plot
    data_main, annotations_id_main, updatemenus = plot_interactions(xyz_df,color,dipole_df=dipole_df, origin=origin,sterimol_params=sterimol_params)

    axis_params = dict(
        showgrid=False,
        showbackground=False,
        showticklabels=False,
        zeroline=False
    )

    scene = dict(
        xaxis=dict(**axis_params, title=dict(text='X', font=dict(color='white'))),
        yaxis=dict(**axis_params, title=dict(text='Y', font=dict(color='white'))),
        zaxis=dict(**axis_params, title=dict(text='Z', font=dict(color='white'))),
        annotations=annotations_id_main
    )

    layout = dict(
        title=dict(
            text=molecule_name,
            x=0.5,
            y=0.9,
            xanchor='center',
            yanchor='top'
        ),
        scene=scene,
        margin=dict(r=0, l=0, b=0, t=0),
        showlegend=False,
        updatemenus=updatemenus
        )
    
    fig = go.Figure(data=data_main, layout=layout)
    html=fig.show()
    run_app(fig)

    
    return html


def run_app(figure):
        # Create a Dash app
    app = dash.Dash(__name__)

    # App layout
    app.layout = html.Div([
        dcc.Graph(id='molecule-plot', figure=figure), # Replace "Water" with your molecule
        html.Div(id='clicked-data', children=[]),
        html.Button('Save Clicked Atom', id='save-button', n_clicks=0),
        html.Div(id='saved-atoms', children=[])
    ])

    # Callback to display clicked data
    @app.callback(
        Output('clicked-data', 'children'),
        Input('molecule-plot', 'clickData'),
        prevent_initial_call=True
    )
    def display_click_data(clickData):
        if clickData:
            return f"Clicked Point: {clickData['points'][0]['pointIndex']}"
        return "Click on an atom."

    # Callback to save clicked atom index
    @app.callback(
        Output('saved-atoms', 'children'),
        Input('save-button', 'n_clicks'),
        State('molecule-plot', 'clickData'),
        State('saved-atoms', 'children'),
        prevent_initial_call=True
    )
    def save_clicked_atom(n_clicks, clickData, saved_atoms):
        if clickData:
            saved_atoms.append(clickData['points'][0]['pointIndex'])
        return f"Saved Atom Indices: {saved_atoms}"
        
import matplotlib.pyplot as plt





# Atom color map (CPK + common metals)
atom_colors = {
    'C': '#404040', 'H': '#AAAAAA', 'O': '#E83030', 'N': '#3060FF',
    'S': '#D4A000', 'Cl': '#44BB44', 'F': '#55DD55', 'Br': '#994400',
    'I': '#773377', 'P': '#FF8800', 'B': '#FFB5B5', 'Si': '#F0C8A0',
    # transition metals
    'Sc': '#E6E6E6', 'Ti': '#BFC2C7', 'V': '#A6A6AB', 'Cr': '#8A99C7',
    'Mn': '#9C7AC7', 'Fe': '#E06633', 'Co': '#F090A0', 'Ni': '#50D050',
    'Cu': '#C88033', 'Zn': '#7D80B0',
    'Y': '#94FFFF', 'Zr': '#94E0E0', 'Nb': '#73C2C9', 'Mo': '#54B5B5',
    'Tc': '#3B9E9E', 'Ru': '#248F8F', 'Rh': '#0A7D8C', 'Pd': '#006985',
    'Ag': '#C0C0C0', 'Cd': '#FFD98F',
    'Hf': '#4DC2FF', 'Ta': '#4DA6FF', 'W': '#2194D6', 'Re': '#267DAB',
    'Os': '#266696', 'Ir': '#175487', 'Pt': '#D0D0E0', 'Au': '#FFD123',
    'Hg': '#B8B8D0',
    'Li': '#CC80FF', 'Na': '#AB5CF2', 'K': '#8F40D4',
    'Mg': '#8AFF00', 'Ca': '#3DFF00',
    'La': '#70D4FF', 'Ce': '#FFFFC7', 'Lu': '#00AB24',
}
_B1_COLOR = '#2E8B57'
_B5_COLOR = '#CC2222'
_L_COLOR  = '#2266CC'


from matplotlib.patches import ConnectionPatch, FancyArrowPatch
import matplotlib.patches as mpatches

def plot_b1_visualization(rotated_plane, edited_coordinates_df,
                          sterimol_df=None, n_points=100,
                          title="XZ plane — End-on view"):
    """
    End-on view of the XZ plane showing atom VdW circles, B1 and B5 arrows,
    and (if sterimol_df supplied) the 3-D vector-to-plane angle.
    Returns the matplotlib Figure.
    """
    # --- Extremes in the XZ plane ---
    max_x, min_x = rotated_plane[:, 0].max(), rotated_plane[:, 0].min()
    max_z, min_z = rotated_plane[:, 1].max(), rotated_plane[:, 1].min()
    avs = np.abs([max_x, min_x, max_z, min_z])
    b1_val   = float(avs.min())
    b1_index = int(avs.argmin())
    b1_tip   = [(max_x, 0), (min_x, 0), (0, max_z), (0, min_z)][b1_index]

    # B5: atom with largest XZ radius (from edited_coordinates_df)
    b5_row  = edited_coordinates_df.loc[edited_coordinates_df['B5'].idxmax()]
    b5_x    = float(b5_row['x'])
    b5_z    = float(b5_row['z'])
    b5_val  = float(b5_row['B5'])

    # angle label
    angle_label = None
    if sterimol_df is not None and 'B1_B5_angle' in sterimol_df.columns:
        angle_label = float(sterimol_df['B1_B5_angle'].iloc[0])

    # --- Figure ---
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=13, fontweight='bold', pad=10)
    ax.set_xlabel("X (Å)", fontsize=11)
    ax.set_ylabel("Z (Å)", fontsize=11)
    ax.axhline(0, color='k', lw=0.4, alpha=0.3)
    ax.axvline(0, color='k', lw=0.4, alpha=0.3)

    # --- Atom circles ---
    n_circles = rotated_plane.shape[0] // n_points
    atom_list = list(edited_coordinates_df.iterrows())

    for i in range(n_circles):
        pts = rotated_plane[i*n_points:(i+1)*n_points, :]
        closed = np.vstack([pts, pts[0]])
        atom_idx, row = atom_list[i]
        color = atom_colors.get(row['atom'], '#808080')
        ax.plot(closed[:, 0], closed[:, 1],
                color=color, lw=1.0, alpha=0.85, solid_joinstyle='round')

        cx, cz = pts.mean(axis=0)
        vec = np.array([cx, cz])
        vn  = np.linalg.norm(vec)
        offset = 0.28 * vec / vn if vn > 1e-5 else np.array([0.28, 0.0])
        lx, lz = vec + offset
        ax.text(lx, lz, str(atom_idx),
                ha='center', va='center', fontsize=7, fontweight='bold',
                color='black',
                bbox=dict(boxstyle='circle,pad=0.15', fc='white', ec='none', alpha=0.6))
        ax.add_patch(ConnectionPatch(
            xyA=(cx, cz), xyB=(lx, lz),
            coordsA='data', coordsB='data',
            arrowstyle='-', lw=0.35, color=color, alpha=0.5))

    # --- B1 tangent wall (dashed) + arrow ---
    lim = max(abs(max_x), abs(min_x), abs(max_z), abs(min_z)) * 1.15
    if b1_index in (0, 1):   # vertical wall at x = ±b1_val
        wall_x = b1_tip[0]
        ax.axvline(wall_x, color=_B1_COLOR, ls='--', lw=1.0, alpha=0.55)
    else:                     # horizontal wall at z = ±b1_val
        wall_z = b1_tip[1]
        ax.axhline(wall_z, color=_B1_COLOR, ls='--', lw=1.0, alpha=0.55)

    ax.annotate('', xy=b1_tip, xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color=_B1_COLOR, lw=2.0))
    ax.text(b1_tip[0]*0.52, b1_tip[1]*0.52,
            f"B1\n{b1_val:.2f} Å",
            fontsize=9, ha='center', va='center',
            fontweight='bold', color=_B1_COLOR,
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=_B1_COLOR, lw=0.8, alpha=0.8))

    # --- B5 arrow ---
    ax.annotate('', xy=(b5_x, b5_z), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color=_B5_COLOR, lw=2.0))
    ax.text(b5_x*0.60, b5_z*0.60,
            f"B5\n{b5_val:.2f} Å",
            fontsize=9, ha='center', va='center',
            fontweight='bold', color=_B5_COLOR,
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=_B5_COLOR, lw=0.8, alpha=0.8))

    # --- Angle arc between B1 and B5 vectors ---
    a1 = np.arctan2(b1_tip[1], b1_tip[0])
    a5 = np.arctan2(b5_z, b5_x)
    a_lo, a_hi = sorted([a1, a5])
    if a_hi - a_lo > np.pi:
        a_lo, a_hi = a_hi, a_lo + 2*np.pi
    arc = np.linspace(a_lo, a_hi, 80)
    r_arc = min(b1_val, b5_val) * 0.55
    ax.plot(r_arc*np.cos(arc), r_arc*np.sin(arc), color='gray', lw=1.0, alpha=0.6)
    mid = (a_lo + a_hi) / 2
    disp_angle = angle_label if angle_label is not None else float(np.degrees(abs(a5 - a1)))
    ax.text(r_arc*1.25*np.cos(mid), r_arc*1.25*np.sin(mid),
            f"{disp_angle:.1f}°",
            fontsize=8, ha='center', va='center', color='gray', fontweight='bold')

    # --- Angle annotation box (upper-left corner) ---
    angle_src = "3-D" if angle_label is not None else "2-D proj."
    ax.text(0.03, 0.97,
            f"∠B1–B5 = {disp_angle:.1f}°  ({angle_src})",
            transform=ax.transAxes, ha='left', va='top',
            fontsize=10, fontweight='bold', color='#444444',
            bbox=dict(boxstyle='round,pad=0.35', fc='#F7F7F7', ec='#AAAAAA', lw=1))

    ax.grid(alpha=0.12)
    plt.tight_layout()
    return fig



def generate_circle(center_x, center_y, radius, n_points=20):
    """
    Generate circle coordinates given a center and radius.
    Returns a DataFrame with columns 'x' and 'y'.
    """
    theta = np.linspace(0, 2 * np.pi, n_points)
    x = center_x + radius * np.cos(theta)
    y = center_y + radius * np.sin(theta)
    return np.column_stack((x, y))


def plot_L_B5_plane(edited_coordinates_df, sterimol_df, n_points=100,
                    title="YZ plane — Side view"):
    """
    Side view (YZ plane) showing L, B5, and B1 reference.
    B5 atom is identified by the 'B5' column maximum (same atom as calc_sterimol).
    Returns the matplotlib Figure.
    """
    # --- Sterimol values from df ---
    L_val  = float(sterimol_df['L'].iloc[0])
    B5_val = float(sterimol_df['B5'].iloc[0])
    B1_val = float(sterimol_df['B1'].iloc[0])
    loc_B5 = float(sterimol_df['loc_B5'].iloc[0])
    angle  = float(sterimol_df['B1_B5_angle'].iloc[0]) if 'B1_B5_angle' in sterimol_df.columns else None

    # --- B5 atom (same selection as calc_sterimol) ---
    b5_idx = edited_coordinates_df['B5'].idxmax()
    b5_row = edited_coordinates_df.loc[b5_idx]
    b5_y, b5_z = float(b5_row['y']), float(b5_row['z'])

    # --- Build YZ circles ---
    circles, centers = [], []
    seen_types = {}
    for atom_idx, r in edited_coordinates_df.iterrows():
        pts = generate_circle(r['y'], r['z'], r['radius'], n_points=n_points)
        circles.append(pts)
        centers.append((atom_idx, r['atom'], float(r['y']), float(r['z'])))
        seen_types[r['atom']] = atom_colors.get(r['atom'], '#808080')

    # --- Figure ---
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=13, fontweight='bold', pad=10)
    ax.set_xlabel("Y (Å) — substituent axis", fontsize=11)
    ax.set_ylabel("Z (Å)", fontsize=11)
    ax.axhline(0, color='k', lw=0.4, alpha=0.3)
    ax.axvline(0, color='k', lw=0.4, alpha=0.3)

    # --- Atom circles + index labels ---
    for i, pts in enumerate(circles):
        closed = np.vstack([pts, pts[0]])
        atom_idx, atom_type, y0, z0 = centers[i]
        color = seen_types[atom_type]
        ax.plot(closed[:, 0], closed[:, 1], color=color, lw=1.0, alpha=0.85)
        ax.text(y0, z0, str(atom_idx),
                ha='center', va='center', fontsize=7, color='black',
                bbox=dict(boxstyle='circle,pad=0.1', fc='white', ec='none', alpha=0.5))

    # --- L arrow (Y axis, length = L_val) ---
    ax.annotate('', xy=(L_val, 0), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color=_L_COLOR, lw=2.2))
    ax.text(L_val * 0.55, 0.12,
            f"L = {L_val:.2f} Å",
            fontsize=9, color=_L_COLOR, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=_L_COLOR, lw=0.8, alpha=0.8))

    # --- loc_B5 dashed vertical line ---
    ax.axvline(loc_B5, color=_B5_COLOR, ls=':', lw=1.0, alpha=0.55)
    ax.text(loc_B5, ax.get_ylim()[0] if ax.get_ylim()[0] != 0 else -0.3,
            f"loc_B5\n{loc_B5:.2f}", fontsize=7, color=_B5_COLOR,
            ha='center', va='top', alpha=0.7)

    # --- B5 arrow to actual B5 atom ---
    ax.annotate('', xy=(b5_y, b5_z), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color=_B5_COLOR, lw=2.2))
    ax.text(b5_y*0.58, b5_z*0.58,
            f"B5 = {B5_val:.2f} Å",
            fontsize=9, color=_B5_COLOR, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=_B5_COLOR, lw=0.8, alpha=0.8))

    # --- B1 reference: double-headed arrow ±B1/2 in Z at Y=0 ---
    ax.annotate('', xy=(0,  B1_val), xytext=(0, -B1_val),
                arrowprops=dict(arrowstyle='<->', color=_B1_COLOR, lw=1.8))
    ax.text(0.08, 0,
            f"B1 = {B1_val:.2f} Å\n(XZ property)",
            fontsize=8, color=_B1_COLOR, fontweight='bold', va='center',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=_B1_COLOR, lw=0.8, alpha=0.8))

    # --- Atom type legend ---
    legend_handles = [mpatches.Patch(color=c, label=a) for a, c in seen_types.items()]
    ax.legend(handles=legend_handles, title="Atom types",
              fontsize=8, title_fontsize=8,
              loc='upper left', bbox_to_anchor=(1.02, 1.0),
              framealpha=0.9)

    # --- Summary text box ---
    summary_lines = [
        f"B1 = {B1_val:.3f} Å",
        f"B5 = {B5_val:.3f} Å",
        f"L  = {L_val:.3f} Å",
    ]
    if angle is not None:
        summary_lines.append(f"∠B1–B5 = {angle:.1f}° (3-D)")
    ax.text(1.02, 0.48, "\n".join(summary_lines),
            transform=ax.transAxes, ha='left', va='top',
            fontsize=9, family='monospace',
            bbox=dict(boxstyle='round,pad=0.4', fc='#F5F5F5', ec='#AAAAAA', lw=1))

    ax.grid(alpha=0.12)
    plt.tight_layout()
    return fig


if __name__ == '__main__':
    
    pass
