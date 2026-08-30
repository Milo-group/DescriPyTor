"""
DescriPytor descriptor-extraction web app (Streamlit)
=====================================================

A browser front-end over the toolkit so lab mates can, without touching the
command line:

  1. provide a molecule set (upload .feather/.xyz, or point at a server folder),
  2. pick atoms visually in the embedded 3D picker and download a run_config.json,
  3. upload that config (or edit the template), choose engines,
  4. run the extraction and download the merged feature CSV.

Run it:
    pip install -r requirements.txt
    streamlit run app.py

Then share the URL it prints (use --server.address 0.0.0.0 to expose on the LAN).
See README_webapp.md for deployment + the dependency reality (DescriPytor, xTB...).
"""

from __future__ import annotations

import contextlib
import io
import json
import os
import sys
import tempfile
import threading
import time
from pathlib import Path

import pandas as pd
import streamlit as st

# --------------------------------------------------------------------------- #
# Make the toolkit importable (this file lives in .../descriptor_extraction_toolkit/webapp)
# --------------------------------------------------------------------------- #
TOOLKIT_DIR = Path(__file__).resolve().parent.parent
if str(TOOLKIT_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLKIT_DIR))

DEFAULT_ROOT = os.environ.get(
    "DESCRIPYTOR_ROOT",
    str(Path(__file__).resolve().parents[3]),
)

OPTIONAL_PACKAGES = [
    ("rdkit", "rdkit", "2D/3D structures, fingerprints, SMILES perception"),
    ("python-igraph", "igraph", "2D bond-graph view in the atom picker"),
    ("morfeus-ml", "morfeus", "Morfeus suite + Sterimol reference engine"),
    ("mordred", "mordred", "Mordred 2D/3D descriptors"),
    ("deepchem", "deepchem", "DeepChem descriptors"),
    ("ase", "ase", "xyz parsing for SMILES perception"),
    ("bassa-reg", "bassa_reg", "BASSA Bayesian spike-and-slab modeling"),
    ("autoqchem", "autoqchem", "qm engine (Gaussian log parsing)"),
    ("scikit-learn", "sklearn", "PCA projection, M3 modeling"),
    ("umap-learn", "umap", "UMAP projection (Results tab)"),
    ("plotly", "plotly", "interactive PCA/UMAP scatter (Results tab)"),
]


def _check_optional_deps():
    """Availability of every optional package a webapp feature depends on."""
    import importlib.util
    rows = []
    for pip_name, import_name, purpose in OPTIONAL_PACKAGES:
        try:
            available = importlib.util.find_spec(import_name) is not None
        except Exception:
            available = False
        rows.append({"package": pip_name, "available": available, "purpose": purpose})
    return rows


st.set_page_config(page_title="DescriPytor Descriptor Extractor", layout="wide")
ss = st.session_state
ss.setdefault("work_dir", tempfile.mkdtemp(prefix="descripytor_web_"))
ss.setdefault("feather_dir", "")
ss.setdefault("xyz_dir", "")
ss.setdefault("result_df", None)
ss.setdefault("run_logs", "")
ss.setdefault("last_config", None)
ss.setdefault("last_output_csv", "")
ss.setdefault("model_df", None)
ss.setdefault("model_df_source", "")
ss.setdefault("m3_results", None)
ss.setdefault("m3_model", None)
ss.setdefault("m3_report", None)
ss.setdefault("bassa_run_dir", "")
ss.setdefault("dataset_preview", None)
ss.setdefault("fp_mol_b64", None)
ss.setdefault("fp_on_bits", None)
ss.setdefault("fp_bit_info", None)
ss.setdefault("fp_mol_obj", None)
ss.setdefault("scope_page_png", None)
ss.setdefault("extraction_thread", None)
ss.setdefault("extraction_buf", None)
ss.setdefault("extraction_outcome", None)
ss.setdefault("extraction_cancel_event", None)
ss.setdefault("extraction_cfg", None)
ss.setdefault("export_zip_bytes", None)
ss.setdefault("proj_fig", None)
ss.setdefault("proj_info", None)
ss.setdefault("proj_scree_fig", None)
ss.setdefault("proj_loadings_fig", None)


# --------------------------------------------------------------------------- #
# Sidebar — environment
# --------------------------------------------------------------------------- #
st.sidebar.header("Environment")
root_dir = st.sidebar.text_input(
    "DescriPyTor root", value=DEFAULT_ROOT,
    help="Folder that contains M2_data_extractor etc. Must exist on THIS server.")
os.environ["DESCRIPYTOR_ROOT"] = root_dir
st.sidebar.caption(f"Working dir: `{ss['work_dir']}`")
st.sidebar.caption("Engines whose Python deps (rdkit, morfeus-ml, xTB, aqme…) "
                   "aren't installed here will skip themselves with a note.")

with st.sidebar.expander("Environment check (optional packages)"):
    for row in _check_optional_deps():
        icon = "✅" if row["available"] else "⚠️"
        st.caption(f"{icon} **{row['package']}** — {row['purpose']}")

if st.sidebar.button("Reset everything", help="Clear loaded data, run results, and model results "
                     "(does not touch files already downloaded/saved on disk)."):
    _keep = {"work_dir"}
    for _k in list(ss.keys()):
        if _k not in _keep:
            del ss[_k]
    st.rerun()

st.title("DescriPytor — Descriptor Extraction")
tab_data, tab_pick, tab_run, tab_res, tab_model = st.tabs(
    ["1 · Data", "2 · Pick atoms", "3 · Configure & run", "4 · Results", "5 · Modeling"])


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _dir_signature(dirpath: str) -> str:
    """
    Cheap signature of a directory's contents (file count + latest mtime +
    total size), used to bust the feather-load cache below when the folder's
    contents change without the path itself changing (e.g. re-uploading).
    """
    try:
        entries = [f for f in Path(dirpath).iterdir() if f.is_file()]
        if not entries:
            return "0:0:0"
        return (f"{len(entries)}:{max(f.stat().st_mtime for f in entries):.1f}:"
               f"{sum(f.stat().st_size for f in entries)}")
    except Exception:
        return "0:0:0"


@st.cache_resource(show_spinner="Loading feather set…")
def _load_feather_molecules_cached(feather_dir: str, root_dir: str, sig: str):
    """
    Load a DescriPytor feather set once and keep it cached across reruns.

    Every widget interaction anywhere in the app reruns this whole script, so
    without caching, just rendering the fingerprint-viewer's molecule dropdown
    would reparse every feather file on disk each time. `sig` (a cheap
    directory signature) is part of the cache key so re-uploading different
    files to the same folder still busts the cache.
    """
    import descriptor_extractor as dx
    dx._add_descripytor_to_path(root_dir)
    from data_extractor import Molecules
    return Molecules(feather_dir)


def _save_uploads(files, subdir):
    """Save uploaded files into work_dir/subdir and return that dir (or '')."""
    if not files:
        return ""
    dest = Path(ss["work_dir"]) / subdir
    dest.mkdir(parents=True, exist_ok=True)
    for f in files:
        (dest / f.name).write_bytes(f.getbuffer())
    return str(dest)


def _first_feather_molecule_xyz(feather_dir):
    """Return (xyz_text, mol_payload) for the first molecule, or (None, None)."""
    try:
        import make_picker as mp
        mols = _load_feather_molecules_cached(feather_dir, root_dir, _dir_signature(feather_dir))
        if not getattr(mols, "molecules", None):
            st.info("No molecules loaded from that feather folder.")
            return None, None
        mol = mols.molecules[0]
        return mp.molecule_to_xyz(mol), mp.molecule_payload(mol)
    except Exception as e:  # noqa
        st.info(f"Could not load a feather molecule for the picker ({e}). "
                "You can still load an .xyz in the picker.")
        return None, None


def _build_picker_html(xyz_text=None, mol_data=None):
    """Inject (optional) molecule into atom_picker.html and return the HTML string."""
    import make_picker as mp
    template = (TOOLKIT_DIR / "atom_picker.html").read_text(encoding="utf-8")
    html = template.replace("__XYZ_DATA__", mp._js_escape(xyz_text) if xyz_text else "__XYZ_DATA__")
    html = html.replace("__MOL_DATA__",
                        mp._js_escape(json.dumps(mol_data)) if mol_data else
                        mp._js_escape(json.dumps({"dipole": None, "modes": []})))
    if xyz_text:
        graph2d_b64 = mp.graph2d_png_base64(xyz_text)
        html = html.replace("__GRAPH2D_IMG__", graph2d_b64)
        html = html.replace("__GRAPH2D_DISPLAY__", "" if graph2d_b64 else "none")
    # paths are filled in the Configure tab instead, leave placeholders -> JS defaults
    return html


def _template_config():
    """Return the toolkit's full starter config, or a small fallback."""
    try:
        from descriptor_extractor import build_template_config
        return build_template_config()
    except Exception:
        return {
            "root_dir": root_dir,
            "feather_dir": "",
            "xyz_dir": "",
            "derive_xyz_from_feathers": False,
            "output_csv": "merged_features.csv",
            "engines": {},
        }


def _overlay_run_paths(cfg, derive_xyz=False):
    """Overlay the server-side paths selected in the app onto a run config."""
    out = dict(cfg)
    out["root_dir"] = root_dir
    if ss["feather_dir"]:
        out["feather_dir"] = ss["feather_dir"]
    if ss["xyz_dir"]:
        out["xyz_dir"] = ss["xyz_dir"]
    out["derive_xyz_from_feathers"] = bool(derive_xyz)
    out["output_csv"] = str(Path(ss["work_dir"]) / "merged_features.csv")
    return out


def _json_area(label, value, key, height=220):
    """Text-area JSON editor with a compact error message."""
    txt = st.text_area(label, value=json.dumps(value, indent=2), height=height, key=key)
    try:
        return json.loads(txt), None
    except Exception as e:  # noqa
        return None, e


def _parse_json_or_pyliteral(text):
    """
    Parse text as JSON first, falling back to a Python literal.

    The picker's "Python answers_dict" export tab produces Python syntax
    (True/False/None, single quotes) rather than strict JSON, so a straight
    paste of that would fail json.loads - this accepts either.
    """
    try:
        return json.loads(text), None
    except Exception:
        pass
    try:
        import ast
        return ast.literal_eval(text), None
    except Exception as e:  # noqa
        return None, e


def _enabled_engines(cfg):
    """Return the enabled engine names in a config."""
    return [
        name for name, engine_cfg in cfg.get("engines", {}).items()
        if engine_cfg.get("enabled")
    ]


def _run_readiness(cfg):
    """Return user-facing blockers before extraction can run."""
    blockers = []
    if not ss["feather_dir"] and not ss["xyz_dir"]:
        blockers.append("Load `.feather` and/or `.xyz` data in tab 1.")
    if not _enabled_engines(cfg):
        blockers.append("Enable at least one descriptor engine.")
    return blockers


# Scalar (non atom-index) fields in descripytor_full's "atoms" dict - everything
# else in that dict is atom indices/pairs/triplets.
_THRESHOLD_LIKE_KEYS = {"stretch_threshold", "bend_threshold"}


def _collect_atom_indices(obj):
    """Recursively collect every int found in a nested list/tuple (atom-pair configs)."""
    out = set()
    if isinstance(obj, bool):
        return out
    if isinstance(obj, int):
        out.add(obj)
    elif isinstance(obj, (list, tuple)):
        for item in obj:
            out |= _collect_atom_indices(item)
    return out


def _atom_indices_from_cfg(cfg):
    """Pull every referenced atom index out of all enabled engine configs."""
    indices = set()
    for name, ecfg in cfg.get("engines", {}).items():
        if not ecfg.get("enabled"):
            continue
        if name == "descripytor_full":
            for key, val in (ecfg.get("atoms") or {}).items():
                if key in _THRESHOLD_LIKE_KEYS:
                    continue
                indices |= _collect_atom_indices(val)
        for key in ("pairs", "angles", "bond_lengths", "sterimol_pairs", "cone_atoms", "pyramid_atoms"):
            if key in ecfg:
                indices |= _collect_atom_indices(ecfg[key])
        metal_index = ecfg.get("metal_index")
        if isinstance(metal_index, int):
            indices.add(metal_index)
        elif isinstance(metal_index, list):
            indices |= _collect_atom_indices(metal_index)
    return indices


def _check_atom_indices_against_feather(cfg, feather_dir):
    """
    Return [(molecule_name, atom_count, max_referenced_index), ...] for every
    molecule whose atom count is smaller than the largest atom index referenced
    anywhere in the enabled engine configs. Best-effort static check - it only
    catches "index bigger than this molecule has", not other misconfigurations.
    """
    indices = _atom_indices_from_cfg(cfg)
    if not indices or not feather_dir:
        return []
    max_idx = max(indices)
    problems = []
    _, by_name = _molecule_names_and_objs(feather_dir)
    for name, mol in by_name.items():
        n_atoms = len(getattr(mol, "xyz_df", []))
        if n_atoms and max_idx > n_atoms:
            problems.append((name, n_atoms, max_idx))
    return problems


def _dataset_preview(feather_dir, xyz_dir):
    """Best-effort preview of what's actually loaded: molecule counts/names (not just paths)."""
    info = {}
    if feather_dir:
        try:
            mols = _load_feather_molecules_cached(feather_dir, root_dir, _dir_signature(feather_dir))
            names = list(getattr(mols, "molecules_names", None) or [
                getattr(m, "molecule_name", "?") for m in getattr(mols, "molecules", [])
            ])
            entry = {"count": len(names), "names": names}
            failed = getattr(mols, "failed_molecules", None)
            if failed:
                entry["failed"] = list(failed)
            info["feather"] = entry
        except Exception as e:  # noqa
            info["feather"] = {"error": str(e)}
    if xyz_dir:
        try:
            xyz_files = sorted(Path(xyz_dir).glob("*.xyz"))
            info["xyz"] = {"count": len(xyz_files), "names": [p.stem for p in xyz_files]}
        except Exception as e:  # noqa
            info["xyz"] = {"error": str(e)}
    return info


# --------------------------------------------------------------------------- #
# Modeling (tab 4) helpers
# --------------------------------------------------------------------------- #
def _parse_pasted_outputs(text):
    """
    Parse the "Output values to merge" textarea into a list of (name_or_None, value).

    Accepts either:
      - one number per line (or comma-separated on one line) -> matched by row order
      - `name,value` pairs, one per line -> matched by name
    """
    lines = [ln.strip() for ln in text.replace(",", "\n").splitlines() if ln.strip()]
    pairs = []
    for ln in lines:
        parts = [p.strip() for p in ln.split(",") if p.strip() != ""] if "," in ln else [ln]
        if len(parts) >= 2:
            name, val = parts[0], parts[1]
            try:
                pairs.append((name, float(val)))
            except ValueError:
                continue
        else:
            try:
                pairs.append((None, float(parts[0])))
            except ValueError:
                continue
    return pairs


def _merge_outputs_into_df(df, pairs, new_col_name, name_col=None):
    """Return a copy of df with new_col_name populated from pairs (name,value) or (None,value)."""
    out = df.copy()
    if not pairs:
        return out
    by_name = all(p[0] is not None for p in pairs)
    if by_name:
        name_col = name_col or next(
            (c for c in out.columns if c.lower() in ("name", "names", "molecule", "molecule_name", "mol", "compound")),
            out.columns[0],
        )
        value_map = {str(n): v for n, v in pairs}
        out[new_col_name] = out[name_col].astype(str).map(value_map)
    else:
        values = [v for _, v in pairs]
        out[new_col_name] = pd.Series(values[: len(out)]).reindex(range(len(out))).values
    return out


def _add_m3_to_path():
    """Make `modeling` (M3_modeler) importable, same convention as descriptor_extractor."""
    try:
        import descriptor_extractor as dx
        dx._add_descripytor_to_path(root_dir)
    except Exception as e:  # noqa
        st.error(f"Could not add M3_modeler to sys.path (check DescriPyTor root): {e}")
        raise


def _build_export_zip() -> bytes:
    """
    Bundle everything generated so far in this session - run_config.json, the
    merged features CSV, the run log, and any M3/BASSA modeling outputs - into
    one zip, so users aren't stuck downloading each piece separately.
    """
    import zipfile

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        cfg = ss.get("last_config")
        if cfg:
            zf.writestr("run_config.json", json.dumps(cfg, indent=2))

        df = ss.get("result_df")
        csv_path = ss.get("last_output_csv")
        if csv_path and Path(csv_path).exists():
            zf.write(csv_path, "merged_features.csv")
        elif df is not None and len(df):
            zf.writestr("merged_features.csv", df.to_csv())

        if ss.get("run_logs"):
            zf.writestr("run_log.txt", ss["run_logs"])

        m3_results = ss.get("m3_results")
        if m3_results is not None and len(m3_results):
            zf.writestr("modeling/m3_results.csv", m3_results.to_csv(index=False))

        report = ss.get("m3_report") or {}
        pdf_path = report.get("pdf_path")
        if pdf_path and Path(pdf_path).exists():
            zf.write(pdf_path, f"modeling/{Path(pdf_path).name}")
        for p in report.get("png_files") or []:
            if Path(p).exists():
                zf.write(p, f"modeling/report_images/{Path(p).name}")

        bassa_dir = ss.get("bassa_run_dir")
        if bassa_dir and Path(bassa_dir).exists():
            for f in Path(bassa_dir).iterdir():
                if f.is_file():
                    zf.write(f, f"modeling/bassa/{f.name}")

        scope_png = ss.get("scope_page_png")
        if scope_png:
            zf.writestr("scope_page.png", scope_png)

    buf.seek(0)
    return buf.getvalue()


def _build_projection_figure(df, feature_cols, method="PCA", color_col=None, n_neighbors=15):
    """
    Reduce the chosen numeric columns to 2D (PCA or UMAP) and return a plotly
    scatter figure plus a short info caption.

    For PCA, the caption (and the axis labels) report the % of total variance
    each of the two components explains - UMAP is nonlinear and has no
    equivalent variance-explained figure, so its axes are just labeled UMAP-1/2.
    For PCA only, also returns a scree plot (variance per component, up to 10)
    and a loadings chart (top features driving PC1/PC2).

    Returns (fig, info, error, scree_fig, loadings_fig) - exactly one of
    (fig, error) is None; scree_fig/loadings_fig are None for UMAP.
    """
    from sklearn.impute import SimpleImputer
    from sklearn.preprocessing import StandardScaler

    sub = df[feature_cols].apply(pd.to_numeric, errors="coerce")
    # Constant or all-NaN columns can't be scaled meaningfully - drop them.
    keep_cols = [c for c in sub.columns if sub[c].notna().any() and sub[c].nunique(dropna=True) > 1]
    dropped = [c for c in feature_cols if c not in keep_cols]
    sub = sub[keep_cols]

    if sub.shape[1] < 2:
        return None, None, "Need at least 2 usable numeric columns (non-constant, not all-NaN).", None, None
    if len(sub) < 3:
        return None, None, "Need at least 3 molecules to compute a 2D projection.", None, None

    x_imputed = SimpleImputer(strategy="mean").fit_transform(sub.to_numpy())
    x_scaled = StandardScaler().fit_transform(x_imputed)

    scree_fig = None
    loadings_fig = None

    if method == "PCA":
        from sklearn.decomposition import PCA
        if min(x_scaled.shape[0] - 1, x_scaled.shape[1]) < 2:
            return None, None, "Not enough columns/rows left for a 2D PCA projection.", None, None
        pca = PCA(n_components=2, random_state=42)
        coords = pca.fit_transform(x_scaled)
        var_pct = pca.explained_variance_ratio_ * 100
        x_label = f"PC1 ({var_pct[0]:.1f}% variance)"
        y_label = f"PC2 ({var_pct[1]:.1f}% variance)"
        info = (f"PC1 + PC2 explain {var_pct.sum():.1f}% of total variance "
                f"across {len(keep_cols)} feature(s).")

        try:
            import plotly.express as px

            # Scree plot: refit with more components (capped at 10) to show
            # how much a 2D projection is actually leaving on the table.
            n_full = min(10, x_scaled.shape[0] - 1, x_scaled.shape[1])
            pca_full = PCA(n_components=n_full, random_state=42).fit(x_scaled)
            scree_df = pd.DataFrame({
                "component": [f"PC{i + 1}" for i in range(n_full)],
                "variance_pct": pca_full.explained_variance_ratio_ * 100,
            })
            scree_fig = px.bar(
                scree_df, x="component", y="variance_pct",
                labels={"variance_pct": "% variance explained"},
                title="Scree plot — variance explained per component",
            )

            # Loadings: which original features drive PC1/PC2 the most.
            loadings = pd.DataFrame(pca.components_.T, index=keep_cols, columns=["PC1", "PC2"])
            magnitude = (loadings["PC1"] ** 2 + loadings["PC2"] ** 2) ** 0.5
            top_loadings = loadings.loc[magnitude.sort_values(ascending=False).head(15).index]
            loadings_long = top_loadings.reset_index().rename(columns={"index": "feature"}).melt(
                id_vars="feature", var_name="component", value_name="loading")
            loadings_fig = px.bar(
                loadings_long, x="feature", y="loading", color="component", barmode="group",
                title="Top feature loadings on PC1/PC2",
            )
        except Exception as e:  # noqa
            info += f" (scree/loadings charts unavailable: {e})"
    else:
        try:
            import umap
        except Exception as e:  # noqa
            return None, None, f"umap-learn is not installed/importable: {e}", None, None
        n_neighbors_eff = max(2, min(n_neighbors, len(sub) - 1))
        reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=n_neighbors_eff)
        coords = reducer.fit_transform(x_scaled)
        x_label, y_label = "UMAP-1", "UMAP-2"
        info = (f"UMAP is a nonlinear projection - unlike PCA, its axes have no "
                f"variance-explained percentage. Used {len(keep_cols)} feature(s), "
                f"n_neighbors={n_neighbors_eff}.")

    if dropped:
        shown = ", ".join(dropped[:8]) + (" …" if len(dropped) > 8 else "")
        info += f" Dropped {len(dropped)} constant/all-NaN column(s): {shown}."

    try:
        import plotly.express as px
    except Exception as e:  # noqa
        return None, None, f"plotly is not installed/importable: {e}", None, None

    plot_df = pd.DataFrame({x_label: coords[:, 0], y_label: coords[:, 1]}, index=sub.index)
    plot_df["molecule"] = sub.index.astype(str)
    if color_col and color_col in df.columns:
        plot_df[color_col] = df.loc[sub.index, color_col].values

    fig = px.scatter(
        plot_df, x=x_label, y=y_label,
        color=color_col if (color_col and color_col in plot_df.columns) else None,
        hover_name="molecule",
        title=f"{method} projection of {len(keep_cols)} feature(s)",
    )
    return fig, info, None, scree_fig, loadings_fig


def _coef_bar_chart(coef_df):
    """Signed bar chart of model coefficients (excludes the intercept), sorted by magnitude."""
    try:
        import plotly.express as px
    except Exception:
        return None
    if coef_df is None or "Estimate" not in coef_df.columns:
        return None
    work = coef_df.drop(index="(Intercept)", errors="ignore")
    if work.empty:
        return None
    plot_df = work[["Estimate"]].reset_index().rename(columns={"index": "feature"})
    plot_df = plot_df.reindex(plot_df["Estimate"].abs().sort_values(ascending=False).index)
    return px.bar(
        plot_df, x="feature", y="Estimate", color="Estimate",
        color_continuous_scale="RdBu", color_continuous_midpoint=0,
        title="Model coefficients (excl. intercept)",
    )


def _best_combo_plot(m3_model, combo):
    """Refit one feature combination and return a matplotlib Figure of actual vs. predicted."""
    import matplotlib.pyplot as plt
    from modeling import fit_and_evaluate_single_combination_regression

    res = fit_and_evaluate_single_combination_regression(m3_model, combo, persist_result=False)
    y_pred = res.get("predictions")
    y_true = m3_model.target_vector.to_numpy()
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    ax.scatter(y_true, y_pred, s=28, alpha=0.8)
    lo, hi = min(y_true.min(), y_pred.min()), max(y_true.max(), y_pred.max())
    ax.plot([lo, hi], [lo, hi], "--", color="gray", linewidth=1)
    ax.set_xlabel("actual")
    ax.set_ylabel("predicted (in-sample)")
    ax.set_title(" + ".join(combo), fontsize=9)
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------- #
# Tab 1 (structure/fingerprint viewer + scope page) helpers
# --------------------------------------------------------------------------- #
HARTREE_TO_KCAL = 627.5094740631


def _molecule_names_and_objs(feather_dir):
    """Load a (cached) feather Molecules set; return (mols_obj, {name: Molecule})."""
    try:
        mols = _load_feather_molecules_cached(feather_dir, root_dir, _dir_signature(feather_dir))
        by_name = {getattr(m, "molecule_name", f"mol_{i}"): m
                  for i, m in enumerate(mols.molecules)}
        return mols, by_name
    except Exception as e:  # noqa
        st.error(f"Could not load feather set: {e}")
        return None, {}


def _molecule_energy(mol):
    """Pull a single float energy value off a Molecule's energy_value df, or None."""
    try:
        ev = getattr(mol, "energy_value", None)
        if ev is None or len(ev) == 0:
            return None
        col = "energy" if "energy" in ev.columns else ev.columns[0]
        val = float(ev.iloc[0][col])
        return val if val == val else None  # filters NaN
    except Exception:
        return None


# --------------------------------------------------------------------------- #
# Tab 1 — Data
# --------------------------------------------------------------------------- #
with tab_data:
    st.subheader("Provide the molecule set")
    mode = st.radio("Source", ["Upload files", "Server folder paths"], horizontal=True)

    if mode == "Upload files":
        c1, c2 = st.columns(2)
        with c1:
            fe = st.file_uploader("Feather files (.feather)", type=["feather"],
                                  accept_multiple_files=True)
            if st.button("Use uploaded feathers"):
                ss["feather_dir"] = _save_uploads(fe, "feathers")
                st.success(f"{len(fe or [])} feathers saved")
        with c2:
            xy = st.file_uploader("XYZ files (.xyz)", type=["xyz"],
                                  accept_multiple_files=True)
            if st.button("Use uploaded xyz"):
                ss["xyz_dir"] = _save_uploads(xy, "xyz")
                st.success(f"{len(xy or [])} xyz saved")
    else:
        ss["feather_dir"] = st.text_input("feather_dir (on server)", value=ss["feather_dir"])
        ss["xyz_dir"] = st.text_input("xyz_dir (on server)", value=ss["xyz_dir"])

    st.write("Current set:")
    st.json({"feather_dir": ss["feather_dir"] or None, "xyz_dir": ss["xyz_dir"] or None})

    if ss["feather_dir"] or ss["xyz_dir"]:
        if st.button("Preview loaded molecules"):
            with st.spinner("Checking what's actually in there…"):
                ss["dataset_preview"] = _dataset_preview(ss["feather_dir"], ss["xyz_dir"])

        preview = ss.get("dataset_preview")
        if preview:
            name_filter = st.text_input(
                "Filter names (optional)", key="preview_name_filter",
                help="Narrows the lists below to names containing this text.").strip().lower()

            feather_info = preview.get("feather")
            if feather_info:
                if "error" in feather_info:
                    st.warning(f"Couldn't preview the feather set: {feather_info['error']}")
                else:
                    names = feather_info["names"]
                    shown_names = [n for n in names if name_filter in str(n).lower()] if name_filter else names
                    shown = ", ".join(shown_names[:10]) + (" …" if len(shown_names) > 10 else "")
                    count_label = (f"{len(shown_names)} of {feather_info['count']}" if name_filter
                                  else str(feather_info["count"]))
                    st.success(f"Feather set: {count_label} molecule(s) — {shown or '(no matches)'}")
                    if feather_info.get("failed"):
                        st.warning(f"Failed to load: {', '.join(feather_info['failed'])}")
            xyz_info = preview.get("xyz")
            if xyz_info:
                if "error" in xyz_info:
                    st.warning(f"Couldn't preview the xyz set: {xyz_info['error']}")
                else:
                    names = xyz_info["names"]
                    shown_names = [n for n in names if name_filter in str(n).lower()] if name_filter else names
                    shown = ", ".join(shown_names[:10]) + (" …" if len(shown_names) > 10 else "")
                    count_label = (f"{len(shown_names)} of {xyz_info['count']}" if name_filter
                                  else str(xyz_info["count"]))
                    st.success(f"XYZ set: {count_label} file(s) — {shown or '(no matches)'}")

    st.caption("Feather-only is fine — enable 'derive_xyz_from_feathers' in the config "
               "to run the xyz/SMILES engines.")

    st.divider()
    with st.expander("Energy overview"):
        st.caption(
            "Relative energies (kcal/mol) of every loaded molecule vs. the lowest-energy one — "
            "a quick sanity check before running descriptors, e.g. to spot outlier conformers."
        )
        if not ss["feather_dir"]:
            st.info("Load a feather set above first.")
        elif st.button("Compute energy overview"):
            with st.spinner("Reading energies…"):
                _, by_name_e = _molecule_names_and_objs(ss["feather_dir"])
                rows = []
                n_missing = 0
                for name, mol in by_name_e.items():
                    e = _molecule_energy(mol)
                    if e is None:
                        n_missing += 1
                    else:
                        rows.append({"molecule": name, "energy_hartree": e})
                if not rows:
                    st.warning("No molecules with usable energy values found.")
                else:
                    edf = pd.DataFrame(rows)
                    e_min = edf["energy_hartree"].min()
                    edf["dE_kcal_mol"] = (edf["energy_hartree"] - e_min) * HARTREE_TO_KCAL
                    edf = edf.sort_values("dE_kcal_mol")
                    try:
                        import plotly.express as px
                        efig = px.bar(
                            edf, x="molecule", y="dE_kcal_mol",
                            title="Relative energy vs. lowest-energy molecule",
                            labels={"dE_kcal_mol": "ΔE (kcal/mol)"},
                        )
                        st.plotly_chart(efig, use_container_width=True)
                    except Exception as e:  # noqa
                        st.error(f"Couldn't render the energy chart: {e}")
                        st.dataframe(edf, use_container_width=True)
                    msg = f"{len(rows)} molecule(s) with energy data"
                    if n_missing:
                        msg += f", {n_missing} missing energy"
                    st.caption(msg + ".")

    st.divider()
    with st.expander("RDKit structure & fingerprint viewer"):
        st.caption(
            "Pick a loaded molecule to see its 2D structure and which substructure "
            "lights up a given Morgan fingerprint bit — the same ECFP4-style bits "
            "the `rdkit_fp` engine extracts in tab 3, just visualized."
        )
        fp_names, fp_xyz_lookup = [], {}
        if ss["feather_dir"]:
            _, by_name = _molecule_names_and_objs(ss["feather_dir"])
            import make_picker as mp
            for name, mol_obj in by_name.items():
                fp_names.append(name)
                fp_xyz_lookup[name] = (lambda m=mol_obj: mp.molecule_to_xyz(m))
        elif ss["xyz_dir"]:
            for p in sorted(Path(ss["xyz_dir"]).glob("*.xyz")):
                fp_names.append(p.stem)
                fp_xyz_lookup[p.stem] = (lambda path=p: path.read_text(encoding="utf-8"))

        if not fp_names:
            st.info("Load a feather or xyz set above first.")
        else:
            if len(fp_names) > 15:
                fp_filter = st.text_input("Filter molecules by name", key="fp_name_filter").strip().lower()
                fp_names_filtered = [n for n in fp_names if fp_filter in str(n).lower()] if fp_filter else fp_names
                if fp_filter:
                    st.caption(f"{len(fp_names_filtered)} of {len(fp_names)} match.")
            else:
                fp_names_filtered = fp_names
            if not fp_names_filtered:
                st.warning("No molecule names match that filter.")
                picked_name = None
            else:
                picked_name = st.selectbox("Molecule", fp_names_filtered, key="fp_mol_name")
            n_bits = st.slider("Fingerprint size (bits)", 32, 512, 256, step=32, key="fp_n_bits")
            if st.button("Analyze structure & fingerprint", key="fp_analyze_btn",
                        disabled=not fp_names_filtered):
                import make_picker as mp
                xyz_text = fp_xyz_lookup[picked_name]()
                mol = mp.rdkit_mol_from_xyz(xyz_text)
                if mol is None:
                    st.error("RDKit couldn't perceive bonds for this molecule (can happen "
                             "on unusual geometries or charge states).")
                    ss["fp_mol_b64"] = None
                    ss["fp_on_bits"] = None
                else:
                    ss["fp_mol_b64"] = mp.mol_structure_png_base64(mol)
                    on_bits, bit_info = mp.morgan_fingerprint_bits(mol, n_bits=n_bits)
                    ss["fp_on_bits"] = on_bits
                    ss["fp_bit_info"] = bit_info
                    ss["fp_mol_obj"] = mol
                    st.success(f"{len(on_bits)} of {n_bits} bits are 'on' for {picked_name}.")

            if ss.get("fp_mol_b64"):
                st.image(f"data:image/png;base64,{ss['fp_mol_b64']}", width=320)
            if ss.get("fp_on_bits"):
                bit_choice = st.selectbox("Fingerprint bit to visualize", ss["fp_on_bits"],
                                          key="fp_bit_choice")
                if st.button("Show substructure for this bit", key="fp_bit_btn"):
                    mol = ss.get("fp_mol_obj")
                    if mol is None:
                        st.warning("Run 'Analyze structure & fingerprint' again first.")
                    else:
                        import make_picker as mp
                        b64 = mp.fingerprint_bit_png_base64(mol, bit_choice, ss["fp_bit_info"])
                        if b64:
                            st.image(f"data:image/png;base64,{b64}", caption=f"Bit {bit_choice}")
                        else:
                            st.warning("Couldn't render that bit.")

    st.divider()
    with st.expander("Build a scope page (RDKit structures + energies from the feather set)"):
        st.caption(
            "Renders every molecule in the feather set as a 2D RDKit structure, labeled "
            "with its name and the energy stored in the feather (as-is, typically Hartree "
            "from the parsed QM log) plus ΔE vs. the lowest-energy entry in kcal/mol — "
            "laid out like a synthesis scope figure."
        )
        if not ss["feather_dir"]:
            st.info("Load a feather set above first — the scope page needs the per-molecule "
                    "energy that's only stored there.")
        else:
            cols_per_row = st.slider("Structures per row", 2, 6, 4, key="scope_cols")
            _, by_name_all = _molecule_names_and_objs(ss["feather_dir"])
            all_names = list(by_name_all.keys())
            if len(all_names) > 15:
                scope_filter = st.text_input(
                    "Filter molecules by name (optional — narrows which ones get built below)",
                    key="scope_name_filter").strip().lower()
                filtered_names = ([n for n in all_names if scope_filter in str(n).lower()]
                                  if scope_filter else all_names)
                if scope_filter:
                    st.caption(f"{len(filtered_names)} of {len(all_names)} will be included.")
            else:
                filtered_names = all_names

            if st.button("Build scope page", type="primary", key="scope_build_btn",
                        disabled=not filtered_names):
                import make_picker as mp
                by_name = {n: by_name_all[n] for n in filtered_names}
                if not by_name:
                    st.warning("No molecules loaded.")
                else:
                    entries_raw = []
                    with st.spinner(f"Perceiving structures for {len(by_name)} molecule(s)…"):
                        for name, mol_obj in by_name.items():
                            xyz_text = mp.molecule_to_xyz(mol_obj)
                            rdmol = mp.rdkit_mol_from_xyz(xyz_text)
                            energy = _molecule_energy(mol_obj)
                            entries_raw.append((name, rdmol, energy))

                    energies = {n: e for n, _, e in entries_raw if e is not None}
                    e_min = min(energies.values()) if energies else None

                    entries, failed = [], []
                    for name, rdmol, energy in entries_raw:
                        if rdmol is None:
                            failed.append(name)
                            continue
                        if energy is not None:
                            legend = f"{name}\nE = {energy:.6f}"
                            if e_min is not None:
                                legend += f"\nΔE = {(energy - e_min) * HARTREE_TO_KCAL:.2f} kcal/mol"
                        else:
                            legend = f"{name}\n(no energy in feather)"
                        entries.append((rdmol, legend))

                    if not entries:
                        st.error("Couldn't build a structure for any molecule "
                                "(bond perception failed for all of them).")
                    else:
                        png_bytes = mp.scope_grid_png_bytes(entries, mols_per_row=cols_per_row)
                        if not png_bytes:
                            st.error("RDKit couldn't render the grid image.")
                        else:
                            ss["scope_page_png"] = png_bytes
                            if failed:
                                st.warning(f"Skipped {len(failed)} molecule(s) (bond perception "
                                          f"failed): {', '.join(failed)}")
                            st.success(f"Built a scope page with {len(entries)} structure(s).")

            if ss.get("scope_page_png"):
                st.image(ss["scope_page_png"], use_container_width=True)
                st.download_button(
                    "Download scope page (PNG)",
                    ss["scope_page_png"],
                    file_name="scope_page.png",
                    mime="image/png",
                )


# --------------------------------------------------------------------------- #
# Tab 2 — Pick atoms (embedded 3D picker)
# --------------------------------------------------------------------------- #
with tab_pick:
    st.subheader("Pick atoms in 3D, then download run_config.json")
    st.caption("Build your selections + tick engines in the picker, then either "
               "**Copy** the answers_dict/config and paste it into tab 3's "
               "**Paste from picker** option, or **Save run_config.json** and upload it there.")
    load_choice = st.selectbox(
        "Load into the picker",
        ["Demo molecule", "First molecule of the feather set", "Upload an .xyz here"])

    xyz_text, mol_data = None, None
    if load_choice == "First molecule of the feather set" and ss["feather_dir"]:
        xyz_text, mol_data = _first_feather_molecule_xyz(ss["feather_dir"])
    elif load_choice == "Upload an .xyz here":
        up = st.file_uploader("xyz for the picker", type=["xyz"], key="pickxyz")
        if up:
            xyz_text = up.getvalue().decode("utf-8", "ignore")

    try:
        import streamlit.components.v1 as components
        components.html(_build_picker_html(xyz_text, mol_data), height=900, scrolling=True)
    except Exception as e:  # noqa
        st.error(f"Could not embed the picker: {e}")


# --------------------------------------------------------------------------- #
# Tab 3 — Configure & run
# --------------------------------------------------------------------------- #
with tab_run:
    st.subheader("Configure the run")

    with st.expander("Descriptor glossary"):
        st.markdown(
            """
- **Sterimol L / B1 / B5** — Verloop's steric parameters for a substituent: `L` is length
  along the primary bond axis, `B1` is the minimum width, `B5` is the maximum width.
- **%Vbur (buried volume)** — the percentage of a sphere around a metal/central atom that is
  occupied by the ligand's atoms; a common ligand-size descriptor in organometallic chemistry.
- **Cone angle** — the apex angle of the narrowest cone, centered on the metal, that contains
  all ligand atoms (Tolman cone angle and its generalizations).
- **SASA** — solvent-accessible surface area: the surface area of a molecule (or fragment)
  reachable by a solvent probe sphere.
- **Dispersion (P_int)** — a Sterimol-derived descriptor capturing London dispersion
  interaction potential, based on polarizability/volume along the substituent.
- **Pyramidalization** — a measure of how far an atom (e.g. a trigonal center) deviates from
  perfect planarity, often reported in degrees or as an improper-angle-derived index.
- **VIF (variance inflation factor)** — quantifies multicollinearity between features in a
  regression; VIF > ~5-10 for a feature signals it's highly correlated with other features.
- **Adjusted R² / Q²** — adjusted R² penalizes R² for the number of predictors used; Q² is a
  cross-validated analogue of R² (predictive, out-of-sample performance) — a model that fits
  well (high R²) but predicts poorly (low Q²) is likely overfit.
            """
        )

    src = st.radio(
        "Config source",
        ["Build in app", "Paste from picker", "Upload run_config.json", "Edit full JSON"],
        horizontal=True,
    )

    cfg = None
    template = _template_config()
    if src == "Paste from picker":
        st.caption(
            "Paste either the full `run_config.json` text, or just the "
            "**Python answers_dict** / atoms JSON copied from the picker's Export panel "
            "(tab 2) — no need to save/download a file first. Either syntax works."
        )
        pasted_cfg_text = st.text_area("Paste here", height=260, key="pasted_picker_cfg")
        if pasted_cfg_text.strip():
            parsed, perr = _parse_json_or_pyliteral(pasted_cfg_text)
            if perr:
                st.error(f"Couldn't parse that: {perr}")
            elif not isinstance(parsed, dict):
                st.error("Expected a JSON/Python object (a `{...}` dict).")
            elif "engines" in parsed:
                cfg = parsed
                st.success("Parsed as a full run_config.")
            else:
                # Looks like just the atoms/answers_dict -> merge into the template.
                cfg = dict(template)
                cfg["engines"] = dict(template.get("engines", {}))
                full_cfg = dict(cfg["engines"].get("descripytor_full", {}))
                full_cfg["enabled"] = True
                full_cfg["atoms"] = parsed
                cfg["engines"]["descripytor_full"] = full_cfg
                st.success("Parsed as an atoms/answers_dict — merged into the descripytor_full engine.")
    elif src == "Upload run_config.json":
        cf = st.file_uploader("run_config.json (exported by the picker)", type=["json"])
        if cf:
            try:
                cfg = json.loads(cf.getvalue().decode("utf-8"))
            except Exception as e:  # noqa
                st.error(f"Invalid JSON upload: {e}")
    elif src == "Edit full JSON":
        txt = st.text_area("Config JSON", value=json.dumps(template, indent=2), height=420)
        try:
            cfg = json.loads(txt)
        except Exception as e:  # noqa
            st.error(f"Invalid JSON: {e}")
    else:
        st.caption(
            "Grouped by what each engine actually computes, so you're not stuck "
            "picking between several similarly-named options. Where two or three "
            "backends compute the same quantity (e.g. Sterimol from the feather "
            "set vs. straight from xyz vs. Morfeus), enter the atom pairs once and "
            "tick whichever backend(s) you want — tick more than one to "
            "cross-check them against each other."
        )
        engine_templates = template.get("engines", {})
        required = ("descripytor_full", "descripytor_steric", "xyz_sterimol", "morfeus_sterimol",
                   "xyz_geometry", "xyz_buried_volume", "morfeus_suite",
                   "rdkit", "rdkit_fp", "mordred", "deepchem", "rafbl", "qm", "aqme_qdescp")
        if not all(name in engine_templates for name in required):
            st.error("Engine templates unavailable (toolkit import failed) — "
                     "check the DescriPyTor root in the sidebar, or use 'Edit full JSON' instead.")
            cfg = None
        else:
            cfg = {key: val for key, val in template.items() if key != "engines"}
            cfg["engines"] = {}
            any_error = False

            # ---- 1. Core DescriPytor descriptor set ---------------------- #
            use_full = st.checkbox(
                "DescriPytor full descriptor set (IR, dipole, charges, sterimol, bond "
                "length/angle…)", value=True, key="grp_full")
            if use_full:
                with st.expander("Atom selections (descripytor_full)", expanded=True):
                    atoms, err = _json_area(
                        "Atoms JSON", engine_templates["descripytor_full"].get("atoms", {}),
                        key="builder_atoms_json", height=260,
                    )
                    if err:
                        st.error(f"Invalid atoms JSON: {err}")
                        any_error = True
                    else:
                        full_cfg = dict(engine_templates["descripytor_full"])
                        full_cfg["enabled"] = True
                        full_cfg["atoms"] = atoms
                        cfg["engines"]["descripytor_full"] = full_cfg

            st.divider()

            # ---- 2. Sterimol: one shared atom-pair input, pick backend(s) - #
            st.markdown("**Sterimol** — same atom pairs, different backend(s):")
            use_sterimol = st.checkbox("Compute Sterimol", value=False, key="grp_sterimol")
            if use_sterimol:
                c1, c2, c3 = st.columns(3)
                backend_feather = c1.checkbox("DescriPytor (feather)", value=True, key="ster_backend_feather")
                backend_xyz = c2.checkbox("Standalone (xyz)", value=False, key="ster_backend_xyz")
                backend_morfeus = c3.checkbox("Morfeus (reference check)", value=False, key="ster_backend_morfeus")
                pairs, perr = _json_area(
                    "Atom pairs (shared by every backend ticked above)",
                    engine_templates["xyz_sterimol"].get("pairs", [[6, 7], [6, 4]]),
                    key="sterimol_pairs_json", height=90,
                )
                if perr:
                    st.error(f"Invalid pairs JSON: {perr}")
                    any_error = True
                else:
                    if backend_feather:
                        fcfg = dict(engine_templates["descripytor_steric"])
                        fcfg.update(enabled=True, pairs=pairs)
                        cfg["engines"]["descripytor_steric"] = fcfg
                    if backend_xyz:
                        radii = st.selectbox("xyz backend radii", ["CPK", "bondi"], key="ster_radii")
                        xcfg = dict(engine_templates["xyz_sterimol"])
                        xcfg.update(enabled=True, pairs=pairs, radii=radii)
                        cfg["engines"]["xyz_sterimol"] = xcfg
                    if backend_morfeus:
                        mscfg = dict(engine_templates["morfeus_sterimol"])
                        mscfg.update(enabled=True, pairs=pairs, prefix="morfeus_ref_")
                        cfg["engines"]["morfeus_sterimol"] = mscfg
                    if not (backend_feather or backend_xyz or backend_morfeus):
                        st.info("Tick at least one backend above, or turn off 'Compute Sterimol'.")

            st.divider()

            # ---- 3. Geometry (xyz) ---------------------------------------- #
            st.markdown(
                "**Geometry** — angles & bond lengths straight from the xyz set. "
                "*(`descripytor_full`'s own `bond_length`/`bond_angle` atom fields do the "
                "same job on the feather set — use one or the other, not both, to avoid "
                "duplicate columns.)*"
            )
            use_geom = st.checkbox("Compute geometry (xyz_geometry)", value=False, key="grp_geometry")
            if use_geom:
                with st.expander("Geometry settings", expanded=True):
                    gcfg_default = dict(engine_templates["xyz_geometry"])
                    gcfg_default["enabled"] = True
                    gcfg, gerr = _json_area("Engine JSON", gcfg_default,
                                           key="builder_engine_xyz_geometry", height=140)
                    if gerr:
                        st.error(f"Invalid JSON: {gerr}")
                        any_error = True
                    else:
                        gcfg["enabled"] = True
                        cfg["engines"]["xyz_geometry"] = gcfg

            # ---- 4. Buried volume, single anchor (xyz) --------------------- #
            st.markdown(
                "**Buried volume** — quick single-anchor version (xyz). For several "
                "anchor descriptors at once (buried volume + cone angle + SASA…), use "
                "the Morfeus suite below instead."
            )
            use_bv = st.checkbox("Compute buried volume (xyz_buried_volume)", value=False, key="grp_bv")
            if use_bv:
                with st.expander("Buried volume settings", expanded=True):
                    bcfg_default = dict(engine_templates["xyz_buried_volume"])
                    bcfg_default["enabled"] = True
                    bcfg, berr = _json_area("Engine JSON", bcfg_default,
                                           key="builder_engine_xyz_bv", height=100)
                    if berr:
                        st.error(f"Invalid JSON: {berr}")
                        any_error = True
                    else:
                        bcfg["enabled"] = True
                        cfg["engines"]["xyz_buried_volume"] = bcfg

            st.divider()

            # ---- 5. Morfeus, kept separate from the descriptor engines ---- #
            st.markdown(
                "**Morfeus** — independent reference descriptors (Sterimol, buried volume, "
                "cone angle, SASA, dispersion, pyramidalization) computed with the `morfeus` "
                "package, kept in one place rather than mixed in with the descriptor engines "
                "above."
            )
            use_morfeus = st.checkbox("Compute Morfeus descriptors (morfeus_suite)",
                                      value=False, key="grp_morfeus")
            if use_morfeus:
                with st.expander("Morfeus settings", expanded=True):
                    available_desc = ["sterimol", "buried_volume", "cone_angle",
                                      "sasa", "dispersion", "pyramidalization"]
                    default_desc = engine_templates["morfeus_suite"].get(
                        "descriptors", ["sterimol", "buried_volume", "sasa", "dispersion"])
                    chosen_desc = st.multiselect("Descriptors", available_desc,
                                                 default=default_desc, key="morfeus_desc")
                    mcfg_default = dict(engine_templates["morfeus_suite"])
                    mcfg_default["enabled"] = True
                    mcfg_default["descriptors"] = chosen_desc
                    mcfg, merr = _json_area(
                        "Engine JSON (anchor atoms, radii…)", mcfg_default,
                        key="builder_engine_morfeus_suite", height=200)
                    if merr:
                        st.error(f"Invalid JSON: {merr}")
                        any_error = True
                    else:
                        mcfg["enabled"] = True
                        mcfg["descriptors"] = chosen_desc
                        cfg["engines"]["morfeus_suite"] = mcfg

            st.divider()

            # ---- 6. External cheminformatics (SMILES-based) ---------------- #
            st.markdown("**External cheminformatics** — computed from SMILES perceived off the xyz set.")
            ext_names = ["rdkit", "rdkit_fp", "mordred", "deepchem", "rafbl"]
            ext_chosen = st.multiselect("Enable", ext_names, default=[], key="grp_external")
            for name in ext_chosen:
                with st.expander(f"{name} settings", expanded=False):
                    block_default = dict(engine_templates[name])
                    block_default["enabled"] = True
                    block, err = _json_area("Engine JSON", block_default,
                                           key=f"builder_engine_{name}", height=140)
                    if err:
                        st.error(f"Invalid JSON for {name}: {err}")
                        any_error = True
                    else:
                        block["enabled"] = True
                        cfg["engines"][name] = block

            # ---- 7. QM / high-cost ------------------------------------------ #
            st.markdown("**QM / high-cost** — Gaussian log parsing or AQME/xTB (can be slow).")
            qm_names = ["qm", "aqme_qdescp"]
            qm_chosen = st.multiselect("Enable", qm_names, default=[], key="grp_qm")
            for name in qm_chosen:
                with st.expander(f"{name} settings", expanded=False):
                    block_default = dict(engine_templates[name])
                    block_default["enabled"] = True
                    block, err = _json_area("Engine JSON", block_default,
                                           key=f"builder_engine_{name}", height=160)
                    if err:
                        st.error(f"Invalid JSON for {name}: {err}")
                        any_error = True
                    else:
                        block["enabled"] = True
                        cfg["engines"][name] = block

            if any_error:
                cfg = None

    derive = st.checkbox("derive_xyz_from_feathers (feather-only sets that use xyz engines)",
                         value=False)

    if cfg is not None:
        cfg = _overlay_run_paths(cfg, derive_xyz=derive)
        ss["last_config"] = cfg

        enabled = _enabled_engines(cfg)
        st.write("Engines enabled:", ", ".join(enabled) or "none")
        blockers = _run_readiness(cfg)
        for blocker in blockers:
            st.warning(blocker)

        try:
            atom_idx_problems = _check_atom_indices_against_feather(cfg, cfg.get("feather_dir"))
        except Exception:  # noqa
            atom_idx_problems = []
        if atom_idx_problems:
            with st.expander(
                f"⚠️ Atom-index check found {len(atom_idx_problems)} molecule(s) too small "
                "for the configured atom indices",
                expanded=True,
            ):
                st.caption(
                    "Best-effort static check: an engine config references an atom index "
                    "larger than a molecule's atom count. This will likely fail or silently "
                    "misbehave for these molecules — double-check the atom indices in tab 2."
                )
                for name, n_atoms, max_idx in atom_idx_problems:
                    st.write(f"- **{name}**: has {n_atoms} atoms, but max referenced index is {max_idx}")

        st.json({
            "root_dir": cfg.get("root_dir"),
            "feather_dir": cfg.get("feather_dir") or None,
            "xyz_dir": cfg.get("xyz_dir") or None,
            "derive_xyz_from_feathers": cfg.get("derive_xyz_from_feathers"),
            "output_csv": cfg.get("output_csv"),
        })

        st.download_button(
            "Download generated run_config.json",
            json.dumps(cfg, indent=2).encode("utf-8"),
            file_name="run_config.json",
            mime="application/json",
        )

        if ss.get("extraction_thread") is None:
            if st.button("Run extraction", type="primary", disabled=bool(blockers)):
                try:
                    from descriptor_extractor import run_from_config
                except Exception as e:  # noqa
                    st.error(f"Toolkit import failed (check DescriPyTor root): {e}")
                else:
                    buf = io.StringIO()
                    outcome = {}
                    cancel_event = threading.Event()

                    def _run_worker(cfg=cfg, buf=buf, outcome=outcome, cancel_event=cancel_event):
                        try:
                            with contextlib.redirect_stdout(buf):
                                outcome["df"] = run_from_config(cfg, stop_check=cancel_event.is_set)
                        except Exception as e:  # noqa
                            outcome["error"] = e

                    thread = threading.Thread(target=_run_worker, daemon=True)
                    ss["extraction_thread"] = thread
                    ss["extraction_buf"] = buf
                    ss["extraction_outcome"] = outcome
                    ss["extraction_cancel_event"] = cancel_event
                    ss["extraction_cfg"] = cfg
                    thread.start()
                    st.rerun()
        else:
            # A run is in progress (or just finished) - rendered via rerun-driven
            # polling rather than a blocking while-loop, so the Cancel button
            # below is actually clickable mid-run (Streamlit can't react to a
            # click while the script is stuck in a blocking loop).
            ex_thread = ss["extraction_thread"]
            ex_buf = ss["extraction_buf"]
            ex_outcome = ss["extraction_outcome"]
            ex_cancel_event = ss["extraction_cancel_event"]
            ex_cfg = ss["extraction_cfg"]

            text = ex_buf.getvalue()
            engine_lines = [ln for ln in text.splitlines() if ln.startswith("=== engine:")]
            if ex_cancel_event.is_set() and ex_thread.is_alive():
                label = "Cancelling — finishing the current engine…"
            elif engine_lines:
                label = f"Running engines… ({engine_lines[-1].strip('= ')})"
            else:
                label = "Running engines…"

            state = "running" if ex_thread.is_alive() else ("error" if "error" in ex_outcome else "complete")
            status_box = st.status(label, expanded=True, state=state)
            status_box.code(text[-4000:] or "(starting…)", language="text")

            if ex_thread.is_alive():
                if st.button("Cancel run", key="extraction_cancel_btn", disabled=ex_cancel_event.is_set()):
                    ex_cancel_event.set()
                    st.info("Cancelling — will stop before the next engine starts "
                           "(the one currently running still finishes).")
                time.sleep(0.6)
                st.rerun()
            else:
                ss["run_logs"] = text
                if "error" in ex_outcome:
                    ss["result_df"] = None
                    st.error(f"Run failed: {ex_outcome['error']}")
                elif ex_cancel_event.is_set():
                    ss["result_df"] = ex_outcome.get("df")
                    ss["last_output_csv"] = ex_cfg["output_csv"]
                    st.warning("Run was cancelled — Results shows whatever engines finished before that.")
                else:
                    ss["result_df"] = ex_outcome.get("df")
                    ss["last_output_csv"] = ex_cfg["output_csv"]
                    st.success("Done — see the Results tab.")
                # clear the in-progress state so the Run button reappears
                ss["extraction_thread"] = None
                ss["extraction_buf"] = None
                ss["extraction_outcome"] = None
                ss["extraction_cancel_event"] = None
                ss["extraction_cfg"] = None


# --------------------------------------------------------------------------- #
# Tab 4 — Results
# --------------------------------------------------------------------------- #
with tab_res:
    st.subheader("Merged features")
    df = ss.get("result_df")
    if df is not None and len(df):
        st.write(f"{df.shape[0]} molecules × {df.shape[1]} features")
        st.dataframe(df, use_container_width=True)
        csv_path = ss.get("last_output_csv")
        if csv_path and Path(csv_path).exists():
            st.caption(f"Saved on server: `{csv_path}`")
            csv_bytes = Path(csv_path).read_bytes()
        else:
            csv_bytes = df.to_csv().encode("utf-8")
        st.download_button("Download merged_features.csv",
                           csv_bytes,
                           file_name="merged_features.csv", mime="text/csv")
    else:
        st.info("No results yet — run an extraction in tab 3.")
    if ss.get("run_logs"):
        with st.expander("Run log"):
            st.code(ss["run_logs"])

    st.divider()
    st.caption("Bundle everything from this session — run_config.json, the merged CSV, the "
               "run log, and any M3/BASSA modeling outputs generated in tab 5 — into one zip.")
    _has_exportable = bool(
        ss.get("last_config") or (df is not None and len(df)) or
        (ss.get("m3_results") is not None and len(ss["m3_results"])) or
        ss.get("bassa_run_dir") or ss.get("scope_page_png")
    )
    if st.button("Prepare export bundle", disabled=not _has_exportable):
        ss["export_zip_bytes"] = _build_export_zip()
    if ss.get("export_zip_bytes"):
        st.download_button(
            "Download everything (zip)", ss["export_zip_bytes"],
            file_name="descripytor_export.zip", mime="application/zip",
        )

    st.divider()
    st.subheader("Dimensionality-reduction scatter")
    st.caption("Project the merged features to 2D to see whether the molecule set spreads out "
               "or clusters meaningfully in descriptor space.")
    if df is None or not len(df):
        st.info("No results yet — run an extraction in tab 3.")
    else:
        numeric_cols = df.select_dtypes(include="number").columns.tolist()
        if len(numeric_cols) < 2:
            st.info("Need at least 2 numeric feature columns to compute a projection.")
        else:
            proj_method = st.radio("Method", ["PCA", "UMAP"], horizontal=True, key="proj_method")
            proj_cols = st.multiselect(
                "Features to include", options=numeric_cols, default=numeric_cols, key="proj_cols",
                help="Constant or all-NaN columns are dropped automatically.",
            )
            color_choice = st.selectbox("Color by", ["(none)"] + numeric_cols, key="proj_color")
            color_col = None if color_choice == "(none)" else color_choice
            n_neighbors = (st.slider("UMAP n_neighbors", 2, 50, 15, key="proj_umap_neighbors")
                          if proj_method == "UMAP" else 15)

            if st.button("Compute projection", key="proj_compute_btn"):
                if len(proj_cols) < 2:
                    st.error("Pick at least 2 feature columns.")
                else:
                    with st.spinner(f"Computing {proj_method} projection…"):
                        fig, info, error, scree_fig, loadings_fig = _build_projection_figure(
                            df, proj_cols, method=proj_method,
                            color_col=color_col, n_neighbors=n_neighbors)
                    if error:
                        st.error(error)
                        ss["proj_fig"] = None
                        ss["proj_info"] = None
                        ss["proj_scree_fig"] = None
                        ss["proj_loadings_fig"] = None
                    else:
                        ss["proj_fig"] = fig
                        ss["proj_info"] = info
                        ss["proj_scree_fig"] = scree_fig
                        ss["proj_loadings_fig"] = loadings_fig

            if ss.get("proj_fig") is not None:
                st.plotly_chart(ss["proj_fig"], use_container_width=True)
                if ss.get("proj_info"):
                    st.caption(ss["proj_info"])
                if ss.get("proj_scree_fig") is not None or ss.get("proj_loadings_fig") is not None:
                    with st.expander("PCA diagnostics (scree plot + loadings)"):
                        if ss.get("proj_scree_fig") is not None:
                            st.plotly_chart(ss["proj_scree_fig"], use_container_width=True)
                        if ss.get("proj_loadings_fig") is not None:
                            st.plotly_chart(ss["proj_loadings_fig"], use_container_width=True)

    st.divider()
    with st.expander("Feature exploration (correlation, distributions, missing data)"):
        if df is None or not len(df):
            st.info("No results yet — run an extraction in tab 3.")
        else:
            numeric_cols_fe = df.select_dtypes(include="number").columns.tolist()
            if not numeric_cols_fe:
                st.info("No numeric columns to explore.")
            else:
                st.markdown("**Correlation heatmap**")
                corr_cols = st.multiselect(
                    "Columns to include", numeric_cols_fe,
                    default=numeric_cols_fe[:30], key="corr_cols",
                    help="Capped at 30 by default for readability — add/remove as needed.",
                )
                if len(corr_cols) >= 2:
                    try:
                        import plotly.express as px
                        corr = df[corr_cols].corr()
                        heat_fig = px.imshow(
                            corr, text_auto=".2f", aspect="auto",
                            color_continuous_scale="RdBu", zmin=-1, zmax=1,
                            title="Feature correlation",
                        )
                        st.plotly_chart(heat_fig, use_container_width=True)
                    except Exception as e:  # noqa
                        st.error(f"Couldn't render the heatmap: {e}")
                else:
                    st.info("Pick at least 2 columns.")

                st.markdown("**Column distribution**")
                hist_col = st.selectbox("Column", numeric_cols_fe, key="hist_col")
                try:
                    import plotly.express as px
                    hist_fig = px.histogram(df, x=hist_col, title=f"Distribution of {hist_col}")
                    st.plotly_chart(hist_fig, use_container_width=True)
                except Exception as e:  # noqa
                    st.error(f"Couldn't render the histogram: {e}")

                st.markdown("**Missing data**")
                null_counts = df.isna().sum()
                null_counts = null_counts[null_counts > 0].sort_values(ascending=False)
                if null_counts.empty:
                    st.success("No missing values in the merged features.")
                else:
                    total_cells = df.shape[0] * df.shape[1]
                    total_missing = int(null_counts.sum())
                    st.caption(
                        f"{total_missing} of {total_cells} cells missing "
                        f"({100 * total_missing / total_cells:.1f}%), across {len(null_counts)} column(s)."
                    )
                    try:
                        import plotly.express as px
                        miss_fig = px.bar(
                            x=null_counts.index.astype(str), y=null_counts.values,
                            labels={"x": "column", "y": "missing count"},
                            title="Missing values per column",
                        )
                        st.plotly_chart(miss_fig, use_container_width=True)
                    except Exception as e:  # noqa
                        st.error(f"Couldn't render missing-data chart: {e}")

    st.divider()
    with st.expander("Compare two columns (parity plot — e.g. cross-check Sterimol backends)"):
        st.caption(
            "Pick any two numeric columns to scatter against each other with a y = x reference "
            "line — most useful right after running Sterimol through more than one backend "
            "(tab 3) to see how closely they agree."
        )
        if df is None or not len(df):
            st.info("No results yet — run an extraction in tab 3.")
        else:
            numeric_cols_cmp = df.select_dtypes(include="number").columns.tolist()
            if len(numeric_cols_cmp) < 2:
                st.info("Need at least 2 numeric columns.")
            else:
                cc1, cc2 = st.columns(2)
                col_x = cc1.selectbox("X column", numeric_cols_cmp, key="parity_col_x")
                default_y_idx = 1 if len(numeric_cols_cmp) > 1 else 0
                col_y = cc2.selectbox("Y column", numeric_cols_cmp, index=default_y_idx, key="parity_col_y")
                if col_x == col_y:
                    st.info("Pick two different columns.")
                else:
                    sub = df[[col_x, col_y]].apply(pd.to_numeric, errors="coerce").dropna()
                    if len(sub) < 2:
                        st.warning("Not enough overlapping non-null values between these two columns.")
                    else:
                        try:
                            import plotly.express as px
                            import plotly.graph_objects as go
                            corr = sub[col_x].corr(sub[col_y])
                            rmse = ((sub[col_x] - sub[col_y]) ** 2).mean() ** 0.5
                            fig = px.scatter(
                                sub, x=col_x, y=col_y, hover_name=sub.index.astype(str),
                                title=f"{col_x} vs {col_y}",
                            )
                            lo = min(sub[col_x].min(), sub[col_y].min())
                            hi = max(sub[col_x].max(), sub[col_y].max())
                            fig.add_trace(go.Scatter(
                                x=[lo, hi], y=[lo, hi], mode="lines",
                                line=dict(dash="dash", color="gray"), name="y = x",
                            ))
                            st.plotly_chart(fig, use_container_width=True)
                            st.caption(f"n={len(sub)}, Pearson r = {corr:.3f}, RMSE = {rmse:.4g}.")
                        except Exception as e:  # noqa
                            st.error(f"Couldn't render the parity plot: {e}")


# --------------------------------------------------------------------------- #
# Tab 5 — Modeling (comes after Results; models the merged CSV)
# --------------------------------------------------------------------------- #
with tab_model:
    st.subheader("Model the extracted CSV")
    st.caption("Runs on the Results CSV (or one you upload here) — M3 for exhaustive "
               "linear/lasso feature-combination search, BASSA for Bayesian spike-and-slab "
               "variable selection.")

    src = st.radio("Data source", ["Results CSV", "Upload a CSV"], horizontal=True, key="model_src")
    df_model = None
    if src == "Results CSV":
        df_model = ss.get("result_df")
        csv_path = ss.get("last_output_csv")
        if (df_model is None or not len(df_model)) and csv_path and Path(csv_path).exists():
            df_model = pd.read_csv(csv_path)
        if df_model is None or not len(df_model):
            st.info("No results yet — run an extraction (tab 3) or upload a CSV below.")
    else:
        up_model_csv = st.file_uploader("Descriptor CSV", type=["csv"], key="model_csv_upload")
        if up_model_csv:
            df_model = pd.read_csv(up_model_csv)

    if df_model is not None and len(df_model):
        st.write(f"{df_model.shape[0]} rows × {df_model.shape[1]} columns")
        numeric_cols = df_model.select_dtypes(include="number").columns.tolist()

        with st.expander("Target column", expanded=not numeric_cols):
            st.caption("If the descriptor CSV doesn't yet contain the measured target, "
                       "paste the values here — matched to CSV row order, or as "
                       "`name,value` pairs (matched against a name/molecule column). "
                       "This writes a new in-memory dataset; the descriptor CSV on disk "
                       "is not touched.")
            pasted = st.text_area(
                "Output values to merge (optional)",
                placeholder="Row order:\n0.42\n0.55\n0.61\n\nor name,value pairs:\nmol_a,0.42\nmol_b,0.55",
                key="pasted_outputs",
            )
            new_col_name = st.text_input("New column name", value="output", key="pasted_outputs_col")
            if st.button("Merge into a modeling dataset"):
                pairs = _parse_pasted_outputs(pasted)
                if not pairs:
                    st.error("Couldn't parse any values out of that text.")
                else:
                    merged = _merge_outputs_into_df(df_model, pairs, new_col_name)
                    ss["model_df"] = merged
                    ss["model_df_source"] = "merged"
                    st.success(f"Merged {len(pairs)} value(s) into column '{new_col_name}'.")

        active_df = ss["model_df"] if ss.get("model_df_source") == "merged" else df_model
        active_numeric = active_df.select_dtypes(include="number").columns.tolist()

        if not active_numeric:
            st.warning("No numeric columns available to model. Merge an output column above.")
        else:
            target_col = st.selectbox("Target column (y)", options=active_numeric, key="model_target")
            default_feats = [c for c in active_numeric if c != target_col]
            feature_cols = st.multiselect(
                "Feature columns (X)", options=default_feats, default=default_feats, key="model_features"
            )

            if ss.get("model_df_source") == "merged":
                st.download_button(
                    "Download merged modeling CSV",
                    active_df.to_csv(index=False).encode("utf-8"),
                    file_name="modeling_dataset.csv", mime="text/csv",
                )

            if not feature_cols:
                st.info("Pick at least one feature column.")
            else:
                engine = st.radio(
                    "Engine", ["M3 — linear / lasso search", "BASSA — Bayesian spike-and-slab"],
                    horizontal=True, key="model_engine",
                )

                # ----------------------------------------------------------- #
                # M3 feature-combination search
                # ----------------------------------------------------------- #
                if engine.startswith("M3"):
                    c1, c2, c3, c4 = st.columns(4)
                    model_type = c1.selectbox("Model type", ["linear", "lasso"], key="m3_model_type")
                    min_feat = c2.number_input("Min features", min_value=1, max_value=len(feature_cols),
                                               value=1, key="m3_min_feat")
                    max_feat = c3.number_input("Max features", min_value=min_feat, max_value=len(feature_cols),
                                               value=min(3, len(feature_cols)), key="m3_max_feat")
                    top_n = c4.number_input("Top N", min_value=1, max_value=500, value=20, key="m3_top_n")
                    threshold = st.slider("R² threshold to keep a combination", 0.0, 1.0, 0.5, key="m3_threshold")
                    n_combos = sum(
                        1 for k in range(int(min_feat), int(max_feat) + 1)
                        for _ in __import__("itertools").combinations(feature_cols, k)
                    )
                    st.caption(f"~{n_combos:,} feature combinations will be evaluated.")

                    LARGE_COMBO_THRESHOLD = 20_000
                    run_blocked = False
                    if n_combos > LARGE_COMBO_THRESHOLD:
                        st.warning(
                            f"That's {n_combos:,} combinations to fit and evaluate — this can take a "
                            "long time and there's no cancel button once it starts. Consider narrowing "
                            "the feature list or the min/max feature range, or confirm below."
                        )
                        run_blocked = not st.checkbox(f"Yes, run all {n_combos:,} combinations anyway")

                    st.caption("Unlike extraction in tab 3, there's no Cancel button for this once it "
                              "starts — it runs as a single blocking call. The confirmation above is "
                              "the only safety net, so double-check the combination count first.")
                    if st.button("Run M3 search", type="primary", disabled=run_blocked):
                        try:
                            _add_m3_to_path()
                            from modeling import LinearRegressionModel
                        except Exception as e:  # noqa
                            st.error(f"Could not import M3_modeler.modeling: {e}")
                        else:
                            sub_df = active_df[feature_cols + [target_col]].copy()
                            db_dir = Path(ss["work_dir"]) / "m3_runs"
                            db_dir.mkdir(parents=True, exist_ok=True)
                            with st.spinner("M3 is evaluating feature combinations…"):
                                try:
                                    m3_model = LinearRegressionModel(
                                        sub_df,
                                        process_method="one csv",
                                        y_value=target_col,
                                        min_features_num=int(min_feat),
                                        max_features_num=int(max_feat),
                                        model_type=model_type,
                                        db_path=str(db_dir / "results"),
                                    )
                                    results = m3_model.search_models(top_n=int(top_n), threshold=float(threshold))
                                    ss["m3_results"] = results
                                    ss["m3_model"] = m3_model
                                except Exception as e:  # noqa
                                    st.error(f"M3 search failed: {e}")
                                    ss["m3_results"] = None

                    if ss.get("m3_results") is not None and len(ss["m3_results"]):
                        results = ss["m3_results"]
                        st.dataframe(results, use_container_width=True)
                        st.download_button(
                            "Download M3 results CSV",
                            results.to_csv(index=False).encode("utf-8"),
                            file_name="m3_results.csv", mime="text/csv",
                        )
                        if ss.get("m3_model") is not None and "combination" in results.columns:
                            def _format_result_option(i, _results=results):
                                row = _results.iloc[i]
                                bits = [str(row.get("combination", ""))]
                                if "adj_r2" in _results.columns:
                                    bits.append(f"adjR²={row['adj_r2']:.3f}")
                                elif "r2" in _results.columns:
                                    bits.append(f"R²={row['r2']:.3f}")
                                if "q2" in _results.columns:
                                    bits.append(f"Q²={row['q2']:.3f}")
                                return f"[{i}] " + "  ·  ".join(bits)

                            row_idx = st.selectbox(
                                "Result to plot / report on (top = best)",
                                options=list(range(len(results))),
                                format_func=_format_result_option,
                                key="m3_report_row",
                            )
                            combo_str = results.iloc[int(row_idx)]["combination"]
                            combo = None
                            try:
                                from modeling import _normalize_combination_to_columns
                                combo = _normalize_combination_to_columns(combo_str)
                            except Exception as e:  # noqa
                                st.caption(f"(couldn't parse that combination: {e})")

                            if combo:
                                try:
                                    fig = _best_combo_plot(ss["m3_model"], combo)
                                    st.pyplot(fig)
                                except Exception as e:  # noqa
                                    st.caption(f"(couldn't render the quick scatter: {e})")

                                if st.button("Generate full report (scatter, violin, VIF, CV, SHAP, diagnostics…)"):
                                    try:
                                        from plot import run_single_combo_report
                                    except Exception as e:  # noqa
                                        st.error(f"Could not import M3_modeler.plot: {e}")
                                    else:
                                        report_dir = Path(ss["work_dir"]) / "m3_reports"
                                        report_dir.mkdir(parents=True, exist_ok=True)
                                        pdf_path = report_dir / f"report_row{int(row_idx)}.pdf"
                                        with st.spinner("Building the full model report…"):
                                            try:
                                                ss["m3_report"] = run_single_combo_report(
                                                    ss["m3_model"], combo, show=True,
                                                    pdf_name=str(pdf_path),
                                                )
                                            except Exception as e:  # noqa
                                                st.error(f"Report generation failed: {e}")
                                                ss["m3_report"] = None

                                report = ss.get("m3_report")
                                if report:
                                    figs = report.get("figures") or {}
                                    if figs.get("scatter") is not None:
                                        st.pyplot(figs["scatter"])
                                    if figs.get("violin") is not None:
                                        st.pyplot(figs["violin"])
                                    for label, df_key in (
                                        ("Cross-validation metrics", "folds_df"),
                                        ("VIF", "vif_df"),
                                        ("Coefficients", "coef_df"),
                                    ):
                                        df_val = report.get(df_key)
                                        if df_val is not None and len(df_val):
                                            st.write(f"{label}:")
                                            st.dataframe(df_val, use_container_width=True)
                                            if df_key == "coef_df":
                                                coef_fig = _coef_bar_chart(df_val)
                                                if coef_fig is not None:
                                                    st.plotly_chart(coef_fig, use_container_width=True)
                                    leftout = report.get("leftout") or {}
                                    if leftout.get("df") is not None and len(leftout["df"]):
                                        st.write("Held-out predictions:")
                                        st.dataframe(leftout["df"], use_container_width=True)
                                    if report.get("pdf_path") and Path(report["pdf_path"]).exists():
                                        st.download_button(
                                            "Download full PDF report",
                                            Path(report["pdf_path"]).read_bytes(),
                                            file_name=Path(report["pdf_path"]).name,
                                            mime="application/pdf",
                                        )
                                    extra_pngs = report.get("png_files") or []
                                    if extra_pngs:
                                        with st.expander(
                                            f"All report pages ({len(extra_pngs)} images: "
                                            "SHAP, thresholds, diagnostics, sanity checks…)"
                                        ):
                                            for p in extra_pngs:
                                                st.image(p, caption=Path(p).name, use_container_width=True)

                # ----------------------------------------------------------- #
                # BASSA — Bayesian spike-and-slab
                # ----------------------------------------------------------- #
                else:
                    try:
                        from bassa_reg import Bassa
                        from bassa_reg.spike_and_slab.spike_and_slab import (
                            SpikeAndSlabConfigurations, SpikeAndSlabRegression,
                        )
                        from bassa_reg.spike_and_slab.spike_and_slab_util_models import SpikeAndSlabPriors
                    except Exception as e:  # noqa
                        st.error(f"bassa-reg is not installed/importable: {e}. "
                                "It ships in the Docker image's environment.yml (pip: bassa-reg).")
                    else:
                        c1, c2 = st.columns(2)
                        iterations = c1.number_input("Sampler iterations", min_value=100, max_value=100_000,
                                                     value=5000, step=500, key="bassa_iterations")
                        experiment_name = c2.text_input("Experiment name", value="streamlit_run",
                                                        key="bassa_experiment")
                        if st.button("Run BASSA", type="primary"):
                            project_path = Path(ss["work_dir"]) / "bassa_runs"
                            project_path.mkdir(parents=True, exist_ok=True)
                            x = active_df[feature_cols]
                            y = active_df[target_col]
                            with st.spinner("BASSA is sampling models (spike-and-slab)…"):
                                try:
                                    priors = SpikeAndSlabPriors()
                                    config = SpikeAndSlabConfigurations(sampler_iterations=int(iterations))
                                    regression = SpikeAndSlabRegression(
                                        x=x, y=y, priors=priors, config=config,
                                        project_path=str(project_path), experiment_name=experiment_name,
                                    )
                                    regression.run()
                                    bassa = Bassa(model=regression)
                                    bassa.run()
                                    ss["bassa_run_dir"] = str(project_path / regression.full_experiment_name)
                                    st.success("BASSA run complete.")
                                except Exception as e:  # noqa
                                    st.error(f"BASSA run failed: {e}")

                        if ss.get("bassa_run_dir"):
                            run_dir = Path(ss["bassa_run_dir"])
                            if run_dir.exists():
                                for fname, caption in (
                                    ("bassa_plot.png", "Models chosen by BASSA"),
                                    ("markov_chain_visualization.png", "Feature exploration over iterations"),
                                    ("survival_plot.png", "Model survival over iterations"),
                                ):
                                    fpath = run_dir / fname
                                    if fpath.exists():
                                        st.image(str(fpath), caption=caption, use_container_width=True)
                                stats_csv = run_dir / "feature_stats.csv"
                                if stats_csv.exists():
                                    st.write("Feature inclusion frequencies:")
                                    st.dataframe(pd.read_csv(stats_csv), use_container_width=True)
                            else:
                                st.caption(f"Run directory not found on disk: `{run_dir}`")
