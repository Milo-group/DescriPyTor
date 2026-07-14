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
    r"C:\Users\edens\Desktop\DescriPytor\DescriPyTor-main\MolFeatures",
)

st.set_page_config(page_title="DescriPytor Descriptor Extractor", layout="wide")
ss = st.session_state
ss.setdefault("work_dir", tempfile.mkdtemp(prefix="descripytor_web_"))
ss.setdefault("feather_dir", "")
ss.setdefault("xyz_dir", "")
ss.setdefault("result_df", None)
ss.setdefault("run_logs", "")
ss.setdefault("last_config", None)
ss.setdefault("last_output_csv", "")


# --------------------------------------------------------------------------- #
# Sidebar — environment
# --------------------------------------------------------------------------- #
st.sidebar.header("Environment")
root_dir = st.sidebar.text_input(
    "MolFeatures root", value=DEFAULT_ROOT,
    help="Folder that contains M2_data_extractor etc. Must exist on THIS server.")
os.environ["DESCRIPYTOR_ROOT"] = root_dir
st.sidebar.caption(f"Working dir: `{ss['work_dir']}`")
st.sidebar.caption("Engines whose Python deps (rdkit, morfeus-ml, xTB, aqme…) "
                   "aren't installed here will skip themselves with a note.")

st.title("DescriPytor — Descriptor Extraction")
tab_data, tab_pick, tab_run, tab_res = st.tabs(
    ["1 · Data", "2 · Pick atoms", "3 · Configure & run", "Results"])


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
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
        import descriptor_extractor as dx
        dx._add_descripytor_to_path(root_dir)   # adds root, utils, M2/M3, MolAlign
        import make_picker as mp
        from data_extractor import Molecules
        mols = Molecules(feather_dir)
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
    st.caption("Feather-only is fine — enable 'derive_xyz_from_feathers' in the config "
               "to run the xyz/SMILES engines.")


# --------------------------------------------------------------------------- #
# Tab 2 — Pick atoms (embedded 3D picker)
# --------------------------------------------------------------------------- #
with tab_pick:
    st.subheader("Pick atoms in 3D, then download run_config.json")
    st.caption("Build your selections + tick engines in the picker, click "
               "**Save run_config.json**, then upload it in tab 3.")
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
    src = st.radio(
        "Config source",
        ["Build in app", "Upload run_config.json", "Edit full JSON"],
        horizontal=True,
    )

    cfg = None
    template = _template_config()
    if src == "Upload run_config.json":
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
            "Choose engines here, edit their small JSON blocks only if needed, "
            "then run or download the generated config."
        )
        engine_templates = template.get("engines", {})
        default_engines = [
            name for name, ecfg in engine_templates.items()
            if ecfg.get("enabled")
        ] or ["descripytor_full"]
        enabled_names = st.multiselect(
            "Engines to run",
            options=list(engine_templates),
            default=[name for name in default_engines if name in engine_templates],
        )

        cfg = {key: val for key, val in template.items() if key != "engines"}
        cfg["engines"] = {}

        atoms_default = engine_templates.get("descripytor_full", {}).get("atoms", {})
        if "descripytor_full" in enabled_names:
            with st.expander("descripytor_full atom selections", expanded=True):
                atoms, err = _json_area(
                    "Atoms JSON",
                    atoms_default,
                    key="builder_atoms_json",
                    height=260,
                )
                if err:
                    st.error(f"Invalid atoms JSON: {err}")
                    cfg = None
                else:
                    full_cfg = dict(engine_templates["descripytor_full"])
                    full_cfg["enabled"] = True
                    full_cfg["atoms"] = atoms
                    cfg["engines"]["descripytor_full"] = full_cfg

        if cfg is not None:
            for name in enabled_names:
                if name == "descripytor_full":
                    continue
                with st.expander(f"{name} settings", expanded=False):
                    block_default = dict(engine_templates[name])
                    block_default["enabled"] = True
                    block, err = _json_area(
                        "Engine JSON",
                        block_default,
                        key=f"builder_engine_{name}",
                        height=180,
                    )
                    if err:
                        st.error(f"Invalid JSON for {name}: {err}")
                        cfg = None
                        break
                    block["enabled"] = True
                    cfg["engines"][name] = block

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

        if st.button("Run extraction", type="primary", disabled=bool(blockers)):
            try:
                from descriptor_extractor import run_from_config
            except Exception as e:  # noqa
                st.error(f"Toolkit import failed (check MolFeatures root): {e}")
            else:
                buf = io.StringIO()
                run_ok = False
                with st.spinner("Running engines…"):
                    try:
                        with contextlib.redirect_stdout(buf):
                            df = run_from_config(cfg)
                        ss["result_df"] = df
                        ss["last_output_csv"] = cfg["output_csv"]
                        run_ok = True
                    except Exception as e:  # noqa
                        ss["result_df"] = None
                        st.error(f"Run failed: {e}")
                ss["run_logs"] = buf.getvalue()
                if not run_ok:
                    st.stop()
                st.success("Done — see the Results tab.")


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
