import sys
import os
if not os.environ.get("MPLBACKEND"):
    os.environ["MPLBACKEND"] = "Agg"
import gc
import time
import argparse
import inspect
import json
import re
from datetime import datetime
from typing import Optional, List


def _project_root() -> str:
    """Repository root when running from a clone; site-packages after install."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _ensure_project_on_path() -> None:
    """Allow both `python __main__.py` and installed-package execution."""
    root = _project_root()
    if root not in sys.path:
        sys.path.insert(0, root)


def _load_runtime_dependencies(include_modeling: bool = False) -> None:
    """
    Import optional scientific and GUI dependencies only when a command needs them.

    Keeping these imports lazy lets `python __main__.py --help` work even when a
    local scientific environment is temporarily inconsistent.
    """
    _ensure_project_on_path()

    global pd, shutil, subprocess, warnings, webbrowser, traceback, urllib
    import pandas as pd
    import shutil
    import subprocess
    import warnings
    import webbrowser
    import traceback
    import urllib.request

    global help_functions, file_handlers, Molecules, Module1Handler
    global logs_to_feather, cube_many, Molecules_xyz
    import utils.help_functions as help_functions
    import utils.file_handlers as file_handlers
    from M2_data_extractor.data_extractor import Molecules
    from M1_pre_calculations.main import Module1Handler
    from M2_data_extractor.feather_extractor import logs_to_feather
    from M2_data_extractor.cube_sterimol import cube_many
    from M2_data_extractor.sterimol_standalone import Molecules_xyz

    if include_modeling:
        global ClassificationModel, LinearRegressionModel
        from M3_modeler.modeling import ClassificationModel, LinearRegressionModel


def load_answers_json(path = None):

    file_path = path
    if not file_path:
        return {}
    with open(file_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    print(f'printing loaded data from json file: {data}')
    return data


def get_local_setup_version():
    """Return the installed package version."""
    try:
        from descripytor import __version__ as version
        return version
    except Exception:
        return "0.4.7"

__version__ = get_local_setup_version()

# Function to get the current date
def get_current_date():
    return datetime.now().strftime("%Y-%m-%d")

            


def convert_to_custom_nested_list(input_str):
    """
    Converts a comma-separated string into a nested list based on the format:
    - "1,2,3,4,5,6" -> [[1,2,3,4], 5, 6]
    - "1,2,3,4,5,6 2,3,1,5,6,7" -> [[[1,2,3,4], 5, 6], [[2,3,1,5], 6, 7]]
    """
   
    split_by_space = input_str.split(' ')  # Split by space for multiple sections

    def process_sublist(sublist_str):
        """Convert a single comma-separated list to the required nested structure."""
        elements = list(map(int, sublist_str.split(',')))  # Convert to list of integers
        if len(elements) > 2:  # Ensure we separate all but the last two elements
            return [elements[:-2]] + elements[-2:]
        return elements  # If fewer than 3 elements, return as-is

    # Process each segment and decide if it's a nested list or a single flat list
    if len(split_by_space) == 1:
        # Single segment, no spaces
        return process_sublist(split_by_space[0])
    else:
        # Multiple segments separated by spaces
        nested_list = []
        for sublist_str in split_by_space:
            nested_list.append(process_sublist(sublist_str))
        return nested_list


def convert_to_list_or_nested_list(input_str):
    
    split_by_space = input_str.split(' ')

    # If there are no spaces, return a flat list
    if len(split_by_space) == 1:
        try:
            return list(map(int, split_by_space[0].split(',')))
        except:
            return list(map(int, split_by_space[0]))
    
    # Otherwise, return a nested list
    nested_list = []
    for sublist_str in split_by_space:
        sublist = list(map(int, sublist_str.split(',')))
        nested_list.append(sublist)

    if nested_list:
        return nested_list
    else:
        return None
   


def resolve_n_jobs(user_n_jobs) -> int:
    """Pick n_jobs from CLI, NSLOTS, or all logical cores."""
    if user_n_jobs is not None:
        return int(user_n_jobs)
    env_slots = os.environ.get("NSLOTS")
    if env_slots:
        try:
            return int(env_slots)
        except ValueError:
            pass
    import psutil
    return int(psutil.cpu_count(logical=True))


def resolve_effective_n_jobs(requested_n_jobs: int, bool_parallel: bool) -> int:
    """
    Normalize the runtime worker count for CLI model execution.

    `n_jobs` is treated as the primary control. A value greater than 1 enables
    parallel execution even when `--bool-parallel` is omitted, while
    `--bool-parallel` remains supported for backward compatibility.
    """
    requested_n_jobs = int(requested_n_jobs)
    if requested_n_jobs == -1:
        import psutil
        requested_n_jobs = int(psutil.cpu_count(logical=True))

    requested_n_jobs = max(1, requested_n_jobs)
    if bool_parallel or requested_n_jobs > 1:
        return requested_n_jobs
    return 1

def add_modeling_path():
    """Add M3_modeler to sys.path if a clone-layout folder exists."""
    path = os.path.join(_project_root(), "M3_modeler")
    if os.path.isdir(path) and path not in sys.path:
        sys.path.append(path)

def safe_chdir(workdir, features_csv):
    """
    Set working directory:
    - If --workdir given, chdir there.
    - Else chdir to the directory containing features_csv.
    """
    if workdir and os.path.isdir(workdir):
        os.chdir(workdir)
        return
    if features_csv:
        abs_feat = os.path.abspath(features_csv)
        base = os.path.dirname(abs_feat)
        if os.path.isdir(base):
            os.chdir(base)

def run_experiment(
    mode: str,
    features_csv: str,
    target_csv: str = "",
    y_value: str = "output",
    n_jobs: int = 1,
    min_features_num: int = 2,
    max_features_num: int = 4,
    threshold: float = 0.70,
    top_n: int = 50,
    leave_out: Optional[List[str]] = None,
    process_method: str = "one csv",
    bool_parallel: bool = False,
) -> float:
    """
    Run regression or classification and return elapsed seconds.
    - For regression, calls: fit_and_evaluate_combinations(...) if available, else search_models(...)
    - For classification, calls: search_models(...)
    """
    _load_runtime_dependencies(include_modeling=True)


    csv_filepaths = {
        "features_csv_filepath": features_csv,
        "target_csv_filepath": target_csv or ""
    }

    leave_out = list(leave_out or [])
    if target_csv and process_method == "one csv":
        process_method = "two csvs"      # -t given: the target comes from its own file

    # Instantiate correct model class
    if mode.lower() == "regression":
        model = LinearRegressionModel(
            csv_filepaths,
            process_method=process_method,
            y_value=y_value,
            leave_out=leave_out,
            min_features_num=min_features_num,
            max_features_num=max_features_num,
           
        )
    elif mode.lower() == "classification":
        # Assumes you have a ClassificationModel class with similar init signature
        model = ClassificationModel(
            csv_filepaths,
            process_method=process_method,
            y_value=y_value,
            leave_out=leave_out,
            min_features_num=min_features_num,
            max_features_num=max_features_num,
    
        )
    else:
        raise ValueError("mode must be 'regression' or 'classification'")

    # Decide parallel cores
    eff_n_jobs = resolve_effective_n_jobs(n_jobs, bool_parallel)
    print(f"\n[run_experiment] mode={mode} | n_jobs={eff_n_jobs} (requested={n_jobs}, parallel={bool_parallel})")
    print(f"[run_experiment] features={features_csv} | target={target_csv or '(none)'} | y_value={y_value}")
    print(f"[run_experiment] min/max features: {min_features_num}/{max_features_num}\n")

    # Run
    start = time.time()
    if mode.lower() == "regression":
        if hasattr(model, "search_models"):
            model.search_models(
                top_n=top_n,
                n_jobs=eff_n_jobs,
                threshold=threshold,
                bool_parallel=bool_parallel
            )
        else:
            raise AttributeError("Regression model lacks search_models.")
    else:  # classification
        if hasattr(model, "search_models"):
            model.search_models(
                top_n=top_n,
                n_jobs=eff_n_jobs,
                threshold=threshold,
                bool_parallel=bool_parallel
            )
        else:
            raise AttributeError("Classification model lacks search_models().")

    elapsed = time.time() - start
    print(f"\n[run_experiment] Completed in {elapsed:.2f} s\n")
    gc.collect()
    return elapsed

    # Your code to launch the GUI app goes here


def run_visual_app(host=None, port=None, open_browser=True, page="/visual"):
    """Open the 3D atom picker in a browser."""
    _ensure_project_on_path()
    from M2_data_extractor.gui_server import serve

    serve(host=host, port=port, open_browser=open_browser, page=page)


CLI_HELP = """DescriPyTor: molecular descriptors from quantum-chemistry output, and the models built on them.

Pages (each opens in the browser; the three link to each other):
  visual_extractor   pick atoms on a molecule, extract descriptors for a folder of structures, fit a model
  feature_explorer   draw any descriptor on the molecule: Sterimol, cone angle, %Vbur, dipole, charges, lengths, angles
  bassa_explorer     the models a BASSA run visited, beside the molecule
  bassa_setup        once: an environment for BASSA's sampler, which new tables in bassa_explorer need

Batch commands:
  extractor          descriptors for a folder of .feather files, from an input file of atom choices
  model              search regression or classification models over a feature table
  design             does a model find the best compounds, and which candidate to make next
  joint              choose descriptors on related reactions together, each reaction with its own weights
  forms              logistic link for yields, a GP / forest ceiling, threshold analysis
  audit              columns one compound dominates (z-audit); cut a wide table to one column per correlation cluster
  bassa_protocol     a kept BASSA case's validation: 20 deals of 5-fold + LOO, chance line, seed check, drop-one
  pipeline          SMILES -> structures -> chained cluster stages (GOAT / UMA / xTB / ORCA / Gaussian)
  renumber           number every molecule of a folder like a reference (Conformaty)

Converters (they ask for their inputs):
  logs_to_feather    Gaussian .log files -> .feather
  sterimol           Sterimol from .xyz files
  cube               Sterimol from Gaussian cube files"""
CLI_EPILOG = """examples:
  descripytor visual_extractor
  descripytor feature_explorer --mol FL_lig_13 --preset "theta construction"
  descripytor bassa_setup
  descripytor bassa_explorer --cases my_cases
  descripytor extractor -i input.json -o features -f feathers/
  descripytor model -m regression -f features.csv -t targets.csv --min-features 2 --max-features 3
  descripytor design table.csv -y ddG --columns B1 loc_B5 --lower --candidates new.csv --batch 3
  descripytor joint features.csv outcomes.csv --columns q_Ni q_donor_mean B1 L
  descripytor forms table.csv -y yield --columns cone vbur --logit
  descripytor audit wide.csv -y ddG --reduce 50
  descripytor bassa_protocol protocol my_case
  descripytor pipeline status protocol.json
  descripytor renumber feathers/ -r basic

`descripytor <command> --help` lists the options of one command. Older names still work:
visual = visual_extractor, theta = feature_explorer, bassa = bassa_explorer.
Manuals: https://github.com/Milo-group/DescriPyTor/blob/main/docs/WORKFLOW.md, FEATURE_EXPLORER.md, BASSA_EXPLORER.md, PIPELINE.md (same folder)."""


DOCS = "https://github.com/Milo-group/DescriPyTor/blob/main/docs/"


def run_bassa_explorer(cases=None, host=None, port=None, open_browser=True, bassa_python=None):
    """Start the BASSA model explorer (M3_modeler/bassa_explorer/server.py) and open it.

    The server runs under ``bassa_python`` (or DESCRIPYTOR_BASSA_PYTHON), the Python that has ``bassa_reg``; by default
    this one. Without ``bassa_reg`` the page opens the example and every kept case, and says that new tables need it.
    Structures are always read by the Python running this command. Returns the server's exit code.
    """
    import subprocess
    import threading
    import webbrowser
    from pathlib import Path

    server = Path(__file__).resolve().parents[1] / "M3_modeler" / "bassa_explorer" / "server.py"
    env = dict(os.environ)
    env.setdefault("BASSA_EXPLORER_READER", sys.executable)
    if cases:
        env["BASSA_EXPLORER_CASES"] = str(Path(cases).resolve())
        try:                                                 # remembered for starts without --cases (the visual extractor's link)
            _ensure_project_on_path()
            from M3_modeler.bassa_explorer.settings import REMEMBERED
            REMEMBERED.parent.mkdir(parents=True, exist_ok=True)
            REMEMBERED.write_text(env["BASSA_EXPLORER_CASES"], encoding="utf-8")
        except (OSError, ImportError):
            pass
    host, port = host or "127.0.0.1", int(port or os.environ.get("BASSA_EXPLORER_PORT", 8765))
    env["BASSA_EXPLORER_PORT"] = str(port)                   # so the other pages' links find this server
    url = f"http://{'127.0.0.1' if host in ('0.0.0.0', '') else host}:{port}/"
    print(f"\n  BASSA model explorer\n  Open:  {url}\n  Press Ctrl+C to stop\n", flush=True)
    python = _sampler_python(bassa_python)
    from M3_modeler.bassa_explorer import running
    already = running.state(port)
    if already in ("current", running.BUSY):                 # this version's server is up (or busy with a run): use it
        print("  It is already running%s; this process will not start a second server.%s"
              % (" (busy: it answered slowly)" if already == running.BUSY else "", chr(10)), flush=True)
        if open_browser:
            webbrowser.open(url)
        return 0
    if already == "old":                                     # a server from before an update would give a new page on old routes
        print("  " + running.stop_hint(port, "the BASSA model explorer") + chr(10), flush=True)
        return 1

    def open_when_up():                                      # the first import of the sampler can take a while
        import socket
        import time
        for _ in range(240):
            with socket.socket() as s:
                s.settimeout(0.3)
                if s.connect_ex(("127.0.0.1", port)) == 0:
                    return webbrowser.open(url)
            time.sleep(0.25)
    if open_browser:
        threading.Thread(target=open_when_up, daemon=True).start()
    try:
        return subprocess.call([python, str(server), str(port), host], env=env)
    except KeyboardInterrupt:
        return 0


def run_bassa_protocol(args):
    """Run M3_modeler/bassa_explorer/protocol.py (validation protocol, chance line, seed check, drop-one) under the
    sampler's Python, the one bassa_explorer uses. Returns its exit code."""
    import subprocess
    from pathlib import Path

    script = Path(__file__).resolve().parents[1] / "M3_modeler" / "bassa_explorer" / "protocol.py"
    env = dict(os.environ)
    if args.cases_dir:
        env["BASSA_EXPLORER_CASES"] = str(Path(args.cases_dir).resolve())
    python = _sampler_python(args.bassa_python, remember=False)
    rest = [args.check, *args.cases, "--deals", str(args.deals), "--shuffles", str(args.shuffles), "--seed", str(args.seed)]
    rest += ["--no-loo"] * args.no_loo + ["--dry"] * args.dry
    try:
        return subprocess.call([python, str(script), *rest], env=env)
    except KeyboardInterrupt:
        return 1


def _sampler_python(bassa_python=None, remember=True):
    """The Python with bassa_reg, first found wins: ``bassa_python`` (remembered when ``remember``),
    DESCRIPYTOR_BASSA_PYTHON, the one remembered from an earlier --bassa-python or bassa_setup, this Python.
    Prints which, so a remembered one is never a surprise."""
    from pathlib import Path

    _ensure_project_on_path()
    from M3_modeler.bassa_explorer.settings import REMEMBERED_PY
    source = "--bassa-python"
    if bassa_python:                                        # remembered for starts without it (the visual extractor's link)
        if remember:
            try:
                REMEMBERED_PY.parent.mkdir(parents=True, exist_ok=True)
                REMEMBERED_PY.write_text(str(bassa_python), encoding="utf-8")
            except OSError:
                pass
    elif os.environ.get("DESCRIPYTOR_BASSA_PYTHON"):
        bassa_python, source = os.environ["DESCRIPYTOR_BASSA_PYTHON"], "DESCRIPYTOR_BASSA_PYTHON"
    else:
        try:
            kept = REMEMBERED_PY.read_text(encoding="utf-8").strip()
            if Path(kept).is_file():
                bassa_python, source = kept, "remembered in %s" % REMEMBERED_PY
        except OSError:
            pass
    python = bassa_python or sys.executable
    print("  Sampler Python: %s (%s)%s" % (python, source if bassa_python else "this one", chr(10)), flush=True)
    return python


def run_bassa_setup(where=None, python=None):
    """Make the environment BASSA's sampler runs in, once, and remember it for bassa_explorer.

    bassa-reg pins NumPy 2.0 and its other dependencies exactly, while DescriPyTor needs NumPy < 2 (its RDKit and PyArrow
    builds), so the two cannot share one Python. This makes a virtual environment (default ~/.descripytor/bassa_env) with
    bassa-reg alone, checks that the sampler imports, and writes its Python where bassa_explorer looks. Returns 0 on success.
    """
    import subprocess
    import venv
    from pathlib import Path

    _ensure_project_on_path()
    from M3_modeler.bassa_explorer.settings import REMEMBERED_PY
    home = Path(where).expanduser().resolve() if where else REMEMBERED_PY.parent / "bassa_env"
    py = home / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    if not py.is_file():
        print("  Making %s ..." % home, flush=True)
        if python:
            subprocess.check_call([python, "-m", "venv", str(home)])
        else:
            venv.create(str(home), with_pip=True)
    spec = "bassa-reg>=0.1.16,<0.2"
    print("  Installing %s into it (a few minutes: it brings its own NumPy, SciPy and Numba) ..." % spec, flush=True)
    if subprocess.call([str(py), "-m", "pip", "install", "--upgrade", spec]):
        print("  pip failed: see the lines above. Nothing was remembered.", flush=True)
        return 1
    check = subprocess.run([str(py), "-c", "from bassa_reg.spike_and_slab.spike_and_slab import SpikeAndSlabConfigurations, SpikeAndSlabRegression; "
                            "import importlib.metadata as m; print(m.version('bassa-reg'))"], capture_output=True, text=True)
    if check.returncode:
        print("  The sampler did not import:" + chr(10) + check.stderr.strip()[-800:], flush=True)
        return 1
    REMEMBERED_PY.parent.mkdir(parents=True, exist_ok=True)
    REMEMBERED_PY.write_text(str(py), encoding="utf-8")
    print("  Done: bassa-reg %s in %s. bassa_explorer uses it from now on (said at its start)." % (check.stdout.strip(), py), flush=True)
    return 0


def feature_explorer_page(mol=None, preset=None):
    """The explorer's address, with the structure and preset the URL parameters ask for."""
    from urllib.parse import urlencode

    query = urlencode({k: v for k, v in (("mol", mol), ("preset", preset)) if v})
    return "/feature_explorer" + ("?" + query if query else "")


theta_page = feature_explorer_page          # the name before the explorer was renamed

def run_feature_extraction(input_file, output_file = 'features_set', molecules_dir_name=None):
    _load_runtime_dependencies()
    if not molecules_dir_name:
        from descripytor.examples import feather_example_dir
        molecules_dir_name = str(feather_example_dir())
    answers = load_answers_json(input_file)
    mols = load_molecules(molecules_dir_name, renumber=False)
    mols.get_molecules_features_set(answers, save_as=True, csv_file_name=output_file)

def load_molecules(molecules_dir_name, renumber=False):
    _load_runtime_dependencies()
    return Molecules(molecules_dir_name, renumber=renumber)


def _bar(done, n, width=24):
    full = round(width * done / n) if n else width
    return "#" * full + "." * (width - full)


def run_renumber(args):
    """`descripytor renumber`: every molecule of a folder numbered like the reference, written to a new folder with
    renumbering.json (the {old: new} dict per molecule). The originals are never changed."""
    import os
    import time

    from M2_data_extractor.alignment import map_folder, write_renumbered
    from M2_data_extractor.data_extractor import Molecules

    folder = os.path.abspath(args.folder)
    mols = Molecules(folder)
    sizes = {m.molecule_name: len(m.xyz_df) for m in mols.molecules}
    if not sizes:
        raise SystemExit(f"No molecules in {folder}")
    reference = args.reference or min(sizes, key=lambda k: (sizes[k], k))      # the common part: the smallest
    if reference not in sizes:
        raise SystemExit(f"No molecule named {reference!r} in {folder}. Names: {', '.join(sorted(sizes))}")
    anchors = {}
    for text in args.anchor or []:
        name, _, pair = text.partition(":")
        ref_atom, _, atom = pair.partition(",")
        anchors.setdefault(name, []).append((int(ref_atom), int(atom)))
    out = os.path.abspath(args.output or os.path.join(folder, "renumbered"))
    t0 = time.perf_counter()

    def progress(ev):
        if ev["event"] == "start":
            print(f"Numbering {ev['n']} molecules like {ev['reference']} ({sizes[ev['reference']]} atoms), "
                  f"{ev['workers']} at a time.")
            return
        done, n = ev["i"], ev["n"]
        left = (time.perf_counter() - t0) / done * (n - done)
        if "error" in ev:
            what = "FAILED: " + ev["error"]
        else:
            what = "same numbering" if not ev["moved"] else f"{ev['moved']} atoms renumbered"
            if ev["differs"]:
                what += f"; other element at {', '.join(map(str, ev['differs']))}"
            what += f"; RMSD {ev['rmsd']:.2f} A" if ev["rmsd"] is not None else ""
            if ev["rmsd"] is not None and ev.get("rmsd_kabsch") is not None and ev["rmsd"] - ev["rmsd_kabsch"] > 0.5:
                what += f" (the {args.fit} fit stalled: Kabsch gives {ev['rmsd_kabsch']:.2f} A)"
        eta = f"~{left / 60:.0f} min left" if left >= 90 else f"~{left:.0f} s left" if done < n else "done"
        print(f"[{_bar(done, n)}] {done:>{len(str(n))}}/{n}  {ev['name']:<16} {ev.get('seconds', 0):5.1f} s  {what}  ({eta})",
              flush=True)

    results = map_folder(mols, reference, anchors=anchors, stages=args.stages or None,
                         match_charge=args.match_charge, workers=args.workers, progress=progress,
                         fit=args.fit, sigma=args.sigma, fit_atoms=args.fit_atoms)
    written, skipped = write_renumbered(mols, results, out, reference)
    print(f"\nWrote {len(written)} renumbered molecules and renumbering.json to {out} "
          f"in {time.perf_counter() - t0:.0f} s.")
    if skipped:
        print(f"Not written (no numbering): {', '.join(skipped)}")
    return 1 if skipped else 0


def run_pipeline(args):
    import json as _json
    from M1_pre_calculations import pipeline as pl

    p = pl.load(args.protocol)
    if args.action == "build":
        rows = pl.build_all(p, strict=not args.keep_going)
        for r in rows:
            print(f"{r['id']:>8}  {r['status']:<8} {r.get('n_atoms', ''):>4}  {r.get('formula', '')}  {r['name']}")
        print(f"built {sum(r['status'] == 'built' for r in rows)} of {len(rows)} into {p.workdir}")
    elif args.action == "adopt":
        if not args.sources:
            raise SystemExit("adopt needs --sources")
        print(_json.dumps(pl.adopt(p, args.stage, pl.read_sources(args.sources)), indent=2))
    elif args.action == "submit":
        only = [i for i in args.only.split(",") if i.strip()] if args.only is not None else None
        if only == []:      # an empty list from a shell variable used to mean "every molecule"
            raise SystemExit("--only was given but names no molecule; leave it out to submit all")
        if args.only_adopted:
            f = p.workdir / "adopted.txt"
            only = f.read_text(encoding="utf-8").split() if f.exists() else []
        print(_json.dumps(pl.submit(p, dry_run=args.dry_run, force=args.force, from_stage=args.from_stage, only=only), indent=2))
    elif args.action == "status":
        rows, queue, checks = pl.status(p)
        stages = p.stage_dirs()
        print("	".join(["id"] + stages))
        for r in rows:
            print("	".join([r["id"]] + [r.get(s, "") for s in stages]))
        for s, c in pl.summarize(p, rows).items():
            print(f"{s}: {c}")
        if checks:
            print("checks (done* = finished with warnings):")
            print(checks)
        print("queue:", queue or "empty")
    elif args.action == "retry":
        plan = pl.retry(p, dry_run=args.dry_run, all_failures=args.all)
        print(_json.dumps(plan, indent=2) if plan else "nothing to retry")
    elif args.action == "reset":
        if not args.only:
            raise SystemExit("reset needs --only ids and --from-stage")
        pl.reset(p, args.only.split(","), args.from_stage)
        print("reset", args.only, "from stage", args.from_stage)
    elif args.action == "watch":
        print(pl.watch(p, interval=args.interval))
    elif args.action == "fetch":
        print("fetched into", pl.fetch(p, all_stages=args.all_stages))
    return 0


def main():
    # A Windows console with a legacy code page (cp1252, cp1255, ...) cannot print every
    # character a report may contain; replace those instead of stopping the run.
    from utils.console import make_console_safe
    make_console_safe()
    raw = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(prog="descripytor", usage="descripytor [-h] [-V] <command> [options]", formatter_class=raw, description=CLI_HELP, epilog=CLI_EPILOG)
    parser.add_argument("-V", "--version", action="version", version="descripytor " + __version__)
    subparsers = parser.add_subparsers(dest="command", metavar="<command>", help=argparse.SUPPRESS)      # listed, grouped, in the description

    visual_parser = subparsers.add_parser(
        "visual_extractor", aliases=["visual"], help="Open the 3D atom picker in a browser", formatter_class=raw,
        description="The atom picker: load a .xyz or .feather, click the atoms each descriptor is read at (the molecule is drawn by\n"
                    "the feature explorer), extract the feature table for a folder of feathers, and fit a model.\n"
                    "Starts a local server and opens http://127.0.0.1:7432/visual. Ctrl+C stops it.",
    )
    visual_parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Start the server without opening a browser tab",
    )
    visual_parser.add_argument("--host", default=None, help="Bind address (default 127.0.0.1)")
    visual_parser.add_argument("--port", type=int, default=None, help="Port (default 7432)")
    theta_parser = subparsers.add_parser(
        "feature_explorer", aliases=["theta"],
        help="Open the feature explorer (Sterimol, cone angle, %%Vbur, dipoles, charges) in a browser", formatter_class=raw,
        description="The feature explorer (until October 2026: the Sterimol theta explorer). Draws B1, B5, L, loc B5, theta, the cone\n"
                    "angle, %Vbur, dipole frames, charges, lengths and angles on the molecule, stacks aligned structures, and exports\n"
                    "figures. Same server as visual_extractor, at /feature_explorer. Manual: " + DOCS + "FEATURE_EXPLORER.md",
    )
    theta_parser.add_argument("--mol", default=None, help="Open this built-in structure, e.g. FL_lig_13")
    theta_parser.add_argument("--preset", default=None,
                              help="Apply this preset (a number, or the start of its name), e.g. 'theta construction'")
    theta_parser.add_argument("--no-browser", action="store_true", help="Start the server without opening a browser tab")
    theta_parser.add_argument("--host", default=None, help="Bind address (default 127.0.0.1)")
    theta_parser.add_argument("--port", type=int, default=None, help="Port (default 7432)")
    model_parser = subparsers.add_parser("model", help="Run regression or classification",
        description="Search linear (regression) or logistic (classification) models over every combination of min..max features of a\n"
                    "feature table, keep the best, and write the report. Example:\n"
                    "  descripytor model -m regression -f features.csv -t targets.csv --min-features 2 --max-features 3", formatter_class=raw)
    feature_extraction = subparsers.add_parser("extractor", help="Run feature extraction - complete set - from input file",
        description="Extract every descriptor named in an input file (the atom choices; visual_extractor writes one) for all the\n"
                    ".feather files of a folder. Example:\n  descripytor extractor -i input.json -o features -f feathers/", formatter_class=raw)
    conver_parser = subparsers.add_parser("logs_to_feather", help="Convert log files to feather files",
        description="Convert a folder of Gaussian .log files to .feather files. Takes no options: it asks for the folder.")
    cube_parser = subparsers.add_parser("cube", help="Calculates cube sterimol from cube files",
        description="Sterimol from Gaussian cube files, written to cube_sterimol.csv. Takes no options: it asks for the folder and the atoms.")
    sterimol_parser = subparsers.add_parser("sterimol", help="Calculate sterimol values from xyz files",
        description="Sterimol for a folder of .xyz files. Takes no options: it asks for the folder, the atoms and the radii.")
    bassa_parser = subparsers.add_parser(
        "bassa_explorer", aliases=["bassa"],
        help="Open the BASSA model explorer: the sampler's models beside the molecule", formatter_class=raw,
        description="The BASSA model explorer: columns are the models the spike-and-slab sampler visited, rows are descriptors, and the\n"
                    "feature explorer draws the selected descriptors on a compound. Opens http://127.0.0.1:8765/ with a bundled example.\n"
                    "Running a new table needs BASSA's sampler in its own environment: run  descripytor bassa_setup  once.\n"
                    "Manual: " + DOCS + "BASSA_EXPLORER.md",
    )
    bassa_parser.add_argument("--cases", default=None, help="Folder where kept cases live (default: ./bassa_cases if it exists, else the last --cases, else ~/.descripytor/cases)")
    bassa_parser.add_argument("--bassa-python", default=None,
                              help="The Python that has bassa_reg, for new tables; remembered. Default: DESCRIPYTOR_BASSA_PYTHON, "
                                   "else the one remembered (by bassa_setup or an earlier --bassa-python), else this one")
    setup_parser = subparsers.add_parser(
        "bassa_setup", help="Once: make the environment BASSA's sampler runs in, for new tables in bassa_explorer", formatter_class=raw,
        description="bassa-reg pins NumPy 2 while DescriPyTor needs NumPy < 2, so the sampler gets its own virtual environment\n"
                    "(default ~/.descripytor/bassa_env). This makes it, installs bassa-reg, checks the sampler imports and remembers\n"
                    "its Python for bassa_explorer. Takes a few minutes.")
    setup_parser.add_argument("--dir", default=None, help="Where to make the environment (default ~/.descripytor/bassa_env)")
    setup_parser.add_argument("--python", default=None, help="The Python to make it with (default: this one)")
    bassa_parser.add_argument("--no-browser", action="store_true", help="Start the server without opening a browser tab")
    bassa_parser.add_argument("--host", default=None, help="Bind address (default 127.0.0.1)")
    bassa_parser.add_argument("--port", type=int, default=None, help="Port (default 8765, or BASSA_EXPLORER_PORT)")
    renumber_parser = subparsers.add_parser(
        "renumber", help="Number every molecule of a folder like a reference (Conformaty)", formatter_class=raw,
        description="Number the atoms of every .feather in a folder like a reference molecule, so a descriptor picked on\n"
                    "the reference reads the same atom in every molecule. Conformaty matches the common substructure,\n"
                    "fits each molecule onto the reference and assigns the remaining atoms by position (the two H of a\n"
                    "CH2 by where they sit). The renumbered copies and renumbering.json go to a new folder; the\n"
                    "originals are not changed. Needs Conformaty: pip install git+https://github.com/Milo-group/Conformaty\n"
                    "Example:\n  descripytor renumber feathers/ -r basic --anchor m_Br:9,9")
    renumber_parser.add_argument("folder", help="Folder of .feather files")
    renumber_parser.add_argument("-r", "--reference", help="Molecule to number like (its name; default: the smallest)")
    renumber_parser.add_argument("-o", "--output", help="Output folder (default: <folder>/renumbered)")
    renumber_parser.add_argument("--anchor", action="append", metavar="NAME:REF,ATOM",
                                 help="Atom ATOM of molecule NAME must be reference atom REF (1-based); repeatable")
    renumber_parser.add_argument("--stages", nargs="+", choices=["complete_rings", "extend_mcs", "coords_proximity_any"],
                                 help="Conformaty's middle stages, in order (default: all three)")
    renumber_parser.add_argument("--match-charge", action="store_true", help="Match formal charges and radicals in the MCS")
    renumber_parser.add_argument("--fit", default="rmsd", choices=["rmsd", "msd", "mixed", "mad", "sig", "kabsch"],
                                 help="Rigid-body fit of the reported overlay and RMSD (not the numbering): Conformaty's "
                                      "losses (mixed, its own default, favours atoms in near-exact overlap but can stall "
                                      "far away) or kabsch (closed-form least squares). Default rmsd")
    renumber_parser.add_argument("--sigma", type=float, default=30, help="Sigma of the mixed and sig losses (default 30)")
    renumber_parser.add_argument("--fit-atoms", type=int, nargs="+", metavar="ATOM",
                                 help="Reference atoms (1-based) the reported fit uses (default: the heavy atoms with the same element)")
    renumber_parser.add_argument("--workers", type=int, help="Molecules at a time in separate processes (default 1)")
    pipeline_parser = subparsers.add_parser(
        "pipeline", help="SMILES -> structures -> chained cluster stages (GOAT/UMA/xTB/ORCA/Gaussian)", formatter_class=raw,
        description="Run a protocol on the cluster, one action at a time:\n"
                    "  build    SMILES -> numbered xyz          submit   send the stages to the queue\n"
                    "  status   what is running, done, failed   watch    status every --interval seconds\n"
                    "  retry    resubmit failures               reset    forget a run\n"
                    "  fetch    bring the outputs back          adopt    take over calculations that already exist\n"
                    "Manual: docs/PIPELINE.md.")
    design_parser = subparsers.add_parser(
        "design", help="Does a model find the best compounds; rank candidates for the next round", formatter_class=raw,
        description="OLS on the chosen columns of a table (first column: compound names; one column the response).\n"
                    "Reports, leave-one-out: how many of the best quarter are predicted in the best quarter (hypergeometric p),\n"
                    "how far the best are under-predicted, and a leave-the-best-out test against shuffled responses.\n"
                    "With --candidates: ranks new compounds (same columns) by expected improvement, with a prediction\n"
                    "interval and a flag for compounds outside the measured domain (picked last unless --allow-outside).\n"
                    "Read the predictions for rank, not size: the best are usually under-predicted.")
    design_parser.add_argument("table", help="CSV: names, the response, descriptor columns")
    design_parser.add_argument("-y", required=True, help="Response column")
    design_parser.add_argument("--columns", nargs="+", help="Descriptor columns of the model (default: all but -y)")
    design_parser.add_argument("--lower", action="store_true", help="Lower is better (e.g. a ΔΔG‡ that favours the other product)")
    design_parser.add_argument("--k", type=int, help="How many count as the best (default: a quarter, rounded up)")
    design_parser.add_argument("--null", type=int, default=20, help="Shuffles for the leave-the-best-out null (default 20)")
    design_parser.add_argument("--candidates", help="CSV of new compounds to rank (names, the same columns)")
    design_parser.add_argument("--batch", type=int, default=1, help="Candidates to pick (default 1)")
    design_parser.add_argument("--allow-outside", action="store_true", help="Rank by expected improvement alone, inside the domain or not")
    design_parser.add_argument("-o", "--output", help="Ranked candidates CSV (default: <candidates>_ranked.csv)")
    design_parser.add_argument("--json", help="Write the report to this JSON file")
    joint_parser = subparsers.add_parser(
        "joint", help="Choose descriptors on related reactions together", formatter_class=raw,
        description="Separate against joint descriptor selection over reactions run on the same compounds, leave-one-compound-out\n"
                    "with the selection redone in every fold. Joint: one set of 1 to --max-terms columns for all reactions, each\n"
                    "reaction fitted with its own coefficients (a descriptor may change sign between reactions). Responses are\n"
                    "never averaged over reactions. Passes when the mean gain is >= 0.05, p < 0.05 and above shuffled partners.\n"
                    "Keep the column pool small and chosen before looking at the responses: every pair is a candidate.")
    joint_parser.add_argument("features", help="CSV: compound names, descriptor columns")
    joint_parser.add_argument("outcomes", help="CSV, one row per compound and reaction (long format)")
    joint_parser.add_argument("--name-col", default="name", help="Compound column of the outcomes (default name)")
    joint_parser.add_argument("--reaction-col", default="reaction", help="Reaction column of the outcomes (default reaction)")
    joint_parser.add_argument("--y-col", default="y", help="Response column of the outcomes (default y)")
    joint_parser.add_argument("--columns", nargs="+", help="Descriptor pool (default: every column of the features)")
    joint_parser.add_argument("--max-terms", type=int, default=2, choices=[1, 2, 3], help="Columns per model (default 2)")
    joint_parser.add_argument("--null", type=int, default=200, help="Shuffles of the partner reactions (default 200)")
    joint_parser.add_argument("-o", "--output", help="Write the per-reaction coefficients to this CSV")
    forms_parser = subparsers.add_parser(
        "forms", help="Logistic link, flexible ceiling, threshold analysis", formatter_class=raw,
        description="Other model forms, each against OLS on the same compounds, leave-one-out:\n"
                    "  --logit           binomial GLM on y/scale (yields): passes only if it beats OLS and clipped OLS\n"
                    "  --ceiling         a Gaussian process and a random forest on all the columns: is linear enough?\n"
                    "  --threshold VALUE one column and one cut separating y >= VALUE (or <= with --lower), chosen in every fold")
    forms_parser.add_argument("table", help="CSV: names, the response, descriptor columns")
    forms_parser.add_argument("-y", required=True, help="Response column")
    forms_parser.add_argument("--columns", nargs="+", help="Descriptor columns (default: all but -y)")
    forms_parser.add_argument("--logit", action="store_true", help="Logistic link against OLS")
    forms_parser.add_argument("--scale", type=float, default=100.0, help="Upper bound of the response for --logit (default 100)")
    forms_parser.add_argument("--ceiling", action="store_true", help="Gaussian process and random forest against OLS")
    forms_parser.add_argument("--threshold", type=float, help="Response value that separates good from bad")
    forms_parser.add_argument("--lower", action="store_true", help="With --threshold: lower is good")
    forms_parser.add_argument("--null", type=int, default=200, help="Permutations for the nulls (default 200)")
    audit_parser = subparsers.add_parser(
        "audit", help="Columns one compound dominates; cut a wide table by correlation clusters", formatter_class=raw,
        description="Checks on a descriptor table that never read the response:\n"
                    "  z-audit    for each column, how far the most extreme compound sits from the others (in their SD). A column\n"
                    "             beyond --zcut whose SD collapses without that compound (ratio < 0.25) is degenerate: drop it\n"
                    "             before any search (the cc sub_B1_asym retraction). Otherwise the compound is a real outlier: look at it.\n"
                    "  --reduce N clean (incomplete, constant, < 4 levels, duplicate), cluster on 1 - |r|, keep one medoid per cluster")
    audit_parser.add_argument("table", help="CSV: names, descriptor columns (and any response columns, named with -y)")
    audit_parser.add_argument("-y", nargs="+", help="Response or label columns to leave out of the audit (kept in --reduce's table)")
    audit_parser.add_argument("--zcut", type=float, default=6.0, help="|z| that flags a column (default 6)")
    audit_parser.add_argument("--show", type=int, default=20, help="Flagged columns to print (default 20)")
    audit_parser.add_argument("--reduce", type=int, help="Cut to this many columns, one per correlation cluster")
    audit_parser.add_argument("--min-levels", type=int, default=4, help="--reduce: drop columns with fewer distinct values (default 4)")
    audit_parser.add_argument("--drop", help="--reduce: regex of column names to drop first (e.g. 'Stretch|^energy$')")
    audit_parser.add_argument("-o", "--output", help="Prefix of the written files (default: the table's name)")
    protocol_parser = subparsers.add_parser(
        "bassa_protocol", help="Validation of a kept BASSA case: deals of 5-fold + LOO, chance, seed, drop-one", formatter_class=raw,
        description="The validation protocol of a case kept by bassa_explorer, with the sampler rerun in every fold:\n"
                    "  protocol   20 deals of stratified 5-fold and leave-one-out: Q2, MAE, 90 % interval coverage -> protocol_<case>.json\n"
                    "  chance     the response shuffled, one deal each: the chance line and its p -> chance_<case>.json\n"
                    "  seed       the deals again with another sampler seed -> seed_<case>.json\n"
                    "  drop_one   BASSA refitted without each compound: which compounds carry the model -> dropone_<case>.json\n"
                    "Runs under the sampler's Python (as bassa_explorer). Each sampler run takes seconds; protocol is ~100 + n runs.\n"
                    "Manual: " + DOCS + "BASSA_EXPLORER.md")
    protocol_parser.add_argument("check", choices=["protocol", "chance", "seed", "drop_one"])
    protocol_parser.add_argument("cases", nargs="+", help="Case names in the cases folder, or paths to their .js files")
    protocol_parser.add_argument("--cases-dir", help="The cases folder (default: the one bassa_explorer uses)")
    protocol_parser.add_argument("--deals", type=int, default=20, help="protocol: deals of 5-fold (default 20)")
    protocol_parser.add_argument("--no-loo", action="store_true", help="protocol: skip leave-one-out")
    protocol_parser.add_argument("--shuffles", type=int, default=20, help="chance: shuffles (default 20)")
    protocol_parser.add_argument("--seed", type=int, default=11, help="seed: the other sampler seed (default 11)")
    protocol_parser.add_argument("--dry", action="store_true", help="Write nothing; print only")
    protocol_parser.add_argument("--bassa-python", default=None, help="The Python that has bassa_reg (default: as bassa_explorer)")
    pipeline_parser.add_argument("action", choices=["build", "adopt", "submit", "status", "retry", "reset", "watch", "fetch"])
    pipeline_parser.add_argument("protocol", help="Protocol JSON (see docs/PIPELINE.md)")
    pipeline_parser.add_argument("--dry-run", action="store_true", help="submit/retry: show what would be sent, send nothing")
    pipeline_parser.add_argument("--force", action="store_true", help="submit: overwrite a run already on the cluster")
    pipeline_parser.add_argument("--keep-going", action="store_true", help="build: skip molecules that fail to build")
    pipeline_parser.add_argument("--interval", type=int, default=600, help="watch: seconds between reports")
    pipeline_parser.add_argument("--all-stages", action="store_true", help="fetch: every stage's outputs, not only the last")
    pipeline_parser.add_argument("--sources", help="adopt: CSV of id, remote_source, local_start")
    pipeline_parser.add_argument("--stage", type=int, default=1, help="adopt: the stage the sources complete (1-based)")
    pipeline_parser.add_argument("--from-stage", type=int, default=1, help="submit: first stage to submit (1-based)")
    pipeline_parser.add_argument("--only", default=None, help="submit: comma-separated ids to submit")
    pipeline_parser.add_argument("--only-adopted", action="store_true", help="submit: only the ids in adopted.txt")
    pipeline_parser.add_argument("--all", action="store_true", help="retry: also failures the same input would repeat (structure, imaginary)")

    
    model_parser.add_argument("-m", "--mode", choices=["regression", "classification"], required=True,
                              help="Which task to run.")
    model_parser.add_argument("-f", "--features_csv", required=True, help="Path to features CSV.")
    model_parser.add_argument("-t", "--target_csv", default="", help="Path to target/labels CSV.")
    model_parser.add_argument("-y", "--y_value", default="output", help="Base name for output files.")
    model_parser.add_argument("-j", "--n_jobs", type=int,
                              help="Cores to use (-1 = all). If omitted, uses NSLOTS or all logical cores. Values > 1 enable parallel execution.")
    model_parser.add_argument("--min-features", type=int, default=4, help="Minimum features per model.")
    model_parser.add_argument("--max-features", type=int, default=4, help="Maximum features per model.")
    model_parser.add_argument("--top-n", type=int, default=50, help="How many top models to keep/evaluate.")
    model_parser.add_argument("--bool-parallel", action="store_true", help="Enable parallel evaluation.")
    model_parser.add_argument("--threshold", type=float, default=0.70, help="Initial threshold (regression(R2)/classification(mcfadden_R2)).")
    model_parser.add_argument("--leave-out", nargs="*", default=[], help="Space-separated list of samples to leave out.")
    feature_extraction.add_argument("-i", "--input", required=True, help="Path to input file.")
    feature_extraction.add_argument("-o", "--output", required=True, help="Path to output file.")
    feature_extraction.add_argument("-f", "--feather_directory", default=".", help="Directory of feather files set to extract features from.")

    args = parser.parse_args()
    if args.command == "bassa_setup":
        return run_bassa_setup(args.dir, args.python)
    if args.command in ("bassa_explorer", "bassa"):
        return run_bassa_explorer(cases=args.cases, host=args.host, port=args.port,
                                  open_browser=not args.no_browser, bassa_python=args.bassa_python)
    if args.command == "pipeline":
        return run_pipeline(args)
    if args.command == "bassa_protocol":
        return run_bassa_protocol(args)
    if args.command in ("design", "joint", "forms", "audit"):
        from M3_modeler import commands
        return getattr(commands, f"run_{args.command}")(args)
    if args.command == "renumber":
        return run_renumber(args)
    if args.command in ("feature_explorer", "theta"):
        run_visual_app(host=args.host, port=args.port, open_browser=not args.no_browser,
                       page=feature_explorer_page(args.mol, args.preset))
    elif args.command in ("visual_extractor", "visual"):
        run_visual_app(
            host=args.host,
            port=args.port,
            open_browser=not args.no_browser,
        )
    elif args.command == "logs_to_feather":
        _load_runtime_dependencies()
        log_dir = input("Enter the path to the log files directory: ")
        logs_to_feather(log_dir)
        

    elif args.command == "cube":
        _load_runtime_dependencies()
        while True:
            cube_file_path = input("Enter the full path to directory with cube files: ")
            base_atoms = input("Enter the atom indices: ")
            base_atoms=convert_to_list_or_nested_list(base_atoms)
            sterimol=cube_many(cube_file_path, base_atoms)
            ## save to csv
            # df=pd.DataFrame(sterimol.sterimol_df)
            sterimol.sterimol_df.to_csv('cube_sterimol.csv')

            another_dir = input("Do you want to input another directory? (yes/no): ").strip().lower()
            if another_dir != 'yes':
                break
    elif args.command == "sterimol":
        _load_runtime_dependencies()
        while True:
            xyz_dir = input("Enter the path to the xyz files directory: ")
            base_atoms = input("Enter the atom indices: ")
            params = inspect.signature(Molecules_xyz.get_sterimol_df)
            radii=list(params.parameters.values())[2].default
            radius_input = input(f"Enter the radii (default: {radii}): ")
            if radius_input=='':
                radius_input=radii
            base_atoms=convert_to_list_or_nested_list(base_atoms)
            sterimol=Molecules_xyz(xyz_dir)
            sterimol_df=sterimol.get_sterimol_df(base_atoms, radii=radius_input)
            os.chdir(xyz_dir)
            sterimol_df.to_csv('xyz_sterimol.csv')
            print(sterimol_df.head())
            print(f"Sterimol values saved to xyz_sterimol.csv in {xyz_dir}")

            another_dir = input("Do you want to input another directory? (yes/no): ").strip().lower()
            if another_dir != 'yes':
                break

            
    elif args.command == "model":
        # Prepare environment/import path & working dir
        n_jobs = resolve_n_jobs(args.n_jobs)
        if n_jobs == -1:
            import psutil
            n_jobs = int(psutil.cpu_count(logical=True))
        add_modeling_path()
        # safe_chdir(args.workdir, args.features_csv)

        elapsed = run_experiment(
            mode=args.mode,
            features_csv=args.features_csv,
            target_csv=args.target_csv,
            y_value=args.y_value,
            n_jobs=n_jobs,
            min_features_num=args.min_features,
            max_features_num=args.max_features,
            threshold=args.threshold,
            top_n=args.top_n,
            leave_out=args.leave_out,
            bool_parallel=args.bool_parallel,
        )
        print(f"[model] Total execution time: {elapsed:.2f} seconds")

    elif args.command == "extractor":
        # Run the feature extraction
        run_feature_extraction(
            input_file=args.input,
            molecules_dir_name=args.feather_directory,
            output_file=args.output
        )


    else:
        parser.print_help()

if __name__ == "__main__":
    sys.exit(main())     # as the installed command does: a refused start exits with 1


