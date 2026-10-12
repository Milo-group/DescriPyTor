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
        return "0.3.3"

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

Batch commands:
  extractor          descriptors for a folder of .feather files, from an input file of atom choices
  model              search regression or classification models over a feature table
  pipeline           SMILES -> structures -> chained cluster stages (GOAT / UMA / xTB / ORCA / Gaussian)

Converters (they ask for their inputs):
  logs_to_feather    Gaussian .log files -> .feather
  sterimol           Sterimol from .xyz files
  cube               Sterimol from Gaussian cube files"""
CLI_EPILOG = """examples:
  descripytor visual_extractor
  descripytor feature_explorer --mol FL_lig_13 --preset "theta construction"
  descripytor bassa_explorer --cases my_cases --bassa-python C:/envs/bassa/python.exe
  descripytor extractor -i input.json -o features -f feathers/
  descripytor model -m regression -f features.csv -t targets.csv --min-features 2 --max-features 3
  descripytor pipeline status protocol.json

`descripytor <command> --help` lists the options of one command. Older names still work:
visual = visual_extractor, theta = feature_explorer, bassa = bassa_explorer.
Manuals: docs/WORKFLOW.md, docs/FEATURE_EXPLORER.md, docs/BASSA_EXPLORER.md, docs/PIPELINE.md."""


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
    host, port = host or "127.0.0.1", int(port or os.environ.get("BASSA_EXPLORER_PORT", 8765))
    env["BASSA_EXPLORER_PORT"] = str(port)                   # so the other pages' links find this server
    url = f"http://{'127.0.0.1' if host in ('0.0.0.0', '') else host}:{port}/"
    print(f"\n  BASSA model explorer\n  Open:  {url}\n  Press Ctrl+C to stop\n", flush=True)
    python = bassa_python or os.environ.get("DESCRIPYTOR_BASSA_PYTHON") or sys.executable
    _ensure_project_on_path()
    from M3_modeler.bassa_explorer import running
    already = running.state(port)
    if already == "current":                                 # this version's server is up: use it
        print("  It is already running; this process will not start a second server." + chr(10), flush=True)
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
        only = args.only.split(",") if args.only else None
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
    parser = argparse.ArgumentParser(prog="descripytor", formatter_class=raw, description=CLI_HELP, epilog=CLI_EPILOG)
    subparsers = parser.add_subparsers(dest="command", metavar="<command>", title="commands")

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
                    "figures. Same server as visual_extractor, at /feature_explorer. Manual: docs/FEATURE_EXPLORER.md.",
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
                    "Running a new table needs bassa_reg (pip install bassa-reg) in the Python given by --bassa-python.\n"
                    "Manual: docs/BASSA_EXPLORER.md.",
    )
    bassa_parser.add_argument("--cases", default=None, help="Folder where kept cases live (default ./bassa_cases)")
    bassa_parser.add_argument("--bassa-python", default=None,
                              help="The Python that has bassa_reg, for running new tables (default: this one, or DESCRIPYTOR_BASSA_PYTHON)")
    bassa_parser.add_argument("--no-browser", action="store_true", help="Start the server without opening a browser tab")
    bassa_parser.add_argument("--host", default=None, help="Bind address (default 127.0.0.1)")
    bassa_parser.add_argument("--port", type=int, default=None, help="Port (default 8765, or BASSA_EXPLORER_PORT)")
    pipeline_parser = subparsers.add_parser(
        "pipeline", help="SMILES -> structures -> chained cluster stages (GOAT/UMA/xTB/ORCA/Gaussian)", formatter_class=raw,
        description="Run a protocol on the cluster, one action at a time:\n"
                    "  build    SMILES -> numbered xyz          submit   send the stages to the queue\n"
                    "  status   what is running, done, failed   watch    status every --interval seconds\n"
                    "  retry    resubmit failures               reset    forget a run\n"
                    "  fetch    bring the outputs back          adopt    take over calculations that already exist\n"
                    "Manual: docs/PIPELINE.md.")
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
    if args.command in ("bassa_explorer", "bassa"):
        return run_bassa_explorer(cases=args.cases, host=args.host, port=args.port,
                                  open_browser=not args.no_browser, bassa_python=args.bassa_python)
    if args.command == "pipeline":
        return run_pipeline(args)
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


