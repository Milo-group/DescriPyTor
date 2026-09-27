"""
Crash-resistant, resumable Boltzmann-averaged descriptor extraction.

Extracts the same features as the original xyz_best script:

Sterimol:
    [8,1], [1,8], [10,8], [9,6], [4,1], [3,1], [3,5]

Angles:
    [2,1,5], [5,18,6]

Dihedrals:
    [5,1,8,6], [2,3,4,20]

Bond lengths:
    [1,5], [5,7], [18,5], [18,6], [1,18]

Buried volume:
    metal atom 18, sphere radius 3.5 Å

Features:
    1. Uses every conformer retained by CREST.
    2. Boltzmann-averages descriptors at the selected temperature.
    3. Pools specified rerun folders before Boltzmann averaging.
    4. Continues if an individual descriptor fails.
    5. Continues if an individual ligand fails.
    6. Saves a checkpoint after every ligand.
    7. Automatically resumes from the checkpoint.
    8. Reports descriptor time, ligand time, elapsed time, and ETA.
"""

# =============================================================================
# Imports
# =============================================================================

import os
import re
import sys
import time
import traceback
from datetime import datetime

import numpy as np
import pandas as pd


# =============================================================================
# Import local DescripyTor module
# =============================================================================

ROOT_DIR = (
    r"C:\Users\edens\Desktop\DescriPytor"
    r"\DescriPyTor-main\MolFeatures"
)

EXTRACTOR_DIR = os.path.join(
    ROOT_DIR,
    "M2_data_extractor",
)

if EXTRACTOR_DIR not in sys.path:
    sys.path.append(EXTRACTOR_DIR)

os.chdir(EXTRACTOR_DIR)

from sterimol_standalone import (
    ConformerEnsemble,
    CrestLigandSet,
)


# =============================================================================
# Input configuration
# =============================================================================

conformers_dir = (
    r"C:\Users\edens\Desktop\Lab\datasets"
    r"\doyle_conformers"
)

# All atom indices are 1-indexed
sterimol_pairs = [
    [8, 1],
    [1, 8],
    [10, 8],
    [9, 6],
    [4, 1],
    [3, 1],
    [3, 5],
]

angle_groups = [
    [2, 1, 5],
    [5, 18, 6],
]

dihedral_groups = [
    [5, 1, 8, 6],
    [2, 3, 4, 20],
]

bond_length_pairs = [
    [1, 5],
    [5, 7],
    [18, 5],
    [18, 6],
    [1, 18],
]

metal_index = 18
buried_volume_radius = 3.5

radii = "CPK"
temperature = 298.15

# Optional conformer filters
energy_cutoff_kcal = None
min_weight = None


# =============================================================================
# Output configuration
# =============================================================================

save_csv = True

output_dir = os.path.join(
    conformers_dir,
    "boltzmann_descriptor_results",
)

output_csv = (
    "doyle_boltzmann_averaged_features.csv"
)

detail_csv = output_csv.replace(
    ".csv",
    "_per_conformer.csv",
)

warnings_txt = output_csv.replace(
    ".csv",
    "_warnings.txt",
)

timing_csv = output_csv.replace(
    ".csv",
    "_timing.csv",
)


# =============================================================================
# Checkpoint configuration
# =============================================================================

resume_from_checkpoint = True
save_checkpoint_after_each_ligand = True

checkpoint_dir = os.path.join(
    output_dir,
    "checkpoints",
)

checkpoint_summary_csv = os.path.join(
    checkpoint_dir,
    "summary_checkpoint.csv",
)

checkpoint_detail_csv = os.path.join(
    checkpoint_dir,
    "per_conformer_checkpoint.csv",
)

checkpoint_timing_csv = os.path.join(
    checkpoint_dir,
    "timing_checkpoint.csv",
)

checkpoint_warnings_txt = os.path.join(
    checkpoint_dir,
    "warnings_checkpoint.txt",
)


# =============================================================================
# CREST folders to pool
# =============================================================================

MERGE_GROUPS = {
    "L12": ["L12", "L12_T"],
    "L20": ["L20", "L20b"],
    "L22": ["L22", "L22b"],
    "L23": ["L23", "L23b"],
    "L24": ["L24", "L24b"],
}


# =============================================================================
# General helper functions
# =============================================================================

def natural_key(value):
    """Natural sorting, so L2 appears before L10."""
    return [
        int(token) if token.isdigit()
        else token.lower()
        for token in re.split(
            r"(\d+)",
            str(value),
        )
    ]


def format_duration(seconds):
    """Convert seconds into a readable time string."""
    if seconds is None or not np.isfinite(seconds):
        return "unknown"

    seconds = max(0, int(round(seconds)))

    days, remainder = divmod(seconds, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, seconds = divmod(remainder, 60)

    parts = []

    if days:
        parts.append(f"{days}d")

    if hours or days:
        parts.append(f"{hours:02d}h")

    if minutes or hours or days:
        parts.append(f"{minutes:02d}m")

    parts.append(f"{seconds:02d}s")

    return " ".join(parts)


def current_timestamp():
    """Return the current local timestamp."""
    return datetime.now().strftime(
        "%Y-%m-%d %H:%M:%S"
    )


def format_exception(
    ligand_name,
    descriptor_name,
    exception,
):
    """Create a detailed warning with traceback."""
    return (
        "\n"
        + "=" * 80
        + f"\nTime: {current_timestamp()}"
        + f"\nLigand: {ligand_name}"
        + f"\nDescriptor: {descriptor_name}"
        + f"\nError type: "
        + f"{type(exception).__name__}"
        + f"\nError: {exception}"
        + "\n\nTraceback:\n"
        + traceback.format_exc()
    )


# =============================================================================
# Atomic saving functions
# =============================================================================

def atomic_save_csv(
    dataframe,
    path,
    index=True,
):
    """
    Save to a temporary file first, then replace the previous file.

    This helps prevent a partially written CSV if the process is interrupted.
    """
    directory = os.path.dirname(path)

    if directory:
        os.makedirs(
            directory,
            exist_ok=True,
        )

    temporary_path = path + ".tmp"

    dataframe.to_csv(
        temporary_path,
        index=index,
    )

    os.replace(
        temporary_path,
        path,
    )


def atomic_save_text(
    text,
    path,
):
    """Atomically save a text file."""
    directory = os.path.dirname(path)

    if directory:
        os.makedirs(
            directory,
            exist_ok=True,
        )

    temporary_path = path + ".tmp"

    with open(
        temporary_path,
        "w",
        encoding="utf-8",
    ) as output_file:
        output_file.write(text)

    os.replace(
        temporary_path,
        path,
    )


def save_warning_log(
    warnings,
    path,
):
    """Save accumulated warnings."""
    text = (
        "\n".join(map(str, warnings))
        if warnings
        else "No warnings."
    )

    atomic_save_text(
        text,
        path,
    )


# =============================================================================
# Descriptor result helpers
# =============================================================================

def split_method_result(result):
    """
    Validate that an ensemble method returned:

        per_conformer_dataframe, boltzmann_summary
    """
    if not isinstance(result, tuple):
        raise TypeError(
            "Expected a tuple containing "
            "(per_conformer, boltzmann_summary), "
            f"but received {type(result).__name__}."
        )

    if len(result) != 2:
        raise ValueError(
            "Expected exactly two returned values, "
            f"but received {len(result)}."
        )

    return result[0], result[1]


def scalar_from_result(
    result,
    preferred_columns=None,
):
    """
    Extract a scalar from a DataFrame, Series, dictionary, array, or scalar.
    """
    if preferred_columns is None:
        preferred_columns = []

    metadata_columns = {
        "n_conformers",
        "n_conformers_used",
        "temperature",
        "temperature_K",
        "weight_sum",
    }

    if isinstance(result, pd.DataFrame):

        for column in preferred_columns:
            if column in result.columns:
                return result[column].iloc[0]

        candidate_columns = [
            column
            for column in result.columns
            if column not in metadata_columns
        ]

        if len(candidate_columns) == 1:
            return result[
                candidate_columns[0]
            ].iloc[0]

        raise ValueError(
            "Could not determine the value column. "
            f"Available columns: {result.columns.tolist()}"
        )

    if isinstance(result, pd.Series):

        for key in preferred_columns:
            if key in result.index:
                return result[key]

        candidate_keys = [
            key
            for key in result.index
            if key not in metadata_columns
        ]

        if len(candidate_keys) == 1:
            return result[candidate_keys[0]]

        if len(result) == 1:
            return result.iloc[0]

        raise ValueError(
            "Could not determine the value from Series. "
            f"Available keys: {result.index.tolist()}"
        )

    if isinstance(result, dict):

        for key in preferred_columns:
            if key in result:
                return result[key]

        candidate_keys = [
            key
            for key in result
            if key not in metadata_columns
        ]

        if len(candidate_keys) == 1:
            return result[candidate_keys[0]]

        if len(result) == 1:
            return next(iter(result.values()))

        raise ValueError(
            "Could not determine the value from dictionary. "
            f"Available keys: {list(result)}"
        )

    if isinstance(result, np.ndarray):

        if result.size == 1:
            return result.item()

        raise ValueError(
            "Expected one value, but NumPy array "
            f"contains {result.size} values."
        )

    return result


def extract_metadata(
    boltzmann_result,
    row,
    pair_label=None,
):
    """Extract conformer-count metadata when available."""
    if isinstance(boltzmann_result, pd.DataFrame):

        if "n_conformers" in boltzmann_result.columns:
            row["n_conformers"] = int(
                boltzmann_result[
                    "n_conformers"
                ].iloc[0]
            )

        if (
            pair_label is not None
            and "n_conformers_used"
            in boltzmann_result.columns
        ):
            row[
                f"n_conformers_used_{pair_label}"
            ] = int(
                boltzmann_result[
                    "n_conformers_used"
                ].iloc[0]
            )

    elif isinstance(boltzmann_result, pd.Series):

        if "n_conformers" in boltzmann_result.index:
            row["n_conformers"] = int(
                boltzmann_result[
                    "n_conformers"
                ]
            )

        if (
            pair_label is not None
            and "n_conformers_used"
            in boltzmann_result.index
        ):
            row[
                f"n_conformers_used_{pair_label}"
            ] = int(
                boltzmann_result[
                    "n_conformers_used"
                ]
            )

    elif isinstance(boltzmann_result, dict):

        if "n_conformers" in boltzmann_result:
            row["n_conformers"] = int(
                boltzmann_result[
                    "n_conformers"
                ]
            )

        if (
            pair_label is not None
            and "n_conformers_used"
            in boltzmann_result
        ):
            row[
                f"n_conformers_used_{pair_label}"
            ] = int(
                boltzmann_result[
                    "n_conformers_used"
                ]
            )


# =============================================================================
# Timing and safe-execution helpers
# =============================================================================

def add_timing_record(
    timing_records,
    ligand_name,
    descriptor_type,
    descriptor_name,
    status,
    seconds,
    error=None,
):
    """Add one descriptor timing record."""
    timing_records.append(
        {
            "ligand": ligand_name,
            "descriptor_type": descriptor_type,
            "descriptor": descriptor_name,
            "status": status,
            "seconds": float(seconds),
            "formatted_time": format_duration(seconds),
            "error_type": (
                type(error).__name__
                if error is not None
                else ""
            ),
            "error_message": (
                str(error)
                if error is not None
                else ""
            ),
            "timestamp": current_timestamp(),
        }
    )


def safe_descriptor_call(
    ligand_name,
    descriptor_type,
    descriptor_name,
    calculation,
    run_warnings,
    timing_records,
):
    """
    Run one descriptor calculation safely and measure its execution time.
    """
    descriptor_start = time.perf_counter()

    try:
        result = calculation()

        per_conf, boltzmann_summary = (
            split_method_result(result)
        )

        descriptor_seconds = (
            time.perf_counter()
            - descriptor_start
        )

        add_timing_record(
            timing_records=timing_records,
            ligand_name=ligand_name,
            descriptor_type=descriptor_type,
            descriptor_name=descriptor_name,
            status="success",
            seconds=descriptor_seconds,
        )

        print(
            f"    {descriptor_name}: "
            f"completed in "
            f"{format_duration(descriptor_seconds)}"
        )

        return (
            True,
            per_conf,
            boltzmann_summary,
        )

    except Exception as error:

        descriptor_seconds = (
            time.perf_counter()
            - descriptor_start
        )

        run_warnings.append(
            format_exception(
                ligand_name,
                descriptor_name,
                error,
            )
        )

        add_timing_record(
            timing_records=timing_records,
            ligand_name=ligand_name,
            descriptor_type=descriptor_type,
            descriptor_name=descriptor_name,
            status="failed",
            seconds=descriptor_seconds,
            error=error,
        )

        print(
            f"    {descriptor_name}: FAILED after "
            f"{format_duration(descriptor_seconds)} | "
            f"{type(error).__name__}: {error}"
        )

        return False, None, None


# =============================================================================
# Checkpoint helper
# =============================================================================

def save_checkpoint(
    summary_rows,
    detail_frames,
    timing_records,
    run_warnings,
):
    """Save all work completed so far."""
    checkpoint_start = time.perf_counter()

    summary_checkpoint = pd.DataFrame.from_dict(
        summary_rows,
        orient="index",
    )

    summary_checkpoint.index.name = "ligand"

    if not summary_checkpoint.empty:
        summary_checkpoint = summary_checkpoint.reindex(
            sorted(
                summary_checkpoint.index,
                key=natural_key,
            )
        )

    if detail_frames:
        detail_checkpoint = pd.concat(
            detail_frames,
            axis=0,
            sort=False,
        )
    else:
        detail_checkpoint = pd.DataFrame()

    timing_checkpoint = pd.DataFrame(
        timing_records
    )

    atomic_save_csv(
        summary_checkpoint,
        checkpoint_summary_csv,
        index=True,
    )

    atomic_save_csv(
        detail_checkpoint,
        checkpoint_detail_csv,
        index=True,
    )

    atomic_save_csv(
        timing_checkpoint,
        checkpoint_timing_csv,
        index=False,
    )

    save_warning_log(
        run_warnings,
        checkpoint_warnings_txt,
    )

    checkpoint_seconds = (
        time.perf_counter()
        - checkpoint_start
    )

    return checkpoint_seconds


# =============================================================================
# Start total timer
# =============================================================================

script_start_time = time.perf_counter()
script_start_timestamp = current_timestamp()

print("=" * 80)
print("BOLTZMANN-AVERAGED DESCRIPTOR EXTRACTION")
print("=" * 80)
print(f"Started: {script_start_timestamp}")
print(f"Input: {conformers_dir}")
print(f"Temperature: {temperature} K")
print(f"Energy cutoff: {energy_cutoff_kcal}")
print(f"Minimum weight: {min_weight}")
print("=" * 80)


# =============================================================================
# Create output folders
# =============================================================================

os.makedirs(
    output_dir,
    exist_ok=True,
)

os.makedirs(
    checkpoint_dir,
    exist_ok=True,
)


# =============================================================================
# Load CREST ligand folders
# =============================================================================

dataset_load_start = time.perf_counter()

dataset = CrestLigandSet(
    conformers_dir,
    temperature=temperature,
)

dataset_load_seconds = (
    time.perf_counter()
    - dataset_load_start
)

print(
    f"\nLoaded {len(dataset.ensembles)} ligand folders "
    f"in {format_duration(dataset_load_seconds)}"
)

print(
    "Skipped during dataset loading:",
    list(dataset.failed_ligands.keys()),
)


# =============================================================================
# Identify ordinary and merged ligand folders
# =============================================================================

merged_sources = {
    source_name
    for source_names in MERGE_GROUPS.values()
    for source_name in source_names
}

plain_ligands = {
    ligand_name: ensemble
    for ligand_name, ensemble
    in dataset.ensembles.items()
    if ligand_name not in merged_sources
}


# =============================================================================
# Pool rerun folders safely
# =============================================================================

final_ensembles = dict(plain_ligands)

merge_warnings = []

print("\nPooling rerun folders:")

for merged_name, source_names in MERGE_GROUPS.items():

    merge_start = time.perf_counter()

    available_directories = [
        os.path.join(
            conformers_dir,
            source_name,
        )
        for source_name in source_names
        if source_name in dataset.ensembles
    ]

    missing_sources = [
        source_name
        for source_name in source_names
        if source_name not in dataset.ensembles
    ]

    if missing_sources:
        warning = (
            f"[{merged_name}] Missing source folders: "
            f"{missing_sources}"
        )
        merge_warnings.append(warning)
        print(f"  WARNING: {warning}")

    if not available_directories:
        warning = (
            f"[{merged_name}] No available source folders."
        )
        merge_warnings.append(warning)
        print(f"  SKIPPED: {warning}")
        continue

    try:
        final_ensembles[merged_name] = (
            ConformerEnsemble.from_multiple(
                available_directories,
                molecule_name=merged_name,
                temperature=temperature,
            )
        )

        merge_seconds = (
            time.perf_counter()
            - merge_start
        )

        print(
            f"  {merged_name}: pooled "
            f"{len(available_directories)} folders in "
            f"{format_duration(merge_seconds)}"
        )

    except Exception as error:
        merge_seconds = (
            time.perf_counter()
            - merge_start
        )

        warning = (
            f"[{merged_name}] Merge failed after "
            f"{format_duration(merge_seconds)}: "
            f"{type(error).__name__}: {error}"
        )

        merge_warnings.append(warning)
        print(f"  FAILED: {warning}")


print(
    f"\nFinal ligand count after merging: "
    f"{len(final_ensembles)}"
)


# =============================================================================
# Initialize result containers
# =============================================================================

summary_rows = {}
detail_frames = []
timing_records = []

run_warnings = [
    f"[{name}] Dataset loading failed: {error}"
    for name, error in dataset.failed_ligands.items()
]

run_warnings.extend(merge_warnings)

completed_ligands = set()


# =============================================================================
# Resume from an existing checkpoint
# =============================================================================

if (
    resume_from_checkpoint
    and os.path.exists(checkpoint_summary_csv)
):
    try:
        previous_summary = pd.read_csv(
            checkpoint_summary_csv,
            index_col=0,
        )

        previous_summary.index = (
            previous_summary.index.astype(str)
        )

        summary_rows = previous_summary.to_dict(
            orient="index",
        )

        # Only rows explicitly marked completed are skipped
        if "run_status" in previous_summary.columns:
            completed_mask = (
                previous_summary["run_status"]
                .astype(str)
                .str.startswith("completed")
            )

            completed_ligands = set(
                previous_summary.index[
                    completed_mask
                ]
            )
        else:
            completed_ligands = set(
                previous_summary.index
            )

        print(
            f"\nResume checkpoint found: "
            f"{len(completed_ligands)} completed ligands"
        )

    except Exception as error:
        print(
            "\nCould not load summary checkpoint:"
        )
        print(
            f"{type(error).__name__}: {error}"
        )

        summary_rows = {}
        completed_ligands = set()


if (
    resume_from_checkpoint
    and os.path.exists(checkpoint_detail_csv)
):
    try:
        previous_detail = pd.read_csv(
            checkpoint_detail_csv,
            index_col=0,
        )

        if not previous_detail.empty:
            detail_frames.append(
                previous_detail
            )

        print(
            f"Loaded {len(previous_detail)} "
            f"per-conformer checkpoint rows"
        )

    except Exception as error:
        print(
            "Could not load the per-conformer checkpoint:"
        )
        print(
            f"{type(error).__name__}: {error}"
        )


if (
    resume_from_checkpoint
    and os.path.exists(checkpoint_timing_csv)
):
    try:
        previous_timing = pd.read_csv(
            checkpoint_timing_csv,
        )

        if not previous_timing.empty:
            timing_records.extend(
                previous_timing.to_dict(
                    orient="records",
                )
            )

        print(
            f"Loaded {len(previous_timing)} "
            f"timing records"
        )

    except Exception as error:
        print(
            "Could not load timing checkpoint:"
        )
        print(
            f"{type(error).__name__}: {error}"
        )


if (
    resume_from_checkpoint
    and os.path.exists(checkpoint_warnings_txt)
):
    try:
        with open(
            checkpoint_warnings_txt,
            "r",
            encoding="utf-8",
        ) as warnings_file:
            previous_warnings = (
                warnings_file.read().strip()
            )

        if (
            previous_warnings
            and previous_warnings != "No warnings."
        ):
            run_warnings.insert(
                0,
                previous_warnings,
            )

    except Exception as error:
        print(
            "Could not load warning checkpoint:"
        )
        print(
            f"{type(error).__name__}: {error}"
        )


# =============================================================================
# Prepare ligand processing order
# =============================================================================

sorted_ligand_names = sorted(
    final_ensembles,
    key=natural_key,
)

ligands_to_process = [
    ligand_name
    for ligand_name in sorted_ligand_names
    if ligand_name not in completed_ligands
]

print("\n" + "=" * 80)
print(f"Total ligands: {len(sorted_ligand_names)}")
print(f"Already completed: {len(completed_ligands)}")
print(f"Remaining: {len(ligands_to_process)}")
print("=" * 80)


# Times from only the current execution
current_run_ligand_times = []
current_run_completed = 0


# =============================================================================
# Main descriptor extraction loop
# =============================================================================

try:

    for overall_position, ligand_name in enumerate(
        sorted_ligand_names,
        start=1,
    ):

        if ligand_name in completed_ligands:
            print(
                f"\n[{overall_position}/"
                f"{len(sorted_ligand_names)}] "
                f"{ligand_name}: already completed"
            )
            continue

        ensemble = final_ensembles[ligand_name]

        ligand_start = time.perf_counter()
        ligand_started_at = current_timestamp()

        print("\n" + "-" * 80)
        print(
            f"[{overall_position}/"
            f"{len(sorted_ligand_names)}] "
            f"Processing {ligand_name}"
        )
        print(f"Started: {ligand_started_at}")
        print("-" * 80)

        row = {
            "run_status": "incomplete",
            "started_at": ligand_started_at,
        }

        ligand_detail_parts = []
        ligand_failure_count = 0
        ligand_success_count = 0
        unexpected_ligand_failure = False

        try:

            # =================================================================
            # Sterimol descriptors
            # =================================================================

            for origin, attached in sterimol_pairs:

                pair = [origin, attached]
                pair_label = f"{origin}-{attached}"
                descriptor_name = (
                    f"sterimol_{pair_label}"
                )

                success, per_conf, boltz = (
                    safe_descriptor_call(
                        ligand_name=ligand_name,
                        descriptor_type="sterimol",
                        descriptor_name=descriptor_name,
                        calculation=(
                            lambda pair=pair.copy():
                            ensemble.get_sterimol(
                                base_atoms=pair,
                                radii=radii,
                                energy_cutoff_kcal=(
                                    energy_cutoff_kcal
                                ),
                                min_weight=min_weight,
                            )
                        ),
                        run_warnings=run_warnings,
                        timing_records=timing_records,
                    )
                )

                sterimol_features = [
                    "B1",
                    "B5",
                    "L",
                    "loc_B5",
                    "B1_B5_angle",
                ]

                if not success:
                    ligand_failure_count += 1

                    for feature in sterimol_features:
                        row[
                            f"sterimol_"
                            f"{pair_label}_{feature}"
                        ] = np.nan

                    row[
                        f"n_conformers_used_"
                        f"{pair_label}"
                    ] = 0

                    continue

                ligand_success_count += 1

                extract_metadata(
                    boltzmann_result=boltz,
                    row=row,
                    pair_label=pair_label,
                )

                for feature in sterimol_features:

                    output_column = (
                        f"sterimol_"
                        f"{pair_label}_{feature}"
                    )

                    try:
                        row[output_column] = (
                            scalar_from_result(
                                boltz,
                                preferred_columns=[
                                    feature
                                ],
                            )
                        )

                    except Exception as error:
                        row[output_column] = np.nan
                        ligand_failure_count += 1

                        run_warnings.append(
                            format_exception(
                                ligand_name,
                                output_column,
                                error,
                            )
                        )

                if (
                    per_conf is not None
                    and isinstance(
                        per_conf,
                        pd.DataFrame,
                    )
                ):
                    ligand_detail_parts.append(
                        per_conf.add_prefix(
                            f"sterimol_"
                            f"{pair_label}_"
                        )
                    )

            # =================================================================
            # Angles
            # =================================================================

            for group in angle_groups:

                group_copy = group.copy()
                label = "-".join(
                    map(str, group)
                )
                column_name = f"angle_{label}"

                success, per_conf, boltz = (
                    safe_descriptor_call(
                        ligand_name=ligand_name,
                        descriptor_type="angle",
                        descriptor_name=column_name,
                        calculation=(
                            lambda atoms=group_copy:
                            ensemble.get_angle(
                                atoms=atoms,
                                energy_cutoff_kcal=(
                                    energy_cutoff_kcal
                                ),
                                min_weight=min_weight,
                            )
                        ),
                        run_warnings=run_warnings,
                        timing_records=timing_records,
                    )
                )

                if not success:
                    row[column_name] = np.nan
                    ligand_failure_count += 1
                    continue

                ligand_success_count += 1

                try:
                    row[column_name] = (
                        scalar_from_result(
                            boltz,
                            preferred_columns=[
                                "angle",
                                column_name,
                                "value",
                            ],
                        )
                    )

                except Exception as error:
                    row[column_name] = np.nan
                    ligand_failure_count += 1

                    run_warnings.append(
                        format_exception(
                            ligand_name,
                            column_name,
                            error,
                        )
                    )

                if (
                    per_conf is not None
                    and isinstance(
                        per_conf,
                        pd.DataFrame,
                    )
                ):
                    ligand_detail_parts.append(
                        per_conf.add_prefix(
                            f"{column_name}_"
                        )
                    )

            # =================================================================
            # Dihedral angles
            # =================================================================

            for group in dihedral_groups:

                group_copy = group.copy()
                label = "-".join(
                    map(str, group)
                )
                column_name = (
                    f"dihedral_{label}"
                )

                success, per_conf, boltz = (
                    safe_descriptor_call(
                        ligand_name=ligand_name,
                        descriptor_type="dihedral",
                        descriptor_name=column_name,
                        calculation=(
                            lambda atoms=group_copy:
                            ensemble.get_dihedral(
                                atoms=atoms,
                                energy_cutoff_kcal=(
                                    energy_cutoff_kcal
                                ),
                                min_weight=min_weight,
                            )
                        ),
                        run_warnings=run_warnings,
                        timing_records=timing_records,
                    )
                )

                if not success:
                    row[column_name] = np.nan
                    ligand_failure_count += 1
                    continue

                ligand_success_count += 1

                try:
                    row[column_name] = (
                        scalar_from_result(
                            boltz,
                            preferred_columns=[
                                "dihedral",
                                column_name,
                                "value",
                            ],
                        )
                    )

                except Exception as error:
                    row[column_name] = np.nan
                    ligand_failure_count += 1

                    run_warnings.append(
                        format_exception(
                            ligand_name,
                            column_name,
                            error,
                        )
                    )

                if (
                    per_conf is not None
                    and isinstance(
                        per_conf,
                        pd.DataFrame,
                    )
                ):
                    ligand_detail_parts.append(
                        per_conf.add_prefix(
                            f"{column_name}_"
                        )
                    )

            # =================================================================
            # Bond lengths
            # =================================================================

            for atom1, atom2 in bond_length_pairs:

                pair = [atom1, atom2]
                label = f"{atom1}-{atom2}"

                column_name = (
                    f"bond_length_{label}"
                )

                success, per_conf, boltz = (
                    safe_descriptor_call(
                        ligand_name=ligand_name,
                        descriptor_type="bond_length",
                        descriptor_name=column_name,
                        calculation=(
                            lambda atoms=pair.copy():
                            ensemble.get_bond_length(
                                atoms=atoms,
                                energy_cutoff_kcal=(
                                    energy_cutoff_kcal
                                ),
                                min_weight=min_weight,
                            )
                        ),
                        run_warnings=run_warnings,
                        timing_records=timing_records,
                    )
                )

                if not success:
                    row[column_name] = np.nan
                    ligand_failure_count += 1
                    continue

                ligand_success_count += 1

                try:
                    row[column_name] = (
                        scalar_from_result(
                            boltz,
                            preferred_columns=[
                                "bond_length",
                                "distance",
                                column_name,
                                "value",
                            ],
                        )
                    )

                except Exception as error:
                    row[column_name] = np.nan
                    ligand_failure_count += 1

                    run_warnings.append(
                        format_exception(
                            ligand_name,
                            column_name,
                            error,
                        )
                    )

                if (
                    per_conf is not None
                    and isinstance(
                        per_conf,
                        pd.DataFrame,
                    )
                ):
                    ligand_detail_parts.append(
                        per_conf.add_prefix(
                            f"{column_name}_"
                        )
                    )

            # =================================================================
            # Buried volume
            # =================================================================

            buried_volume_column = (
                f"buried_volume_"
                f"{metal_index}_"
                f"{buried_volume_radius:g}A"
            )

            success, per_conf, boltz = (
                safe_descriptor_call(
                    ligand_name=ligand_name,
                    descriptor_type="buried_volume",
                    descriptor_name=(
                        buried_volume_column
                    ),
                    calculation=(
                        lambda:
                        ensemble.get_buried_volume(
                            metal_index=metal_index,
                            radius=(
                                buried_volume_radius
                            ),
                            energy_cutoff_kcal=(
                                energy_cutoff_kcal
                            ),
                            min_weight=min_weight,
                        )
                    ),
                    run_warnings=run_warnings,
                    timing_records=timing_records,
                )
            )

            if not success:
                row[
                    buried_volume_column
                ] = np.nan
                ligand_failure_count += 1

            else:
                ligand_success_count += 1

                try:
                    row[
                        buried_volume_column
                    ] = scalar_from_result(
                        boltz,
                        preferred_columns=[
                            "buried_volume",
                            "percent_buried_volume",
                            "fraction_buried_volume",
                            "V_bur",
                            "%V_bur",
                            "value",
                        ],
                    )

                except Exception as error:
                    row[
                        buried_volume_column
                    ] = np.nan
                    ligand_failure_count += 1

                    run_warnings.append(
                        format_exception(
                            ligand_name,
                            buried_volume_column,
                            error,
                        )
                    )

                if (
                    per_conf is not None
                    and isinstance(
                        per_conf,
                        pd.DataFrame,
                    )
                ):
                    ligand_detail_parts.append(
                        per_conf.add_prefix(
                            f"{buried_volume_column}_"
                        )
                    )

            # =================================================================
            # Combine per-conformer descriptor tables
            # =================================================================

            if ligand_detail_parts:

                detail_combine_start = (
                    time.perf_counter()
                )

                try:
                    ligand_detail = pd.concat(
                        ligand_detail_parts,
                        axis=1,
                    )

                    ligand_detail.insert(
                        0,
                        "ligand",
                        ligand_name,
                    )

                    detail_frames.append(
                        ligand_detail
                    )

                    detail_combine_seconds = (
                        time.perf_counter()
                        - detail_combine_start
                    )

                    add_timing_record(
                        timing_records=timing_records,
                        ligand_name=ligand_name,
                        descriptor_type="internal",
                        descriptor_name=(
                            "combine_per_conformer_tables"
                        ),
                        status="success",
                        seconds=(
                            detail_combine_seconds
                        ),
                    )

                except Exception as error:
                    detail_combine_seconds = (
                        time.perf_counter()
                        - detail_combine_start
                    )

                    ligand_failure_count += 1

                    run_warnings.append(
                        format_exception(
                            ligand_name,
                            (
                                "combine per-conformer "
                                "tables"
                            ),
                            error,
                        )
                    )

                    add_timing_record(
                        timing_records=timing_records,
                        ligand_name=ligand_name,
                        descriptor_type="internal",
                        descriptor_name=(
                            "combine_per_conformer_tables"
                        ),
                        status="failed",
                        seconds=(
                            detail_combine_seconds
                        ),
                        error=error,
                    )

        except Exception as error:

            unexpected_ligand_failure = True
            ligand_failure_count += 1

            run_warnings.append(
                format_exception(
                    ligand_name,
                    "unexpected ligand-level error",
                    error,
                )
            )

            print(
                "  Unexpected ligand-level failure: "
                f"{type(error).__name__}: {error}"
            )

        # =====================================================================
        # Finish ligand timing and status
        # =====================================================================

        ligand_seconds = (
            time.perf_counter()
            - ligand_start
        )

        current_run_ligand_times.append(
            ligand_seconds
        )

        current_run_completed += 1

        row["finished_at"] = current_timestamp()
        row["ligand_seconds"] = ligand_seconds
        row["ligand_time"] = format_duration(
            ligand_seconds
        )
        row["successful_calculations"] = (
            ligand_success_count
        )
        row["failed_calculations"] = (
            ligand_failure_count
        )

        if unexpected_ligand_failure:
            row["run_status"] = "incomplete"

        elif ligand_failure_count > 0:
            row["run_status"] = (
                "completed_with_errors"
            )

        else:
            row["run_status"] = "completed"

        summary_rows[ligand_name] = row

        try:
            ensemble_warnings = getattr(
                ensemble,
                "warnings",
                [],
            )

            if ensemble_warnings:
                run_warnings.extend(
                    [
                        f"[{ligand_name}] {warning}"
                        for warning
                        in ensemble_warnings
                    ]
                )

        except Exception:
            pass

        # Ligand timing record
        add_timing_record(
            timing_records=timing_records,
            ligand_name=ligand_name,
            descriptor_type="ligand_total",
            descriptor_name="all_descriptors",
            status=row["run_status"],
            seconds=ligand_seconds,
        )

        # =====================================================================
        # Save checkpoint
        # =====================================================================

        checkpoint_seconds = 0.0

        if save_checkpoint_after_each_ligand:
            try:
                checkpoint_seconds = (
                    save_checkpoint(
                        summary_rows=summary_rows,
                        detail_frames=detail_frames,
                        timing_records=(
                            timing_records
                        ),
                        run_warnings=run_warnings,
                    )
                )

            except Exception as error:
                run_warnings.append(
                    format_exception(
                        ligand_name,
                        "save checkpoint",
                        error,
                    )
                )

                print(
                    "  WARNING: checkpoint save "
                    f"failed: {error}"
                )

        # =====================================================================
        # Calculate ETA
        # =====================================================================

        mean_ligand_seconds = float(
            np.mean(current_run_ligand_times)
        )

        remaining_ligands = (
            len(ligands_to_process)
            - current_run_completed
        )

        estimated_remaining_seconds = (
            mean_ligand_seconds
            * max(remaining_ligands, 0)
        )

        total_elapsed_seconds = (
            time.perf_counter()
            - script_start_time
        )

        missing_values = sum(
            pd.isna(value)
            for key, value in row.items()
            if key not in {
                "run_status",
                "started_at",
                "finished_at",
                "ligand_time",
            }
        )

        print(
            f"\n  {ligand_name} finished in "
            f"{format_duration(ligand_seconds)}"
        )

        print(
            f"  Status: {row['run_status']} | "
            f"successful={ligand_success_count} | "
            f"failed={ligand_failure_count} | "
            f"missing values={missing_values}"
        )

        if save_checkpoint_after_each_ligand:
            print(
                "  Checkpoint saved in "
                f"{format_duration(checkpoint_seconds)}"
            )

        print(
            "  Total elapsed: "
            f"{format_duration(total_elapsed_seconds)}"
        )

        print(
            "  Mean time per ligand: "
            f"{format_duration(mean_ligand_seconds)}"
        )

        print(
            f"  Remaining ligands: "
            f"{remaining_ligands}"
        )

        print(
            "  Estimated time remaining: "
            f"{format_duration(estimated_remaining_seconds)}"
        )


# =============================================================================
# Save safely when interrupted manually
# =============================================================================

except KeyboardInterrupt:

    print("\n" + "=" * 80)
    print("RUN INTERRUPTED BY USER")
    print("=" * 80)
    print("Saving current checkpoint...")

    try:
        checkpoint_seconds = save_checkpoint(
            summary_rows=summary_rows,
            detail_frames=detail_frames,
            timing_records=timing_records,
            run_warnings=run_warnings,
        )

        print(
            "Checkpoint saved in "
            f"{format_duration(checkpoint_seconds)}"
        )

        print(
            "Restart this cell to continue from "
            "the saved checkpoint."
        )

    except Exception as checkpoint_error:
        print(
            "Checkpoint saving also failed:"
        )
        print(
            f"{type(checkpoint_error).__name__}: "
            f"{checkpoint_error}"
        )

    raise


# =============================================================================
# Save safely after an unexpected global failure
# =============================================================================

except Exception as fatal_error:

    run_warnings.append(
        format_exception(
            "GLOBAL",
            "unexpected fatal error",
            fatal_error,
        )
    )

    print("\n" + "=" * 80)
    print("UNEXPECTED GLOBAL ERROR")
    print("=" * 80)
    print(
        f"{type(fatal_error).__name__}: "
        f"{fatal_error}"
    )
    print("Saving current checkpoint...")

    try:
        checkpoint_seconds = save_checkpoint(
            summary_rows=summary_rows,
            detail_frames=detail_frames,
            timing_records=timing_records,
            run_warnings=run_warnings,
        )

        print(
            "Checkpoint saved in "
            f"{format_duration(checkpoint_seconds)}"
        )

    except Exception as checkpoint_error:
        print(
            "Checkpoint saving also failed:"
        )
        print(
            f"{type(checkpoint_error).__name__}: "
            f"{checkpoint_error}"
        )

    raise


# =============================================================================
# Assemble final output tables
# =============================================================================

summary_df = pd.DataFrame.from_dict(
    summary_rows,
    orient="index",
)

summary_df.index.name = "ligand"

if not summary_df.empty:
    summary_df = summary_df.reindex(
        sorted(
            summary_df.index,
            key=natural_key,
        )
    )


if detail_frames:
    detail_df = pd.concat(
        detail_frames,
        axis=0,
        sort=False,
    )
else:
    detail_df = pd.DataFrame()


timing_df = pd.DataFrame(
    timing_records
)


# =============================================================================
# Print final run summary
# =============================================================================

total_script_seconds = (
    time.perf_counter()
    - script_start_time
)

if "run_status" in summary_df.columns:

    completed_count = int(
        summary_df["run_status"]
        .astype(str)
        .str.startswith("completed")
        .sum()
    )

    error_count = int(
        (
            summary_df["run_status"]
            == "completed_with_errors"
        ).sum()
    )

    incomplete_count = int(
        (
            summary_df["run_status"]
            == "incomplete"
        ).sum()
    )

else:
    completed_count = len(summary_df)
    error_count = 0
    incomplete_count = 0


print("\n" + "=" * 80)
print("EXTRACTION FINISHED")
print("=" * 80)
print(f"Started: {script_start_timestamp}")
print(f"Finished: {current_timestamp()}")
print(
    "Total execution time: "
    f"{format_duration(total_script_seconds)}"
)
print(f"Expected ligands: {len(final_ensembles)}")
print(f"Rows produced: {len(summary_df)}")
print(f"Completed: {completed_count}")
print(f"Completed with errors: {error_count}")
print(f"Incomplete: {incomplete_count}")
print(f"Summary shape: {summary_df.shape}")
print(f"Per-conformer shape: {detail_df.shape}")
print(f"Timing records: {len(timing_df)}")


# =============================================================================
# Save final output files
# =============================================================================

if save_csv:

    summary_path = os.path.join(
        output_dir,
        output_csv,
    )

    detail_path = os.path.join(
        output_dir,
        detail_csv,
    )

    warnings_path = os.path.join(
        output_dir,
        warnings_txt,
    )

    timing_path = os.path.join(
        output_dir,
        timing_csv,
    )

    final_save_start = time.perf_counter()

    atomic_save_csv(
        summary_df,
        summary_path,
        index=True,
    )

    atomic_save_csv(
        detail_df,
        detail_path,
        index=True,
    )

    atomic_save_csv(
        timing_df,
        timing_path,
        index=False,
    )

    save_warning_log(
        run_warnings,
        warnings_path,
    )

    final_save_seconds = (
        time.perf_counter()
        - final_save_start
    )

    print(
        "\nFinal files saved in "
        f"{format_duration(final_save_seconds)}:"
    )
    print(f"  Summary: {summary_path}")
    print(f"  Per-conformer: {detail_path}")
    print(f"  Timing: {timing_path}")
    print(f"  Warnings: {warnings_path}")


# =============================================================================
# Display slowest calculations
# =============================================================================

if (
    not timing_df.empty
    and "seconds" in timing_df.columns
):

    timing_df["seconds"] = pd.to_numeric(
        timing_df["seconds"],
        errors="coerce",
    )

    descriptor_timings = timing_df[
        timing_df["descriptor_type"]
        != "ligand_total"
    ].copy()

    slowest = descriptor_timings.nlargest(
        10,
        "seconds",
    )[
        [
            "ligand",
            "descriptor_type",
            "descriptor",
            "status",
            "seconds",
            "formatted_time",
        ]
    ]

    print("\nTen slowest descriptor calculations:")
    print(
        slowest.to_string(
            index=False,
        )
    )