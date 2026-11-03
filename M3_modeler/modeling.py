try:
    from utils.console import make_console_safe as _mcs; _mcs()
except ImportError:
    pass
# -*- coding: latin-1 -*-
import cProfile
import copy
import os
if not os.environ.get("MPLBACKEND"):
    os.environ["MPLBACKEND"] = "Agg"   # headless, file-only
import pstats
import random
import sqlite3
import sys
import time
import multiprocessing
from itertools import combinations
from tkinter import filedialog, messagebox
from typing import Iterable, List, Sequence, Set, Tuple, Union, Optional

# --- Third-party ---
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import PIL
import scipy.stats as stats
from tqdm import tqdm
from joblib import Parallel, delayed
from scipy.stats import t

# scikit-learn
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import (
    BayesianRidge,
    Lasso,
    LinearRegression,
    SGDRegressor,
)
from sklearn.utils.validation import check_is_fitted
from sklearn.exceptions import NotFittedError 
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    make_scorer,
    mean_absolute_error,
    mean_squared_error,
    precision_score,
    r2_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import (
    KFold,
    LeaveOneOut,
    RepeatedKFold,
    cross_val_predict,
    cross_validate,
    cross_val_score,
    train_test_split,
)
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.outliers_influence import variance_inflation_factor



_this_dir   = os.path.dirname(os.path.abspath(__file__))
_parent_dir = os.path.dirname(_this_dir)
for _p in (_this_dir, _parent_dir):
    if _p not in sys.path:
        sys.path.insert(0, _p)




_NAME_COLUMNS = ("name", "names", "molecule", "molecule_name", "mol", "compound", "sample", "id")


def merge_features_and_target(features_csv, target_csv, y_value: str) -> pd.DataFrame:
    """One table from a features CSV and a separate target CSV, joined on the molecule name.

    Each file's name column is the first of ``name``, ``molecule``, ``id``, ... it has, or else
    its first column (the extractor writes names as an unnamed first column). Rows present in
    only one file are dropped, and the drop is reported.
    """
    feats, target = pd.read_csv(features_csv), pd.read_csv(target_csv)

    def name_col(df):
        return next((c for c in df.columns if str(c).strip().lower() in _NAME_COLUMNS), df.columns[0])

    tcol = next((c for c in target.columns if str(c).strip().lower() == str(y_value).strip().lower()), None)
    if tcol is None:
        raise KeyError(f"Target column {y_value!r} not found in {target_csv}. Available: {list(target.columns)}")
    fn, tn = name_col(feats), name_col(target)
    feats = feats.rename(columns={fn: "name"})
    feats = feats.drop(columns=[c for c in feats.columns if str(c).strip().lower() == str(y_value).strip().lower()])
    t = target[[tn, tcol]].rename(columns={tn: "name", tcol: y_value})
    feats["name"], t["name"] = feats["name"].astype(str), t["name"].astype(str)
    merged = feats.merge(t, on="name", how="inner")
    only_f = sorted(set(feats["name"]) - set(merged["name"]))
    only_t = sorted(set(t["name"]) - set(merged["name"]))
    if only_f or only_t:
        print(f"merge_features_and_target: {len(merged)} molecules in both files; "
              f"dropped {len(only_f)} without a target {only_f[:5]} and {len(only_t)} without features {only_t[:5]}")
    if merged.empty:
        raise ValueError("No molecule name is shared by the features CSV and the target CSV.")
    return merged

def _combo_key(combo: Sequence[str]) -> str:
    """Canonical string key for a combo (order-insensitive)."""
    return ",".join(sorted(map(str, combo)))

def _to_df(maybe_df) -> pd.DataFrame:
    """Coerce list/None/DataFrame to a DataFrame."""
    
    if isinstance(maybe_df, pd.DataFrame):
        return maybe_df
    if maybe_df is None:
        return pd.DataFrame()
    if isinstance(maybe_df, list):
        return pd.DataFrame(maybe_df)
    return pd.DataFrame(maybe_df)

def _ensure_numeric(s: pd.Series) -> pd.Series:
    """Coerce to numeric; invalid -> NaN."""
    return pd.to_numeric(s, errors="coerce")

def _extract_series(df: pd.DataFrame, names: Tuple[str, ...]) -> Optional[pd.Series]:
    """Return the first matching column (case-insensitive), else None."""
    if df.empty:
        return None
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in lower:
            return df[lower[n.lower()]]
    return None

def _extract_q2(df: pd.DataFrame) -> Optional[pd.Series]:
    """Get Q2 (or q2) as numeric; if only in 'scores' dicts, pull from there."""
    if df.empty:
        return None
    s = _extract_series(df, ("q2", "Q2"))
    if s is not None:
        return _ensure_numeric(s)
    # scores column path
    scores = _extract_series(df, ("scores",))
    if scores is not None:
        s = scores.apply(lambda d: (d or {}).get("q2", (d or {}).get("Q2", np.nan)))
        return _ensure_numeric(s)
    return None

def _adjusted_r2(r2: float, n: int, p: int) -> float:
    """Adjusted R²: penalises extra predictors.  Returns NaN when ill-defined."""
    denom = n - p - 1
    if denom <= 0:
        return float("nan")
    return float(1.0 - (1.0 - r2) * (n - 1) / denom)


def _extract_r2(df: pd.DataFrame) -> Optional[pd.Series]:
    # prefer adj_r2 when available, fall back to plain r2
    s = _extract_series(df, ("adj_r2", "adj_R2"))
    if s is None:
        s = _extract_series(df, ("r2", "R2"))
    return _ensure_numeric(s) if s is not None else None


def _resolve_effective_jobs(n_jobs: int, bool_parallel: bool, cpu_count: int) -> int:
    """
    Resolve the worker count for model evaluation.

    `n_jobs` is treated as the primary knob:
    - `-1` means all detected cores.
    - values > 1 enable parallel execution even when `bool_parallel` is False.
    - `bool_parallel=True` keeps backward compatibility for callers that rely on
      the flag and omit `n_jobs`.
    """
    if cpu_count <= 1:
        return 1

    requested_jobs = cpu_count if n_jobs == -1 else max(1, int(n_jobs))
    requested_jobs = min(requested_jobs, cpu_count)

    if bool_parallel or requested_jobs > 1:
        return requested_jobs
    return 1

def _is_all_q2_neg_inf(df: pd.DataFrame) -> bool:
    """True iff Q2 exists and every non-NaN value is exactly -inf (NaN => not -inf)."""
    q2 = _extract_q2(df)
    if q2 is None or q2.dropna().empty:
        return False
    return bool((q2.dropna() == -np.inf).all())

def _best_r2(df: pd.DataFrame) -> Optional[float]:
    r2 = _extract_r2(df)
    if r2 is None or r2.dropna().empty:
        return None
    return float(r2.max())

def _sort_results(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sort model results by Q² (descending) if present, then by adj_r2 / R².
    Falls back to R²-only sorting if Q² not found.

    Parameters
    ----------
    df : pd.DataFrame
        Results dataframe expected to contain columns like:
        ['id', 'threshold', 'model', 'predictions', 'q2', 'adj_r2', 'r2', ...]

    Returns
    -------
    pd.DataFrame
        Sorted dataframe (copy), cleaned of non-metric columns.
    """
    if df is None or df.empty:
        print("WARNING: Empty DataFrame received in _sort_results.")
        return pd.DataFrame()

    # drop non-metric columns if they exist
    drop_cols = [c for c in ['id', 'threshold', 'model', 'predictions'] if c in df.columns]
    df = df.drop(columns=drop_cols, errors='ignore').copy()

    # filter out any rows with inf or NaN values in numeric columns
    numeric_cols = df.select_dtypes(include=[float, int]).columns
    if len(numeric_cols) > 0:
        df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=numeric_cols, how="any")

    if df.empty:
        print("WARNING: No valid rows after cleaning in _sort_results.")
        return df

    # --- your existing logic follows ---
    q2 = _extract_q2(df)
    r2 = _extract_r2(df)

    temp = df.copy()
    if q2 is not None:
        temp["_q2__"] = q2
    if r2 is not None:
        temp["_r2__"] = r2

    by = []
    if "_q2__" in temp.columns:
        by.append("_q2__")
    if "_r2__" in temp.columns:
        by.append("_r2__")

    if not by:
        return df

    temp = temp.sort_values(by=by, ascending=[False] * len(by))
    tmp_re = temp.drop(columns=[c for c in ("_q2__", "_r2__") if c in temp.columns])
    
    return tmp_re 



        # Print the count of combinations for each number of features
def _linear_matrix_fit_predict(
    X_train,
    y_train,
    X_test,
    *,
    alpha: float = 1e-5,
    scale: bool = True,
):
    """
    Fit linear regression with matrix operations and predict on X_test.

    Scaling is fit on X_train only, then applied to X_test, matching the
    no-leakage behavior expected during cross-validation.
    """
    X_train = np.asarray(X_train, dtype=float)
    X_test = np.asarray(X_test, dtype=float)
    y_train = np.asarray(y_train, dtype=float).ravel()

    if scale:
        mean = X_train.mean(axis=0)
        std = X_train.std(axis=0, ddof=0)
        std = np.where(std == 0, 1.0, std)
        X_train = (X_train - mean) / std
        X_test = (X_test - mean) / std

    Xb_train = np.c_[np.ones((X_train.shape[0], 1)), X_train]
    Xb_test = np.c_[np.ones((X_test.shape[0], 1)), X_test]

    gram = Xb_train.T @ Xb_train
    penalty = np.eye(gram.shape[0])
    penalty[0, 0] = 0.0
    rhs = Xb_train.T @ y_train

    try:
        theta = np.linalg.solve(gram + alpha * penalty, rhs)
    except np.linalg.LinAlgError:
        theta = np.linalg.pinv(gram + alpha * penalty) @ rhs

    return Xb_test @ theta


def _persist_regression_result(result_dict, r2_threshold):
    """Persist one regression result to the configured DB/CSV outputs."""
    if 'scores' not in result_dict:
        return result_dict

    scores = result_dict['scores']
    model = result_dict['model']
    combination = result_dict['combination']
    y_pred = result_dict['predictions']

    try:
        r2 = scores.get('r2', float('-inf'))
        adj_r2 = scores.get('adj_r2', float('-inf'))
        q2 = scores.get('Q2', float('-inf'))
        mae = scores.get('MAE', float('inf'))
        rmsd = scores.get('RMSD', float('inf'))
        csv_path = model.db_path.replace('.db', '.csv')

        insert_result_into_db_regression(
            db_path=model.db_path,
            combination=combination,
            r2=r2,
            adj_r2=adj_r2,
            q2=q2,
            mae=mae,
            rmsd=rmsd,
            threshold=r2_threshold,
            csv_path=csv_path,
            model=model,
            predictions=y_pred,
        )
    except Exception as e:
        print(f'failed to insert combination : {combination}')
        print(e)

    return {
        'combination': str(combination),
        'r2': r2,
        'adj_r2': adj_r2,
        'q2': q2,
        'mae': mae,
        'rmsd': rmsd,
        'threshold': r2_threshold,
        'model': model,
        'predictions': y_pred,
    }


def _analytic_loo_linear(X: np.ndarray, y: np.ndarray, alpha: float = 1e-5):
    """
    Exact LOO Q²/MAE/RMSD for ridge-regularised linear regression in O(n·p²).

    Uses the hat-matrix identity  e_loo[i] = e[i] / (1 - h[i])  so no
    re-fitting loop is needed — equivalent to calling calculate_q2_and_mae
    with n_splits=1 but orders-of-magnitude faster for large n.
    """
    from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).ravel()
    n = X.shape[0]

    mean = X.mean(axis=0)
    std  = X.std(axis=0, ddof=0)
    std  = np.where(std == 0, 1.0, std)
    Xs   = (X - mean) / std

    Xb = np.c_[np.ones(n), Xs]
    gram = Xb.T @ Xb
    pen  = np.eye(gram.shape[0])
    pen[0, 0] = 0.0                        # no penalty on intercept
    A_inv = np.linalg.pinv(gram + alpha * pen)

    theta  = A_inv @ (Xb.T @ y)
    y_hat  = Xb @ theta
    e      = y - y_hat

    # Hat-matrix diagonal: h_ii = x_i^T A_inv x_i  (avoid building full n×n H)
    h = np.einsum("ij,jk,ik->i", Xb, A_inv, Xb)
    h = np.clip(h, None, 1.0 - 1e-10)

    y_pred_loo = y - e / (1.0 - h)

    q2   = float(r2_score(y, y_pred_loo))
    mae  = float(mean_absolute_error(y, y_pred_loo))
    rmsd = float(np.sqrt(mean_squared_error(y, y_pred_loo)))
    return q2, mae, rmsd


def fit_and_evaluate_single_combination_regression(model, combination, r2_threshold=0.7, persist_result=True):
    try:
        colms = _normalize_combination_to_columns(combination)
        selected_features = model.features_df[colms]
    except:
        selected_features = model.features_df[_parse_tuple_string(combination)]

    X = selected_features.to_numpy()
    y = model.target_vector.to_numpy()

    # ── Thread-local worker ───────────────────────────────────────────────────
    # search_models() evaluates many combinations concurrently via
    # Parallel(backend="threading"). fit()/evaluate() mutate instance state
    # (self.theta, self.X_b_train, self.model, ...); sharing one `model` object
    # across threads races and corrupts that state — e.g. one thread's 4-feature
    # theta getting matmul'd against another thread's 3-feature design matrix,
    # which surfaces as "matmul: ... size N is different from M" crashes once
    # combos of mixed sizes run in the same parallel batch. Work on a shallow
    # per-call copy so each thread's fit state is isolated; large shared
    # read-only data (features_df, target_vector) is still referenced, not
    # copied, so this stays cheap.
    worker = copy.copy(model)
    if worker.model_type.lower() != "linear":
        worker.model = copy.deepcopy(model.model)

    # ── Fit both estimators ───────────────────────────────────────────────────
    # Ordinary linear regression uses the matrix-based fit/predict path.
    # Nonlinear/penalized sklearn estimators, such as lasso, still need their
    # estimator-specific fit implementation.
    t0 = time.time()
    worker.fit(X, y)
    if worker.model_type.lower() != "linear":
        worker.model.fit(X, y)
    fit_time = time.time() - t0

    # ── In-sample evaluation (R², adj-R²) ────────────────────────────────────
    t1 = time.time()
    worker._trained_features = list(combination)
    evaluation_results, y_pred = worker.evaluate(X, y)
    eval_time = time.time() - t1
    coefficients, intercepts = worker.get_coefficients_from_trained_estimator()

    # ── Q² LOO — always compute, unconditionally ──────────────────────────────
    # Linear models use the analytic hat-matrix LOO (no re-fitting loop).
    t3 = time.time()
    if worker.model_type.lower() == "linear":
        q2, mae, rmsd = _analytic_loo_linear(X, y)
    else:
        q2, mae, rmsd = worker.calculate_q2_and_mae(X, y, n_splits=1)
    evaluation_results['Q2']   = q2
    evaluation_results['MAE']  = mae
    evaluation_results['RMSD'] = rmsd
    q2_time = time.time() - t3

    result = {
        'combination':  combination,
        'scores':       evaluation_results,
        'intercept':    intercepts,
        'coefficients': coefficients,
        'model':        worker,
        'predictions':  y_pred,
    }

    if 'scores' in result and persist_result:
        return _persist_regression_result(result, r2_threshold)

    # Parallel path (persist_result=False): return only scalars.
    # Returning the model object and predictions through the IPC pipe
    # causes MemoryError on large searches with many workers.
    return {
        'combination': combination,
        'scores':      evaluation_results,
    }





class PlotModel:
    def __init__(self, model):
        self.model = model
    
    


class LinearRegressionModel:

    def __init__(
        self,
        csv_filepaths,
        process_method: str = "one csv",
        names_column: Optional[str] = None,
        y_value: str = "output",
        leave_out=None,
        min_features_num: int = 2,
        max_features_num: Optional[int] = None,
        n_splits: int = 5,
        metrics: Optional[List[str]] = None,
        return_coefficients: bool = False,
        model_type: str = "linear",   # "linear" or "lasso"
        alpha: float = 1.0,           # L2 / Lasso regularisation strength
        app=None,
        db_path: str = "results",
        scale: bool = True,
        seed: int = 42,
        drop_columns: Union[str, List[str]] = "id",
    ):
        self.random_state = seed
        self.setup_random_seeds()

        self.csv_filepaths       = csv_filepaths
        self.process_method      = process_method
        self.y_value             = y_value
        self.names_column        = names_column
        self.leave_out           = leave_out
        self.min_features_num    = min_features_num
        self.max_features_num    = max_features_num
        self.metrics             = metrics if metrics is not None else ["r2", "adj_r2", "neg_mean_absolute_error"]
        self.return_coefficients = return_coefficients
        self.model_type          = model_type
        self.alpha               = alpha
        self.app                 = app
        self.n_splits            = n_splits
        self.scale               = scale
        self.drop_columns        = drop_columns

        # Dataset name (used for directory / DB naming)
        try:
            self.name = os.path.splitext(
                os.path.basename(self.csv_filepaths.get("features_csv_filepath"))
            )[0]
        except Exception:
            self.name = "in_memory_dataset"

        # Directory structure & DB
        self.paths = prepare_run_dirs(
            base_dir="runs",
            dataset_name=self.name,
            y_value=self.y_value,
            tag=self.model_type,
        )
        self.db_path = resolve_db_path(db_path, self.name, self.paths.db)
        create_results_table(self.db_path)

        # Sklearn estimator
        model_type_lc = model_type.lower()
        if model_type_lc == "linear":
            self.model = LinearRegression()
        elif model_type_lc == "lasso":
            self.model = Lasso(alpha=alpha)
        else:
            raise ValueError(f"Invalid model_type {model_type!r}. Choose 'linear' or 'lasso'.")

        # Load data
        if csv_filepaths:
            if process_method == "one csv":
                self.process_features_csv(drop_columns=drop_columns)
            elif process_method == "two csvs":
                self.csv_filepaths = merge_features_and_target(
                    csv_filepaths.get("features_csv_filepath"), csv_filepaths.get("target_csv_filepath"), self.y_value)
                self.process_features_csv(drop_columns=drop_columns)

        # Scale features
        self.scaler              = StandardScaler()
        self.original_features_df = self.features_df.copy()
        self.features_df         = pd.DataFrame(
            self.scaler.fit_transform(self.features_df),
            columns=self.features_df.columns,
            index=self.features_df.index,
        )
        self.feature_names = self.features_df.columns.tolist()

        self.leave_out_samples(self.leave_out)
        self.determine_number_of_features()

        # Fit-state (populated by self.fit())
        self.theta              = None
        self.X_b_train          = None
        self.sigma2             = None
        self.XtX_inv            = None
        self.model_covariance_matrix  = None


    # ─────────────────────────────────────────────────────────────────────────
    # Setup / utility
    # ─────────────────────────────────────────────────────────────────────────
    def drop_rows(self, names):
        """
        Remove rows from the training set by molecule name (or list of names).

        Updates features_df, original_features_df, target_vector, and
        molecule_names in place, then re-fits the scaler on the remaining
        data so scaling parameters stay consistent.

        Parameters
        ----------
        names : str or list[str]
            One or more molecule names to drop (matched against molecule_names /
            the DataFrame index).

        Returns
        -------
        self  (for chaining)

        Example
        -------
        model.drop_rows("L23_2")
        model.drop_rows(["mol_3", "mol_7"])
        """
        if isinstance(names, str):
            names = [names]
        names = list(names)

        missing = [n for n in names if n not in self.molecule_names]
        if missing:
            raise ValueError(f"Names not found in training set: {missing}")

        keep = [n for n in self.molecule_names if n not in set(names)]

        # Drop from both scaled and unscaled feature tables
        self.original_features_df = self.original_features_df.loc[keep]
        self.features_df          = self.features_df.loc[keep]
        self.target_vector        = self.target_vector.loc[keep]
        self.molecule_names       = keep

        # Re-fit scaler on the reduced unscaled data so parameters are correct
        self.scaler = StandardScaler()
        self.features_df = pd.DataFrame(
            self.scaler.fit_transform(self.original_features_df),
            columns=self.original_features_df.columns,
            index=self.original_features_df.index,
        )

        # Clear any stale fit state
        self.theta             = None
        self.X_b_train         = None
        self.sigma2            = None
        self.XtX_inv           = None
        self.model_covariance_matrix = None

        print(f"Dropped {len(names)} row(s). Training set now has {len(keep)} molecules.")
        return self

    def load_results(self, sort_key: str = "adj_r2", ascending: bool = False):
        """
        Load all CSV results from the run's db directory into DataFrames,
        optionally sorting them by a metric (default: adjusted R²).

        Parameters
        ----------
        sort_key : str
            Column name to sort by (default 'adj_r2').
        ascending : bool
            Sort order. Default False = best (highest values first).

        Returns
        -------
        results_dict : dict[str, pd.DataFrame]
            Mapping {filename (no extension): DataFrame}.
        """
        results_dict = {}

        _drop_cols = ['threshold', 'model', 'predictions']
        for csv_file in Path(self.paths.db).glob("*.csv"):
            try:
                df = pd.read_csv(csv_file)
                df.drop(columns=_drop_cols, errors='ignore', inplace=True)
                if sort_key in df.columns:
                    df = df.sort_values(by=sort_key, ascending=ascending).reset_index(drop=True)
                else:
                    print(f"WARNING: Column '{sort_key}' not found in {csv_file.name}, returning unsorted.")
                results_dict[csv_file.stem] = df
            except Exception as e:
                print(f"ERROR: Failed to load {csv_file.name}: {e}")
        return results_dict



    def setup_random_seeds(self) -> None:
        """Seed Python and NumPy for reproducibility."""
        random.seed(self.random_state)
        np.random.seed(self.random_state)

    def check_linear_regression_assumptions(self):
        """Check linearity, homoscedasticity, and normality assumptions."""
        if not self.features_df.index.equals(self.target_vector.index):
            raise ValueError(
                "features_df and target_vector have mismatched indices; "
                "align them before calling check_linear_regression_assumptions()."
            )
        return check_linear_regression_assumptions(self.features_df, self.target_vector)

    def compute_multicollinearity(self, vif_threshold: float = 5.0) -> pd.DataFrame:
        """Return VIF table for all features, flagging those above vif_threshold."""
        vif_df = self._compute_vif(self.features_df)
        vif_df["high_vif"] = vif_df["VIF"] > vif_threshold
        return vif_df
    

    def process_features_csv(
        self,
        drop_columns: Optional[Union[str, List[str]]] = None,
        drop_smiles: bool = True,
    ) -> dict:
        """
        Load and process a features dataset.

        Accepts: CSV path (str), dict with 'features_csv_filepath', single-path
        list/tuple, or an in-memory DataFrame.

        Sets:
            self.molecule_names  : list[str]
            self.target_vector   : pd.Series  (indexed by molecule names)
            self.features_df     : pd.DataFrame  (numeric, deduplicated, indexed)
            self.features_list   : list[str]
        """

        # ---------- Load source ----------
        csv_source = getattr(self, "csv_filepaths", None)
        if csv_source is None:
            raise AttributeError(
                "`self.csv_filepaths` must be set to a CSV path (str), a dict, a pandas DataFrame, "
                "or a single-path list/tuple."
            )

        if isinstance(csv_source, (list, tuple)):
            if len(csv_source) == 1 and isinstance(csv_source[0], str):
                csv_source = csv_source[0]
            else:
                raise TypeError(
                    "`csv_filepaths` list/tuple must contain exactly one string path; "
                    f"got {csv_source!r}"
                )

        path_used = None

        if isinstance(csv_source, pd.DataFrame):
            df = csv_source.copy()
        elif isinstance(csv_source, dict):
            csv_path = csv_source.get("features_csv_filepath")
            if not csv_path or not os.path.exists(csv_path):
                raise FileNotFoundError(f"CSV not found at: {csv_path!r}")
            df = pd.read_csv(csv_path)
            path_used = os.path.abspath(csv_path)
        elif isinstance(csv_source, str):
            if not os.path.exists(csv_source):
                raise FileNotFoundError(f"CSV not found at: {csv_source!r}")
            df = pd.read_csv(csv_source)
            path_used = os.path.abspath(csv_source)
        else:
            raise TypeError("`self.csv_filepaths` must be a str, dict, DataFrame, or single-path list/tuple.")

        if df.empty:
            raise ValueError("The provided dataset is empty.")

        # ---------- Helpers ----------
        def _case_insensitive_find(col_name: str, columns: pd.Index):
            lc = str(col_name).lower()
            for c in columns:
                if str(c).lower() == lc:
                    return c
            return None

        def _resolve_many(names: List[str], columns: pd.Index) -> List[str]:
            resolved = []
            for n in names:
                hit = _case_insensitive_find(n, columns)
                if hit is not None:
                    resolved.append(hit)
            return list(dict.fromkeys(resolved))  # dedupe

        def _first_match(candidates: List[str], columns: pd.Index):
            for c in candidates:
                hit = _case_insensitive_find(c, columns)
                if hit is not None:
                    return hit
            return None

        def _clean_and_make_unique(names: List[str]) -> List[str]:
            # strip, normalize 'nan', and make unique with suffixes: name, name#2, name#3, ...
            base = []
            for x in names:
                s = "" if x is None else str(x)
                s = s.strip()
                if s.lower() == "nan":
                    s = ""
                base.append(s)

            result = []
            seen = {}
            for s in base:
                key = s if s else ""  # allow empty then disambiguate
                if key not in seen:
                    seen[key] = 1
                    result.append(key if key else "row")
                else:
                    seen[key] += 1
                    # append a simple counter
                    tag = f"{key if key else 'row'}#{seen[key]}"
                    result.append(tag)
            return result

        # ---------- Target & name resolution ----------
        y_value = getattr(self, "y_value", None)
        if y_value is None:
            raise AttributeError("`self.y_value` must be set to the target column name.")

        if drop_columns is None and hasattr(self, "drop_columns"):
            drop_columns = getattr(self, "drop_columns")

        if isinstance(drop_columns, str):
            drop_columns = [drop_columns]
        elif drop_columns is None:
            drop_columns = []

        provided_names_col = getattr(self, "names_column", None)
        name_candidates = ["name", "names", "molecule", "molecule_name", "mol", "compound", "sample", "id"]

        # Choose names column (honor provided_names_col first)
        names_col = None
        if provided_names_col:
            names_col = _case_insensitive_find(provided_names_col, df.columns)
            if names_col is None:
                raise KeyError(
                    f"names_column {provided_names_col!r} not found. Available: {list(df.columns)}"
                )
        else:
            if not isinstance(df.index, pd.RangeIndex):
                names_col = None  # use index
            else:
                names_col = _first_match(name_candidates, df.columns) or df.columns[0]

        target_col = _case_insensitive_find(y_value, df.columns)
        if target_col is None:
            raise KeyError(f"Target column {y_value!r} not found. Available: {list(df.columns)}")

        # ---------- Build features frame ----------
        cols_to_drop = {c for c in [names_col, target_col] if c is not None}
        features_df_raw = df.drop(columns=[c for c in cols_to_drop if c in df.columns])
        n_features_total = features_df_raw.shape[1]

        dropped_smiles_like_cols, dropped_user_columns = [], []

        # Drop SMILES-like columns
        if drop_smiles:
            smiles_candidates = ["smiles", "canonical_smiles", "random_smiles"]
            smiles_to_drop = _resolve_many(smiles_candidates, features_df_raw.columns)
            if smiles_to_drop:
                features_df_raw = features_df_raw.drop(columns=smiles_to_drop)
                dropped_smiles_like_cols.extend(smiles_to_drop)

        # Drop user-specified columns (support names or integer positions)
        if drop_columns:
            dc = list(drop_columns)
            by_pos = [c for c in dc if isinstance(c, (int, np.integer))]
            by_name = [c for c in dc if not isinstance(c, (int, np.integer))]

            resolved_user_cols = _resolve_many(by_name, features_df_raw.columns)

            if by_pos:
                pos_names = [
                    features_df_raw.columns[int(i)] for i in by_pos
                    if 0 <= int(i) < features_df_raw.shape[1]
                ]
                resolved_user_cols.extend(pos_names)

            resolved_user_cols = list(dict.fromkeys(resolved_user_cols))  # dedupe
            if resolved_user_cols:
                features_df_raw = features_df_raw.drop(columns=resolved_user_cols)
                dropped_user_columns.extend(resolved_user_cols)

        # Deduplicate columns
        features_df_raw = features_df_raw.loc[:, ~features_df_raw.columns.duplicated()]

        # ---------- Coerce to numeric ----------
        coerced = features_df_raw.apply(pd.to_numeric, errors="coerce")
        all_nan_cols = [c for c in coerced.columns if coerced[c].isna().all()]
        numeric_df = coerced.drop(columns=all_nan_cols, errors="ignore").select_dtypes(include=[np.number])
        dropped_non_numeric = [
            c for c in features_df_raw.columns if c not in numeric_df.columns and c not in all_nan_cols
        ]

        # ---------- Names / index handling (MAIN FIX) ----------
        if names_col is None:
            names_series = pd.Series(df.index, copy=False)
            names_label = "<index>"
        else:
            names_series = pd.Series(df[names_col], copy=False)
            names_label = names_col

        clean_unique_names = _clean_and_make_unique(names_series.tolist())

        # Apply the same name index to both X and y
        numeric_df = numeric_df.copy()
        numeric_df.index = clean_unique_names

        target_series_raw = pd.to_numeric(df[target_col], errors="coerce")
        # If everything became NaN, fallback to original dtype
        target_series = target_series_raw if not target_series_raw.isna().all() else df[target_col]
        target_series = pd.Series(target_series.values, index=clean_unique_names, name=target_col)

        # Optionally drop rows with NaN targets (keeps model flows sane)
        nan_target_mask = target_series.isna()
        n_nan_targets = int(nan_target_mask.sum())
        if n_nan_targets:
            target_series = target_series[~nan_target_mask]
            numeric_df = numeric_df.loc[target_series.index]

        # Final align & checks
        shared = numeric_df.index.intersection(target_series.index)
        numeric_df = numeric_df.loc[shared]
        target_series = target_series.loc[shared]

        # Sync public attributes
        self.molecule_names = list(numeric_df.index)
        self.target_vector = target_series
        self.features_df = numeric_df
        self.features_list = numeric_df.columns.tolist()

        # ---------- Summary ----------
        name_src_note = "from index" if names_col is None else f"from column '{names_col}'"
        summary = {
            "path": path_used,
            "names_column": names_label,
            "target_column": target_col,
            "n_rows": df.shape[0],
            "n_features_total": n_features_total,
            "n_features_kept_numeric": numeric_df.shape[1],
            "dropped_non_numeric_cols": dropped_non_numeric,
            "dropped_all_nan_cols": all_nan_cols,
            "dropped_smiles_like_cols": dropped_smiles_like_cols,
            "dropped_user_columns": dropped_user_columns,
            "dropped_rows_nan_target": n_nan_targets,
            "n_final_rows": numeric_df.shape[0],
        }

        # ---------- Output/log ----------
        msg = (
            f"Processed {'DataFrame' if path_used is None else os.path.basename(path_used)}\n"
            f"Molecule names: {name_src_note} | Target: {target_col}\n"
            f"Rows: {df.shape[0]} -> {summary['n_final_rows']} (dropped NaN targets: {n_nan_targets}) | "
            f"Features: {n_features_total} -> {numeric_df.shape[1]}"
        )
        for key, lst in {
            "SMILES-like": dropped_smiles_like_cols,
            "User-specified": dropped_user_columns,
            "Non-numeric": dropped_non_numeric,
            "All-NaN": all_nan_cols,
        }.items():
            if lst:
                msg += f"\nDropped {key} cols: {lst}"

        app = getattr(self, "app", None)
        if app:
            try:
                app.show_result(msg)
            except Exception:
                print(msg)
        else:
            print(msg)

        # Final sanity
        assert self.features_df.index.equals(self.target_vector.index), \
            "X and y must share identical indices after processing."

        return summary



    def compute_correlation(
        self,
        correlation_threshold: float = 0.80,
        vif_threshold: Optional[float] = None,
    ) -> dict:
        """
        Analyse multicollinearity via correlation and (optionally) greedy VIF pruning.

        Parameters
        ----------
        correlation_threshold : float
            Absolute correlation cutoff to flag feature pairs.
        vif_threshold : float or None
            If set, greedily drop the highest-VIF feature until all VIF ≤ threshold.

        Returns
        -------
        dict with keys: corr_matrix, high_corr_features, sub_corr,
                        vif_table_before, vif_table_after,
                        dropped_features, kept_features, message
        """
        app = getattr(self, "app", None)

        num_df = self.features_df.select_dtypes(include=[np.number])
        self.corr_matrix = num_df.corr()
        high_corr_features = self._get_highly_correlated_features(
            self.corr_matrix, threshold=correlation_threshold
        )
        msgs: List[str] = []
        dropped_features: List[str] = []
        sub_corr = vif_table_before = vif_table_after = None

        if high_corr_features:
            msgs.append(
                f"Correlation Report: features with |r| > {correlation_threshold}:\n"
                f"{sorted(high_corr_features)}"
            )

            # Optional heatmap (requires GUI app)
            if app is not None:
                from tkinter import messagebox
                if messagebox.askyesno("Visualize?",
                                       "Show heatmap of correlated features?"):
                    sub_corr = self.corr_matrix.loc[
                        sorted(high_corr_features), sorted(high_corr_features)
                    ]
                    interactive_corr_heatmap(
                        sub_corr,
                        title=f"Correlation Heatmap (|r| > {correlation_threshold})"
                    )

            # Greedy VIF pruning
            apply_vif = vif_threshold is not None
            if apply_vif and app is not None:
                from tkinter import messagebox
                apply_vif = messagebox.askyesno(
                    "VIF Pruning?",
                    f"Apply VIF pruning (threshold={vif_threshold})?"
                )

            if apply_vif:
                vif_cols = [c for c in num_df.columns if c in high_corr_features]
                if len(vif_cols) >= 2:
                    vif_table_before = self._compute_vif(num_df[vif_cols])
                    msgs.append(f"VIF before:\n{vif_table_before.head(10).to_string(index=False)}")
                    keep_cols = vif_cols[:]
                    while len(keep_cols) > 1:
                        vif_tab = self._compute_vif(num_df[keep_cols])
                        max_vif = float(vif_tab["VIF"].max())
                        if np.isnan(max_vif) or max_vif <= vif_threshold:
                            vif_table_after = vif_tab
                            break
                        to_drop = str(vif_tab.iloc[0]["variable"])
                        keep_cols.remove(to_drop)
                        dropped_features.append(to_drop)
                else:
                    msgs.append("VIF pruning skipped — fewer than 2 correlated numeric features.")

        else:
            msgs.append("No feature pairs exceeded the correlation threshold.")

        final_msg = "\n".join(msgs)
        if app is not None:
            app.show_result(final_msg)
        else:
            print(final_msg)

        return {
            "corr_matrix":        self.corr_matrix,
            "high_corr_features": high_corr_features,
            "sub_corr":           sub_corr,
            "vif_table_before":   vif_table_before,
            "vif_table_after":    vif_table_after,
            "dropped_features":   dropped_features,
            "kept_features":      getattr(self, "features_list", list(self.features_df.columns)),
            "message":            final_msg,
        }



    def _compute_vif(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute VIF for each numeric column. Returns DataFrame[variable, VIF]."""
        df_num = (
            df.select_dtypes(include=[np.number])
            .replace([np.inf, -np.inf], np.nan)
            .dropna(axis=1, how="any")
            .dropna(axis=0, how="any")
        )
        if df_num.shape[1] < 2:
            return pd.DataFrame({"variable": df_num.columns, "VIF": [np.nan] * df_num.shape[1]})

        rank = np.linalg.matrix_rank(df_num.values)
        if rank < df_num.shape[1]:
            print(f"Warning: linear dependence detected (rank {rank} < {df_num.shape[1]} features).")

        X_scaled = StandardScaler().fit_transform(df_num)
        vif_vals = [variance_inflation_factor(X_scaled, i) for i in range(df_num.shape[1])]
        return pd.DataFrame({"variable": df_num.columns, "VIF": vif_vals})

    def _flag_high_vif(self, results: pd.DataFrame, vif_threshold: float = 10.0) -> pd.DataFrame:
        """
        For each row in `results` where Q2 was actually calculated (not NaN/-inf),
        compute the max VIF among that combination's features and flag high_vif
        when it exceeds `vif_threshold`. Rows without a calculated Q2 get NaN/False.
        """
        if results is None or results.empty or "combination" not in results.columns:
            return results

        results = results.copy()
        q2 = _extract_q2(results)
        has_q2 = q2.notna() & (q2 != float("-inf")) if q2 is not None else pd.Series(False, index=results.index)

        max_vifs = []
        high_flags = []
        for idx, row in results.iterrows():
            if not has_q2.loc[idx]:
                max_vifs.append(np.nan)
                high_flags.append(False)
                continue
            try:
                cols = _parse_tuple_string(row["combination"])
                cols = [c for c in cols if c in self.features_df.columns]
                if len(cols) < 2:
                    max_vifs.append(np.nan)
                    high_flags.append(False)
                    continue
                vif_df = self._compute_vif(self.features_df[cols])
                max_vif = float(vif_df["VIF"].max()) if not vif_df.empty else np.nan
            except Exception:
                max_vif = np.nan
            max_vifs.append(max_vif)
            high_flags.append(bool(max_vif > vif_threshold) if pd.notna(max_vif) else False)

        results["max_vif"] = max_vifs
        results["high_vif"] = high_flags
        return results

    def _get_highly_correlated_features(self, corr_matrix: pd.DataFrame, threshold: float = 0.8) -> set:
        """Return the set of feature names involved in any pair with |r| > threshold."""
        arr  = corr_matrix.abs().values
        cols = corr_matrix.columns
        mask = np.triu(arr > threshold, k=1)           # upper triangle, exclude diagonal
        i_idx, j_idx = np.where(mask)
        return {cols[i] for i in i_idx} | {cols[j] for j in j_idx}


    def process_target_csv(self, csv_filepath):
        target_vector_unordered = pd.read_csv(csv_filepath)[self.y_value]
        
        self.target_vector = target_vector_unordered.loc[self.molecule_names]
        
    def leave_out_samples(self, leave_out=None, keep_only=False):
        """
        Remove or retain specific rows from features_df/target_vector based on indices or molecule names.
        Stores the removed (or retained) samples into:
            - self.leftout_samples
            - self.leftout_target_vector
            - self.molecule_names_predict

        Parameters:
            leave_out (int, str, list[int|str]): Indices or molecule names to leave out (or to keep if keep_only=True).
            keep_only (bool): If True, only the given indices/names are kept, all others are left out.
        """
       
        if leave_out is None:
            return

        # Normalize input to a list
        if isinstance(leave_out, (int, np.integer, str)):
            leave_out = [leave_out]
        else:
            leave_out = list(leave_out)

        # Determine indices from mixed inputs (names or indices)
        indices = []
        name_to_index = {name: idx for idx, name in enumerate(self.molecule_names)}

        for item in leave_out:
            if isinstance(item, (int, np.integer)):
                indices.append(int(item))

            elif isinstance(item, str):
                # try exact match first
                if item in name_to_index:
                    indices.append(name_to_index[item])
                else:
                    # try to parse as int
                    try:
                        idx_val = int(item)
                        indices.append(idx_val)
                    except ValueError:
                        raise ValueError(
                            f"Leave-out entry '{item}' is not a known molecule name "
                            f"and cannot be parsed as an integer index."
                        )

            else:
                raise ValueError("Items in leave_out must be integers or strings.")

        selected_features = self.features_df.iloc[indices].copy()
        selected_target = self.target_vector.iloc[indices].copy()
        selected_names = [self.molecule_names[i] for i in indices]

        # Ensure we have the full set of columns in the right order
        if hasattr(self, 'feature_names'):
            selected_features = selected_features.reindex(
                columns=self.feature_names,
                fill_value=0
            )

        # locate index labels for dropping
        feat_labels = self.features_df.index.take(indices)

        # --- decide how to address the target_vector (labels vs positions) ---
        use_label_drop_for_target = False
        targ_labels = None
        if self.target_vector.index.equals(self.features_df.index):
            # same index as features → safe to drop by the same labels
            targ_labels = feat_labels
            use_label_drop_for_target = True
        elif hasattr(self, "molecule_names") and len(self.molecule_names) == len(self.target_vector):
            # try mapping via molecule_names if they mirror the target order
            cand = pd.Index(self.molecule_names).take(indices)
            if pd.Index(cand).isin(self.target_vector.index).all():
                targ_labels = cand
                use_label_drop_for_target = True
        # else: we will drop target by positions

        # cache for downstream / debugging
        self.idx_labels = feat_labels
        self.indices_predict = indices

        if keep_only:
            # Everything NOT in indices becomes the "left out" set
            self.leftout_samples = self.features_df.drop(index=feat_labels).copy()
            if use_label_drop_for_target:
                self.leftout_target_vector = self.target_vector.drop(index=targ_labels).copy()
            else:
                self.leftout_target_vector = self.target_vector.drop(self.target_vector.index[indices]).copy()
            self.molecule_names_predict = [n for i, n in enumerate(self.molecule_names) if i not in indices]

            # Keep only the selected rows as the main set (already sliced earlier)
            self.features_df   = selected_features
            self.target_vector = selected_target
            self.molecule_names = selected_names
        else:
            # The selected rows become the "left out" set
            self.leftout_samples        = selected_features
            self.leftout_target_vector  = selected_target
            self.molecule_names_predict = selected_names

            # Drop them from the main set
            self.features_df = self.features_df.drop(index=feat_labels)
            if use_label_drop_for_target:
                self.target_vector = self.target_vector.drop(index=targ_labels)
            else:
                self.target_vector = self.target_vector.drop(self.target_vector.index[indices])
            self.molecule_names = [n for i, n in enumerate(self.molecule_names) if i not in indices]

    def determine_number_of_features(self):
        total_features_num = self.features_df.shape[0]
        self.max_features_num = set_max_features_limit(total_features_num, self.max_features_num)

    def get_feature_combinations(self):
        self.features_combinations = list(get_feature_combinations(self.features_list, self.min_features_num, self.max_features_num))


    def _add_leftout_mae_to_results(self, results: pd.DataFrame) -> pd.DataFrame:
        """
        Compute MAE of left-out predictions for each model in results.
        Fits a fresh model on training data for each combination and evaluates on left-out samples.
        
        Parameters
        ----------
        results : pd.DataFrame
            Results dataframe with 'combination' column.
            
        Returns
        -------
        pd.DataFrame
            Updated results with 'leftout_mae' column added.
        """
        if results.empty or not hasattr(self, 'leftout_samples') or self.leftout_samples is None:
            return results
        
        results = results.copy()
        leftout_mae_values = []
        
        for idx, row in results.iterrows():
            try:
                combination = row.get('combination')
                
                if combination is None:
                    leftout_mae_values.append(np.nan)
                    continue
                
                # Extract features for both training and left-out samples
                try:
                    colms = _normalize_combination_to_columns(combination)
                    train_X = self.features_df[colms]
                    leftout_X = self.leftout_samples[colms]
                except:
                    colms = _parse_tuple_string(combination)
                    train_X = self.features_df[colms]
                    leftout_X = self.leftout_samples[colms]
                
                # Get target values (ensure they're numpy arrays)
                train_y = self.target_vector.values if hasattr(self.target_vector, 'values') else np.asarray(self.target_vector)
                leftout_y = self.leftout_target_vector.values if hasattr(self.leftout_target_vector, 'values') else np.asarray(self.leftout_target_vector)
                
                # Ensure X data is numpy arrays with correct shape
                train_X_arr = train_X.values if hasattr(train_X, 'values') else np.asarray(train_X)
                leftout_X_arr = leftout_X.values if hasattr(leftout_X, 'values') else np.asarray(leftout_X)
                
                # Validation: check that shapes align
                if train_X_arr.shape[0] != train_y.shape[0]:
                    raise ValueError(
                        f"Training data mismatch: X has {train_X_arr.shape[0]} samples, "
                        f"y has {train_y.shape[0]} samples"
                    )
                if leftout_X_arr.shape[0] != leftout_y.shape[0]:
                    raise ValueError(
                        f"Left-out data mismatch: X has {leftout_X_arr.shape[0]} samples, "
                        f"y has {leftout_y.shape[0]} samples"
                    )
                
                if self.model_type.lower() == 'linear':
                    leftout_predictions = _linear_matrix_fit_predict(
                        train_X_arr,
                        train_y,
                        leftout_X_arr,
                        alpha=1e-5,
                        scale=False,
                    )
                elif self.model_type.lower() == 'lasso':
                    temp_model = Lasso(alpha=self.alpha)
                    temp_model.fit(train_X_arr, train_y)
                    leftout_predictions = temp_model.predict(leftout_X_arr)
                else:
                    raise ValueError(f"Unknown model_type: {self.model_type}")
                
                # Calculate MAE
                mae = mean_absolute_error(leftout_y, leftout_predictions)
                leftout_mae_values.append(mae)
                
            except Exception as e:
                print(f"WARNING: Failed to compute left-out MAE for combination {row.get('combination')}: {e}")
                leftout_mae_values.append(np.nan)
        
        results['leftout_mae'] = leftout_mae_values
        return results

    def plot_shap_values(self, X, feature_names=None, target_name="output",
                         n_top_features: int = 10, plot: bool = True, dir=None):
        """Plot SHAP feature-importance values for the fitted model."""
        sample_names = getattr(self, "molecule_names", None)
        save_dir     = getattr(getattr(self, "paths", None), "png_dir", dir)
        analyze_shap_values(
            self.model, X,
            feature_names=feature_names,
            sample_names=sample_names,
            target_name=target_name,
            n_top_features=n_top_features,
            plot=plot,
            dir=save_dir,
        )

    def calculate_q2_and_mae(
        self,
        X,
        y,
        n_splits=None,
        test_size=0.1,
        random_state=42,
        n_iterations=5,
        exclude_leftout: bool = True,
        plot: bool = False,
    ):
        """
        Cross-validation metrics: Q², MAE, RMSD.
        Optionally plots out-of-fold predictions instead of printing per-fold stats.
        """
        import numpy as np
        import matplotlib.pyplot as plt
        from sklearn.model_selection import train_test_split, LeaveOneOut, KFold
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
        from sklearn.base import clone
        from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

        X = np.asarray(X)
        y = np.asarray(y).ravel()

        # Exclude reserved molecules
        if exclude_leftout and hasattr(self, "molecule_names_predict") and hasattr(self, "molecule_names"):
            try:
                keep_mask = np.array(
                    [name not in set(self.molecule_names_predict) for name in self.molecule_names],
                    dtype=bool,
                )
                if keep_mask.shape[0] == X.shape[0]:
                    X, y = X[keep_mask], y[keep_mask]
            except Exception:
                pass

        n = X.shape[0]
        if X.ndim != 2 or y.ndim != 1 or y.shape[0] != n:
            raise ValueError("X must be 2D and y must be 1D with matching rows")

        def _fit_predict(train_idx, test_idx):
            if self.model_type.lower() == "linear":
                return _linear_matrix_fit_predict(
                    X[train_idx],
                    y[train_idx],
                    X[test_idx],
                    alpha=1e-5,
                    scale=True,
                )
            pipe = make_pipeline(StandardScaler(), clone(self.model))
            pipe.fit(X[train_idx], y[train_idx])
            return pipe.predict(X[test_idx])

        q2_list, mae_list, rmsd_list = [], [], []
        all_preds, all_true, all_folds = [], [], []
        n_iter = max(1, int(n_iterations))
        base_seed = None if random_state is None else int(random_state)

        # -----------------
        # HOLDOUT
        # -----------------
        if n_splits is None or n_splits == 0:
            for it in range(n_iter):
                rs = None if base_seed is None else base_seed + it
                tr, te = train_test_split(np.arange(n), test_size=test_size, shuffle=True, random_state=rs)
                y_hat = _fit_predict(tr, te)
                q2 = float(r2_score(y[te], y_hat)) if len(te) >= 2 else float("nan")
                mae = float(mean_absolute_error(y[te], y_hat))
                rmsd = float(np.sqrt(mean_squared_error(y[te], y_hat)))
                q2_list.append(q2); mae_list.append(mae); rmsd_list.append(rmsd)
                all_preds.extend(y_hat); all_true.extend(y[te]); all_folds.extend([it]*len(te))

        # -----------------
        # LOO
        # -----------------
        elif n_splits == 1:
            loo = LeaveOneOut()
            y_pred_oof = np.empty_like(y, dtype=float)
            for tr, te in loo.split(X):
                y_pred_oof[te] = _fit_predict(tr, te)
            q2 = float(r2_score(y, y_pred_oof))
            mae = float(mean_absolute_error(y, y_pred_oof))
            rmsd = float(np.sqrt(mean_squared_error(y, y_pred_oof)))
            q2_list = [q2]; mae_list = [mae]; rmsd_list = [rmsd]
            all_preds = list(y_pred_oof); all_true = list(y); all_folds = list(range(n))

        # -----------------
        # K-FOLD
        # -----------------
        else:
            k = int(min(max(2, n // 2), int(n_splits)))
            for it in range(n_iter):
                rs = None if base_seed is None else base_seed + it
                kf = KFold(n_splits=k, shuffle=True, random_state=rs)
                y_pred_oof = np.empty_like(y, dtype=float)
                for fold, (tr, te) in enumerate(kf.split(X)):
                    y_hat = _fit_predict(tr, te)
                    y_pred_oof[te] = y_hat
                    all_preds.extend(y_hat); all_true.extend(y[te]); all_folds.extend([fold]*len(te))
                q2 = float(r2_score(y, y_pred_oof))
                mae = float(mean_absolute_error(y, y_pred_oof))
                rmsd = float(np.sqrt(mean_squared_error(y, y_pred_oof)))
                q2_list.append(q2); mae_list.append(mae); rmsd_list.append(rmsd)

        # --- Final metrics ---
        q2 = float(np.nanmean(q2_list))
        mae = float(np.mean(mae_list))
        rmsd = float(np.mean(rmsd_list))

        # -----------------
        # Plot diagnostics
        # -----------------
        if plot and len(all_preds) > 0:
            all_preds = np.array(all_preds)
            all_true = np.array(all_true)
            all_folds = np.array(all_folds)
            plt.figure(figsize=(6,6))
            scatter = plt.scatter(all_true, all_preds, c=all_folds, cmap="tab10", s=60, alpha=0.8, edgecolor="k")
            lims = [min(all_true.min(), all_preds.min())-0.5, max(all_true.max(), all_preds.max())+0.5]
            plt.plot(lims, lims, "r--", lw=2, label="Ideal")
            plt.xlabel("True values")
            plt.ylabel("Predicted values (OOF)")
            plt.title(f"CV predictions\nQ²={q2:.3f}, MAE={mae:.3f}, RMSD={rmsd:.3f}")
            plt.legend()
            plt.colorbar(scatter, label="Fold")
            plt.tight_layout()
            plt.show()

        return q2, mae, rmsd




    def fit(self, X, y, alpha: float = 1e-5):
        """
        Fit with bias term and ridge regularisation (no penalty on bias).
        Precomputes σ², (XᵀX)⁻¹, and the coefficient covariance matrix.
        """
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float).ravel()
        X_b = np.c_[np.ones((X.shape[0], 1)), X]
        self.X_b_train = X_b

        # Ridge normal equation
        A = X_b.T @ X_b
        I = np.eye(A.shape[0])
        I[0, 0] = 0.0   # no penalty on bias
        rhs = X_b.T @ y
        try:
            self.theta = np.linalg.solve(A + alpha * I, rhs)
        except np.linalg.LinAlgError:
            self.theta = np.linalg.pinv(A + alpha * I) @ rhs

        # Residual variance σ² = SSE / (n − p)
        residuals   = y - X_b @ self.theta
        n, p        = X_b.shape
        self.sigma2 = float((residuals ** 2).sum() / max(n - p, 1))

        # (XᵀX)⁻¹ used for prediction intervals and covariance
        self.XtX_inv = np.linalg.pinv(X_b.T @ X_b)

        # Coefficient covariance matrix (σ² · (XᵀX)⁻¹)
        self.model_covariance_matrix = self.sigma2 * self.XtX_inv

        return self

    def predict(self, X, return_interval: bool = False, cl: float = 0.95):
        """
        Point predictions, and optionally (lower, upper) prediction intervals.

        Parameters
        ----------
        X : array-like, shape (n, p)
        return_interval : bool
        cl : float  confidence level for the interval
        """
        X_b  = np.c_[np.ones((X.shape[0], 1)), X]
        yhat = X_b @ self.theta

        if not return_interval:
            return yhat

        dof  = self.X_b_train.shape[0] - self.X_b_train.shape[1]
        tval = t.ppf(1 - (1 - cl) / 2, dof)

        # Prediction-interval SE: sqrt(σ² · (1 + xᵀ(XᵀX)⁻¹x))
        inside  = 1 + np.sum((X_b @ self.XtX_inv) * X_b, axis=1)
        se_pred = np.sqrt(self.sigma2 * inside)

        return yhat, yhat - tval * se_pred, yhat + tval * se_pred



    def predict_for_leftout(
    self,
    X: pd.DataFrame,
    y: Optional[object] = None,
    X_train: Optional[pd.DataFrame] = None,
    y_train: Optional[object] = None,
    *,
    calc_interval: bool = False,
    confidence_level: float = 0.95,
):
        """
        Predict on left-out samples using an already-fitted or freshly refitted model.

        Parameters
        ----------
        X : pd.DataFrame or np.ndarray
            Left-out feature matrix.
        y : array-like, optional
            Left-out target values.
        X_train, y_train : array-like, optional
            If provided, the model will be (re)fit on this data before predicting.
        calc_interval : bool, default False
            If True and the model/predict wrapper supports it, compute prediction intervals.
        confidence_level : float, default 0.95
            Confidence level for prediction intervals (if supported).

        Returns
        -------
        preds or (preds, errors) or (preds, lower, upper) or (preds, lower, upper, errors)
        """

        # ----------------------- 1) Ensure fitted (or fit now) -----------------------
        # NOTE: self.model (the plain sklearn estimator) and self.theta (this class's
        # own ridge-regression coefficients, set by self.fit() and consumed by
        # self.predict() below) are two separate fitted states. Refitting only
        # self.model here used to leave self.theta stale/unset, which made the
        # self.predict() call further down crash with "matmul: ... has 0 dimensions"
        # whenever self.theta was still None - so both must be kept in sync.
        must_refit = (X_train is not None) and (y_train is not None)
        if must_refit:
            Xt = X_train.to_numpy() if hasattr(X_train, "to_numpy") else np.asarray(X_train)
            yt = np.asarray(y_train).ravel()
            self.fit(Xt, yt)                # refit theta (used by self.predict())
            if self.model_type.lower() != "linear":
                self.model.fit(Xt, yt)      # keep the sklearn estimator (e.g. Lasso) in sync too
            # Optionally keep for future checks
            self.X_, self.y_ = Xt, yt
        else:
            try:
                check_is_fitted(self.model)
            except Exception:
                # try cached training data, else error
                if hasattr(self, "X_") and hasattr(self, "y_"):
                    self.model.fit(self.X_, self.y_)
                else:
                    raise NotFittedError(
                        "This estimator is not fitted yet. Provide X_train,y_train or fit before calling predict_for_leftout."
                    )
            if self.theta is None:
                # self.model was fitted above (or already fit), but self.predict()
                # below needs self.theta specifically - refit it too if possible.
                if hasattr(self, "X_") and hasattr(self, "y_"):
                    self.fit(self.X_, self.y_)
                else:
                    raise NotFittedError(
                        "This estimator's ridge coefficients (theta) are not fitted yet. "
                        "Provide X_train,y_train or call fit() before calling predict_for_leftout."
                    )

        # ----------------------- 2) Align columns/schema -----------------------
        X_is_df = isinstance(X, pd.DataFrame)
        if X_is_df and isinstance(X_train, pd.DataFrame):
            # reindex to training columns; missing -> NaN
            X_aligned = X.reindex(columns=X_train.columns)
        elif X_is_df and hasattr(self, "features_df") and isinstance(self.features_df, pd.DataFrame):
            X_aligned = X.reindex(columns=self.features_df.columns)
        else:
            X_aligned = X

        # handle missing columns (NaNs after reindex)
        if isinstance(X_aligned, pd.DataFrame):
            if X_aligned.isna().any().any():
                # Simple, safe default: fill NaNs with 0.0
                # (swap to column means if you store them, e.g., self.feature_means_)
                X_aligned = X_aligned.fillna(0.0)
            X_arr = X_aligned.to_numpy()
        else:
            X_arr = np.asarray(X_aligned)

        # Ensure 2D
        if X_arr.ndim == 1:
            X_arr = X_arr.reshape(-1, 1)

        # ----------------------- 3) Apply scaler if present -----------------------
        if hasattr(self, "scaler_") and self.scaler_ is not None:
            try:
                X_arr = self.scaler_.transform(X_arr)
            except Exception:
                # if scaler fails (shape mismatch), skip silently or raise as you prefer
                pass

        # ----------------------- 4) Predict (with or without intervals) -----------------------
        # If your model class has a wrapper predict(X, return_interval=..., cl=...),
        # prefer that; else fall back to sklearn's .predict
        preds = lower = upper = None

        if calc_interval and hasattr(self, "predict"):
            # Your wrapper signature uses 'cl' in earlier code; keep both for safety
            try:
                preds, lower, upper = self.predict(X_arr, return_interval=True, cl=confidence_level)
            except TypeError:
                preds, lower, upper = self.predict(X_arr, return_interval=True, confidence_level=confidence_level)
        elif hasattr(self, "predict"):
            # If your wrapper exists but intervals not requested
            out = self.predict(X_arr, return_interval=False)
            preds = out[0] if isinstance(out, (tuple, list)) else out
        else:
            # plain sklearn estimator
            preds = self.model.predict(X_arr)

        preds = np.asarray(preds).ravel()

        # ----------------------- 5) Errors (if y given) -----------------------
        errors = None
        if y is not None:
            y_arr = y.values.ravel() if isinstance(y, (pd.Series, pd.DataFrame)) else np.asarray(y).ravel()
            if preds.shape != y_arr.shape:
                raise ValueError(f"Prediction shape {preds.shape} != y shape {y_arr.shape}")
            errors = y_arr - preds

        # ----------------------- 6) Return -----------------------
        if calc_interval and (lower is not None) and (upper is not None):
            lower = np.asarray(lower).ravel()
            upper = np.asarray(upper).ravel()
            return (preds, lower, upper, errors) if errors is not None else (preds, lower, upper)
        else:
            return (preds, errors) if errors is not None else preds

        
    def evaluate(self, X, y) -> Tuple[dict, np.ndarray]:
        """Evaluate on (X, y) using all configured metrics.

        Always computes both R² and adjusted R² regardless of self.metrics so
        downstream code that keys on either is never surprised.
        """
        if self.model_type.lower() == "linear" and self.theta is not None:
            y_pred = self.predict(X)
        else:
            y_pred = self.model.predict(X)
        n, p   = X.shape[0], X.shape[1]
        results = {}

        r2 = float(r2_score(y, y_pred))
        results["r2"]     = r2
        results["adj_r2"] = _adjusted_r2(r2, n, p)

        if "neg_mean_absolute_error" in self.metrics or "mae" in self.metrics:
            results["mae"] = float(mean_absolute_error(y, y_pred))
        if "rmse" in self.metrics:
            results["rmse"] = float(np.sqrt(mean_squared_error(y, y_pred)))
        return results, y_pred

    def cross_validate(self, X, y) -> Tuple[np.ndarray, dict]:
        """K-fold cross-validation; returns (oof_predictions, scores_dict)."""
        cv     = KFold(n_splits=self.n_splits, shuffle=True, random_state=self.random_state)
        y_pred = cross_val_predict(self.model, X, y, cv=cv)
        scores = cross_validate(self.model, X, y, cv=cv, scoring=self.metrics)
        return y_pred, scores

    def get_coefficients_from_trained_estimator(self) -> Tuple[float, np.ndarray]:
        """Return (intercept, coefficients) from the fitted theta vector."""
        return float(self.theta[0]), self.theta[1:]

    def get_covariance_matrix(self, features: List[str]) -> pd.DataFrame:
        """
        Return a coefficient summary table with Estimate, Std. Error, t value, p value.
        Requires self.fit() to have been called first.
        """
        if self.model_covariance_matrix is None:
            raise RuntimeError("Call fit() before get_covariance_matrix().")

        intercept    = float(self.theta[0])
        coefficients = self.theta[1:]
        std_errors   = np.sqrt(np.diag(self.model_covariance_matrix))
        dof          = len(self.target_vector) - len(coefficients) - 1

        t_vals_coef  = coefficients / std_errors[1:]
        p_vals_coef  = 2 * (1 - stats.t.cdf(np.abs(t_vals_coef), df=dof))

        t_val_int    = intercept / std_errors[0]
        p_val_int    = 2 * (1 - stats.t.cdf(np.abs(t_val_int), df=dof))

        coef_df = pd.DataFrame({
            "Estimate":   np.concatenate([[intercept],   coefficients]),
            "Std. Error": std_errors,
            "t value":    np.concatenate([[t_val_int],   t_vals_coef]),
            "p value":    np.concatenate([[p_val_int],   p_vals_coef]),
        }, index=["(Intercept)"] + features)

        return coef_df



    def search_models(
        self,
        top_n: int = 20,
        n_jobs: int = -1,
        threshold: float = 0.7,
        bool_parallel: bool = False,
        required_features: Optional[Iterable[str]] = None,
        min_models_to_keep: int = 5,
        threshold_step: float = 0.05,
        min_threshold: float = 0.2,
        vif_threshold: float = 10.0,
    ):
        """
        Fit and evaluate feature combinations with an optional set of required features.
        Threshold is relaxed automatically if all Q2 are -inf.

        Adds a `high_vif` flag (and `max_vif` value) for every model whose Q2 was
        actually calculated, marking combinations whose max feature VIF exceeds
        `vif_threshold` (default 10).

        Returns
        -------
        pd.DataFrame
            Top results, sorted by q2 then r2, limited to top_n.
        """
        app = self.app

        # Generate combinations up-front
        self.get_feature_combinations()
        all_combos: List[Tuple[str, ...]] = list(self.features_combinations or [])
        req: Set[str] = set(required_features or [])

        # Decide parallel
        cpu_count = multiprocessing.cpu_count()
        effective_jobs = _resolve_effective_jobs(
            n_jobs=n_jobs,
            bool_parallel=bool_parallel,
            cpu_count=cpu_count,
        )
        print(f"Using {effective_jobs} jobs for evaluation. Found {cpu_count} cores.")

        # Load existing results (avoid re-doing work)
        try:
            existing_results = _to_df(load_results_from_db(self.db_path))
            print(f"Loaded {len(existing_results)} existing results from DB.")
            # normalize combination keys
            if "combination" in existing_results.columns:
                done_combos = set(map(str, existing_results["combination"].tolist()))
            else:
                done_combos = set()
        except Exception:
            existing_results = pd.DataFrame()
            done_combos = set()

        def _filter_new_combos() -> List[Tuple[str, ...]]:
            """Not in DB and contains all required features."""
            out = []
            for combo in all_combos:
                key = str(combo)
                if key in done_combos:
                    continue
                if req and not req.issubset(set(combo)):
                    continue
                out.append(combo)
            return out

        def _evaluate_block(threshold: float) -> pd.DataFrame:
            """
            Evaluate only combos not yet done & satisfying required features.
            If parallel, assume each call writes to DB; then reload.
            """
            results_df = existing_results.copy()

            combos_to_run = _filter_new_combos()
            print(f"Combos to run: {len(combos_to_run)}, done_combos: {len(done_combos)}")
            if not combos_to_run:
                print(f"No new combinations to evaluate at threshold {threshold:.3f}.")
                return results_df

            print(f"Evaluating {len(combos_to_run)} new combos with R2 >= {threshold:.3f}...")
            if effective_jobs == 1:
                new_results = []
                for combo in tqdm(combos_to_run, desc=f"Threshold {threshold:.3f} (single-core)"):
                    try:
                        res = fit_and_evaluate_single_combination_regression(self, combo, threshold)
                        new_results.append(res)
                    except Exception:
                        # swallow or log; keep going
                        pass
                new_df = _to_df(new_results)
                results_df = pd.concat([results_df, new_df], ignore_index=True) if not new_df.empty else results_df
               
            else:
                # Always use threading:
                #   - NumPy/SciPy release the GIL during matrix ops, so
                #     threads achieve real CPU parallelism.
                #   - Zero process-spawn cost (no 2-10 s re-import per worker).
                #   - Results stay in shared memory — no IPC pickle, no
                #     MemoryError even with millions of combinations.
                #   - loky (process pool) is only better when the worker
                #     function runs pure-Python loops that hold the GIL, which
                #     is not the case here.
                print(f"Parallel backend: threading  "
                      f"(NumPy releases GIL — {effective_jobs} threads)")
                new_results = Parallel(
                    n_jobs=effective_jobs,
                    backend="threading",
                )(
                    delayed(fit_and_evaluate_single_combination_regression)(
                        self,
                        combo,
                        threshold,
                        False,
                    )
                    for combo in tqdm(combos_to_run, desc=f"Threshold {threshold:.3f} (threads)")
                )
                # Build flat result dicts and batch-write to DB/CSV in one transaction.
                # Workers return only {combination, scores} — no model object or predictions.
                csv_path = self.db_path.replace(".db", ".csv")
                flat_results = []
                for result in new_results:
                    if not isinstance(result, dict) or "scores" not in result:
                        continue
                    scores = result["scores"]
                    flat = {
                        "combination": str(result["combination"]),
                        "r2":          scores.get("r2",    float("-inf")),
                        "adj_r2":      scores.get("adj_r2", float("-inf")),
                        "q2":          scores.get("Q2",    float("-inf")),
                        "mae":         scores.get("MAE",   float("inf")),
                        "rmsd":        scores.get("RMSD",  float("inf")),
                        "threshold":   threshold,
                        "model":       "",
                        "predictions": "",
                    }
                    flat_results.append(flat)

                if flat_results:
                    try:
                        batch_insert_results_into_db_regression(self.db_path, flat_results, csv_path)
                    except Exception as _e:
                        print(f"Batch DB insert failed ({_e}); falling back to per-row inserts.")
                        for r in flat_results:
                            try:
                                insert_result_into_db_regression(
                                    db_path=self.db_path,
                                    combination=r["combination"],
                                    r2=r["r2"], adj_r2=r["adj_r2"],
                                    q2=r["q2"], mae=r["mae"], rmsd=r["rmsd"],
                                    threshold=threshold,
                                    csv_path=csv_path,
                                    model=r["model"],
                                    predictions=r["predictions"],
                                )
                            except Exception:
                                pass

                new_df = _to_df(flat_results)
                results_df = pd.concat([results_df, new_df], ignore_index=True) if not new_df.empty else results_df
          
            return results_df

        # ---- initial pass -------------------------------------------------------
        results = _evaluate_block(threshold)
        results = _sort_results(results)
      

        # If everything's -inf on Q2, relax threshold and retry (bounded loop)
        attempts = 0
        threshold = float(threshold)
        while _is_all_q2_neg_inf(results) and threshold > min_threshold:
            print("All Q2 values are -inf, lowering R2 threshold and retrying...")
            best = _best_r2(results)
            # choose the next threshold: either step down from best or by a fixed step
            if best is not None:
                threshold = max(min_threshold, best - 0.15)
            else:
                threshold = max(min_threshold, threshold - threshold_step)
            print(f"New threshold: {threshold:.3f}")
            results = _evaluate_block(threshold)
            results = _sort_results(results)
            attempts += 1
            if attempts > 5:
                # avoid infinite loops
                break

        # keep at least min_models_to_keep even if top_n is small
        results = results.head(max(top_n, min_models_to_keep))

        # Flag high multicollinearity (VIF > vif_threshold) for models with a calculated Q2
        results = self._flag_high_vif(results, vif_threshold=vif_threshold)

        # Show table when the optional plot module loaded
        printer = globals().get("print_models_regression_table")
        if not results.empty and callable(printer):
            printer(results, app, self)

        self.results=results

        return results



# ---------------------------------------------------------------------------
# Deferred cross-module imports — placed here so that all functions in this
# module are fully defined before plot.py tries to import them, avoiding the
# circular-import failure that occurs when these lines run at the top.
# Plotting extras (adjustText, shap, plotly) must not block LinearRegressionModel.
# ---------------------------------------------------------------------------
try:
    from modeling_utils import (
        simi_sampler, stratified_sampling_with_plots,
        _normalize_combination_to_columns, _parse_tuple_string,
    )
    from modeling_utils import *
except ImportError:
    from M3_modeler.modeling_utils import (
        simi_sampler, stratified_sampling_with_plots,
        _normalize_combination_to_columns, _parse_tuple_string,
    )
    from M3_modeler.modeling_utils import *

try:
    from plot import *
except ImportError:
    try:
        from M3_modeler.plot import *
    except ImportError:
        pass
