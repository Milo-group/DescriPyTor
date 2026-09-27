"""
Nested cross-validation utilities for comparative descriptor-set modeling.

Mirrors the repeated stratified nested-CV protocol described in Lau, Borden,
Steiman, Wang, Parasram & Doyle, "Ni/Photoredox-Catalyzed Enantioselective
Cross-Electrophile Coupling of Styrene Oxides with Aryl Iodides", J. Am.
Chem. Soc. 2021, 143, 15873-15881 (SI Table S17): an outer loop estimates
generalization performance on data that never participated in feature
selection, while an inner loop performs exhaustive-subset feature selection.

Why this exists (as opposed to just calling `LinearRegressionModel.search_models`):
`search_models` does a single-pass exhaustive search over the *entire* dataset
and reports each combination's LOOCV/k-fold Q2. That protects against
overfitting the coefficients, but not against "search bias" -- if you try
thousands of combinations and report the CV score of whichever one happened
to win, that reported score is optimistic, because the same data used to
*pick* the winner was also used to *score* it. Nested CV fixes this by moving
feature selection inside the training fold of an outer, held-out evaluation
loop, so the reported test performance is for data the selection process
never saw.

Design, matching the reference paper unless noted:
  * Outer loop:  k_out-fold CV, repeated r_out times (default 4x5 -> 20 outer
    splits). Each outer test fold is scored by a model chosen without seeing it.
  * Inner loop:  within each outer training set, k_in-fold CV repeated r_in
    times (default 5x10 -> 50 inner folds). For each inner fold: rank all
    candidate feature combinations by training-set RSS, take the
    `top_k_by_rss` lowest-RSS combos (paper uses 5), evaluate those on the
    inner validation fold, and keep the one with lowest validation RMSE.
  * The most frequently selected combination across the 50 inner folds is the
    "champion" for that outer split (the paper's own wording).
  * Across all outer splits, tally how often each combination won -> the
    `count` column in the summary table, directly comparable to Table S17's
    Model / Average Training RMSE / Average Test RMSE / Adj. R2 / Count.

Stratification: pass `groups` (array-like, aligned to model row order) to use
StratifiedKFold in both loops -- e.g. ligand class, matching the paper.
Without it, plain shuffled KFold is used.

Scalability note: the reference paper pre-filtered its descriptor pool down to
20 candidates before nested CV (exhaustive C(20,3) = 1140 combinations, cheap
to loop over in R). DescriPyTor feature sets are often much larger
(C(68,3) ~= 50k), so a naive per-combination Python loop repeated for every
one of the ~1000 inner folds would mean tens of millions of individual
regressions. `rank_combos_by_train_rss` avoids that by batching every
combination of a given width into one vectorized linear-algebra call (no
per-combination Python loop, chunked to bound memory), and
`shortlist_candidate_features` gives an automatic, defensible way to prune a
large pool down to a Doyle-like ~20 candidates first.

Verified against the winning_features_lowest.csv case study (29 ligands, 68
features, leave_out=['L10','L27','L21','L18']): a reduced-scale run
(outer 4x1, inner 5x2) recovers the same ('nbo_diff_1-5', '0_bond_length_1-2',
'dihedral_[5, 1, 8, 6]') combination discussed throughout this notebook as one
of the outer-split champions.
"""

from __future__ import annotations

import itertools
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold, StratifiedKFold

try:
    from modeling import _linear_matrix_fit_predict
except ImportError:  # pragma: no cover - depends on caller's sys.path
    from M3_modeler.modeling import _linear_matrix_fit_predict

try:
    from modeling_utils import _normalize_combination_to_columns
except ImportError:  # pragma: no cover
    from M3_modeler.modeling_utils import _normalize_combination_to_columns


# ---------------------------------------------------------------------------
# small numeric helpers
# ---------------------------------------------------------------------------

def _rmse(y_true, y_pred) -> float:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    return float(np.sqrt(np.mean((y_true - y_pred) ** 2)))


def _r2(y_true, y_pred) -> float:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - y_true.mean()) ** 2)
    return float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")


def _adj_r2(r2: float, n: int, p: int) -> float:
    denom = n - p - 1
    if denom <= 0:
        return float("nan")
    return 1.0 - (1.0 - r2) * (n - 1) / denom


def _make_splitter(k: int, groups: Optional[np.ndarray], random_state: int):
    if groups is not None:
        return StratifiedKFold(n_splits=k, shuffle=True, random_state=random_state)
    return KFold(n_splits=k, shuffle=True, random_state=random_state)


def _split(splitter, X: np.ndarray, groups: Optional[np.ndarray]):
    if groups is not None:
        return splitter.split(X, groups)
    return splitter.split(X)


# ---------------------------------------------------------------------------
# vectorized exhaustive-subset RSS ranking (the expensive step, batched)
# ---------------------------------------------------------------------------

def _batched_train_rss(X_train: np.ndarray, y_train: np.ndarray,
                        combo_idx: np.ndarray, alpha: float = 1e-5) -> np.ndarray:
    """
    Training-set RSS for every combination in `combo_idx`, in one batched
    linear-algebra call (no per-combination Python loop).

    combo_idx : (n_combos, width) int array of column indices into X_train.
                All combos here must share the same `width`.
    Returns   : (n_combos,) array of RSS values.
    """
    m = X_train.shape[0]
    Xb = X_train[:, combo_idx]                # (m, n_combos, width) fancy index
    Xb = np.transpose(Xb, (1, 0, 2))          # (n_combos, m, width)

    mean = Xb.mean(axis=1, keepdims=True)
    std = Xb.std(axis=1, keepdims=True)
    std = np.where(std == 0, 1.0, std)
    Xb = (Xb - mean) / std

    ones = np.ones((Xb.shape[0], m, 1))
    Xb = np.concatenate([ones, Xb], axis=2)   # (n_combos, m, width+1)

    y = np.broadcast_to(y_train, (Xb.shape[0], m))  # (n_combos, m)

    gram = np.einsum("nip,niq->npq", Xb, Xb)
    pen = np.eye(gram.shape[-1]) * alpha
    pen[0, 0] = 0.0
    rhs = np.einsum("nip,ni->np", Xb, y)

    # np.linalg.solve's gufunc signature is (m,m),(m,n)->(m,n); with rhs as a
    # plain 2D (n_combos, width+1) array, numpy treats the *whole* array as a
    # single (m,n) rhs instead of batching over n_combos. Giving rhs an
    # explicit trailing axis makes each combo's target a (width+1, 1) matrix,
    # so solve correctly batches over the leading n_combos dimension.
    theta = np.linalg.solve(gram + pen, rhs[..., None])[..., 0]  # (n_combos, width+1)
    fitted = np.einsum("nip,np->ni", Xb, theta)  # (n_combos, m)
    rss = np.sum((y - fitted) ** 2, axis=1)
    return rss


def rank_combos_by_train_rss(X_train: np.ndarray, y_train: np.ndarray,
                              combos: Sequence[Tuple[str, ...]],
                              feature_index: Dict[str, int],
                              alpha: float = 1e-5,
                              chunk_size: int = 20_000) -> List[Tuple[Tuple[str, ...], float]]:
    """Rank every combo in `combos` by training RSS, batched per combo width
    and chunked so memory stays bounded regardless of how many combos there are."""
    by_width: Dict[int, List[Tuple[str, ...]]] = defaultdict(list)
    for c in combos:
        by_width[len(c)].append(c)

    out: List[Tuple[Tuple[str, ...], float]] = []
    for width, group in by_width.items():
        for start in range(0, len(group), chunk_size):
            chunk = group[start:start + chunk_size]
            idx = np.array([[feature_index[f] for f in c] for c in chunk], dtype=int)
            rss = _batched_train_rss(X_train, y_train, idx, alpha=alpha)
            out.extend(zip(chunk, rss.tolist()))

    out.sort(key=lambda t: t[1])
    return out


def select_best_combo_on_validation(X: np.ndarray, y: np.ndarray,
                                     train_idx: np.ndarray, val_idx: np.ndarray,
                                     shortlist: Sequence[Tuple[str, ...]],
                                     feature_index: Dict[str, int],
                                     alpha: float = 1e-5) -> Tuple[Tuple[str, ...], float]:
    """Among `shortlist`, keep the combo with lowest RMSE on (X[val_idx], y[val_idx])."""
    best_combo, best_rmse = None, np.inf
    for combo in shortlist:
        idx = [feature_index[f] for f in combo]
        y_hat = _linear_matrix_fit_predict(
            X[np.ix_(train_idx, idx)], y[train_idx], X[np.ix_(val_idx, idx)],
            alpha=alpha, scale=True,
        )
        rmse = _rmse(y[val_idx], y_hat)
        if rmse < best_rmse:
            best_combo, best_rmse = combo, rmse
    return best_combo, best_rmse


# ---------------------------------------------------------------------------
# main nested-CV driver
# ---------------------------------------------------------------------------

def nested_cv_search(
    model,
    candidate_features: Optional[Sequence[str]] = None,
    candidate_combos: Optional[Sequence[Tuple[str, ...]]] = None,
    min_features: int = 3,
    max_features: int = 3,
    outer_k: int = 4,
    outer_repeats: int = 5,
    inner_k: int = 5,
    inner_repeats: int = 10,
    top_k_by_rss: int = 5,
    groups: Optional[Sequence] = None,
    random_state: int = 42,
    alpha: float = 1e-5,
    verbose: bool = True,
    return_outer_predictions: bool = False,
):
    """
    Repeated stratified nested cross-validation, replicating the protocol in
    Lau et al. JACS 2021 SI Table S17.

    Parameters
    ----------
    model : modeling.LinearRegressionModel
        Must already have `.original_features_df`, `.target_vector`,
        `.molecule_names` populated (i.e. after `process_features_csv`, which
        `LinearRegressionModel.__init__` already runs).
    candidate_features : column names to search over. Defaults to
        `model.features_list`. For large feature sets (DescriPyTor commonly
        extracts 60-70 descriptors), pre-filter first with
        `shortlist_candidate_features` -- see module docstring. Ignored if
        `candidate_combos` is given.
    candidate_combos : an explicit list of feature-name tuples to use as the
        *entire* search space inside the inner loop, instead of generating
        every combination of `candidate_features` at widths
        `min_features..max_features`. Use this to run nested CV over a
        specific shortlist (e.g. the top combinations from a prior
        `search_models()` / brute-force run) rather than a fresh exhaustive
        search -- this is what makes DescriPyTor-scale feature pools tractable
        at Doyle's full 4x5/5x10 CV design. Any combo referencing a column not
        present in `model.original_features_df` is dropped with a warning
        (e.g. a feature name from a different feature-CSV/run than this
        model's).
    min_features, max_features : combination width to search when
        `candidate_combos` is not given. The reference paper fixed this at
        exactly 3 (default here); widen if you want to compare against your
        existing 3-4 feature convention.
    outer_k, outer_repeats : outer CV design (default 4-fold x 5 repeats = 20
        outer splits, matching the paper's k_out=4, 5 repeats, n=29).
    inner_k, inner_repeats : inner CV design used for feature selection only
        (default 5-fold x 10 repeats = 50 inner folds per outer split).
    top_k_by_rss : shortlist size before validation scoring (paper uses 5).
    groups : optional array-like aligned to model.molecule_names giving a
        stratification label (e.g. ligand class / coordination class). If
        omitted, plain shuffled KFold is used in both loops.
    return_outer_predictions : if True, also return a second DataFrame with
        one row per outer-test-fold prediction (sample_id, outer_repeat,
        outer_fold, combination, y_true, y_pred) -- everything
        `compute_nested_q2` needs to turn the champion table's
        `avg_test_rmse` into an exact nested Q2 (pooled PRESS / SS_tot over
        every held-out prediction), instead of only an RMSE average. Off by
        default to keep the return type unchanged for existing callers
        (e.g. `compare_representations`).

    Returns
    -------
    pd.DataFrame with columns
        combination, avg_train_rmse, avg_test_rmse, adj_r2_full_fit, count,
        of_total_outer_splits
    sorted by `count` descending -- directly comparable to Table S17.

    If `return_outer_predictions=True`, returns instead the tuple
    `(summary_df, outer_predictions_df)`, where `outer_predictions_df` has
    columns `[sample_id, outer_repeat, outer_fold, combination, y_true,
    y_pred]` -- one row per molecule per outer test fold it appeared in
    (each molecule appears once per outer repeat). Pass this straight to
    `compute_nested_q2` or `champion_report`.

    Runtime note: with the full ~60-70 feature DescriPyTor pool (min_features=
    max_features=3 -> C(68,3) ~= 50k combos), the default protocol
    (4x5 outer, 5x10 inner = 1000 inner fits) takes several minutes -- run it
    in a background thread (see the `search_models` cells earlier in this
    notebook for the pattern). Passing a pre-filtered `candidate_combos` (e.g.
    ~10-30 shortlisted combos) cuts the per-fold RSS-ranking cost
    proportionally and can make the full protocol finish in seconds.
    """
    # `original_features_df` retains every row from the source CSV, including
    # any `leave_out` molecules; `target_vector`/`molecule_names` are already
    # filtered to the training set. Re-align explicitly so X and y correspond
    # 1:1 -- silently using the unfiltered frame would leak held-out
    # molecules' features into the nested CV.
    raw_df = model.original_features_df.loc[model.molecule_names]
    feats = list(candidate_features) if candidate_features is not None else list(model.features_list)
    feature_index = {f: i for i, f in enumerate(raw_df.columns)}
    X_full = raw_df.to_numpy(dtype=float)
    y_full = np.asarray(model.target_vector, dtype=float).ravel()
    n = X_full.shape[0]

    if groups is not None:
        groups = np.asarray(groups)
        if groups.shape[0] != n:
            raise ValueError(f"groups has {groups.shape[0]} entries, expected {n}")

    if candidate_combos is not None:
        combos = []
        dropped = []
        for c in candidate_combos:
            c = tuple(c)
            if all(f in feature_index for f in c):
                combos.append(c)
            else:
                dropped.append(c)
        if not combos:
            raise ValueError(
                "None of the provided candidate_combos reference columns present "
                "in this model's feature set -- check they came from the same CSV."
            )
        if verbose:
            print(f"Using {len(combos)} explicit candidate combination(s)"
                  + (f" ({len(dropped)} dropped -- feature(s) not in this dataset: {dropped})"
                     if dropped else ""))
    else:
        combos = list(itertools.chain.from_iterable(
            itertools.combinations(feats, w) for w in range(min_features, max_features + 1)
        ))
        if verbose:
            print(f"Candidate features: {len(feats)}  ->  {len(combos):,} combinations "
                  f"(width {min_features}-{max_features})")

    top_k_by_rss = min(top_k_by_rss, len(combos))

    total_outer_splits = outer_k * outer_repeats
    if verbose:
        print(f"Outer: {outer_k}-fold x {outer_repeats} repeats = {total_outer_splits} splits")
        print(f"Inner: {inner_k}-fold x {inner_repeats} repeats = "
              f"{inner_k * inner_repeats} folds per outer split")

    champion_counter: Counter = Counter()
    train_rmse: Dict[tuple, List[float]] = defaultdict(list)
    test_rmse: Dict[tuple, List[float]] = defaultdict(list)
    outer_prediction_rows: List[dict] = []
    molecule_names_arr = np.asarray(model.molecule_names)

    outer_split_id = 0
    for rep in range(outer_repeats):
        outer_splitter = _make_splitter(outer_k, groups, random_state + rep)
        fold_in_rep = 0
        for outer_train_idx, outer_test_idx in _split(outer_splitter, X_full, groups):
            outer_split_id += 1
            fold_in_rep += 1
            inner_groups = groups[outer_train_idx] if groups is not None else None
            inner_counter: Counter = Counter()

            for inner_rep in range(inner_repeats):
                seed = random_state + 1000 * rep + 17 * inner_rep + outer_split_id
                inner_splitter = _make_splitter(inner_k, inner_groups, seed)
                sub_X = X_full[outer_train_idx]
                for tr_local, val_local in _split(inner_splitter, sub_X, inner_groups):
                    tr_idx = outer_train_idx[tr_local]
                    val_idx = outer_train_idx[val_local]
                    ranked = rank_combos_by_train_rss(
                        X_full[tr_idx], y_full[tr_idx], combos, feature_index, alpha=alpha
                    )
                    shortlist = [c for c, _ in ranked[:top_k_by_rss]]
                    best_combo, _ = select_best_combo_on_validation(
                        X_full, y_full, tr_idx, val_idx, shortlist, feature_index, alpha=alpha
                    )
                    inner_counter[best_combo] += 1

            champion_combo = inner_counter.most_common(1)[0][0]
            champion_counter[champion_combo] += 1

            idx = [feature_index[f] for f in champion_combo]
            y_hat_tr = _linear_matrix_fit_predict(
                X_full[np.ix_(outer_train_idx, idx)], y_full[outer_train_idx],
                X_full[np.ix_(outer_train_idx, idx)], alpha=alpha, scale=True,
            )
            y_hat_te = _linear_matrix_fit_predict(
                X_full[np.ix_(outer_train_idx, idx)], y_full[outer_train_idx],
                X_full[np.ix_(outer_test_idx, idx)], alpha=alpha, scale=True,
            )
            train_rmse[champion_combo].append(_rmse(y_full[outer_train_idx], y_hat_tr))
            test_rmse[champion_combo].append(_rmse(y_full[outer_test_idx], y_hat_te))

            if return_outer_predictions:
                for sample_id, y_true_i, y_pred_i in zip(
                    molecule_names_arr[outer_test_idx], y_full[outer_test_idx], y_hat_te
                ):
                    outer_prediction_rows.append({
                        "sample_id": sample_id,
                        "outer_repeat": rep,
                        "outer_fold": fold_in_rep,
                        "combination": champion_combo,
                        "y_true": float(y_true_i),
                        "y_pred": float(y_pred_i),
                    })

            if verbose:
                print(f"  outer split {outer_split_id}/{total_outer_splits}: "
                      f"champion = {champion_combo}")

    rows = []
    for combo, count in champion_counter.most_common():
        idx = [feature_index[f] for f in combo]
        y_hat_full = _linear_matrix_fit_predict(
            X_full[:, idx], y_full, X_full[:, idx], alpha=alpha, scale=True
        )
        r2_full = _r2(y_full, y_hat_full)
        rows.append({
            "combination": combo,
            "avg_train_rmse": float(np.mean(train_rmse[combo])),
            "avg_test_rmse": float(np.mean(test_rmse[combo])),
            "adj_r2_full_fit": _adj_r2(r2_full, n, len(combo)),
            "count": count,
            "of_total_outer_splits": total_outer_splits,
        })

    summary_df = pd.DataFrame(rows).sort_values("count", ascending=False).reset_index(drop=True)

    if return_outer_predictions:
        outer_predictions_df = pd.DataFrame(outer_prediction_rows)
        return summary_df, outer_predictions_df

    return summary_df


# ---------------------------------------------------------------------------
# exact nested Q2 from stored outer-test predictions
# ---------------------------------------------------------------------------

def compute_nested_q2(outer_predictions: pd.DataFrame,
                       combination: Optional[Tuple[str, ...]] = None) -> Dict[str, float]:
    """
    Exact nested Q2 = 1 - PRESS/SS_tot, computed directly from stored
    outer-test predictions (the `outer_predictions` DataFrame returned by
    `nested_cv_search(..., return_outer_predictions=True)`), instead of
    approximating it from an RMSE average.

    Parameters
    ----------
    outer_predictions : DataFrame with columns
        [sample_id, outer_repeat, outer_fold, combination, y_true, y_pred].
    combination : if given, restrict to the outer-test predictions made when
        this exact combination was the outer split's champion (i.e. the same
        subset `avg_test_rmse`/`count` in the summary table are computed
        over). If None, uses every outer-test prediction regardless of which
        combination won that split -- the single overall "procedure" nested
        Q2 for the whole nested-CV run.

    Returns
    -------
    dict with:
        q2_pooled          : PRESS/SS_tot computed once over every matching
                              prediction, pooled across repeats.
        q2_per_repeat_mean  : mean of per-repeat Q2 (each repeat's outer folds
                              cover every sample exactly once, so this is the
                              more paper-faithful "repeated CV" Q2).
        q2_per_repeat_std   : std-dev across repeats (0.0 if only one repeat).
        n_predictions       : how many outer-test predictions went into
                              q2_pooled (sanity-check against `count *
                              outer_k` from the summary table).
    """
    df = outer_predictions
    if combination is not None:
        df = df[df["combination"] == tuple(combination)]

    if df.empty:
        return {
            "q2_pooled": float("nan"),
            "q2_per_repeat_mean": float("nan"),
            "q2_per_repeat_std": float("nan"),
            "n_predictions": 0,
        }

    def _q2(sub: pd.DataFrame) -> float:
        y = sub["y_true"].to_numpy(dtype=float)
        yhat = sub["y_pred"].to_numpy(dtype=float)
        press = float(np.sum((y - yhat) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        return 1.0 - press / ss_tot if ss_tot > 0 else float("nan")

    q2_pooled = _q2(df)
    per_repeat = df.groupby("outer_repeat").apply(_q2)

    return {
        "q2_pooled": float(q2_pooled),
        "q2_per_repeat_mean": float(per_repeat.mean()),
        "q2_per_repeat_std": float(per_repeat.std(ddof=1)) if len(per_repeat) > 1 else 0.0,
        "n_predictions": int(len(df)),
    }


# ---------------------------------------------------------------------------
# shortlist helper for large descriptor pools
# ---------------------------------------------------------------------------

def shortlist_candidate_features(model, top_n_combos: int = 200, top_k_features: int = 20,
                                  n_jobs: int = -1, threshold: float = 0.5) -> List[str]:
    """
    Reduce a large feature pool (e.g. DescriPyTor's ~60-70 descriptors) down to
    a Doyle-like ~20 candidates before nested CV, by frequency of appearance
    among the best combinations of a broad single-pass search.

    This mirrors what the reference paper's authors did by hand (picking 20
    candidate descriptors up front); it is *not* a substitute for nested CV --
    it only controls the combinatorial cost of the exhaustive step inside it.
    Uses the model's own `search_models`, so results also land in the model's
    usual results DB as a side effect.
    """
    results = model.search_models(
        top_n=top_n_combos, n_jobs=n_jobs, threshold=threshold,
        min_models_to_keep=top_n_combos,
    )
    counts: Counter = Counter()
    for combo in results["combination"]:
        for f in _normalize_combination_to_columns(combo):
            counts[f] += 1
    return [f for f, _ in counts.most_common(top_k_features)]


# ---------------------------------------------------------------------------
# decorrelation: prune a candidate pool BEFORE any search sees the response
# ---------------------------------------------------------------------------

def decorrelate_features(
    model,
    candidate_features: Sequence[str],
    corr_threshold: float = 0.9,
) -> Tuple[List[str], pd.DataFrame]:
    """
    Cluster `candidate_features` by pairwise |correlation| > `corr_threshold`
    (single-linkage: any chain of pairwise-correlated features ends up in one
    cluster) and keep exactly one representative per cluster -- whichever
    member has the strongest |correlation| with the response.

    Why this matters more than any CV setting: RSS/RMSE-ranked search (brute
    force or the inner loop of `nested_cv_search`) will always favor a pair of
    near-duplicate features over a single clean one, because together they
    give near-free extra fitting power -- that's a property of least squares,
    not something nested CV corrects for. Nested CV only guards against
    overfitting the *choice of combination*; it does nothing about
    multicollinearity *within* the winning combination. The fix has to happen
    in the candidate pool, before search.

    This also happens to close a previously disclosed gap: because clustering
    here runs on feature-feature correlation only (not feature-response Q2),
    the resulting pool is response-blind, unlike a Q2-frequency shortlist.

    Parameters
    ----------
    model : modeling.LinearRegressionModel
        Must already have `.original_features_df`, `.target_vector`,
        `.molecule_names` populated.
    candidate_features : the pool to prune (e.g. a prior
        `shortlist_candidate_features` output, or the full feature list).
    corr_threshold : |r| above which two features are considered redundant
        and merged into the same cluster. 0.9 is a conservative default;
        the correlation screen run earlier in this notebook found several
        clusters well above this (0.995, 0.998, 0.96-0.98).

    Returns
    -------
    (pruned_features, cluster_report)
        pruned_features : one representative feature name per cluster --
            pass this straight into `search_models`/`nested_cv_search` in
            place of the original pool.
        cluster_report : DataFrame with one row per cluster --
            [cluster_representative, cluster_size, cluster_members,
            target_corr_of_rep] -- so you can see what got merged and why a
            given representative was chosen. Clusters of size 1 are features
            that weren't correlated with anything else above threshold.
    """
    candidate_features = list(candidate_features)
    raw_df = model.original_features_df.loc[model.molecule_names]
    X = raw_df[candidate_features].astype(float)
    y = np.asarray(model.target_vector, dtype=float).ravel()

    # pairwise |correlation| among candidates only
    corr = X.corr().abs()

    # union-find clustering on |r| > threshold
    parent = {f: f for f in candidate_features}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i, f1 in enumerate(candidate_features):
        for f2 in candidate_features[i + 1:]:
            if corr.loc[f1, f2] > corr_threshold:
                union(f1, f2)

    clusters: Dict[str, List[str]] = defaultdict(list)
    for f in candidate_features:
        clusters[find(f)].append(f)

    # univariate |corr| with the response, used only to pick a representative
    # within an already-formed cluster -- clustering itself never looks at y
    target_corr = {
        f: float(abs(np.corrcoef(X[f].to_numpy(), y)[0, 1])) for f in candidate_features
    }

    pruned: List[str] = []
    rows = []
    for members in clusters.values():
        best = max(members, key=lambda f: target_corr[f])
        pruned.append(best)
        rows.append({
            "cluster_representative": best,
            "cluster_size": len(members),
            "cluster_members": members,
            "target_corr_of_rep": target_corr[best],
        })

    report = (
        pd.DataFrame(rows)
        .sort_values(["cluster_size", "target_corr_of_rep"], ascending=[False, False])
        .reset_index(drop=True)
    )
    return pruned, report


def vif_prune(
    model,
    candidate_features: Sequence[str],
    vif_threshold: float = 5.0,
    min_features: int = 6,
    verbose: bool = True,
) -> Tuple[List[str], pd.DataFrame]:
    """
    Iteratively drop the single highest-VIF feature from `candidate_features`,
    recomputing VIF each time, until every remaining feature is at or below
    `vif_threshold` (or only `min_features` remain).

    Necessary *in addition to* `decorrelate_features`, not instead of it:
    pairwise correlation clustering only catches 2-way redundancy (feature A
    duplicates feature B). VIF_i = 1/(1-R_i^2), where R_i^2 comes from
    regressing feature i on *every other* feature in the pool -- it catches
    N-way linear dependency even when no single pairwise |r| is high. This is
    exactly what showed up after pairwise pruning here: a 14-feature pool
    with no remaining pair above |r|=0.9 still had VIFs in the hundreds of
    thousands, meaning several features are jointly, almost exactly
    collinear as a group.

    Parameters
    ----------
    vif_threshold : conventional cutoffs are 5 (conservative) or 10
        (permissive); 5 is the default here.
    min_features : stop early if the pool would shrink below this, so you
        don't strip a search down to nothing on a very collinear pool --
        inspect the returned history and consider raising `vif_threshold`
        or reworking the candidate list if this triggers.

    Returns
    -------
    (pool, history)
        pool : surviving feature names.
        history : DataFrame, one row per drop -- [dropped, vif_at_drop,
            pool_size_before] -- so you can see exactly what was removed and
            at what VIF, in order.
    """
    from statsmodels.stats.outliers_influence import variance_inflation_factor

    pool = list(candidate_features)
    raw_df = model.original_features_df.loc[model.molecule_names]
    history = []

    while len(pool) > min_features:
        X = raw_df[pool].astype(float).to_numpy()
        vifs = [variance_inflation_factor(X, i) for i in range(X.shape[1])]
        worst_i = int(np.argmax(vifs))
        worst_feat, worst_vif = pool[worst_i], float(vifs[worst_i])
        if worst_vif <= vif_threshold:
            break
        history.append({
            "dropped": worst_feat,
            "vif_at_drop": worst_vif,
            "pool_size_before": len(pool),
        })
        if verbose:
            print(f"  dropping {worst_feat!r} (VIF={worst_vif:,.1f}), "
                  f"{len(pool) - 1} feature(s) remain")
        pool.pop(worst_i)

    return pool, pd.DataFrame(history)


def pool_vif(model, candidate_features: Sequence[str]) -> pd.DataFrame:
    """
    VIF of every feature in `candidate_features` against every other feature
    in the same pool (not just within one winning combination) -- a quick
    sanity check that `decorrelate_features` actually removed the
    problem pairs. Requires statsmodels.
    """
    from statsmodels.stats.outliers_influence import variance_inflation_factor

    candidate_features = list(candidate_features)
    raw_df = model.original_features_df.loc[model.molecule_names]
    X = raw_df[candidate_features].astype(float).to_numpy()
    vifs = [variance_inflation_factor(X, i) for i in range(X.shape[1])]
    return (
        pd.DataFrame({"feature": candidate_features, "vif": vifs})
        .sort_values("vif", ascending=False)
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# champion report: top-N candidates, exact nested Q2, VIF, pick-one
# ---------------------------------------------------------------------------

def champion_report(
    table: pd.DataFrame,
    outer_predictions: Optional[pd.DataFrame] = None,
    model: Optional["LinearRegressionModel"] = None,
    top_n: int = 3,
) -> List[Tuple[str, ...]]:
    """
    Print a champion-model report for the top `top_n` rows of a
    `nested_cv_search` summary table (by default just the single champion --
    `top_n=3` shows the top 3 competing combinations instead).

    For each of the top `top_n` combinations, prints: features, outer-fold
    win count/rate, average train/test RMSE, full-fit adjusted R2, exact
    nested Q2 (if `outer_predictions` is supplied -- see
    `nested_cv_search(..., return_outer_predictions=True)`), and a VIF table
    (if `model` is supplied). Also prints the single overall "procedure"
    nested Q2 -- pooled over every outer-test prediction regardless of which
    combination won that split -- as the one number that best represents the
    whole nested-CV run's generalization performance.

    Parameters
    ----------
    table : the DataFrame returned by `nested_cv_search` (sorted by `count`
        descending; row 0 is the single most-frequent champion).
    outer_predictions : optional, the second element of the tuple returned by
        `nested_cv_search(..., return_outer_predictions=True)`. Without it,
        the report falls back to the RMSE-average comparison only and prints
        a note that exact nested Q2 isn't available.
    model : optional `LinearRegressionModel`, needed only for the VIF section
        (uses `model.original_features_df` / `model.molecule_names`, same as
        `pool_vif`). Skipped if not given.
    top_n : how many rows of `table` to report (default 3).

    Returns
    -------
    List of the top `top_n` combinations (as feature-name tuples), in the
    same order printed -- i.e. index 0 is `table.iloc[0]['combination']`.
    Use this to let the user pick one afterward, e.g.::

        top_combos = ncv.champion_report(bg_nested['table'], bg_nested['outer_predictions'], MODEL, top_n=3)
        CHOICE = 0  # index into the printed list -- change and re-run to compare
        chosen_combo = top_combos[CHOICE]
    """
    top_n = min(top_n, len(table))
    top_rows = table.head(top_n).reset_index(drop=True)
    total_outer_splits = int(top_rows.iloc[0]["of_total_outer_splits"]) if top_n else 0

    print("=" * 80)
    print(f"TOP {top_n} CHAMPION MODEL{'S' if top_n != 1 else ''}")
    print("=" * 80)

    top_combos: List[Tuple[str, ...]] = []
    for i, row in top_rows.iterrows():
        combo = tuple(row["combination"])
        top_combos.append(combo)
        win_rate = 100.0 * row["count"] / row["of_total_outer_splits"]

        print(f"\n--- #{i}  {'(current champion)' if i == 0 else ''} ---".rstrip())
        print("Features:")
        for f in combo:
            print(f"   {f}")
        print(f"Outer-fold wins       : {int(row['count'])}/{int(row['of_total_outer_splits'])} "
              f"({win_rate:.1f}%)")
        print(f"Average training RMSE : {row['avg_train_rmse']:.4f}")
        print(f"Average test RMSE     : {row['avg_test_rmse']:.4f}")
        print(f"Full-fit adjusted R2  : {row['adj_r2_full_fit']:.4f}")

        if outer_predictions is not None and not outer_predictions.empty:
            q2 = compute_nested_q2(outer_predictions, combination=combo)
            n_repeats_seen = outer_predictions.loc[
                outer_predictions["combination"] == combo, "outer_repeat"
            ].nunique()
            print(f"Nested Q2 (this combo's own {q2['n_predictions']} outer-test "
                  f"prediction(s), pooled)      : {q2['q2_pooled']:.4f}")
            print(f"Nested Q2 (per-repeat mean +/- SD, {n_repeats_seen} repeat(s)) : "
                  f"{q2['q2_per_repeat_mean']:.4f} +/- {q2['q2_per_repeat_std']:.4f}")
        else:
            print("Nested Q2              : n/a -- rerun nested_cv_search(..., "
                  "return_outer_predictions=True) and pass the result here.")

        if model is not None:
            try:
                vif_table = pool_vif(model, list(combo))
                print("VIF (this combination's own features against each other):")
                print(vif_table.to_string(index=False))
            except Exception as e:
                print(f"VIF could not be computed: {e}")

    if outer_predictions is not None and not outer_predictions.empty:
        overall_q2 = compute_nested_q2(outer_predictions, combination=None)
        print("\n" + "=" * 80)
        print("OVERALL NESTED Q2 (whole procedure, every outer-test prediction pooled "
              "regardless of which combo won that split)")
        print("=" * 80)
        print(f"Q2 (pooled)                 : {overall_q2['q2_pooled']:.4f}  "
              f"(n={overall_q2['n_predictions']} predictions)")
        print(f"Q2 (mean +/- SD per repeat)  : {overall_q2['q2_per_repeat_mean']:.4f} "
              f"+/- {overall_q2['q2_per_repeat_std']:.4f}")
    else:
        print("\n" + "=" * 80)
        print("NESTED Q2 NOT AVAILABLE")
        print("=" * 80)
        print("Pass `outer_predictions` (rerun nested_cv_search with "
              "return_outer_predictions=True) to get an exact nested Q2 instead of "
              "just the average test RMSE.")

    if len(top_combos) > 1:
        counts = [int(r["count"]) for _, r in top_rows.iterrows()]
        if counts[0] < 2 * counts[-1]:
            print("\nNo strongly dominant combination -- interpret the top "
                  f"{top_n} models as competing alternatives.")

    return top_combos


# ---------------------------------------------------------------------------
# comparing multiple representations (Doyle Table S17 (a)-(d) + 5x2cv bottom line)
# ---------------------------------------------------------------------------

def compare_representations(
    models: Dict[str, "LinearRegressionModel"],
    candidate_features: Optional[Dict[str, Sequence[str]]] = None,
    groups: Optional[Dict[str, Sequence]] = None,
    n_5x2_repeats: int = 5,
    random_state: int = 42,
    **nested_cv_kwargs,
) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
    """
    Run `nested_cv_search` for each named representation (e.g.
    {'L*NiArCl': model_a, 'L*NiF2': model_b, 'Free ligand': model_c}), then
    compare their overall champion models with repeated 5x2 cross-validation
    (Dietterich, 1998) -- the same final comparison step the reference paper
    uses to rank molecular representations against each other and against a
    random-descriptor control.

    Returns
    -------
    (per_representation_tables, summary_table)
        per_representation_tables : dict of Table-S17-style DataFrames, one
            per representation (sorted by `count`, champion is row 0).
        summary_table : DataFrame with [representation, champion,
            test_rmse_mean, test_rmse_std, n_5x2_folds] -- directly comparable
            to the paper's final "Test RMSE, mean +/- SD" table.
    """
    per_rep_tables: Dict[str, pd.DataFrame] = {}
    champions: Dict[str, Tuple[str, ...]] = {}

    for name, model in models.items():
        feats = (candidate_features or {}).get(name)
        grp = (groups or {}).get(name)
        if nested_cv_kwargs.get("verbose", True):
            print(f"\n=== Nested CV: {name} ===")
        table = nested_cv_search(model, candidate_features=feats, groups=grp,
                                  random_state=random_state, **nested_cv_kwargs)
        per_rep_tables[name] = table
        champions[name] = tuple(table.iloc[0]["combination"])

    summary_rows = []
    for name, model in models.items():
        raw_df = model.original_features_df.loc[model.molecule_names]
        feature_index = {f: i for i, f in enumerate(raw_df.columns)}
        X_full = raw_df.to_numpy(dtype=float)
        y_full = np.asarray(model.target_vector, dtype=float).ravel()
        idx = [feature_index[f] for f in champions[name]]

        rmses = []
        for rep in range(n_5x2_repeats):
            kf = KFold(n_splits=2, shuffle=True, random_state=random_state + rep)
            for train_idx, test_idx in kf.split(X_full):
                y_hat = _linear_matrix_fit_predict(
                    X_full[np.ix_(train_idx, idx)], y_full[train_idx],
                    X_full[np.ix_(test_idx, idx)], alpha=1e-5, scale=True,
                )
                rmses.append(_rmse(y_full[test_idx], y_hat))

        summary_rows.append({
            "representation": name,
            "champion": champions[name],
            "test_rmse_mean": float(np.mean(rmses)),
            "test_rmse_std": float(np.std(rmses, ddof=1)),
            "n_5x2_folds": len(rmses),
        })

    summary_table = pd.DataFrame(summary_rows).sort_values("test_rmse_mean").reset_index(drop=True)
    return per_rep_tables, summary_table


def dietterich_5x2cv_test(model_a, features_a: Sequence[str],
                           model_b, features_b: Sequence[str],
                           n_repeats: int = 5, random_state: int = 42) -> Tuple[float, float]:
    """
    Dietterich's (1998) paired 5x2cv test comparing two champion models
    (which may come from different representations/datasets, as long as the
    row order/molecule set lines up -- pass the same `model` twice with two
    different `features_*` lists to compare two combinations on one dataset).

    Returns (t_statistic, p_value) for the null hypothesis that both models
    have the same generalization RMSE. This is the test the reference paper
    cites when comparing molecular representations.
    """
    from scipy import stats as _stats

    raw_a = model_a.original_features_df.loc[model_a.molecule_names]
    y_a = np.asarray(model_a.target_vector, dtype=float).ravel()
    idx_a = [raw_a.columns.get_loc(f) for f in features_a]
    X_a = raw_a.to_numpy(dtype=float)

    raw_b = model_b.original_features_df.loc[model_b.molecule_names]
    y_b = np.asarray(model_b.target_vector, dtype=float).ravel()
    idx_b = [raw_b.columns.get_loc(f) for f in features_b]
    X_b = raw_b.to_numpy(dtype=float)

    if X_a.shape[0] != X_b.shape[0]:
        raise ValueError(
            "model_a and model_b must have the same number of rows (same "
            "molecule set) for a paired 5x2cv test."
        )

    p1_first_rep = None
    variances = []
    for rep in range(n_repeats):
        kf = KFold(n_splits=2, shuffle=True, random_state=random_state + rep)
        splits = list(kf.split(X_a))
        p_this_rep = []
        for train_idx, test_idx in splits:
            yhat_a = _linear_matrix_fit_predict(
                X_a[np.ix_(train_idx, idx_a)], y_a[train_idx], X_a[np.ix_(test_idx, idx_a)],
                alpha=1e-5, scale=True,
            )
            yhat_b = _linear_matrix_fit_predict(
                X_b[np.ix_(train_idx, idx_b)], y_b[train_idx], X_b[np.ix_(test_idx, idx_b)],
                alpha=1e-5, scale=True,
            )
            rmse_a = _rmse(y_a[test_idx], yhat_a)
            rmse_b = _rmse(y_b[test_idx], yhat_b)
            p_this_rep.append(rmse_a - rmse_b)

        p1, p2 = p_this_rep
        pbar = (p1 + p2) / 2.0
        s2 = (p1 - pbar) ** 2 + (p2 - pbar) ** 2
        variances.append(s2)
        if rep == 0:
            p1_first_rep = p1

    t_stat = p1_first_rep / np.sqrt(np.mean(variances))
    df = n_repeats
    p_value = float(2 * (1 - _stats.t.cdf(abs(t_stat), df)))
    return float(t_stat), p_value
