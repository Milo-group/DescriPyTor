"""How many torsion wells an ensemble supports.

Angles are embedded as sin and cos, so -179 degrees and +179 degrees are
neighbours. Four votes are reported and none is hidden when they disagree:

- the knee of k-means inertia (Kneedle)
- the largest jump in the silhouette, and the silhouette maximum
- density peaks on the torus, counted when a peak is at least a tenth of the tallest
- the largest k that still puts every file in every well (needs at least two files)
"""
from __future__ import annotations

import numpy as np

K_MAX = 8


def _embed(angles):
    rad = np.radians(np.asarray(angles, float))
    if rad.ndim == 1:
        rad = rad[:, None]
    return np.concatenate([np.sin(rad), np.cos(rad)], axis=1)


def _kmeans(x, k, restarts=6, iters=40):
    rng = np.random.default_rng(0)
    n = len(x)
    best = None
    for _restart in range(restarts):
        centres = x[rng.choice(n, size=k, replace=False)].copy()
        labels = np.zeros(n, int)
        for _it in range(iters):
            dist = ((x[:, None, :] - centres[None, :, :]) ** 2).sum(axis=2)
            labels = dist.argmin(axis=1)
            for j in range(k):
                members = x[labels == j]
                if len(members) == 0:
                    centres[j] = x[int(rng.integers(0, n))]
                else:
                    centres[j] = members.mean(axis=0)
        dist = ((x[:, None, :] - centres[None, :, :]) ** 2).sum(axis=2)
        labels = dist.argmin(axis=1)
        inertia = float(dist[np.arange(n), labels].sum())
        if best is None or inertia < best[0]:
            best = (inertia, labels.copy())
    return best


def _silhouette(x, labels):
    n = len(x)
    use = np.arange(n)
    if n > 400:
        use = np.random.default_rng(1).choice(n, size=400, replace=False)
    sample = x[use]
    lab = labels[use]
    d = np.sqrt(((sample[:, None, :] - sample[None, :, :]) ** 2).sum(axis=2))
    scores = []
    for i in range(len(sample)):
        own = lab == lab[i]
        own[i] = False
        a = float(d[i, own].mean()) if own.any() else 0.0
        b = None
        for cluster in set(lab.tolist()):
            if cluster == lab[i]:
                continue
            mask = lab == cluster
            if not mask.any():
                continue
            mean = float(d[i, mask].mean())
            b = mean if b is None else min(b, mean)
        scores.append(0.0 if b is None or max(a, b) == 0 else (b - a) / max(a, b))
    return float(np.mean(scores)) if scores else 0.0


def _kneedle(ks, values):
    x = np.asarray(ks, float)
    y = np.asarray(values, float)
    if len(x) < 3 or x[-1] == x[0]:
        return int(x[0])
    xn = (x - x[0]) / (x[-1] - x[0])
    span = float(y.max() - y.min()) or 1.0
    yn = (y - y.min()) / span
    chord = yn[0] + (yn[-1] - yn[0]) * xn
    return int(x[int(np.argmax(chord - yn))])


def _density_count(angles, kappa=8.0, floor=0.10, neighbourhood=25.0):
    rad = np.radians(np.asarray(angles, float))
    if rad.ndim == 1:
        rad = rad[:, None]
    n = len(rad)
    dens = np.zeros(n)
    for i in range(n):
        delta = rad[i] - rad
        dens[i] = np.exp(kappa * np.cos(delta).sum(axis=1)).sum()
    peaks = []
    limit = np.radians(neighbourhood)
    for i in np.argsort(-dens):
        if dens[i] < floor * dens.max():
            break
        delta = np.abs((rad[i] - rad[peaks] + np.pi) % (2 * np.pi) - np.pi) if peaks else None
        if delta is not None and np.any(np.all(delta < limit, axis=1)):
            continue
        peaks.append(i)
    return len(peaks)


def _covers(labels, groups, k):
    present = {}
    for lab, group in zip(labels, groups):
        present.setdefault(str(group), set()).add(int(lab))
    if len(present) < 2:
        return None
    return all(len(seen) >= k for seen in present.values())


def vote_wells(angles, groups, k_max=K_MAX):
    angles = np.asarray(angles, float)
    if angles.ndim == 1:
        angles = angles[:, None]
    n, _p = angles.shape
    if n < 4:
        return {"error": "at least 4 conformers"}
    x = _embed(angles)
    hi = min(int(k_max), n - 1)
    rows, labels_by_k = [], {}
    sil_by_k = {}
    for k in range(2, hi + 1):
        inertia, labels = _kmeans(x, k)
        sil = _silhouette(x, labels)
        cover = _covers(labels, groups, k)
        rows.append({
            "k": k,
            "inertia": round(inertia, 4),
            "silhouette": round(sil, 4),
            "every_file": cover,
        })
        labels_by_k[str(k)] = [int(v) for v in labels]
        sil_by_k[k] = sil
    ks = [r["k"] for r in rows]
    inertias = [r["inertia"] for r in rows]
    knee = _kneedle(ks, inertias)
    jump_k, jump = ks[0], -1.0
    for prev, nxt in zip(ks, ks[1:]):
        gain = sil_by_k[nxt] - sil_by_k[prev]
        if gain > jump:
            jump, jump_k = gain, nxt
    sil_max = max(ks, key=lambda k: sil_by_k[k])
    density = _density_count(angles)
    coverage = None
    for row in rows:
        if row["every_file"]:
            coverage = row["k"]
        elif coverage is not None:
            break
    votes = [knee, jump_k, density]
    counts = {v: votes.count(v) for v in votes}
    agreed = [v for v, nvote in counts.items() if nvote >= 2 and 2 <= v <= hi]
    suggested = agreed[0] if agreed else (density if 2 <= density <= hi else knee)
    return {
        "criteria": {
            "knee": knee,
            "silhouette_jump": jump_k,
            "silhouette_max": sil_max,
            "density": density,
            "coverage": coverage,
        },
        "suggested": int(suggested),
        "rows": rows,
        "labels": labels_by_k,
    }
