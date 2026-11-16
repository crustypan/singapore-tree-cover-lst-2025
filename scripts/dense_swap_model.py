"""Check and fit the local five-point pair-difference design."""
import numpy as np

from analyze_local_comparison import DRAWS, SEED
from prepare_local_comparison import BLOCK


COEFFICIENTS = ["tree_gap_pp", "built_gap_pp", "other_gap_pp", "land_gap_pp"]
CHANGE_PP = 5
MIN_RESAMPLE_RETENTION = 0.95


def design(pairs):
    x = pairs[COEFFICIENTS].to_numpy() / CHANGE_PP
    # An identically zero difference needs no adjustment.
    varies = np.any(x != 0, axis=0)
    _, dates = np.unique(pairs.date, return_inverse=True)
    w = 1 / np.bincount(dates)[dates]
    blocks, ids = np.unique(pairs[BLOCK], axis=0, return_inverse=True)
    return x[:, varies], w, blocks, ids, np.asarray(COEFFICIENTS)[varies]


def spatial_weights(w, ids, blocks, draws=DRAWS):
    rng = np.random.default_rng(SEED)
    multiplicity = rng.multinomial(blocks, np.full(blocks, 1 / blocks), size=draws)
    return w * multiplicity[:, ids]


def full_rank(x, weights):
    weighted = x * np.sqrt(weights[..., None])
    return np.linalg.matrix_rank(weighted) == x.shape[1]


def assess(pairs):
    distinct = len(pairs[["lower_tree_cell", "higher_tree_cell"]].drop_duplicates())
    summary = {
        "pairs": len(pairs), "spatial_pairs": distinct, "dates": int(pairs.date.nunique()),
        "blocks": len(pairs[BLOCK].drop_duplicates()), "matched_ready": False, "adjusted_ready": False,
    }
    if summary["blocks"] < 3 or summary["dates"] < 2:
        return {**summary, "reason": "fewer than three spatial blocks or two overpasses"}
    x, w, blocks, ids, terms = design(pairs)
    block_weights = np.bincount(ids, weights=w) / w.sum()
    summary.update(
        matched_ready=True, max_block_weight=float(block_weights.max()),
        effective_blocks=float(1 / (block_weights @ block_weights)),
        terms=terms.tolist(), draws_requested=DRAWS,
    )
    if "tree_gap_pp" not in terms or not full_rank(x, w):
        return {**summary, "reason": "rank-deficient adjustment design"}
    leave_block = w * (ids != np.arange(len(blocks))[:, None])
    retained = int(full_rank(x, spatial_weights(w, ids, len(blocks))).sum())
    summary.update(adjustment_draws_retained=retained,
                   rank_after_each_block_omission=bool(full_rank(x, leave_block).all()))
    if not summary["rank_after_each_block_omission"]:
        return {**summary, "reason": "adjustment loses rank after omitting a spatial block"}
    if retained < MIN_RESAMPLE_RETENTION * DRAWS:
        return {**summary, "reason": "adjustment has fewer than 95% full-rank draws"}
    return {**summary, "adjusted_ready": True, "reason": "ready"}


def adjusted_fits(x, y, weights):
    root = np.sqrt(weights)
    return np.array([
        np.linalg.lstsq(x * row[:, None], y * row, rcond=None)[0][0]
        for row in root
    ])
