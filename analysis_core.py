from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor


GroupMethod = Literal["rank", "quantile"]
TreeMethod = Literal["RF", "ET"]
TransformMethod = Literal["log2p1", "none"]


@dataclass(frozen=True)
class StratificationResult:
    target_gene: str
    expression: pd.Series
    groups: pd.Series
    high_samples: list[str]
    low_samples: list[str]
    excluded_samples: list[str]
    low_cutoff: float
    high_cutoff: float
    tie_warning: str | None


@dataclass(frozen=True)
class Genie3Result:
    vim: pd.DataFrame
    edges: pd.DataFrame
    expression_used: pd.DataFrame  # samples x genes


def _clean_axis(labels: Sequence[object]) -> list[str]:
    return [str(x).strip() for x in labels]


def clean_expression_matrix(
    matrix: pd.DataFrame,
    *,
    matrix_name: str,
    duplicate_gene_rule: Literal["mean", "sum"],
) -> pd.DataFrame:
    """Validate a genes x samples matrix and collapse duplicate gene IDs."""
    if matrix.empty:
        raise ValueError(f"{matrix_name} is empty.")

    out = matrix.copy()
    out.index = _clean_axis(out.index)
    out.columns = _clean_axis(out.columns)

    if any(x == "" for x in out.index):
        raise ValueError(f"{matrix_name} contains empty gene IDs.")
    if any(x == "" for x in out.columns):
        raise ValueError(f"{matrix_name} contains empty sample IDs.")
    if out.columns.duplicated().any():
        dup = out.columns[out.columns.duplicated()].unique().tolist()
        raise ValueError(f"{matrix_name} contains duplicated sample IDs: {dup[:8]}")

    out = out.apply(pd.to_numeric, errors="coerce")
    if out.isna().all(axis=None):
        raise ValueError(f"{matrix_name} contains no numeric values.")

    if out.index.duplicated().any():
        if duplicate_gene_rule == "sum":
            out = out.groupby(level=0, sort=False).sum(min_count=1)
        else:
            out = out.groupby(level=0, sort=False).mean()

    return out


def resolve_gene_name(gene: str, genes: Sequence[str]) -> str:
    """Resolve an exact or case-insensitive gene name without fuzzy matching."""
    query = str(gene).strip()
    if not query:
        raise ValueError("Target gene is empty.")
    genes_list = list(genes)
    if query in genes_list:
        return query

    matches = [g for g in genes_list if g.upper() == query.upper()]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(f"Gene name '{query}' is ambiguous after case-insensitive matching: {matches}")
    raise ValueError(f"Gene '{query}' was not found in the matrix.")


def align_samples(
    normalized: pd.DataFrame,
    raw_counts: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str], list[str]]:
    """Align matrices to shared samples while retaining all genes in each matrix."""
    shared = [s for s in normalized.columns if s in raw_counts.columns]
    if len(shared) < 8:
        raise ValueError(
            f"Only {len(shared)} shared samples were found. At least 8 are required, "
            "and 12+ are strongly recommended for quartile stratification."
        )
    norm_only = [s for s in normalized.columns if s not in raw_counts.columns]
    count_only = [s for s in raw_counts.columns if s not in normalized.columns]
    return normalized.loc[:, shared], raw_counts.loc[:, shared], shared, norm_only, count_only


def validate_raw_counts(raw_counts: pd.DataFrame) -> pd.DataFrame:
    """Require finite, non-negative integer raw counts."""
    if raw_counts.isna().any(axis=None):
        locations = int(raw_counts.isna().sum().sum())
        raise ValueError(f"Raw count matrix contains {locations} missing/non-numeric values.")
    values = raw_counts.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Raw count matrix contains infinite values.")
    if (values < 0).any():
        raise ValueError("Raw count matrix contains negative values.")
    if not np.allclose(values, np.rint(values), atol=1e-8):
        raise ValueError(
            "Raw count matrix contains non-integer values. DESeq2 requires unnormalized integer counts."
        )
    return pd.DataFrame(
        np.rint(values).astype(np.int64),
        index=raw_counts.index,
        columns=raw_counts.columns,
    )


def stratify_samples_by_target(
    normalized: pd.DataFrame,
    target_gene: str,
    fraction: float = 0.25,
    method: GroupMethod = "rank",
    min_group_size: int = 3,
) -> StratificationResult:
    """Create target-high and target-low sample groups from normalized expression."""
    if not 0 < fraction < 0.5:
        raise ValueError("fraction must be between 0 and 0.5.")

    target = resolve_gene_name(target_gene, normalized.index)
    expr = pd.to_numeric(normalized.loc[target], errors="coerce").dropna().astype(float)
    expr.index = _clean_axis(expr.index)

    if not np.isfinite(expr.to_numpy()).all():
        raise ValueError(f"{target} expression contains infinite values.")
    if expr.nunique(dropna=True) < 3:
        raise ValueError(
            f"{target} has fewer than three distinct expression values; robust high/low stratification is not possible."
        )

    low_q = float(expr.quantile(fraction))
    high_q = float(expr.quantile(1.0 - fraction))
    tie_warning: str | None = None

    if method == "quantile":
        if low_q >= high_q:
            raise ValueError(
                "Lower and upper quantile cutoffs overlap because the target expression contains too many ties. "
                "Use rank-based grouping or choose another target."
            )
        low_samples = expr.index[expr <= low_q].tolist()
        high_samples = expr.index[expr >= high_q].tolist()
    elif method == "rank":
        n_each = int(np.floor(len(expr) * fraction))
        if n_each < min_group_size:
            raise ValueError(
                f"Only {n_each} samples would enter each group. Increase cohort size or grouping fraction."
            )
        ordered = (
            expr.rename("expression")
            .to_frame()
            .assign(sample=lambda x: x.index.astype(str))
            .sort_values(["expression", "sample"], kind="mergesort")
        )
        low_samples = ordered.index[:n_each].tolist()
        high_samples = ordered.index[-n_each:].tolist()

        low_boundary = float(ordered.iloc[n_each - 1]["expression"])
        high_boundary = float(ordered.iloc[-n_each]["expression"])
        low_tied = int((expr == low_boundary).sum()) > 1
        high_tied = int((expr == high_boundary).sum()) > 1
        if low_tied or high_tied:
            tie_warning = (
                "The target expression has ties at a rank boundary. Exact group sizes were retained, "
                "so tied samples were deterministically ordered by sample ID."
            )
    else:
        raise ValueError(f"Unsupported grouping method: {method}")

    overlap = set(low_samples).intersection(high_samples)
    if overlap:
        raise ValueError(f"High and low groups overlap: {sorted(overlap)[:8]}")
    if min(len(low_samples), len(high_samples)) < min_group_size:
        raise ValueError(
            f"Group sizes are Low={len(low_samples)} and High={len(high_samples)}; "
            f"at least {min_group_size} per group are required."
        )

    groups = pd.Series(index=expr.index, data="Excluded", name="group", dtype="object")
    groups.loc[low_samples] = "Low"
    groups.loc[high_samples] = "High"
    excluded = groups.index[groups == "Excluded"].tolist()

    return StratificationResult(
        target_gene=target,
        expression=expr,
        groups=groups,
        high_samples=high_samples,
        low_samples=low_samples,
        excluded_samples=excluded,
        low_cutoff=low_q,
        high_cutoff=high_q,
        tie_warning=tie_warning,
    )


def run_pydeseq2(
    raw_counts: pd.DataFrame,
    stratification: StratificationResult,
    *,
    alpha: float = 0.05,
    lfc_threshold: float = 1.0,
    min_total_count: int = 10,
    n_cpus: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run High-versus-Low PyDESeq2 using raw counts only."""
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:
        raise ImportError(
            "PyDESeq2 is not installed. Install dependencies with: pip install -r requirements.txt"
        ) from exc

    selected = stratification.low_samples + stratification.high_samples
    missing = [s for s in selected if s not in raw_counts.columns]
    if missing:
        raise ValueError(f"Selected samples missing from raw counts: {missing[:8]}")

    counts_sg = raw_counts.loc[:, selected].T.copy()  # samples x genes
    keep = counts_sg.sum(axis=0) >= int(min_total_count)
    counts_sg = counts_sg.loc[:, keep]
    if counts_sg.shape[1] < 2:
        raise ValueError("Fewer than two genes remain after low-count filtering.")

    metadata = pd.DataFrame(
        {"group": ["Low"] * len(stratification.low_samples) + ["High"] * len(stratification.high_samples)},
        index=selected,
    )
    metadata["group"] = pd.Categorical(metadata["group"], categories=["Low", "High"])

    dds = DeseqDataSet(
        counts=counts_sg,
        metadata=metadata,
        design="~group",
        refit_cooks=True,
        n_cpus=max(1, int(n_cpus)),
        quiet=True,
    )
    dds.deseq2()

    stats = DeseqStats(
        dds,
        contrast=["group", "High", "Low"],
        alpha=float(alpha),
        n_cpus=max(1, int(n_cpus)),
        quiet=True,
    )
    stats.summary()
    result = stats.results_df.copy()
    result.index = result.index.astype(str)
    result.index.name = "gene"
    result = result.reset_index()

    numeric_cols = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    for col in numeric_cols:
        if col in result.columns:
            result[col] = pd.to_numeric(result[col], errors="coerce")

    result["significant"] = (
        result["padj"].notna()
        & (result["padj"] < float(alpha))
        & (result["log2FoldChange"].abs() >= float(lfc_threshold))
    )
    result = result.sort_values(
        ["significant", "padj", "log2FoldChange"],
        ascending=[False, True, False],
        na_position="last",
    ).reset_index(drop=True)
    return result, metadata


def choose_network_genes(
    deseq_results: pd.DataFrame,
    normalized: pd.DataFrame,
    target_gene: str,
    *,
    max_genes: int = 100,
    forced_partner: str | None = None,
) -> list[str]:
    """Select target A plus significant genes available in normalized expression."""
    target = resolve_gene_name(target_gene, normalized.index)
    significant = deseq_results.loc[deseq_results["significant"], :].copy()
    if significant.empty:
        raise ValueError("No significant genes are available for GENIE3.")

    significant["abs_lfc"] = significant["log2FoldChange"].abs()
    significant = significant.sort_values(["padj", "abs_lfc"], ascending=[True, False])
    available = [g for g in significant["gene"].astype(str) if g in normalized.index and g != target]

    partner: str | None = None
    if forced_partner is not None:
        partner = resolve_gene_name(forced_partner, normalized.index)
        significant_set = set(significant["gene"].astype(str))
        if partner not in significant_set:
            raise ValueError(f"Partner gene '{partner}' is not significant under the current thresholds.")
        if partner == target:
            raise ValueError("Partner gene must differ from target gene A.")

    max_genes = max(3, int(max_genes))
    chosen: list[str] = [target]
    if partner is not None:
        chosen.append(partner)
    for gene in available:
        if gene not in chosen:
            chosen.append(gene)
        if len(chosen) >= max_genes:
            break

    if len(chosen) < 3:
        raise ValueError(
            "GENIE3 requires a multigene context. Fewer than three usable genes were selected."
        )
    return chosen


def prepare_network_expression(
    normalized: pd.DataFrame,
    genes: Sequence[str],
    samples: Sequence[str],
    *,
    transform: TransformMethod = "log2p1",
) -> pd.DataFrame:
    """Prepare samples x genes normalized expression for GENIE3."""
    missing_genes = [g for g in genes if g not in normalized.index]
    missing_samples = [s for s in samples if s not in normalized.columns]
    if missing_genes:
        raise ValueError(f"Genes missing from normalized matrix: {missing_genes[:8]}")
    if missing_samples:
        raise ValueError(f"Samples missing from normalized matrix: {missing_samples[:8]}")

    expr = normalized.loc[list(genes), list(samples)].T.copy()
    expr = expr.apply(pd.to_numeric, errors="coerce")
    if expr.isna().any(axis=None):
        bad = expr.columns[expr.isna().any(axis=0)].tolist()
        raise ValueError(f"Normalized matrix contains missing/non-numeric values for network genes: {bad[:8]}")
    if not np.isfinite(expr.to_numpy(dtype=float)).all():
        raise ValueError("Normalized network expression contains infinite values.")

    if transform == "log2p1":
        if (expr.to_numpy(dtype=float) < 0).any():
            raise ValueError("log2(x+1) cannot be applied because normalized expression contains negative values.")
        expr = np.log2(expr.astype(float) + 1.0)
    elif transform != "none":
        raise ValueError(f"Unsupported transform: {transform}")

    variances = expr.var(axis=0, ddof=0)
    zero_var = variances.index[variances <= 0].tolist()
    if zero_var:
        expr = expr.drop(columns=zero_var)
    if expr.shape[1] < 3:
        raise ValueError("Fewer than three non-constant genes remain for GENIE3.")
    return expr


def infer_genie3(
    expression_sg: pd.DataFrame,
    *,
    tree_method: TreeMethod = "RF",
    n_trees: int = 1000,
    max_features: str | int = "sqrt",
    n_jobs: int = -1,
    random_state: int = 1234,
) -> Genie3Result:
    """Infer a GENIE3-style weighted directed network using tree ensembles.

    Rows of expression_sg are samples and columns are genes. For each target gene,
    all other genes are used as predictors. Feature importances form incoming edge weights.
    """
    if expression_sg.shape[0] < 6:
        raise ValueError("At least six samples are required for GENIE3; more are strongly recommended.")
    if expression_sg.shape[1] < 3:
        raise ValueError("At least three genes are required for GENIE3.")

    values = expression_sg.to_numpy(dtype=float)
    genes = expression_sg.columns.astype(str).tolist()
    p = len(genes)
    vim = np.zeros((p, p), dtype=float)

    if isinstance(max_features, str):
        if max_features not in {"sqrt", "all"}:
            raise ValueError("max_features must be 'sqrt', 'all', or a positive integer.")
        sklearn_max_features: str | float | int = "sqrt" if max_features == "sqrt" else 1.0
    else:
        sklearn_max_features = max(1, int(max_features))

    model_cls = RandomForestRegressor if tree_method == "RF" else ExtraTreesRegressor
    if tree_method not in {"RF", "ET"}:
        raise ValueError("tree_method must be 'RF' or 'ET'.")

    all_idx = np.arange(p)
    for target_idx in range(p):
        regulator_idx = all_idx[all_idx != target_idx]
        x = values[:, regulator_idx]
        y = values[:, target_idx]

        model = model_cls(
            n_estimators=max(10, int(n_trees)),
            max_features=sklearn_max_features,
            random_state=int(random_state) + target_idx,
            n_jobs=int(n_jobs),
        )
        model.fit(x, y)
        importance = np.asarray(model.feature_importances_, dtype=float)
        total = float(importance.sum())
        if total > 0:
            importance = importance / total
        vim[regulator_idx, target_idx] = importance

    vim_df = pd.DataFrame(vim, index=genes, columns=genes)
    records: list[dict[str, object]] = []
    for regulator_idx, regulator in enumerate(genes):
        for target_idx, target in enumerate(genes):
            if regulator_idx == target_idx:
                continue
            records.append(
                {
                    "regulator": regulator,
                    "target": target,
                    "weight": float(vim[regulator_idx, target_idx]),
                }
            )
    edges = pd.DataFrame.from_records(records).sort_values("weight", ascending=False).reset_index(drop=True)
    edges["global_rank"] = np.arange(1, len(edges) + 1)
    edges["target_rank"] = edges.groupby("target")["weight"].rank(method="min", ascending=False).astype(int)
    denom = max(1, p - 2)
    edges["target_rank_percentile"] = 1.0 - (edges["target_rank"] - 1) / denom

    return Genie3Result(vim=vim_df, edges=edges, expression_used=expression_sg.copy())


def extract_pair_result(
    result: Genie3Result,
    gene_a: str,
    gene_b: str,
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Extract A→B and B→A weights and an unsigned association diagnostic."""
    genes = result.expression_used.columns.tolist()
    a = resolve_gene_name(gene_a, genes)
    b = resolve_gene_name(gene_b, genes)
    if a == b:
        raise ValueError("gene_a and gene_b must differ.")

    pair_edges = result.edges.loc[
        ((result.edges["regulator"] == a) & (result.edges["target"] == b))
        | ((result.edges["regulator"] == b) & (result.edges["target"] == a))
    ].copy()
    pair_edges.insert(0, "edge", pair_edges["regulator"] + " → " + pair_edges["target"])
    pair_edges = pair_edges.sort_values("weight", ascending=False).reset_index(drop=True)

    rho, pvalue = spearmanr(
        result.expression_used[a].to_numpy(dtype=float),
        result.expression_used[b].to_numpy(dtype=float),
    )
    diagnostics = {
        "spearman_rho": float(rho),
        "spearman_pvalue": float(pvalue),
        "n_samples": float(result.expression_used.shape[0]),
        "n_network_genes": float(result.expression_used.shape[1]),
    }
    return pair_edges, diagnostics


def target_neighbor_tables(
    result: Genie3Result,
    target_gene: str,
    top_n: int = 15,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    target = resolve_gene_name(target_gene, result.expression_used.columns)
    downstream = result.edges.loc[result.edges["regulator"] == target].head(int(top_n)).copy()
    upstream = result.edges.loc[result.edges["target"] == target].head(int(top_n)).copy()
    return downstream, upstream
        
        
        
