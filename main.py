from __future__ import annotations

import io

import pandas as pd
import streamlit as st

from analysis_core import (
    align_samples,
    choose_network_genes,
    clean_expression_matrix,
    extract_pair_result,
    infer_genie3,
    prepare_network_expression,
    resolve_gene_name,
    run_pydeseq2,
    stratify_samples_by_target,
    target_neighbor_tables,
    validate_raw_counts,
)
from plotting import plot_pair_scatter, plot_stratification, plot_target_network


st.set_page_config(page_title="Target-gene stratified DESeq2 + GENIE3", layout="wide")

# ====================== 【核心修复：防止下载刷新丢失session】 ======================
# 预初始化session，杜绝None刷新清空
if "de_bundle" not in st.session_state:
    st.session_state["de_bundle"] = None
if "genie3_bundle" not in st.session_state:
    st.session_state["genie3_bundle"] = None

# JS拦截下载按钮刷新行为（根治跳转首页）
st.markdown("""
<script>
document.addEventListener("click", function(e){
    const btn = e.target.closest('button[data-testid="stDownloadButton"]');
    if(btn){
        e.stopImmediatePropagation();
        e.preventDefault();
    }
})
</script>
""", unsafe_allow_html=True)
# =================================================================================

st.title("Target-gene stratified differential analysis and GENIE3")
st.caption(
    "Normalized expression defines target-high/target-low samples; raw counts are used only for DESeq2; "
    "normalized expression is used for GENIE3."
)

with st.expander("Required input format", expanded=False):
    st.markdown(
        """
Both files must be CSV matrices with **genes in rows**, **samples in columns**, and the first column containing gene IDs.

- Normalized matrix: TPM, FPKM, or another continuous normalized expression matrix.
- Raw count matrix: non-negative integer read counts.
- Sample IDs must match between the two files.
- Gene IDs should use the same identifier system whenever possible.
        """
    )

left, right = st.columns(2)
with left:
    normalized_file = st.file_uploader("1. Upload normalized expression matrix", type=["csv"], key="norm")
with right:
    counts_file = st.file_uploader("2. Upload raw count matrix", type=["csv"], key="counts")


def read_uploaded_csv(uploaded) -> pd.DataFrame:
    uploaded.seek(0)
    return pd.read_csv(uploaded, index_col=0)


if normalized_file is None or counts_file is None:
    st.info("Upload both matrices to begin.")
    st.stop()

try:
    normalized = clean_expression_matrix(
        read_uploaded_csv(normalized_file),
        matrix_name="Normalized expression matrix",
        duplicate_gene_rule="mean",
    )
    raw_counts = clean_expression_matrix(
        read_uploaded_csv(counts_file),
        matrix_name="Raw count matrix",
        duplicate_gene_rule="sum",
    )
    raw_counts = validate_raw_counts(raw_counts)
    normalized, raw_counts, shared_samples, norm_only, count_only = align_samples(normalized, raw_counts)
except Exception as exc:
    st.error(str(exc))
    st.stop()

st.success(
    f"Loaded {normalized.shape[0]:,} normalized-expression genes, {raw_counts.shape[0]:,} count genes, "
    f"and {len(shared_samples)} shared samples."
)
if norm_only or count_only:
    st.warning(
        f"Non-shared samples were excluded: normalized-only={len(norm_only)}, count-only={len(count_only)}."
    )

st.sidebar.header("Stratification and DESeq2")
target_gene_input = st.sidebar.text_input("Target gene A", value="")
fraction = st.sidebar.slider("High/low fraction per group", 0.10, 0.40, 0.25, 0.05, help=(
        "The default is the top and bottom 25%. "
        "Adjust the proportion according to cohort size and "
        "the expression distribution of target gene A."),
    )
group_method_label = st.sidebar.selectbox(
    "Grouping method",
    ["Exact rank groups", "Quantile thresholds"],
    index=0,
)
group_method = "rank" if group_method_label == "Exact rank groups" else "quantile"
alpha = st.sidebar.number_input("Adjusted P-value threshold", 0.001, 0.20, 0.05, 0.005, format="%.3f")
lfc_threshold = st.sidebar.number_input("|log2FC| threshold", 0.0, 5.0, 1.0, 0.1)
min_total_count = st.sidebar.number_input("Minimum total raw count", 0, 10000, 10, 1)
n_cpus = st.sidebar.number_input("CPU workers for DESeq2", 1, 32, 1, 1)

run_de = st.sidebar.button("Run stratification and DESeq2", type="primary")

if run_de:
    try:
        target_gene = resolve_gene_name(target_gene_input, normalized.index)
        if target_gene not in raw_counts.index:
            raise ValueError(
                f"Target gene '{target_gene}' is present in the normalized matrix but absent from raw counts."
            )
        stratification = stratify_samples_by_target(
            normalized,
            target_gene,
            fraction=float(fraction),
            method=group_method,
        )
        de_results, metadata = run_pydeseq2(
            raw_counts,
            stratification,
            alpha=float(alpha),
            lfc_threshold=float(lfc_threshold),
            min_total_count=int(min_total_count),
            n_cpus=int(n_cpus),
        )
        st.session_state["de_bundle"] = {
            "target_gene": target_gene,
            "stratification": stratification,
            "de_results": de_results,
            "metadata": metadata,
            "normalized": normalized,
            "raw_counts": raw_counts,
            "shared_samples": shared_samples,
            "alpha": float(alpha),
            "lfc_threshold": float(lfc_threshold),
        }
        st.session_state.pop("genie3_bundle", None)
    except Exception as exc:
        st.error(str(exc))

bundle = st.session_state.get("de_bundle")
if bundle is None:
    st.stop()

stratification = bundle["stratification"]
de_results = bundle["de_results"]
target_gene = bundle["target_gene"]
normalized = bundle["normalized"]
shared_samples = bundle["shared_samples"]

st.header("1. Target-gene sample stratification")
col1, col2, col3 = st.columns(3)
col1.metric("A-low samples", len(stratification.low_samples))
col2.metric("A-high samples", len(stratification.high_samples))
col3.metric("Excluded middle samples", len(stratification.excluded_samples))
if stratification.tie_warning:
    st.warning(stratification.tie_warning)

st.pyplot(plot_stratification(stratification), clear_figure=True)

group_export = pd.DataFrame(
    {
        "sample": stratification.expression.index,
        f"{target_gene}_expression": stratification.expression.values,
        "group": stratification.groups.reindex(stratification.expression.index).values,
    }
)
# 【修复】下载防刷新
st.download_button(
    "Download sample stratification",
    data=group_export.to_csv(index=False).encode("utf-8"),
    file_name=f"{target_gene}_sample_stratification.csv",
    mime="text/csv",
    on_click=lambda: None
)

st.header("2. DESeq2: A-high versus A-low")
st.caption("Positive log2FC means higher expression in the A-high group.")
sig_count = int(de_results["significant"].sum())
st.metric("Significant genes", sig_count)
st.dataframe(de_results, use_container_width=True, height=420)

# 【修复】下载防刷新
st.download_button(
    "Download complete DESeq2 results",
    data=de_results.to_csv(index=False).encode("utf-8"),
    file_name=f"{target_gene}_high_vs_low_pydeseq2.csv",
    mime="text/csv",
    on_click=lambda: None
)

significant_available = [
    gene
    for gene in de_results.loc[de_results["significant"], "gene"].astype(str)
    if gene in normalized.index and gene != target_gene
]
if not significant_available:
    st.warning("No significant non-target genes are shared with the normalized matrix; GENIE3 cannot proceed.")
    st.stop()

st.header("3. GENIE3 and target-pair extraction")
st.info(
    "GENIE3 is fitted in a multigene context (A plus significant genes). The app then extracts A→B and B→A. "
    "A two-gene-only GENIE3 model is intentionally not used because its feature importance would be uninformative."
)

c1, c2, c3 = st.columns(3)
with c1:
    partner_gene = st.selectbox("Select significant partner gene B", significant_available)
    max_genes = st.number_input("Maximum network genes (including A and B)", 3, 500, 100, 10)
with c2:
    network_sample_mode = st.selectbox(
        "Samples used for GENIE3",
        [
            "Only A-high and A-low samples (requested workflow)",
            "All shared samples (recommended sensitivity analysis)",
        ],
    )
    transform_label = st.selectbox(
        "Normalized-expression transform for GENIE3",
        ["log2(x + 1) — use for unlogged TPM/FPKM", "No transform — use for already log-transformed data"],
    )
with c3:
    tree_method = st.selectbox("Tree ensemble", ["RF", "ET"], format_func=lambda x: "Random Forest" if x == "RF" else "Extra Trees")
    n_trees = st.number_input("Number of trees per target", 50, 5000, 1000, 50)
    random_state = st.number_input("Random seed", 0, 1000000, 1234, 1)

extra_regulator_text = st.text_area(
    "Optional expressed TFs / pathway regulators",
    value="",
    help=(
        "Optional gene symbols separated by commas. "
        "Only genes present in the normalized matrix will be retained."
    ),
)

extra_regulators = [
    gene.strip()
    for gene in extra_regulator_text.replace("\n", ",").split(",")
    if gene.strip()
]

run_genie3 = st.button("Run GENIE3", type="primary")
if run_genie3:
    try:
        network_genes = choose_network_genes(
            de_results,
            normalized,
            target_gene,
            max_genes=int(max_genes),
            forced_partner=partner_gene,
            extra_regulators=extra_regulators,
        )
        if network_sample_mode.startswith("Only"):
            network_samples = stratification.low_samples + stratification.high_samples
        else:
            network_samples = shared_samples
        transform = "log2p1" if transform_label.startswith("log2") else "none"
        network_expression = prepare_network_expression(
            normalized,
            network_genes,
            network_samples,
            transform=transform,
        )
        if target_gene not in network_expression.columns or partner_gene not in network_expression.columns:
            raise ValueError(
                "A or B had zero variance or invalid values after network preprocessing and was removed."
            )
        genie3_result = infer_genie3(
            network_expression,
            tree_method=tree_method,
            n_trees=int(n_trees),
            n_jobs=-1,
            random_state=int(random_state),
        )
        pair_edges, diagnostics = extract_pair_result(genie3_result, target_gene, partner_gene)
        downstream, upstream = target_neighbor_tables(genie3_result, target_gene, top_n=15)
        st.session_state["genie3_bundle"] = {
            "result": genie3_result,
            "pair_edges": pair_edges,
            "diagnostics": diagnostics,
            "downstream": downstream,
            "upstream": upstream,
            "target_gene": target_gene,
            "partner_gene": partner_gene,
            "sample_mode": network_sample_mode,
            "network_genes": network_genes,
        }
    except Exception as exc:
        st.error(str(exc))

g_bundle = st.session_state.get("genie3_bundle")
if g_bundle is None:
    st.stop()

result = g_bundle["result"]
pair_edges = g_bundle["pair_edges"]
diagnostics = g_bundle["diagnostics"]
partner_gene = g_bundle["partner_gene"]

st.subheader(f"Predicted relationship: {target_gene} and {partner_gene}")
metric1, metric2, metric3 = st.columns(3)
metric1.metric("Spearman ρ", f"{diagnostics['spearman_rho']:.3f}")
metric2.metric("Spearman P", f"{diagnostics['spearman_pvalue']:.3g}")
metric3.metric("GENIE3 samples", int(diagnostics["n_samples"]))

st.dataframe(
    pair_edges[["edge", "weight", "target_rank", "target_rank_percentile", "global_rank"]],
    use_container_width=True,
)
st.caption(
    "GENIE3 weights are non-negative predictive-importance scores. They do not by themselves establish causality, "
    "direct binding, activation, or repression. Spearman correlation is shown separately only as a sign/association diagnostic."
)

p1, p2 = st.columns(2)
with p1:
    st.pyplot(plot_pair_scatter(result.expression_used, target_gene, partner_gene), clear_figure=True)
with p2:
    st.pyplot(plot_target_network(result, target_gene), clear_figure=True)

st.subheader(f"Top predicted outgoing edges from {target_gene}")
st.dataframe(g_bundle["downstream"], use_container_width=True)
st.subheader(f"Top predicted incoming edges to {target_gene}")
st.dataframe(g_bundle["upstream"], use_container_width=True)

# 【修复】全部GENIE3下载防刷新
st.download_button(
    "Download all GENIE3 edges",
    data=result.edges.to_csv(index=False).encode("utf-8"),
    file_name=f"{target_gene}_GENIE3_edges.csv",
    mime="text/csv",
    on_click=lambda: None
)
st.download_button(
    "Download A–B pair result",
    data=pair_edges.to_csv(index=False).encode("utf-8"),
    file_name=f"{target_gene}_{partner_gene}_GENIE3_pair.csv",
    mime="text/csv",
    on_click=lambda: None
)
st.download_button(
    "Download GENIE3 weight matrix",
    data=result.vim.to_csv().encode("utf-8"),
    file_name=f"{target_gene}_GENIE3_weight_matrix.csv",
    mime="text/csv",
    on_click=lambda: None
)
