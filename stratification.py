import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import streamlit as st
from GENIE3 import GENIE3
from pydeseq2.dds import DeseqDataSet as Deseq2DataSet
from pydeseq2.ds import DeseqStats as DeseqStats

# Force global UTF-8 encoding
os.environ['PYTHONIOENCODING'] = 'utf-8'
os.environ['LC_ALL'] = 'en_US.UTF-8'

if sys.platform == "win32":
    import multiprocessing
    multiprocessing.set_start_method("spawn", force=True)


# ========== Temporary Folder Permission Check ==========
TEMP_FOLDER = "E:/cell_definition_temp"
if not os.path.exists(TEMP_FOLDER):
    os.makedirs(TEMP_FOLDER, exist_ok=True)
if not os.access(TEMP_FOLDER, os.W_OK):
    st.error(f"""
        ❌ Error: No write permission for the temporary folder **{TEMP_FOLDER}**!  
        Solutions:  
        1. Manually create the folder and ensure the current user has write permission;  
        2. Change the path to an English-only directory without spaces (e.g., `E:/temp`).  
    """)
    st.stop()
os.environ["JOBLIB_TEMP_FOLDER"] = TEMP_FOLDER


# ========== Gene Stratification Function ==========
def stratify_genes(counts, high_quantile=0.75, low_quantile=0.25, save_path="gene_stratification_plot.png"):
    gene_max = counts.max(axis=1)
    threshold_high = gene_max.quantile(high_quantile)
    threshold_low = gene_max.quantile(low_quantile)

    high_expr = gene_max[gene_max > threshold_high].index.tolist()
    mid_expr = gene_max[(gene_max > threshold_low) & (gene_max <= threshold_high)].index.tolist()
    low_expr = gene_max[gene_max <= threshold_low].index.tolist()

    plt.figure(figsize=(10, 6))
    sns.histplot(gene_max, kde=True, bins=30, color='#4CAF50')
    plt.axvline(threshold_high, color='r', linestyle='--', label=f'High Expression Threshold ({high_quantile * 100}% Quantile)')
    plt.axvline(threshold_low, color='b', linestyle='--', label=f'Low Expression Threshold ({low_quantile * 100}% Quantile)')
    plt.axvspan(threshold_high, gene_max.max(), alpha=0.2, color='red', label='High Expression')
    plt.axvspan(threshold_low, threshold_high, alpha=0.2, color='yellow', label='Medium Expression')
    plt.axvspan(gene_max.min(), threshold_low, alpha=0.2, color='blue', label='Low Expression')
    plt.title('Gene Expression Stratification (High, Medium, Low Layers)')
    plt.xlabel('Maximum Gene Expression (count value)')
    plt.ylabel('Number of Genes')
    plt.legend(loc='upper right')
    plt.tight_layout()
    st.pyplot(plt)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    st.success(f"Gene stratification plot saved to: {save_path}")

    st.info(f"Gene Stratification Results:\n"
            f"- High-expression genes: {len(high_expr)} (> {threshold_high:.2f})\n"
            f"- Medium-expression genes: {len(mid_expr)} ({threshold_low:.2f}-{threshold_high:.2f})\n"
            f"- Low-expression genes: {len(low_expr)} (< {threshold_low:.2f})")

    return high_expr, mid_expr, low_expr, threshold_high, threshold_low


# ========== Differential Analysis Function ==========
def run_pydeseq2_analysis(counts, metadata, group_col='group',
                          log2fc_thresh=0.5,
                          p_thresh=0.1,
                          n_jobs=4):

    counts = counts.astype(int)

    # Data Quality Preview (aids in troubleshooting)
    st.subheader("Data Preview (First 5 Rows)")
    st.dataframe(counts.head())
    st.dataframe(metadata.head())

    # Check if group column exists
    if group_col not in metadata.columns:
        st.error(f"❌ Error: Group column `{group_col}` not found in metadata!")
        return None, []

    # Check number of group levels (must be exactly 2: treatment group + control group)
    group_levels = metadata[group_col].unique()
    if len(group_levels) != 2:
        st.error(f"""
            ❌ Error: Group column `{group_col}` must contain **2 levels** (treatment group + control group), but {len(group_levels)} levels found:  
            {group_levels}  
            Solution: Ensure this column in metadata has only two unique values (e.g., 'Treat' and 'Ctrl').
        """)
        return None, []

    # Automatically infer treatment and control groups (logic can be adjusted based on scenarios)
    group_levels_sorted = sorted(group_levels)
    ctrl_level, treat_level = group_levels_sorted  # Example: ['Ctrl', 'Treat'] → Ctrl = control, Treat = treatment
    # (Optional: Manual specification is more reliable, e.g., treat_level = 'Treat', ctrl_level = 'Ctrl')

    # Construct DESeq2 dataset
    dds = Deseq2DataSet(
        counts=counts.T,  # Transpose to samples × genes
        metadata=metadata,
        design_factors=group_col,  # Name of the group column (e.g., 'group')
        refit_cooks=True  # Handle outliers
    )

    # Execute core DESeq2 analysis
    try:
        dds.deseq2()
    except Exception as e:
        st.error(f"DESeq2 Analysis Failed: {str(e)}")
        st.warning("Troubleshooting Suggestions:\n"
                   "1. Check if gene names/paths contain Chinese characters\n"
                   "2. Run `pip check` to verify dependency conflicts\n"
                   "3. Upgrade dependencies: `pip install --upgrade pydeseq2 numpy scipy`")
        return None, []

    # ========== Key Fix: Explicitly Define contrast Parameter ==========
    try:
        # contrast format: [group column name, treatment group level, control group level]
        contrast = [group_col, treat_level, ctrl_level]
        stat_res = DeseqStats(dds, contrast=contrast)
        stat_res.summary()  # Calculate p-values and padj
    except Exception as e:
        st.error(f"DeseqStats Initialization Failed: {str(e)}")
        st.warning(f"""
            Possible Causes:  
            1. Incorrect contrast format (required: `[ '{group_col}', '{treat_level}', '{ctrl_level}' ]`)  
            2. Group levels do not match actual values in metadata (case-sensitive!)  
            Currently inferred contrast: {contrast}  
            Actual group levels: {group_levels}
        """)
        return None, []
    # ========== Fix Completed ==========

    # Extract differential analysis results (fix column names)
    res_df = stat_res.results_df.reset_index()

    # Fix abnormal column names: convert 0 or 'index' to 'gene'
    if 0 in res_df.columns:
        res_df = res_df.rename(columns={0: 'gene'})
    else:
        res_df = res_df.rename(columns={'index': 'gene'})
    st.write("✅ Gene Column Name:", res_df.columns)

    # ========== Deep Cleaning of Gene Column (Adapt to New Arrow Version) ==========
    res_df['gene'] = res_df['gene'].astype(str)  # Force string type
    res_df['gene'] = res_df['gene'].str.strip()  # Remove extra spaces
    res_df['gene'] = res_df['gene'].replace(
        r'[^\w\-_.]', '', regex=True  # Remove special characters (retain letters, numbers, -, _, .)
    )
    res_df.loc[res_df['gene'] == '', 'gene'] = 'unknown_gene'  # Fill empty values
    # ========== Cleaning Completed ==========

    # Fix numeric columns (continue previous type conversion code)
    numeric_cols = ['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
    for col in numeric_cols:
        if col in res_df.columns and res_df[col].dtype == 'object':
            res_df[col] = pd.to_numeric(res_df[col], errors='coerce').fillna(0)

    # Verify Arrow Compatibility
    try:
        import pyarrow as pa
        pa.Table.from_pandas(res_df)
        st.success("✅ Data can be normally converted to Arrow table; Streamlit can display it properly!")
    except Exception as e:
        st.error(f"❌ Arrow Conversion Still Fails: {str(e)}")
        st.warning("Suggestions:\n"
                   "1. Check if the gene column contains special characters (e.g., \\, |, tabs)\n"
                   "2. Downgrade Streamlit to version 1.27 (run `pip install streamlit==1.27`)")

    # Filter Significant Differential Genes
    significant_degs = res_df[
        (abs(res_df['log2FoldChange']) > log2fc_thresh) &
        (res_df['padj'] < p_thresh)
    ].sort_values('padj')['gene'].tolist()

    # Display Results (remove use_arrow to adapt to new Streamlit version)
    st.subheader("Differential Analysis Results Preview (First 10 Rows)")
    st.dataframe(
        res_df[['gene', 'baseMean', 'log2FoldChange', 'padj']]
        .sort_values('padj')
        .head(10)
    )

# Visualize Distribution of Differential Indicators
    st.subheader("Distribution of Differential Indicators (pyDESeq2 Results)")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    # Distribution of log2FoldChange
    sns.histplot(res_df['log2FoldChange'], kde=True, ax=ax1, color='blue')
    ax1.axvline(log2fc_thresh, color='r', linestyle='--', label=f'Threshold={log2fc_thresh}')
    ax1.axvline(-log2fc_thresh, color='r', linestyle='--')
    ax1.set_title('Distribution of log2FoldChange')
    ax1.legend()
    # Distribution of padj
    sns.histplot(res_df['padj'], kde=True, ax=ax2, color='green')
    ax2.axvline(p_thresh, color='r', linestyle='--', label=f'Threshold={p_thresh}')
    ax2.set_title('Distribution of Adjusted p-value (padj)')
    ax2.legend()
    plt.tight_layout()
    st.pyplot(plt)
    plt.close()

    # Display Results Preview
    st.subheader("Differential Analysis Results Preview (First 10 Rows)")
    st.dataframe(
        res_df[['gene', 'baseMean', 'log2FoldChange', 'padj']]
        .sort_values('padj')
        .head(10)
    )

    if not significant_degs:
        st.warning("No significant differential genes found! You can try:\n"
                   "1. Lower the log2FC threshold (e.g., 0.3)\n"
                   "2. Increase the padj threshold (e.g., 0.15)\n"
                   "3. Check if group levels are correct (treatment group vs control group)")
        return res_df, []

    elif len(significant_degs) < 100:
        st.warning(f"⚠️ Total number of significant differential genes is only {len(significant_degs)}, less than 100! Cannot build regulatory network!")

    res_df = res_df.sort_values('pvalue', ascending=True)
    st.success(f"pyDESeq2 Analysis Completed: A total of {len(significant_degs)} significant differential genes found (sorted by p-value)")
    res_df.to_csv('pydeseq2_results.csv', index=False, encoding='utf-8')
    return res_df, significant_degs


# ========== Regulatory Network Construction ==========
def build_regulatory_network(counts, significant_degs, max_genes, top_edges=300,
                             edge_info_path="gene_regulatory_edges.txt"):
    if not significant_degs:
        st.warning("No significant differential genes; cannot build network!")
        return None
    analysis_genes = significant_degs[:max_genes]
    st.info(f"Using top {len(analysis_genes)} significant differential genes to build global regulatory network")
    # Filter genes not present in the count matrix
    analysis_genes = [g for g in analysis_genes if g in counts.index]
    if len(analysis_genes) < 10:
        st.error("Fewer than 10 valid analysis genes; cannot build network!")
        return None
    expr_matrix = counts.loc[analysis_genes].T.astype(float).values
    try:
        st.info("GENIE3 is predicting global regulatory relationships (about 5-10 minutes, please be patient...)")
        vim_matrix = GENIE3(expr_matrix, gene_names=analysis_genes)
    except Exception as e:
        st.error(f"GENIE3 Run Failed: {str(e)}")
        return None
    edges = []
    for i, target in enumerate(analysis_genes):
        for j, regulator in enumerate(analysis_genes):
            if i == j:
                continue
            weight = vim_matrix[i, j]
            edges.append((regulator, target, weight))
    edges.sort(key=lambda x: x[2], reverse=True)
    top_edges = edges[:top_edges]
    # Save global regulatory edge information
    with open(edge_info_path, 'w', encoding='utf-8') as f:
        f.write("Regulator\tTarget Gene\tWeight\n")
        for reg, target, weight in top_edges:
            f.write(f"{reg}\t{target}\t{weight:.6f}\n")
    st.success(f"Global regulatory edge information saved to: {edge_info_path}")
    G = nx.DiGraph()
    for reg, target, weight in top_edges:
        G.add_edge(reg, target, weight=round(weight, 4))
    st.success(f"Global regulatory network constructed: {len(G.nodes)} nodes, {len(G.edges)} edges")
    return G


# ========== Hub Gene Report Function ==========
def generate_hub_report(G, counts, save_path="hub_gene_report.txt"):
    if not G:
        st.warning("Regulatory network is empty; cannot generate hub gene report!")
        return "", ""
    hub_gene = max(G.nodes, key=lambda x: G.degree(x))
    upstream_perturbed, downstream_perturbed = simulate_perturbation(counts, hub_gene, direction='up')
    edge_data = G.edges(hub_gene, data=True)
    total_edges = len(edge_data)
    activate_count = sum(1 for _, _, d in edge_data if d['weight'] > 0)
    activate_ratio = (activate_count / total_edges * 100) if total_edges > 0 else 0.0
    regulatory_mode = f"Activation Ratio {activate_ratio:.1f}%, Inhibition Ratio {100 - activate_ratio:.1f}%"
    upstream_lines = []
    for gene in upstream_perturbed[:5]:
        truncated_gene = gene[:15].ljust(15)
        upstream_lines.append(f"│ • {truncated_gene} │")
    downstream_lines = []
    for gene in downstream_perturbed[:5]:
        truncated_gene = gene[:15].ljust(15)
        downstream_lines.append(f"│ • {truncated_gene} │")
    report_parts = [
        "┌───────────────────────┐",
        f"│ Global Hub Gene: {hub_gene.ljust(10)} │",
        "├───────────────────────┤",
        f"│ Upstream Regulators({len(upstream_perturbed)}) │",
        *upstream_lines,
        "├───────────────────────┤",
        f"│ Downstream Targets({len(downstream_perturbed)}) │",
        *downstream_lines,
        "├───────────────────────┤",
        "│ Main Regulatory Mode: │",
        f"│ {regulatory_mode.ljust(30)} │",
        "└───────────────────────┘"
    ]
    report = "\n".join(report_parts)
    st.markdown(f"```\n{report}\n```")
    with open(save_path, 'w', encoding='utf-8') as f:
        f.write(report)
    st.success(f"Global hub gene report saved to: {save_path} (Gene: {hub_gene})")
    return hub_gene, report


# ========== Virtual Perturbation Function ==========
def simulate_perturbation(counts, gene, direction='up', factor=2):
    perturbed_counts = counts.copy()
    if direction == 'up':
        perturbed_counts.loc[gene] *= factor
    else:
        perturbed_counts.loc[gene] /= factor
    corr = perturbed_counts.corrwith(perturbed_counts.loc[gene], axis=1).drop(gene)
    upstream = [f"{g} (+)" if c > 0 else f"{g} (-)" for g, c in corr[corr > 0.7].items()]
    downstream = [f"{g} (+)" if c > 0 else f"{g} (-)" for g, c in corr[corr < -0.7].items()]
    return upstream, downstream


# ========== Upstream and Downstream Analysis Function for Specific Gene ==========
def get_specific_gene_relations(G, gene):
    """Extract upstream regulators and downstream target genes of a specific gene"""
    if gene not in G.nodes:
        return [], []

    upstream = [(reg, data['weight']) for reg, target, data in G.edges(data=True) if target == gene]
    downstream = [(target, data['weight']) for reg, target, data in G.edges(data=True) if reg == gene]

    upstream.sort(key=lambda x: x[1], reverse=True)
    downstream.sort(key=lambda x: x[1], reverse=True)
    return upstream, downstream


def generate_specific_report(upstream, downstream, gene, save_path="specific_gene_report.txt"):
    """Generate upstream and downstream regulatory report for a specific gene"""
    upstream_lines = [f"│ • {r[:15].ljust(15)} (Weight: {w:.4f}) │" for r, w in upstream[:10]]
    downstream_lines = [f"│ • {t[:15].ljust(15)} (Weight: {w:.4f}) │" for t, w in downstream[:10]]
    report_parts = [
        "┌───────────────────────┐",
        f"│ Specific Gene of Interest: {gene.ljust(10)} │",
        "├───────────────────────┤",
        f"│ Upstream Regulators({len(upstream)}) │",
        *upstream_lines,
        "├───────────────────────┤",
        f"│ Downstream Targets({len(downstream)}) │",
        *downstream_lines,
        "└───────────────────────┘"
    ]
    report = "\n".join(report_parts)
    st.markdown(f"```\n{report}\n```")
    with open(save_path, 'w', encoding='utf-8') as f:
        f.write(report)
    st.success(f"Specific gene report saved to: {save_path}")

# ========== Main Workflow ==========
def main():
    st.set_page_config(page_title="Gene Expression Analysis Workflow (Supports 500 Genes)", layout="wide")
    st.title("Count Matrix Analysis (pyDESeq2 + GENIE3)")

    # Initialize session state (record last max_genes to detect parameter changes)
    if 'network' not in st.session_state:
        st.session_state['network'] = None  # Store current regulatory network
    if 'hub_report_saved' not in st.session_state:
        st.session_state['hub_report_saved'] = False  # Mark if hub report has been generated
    if 'strat_plot_path' not in st.session_state:
        st.session_state['strat_plot_path'] = "gene_stratification_plot.png"
    if 'last_max_genes' not in st.session_state:
        st.session_state['last_max_genes'] = 500  # Record last number of analysis genes, default 500

    # 1. Data Upload
    st.sidebar.header("Data Upload")
    treat_file = st.sidebar.file_uploader("Treatment Group Count Matrix (CSV)", type="csv")
    ctrl_file = st.sidebar.file_uploader("Control Group Count Matrix (CSV)", type="csv")

    if not treat_file or not ctrl_file:
        st.warning("Please first upload count matrices for the treatment and control groups (when there is no header, the first column is gene names)")
        return

    # 2. Data Reading and Merging
    try:
        has_header = st.sidebar.checkbox("CSV Contains Header", value=False)
        treat_counts = pd.read_csv(treat_file, header=0 if has_header else None, index_col=0)
        ctrl_counts = pd.read_csv(ctrl_file, header=0 if has_header else None, index_col=0)

        treat_samples = [f"Treat_{i + 1}" for i in range(treat_counts.shape[1])]
        ctrl_samples = [f"Ctrl_{i + 1}" for i in range(ctrl_counts.shape[1])]
        treat_counts.columns = treat_samples
        ctrl_counts.columns = ctrl_samples

        counts = pd.merge(treat_counts, ctrl_counts, left_index=True, right_index=True, how='inner')
        st.success(
            f"Data Merged Successfully:\n"
            f"- Total Genes: {counts.shape[0]}\n"
            f"- Treatment Group Samples: {len(treat_samples)}\n"
            f"- Control Group Samples: {len(ctrl_samples)}"
        )
    except Exception as e:
        st.error(f"Data Reading Failed: {str(e)}")
        st.error("Please confirm the CSV format: when there is no header, the first column is gene names, and subsequent columns are sample expression values (pure numbers)")
        return

    # 3. Construct Group Metadata
    group_col = 'group'
    metadata = pd.DataFrame({
        'sample': treat_samples + ctrl_samples,
        group_col: ['Treat'] * len(treat_samples) + ['Ctrl'] * len(ctrl_samples)
    }).set_index('sample')

    # 4. Gene Stratification + Stratification Plot Download
    st.header("1. Gene Expression Stratification")
    col1, col2 = st.columns(2)
    high_quantile = col1.slider("High-Expression Gene Quantile", 0.5, 0.95, 0.75, 0.05)
    low_quantile = col2.slider("Low-Expression Gene Quantile", 0.05, 0.5, 0.25, 0.05)

    # Perform stratification and save plot (store path in session state)
    high_expr, mid_expr, low_expr, _, _ = stratify_genes(
        counts, high_quantile, low_quantile, save_path=st.session_state['strat_plot_path']
    )

    # Button to download gene stratification plot
    if os.path.exists(st.session_state['strat_plot_path']):
        with open(st.session_state['strat_plot_path'], "rb") as f:
            st.download_button(
                label="Download Gene Stratification Plot (PNG)",
                data=f,
                file_name="gene_stratification_plot.png",
                mime="image/png"
            )

    # 5. Differential Expression Analysis + Results Download
    st.header("2. Differential Expression Analysis (pyDESeq2)")
    hl_counts = counts.loc[high_expr + low_expr]
    st.info(f"Differential Analysis Scope: High-Expression({len(high_expr)}) + Low-Expression({len(low_expr)}) = Total {len(hl_counts)} genes")

    st.sidebar.header("Differential Analysis Parameters")
    log2fc_thresh = st.sidebar.slider("log2FC Threshold (smaller = more significant genes)", 0.1, 2.0, 0.5)
    p_thresh = st.sidebar.slider("padj Threshold (larger = more significant genes)", 0.01, 0.2, 0.05)

    # Perform differential analysis (generate pydeseq2_results.csv)
    deg_df, significant_degs = run_pydeseq2_analysis(
        hl_counts,
        metadata,
        group_col=group_col,
        log2fc_thresh=log2fc_thresh,
        p_thresh=p_thresh
    )

    if not significant_degs:
        return

    # Button to download differential analysis results (CSV file)
    if os.path.exists("pydeseq2_results.csv"):
        with open("pydeseq2_results.csv", "rb") as f:
            st.download_button(
                label="Download Differential Analysis Results (CSV)",
                data=f,
                file_name="pydeseq2_results.csv",
                mime="text/csv"
            )

    # 6. Regulatory Network Construction
    st.header("3. Regulatory Network Construction (GENIE3)")
    st.sidebar.header("GENIE3 Network Parameters")
    # Slider to select new max_genes
    current_max_genes = st.sidebar.slider("Total Number of Genes for Analysis", 100, 2000, 500)
    top_edges = st.sidebar.slider("Number of High-Weight Regulatory Edges to Retain", 100, 500, 300)

    # Detect if current max_genes matches the last one; clear old network if not
    if current_max_genes != st.session_state['last_max_genes']:
        st.warning(f"⚠️ The number of genes for analysis has changed from {st.session_state['last_max_genes']} to {current_max_genes}; the network will be reconstructed")
        st.session_state['network'] = None  # Clear old network
        st.session_state['hub_report_saved'] = False  # Reset hub report marker
        st.session_state['last_max_genes'] = current_max_genes  # Update last max_genes

    if st.session_state['network'] is None:
        st.info(f"Constructing regulatory network for {current_max_genes} genes...")
        st.session_state['network'] = build_regulatory_network(
            counts,
            significant_degs,
            max_genes=current_max_genes,  # Use currently selected number of genes
            top_edges=top_edges
        )
    else:
        st.success(f"✅ Reused regulatory network for {current_max_genes} genes (no need to reconstruct)")

    if not st.session_state['network']:
        return

    # 7. Global Hub Gene Report
    st.header("4. Global Hub Genes and Virtual Perturbation")
    if not st.session_state['hub_report_saved']:
        generate_hub_report(st.session_state['network'], counts)
        st.session_state['hub_report_saved'] = True
    else:
        st.success("✅ Hub gene report has been generated (based on the current network)")
        # Display the report of the current network (avoid residual old reports)
        with open("hub_gene_report.txt", "r", encoding='utf-8') as f:
            hub_report = f.read()
        st.markdown(f"```\n{hub_report}\n```")

    # 8. Specific Gene Analysis
    st.header("5. Upstream and Downstream Regulatory Research for Specific Genes (ensure the gene is within the analysis scope)")
    specific_gene = st.text_input("Please enter the name of the specific gene to study", "")
    specific_report_path = "specific_gene_report.txt"

    if specific_gene:
        # Validation
        if specific_gene in significant_degs[:current_max_genes]:
            if specific_gene in st.session_state['network'].nodes:
                upstream, downstream = get_specific_gene_relations(st.session_state['network'], specific_gene)
                generate_specific_report(upstream, downstream, specific_gene, save_path=specific_report_path)
                # Download button for specific gene report (based on current network)
                if os.path.exists(specific_report_path):
                    with open(specific_report_path, "rb") as f:
                        st.download_button(
                            label=f"Download {specific_gene} Regulatory Report (TXT)",
                            data=f,
                            file_name=f"{specific_gene}_regulatory_report.txt",
                            mime="text/plain"
                        )
            else:
                st.error(f"Error: Gene `{specific_gene}` is not in the regulatory network of {current_max_genes} genes; please expand the analysis scope")
        else:
            st.error(f"Error: Gene `{specific_gene}` is not among the top {current_max_genes} significant genes; please expand the analysis scope")

    # 9. Additional Download Buttons
    with open("gene_regulatory_edges.txt", "rb") as f:
        st.download_button(
            label="Download All Regulatory Edge Information (TXT)",
            data=f,
            file_name="gene_regulatory_edges.txt",
            mime="text/plain"
        )

    with open("hub_gene_report.txt", "rb") as f:
        st.download_button(
            label="Download Global Hub Gene Report",
            data=f,
            file_name="hub_gene_report.txt",
            mime="text/plain"
        )


if __name__ == "__main__":
    main()