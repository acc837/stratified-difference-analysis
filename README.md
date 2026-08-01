# Target-gene Stratified DESeq2 + GENIE3 Regulatory Network Workflow
This is a corrected and robust implementation for **target-gene stratified differential analysis and GENIE3-based gene regulatory network inference**, optimized for RNA-seq cohort analysis. The pipeline is inspired by and structurally improved from stratified difference analysis workflow.

---

## 🏛 Affiliated Institutions
- **Department of Otolaryngology, Head and Neck Surgery, Beijing Tongren Hospital, Capital Medical University, Beijing, China**
- **Department of Biochemistry and Molecular Biology, Capital Medical University, Beijing, China**

## 👥 Team Members
- **Developer**: Tian Ye
- **Project Leaders**: Lu Kong, Xiaohong Chen
- **Project Members**: Guoliang Yang, Xudong Wang, Tingyao Ma, Jiaxin Chen, Fang Nan, Qian Chen

## 📚 Research Background
This tool was developed based on the practical needs in our research on molecular mechanisms of otolaryngology-head and neck surgery diseases, aiming to provide an efficient and accurate gene expression analysis platform for researchers in related fields.

---

## 🔬 Intended Workflow
This web application performs a complete stratified transcriptomic analysis pipeline:

1. Upload a **normalized expression matrix** (TPM / FPKM / log-normalized values).
2. Upload the matched **raw integer count matrix** for differential analysis.
3. Define a **target gene A** for sample stratification.
4. Rank all samples by the normalized expression level of target gene A.
5. Group samples into **A-high** and **A-low** subgroups (The default setting uses the top and bottom 25% of samples, while the stratification proportion should be adjusted according to the cohort size and the expression distribution of target gene A).
6. Run **PyDESeq2** differential analysis: `A-high vs A-low`.
7. Filter statistically significant DEGs by `padj` and `|log2FC|` thresholds.
8. Construct a multigene GENIE3 regulatory network using target gene A and significant DEGs, with optional inclusion of expressed known transcription factors and signaling regulators related to the pathway of interest.
9. Select a partner gene B and extract **bidirectional regulatory weights** (`A→B` and `B→A`).

---

## 💡 Multigene Network Design Principle
**Two-gene-only GENIE3 modeling is intentionally avoided.**

A pairwise A/B-only network produces meaningless feature importance, because there is no predictive competition between features. Without multigene context, regulatory weights cannot represent true biological predictability.

This pipeline builds a **context-rich multigene network** comprising target gene A and significant DEGs. Depending on the biological question and dataset, expressed known transcription factors and pathway-related signaling regulators may also be included as candidate genes. The A–B regulatory relationship is then extracted from the resulting multigene network, allowing each edge weight to be evaluated within a broader predictive context.

---

## 📄 Input File Format Requirements
Both CSV matrices follow identical structure: **genes in rows, samples in columns**. First column = gene ID.

```
Gene,S1,S2,S3,S4,...
GENE1,12.3,8.1,9.4,15.6,...
GENE2,0.0,2.5,1.1,0.8,...
```

### Normalized Expression Matrix
- Accepts: TPM, FPKM, or pre-log-transformed values
- Duplicate genes: averaged

### Raw Count Matrix
- Requires: **non-negative integer raw reads**
- Duplicate genes: summed
- Used exclusively for DESeq2 differential testing

### Matching Rules
- Sample IDs must be identical between two matrices
- Gene annotation systems should be consistent

---

## ⚙️ Installation & Running
Python 3.10 – 3.12 recommended

### 1. Create virtual environment
```bash
python -m venv gene_network_env
```

### 2. Activate environment
**Windows PowerShell**
```bash
gene_network_env\Scripts\Activate.ps1
```

**Linux / macOS**
```bash
source gene_network_env/bin/activate
```

### 3. Install dependencies
```bash
pip install streamlit pandas numpy scipy scikit-learn pydeseq2 matplotlib networkx seaborn
```

### 4. Launch web app
```bash
streamlit run main.py
```

---

## 🧭 Key Methodological Choices
### Target-Gene Stratification
- The default setting defines the top 25% of samples as A-high and the bottom 25% as A-low
- The grouping proportion should be adjusted according to sample size and the expression distribution of target gene A
### Differential Analysis (DESeq2)
- Only **raw integer counts** are used for DE modeling
- Normalized TPM/FPKM data are never used for statistics

### GENIE3 Network Construction
- Network inference uses **normalized expression matrix**
- Auto `log2(x+1)` transformation for unlogged TPM/FPKM
- Turn off transform if input data are already log-normalized
- Candidate genes include target gene A and significant DEGs; expressed known transcription factors and pathway-related signaling regulators may be additionally included according to the biological context

### Sample Modes for Network Training
- **Default mode**: Only A-high / A-low extreme samples (strict stratified workflow)
- **Sensitivity mode**: All shared samples, reducing selection bias and circularity

---

## 📌 Result Interpretation
A GENIE3 directed edge represents **predictive gene expression dependency** from tree-ensemble modeling.

It does **not** directly prove:
- direct molecular binding
- definite causality
- transcriptional activation or repression

**Spearman correlation** is provided separately to evaluate co-expression direction and statistical association. Experimental validation is required for causal claims.

---

## 📁 Pipeline Outputs
- Sample stratification table (group / expression value)
- Full PyDESeq2 differential result table
- Filtered significant DEG list
- Complete GENIE3 edge table with weight/rank/percentile
- GENIE3 gene regulatory weight matrix (VIM matrix)
- Bidirectional A↔B regulatory comparison results
- Top upstream / downstream target genes of the target gene
- Visualization plots: stratification barplot, co-expression scatter plot, regulatory network graph

---

## 🗂 Project Structure
```
.
├── main.py              # Streamlit web UI & full pipeline controller
├── analysis_core.py    # Core analysis functions (stratification, DESeq2, self-implemented GENIE3)
├── plotting.py         # All visualization functions
└── README.md           # Project documentation

