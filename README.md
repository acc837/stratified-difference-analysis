Target-gene stratified DESeq2 + GENIE3
This is a corrected implementation inspired by the workflow style of acc837/stratified-difference-analysis.

Intended workflow
Upload a normalized expression matrix, such as TPM or FPKM.
Upload the matching raw integer count matrix.
Specify target gene A.
Rank samples by normalized expression of A.
Define A-high and A-low samples, defaulting to the upper and lower 25%.
Select those samples in the raw count matrix and run PyDESeq2: A-high versus A-low.
Select significant DEGs using padj and |log2FC| thresholds.
Build a GENIE3-style multigene network using normalized expression of A + significant genes.
Select a significant gene B and extract both A → B and B → A weights.
Why the network is multigene
Do not fit GENIE3 using only A and B. With only one predictor for each target, feature importance has no meaningful competition and the resulting weights are not informative. This app first builds a network containing A, B, and other significant genes, then extracts the A–B pair.

Input format
Both CSV files must have genes in rows and samples in columns. The first column contains gene IDs.

Gene,S1,S2,S3,...
GENE1,12.3,8.1,9.4,...
GENE2,0.0,2.5,1.1,...
The normalized matrix may contain TPM, FPKM, or already log-transformed values.
The raw-count matrix must contain non-negative integers.
Sample IDs must match between files.
Duplicate normalized-expression gene IDs are averaged.
Duplicate raw-count gene IDs are summed.
Installation
Python 3.10–3.12 is recommended.

python -m venv gene_network_env

# Windows PowerShell
gene_network_env\Scripts\Activate.ps1

# Linux/macOS
source gene_network_env/bin/activate

pip install -r requirements.txt
streamlit run app.py
Important methodological choices
DESeq2 input
Only raw integer counts are passed to PyDESeq2. TPM/FPKM values are never used for differential-expression modeling.

GENIE3 input
GENIE3 uses the normalized expression matrix. For unlogged TPM/FPKM, the default is log2(x + 1). Disable that option when the uploaded matrix is already log-transformed.

Samples used for GENIE3
The default reproduces the requested workflow and uses only the A-high/A-low samples. The app also provides an all shared samples sensitivity mode, which retains continuous variation and helps assess circularity caused by discovering DEGs and estimating network edges in the same extreme samples.

Interpretation
A GENIE3 edge means that one gene's expression has predictive importance for another target gene within a tree-ensemble model. It does not prove:

direct molecular regulation;
causal direction;
transcription-factor binding;
activation versus repression.
The app reports Spearman correlation separately as an association/sign diagnostic. Experimental or independent-cohort validation is still required.

Main outputs
sample stratification table;
full PyDESeq2 result table;
significant DEG list;
complete GENIE3 edge list;
GENIE3 weighted adjacency matrix;
A→B and B→A edge weights and ranks;
top predicted incoming and outgoing edges for A.
