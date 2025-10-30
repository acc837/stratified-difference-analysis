# Gene Expression Analysis Workflow

An interactive gene expression analysis platform based on Streamlit, integrating complete workflow of gene stratification, differential expression analysis, and regulatory network construction.

## Research Team Information

### Affiliated Institutions
1. **Department of Otolaryngology, Head and Neck Surgery, Beijing Tongren Hospital, Capital Medical University, Beijing, China**
2. **Department of Biochemistry and Molecular Biology, Capital Medical University, Beijing, China**

### Team Members
- **Developer**: Tian Ye
- **Project Leaders**: Lu Kong, Xiaohong Chen
- **Project members**: Guoliang Yang, Xudong Wang, Tingyao Ma, JiaXin Chen, Fang Nan, Qian Chen

### Research Background
This tool was developed based on the practical needs in our research on molecular mechanisms of otolaryngology-head and neck surgery diseases, aiming to provide an efficient and accurate gene expression analysis platform for researchers in related fields.

## Features

### Core Analysis Modules
- **Gene Expression Stratification** - Automatically classify genes into high, medium, and low expression layers based on expression levels
- **Differential Expression Analysis** - Professional-level differential analysis using pyDESeq2
- **Regulatory Network Construction** - Infer gene regulatory relationships based on GENIE3 algorithm
- **Hub Gene Identification** - Automatically discover global key regulatory genes
- **Specific Gene Analysis** - In-depth study of upstream regulation and downstream targets of specific gene

### Visualization & Output
- Interactive gene stratification distribution plots
- Differential expression indicator distribution visualization
- Regulatory network edge information export
- Multi-format report downloads (CSV, TXT, PNG)

## Environment Requirements

### System Requirements
- Python 3.8+
- Windows/Linux/macOS
- Minimum 8GB RAM (recommended 16GB+)
- UTF-8 encoding supported environment

### Dependencies
See `requirements.txt` for details, main dependencies include:
- streamlit==1.28.0
- pandas==1.5.3
- numpy==1.24.3
- scikit-learn==1.3.0
- pydeseq2==0.4.3
- networkx==2.8.8
- matplotlib==3.7.2
- seaborn==0.12.2

## Installation Steps

### 1. Clone or Download Project Files
Ensure the following files are included:
- `stratification.py` - Main program file
- `GENIE3.py` - Regulatory network algorithm
- `requirements.txt` - Dependency package list

### 2. Create Python Environment (Recommended)
```bash
# Create conda environment
conda create -n gene_analysis python=3.10
conda activate gene_analysis

# Or create venv environment
python -m venv gene_analysis
source gene_analysis/bin/activate  # Linux/macOS
gene_analysis\Scripts\activate    # Windows
