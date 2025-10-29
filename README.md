# 基因表达分析工作流 (Gene Expression Analysis Workflow)

一个基于Streamlit的交互式基因表达分析平台，集成了基因分层、差异表达分析和调控网络构建的完整工作流。

## 功能特性

### 核心分析模块
- **基因表达分层** - 基于表达量自动划分高、中、低表达基因
- **差异表达分析** - 使用pyDESeq2进行专业级差异分析
- **调控网络构建** - 基于GENIE3算法推断基因调控关系
- **枢纽基因识别** - 自动发现全局关键调控基因
- **特定基因分析** - 深入研究任意基因的上游调控和下游靶标

### 可视化与输出
- 交互式基因分层分布图
- 差异表达指标分布可视化
- 调控网络边信息导出
- 多格式报告下载(CSV, TXT, PNG)

## 环境要求

### 系统要求
- Python 3.8+
- Windows/Linux/macOS
- 至少8GB内存（推荐16GB+）
- 支持UTF-8编码的环境

### 依赖包
详见 `requirements.txt`，主要依赖包括：
- streamlit==1.28.0
- pandas==1.5.3
- numpy==1.24.3
- scikit-learn==1.3.0
- pydeseq2==0.4.3
- networkx==2.8.8
- matplotlib==3.7.2
- seaborn==0.12.2

## 安装步骤

### 1. 克隆或下载项目文件
确保包含以下文件：
- `stratification.py` - 主程序文件
- `GENIE3.py` - 调控网络算法
- `requirements.txt` - 依赖包列表

### 2. 创建Python环境（推荐）
```bash
# 创建conda环境
conda create -n gene_analysis python=3.8
conda activate gene_analysis

# 或创建venv环境
python -m venv gene_analysis
source gene_analysis/bin/activate  # Linux/macOS
gene_analysis\Scripts\activate    # Windows
