# DESeq2-ANALYSIS_TRIALDATA
Bioinformatics project using DESeq2 for RNA-seq differential expression. Includes CSV results, Shiny app for interactive DNA sequencing data exploration, PDF/HTML reports, and documentation. Organized for reproducibility, visualization, and collaborative analysis.
## 📖 Overview
This project contains a differential gene expression analysis performed using **DESeq2** in R.  
It includes:
- R scripts for data simulation and analysis
- CSV results of significant genes
- A Shiny app for interactive exploration
- A PDF report summarizing the workflow and findings

---

## 📂 Project Structure
DESeq2-Differential-Expression/ │ ├── scripts/ │ └── analysis.R # R script for DESeq2 analysis │ ├── results/ │ └── DESeq2_CSVResults/ # CSV output files │ ├── shiny_app/ │ └── DNA SEQ SHINY APP/ # Shiny app for interactive visualization │ ├── figures/ │ └── volcano.png # Example plots (Volcano, Heatmap, PCA) │ └── heatmap.png │ ├── DESEQ_REPORT.pdf # Final report ├── README.md # Project documentation ├── LICENSE # License file └── .gitignore # Ignore unnecessary files

---

## ⚙️ Requirements
- **R (≥ 4.0)**
- Packages:
  - `DESeq2`
  - `ggplot2`
  - `pheatmap`
  - `dplyr`
  - `shiny`

Install required packages:
```R
install.packages(c("ggplot2", "pheatmap", "dplyr"))
BiocManager::install("DESeq2")
________________________________________
🚀 Usage
1. Run DESeq2 Analysis
Execute the R script:
source("DESeq2_TRIALDATA.R")
This will generate results in results/DESeq2_CSVResults/.
2. Explore with Shiny App
Launch the interactive app:
shiny::runApp("DNA SEQ SHINY APP")
3. View Report
Open DESEQ_REPORT.pdf for a detailed summary of the analysis.
________________________________________
📊 Outputs
CSV files: Differentially expressed genes with fold changes and p-values.
Figures: Volcano plot, PCA plot, heatmap of top DEGs.
Shiny App: Interactive visualization of gene expression data.
________________________________________
📜 License
This project is licensed under the MIT License.
Feel free to use, modify, and share with attribution.
________________________________________
👨‍🔬 Author
PH. Fouad Mohsen
Bioinformatics Analysis TRIAL DATA – January 2026

---
