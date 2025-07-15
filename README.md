# NEN-methylation

DNA methylation patterns facilitate tracing the origin of neuroendocrine neoplasms (NEN) through a reproducible end-to-end pipeline that covers:

**1. CpG site selection**  
**2. Copy-number & deconvolution** (MeDeCom)  
**3. Latent methylation component** (LMC) analysis & plotting  
**4. Feature engineering** & **residualization** (R)  
**5. Random Forest** tissue-of-origin classifier (Python)  
**6. XGBoost** feature extraction & RF re-training (Python)  
**7. Validation** on independent liver-metastasis cohort  

---

## 🚀 Features

* **R-based preprocessing** of Illumina EPIC methylation arrays (probe filtering, normalization, CNA calling)
* **MeDeCom** deconvolution to identify latent methylation components
* **Automated heatmap** and **correlation** plotting of LMC proportions
* **Machine-learning classifiers** (RF & XGB) trained on discovery cohort and validated on liver‐metastasis samples
* **Configurable**: change parameters in one place, run each step as an R script
* **Reproducible**: conda environment spec, clear directory structure, version control

---

## 📦 Repository Structure

```
NEN-methylation/
├── .gitignore                 # ignore R/Python/Jupyter/data artifacts
├── LICENSE
├── environment.yml            # conda environment spec
├── README.md                  # this file
├── scripts/                   # R scripts, run in order
│   ├── 1_site_selection_update.R
│   ├── 2_medecom_validation.R
│   ├── 3_lmc_analysis.R
│   ├── 4_rf_model.R
│   ├── 5_ml_5k_validation.R
│   └── 6_xgb_model.R
├── notebooks/                 # Jupyter notebooks for exploration & figures
│   ├── xgb_training_RF.ipynb
│   ├── ML_val_lmc3.ipynb
│   ├── confusion_matrices_merged.ipynb
└────────────────────────
```

---

## 🔧 Prerequisites

* [Conda](https://docs.conda.io/) (Miniconda or Anaconda)
* R (≥ 4.1) with Bioconductor support
* Python (3.8)

---

## 🛠️ Setup

1. **Clone the repo**

   ```bash
   git clone https://github.com/Danialbt/NEN-methylation.git
   cd NEN-methylation
   ```

2. **Create and activate the conda environment**

   ```bash
   conda env create -f environment.yml
   conda activate nen-pipeline
   ```

---

## ▶️ Running the Pipeline

Each step is an independent R script that reads from `data/` or `data/processed/` and writes outputs to `data/processed/`, `results/`, or `notebooks/`. Run them in order:

**1. CpG site selection**

   * Differential‐methylation site extraction across anatomical sites
   * Outputs: `selected12k_sites_update` RDS in `data/processed/`

   ```bash
   Rscript scripts/1_site_selection_update.R
   ```

**2. MeDeCom deconvolution**

   * Runs MeDeCom on discovery + validation CpG sets

   ```bash
   Rscript scripts/2_medecom_validation.R
   ```

**3. LMC analysis & plotting**

   * Generates heatmaps & boxplots for LMC proportions (with/without LMC2 or LMC3)
   * Correlates LMCs with LUMP score
   * Saves PDF figures to `results/lmc/`

   ```bash
   Rscript scripts/3_lmc_analysis.R
   ```

**4. Feature engineering & residualization**

   * 1. Load MeDeCom set (md.res_...rds) and get LMC proportions
   * 2. From RnBeads object, select 5 000 most variable CpGs
   * 3. Intersect with validation CpGs
   * 4. Pick top 5 000 by variance → meth.data_val_top5000
   * 5. For each CpG, fit β ~ LMC3:
        * If p < 0.01, replace β with residuals
        * Else center β at zero
   * 6. Save:
        * NEN_ml_5k_org.csv (raw 5 000 features)
        * NEN_ml_5k_res_test_update.csv (residualized) in data/
    

   ```bash
   Rscript scripts/4_rf_model.R
   ```

**5. Random Forest classifier**

   * Execute notebooks/ML_val_lmc3.ipynb, or headless:
   * Load NEN_ml_5k_res_test_update.csv + annotations
   * Exclude known metastases & CUP
   * Train RandomForestClassifier with:
        * n_estimators=1500
        * max_features='sqrt'
        * max_depth=12
        * oob_score=True
        * random_state=3
   * Stratified 3-fold CV, compute OOB score, ROC AUC per class
   * Save ROC plots, feature importance, confusion matrices to results/rf/
  

   ```bash
   jupyter nbconvert --execute notebooks/ML_val_lmc3.ipynb \
     --to notebook --output results/xgb/ML_val_lmc3_output.ipynb
   ```

**6. XGBoost feature extraction & RF re-training**

   * Open & execute notebooks/xgb_training_RF.ipynb, or run headless:
   * Train XGBoost multiclass softprob (nrounds = 800)
   * Extract top 404 CpGs by importance
   * Subset methylation matrix → train new RF with same hyperparameters
   * Stratified 3-fold CV, ROC, confusion matrix
   * Save:
        * NEN_ml_xgboosting_update2.csv (full-feature table)
        * NEN_ml_xgb_val_update.csv (validation predictions)
        * Plots to results/xgb/

          
   ```bash
   jupyter nbconvert --execute notebooks/xgb_training_RF.ipynb \
     --to notebook --output results/rf/xgb_training_RF_output.ipynb
   ```

**7. Validation & final figures**  

Open notebooks/confusion_matrices_merged.ipynb to compare RF vs XGB, inspect class-level performance, and produce final summary figures.

---

## 📝 License

This project is licensed under the [MIT License](LICENSE).

---

## 🤝 Contributing

Feel free to submit issues or pull requests. For major changes, please open an issue first to discuss what you’d like to change.

---

## 📧 Contact



