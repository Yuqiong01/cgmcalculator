# cgmcalculator

[![GitHub release](https://img.shields.io/github/v/release/Yuqiong01/cgmcalculator)](https://github.com/Yuqiong01/cgmcalculator/releases)
[![License](https://img.shields.io/github/license/Yuqiong01/cgmcalculator)](LICENSE)

An R package for **Continuous Glucose Monitoring (CGM)** data processing,  
comprehensive metric calculation, AGP-style visualization, and  
function-on-scalar regression using **functional data analysis (FDA)**.

---

## ✨ Key Features

- 📥 Import raw CGM CSV files (one file per subject)
- 🧹 Clean and standardize data via `tidydata()`
- ⏱ Aggregate glucose to an evenly spaced time grid
- 📊 Compute a **comprehensive suite of CGM metrics**
- 🌗 Automatic **daytime / nighttime** stratified metrics
- 📈 Generate **24-hour AGP-style glucose profiles**
- 🧠 Function-on-scalar regression via `refund::pffr`
- 💾 Automatic export of plots, coefficients, and model summaries

---

## 📦 Installation

```r
# install.packages("remotes")
remotes::install_github("Yuqiong01/cgmcalculator")
```

---

## 📁 Required Input Format

All input files must:

- Be **CSV files**
- Contain exactly three columns:

| column name   | description |
|--------------|------------|
| subjectid     | subject identifier |
| timestamp     | date-time of CGM reading |
| sensorglucose | glucose value (mg/dL) |

Use:

```r
tidydata()
```

to standardize raw CGM exports.

---

## 🚀 Workflow

The package uses a **directory-based batch processing pipeline**.

Each function takes a folder containing CGM CSV files as input and writes results to an output folder.

---

### 1️⃣ Tidy raw CGM data

Convert raw CGM exports into standardized files.

```r
library(cgmcalculator)

inputdir  <- "path/to/rawdata_directory"
outputdir <- "path/to/tidydata_directory"

tidydata(
  inputdir  = inputdir,
  outputdir = outputdir,
  skiphours = 0,
  unit = "mg/dL"
)
```

---

### 2️⃣ Calculate CGM metrics

```r
metricsstats(
  inputdir  = "path/to/tidydata_directory",
  outputdir = "path/to/results_directory",
  outputname = "summary_metrics",
  conga_hrs = 1,
  mage_sd = 1,
  calc_day_night = c(6, 22),
  tbr_seq = c(40, 70, 1),
  tar_seq = c(140, 250, 10),
  format = "long",
  unit = "mg/dL"
)
```

Output:

```
summary_metrics.csv
```

---

### 3️⃣ Generate AGP report

```r
agpreport(
  inputdir  = "path/to/tidydata_directory",
  outputdir = "path/to/agp_directory",
  tz = "UTC",
  yaxis = c(0, 500),
  agg_smooth = "loess",
  ptype = "dm",
  mapping_file = "path/to/mapping.csv",
  group_var = "group",
  bin_mins = 5
)
```

Output:

```
AGP_Report_*.pdf
```

---

### 4️⃣ Functional data analysis (FDA)

```r
res <- agpanalyze(
  inputdir       = "path/to/tidydata_directory",
  outputdir      = "path/to/fda_directory",
  covariate_file = "path/to/covariates.csv",
  run_pffr       = TRUE,
  pffr_group     = "micro",
  pffr_covars    = c(
    "age","male","education",
    "dm_duration","egfr","tc","tg"
  ),
  pffr_id_re     = "subjectid"
)
```

Access model results:

```r
res$pffr_plot
summary(res$pffr_fit)
```

---

# 📊 CGM Metrics

The package provides a **comprehensive and extensible CGM metric system**.

---

## 1️⃣ Basic Data & Statistics

- CGM placement and removal time
- Active monitoring time (%)
- Total valid monitoring days

Glucose summary statistics:

- Mean, median
- Minimum, maximum, range
- Standard deviation (SD)
- Interquartile range (IQR)
- Coefficient of variation (CV)
- ...

---

## 2️⃣ Time in Range Metrics (TIR / TAR / TBR)

### ✅ TIR ranges

- 63–140 mg/dL  
- 70–140 mg/dL  
- 70–180 mg/dL  

### 🔺 TAR thresholds

- >140 mg/dL
- >180 mg/dL
- >250 mg/dL
- >400 mg/dL

### 🔻 TBR thresholds

- <70 mg/dL
- <63 mg/dL
- <54 mg/dL
- <40 mg/dL

For each range:

- Time (%)
- AUC
- Event count
- Event duration
- Mean glucose in range

### ⚙ Custom scanning

User-defined:

```r
tar_seq
tbr_seq
```

---

## 3️⃣ Glycemic Variability & Fluctuation Metrics

Includes:

- MAGE
- MODD
- CONGA
- LBGI
- HBGI
- M-value
- J-index
- ADRR
- MAG
- GVP
- ...

---

## 4️⃣ Autocorrelation & Rate of Change

Rate-of-change (ROC) distribution:

- Mean
- Variance
- SD
- Percentiles

(ACF-based metrics supported depending on analysis pipeline)

---

## 5️⃣ GRADE & Glycemic Risk

- GRADE score
- Proportion of:
  - Euglycemia
  - Hypoglycemia
  - Hyperglycemia
- Glycemic Risk Index (GRI)
- Continuous Glucose Monitoring Index (COGI)

---

## 6️⃣ Daytime / Nighttime Metrics 🌗

All major metrics are also computed for:

- Daytime → `*_day`
- Nighttime → `*_night`

This enables:

- Diurnal pattern analysis
- Context-specific glycemic control assessment

---

# 📈 Time-of-Day Features

- Hourly glucose curves
- Smoothed 24-hour glucose profiles
- Functional data objects for FDA

---

# 🧠 Functional Data Analysis

Compatible with:

```r
refund::pffr
```

Supports:

- Function-on-scalar regression
- Subject-specific random effects
- Coefficient curve export
- Confidence intervals

---

# 📤 Typical Outputs

- CGM metric tables
- Day/night stratified metrics
- AGP plots
- FDA coefficient curves
- Model summaries

---

# 🎯 Package Scope

`cgmcalculator` is designed for:

- Cohort-scale CGM studies
- Automated clinical CGM metric pipelines
- AGP-based visualization
- FDA-based inference

---

# 📜 License

MIT License

---

# 👤 Author

Yuqiong Li.

---

# ⭐ Citation

A formal citation will be provided after the associated methods paper is published.