# Biostatistical Significance Analysis & Clinical Outcome Modeling on Regional NCD Datasets

## 📌 Project Overview
This repository contains the end-to-end quantitative biostatistical pipeline engineered as part of the technical selection process for the **Research Engineer (Data Science)** position at the **Institute for Advanced Research (IAR) / International Research Center (IR_IRC), United International University (UIU)**.

The primary objective of this analysis is to evaluate an empirical public health dataset ($N = 29,999$ unique clinical records) across rural geographic cohorts to identify, isolate, and validate key demographic, socioeconomic, and physiological risk factors underlying major Non-Communicable Diseases (NCDs): **Stroke (`had_stroke`), Cardiovascular Disease (`has_cardiovascular_disease`), Diabetes (`diabetic`), and Hypertension (`profile_hypertensive`)**.

---

## 🧬 Data Architecture & Feature Engineering

### 1. Raw Dataset Profiles
The base cohort tracking system spans 31 distinct variables capturing extensive participant metrics:
* **Demographics:** `age`, `gender` (Female: 23,236; Male: 6,763 — indicating a highly skewed gender distribution).
* **Socioeconomic Metrics:** `total_income`, `is_poor` (homogenous across sample array; filtered out).
* **Clinical Telemetry:** `SYSTOLIC`, `DIASTOLIC`, `PULSE_RATE`, `HEIGHT`, `WEIGHT`, `BMI`, `SUGAR`, `SPO2`, `MUAC`.
* **Diagnostic Categorization:** `RESULT_STAT_BP`, `RESULT_STAT_BMI`, `RESULT_STAT_SUGAR`, `RESULT_STAT_PR`, `RESULT_STAT_SPO2`, `TAG_NAME`.
* **Target Vector Labels:** `had_stroke`, `has_cardiovascular_disease`, `diabetic`, `profile_hypertensive`.

### 2. Missing Value & Data Manipulation Strategy
A robust, adaptive data manipulation pipeline was implemented to handle systemic data sparseness:
* **High-Missingness Filtering ($>85\%$ Nulls):** The continuous variables `HEIGHT`, `WEIGHT`, `BMI`, `SUGAR`, and `SPO2` contained massive data gaps ($25,654$ to $28,871$ empty entries). To preserve this structural missingness without corrupting the dataset, **Missingness Indicator Flags** were generated (e.g., `HEIGHT_measured`), mapping availability to discrete binary indicators ($0$ or $1$).
* **High-Missingness Categorical Alignment:** Empty cells within related diagnostic classification strings (`RESULT_STAT_BMI`, `TAG_NAME`, `RESULT_STAT_SUGAR`, `RESULT_STAT_SPO2`) were filled with the explicit category `'Not Measured'`.
* **Low-Missingness Stratified Imputation ($<10\%$ Nulls):** Essential vital metrics (`SYSTOLIC`, `DIASTOLIC`, and `PULSE_RATE`) had sparse missing data (~$8\%$). Missing values were filled via **Stratified Group-Median Imputation**, computed dynamically across corresponding matching combinations of `gender` and structural `age_group` brackets (Binned: `<18`, `18-35`, `36-50`, `51-65`, `65+`). Residual missing cells defaulted to global column medians.
* **Low-Missingness Categorical Mapping:** Missing values in blood pressure and pulse rate statuses (`RESULT_STAT_BP`, `RESULT_STAT_PR`) were mapped to the descriptive fallback category `'Unknown'`.

---

## 🛠️ Statistical Methodology & Mathematical Rigor

### 1. Categorical Feature Dependency Mapping ($\chi^2$ Test)
To map dependencies among categorical risk factors, a Chi-Square Test of Independence was used to evaluate observed vs. expected cell frequencies across constructed contingency matrices.
* **Mathematical Equation:**
    $$\chi^2 = \sum_{i=1}^{r} \sum_{j=1}^{c} \frac{(O_{ij} - E_{ij})^2}{E_{ij}}, \quad E_{ij} = \frac{R_i \times C_j}{N}, \quad df = (r-1)(c-1)$$
* **Hypothesis Formulation:**
    * $H_0$: The independent predictor feature and the clinical target outcome are independent.
    * $H_1$: The independent predictor feature and the clinical target outcome are statistically dependent.

### 2. Continuous Metric Variational Testing (Welch's T-Test)
To identify shifts in continuous biological features, independent group differences were evaluated. Welch's $t$-test was selected as the primary parametric engine to handle unequal group variance ($\sigma_0^2 \neq \sigma_1^2$) and severe sample class imbalances ($n_0 \neq n_1$). It remains highly robust under the Central Limit Theorem ($CLT$) given $N = 29,999$.
* **Mathematical Equation:**
    $$t = \frac{\bar{X}_1 - \bar{X}_2}{\sqrt{\frac{s_1^2}{n_1} + \frac{s_2^2}{n_2}}}, \quad df = \frac{\left(\frac{s_1^2}{n_1} + \frac{s_2^2}{n_2}\right)^2}{\frac{\left(\frac{s_1^2}{n_1}\right)^2}{n_1 - 1} + \frac{\left(\frac{s_2^2}{n_2}\right)^2}{n_2 - 1}}$$
* **Hypothesis Formulation:**
    * $H_0: \mu_1 = \mu_2$ (The population means of healthy and affected diagnostic cohorts are equal).
    * $H_1: \mu_1 \neq \mu_2$ (The population means of both cohorts are significantly different).

### 3. Continuous Distribution Verification (Mann-Whitney U Safeguard)
Formal **Kolmogorov-Smirnov (K-S)** screening significantly rejected absolute data normality across vital metrics ($p < 0.05$ across all tests; `SYSTOLIC` $p = 4.4446 \times 10^{-245}$, `age` $p = 4.9029 \times 10^{-115}$), capturing notable right-skewness and rounding bias. To ensure structural protection against outlier inflation, a non-parametric **Mann-Whitney U Test** was executed concurrently to cross-validate continuous variables based on their ordinal ranks.
* **Mathematical Equation:**
    $$U_i = n_1 n_2 + \frac{n_i(n_i + 1)}{2} - R_i, \quad U = \min(U_1, U_2)$$
* **Hypothesis Formulation:**
    * $H_0: P(X > Y) = P(Y > X)$ (The distribution profiles of both groups are stochastically equal).
    * $H_1: P(X > Y) \neq P(Y > X)$ (One target cohort systematically exhibits greater values than the other).

---

## 📈 Empirical Insights & Core Analytical Outcomes

### 1. Categorical Contingency Matrix Insights ($\alpha = 0.05$)
* **Geographic Regionalism (`union_name`):** Showed extreme statistical significance ($p < 0.001$) across *all four target clinical outcomes*. This confirms that chronic disease distributions are non-randomly clustered within specific regional Union Parishads (e.g., `union_name` vs `diabetic` $p = 2.5741 \times 10^{-129}$).
* **Socioeconomic Class (`total_income`):** Shares a highly significant relationship with underlying metabolic disorders, specifically `diabetic` ($p = 2.2885 \times 10^{-08}$) and `profile_hypertensive` ($p = 1.8919 \times 10^{-02}$), while remaining independent ($p > 0.05$) of acute, sudden stroke occurrences (`had_stroke` $p = 0.8253$).
* **Biological Sex (`gender`):** Demonstrated a statistically significant association exclusively with `had_stroke` incidence rates ($p = 0.0312$), but was non-significant across cardiovascular and diabetes categories ($p > 0.05$).

### 2. Continuous Cohort Variation Analysis
* **Age Profile (`age`):** Confirmed as the most dominant covariate across all pathology profiles ($p < 10^{-5}$); mean age scales up significantly within impacted diagnostic groups.
* **Hemodynamic Markers (`SYSTOLIC` / `DIASTOLIC`):** Showed clear significance ($p < 10^{-20}$) when analyzing risk variance for both `diabetic` and `profile_hypertensive` cohorts.
* **Resting Pulse (`PULSE_RATE`):** Exhibited a localized, statistically isolated significance threshold ($p = 8.7403 \times 10^{-15}$) uniquely linked to the `profile_hypertensive` target configuration.

### 3. Extracted Feature-Dependency Maps for Modeling
Based on the convergence of the statistical tests, the features below have been isolated as highly relevant for subsequent model development:
* **Stroke Matrix (`had_stroke`):** `gender`, `age`, `SYSTOLIC`, `DIASTOLIC`, `union_name`, `has_cardiovascular_disease`, `diabetic`, `profile_hypertensive`.
* **Cardiovascular Disease Matrix (`has_cardiovascular_disease`):** `age`, `SYSTOLIC`, `DIASTOLIC`, `had_stroke`, `union_name`, `diabetic`, `profile_hypertensive`.
* **Diabetes Prediction Matrix (`diabetic`):** `age`, `total_income`, `SYSTOLIC`, `DIASTOLIC`, `union_name`, `has_cardiovascular_disease`, `profile_hypertensive`.
* **Hypertension Diagnostic Matrix (`profile_hypertensive`):** `age`, `total_income`, `SYSTOLIC`, `DIASTOLIC`, `PULSE_RATE`, `union_name`, `has_cardiovascular_disease`, `diabetic`.

---

## 🧑‍🔬 Author & Candidate Profile
* **Name:** Abu Omayed
* **Position Applied For:** Research Engineer (Data Science)
* **Target Institution:** International Research Center (IR_IRC), United International University (UIU)
* **ORCID:** [0009-0006-1607-7979](https://orcid.org/0009-0006-1607-7979)
* **LinkedIn:** [linkedin.com/in/abu-omayed](https://www.linkedin.com/in/abu-omayed/)
* **GitHub:** [github.com/AbuOmayed](https://github.com/AbuOmayed)
* **Portfolio Website:** [abuomayed.netlify.app](https://abuomayed.netlify.app/)
