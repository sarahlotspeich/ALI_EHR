# Overcoming data challenges with enriched validation and targeted sampling to measure whole-person health in electronic health records

This repository contains `R` code and simulation data to reproduce results from the manuscript by [Lotspeich et al. (2025+)](https://arxiv.org/abs/2502.05380). *Please note that due to the sensitive information of our application to the learning health system, we cannot share patient data to replicate the real-data analysis.*

This code relies on two `R` packages developed to accompany the paper, among other standard packages available from CRAN.  
  
  1.  `logiSieve`: implements the sieve maximum likelihood estimator (SMLE) for logistic regression with covariate measurement error from the paper. The package can be found in its own repo [here](https://github.com/sarahlotspeich/logiSieve) and installed in `R` as follows:

``` r
# Run once
# install.packages("devtools")
devtools::install_github("sarahlotspeich/logiSieve", ref = "main")
```

  2.  `auditDesignR`: implements various sampling strategies to select the validation subset, including the residual sampling design from the paper. The package can be found in its own repo [here](https://github.com/sarahlotspeich/auditDesignR) and installed in `R` as follows:

``` r
devtools::install_github("sarahlotspeich/auditDesignR", ref = "main")
```

## Graphical Abstract 

## Figures 

### Made in Canva - No data or code

**Figure 1.** Diagram of potential errors and missingness in a patient's hemoglobin A1c (HbA1c) measurement in the electronic health record (EHR) before and after validation (chart review). In the top scenario, the patient had HbA1c measured but with error, and the correct value was recovered through validation. In the bottom scenario, the patient did not have HbA1c measured, but auxiliary information was located through validation to replace the missing component.

  - [Figure](figures/Fig1_ALI_EHR_Flowchart.png)

**Figure 2.** After extracting the electronic health records (EHR) data, partial validation studies break down into three key steps: analysis, design, and protocol. These steps are intertwined, with design choices depending on the planned analysis, and each step provides opportunities to gain information. For added flexibility, a multi-wave validation framework allows us to revisit these key steps and make modifications as information accumulates.

  - [Figure](figures/Fig2_Full_ALI_EHR_Flowchart_No_Background.png)

**Figure 3.** Diagram of partial validation studies with $N$ patients extracted from the electronic health records (EHR) and $n$ patients chosen for validation via chart review. The least informative strategy to select the validation study is **A)** simple random sampling. More informative targeted strategies include **B)** case-control sampling, **C)** balanced case-control sampling using discretized ALI, **D)** optimal stratified sampling using discretized ALI, and **E)** residual sampling based on the naive model. Within strata, the darker shaded areas indicate which patients would be selected.

  -  [Figure](figures/Fig3_Combined_ALI_Designs.png)

**Figure 4.** Flow chart of instructions given to auditors for the chart reviews (Steps 1-3) and the corresponding auditor finding assigned by the analysts (Step 4). Column names refer to those in Supplemental Figure S2.

  -  [Figure](figures/Fig4_Protocol_Flowchart.png)

**Figure S1.** Ten component stressors of the allostatic load index were taken across three body systems (cardiovascular, metabolic, and inflammation), discretized at clinically meaningful thresholds, and combined to create a whole-person health score.

  -  [Top Figure](figures/FigS1A_CalcALI_Original.png)
  -  [Bottom Figure](figures/FigS1B_CalcALI_Disc.png)

**Figure S2.** Illustrative example of how the original extracted EHR data for patients' vitals were transformed from a wide format (with one row per patient encounter and columns for different variables) to a long format (with one row per patient encounter per variable), merged with the audit roadmap, and augmented with columns where auditors entered their findings. The same was done with patients' lab data.

  -  [Figure](figures/FigS2_Example_Transform_Data.png)

### Made in R - Data (where shareable) and code

**Figure 5.** Auditors' findings from the validation of **A)** $52$ patients' data in the Pilot + Wave I, **B)** $48$ patients' data in  Wave II, and **C)** $100$ patients' data in all waves combined. These findings refer to the original numeric measurements (before discretizing them into allostatic load index components), and there could be multiple per patient per variable. Auditor findings were mutually exclusive and specific to whether the measurements were originally *non-missing* ("extracted value correct," "extracted value incorrect," or "extracted value not found") or *missing* ("no auxiliary information found" or "auxiliary information found") in the extracted EHR data.

  -  [Figure](figures/Fig5_Findings_Bargraph.png)
  -  [Script (Make Figure)](figure-scripts/Fig5_Findings_Bargraph.R)

**Figure 6.** Comparison of the discretized allostatic load index components from the unvalidated EHR data to **A)** $52$ patients' data in the Pilot + Wave I, **B)** $48$ patients' data in  Wave II, and **C)** $100$ patients' data in all waves combined.

  -  [Figure](figures/Fig6_Heatmap.png)
  -  [Script (Make Figure)](figure-scripts/Fig6_Heatmap.R)

**Figure 7.** Among patients with chart reviews (triangles), validated allostatic load index (ALI) could differ noticeably from the unvalidated version using the EHR data. For patients without chart reviews (squares), validated ALI could be robustly predicted using the estimated exposure error mechanism. The dashed line indicates the line of equality (i.e., where unvalidated and validated ALI are the same). The solid line is a loess smoother.

  -  [Figure](figures/Fig7_Scatterplot_PredictedALI.png)
  -  [Script (Make Figure)](figure-scripts/Fig7_Scatterplot_PredictedALI.R)

**Figure S3.** Based on the original extracted EHR data (before validation), the distribution of the error-prone allostatic load index (ALI) was right-skewed with a median of $0.33$ (denoted by the vertical dashed line) and an interquartile range of $(0.17, 0.50)$.

  -  [Figure](figures/FigS3_Histogram_Unvalidated_ALI.png)
  -  [Script (Make Figure)](figure-scripts/FigS3_Histogram_Unvalidated_ALI.R)

**Figure S4.** Missingness in the discretized allostatic load index (ALI) components in the extracted EHR data (before validation), after the Pilot + Wave I validation, and after all waves of validation. In **A)**, missingness is broken down by the ALI component, and in **B)**, it is broken down by the number of missing ALI components per patient.

  -  [Top Figure](figures/FigS4A_missing_by_component.png)
  -  [Script (Make Top Figure)](figure-scripts/FigS4A_missing_by_component.R)
  -  [Bottom Figure](figures/FigS4B_missing_by_patient.png)
  -  [Script (Make Bottom Figure)](figure-scripts/FigS4B_missing_by_patient.R)

**Figure S5.** Most common patterns of missing data in the components of the allostatic load index across $N = 1000$ patients in the extracted EHR data. For example, $432$ patients were missing C-reactive protein, homocysteine, and creatinine clearance. This plot was created using the `naniar` package.

  -  [Figure](figures/FigS5_Missingness_Patterns.png)
  -  [Script (Make Figure)](figure-scripts/FigS5_Missingness_Patterns.R)

**Figure S6.** In simulations, the SMLE for the adjusted log odds ratio on ALI was empirically unbiased (i.e., close to the truth of $\beta_1 = 1.88$, denoted by the dashed line) under varied error rates and different validation study designs. Five designs were considered: (i) simple random sampling (SRS), (ii) case-control sampling (CC) on $Y$, (iii) balanced case-control (BCC) and (iv) optimal (OPT) sampling on $Y$ and $X^*$ (discretized at the median), and (v) residual sampling (RS).

  - [Figure](figures/FigS6_Errors_logOR_Boxplot.png) 
  - [Script (Run Simulations)](sim-scripts/S5.2-Sims-Vary-Error-Rates/)
  - [Data (Simulation Results)](sims-data/S5.2-Sims-Vary-Error-Rates/)
  - [Script (Make Figure)](figure-scripts/FigS6-Errors-Bias-Boxplot.R)

**Figure S7.** In simulations, the targeted validation study designs offered higher efficiency (i.e., smaller variance) for the SMLE log odds ratio estimates on ALI than simple random sampling (SRS) under varied error rates. Four targeted designs were considered: (i) case-control sampling (CC) on $Y$, (ii) balanced case-control (BCC) and (iii) optimal (OPT) sampling on $Y$ and $X^*$ (discretized at the median), and (iv) residual sampling (RS).

  - [Figure](figures/FigS7_Errors_RE_Linegraph.png) 
  - [Script (Run Simulations)](sim-scripts/S5.2-Sims-Vary-Error-Rates/)
  - [Data (Simulation Results)](sims-data/S5.2-Sims-Vary-Error-Rates/)
  - [Script (Make Figure)](figure-scripts/FigS7-Errors-RE-Lineplot.R)

**Figure S8.** In simulations, the SMLE for the adjusted log odds ratio on $X$ was empirically unbiased (i.e., close to the truth of $\beta_1 = 1.88$, denoted by the dashed line) under full/high data recovery rates and different validation study designs. Still, it was closer to the truth for moderate recovery or lower than the naive analysis (i.e., based on the unvalidated EHR data, denoted by the dotted line). Five designs were considered: (i) simple random sampling (SRS), (ii) case-control sampling (CC) on $Y$, (iii) balanced case-control (BCC) and (iv) optimal (OPT) sampling on $Y$ and $X^*$ (discretized at the median), and (v) residual sampling (RS).

  - [Figure](figures/FigS8_Recovery_logOR_Boxplot.png) 
  - [Script (Run Simulations)](sim-scripts/S5.3-Sims-Vary-Data-Recovery/)
  - [Data (Simulation Results)](sims-data/S5.3-Sims-Vary-Data-Recovery/)
  - [Script (Make Figure)](figure-scripts/FigS8-Recovery-Bias-Boxplot.R)

**Figure S9.** In simulations, the targeted validation study designs offered higher efficiency (i.e., smaller variance) for the SMLE log odds ratio estimates than simple random sampling (SRS) under varied data recovery rates. Four targeted designs were considered: (i) case-control sampling (CC) on $Y$, (ii) balanced case-control (BCC) and (iii) optimal (OPT) sampling on $Y$ and $X^*$ (discretized at the median), and (iv) residual sampling (RS).

  - [Figure](figures/FigS9-Recovery-RE-Lineplot.png) 
  - [Script (Run Simulations)](sim-scripts/S5.3-Sims-Vary-Data-Recovery/)
  - [Data (Simulation Results)](sims-data/S5.3-Sims-Vary-Data-Recovery/)
  - [Script (Make Figure)](figure-scripts/FigS9-Recovery-RE-Lineplot.R)

## Tables 

**Table 1.** Simulation results under increasing severity of errors in straight-line proximity to healthy foods, as controlled by the error standard deviation $\sigma_U$.

  - [Script (Run Simulations Locally)](sims-scripts/sims_vary_sigmaU.R)
  - [Script (Make Table)](table-scripts/table1_vary_sigmaU.R)
  - [Data (Simulation Results)](sims-data/vary_sigmaU_sims_combined.csv)
