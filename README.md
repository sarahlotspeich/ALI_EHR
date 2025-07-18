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

## Tables 

**Table 1.** Simulation results under increasing severity of errors in straight-line proximity to healthy foods, as controlled by the error standard deviation $\sigma_U$.

  - [Script (Run Simulations Locally)](sims-scripts/sims_vary_sigmaU.R)
  - [Script (Make Table)](table-scripts/table1_vary_sigmaU.R)
  - [Data (Simulation Results)](sims-data/vary_sigmaU_sims_combined.csv)

**Table 2.** Simulation results under increasing proportion of neighborhoods queried to obtain map-based proximity to healthy foods, as controlled by $p_V$.

  - [Script (Run Simulations Locally)](sims-scripts/sims_vary_pV.R)
  - [Script (Make Table)](table-scripts/table2_vary_pV.R)
  - [Data (Simulation Results)](sims-data/vary_pV_sims_combined.csv)

**Table 3.** Simulation results under higher disease prevalence and prevalence ratios for map-based proximity to healthy foods, as controlled by the coefficients $\beta_0$ and $\beta_1$, respectively.

  - [Script (Run Simulations Locally)](sims-scripts/sims_vary_prev.R)
  - [Script (Make Table)](table-scripts/table3_vary_prev.R)
  - [Data (Simulation Results)](sims-data/vary_prev_sims_combined.csv)

**Table 4.** Simulation results under increasingly severe multiplicative errors in straight-line proximity to healthy foods, as controlled by the max of the error distribution $\tau_W$. 

  - [Script (Run Simulations Locally)](sims-scripts/sims_mult_error.R)
  - [Script (Make Table)](table-scripts/table4_mult_error.R)
  - [Data (Simulation Results)](sims-data/mult_error_sims_combined.csv)

**Table S1.** Simulation results under varied additive errors in straight-line proximity to healthy foods, as controlled by the mean $\mu_U$ of the errors $U$. The standard deviation $\sigma_U = 0.8$ of the errors was fixed.

  - [Script (Run Simulations Locally)](sims-scripts/sims_vary_muU.R)
  - [Script (Make Table)](table-scripts/tableS1_vary_muU.R)
  - [Data (Simulation Results)](sims-data/vary_muU/)

**Table S2.** Simulation results under varied multiplicative errors in straight-line proximity to healthy foods, as controlled by the mean $\mu_W$ of the errors $W$. The standard deviation $\sigma_W = 0.15$ of the errors was fixed.

  - [Script (Run Simulations Locally)](sims-scripts/sims_mult_error2.R)
  - [Script (Make Table)](table-scripts/tableS2_mult_errors2.R)
  - [Data (Simulation Results)](sims-data/mult_error2/)

**Table S3.** Descriptive statistics of the $N = 387$ census tracts in the Piedmont Triad, North Carolina.

  - [Script (Make Table)](table-scripts/tableS3_piedmont.R)
  - [Data (Food Access + Health + RUCA)](piedmont-triad-data/analysis_data.csv)

**Table S4.** Simulation results for the mixed-effects model under additive errors in straight-line proximity to healthy foods, with fixed $\mu_U = -0.7$ and $\sigma_U = 0.8$, and different amounts of variability in the spatial random effect, as controlled by $\tau^2$. 

  - [Script (Run Simulations Locally)](sims-scripts/sims_spatial.R)
  - [Script (Make Table)](table-scripts/tableS4_spatial_vary_tau2.R)
  - [Data (Simulation Results)](sims-data/spatial/)

## Figures 

**Figure 2.** Choropleth maps of the crude prevalence of diagnosed diabetes and obesity for census tracts in the Piedmont Triad, North Carolina.

  - [Figure](figures/fig2_map_piedmont_triad_health_outcomes.png)
  - [Data (PLACES)](piedmont-triad-data/disease_prevalences_2022.csv)
  - [Script (Make Figure)](figure-scripts/fig2_map_piedmont_triad_health_outcomes.R)

**Figure 3.** Scatter plot of straight-line versus map-based proximity to healthy foods store for neighborhoods in the Piedmont Triad, North Carolina using the fully-queried data ($N = 387$).

  - [Figure](figures/fig3_scatterplot_proximity_piedmont.png)
  - [Data (Food Access + Health + RUCA)](piedmont-triad-data/analysis_data.csv)
  - [Script (Make Figure)](figure-scripts/fig3_scatterplot_proximity_piedmont.R)

**Figure 4.** Estimated prevalence ratios (with 95\% confidence intervals) for proximity to healthy foods, metropolitan status, and their interaction with the two health outcomes in the Piedmont Triad, North Carolina using four different analysis methods. Within each health outcome and method, estimates on the right with the dashed error bars came from the mixed effects model allowing for spatial autocorrelation between neighboring census tracts; estimates on the left with the solid error bars came from the non-spatial model assuming independence between tracts.

  - [Figure](figures/fig4_forest_plot_piedmont.png)
  - [Data (Food Access + Health + RUCA)](piedmont-triad-data/analysis_data.csv)
  - [Data (Fitted Models)](piedmont-triad-data/forest_plot_data.csv) 
  - [Script (Make Figure)](figure-scripts/fig4_forest_plot_piedmont.R)

**Figure S1.** Straight-line and map-based distances from Reynolda House (square symbol) to a nearby Food Lion grocery store (triangle symbol) in Winston-Salem, North Carolina.

  - [Figure](figures/figS1_map_comparing_distances.png)
  - [Script (Make Figure)](figure-scripts/figS1_map_comparing_distances.R)

**Figure S2.** Line graph of the cumulative computing time (in seconds) for the map-based versus straight-line distance calculations in the Piedmont Triad data.

  - [Figure](figures/figS2_cum_comp_time_line.png) 
  - [Data (Food Access + Health + RUCA)](piedmont-triad-data/analysis_data.csv)
  - [Script (Make Figure)](figure-scripts/figS2_cum_comp_time_line.R)

**Figure S3.** Estimated prevalence ratios for map-based food access $X$ on health using multiple imputation. The five possible ways to include the analysis model outcome $Y$ (with or without the model offset $Pop$) in the imputation model for $X$ were considered.

  - [Figure](figures/figS3_incl_in_imputation_model_PR.png) 
  - [Script (Run Simulations Locally)](sims-scripts/sims_incl_outcome.R)
  - [Data (Simulation Results)](sims-data/include_outcome/)
  - [Script (Make Figure)](figure-scripts/figS3_inclY_in_imputation_model.R)

**Figure S4.** Map of the $100$ North Carolina counties, colored by whether they belong to the Piedmont Triad ($n = 12$), border the Piedmont Triad ($n = 12$), or fall into the rest of the state ($n = 76$).

  - [Figure](figures/figS4_map_piedmont_triad.pdf)
  - [Script (Make Figure)](figure-scripts/figS4_map_piedmont_triad.R)

**Figure S5.** Choropleth maps of socioeconomic factors across the census tracts of the Piedmont Triad, North Carolina. Data were taken from the 2015 American Community Survey. The gradient for each map is centered at the state median. Two census tracts were missing median family income. All other maps are based on $N = 387$ census tracts. The one census tract with zero population was excluded.

  - [Figure](figures/figS5_map_piedmont_triad_acs_data.png)
  - [Data (RUCA)](piedmont-triad-data/ruca2010revised.csv)
  - [Data (ACS)](piedmont-triad-data/piedmont_triad_acs_data.csv)
  - [Script (Make Figure)](figure-scripts/figS5_map_piedmont_triad_acs_data.R)

**Figure S6.** Choropleth maps of percents of census tract population by self-reported race in the Piedmont Triad, North Carolina. Data were taken from the 2015 American Community Survey. The gradient for each map is centered at the state median. All maps are based on $N = 387$ census tracts. The one census tract with zero population was excluded.

  - [Figure](figures/figS4_map_piedmont_triad_acs_data.png)
  - [Data (ACS)](piedmont-triad-data/piedmont_triad_acs_data.csv)
  - [Script (Make Figure)](figure-scripts/figS6_map_piedmont_triad_acs_race_data.R)

**Figure S7.** Map of $M = 701$ authorized SNAP retailers in the Piedmont Triad, North Carolina, broken down by store type. Data were taken from the 2022 Historical SNAP Retail Locator Data.

  - [Figure](figures/figS7_map_piedmont_triad_SNAP_wide.png)
  - [Data (SNAP)](piedmont-triad-data/healthy_foods_stores_2022.csv)
  - [Script (Make Figure)](figures/figS7_map_SNAP.R)

**Figure S8.** Map of census tracts in the Piedmont Triad, North Carolina, colored according to whether it was treated as queried in the partially queried analysis. For the partially queried analysis, $n = 48$ tracts were chosen to be queried via county-stratified random sampling. The thicker boundaries denote outline the $12$ counties.

  - [Figure](figures/figS8_map_piedmont_queried.png) 
  - [Data (Food Access + Health + RUCA)](piedmont-triad-data/analysis_data.csv)
  - [Script (Make Figure)](figure-scripts/figS8_map_piedmont_queried.R)

**Figure S9.** Choropleth maps of the crude prevalence of adverse health outcomes for census tracts in Forsyth County (top row) and Guilford County (bottom row), North Carolina. The triangles denote "downtown"  Winston-Salem and Greensboro (represented by their City Halls) in the top and bottom row, respectively. The gradient for each map is centered at the state median. Data were taken from the 2022 PLACES dataset.

  - [Figure](figures/figS9_map_forsyth_guilford_health_outcomes.png)
  - [Data (PLACES)](piedmont-triad-data/disease_prevalences_2022.csv)
  - [Script (Make Figure)](figure-scripts/figS7_map_forsyth_guilford_health_outcomes.R)

**Figure S10.** Scatter plot of straight-line versus map-based proximity to healthy food store for neighborhoods (census tracts) in the Piedmont Triad, North Carolina using the fully queried data ($N = 387$) or the partially queried data ($n = 48$). The top row is among only metropolitan census tracts, and the bottom row is only among non-metropolitan census tracts. The solid line follows the fitted least-squares linear regression fit between $X$ and $X ^ {\ast}$ among those tracts, while the dashed line denotes the hypothetical $X = X^{\ast}$ if there had been no errors in $X^{\ast}$.

  - [Figure](figS10_scatterplot_proximity_piedmont_metro.pdf)
  - [Data (Food Access + Health + RUCA)](piedmont-triad-data/analysis_data.csv)
  - [Script (Make Figure)](figure-scripts/figS10_scatterplot_proximity_piedmont_metro.R)

**Figure S11.** Histogram of additive errors ($U$) and multiplicative errors ($W$) in straight-line proximity to healthy foods ($X^*$) from the fully queried data ($N = 387$) for the Piedmont Triad, North Carolina.

  - [Figure](figures/figS11_histogram_errors_piedmont.png)
  - [Data (Food Access + Health + RUCA)](piedmont-triad-data/analysis_data.csv)
  - [Script (Make Figure)](figure-scripts/figS11_histogram_errors_proximity.R)

**Figure S12.** Estimated baseline prevalence (with 95\% confidence intervals) for health outcomes in the Piedmont Triad, North Carolina using different analysis methods. Within each health outcome and method, estimates on the right came from the mixed effects model allowing for spatial autocorrelation between bordering neighborhoods (census tracts); estimates on the left came from the model assuming independence between neighborhoods.

  - [Figure](figures/figS12_forest_plot_intercept_piedmont.png)
  - [Data (Food Access + Health)](piedmont-triad-data/analysis_data.csv)
  - [Data (Fitted Models)](piedmont-triad-data/forest_plot_data.csv) 
  - [Script (Make Figure)](figure-scripts/figS12_forest_plot_intercept_piedmont.R)
