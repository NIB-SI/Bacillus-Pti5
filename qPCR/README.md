# 🧪 qPCR and qPCR-related assays 


## 📁 structure

```text
├── input/              # Processed qPCR data and other input tables
│   └── Standardised qPCR gene expression data
│   └── The amount of biofilm formation on potato roots
│   └── B. subtilis abundance in potato plants
│   └── Mycorrhizal colonization
│   └── Flagellin treatment
│   └── ...
├── other/              # Miscellaneous files
├── output/             # Output tables from scripts as .tsv and restructured analogues in .xlsx
├── reports/            # Separate `.pdf` and `.svg` plots
├── scripts/            # `.Rmd` scripts and `.html` reports
├── README.md           # Project overview
```   


## 👩‍🔬 about qPCR data
 - [MIQE 2.0: Revision of the Minimum Information for Publication of Quantitative Real-Time PCR Experiments Guidelines](https://academic.oup.com/clinchem/article/71/6/634/8119148)
 - [quantGenius - quantification of qPCR data using standard curve](http://quantgenius.nib.si)
   - [Equations used in quantGenius workflow](https://static-content.springer.com/esm/art%3A10.1186%2Fs12859-017-1688-7/MediaObjects/12859_2017_1688_MOESM2_ESM.pdf)
 - quantification based on a standard curve inherently involves a log transformation of the input data
 - LOQ values skew distributions

## 🖥️ R packages of interest    
- [MKinfer](https://stamats.r-universe.dev/MKinfer)
- [exactRankTests](https://thothorn.r-universe.dev/exactRankTests)
- [rstatix](https://kassambara.r-universe.dev/rstatix)
- [Rfast](https://rfastofficial.r-universe.dev/Rfast)
- [emmeans](https://rvlenth.r-universe.dev/emmeans)
- [multcomp](https://r-forge.r-universe.dev/multcomp)
- [ggplot2](https://tidyverse.r-universe.dev/ggplot2)
- [ggplubr](https://kassambara.r-universe.dev/ggpubr)
- [egg](https://baptiste.r-universe.dev/egg)
- ...
See sessionInfo output for package version.

## ⚓ concepts
### distribution, variance and effects
#### 📙 Distribution
  - various robustness under varying skewness and kurtosis
  - can have low power for small sample size
  - avoid overinterpreting p-values — focus on effect size and distribution shape

| Test                 |  Type                          | Sensitivity to Tails | R |
|----------------------|-------------------------------|----------------------|-|
| Shapiro–Wilk         | Parametric                    | Moderate             | `shapiro.test()` |
| Lilliefors (KS)      | Non-parametric (KS variant)   | High                 | `nortest::lillie.test()`|
| Anderson–Darling     | Empirical distribution-based  | Very High           | `nortest::ad.test()` |
| Jarque–Bera          | Parametric                    | Low–Moderate        | `tseries::jarque.bera.test()`|
| D’Agostino Skewness  | Parametric                    | Focused on skewness | `moments::agostino.test()`|

  - 📚 **Q-Q (Quantile-Quantile) plots and a Residual plots ** `ggqqplot` 🔎
      - if most or all points fall inside the shaded confidence band, the sample’s distribution does not show strong evidence of departure from normality at the plotted sample size and confidence level
      - a few isolated points outside the band at the extremes are common with small samples and do not necessarily indicate a severe problem

#### 📘 Heteroscedasticity
- some tests are meant to be used with normally distributed data, but can tolerate relatively low deviation from normality

| Test               | Assumptions                     | Robustness to Outliers | R 
|--------------------|----------------------------------|-------------------------|-|
| Levene’s Test      | Normality (mean-based)           | Moderate                |`levene_test`|
| Brown-Forsythe     | Normality (median-based)         | High                    |`levene_test(center = median)`|
| **Fligner-Killeen**    | Non-parametric                   | Very high               |`fligner.test` 🔎|

#### 📗 Effect Size

| Metric               | Type           | Notes                        |R|
|----------------------|----------------|------------------------------|-|
| Cohen’s d            | Parametric     | Sensitive to outliers        |`cohens_d`|
| Hedges’ g            | Parametric     | Smaller sample size        |`cohen.d(hedges.correction = TRUE)`|
| **Wilcoxon Effect Size** | Non-parametric | More robust alternative      |`wilcox_effsize` 🔎|

- Cohen’s d or Hedges’ g are not used when reporting results from a non-parametric test
- 
### 📔 permutation t-test
- Non-parametric alternative to traditional t-tests
- Suitable for:
  - Small sample sizes
  - Unequal variances
  - Non-i.i.d. data
- Robust and flexible for exploratory comparisons
- `MKinfer::perm.t.test`, `exactRankTests::perm.test`
- `stats::p.adjust(method = "BH")`

### 📜 ANOVA 
- Parametric, non-parametric, and permutation-based variants
- Running ANOVA assumes you're testing whether group means differ significantly, under the null hypothesis that all groups are equal
- When experimental design intentionally introduces differences violates the null hypothesis from the start
- Can obscure meaningful planned contrasts or targeted comparisons
- Post-hoc tests (e.g., Dunnett, Games-Howell, emmeans) are often more informative and aligned with design

### 📔 Games-Howell Post-hoc Test
- `rstatix::games_howell_test`
- Compares all group pairs when variance homogeneity is violated
- Based on:
  - Welch’s degrees of freedom correction
  - Tukey’s studentized range distribution
- Features:
  - Confidence intervals for mean differences
  - No need for additional p-value corrections
  - Assumes approximate normality of residuals 

### 🗞️ Hypothesis Test for Two Means of Percentages
- Use `Rfast::percent.ttest()` for pairwise proportion comparisons
- Assumes beta-distributed data
- Must be run manually for each pair

### 📔 Beta Regression & EMMs
- Use `betareg::betareg()` for modeling data bounded between 0 and 1 (i.e., proportions or percentages)
- Use `emmeans::emmeans()` for estimated marginal means (EMMs)
- Use `emmeans::joint_tests()` to test overall factor effects
- EMMs are interpretable as **model-adjusted group means**
### Multiple Comparison Adjustments
| Method     | Best For                      | Notes                                  |
|------------|-------------------------------|----------------------------------------|
| Šidák      | All pairwise comparisons      | Controls family-wise error rate        |
| Dunnettx   | One-vs-control comparisons    | More powerful for targeted contrasts   |

### 📔 Linear mixed-effects models
- For nested structures
- ```lme4::lmer```, followed by ```MuMIn::r.squaredGLMM```, diagnostics plots, Estimated marginal means (EMMs) ```emmeans::emmeans```, ```emmeans::contras``` and ```multcomp::cld```

### 📔 Bayesian mixed-effects logistic regression models
- ```rstanarm::stan_glmer```, followed by diagnostics (various), and Estimated marginal means (EMMs)
- Because Bayesian models do not yield traditional p-values (```bayestestR::describe_posterior```), try Fisher’s Exact Test for count data to facilitate interpretation for broader audiences


___
🦘 Suggestions:
- [DHARMa package](https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html)
- [sandwich](https://cran.r-project.org/web/packages/sandwich/refman/sandwich.html)
