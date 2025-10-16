🧬 qPCR assays

# qPCR data
 - [MIQE 2.0: Revision of the Minimum Information for Publication of Quantitative Real-Time PCR Experiments Guidelines](https://academic.oup.com/clinchem/article/71/6/634/8119148)
 - [quantGenius - quantification of qPCR data using standard curve](http://quantgenius.nib.si)
   - [Equations used in quantGenius workflow](https://static-content.springer.com/esm/art%3A10.1186%2Fs12859-017-1688-7/MediaObjects/12859_2017_1688_MOESM2_ESM.pdf)
 - quantification based on a standard curve inherently involves a log transformation of the input data
 - LOQ values skew distributions

# R packages    
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

# concepts
## distribution, variance and effects
### Distribution
  - various robustness under varying skewness and kurtosis
  - can have low power for small sample size

| Test Name            | R                  | Type                          | Sensitivity to Tails | 
|----------------------|-----------------------------|-------------------------------|----------------------|
| Shapiro–Wilk     | `shapiro.test()`            | Parametric                    | Moderate             | 
| Lilliefors (KS)  | `nortest::lillie.test()`    | Non-parametric (KS variant)   | High                 | 
| Anderson–Darling | `nortest::ad.test()`        | Empirical distribution-based  | Very High           | 
| Jarque–Bera      | `tseries::jarque.bera.test()` | Parametric                    | Low–Moderate        | 
| D’Agostino Skewness | `moments::agostino.test()` | Parametric                    | Focused on skewness | 

  - Q-Q (Quantile-Quantile) plots and a Residual plots `ggqqplot`
      - if most or all points fall inside the shaded confidence band, the sample’s distribution does not show strong evidence of departure from normality at the plotted sample size and confidence level
      - a few isolated points outside the band at the extremes are common with small samples and do not necessarily indicate a severe problem

### Heteroscedasticity
- some tests are meant to be used with normally distributed data, but can tolerate relatively low deviation from normality

| Test               | Assumptions                     | Robustness to Outliers | R 
|--------------------|----------------------------------|-------------------------|-|
| Levene’s Test      | Normality (mean-based)           | Moderate                |`levene_test`|
| Brown-Forsythe     | Normality (median-based)         | High                    |`levene_test(center = median)`|
| Fligner-Killeen    | Non-parametric                   | Very high               |`fligner.test`|

### Effect Size

| Metric               | Type           | Notes                        |R|
|----------------------|----------------|------------------------------|-|
| Cohen’s d            | Parametric     | Sensitive to outliers        |`cohens_d`|
| Wilcoxon Effect Size | Non-parametric | More robust alternative      |`wilcox_effsize`|

## permutation t-test
- Non-parametric alternative to traditional t-tests
- Suitable for:
  - Small sample sizes
  - Unequal variances
  - Non-i.i.d. data
- Robust and flexible for exploratory comparisons
- `MKinfer::perm.t.test`, `exactRankTests::perm.test`

## Games-Howell Post-hoc Test
- `games_howell_test`
- Compares all group pairs when variance homogeneity is violated
- Based on:
  - Welch’s degrees of freedom correction
  - Tukey’s studentized range distribution
- Features:
  - Confidence intervals for mean differences
  - No need for additional p-value corrections
  - Assumes approximate normality of residuals 

## Hypothesis Test for Two Means of Percentages
- Use `Rfast::percent.ttest()` for pairwise proportion comparisons
- Assumes beta-distributed data
- Must be run manually for each pair

## Beta Regression & EMMs
- Use `betareg::betareg()` for modeling data bounded between 0 and 1 (i.e., proportions or percentages)
- Use `emmeans::emmeans()` for estimated marginal means (EMMs)
- Use `emmeans::joint_tests()` to test overall factor effects
- EMMs are interpretable as **model-adjusted group means**
### Multiple Comparison Adjustments
| Method     | Best For                      | Notes                                  |
|------------|-------------------------------|----------------------------------------|
| Šidák      | All pairwise comparisons      | Controls family-wise error rate        |
| Dunnettx   | One-vs-control comparisons    | More powerful for targeted contrasts   |
  

