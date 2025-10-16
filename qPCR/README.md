# qPCR

# qPCR data
 - [quantGenius - quantification of qPCR data using standard curve](http://quantgenius.nib.si)
   - [Equations used in quantGenius workflow](https://static-content.springer.com/esm/art%3A10.1186%2Fs12859-017-1688-7/MediaObjects/12859_2017_1688_MOESM2_ESM.pdf)
 - [MIQE 2.0: Revision of the Minimum Information for Publication of Quantitative Real-Time PCR Experiments Guidelines](https://academic.oup.com/clinchem/article/71/6/634/8119148)
 - quantification based on a standard curve inherently involves a log transformation of the input data
 - LOQ values skew distributions

# R packages    
- [MKinfer](https://stamats.r-universe.dev/MKinfer)
- [exactRankTests](https://thothorn.r-universe.dev/exactRankTests)
- [rstatix]([https://rpkgs.datanovia.com/rstatix/](https://kassambara.r-universe.dev/rstatix))
- [Rfast](https://rfastofficial.r-universe.dev/Rfast)
- [emmeans](https://rvlenth.r-universe.dev/emmeans)
- [multcomp](https://r-forge.r-universe.dev/multcomp)
- [ggplot2](https://tidyverse.r-universe.dev/ggplot2)
- [ggplubr](https://kassambara.r-universe.dev/ggpubr)
- ...

# concepts
## distribution, variance and effects
### Distribution
  - various robustness under varying skewness and kurtosis
    - Shapiro–Wilk
    - Lilliefors
    - Anderson–Darling (Empirical distribution-based test)
    - Jarque–Bera
    - D’Agostino Skewness 
  - low power for small sample size,
  - Q-Q (Quantile-Quantile) plots and a Residual plots
      - if most or all points fall inside the shaded confidence band, the sample’s distribution does not show strong evidence of departure from normality at the plotted sample size and confidence level
      - a few isolated points outside the band at the extremes are common with small samples and do not necessarily indicate a severe problem
### Heteroscedasticity: some tests are meant to be used with normally distributed data, but can tolerate relatively low deviation from normality
  - Levene’s test with mean, under asumptions
  - more robust test Brown-Forsythe with median, under asumptions
  - Fligner-Killeen when data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved
### Effect Size
  - Cohen’s d - parametric, under assumptions, sensitive to outliers
  - Wilcoxon Effect Size - the non-parametric alternative; different sensitivity to outliers

## permutation t-test
- Non-parametric alternative to traditional t-tests
- Suitable for:
  - Small sample sizes
  - Unequal variances
  - Non-i.i.d. data
- Robust and flexible for exploratory comparisons

## Games-Howell Post-hoc Test
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
  

