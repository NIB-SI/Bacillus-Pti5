# qPCR

- distribution, variance and effects
  - on qPCR data
    - [quantGenius - quantification of qPCR data using standard curve](http://quantgenius.nib.si/user/login)
    - [MIQE 2.0: Revision of the Minimum Information for Publication of Quantitative Real-Time PCR Experiments Guidelines](https://academic.oup.com/clinchem/article/71/6/634/8119148)
    - quantification based on a standard curve inherently involves a log transformation of the input data
  - *_Denote_*: LOQ values skew distributions
  - *_Denote_*: distribution tests: low power for small sample size, see QQ of residuals
    - if most or all points fall inside the shaded confidence band, the sample’s distribution does not show strong evidence of departure from normality at the plotted sample size and confidence level
    - a few isolated points outside the band at the extremes are common with small samples and do not necessarily indicate a severe problem
  - *_Denote_*: heteroscedasticity: some tests are meant to be used with normally distributed data, but can tolerate relatively low deviation from normality; Levene’s test with mean, more robust test Brown-Forsythe with median, Fligner-Killeen when data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved
  - *_Denote_*: effect size: Cohen’s d - parametric with assumptions, Wilcoxon Effect Size - the non-parametric alternative; different sensitivity to outliers

- permutation t-test
  - for non-i.i.d, hypotheses
  - a non-parametric robust alternative to traditional t-tests
  - small sample sizes
  - robust to unequal variances

- Games Howell Post-hoc Tests
  - compare all possible combinations of group differences when the assumption of homogeneity of variances is violated
  - provides confidence intervals for the differences between group means and shows whether the differences are statistically significant
  - based on Welch’s degrees of freedom correction and uses Tukey’s studentized range distribution for computing the p-values
  - compares the difference between each pair of means with appropriate adjustment for the multiple testing, so there is no need to apply additional p-value corrections
  - assumes approximate normality of residuals

- hypothesis test for two means of percentages
  - performs a two-sample test for proportions (like a t-test for percentages), assuming beta-distributed data
  - must manually run for each pair

- beta regression & emmeans
  - analyzing percentage data using beta regression followed by estimated marginal means (EMMs) and pairwise comparisons
  - Beta Regression is well-suited for modeling data bounded between 0 and 1 (i.e., proportions or percentages)
  - Joint Tests is useful for testing whether the factor Genotype has a significant overall effect
  - Estimated Marginal Means – EMMs are interpretable as model-adjusted group means.
  - Sidak adjustment for controlling family-wise error rate - sdjusts p-values or alpha levels for multiple comparisons; Dunnettx is designed for one-vs-control comparisons and is more statistically powerful in that context

