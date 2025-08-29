# qPCR

# general concepts

descriptive statistics
- test of normality
  - [Shapiro-Wilk test](https://rpkgs.datanovia.com/rstatix/reference/shapiro_test.html)
- test for homogeneity of variance
- [Levene test](https://rpkgs.datanovia.com/rstatix/reference/levene_test.html)

effect size
- [Wilcoxon effect size](https://rpkgs.datanovia.com/rstatix/reference/wilcox_effsize.html)
- [Cohen's d](https://rpkgs.datanovia.com/rstatix/reference/cohens_d.html)

t-test
- Welch's t-test -- unequal variances t-test
- [Permutation test](https://cran.r-project.org/web//packages/MKinfer/vignettes/MKinfer.html#permutation-t-test) -- a subset of non-parametric statistics

Multiple testing correction
- Family-wise error rate
  - Holm's procedure
  - Dunnett's correction
  - ...
- False discovery rate
  - Benjamini–Hochberg procedure
  - Benjamini–Yekutieli procedure
  - ...
