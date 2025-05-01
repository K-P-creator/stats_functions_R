# stats_functions_R

My own function library for the STAT411 course with R 
Theres a lot of functions in here to assist with the coursework and the final exam.
My preffered method for using these functions is to manually source the file.

Use:

```ps
setwd("C:/directory_with_functions.R")
source("functions.r")
```

to add the functions to your current instance of R studio. 

The following functions summary is AI generated so proceed with skepticism

## Functions

- **t_stat**: Computes the t-statistic for a sample mean versus a population mean.
- **t_stat_p**: Computes the z-score for a sample proportion versus a population proportion.
- **sse**: Calculates the sum of squared errors (SSE) for a linear regression model.
- **cv_xbar**: Returns the critical t-value for a given confidence level and sample size.
- **ci_xbar**: Computes the confidence interval for a sample mean using the t-distribution.
- **cv_phat**: Returns the critical z-value for a confidence level in proportion testing.
- **ci_phat**: Computes the confidence interval for a population proportion.
- **perm**: Calculates the number of permutations of r items from a set of n.
- **comb**: Calculates the number of combinations of r items from a set of n.
- **drvmean**: Computes the expected value of a discrete random variable.
- **drvvar**: Computes the variance of a discrete random variable.
- **drvdev**: Computes the standard deviation of a discrete random variable.
- **z_score**: Calculates the z-score of a given value.
- **prob_phat**: Calculates tail probability for a sample proportion under normal approximation.
- **me_phat**: Calculates the probability that a sample proportion is within or beyond a margin of error.
- **prob_xbar**: Calculates tail probability for a sample mean under CLT normal approximation.
- **me_xbar**: Computes the probability that a sample mean is within or beyond a margin of error.
- **cv_sd**: Returns chi-square critical values for confidence intervals on standard deviation or variance.
- **ci_sd**: Computes the confidence interval for a population standard deviation.
- **ci_var**: Computes the confidence interval for a population variance.
- **sd_binom**: Calculates the standard deviation of a binomial distribution.
- **binom_var**: Calculates the variance of a binomial distribution.
- **t_stat_two_means_known_pop_sd**: Computes z-test statistic for two means with known population std devs.
- **se_two_means_known_pop_sd**: Computes standard error for two means with known population std devs.
- **ci_two_means_known_pop_sd**: Computes confidence interval for difference in means with known population std devs.
- **sd_two_means_equal_var**: Computes pooled standard deviation assuming equal variances.
- **t_stat_two_means_equal_var**: Computes t-statistic for two means assuming equal variances.
- **se_two_means_equal_var**: Computes standard error for two means assuming equal variances.
- **ci__two_means_equal_var**: Computes confidence interval for difference in means assuming equal variances.
- **p_pooled**: Computes pooled proportion from two independent samples.
- **t_stat_p_pooled**: Computes test statistic for difference in proportions using pooled estimate.
- **se_p_pooled**: Computes standard error for difference in proportions using pooled estimate.
- **ci_p_pooled**: Computes confidence interval for difference in proportions using pooled estimate.
- **t_stat_paired_means**: Computes t-statistic for paired sample means.
- **se_paired_means**: Computes standard error for paired sample means.
- **ci_paired_means**: Computes confidence interval for paired sample means.
- **se_two_means_unknown_pop_sd**: Computes standard error for two means with unknown population std devs.
- **t_stat_two_means_unknown_pop_sd**: Computes t-statistic for two means with unknown population std devs.
- **ci_two_means_unknown_pop_sd**: Computes confidence interval for difference in means with unknown population std devs.
