t_stat <- function(xbar, mu, sd, n) {
    return ((xbar-mu)/(sd/sqrt(n)))
}

t_stat_p <- function(phat, p, n){
    return ((phat - p)/sqrt((p*(1-p))/n))
}

#' Compute Sum of Squared Errors (SSE) for a given linear model
#'
#' This function calculates the sum of squared errors between actual values and predicted values
#' from a linear regression model of the form: y_hat = intercept + slope * x
#'
#' @param x A numeric vector of predictor (independent variable) values
#' @param y A numeric vector of actual outcome (dependent variable) values
#' @param intercept A numeric value representing the intercept (β₀) of the regression model
#' @param slope A numeric value representing the slope (β₁) of the regression model
#'
#' @return The Sum of Squared Errors (SSE), rounded to 4 decimal places
#' @examples
#' x <- c(3.34, 2.14, 2.09, 2.93, 2.26, 3.66)
#' y <- c(2.78, 3.70, 2.27, 3.47, 3.14, 3.95)
#' sse(x, y, intercept, slope)
sse <- function(x, y, intercept, slope) {
  predicted <- intercept + slope * x
  sse <- sum((y - predicted)^2)
  return(sse)
}

#' Critical Value for a confidence level and n with t-dist
#'
#' @param n Integer. The sample size (must be greater than 0).
#' @param conf_level Numeric. The confidence level (must be > 0 and <= 1, e.g., 0.80 for 80% confidence).
#'
#' @return CV for the given input
#'
#' @export
cv_xbar <- function (n, conf_level) {
    df <- n - 1
    t_crit <- qt(1 - (1 - conf_level) / 2, df)
    return (t_crit)
}


#' Confidence Interval for the Mean using t-distribution
#'
#' Calculates the confidence interval for a population mean
#'
#' @param xbar Numeric. The sample mean.
#' @param sd Numeric. The sample standard deviation.
#' @param n Integer. The sample size (must be greater than 0).
#' @param conf_level Numeric. The confidence level (must be > 0 and <= 1, e.g., 0.80 for 80% confidence).
#'
#' @return A vector with the lower and upper bounds of the confidence interval.
#'
#' @export
ci_xbar <- function(xbar, sd, n, conf_level) {
  if (n <= 0) stop("Sample size 'n' must be greater than 0.")
  if (conf_level <= 0 || conf_level > 1) stop("Confidence level must be between 0 and 1.")
  
  df <- n - 1
  t_crit <- qt(1 - (1 - conf_level) / 2, df)
  margin_error <- t_crit * sd / sqrt(n)
  lowerbound <- xbar - margin_error
  upperbound <- xbar + margin_error
  return(c(lowerbound, upperbound))
}


#' Critical Value for a confidence level and n with t-dist
#'
#' @param n Integer. The sample size (must be greater than 0).
#' @param conf_level Numeric. The confidence level (must be > 0 and <= 1, e.g., 0.80 for 80% confidence).
#'
#' @return CV for the given input
#'
#' @export
cv_phat <- function(n, conf_level) {  
  z_crit <- qnorm(1 - (1 - conf_level) / 2)
  return(z_crit)
}


#' Confidence Interval for a Population Proportion
#'
#' Calculates the confidence interval for a population proportion using the normal approximation.
#'
#' @param phat Numeric. The sample proportion (e.g., 0.1 for 10%).
#' @param n Integer. The sample size.
#' @param conf_level Numeric. The confidence level (e.g., 0.98 for 98% confidence).
#'
#' @return A vector with the lower and upper bounds of the confidence interval.
#' 
#' @export
ci_phat <- function(phat, n, conf_level) {
  if (n <= 0) stop("Sample size 'n' must be greater than 0.")
  if (conf_level <= 0 || conf_level > 1) stop("Confidence level must be between 0 and 1.")
  
  z <- cv_phat(n, conf_level)


  se <- sqrt(phat * (1 - phat) / n)

  lower <- phat - z * se
  upper <- phat + z * se

  return(c(lower, upper))
}

perm <- function (n,r){
  if (n < r){
    return ("r cannot be greater than n!")
  }
  return (factorial(n)/factorial(n-r))
}

comb <- function (n,r){
  return(factorial(n)/(factorial(r)*factorial(n-r)))
}

#' Calculate expected value for a discreet random variable given a PDF
#'
#' NOTE - Length of x must equal lenght of prob, and prob must sum to 1
#'
#' @param x The set of all possilbe values
#' @param r The set of probabilities for all expected values
#' @return The expected value (mean)
#' @export
drvmean <- function (x, prob){
    if (length(x) != length(prob)){
        return("X set and Probability sets must have the same length")
    }

    if (sum(prob) != 1){
        return("Sum of Probabilty set must equal 1")
    }

    mean_val <- sum(x * prob)
    return(mean_val)
}


#' Calculate variance for a discreet random variable given a PDF
#'
#' NOTE - Length of x must equal lenght of prob, and prob must sum to 1
#'
#' @param x The set of all possilbe values
#' @param r The set of probabilities for all expected values
#' @return The variance for the sets
#' @export
drvvar <- function (x, prob){
    if (length(x) != length(prob)){
        return("X set and Probability sets must have the same length")
    }

    if (sum(prob) != 1){
        return("Sum of Probabilty set must equal 1")
    }

    mean_val <- drvmean(x,prob)
    var_val <- sum ((x^2) * prob)
    var_val = var_val - mean_val^2

    return(var_val)
}


#' Calculate standard deviation for a discreet random variable given a PDF
#'
#' NOTE - Length of x must equal lenght of prob, and prob must sum to 1
#'
#' @param x The set of all possilbe values
#' @param r The set of probabilities for all expected values
#' @return The standard deviation for the sets
#' @export
drvdev <- function (x,prob){
    dev <- drvvar(x,prob)
    dev <- sqrt(dev)

    return(dev)
}


#' Calculate z score given x, mean and std dev
#'
#'
#' @param x the given x
#' @param mean the mean of the dist
#' @param sd the standard deviation of the dist
#' @return the z score for the given inputs
#' @export
z_score <- function(x, mean, sd) {
  (x - mean) / sd
}


#' Calculate probability that sample proportion is greater or less than given value
#'
#' Uses the normal approximation to compute P(phat > phat_val) or P(phat < phat_val)
#'
#' @param p the population proportion
#' @param phat the sample proportion value to test against
#' @param n the sample size
#' @param lower.tail logical; if TRUE, returns P(phat < phat_val), otherwise P(phat > phat_val)
#' @return the probability that the sample proportion is on the specified tail of phat (rounded to 4 dig)
#' @export
prob_phat <- function(p, phat, n, lower.tail = TRUE) {
  se <- sqrt(p * (1 - p) / n)                    # Standard error AKA sd phat
  z <- z_score(phat, mean = p, sd = se)          
  prob <- pnorm(z, lower.tail = lower.tail)     

  return(round(prob,4))
}


#' Calculate probability that sample proportion is within or greater than a given margin of error
#'
#' Uses the normal approximation to compute P(phat > phat_val) or P(phat < phat_val) and find the difference or sum
#'
#' @param p the population proportion
#' @param me the margin of error
#' @param n the sample size
#' @param lower.tail logical; if TRUE, find the prob that phat is within the ME, otherwise
#'                    prob that phat is outside of it
#' @return the probability that the sample proportion is on the specified side of the ME (rounded to 4 dig)
#' @export
me_phat <- function(p, me, n, lower.tail = TRUE) {
  se <- sqrt(p * (1 - p) / n)

  z_upper <- z_score(p + me, mean = p, sd = se)
  z_lower <- z_score(p - me, mean = p, sd = se)

  if (lower.tail) {
    # Probability that phat is within margin (middle area)
    prob <- pnorm(z_upper, lower.tail = TRUE) - pnorm(z_lower, lower.tail = TRUE)
  } else {
    # Probability that phat is outside margin (two tails)
    prob <- pnorm(z_lower, lower.tail = TRUE) + pnorm(z_upper, lower.tail = FALSE)
  }

  return(round(prob, 4))
}


#' Calculate probability that sample mean is greater or less than a given value
#'
#' Uses the normal approximation based on the Central Limit Theorem to compute
#' P(sample mean > xbar) or P(sample mean < xbar)
#'
#' @param mu the population mean
#' @param sd the population standard deviation
#' @param xbar the sample mean value to test against
#' @param n the sample size
#' @param lower.tail logical; if TRUE, returns P(mean < xbar), otherwise P(mean > xbar)
#' @return the probability that the sample mean is on the specified tail, rounded to 4 digits
#' @export
prob_xbar <- function(mu, sd, xbar, n, lower.tail = TRUE) {
  se <- sd / sqrt(n)
  z <- z_score(xbar, mean = mu, sd = se)
  prob <- pnorm(z, lower.tail = lower.tail)
  return(round(prob, 4))
}


#' Calculate probability that sample mean is within or greater than a given margin of error
#'
#' Uses the normal approximation to compute P(mean > mu + me), P(mean < mu - me),
#' or the probability that the sample mean is within the margin.
#'
#' @param mu the population mean
#' @param sd the population standard deviation
#' @param me the margin of error
#' @param n the sample size
#' @param lower.tail logical; if TRUE, find the probability that the sample mean is within the ME,
#'                   otherwise probability that it is outside
#' @return the probability that the sample mean is on the specified side of the ME, rounded to 4 digits
#' @export
me_xbar <- function(mu, sd, me, n, lower.tail = TRUE) {
  se <- sd / sqrt(n)

  z_upper <- z_score(mu + me, mean = mu, sd = se)
  z_lower <- z_score(mu - me, mean = mu, sd = se)

  if (lower.tail) {
    # Probability that sample mean is within margin (middle area)
    prob <- pnorm(z_upper, lower.tail = TRUE) - pnorm(z_lower, lower.tail = TRUE)
  } else {
    # Probability that sample mean is outside margin (two tails)
    prob <- pnorm(z_lower, lower.tail = TRUE) + pnorm(z_upper, lower.tail = FALSE)
  }

  return(round(prob, 4))
}


#' cv_sd
#'
#' Computes the critical values from the Chi-Square distribution for constructing a 
#' confidence interval for the population variance or standard deviation.
#'
#'
#' @param n Integer. The sample size (must be greater than 1).
#' @param alpha Numeric. 1 - the confidence level for the interval (must be between 0 and 1).
#'
#' @return A numeric vector of length 2, containing the lower and upper critical values 
#' from the Chi-Square distribution, in that order: \code{c(chi_lower, chi_upper)}.
#'
#' @examples
#' # Get 90% confidence interval critical values with a sample size of 11
#' cv_sd(n = 11, conf_level = 0.90)
#'
#' # Get 95% confidence interval critical values with a sample size of 30
#' cv_sd(n = 30, conf_level = 0.95)
#'
#' @export
cv_sd <- function(n, alpha) {
    if (n <= 1) {
        stop("Sample size must be greater than 1.")
    }
    if (alpha <= 0 || alpha >= 1) {
        stop("alpha must be between 0 and 1.")
    }

    df <- n - 1

    # Chi-square critical values
    chi_upper <- qchisq(1 - alpha / 2, df)
    chi_lower <- qchisq(alpha / 2, df)

    return(c(chi_lower, chi_upper))
}


#' ci_sd
#'
#' Computes the confidence interval for the population standard deviation
#' based on a sample standard deviation, sample size, and confidence level.
#'
#' @param sd Numeric. The sample standard deviation.
#' @param n Integer. The sample size.
#' @param conf_level Numeric. The confidence level (e.g., 0.80, 0.95).
#'
#' @return A numeric vector of length 2 containing the lower and upper bounds of the confidence interval.
#' @examples
#' ci_sd(sd = 1.79, n = 83, conf_level = 0.80)
#' ci_sd(sd = 2.5, n = 30, conf_level = 0.95)
ci_sd <- function(sd, n, conf_level) {
  if (n <= 1) {
    stop("Sample size must be greater than 1.")
  }
  if (sd <= 0) {
    stop("Sample standard deviation must be greater than 0.")
  }
  if (conf_level <= 0 || conf_level >= 1) {
    stop("Confidence level must be between 0 and 1.")
  }

  df <- n - 1
  alpha <- 1 - conf_level

    # Critical values from the Chi-Square distribution
  chi_vals <- cv_sd(n, alpha)
  chi_lower <- chi_vals[2]
  chi_upper <- chi_vals[1]

  # Confidence interval for population standard deviation
  lower_bound <- sqrt((df * sd^2) / chi_lower)
  upper_bound <- sqrt((df * sd^2) / chi_upper)

  return(c(lower_bound, upper_bound))
}


#' ci_var
#'  
#' Computes the confidence interval for the population variance based on a sample variance,
#' sample size, and confidence level.
#' 
#' @param var Numeric. The sample variance.
#' @param n Integer. The sample size.   
#' @param conf_level Numeric. The confidence level (e.g., 0.80, 0.95).
#'  
#' @return A numeric vector of length 2 containing the lower and upper bounds of the confidence interval.
#' 
#' 
ci_var <- function (var, n, conf_level) {
  if (n <= 1) {
    stop("Sample size must be greater than 1.")
  }
  if (var <= 0) {
    stop("Sample variance must be greater than 0.")
  }
  if (conf_level <= 0 || conf_level >= 1) {
    stop("Confidence level must be between 0 and 1.")
  }

  df <- n - 1
  alpha <- 1 - conf_level

    # Critical values from the Chi-Square distribution
  chi_vals <- cv_sd(n, alpha)
  chi_lower <- chi_vals[1]
  chi_upper <- chi_vals[2]

  # Confidence interval for population variance
  lower_bound <- (df * var) / chi_upper
  upper_bound <- (df * var) / chi_lower

  return(c(lower_bound, upper_bound))
}


#' Calculate standard deviation for a binomial distribution DRV
#'
#'
#' @param n The number of trials
#' @param p the probability of success
#' @return The standard deviation of the DRV
#' @export
sd_binom <- function(x, prob) {

  sqrt(x*prob*(1-prob))

}


#' Calculate variance for a binomial distribution DRV
#'
#'
#' @param n The number of trials
#' @param p the probability of success
#' @return The variance of the DRV
#' @export
binom_var <- function(x, prob) {

  x*prob*(1-prob)
  
}


# Function to calculate the test statistic for two independent means with known population std deviations
#' @param x1 Mean of sample 1
#' @param x2 Mean of sample 2
#' @param s1 Population standard deviation for sample 1
#' @param s2 Population standard deviation for sample 2
#' @param n1 Sample size for sample 1
#' @param n2 Sample size for sample 2
#' @return z-test statistic for the difference in means
#'
t_stat_two_means_known_pop_sd <- function(x1, x2, s1, s2, n1, n2) {
  
  # Calculate z-test statistic
  z <- (x1 - x2) / sqrt((s1^2 / n1) + (s2^2 / n2))
  
  # Return result
  return(z)
}


# Function to calculate the standard error for two independent means with known population std deviations
#' @param x1 Mean of sample 1
#' @param x2 Mean of sample 2
#' @param s1 Population standard deviation for sample 1
#' @param s2 Population standard deviation for sample 2
#' @param n1 Sample size for sample 1
#' @param n2 Sample size for sample 2
#' @return Standard error for the difference between two means
#' 
se_two_means_known_pop_sd <- function(x1, x2, s1, s2, n1, n2) {

  se <- sqrt((s1^2 / n1) + (s2^2 / n2))
  return(se)
}



#' Compute a confidence interval for the difference between two means
#'
#' Assumes independent samples and known population standard deviations.
#'
#' @param x1 Mean of sample 1
#' @param x2 Mean of sample 2
#' @param s1 Population standard deviation for sample 1
#' @param s2 Population standard deviation for sample 2
#' @param n1 Sample size for sample 1
#' @param n2 Sample size for sample 2
#' @param conf_level Confidence level (e.g., 0.95 for 95%)
#'
#' @return A numeric vector with lower and upper bounds of the confidence interval
#' @examples
#' ci_two_means(65.9, 80, 12.47, 12.84, 157, 129, 0.95)
#' # Returns: c(-16.7, -10.9)
ci_two_means_known_pop_sd <- function(x1, x2, s1, s2, n1, n2, conf_level) {
  # Difference in means
  diff <- x1 - x2
  
  # Standard error
  se <- se_two_means_known_pop_sd(x1, x2, s1, s2, n1, n2)
  
  # z* critical value
  alpha <- 1 - conf_level
  z_star <- qnorm(1 - alpha / 2)
  
  # Confidence interval
  lower <- diff - z_star * se
  upper <- diff + z_star * se
  
  # Return rounded result
  return(c(lower, upper))
}



# Function to calculate the pooled standard deviation for two independent samples
#' @param x1 Mean of sample 1
#' @param x2 Mean of sample 2
#' @param s1 Sample standard deviation for sample 1
#' @param s2 Sample standard deviation for sample 2
#' @param n1 Sample size for sample 1
#' @param n2 Sample size for sample 2
#' @return Pooled standard deviation
#' 
sd_two_means_equal_var <- function (x1,x2,s1,s2,n1,n2){
    return (sqrt(((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2)))
}


#' Compute t-statistic for two-sample t-test with equal variances
#'
#' @param x1 Sample mean from group 1
#' @param x2 Sample mean from group 2
#' @param s1 Sample standard deviation from group 1
#' @param s2 Sample standard deviation from group 2
#' @param n1 Sample size from group 1
#' @param n2 Sample size from group 2
#'
#' @return t-stat
t_stat_two_means_equal_var <- function(x1, x2, s1, s2, n1, n2) {
  sd <- sd_two_means_equal_var(x1, x2, s1, s2, n1, n2)  # Pooled standard deviation
  t <- (x1 - x2) / (sd * sqrt(1/n1 + 1/n2))
  return(t)
}


# Function to calculate the standard error for the difference between two means with equal variances
#' @param x1 Mean of sample 1
#' @param x2 Mean of sample 2
#' @param s1 Sample standard deviation for sample 1
#' @param s2 Sample standard deviation for sample 2
#' @param n1 Sample size for sample 1
#' @param n2 Sample size for sample 2
#' @return Standard error for the difference between two means
#' 
se_two_means_equal_var <- function (x1,x2,s1,s2,n1,n2){
    sp <- sd_two_means_equal_var(x1,x2,s1,s2,n1,n2)
    return(sp * sqrt(1/n1 + 1/n2))
}


# Function to calculate the confidence interval for the difference between two means with equal variances
#' @param x1 Mean of sample 1
#' @param x2 Mean of sample 2
#' @param s1 Sample standard deviation for sample 1
#' @param s2 Sample standard deviation for sample 2
#' @param n1 Sample size for sample 1
#' @param n2 Sample size for sample 2
#' @param conf_level Confidence level (e.g., 0.95 for 95%)
#' @return A numeric vector with lower and upper bounds of the confidence interval
#' 
ci__two_means_equal_var <- function(x1, x2, s1, s2, n1, n2, conf_level) {
  # Point estimate
  point_estimate <- x1 - x2
  
  # Degrees of freedom
  df <- n1 + n2 - 2
  
  # Critical t value
  alpha <- 1 - conf_level
  t_crit <- qt(1 - alpha / 2, df)
  
  # Standard error
  se <- se_sp(x1, x2, s1, s2, n1, n2)
  
  # Confidence interval
  lower <- point_estimate - t_crit * se
  upper <- point_estimate + t_crit * se
  
  return(c(lower, upper))
}


# Function to calculate the pooled proportion for two independent samples
#' @param p1 Proportion from sample 1
#' @param p2 Proportion from sample 2
#' @param n1 Sample size of sample 1
#' @param n2 Sample size of sample 2
#' @return Pooled proportion
#' 
p_pooled <- function (p1,p2,n1,n2){
  return ((n1 * p1 + n2 * p2) / (n1 + n2))
}


# Function to calculate the test statistic for two independent proportions
#' @param p1 Proportion from sample 1
#' @param p2 Proportion from sample 2
#' @param n1 Sample size of sample 1
#' @param n2 Sample size of sample 2
#' @return Test statistic for the difference in proportions
#' 
t_stat_p_pooled <- function (p1,p2,n1,n2){
  if (any(c(p1, p2) < 0 | c(p1, p2) > 1)) stop("Proportions must be between 0 and 1")
  if (any(c(n1, n2) <= 0)) stop("Sample sizes must be positive")
  p_pool <- p_pooled(p1,p2,n1,n2)
  q = 1 - p_pool
  return ((p1 - p2) / sqrt((p_pool * q) / n1 + (p_pool * q) / n2))
}

# Function to calculate the standard error for the difference in proportions
#' @param p1 Proportion from sample 1 
#' @param p2 Proportion from sample 2
#' @param n1 Sample size of sample 1
#' @param n2 Sample size of sample 2
#' @return Standard error for the difference in proportions
#'
se_p_pooled <- function (p1,p2,n1,n2){
  p_pool <- p_pooled(p1,p2,n1,n2)
  q = 1 - p_pool
  return (sqrt((p_pool * q) / n1 + (p_pool * q) / n2))
}


# Function to calculate the confidence interval for the difference in proportions
#' @param p1 Proportion from sample 1
#' @param p2 Proportion from sample 2
#' @param n1 Sample size of sample 1
#' @param n2 Sample size of sample 2
#' @param conf_level Confidence level (e.g., 0.95 for 95% confidence)
#' @return Confidence interval for the difference in proportions
#' 
ci_p_pooled <- function (p1,p2,n1,n2,conf_level){
  p_pool <- p_pooled(p1,p2,n1,n2)
  q = 1 - p_pool
  z_crit <- qnorm(1 - (1 - conf_level) / 2)
  
  # Standard error for the difference in proportions
  se_diff <- p_pooled_se(p1, p2, n1, n2)
  
  # Confidence interval for the difference in proportions
  lower_bound <- (p1 - p2) - z_crit * se_diff
  upper_bound <- (p1 - p2) + z_crit * se_diff
  
  return(c(lower_bound, upper_bound))
}


# Function to calculate the test statistic for paired means
#' @param xbar_diff Difference in sample means (x1 - x2)
#' @param sd_diff Sample standard deviation of the differences
#' @param n Sample size (number of paired observations)
#' @return t-statistic for the paired means test
#'
t_stat_paired_means <- function(xbar_diff, sd_diff, n) {
  
  # Calculate t-statistic
  t_stat <- (xbar_diff) / (sd_diff / sqrt(n))
  
  return(t_stat)
}


# Function to calculate the standard error for paired means
#' @param sd_diff Sample standard deviation of the differences
#' @param n Sample size (number of paired observations)
#' @return Standard error for the paired means test
#' 
se_paired_means <- function(sd_diff, n) {
  
  # Calculate standard error for paired means
  se <- sd_diff / sqrt(n)
  
  return(se)
}


# Function to calculate the confidence interval for paired means
#' @param xbar_diff Difference in sample means (x1 - x2)
#' @param sd_diff Sample standard deviation of the differences
#' @param n Sample size (number of paired observations)
#' @param conf_level Confidence level (e.g., 0.95 for 95% confidence)
#' @return Confidence interval for the paired means test
#' 
ci_paired_means <- function(xbar_diff, sd_diff, n, conf_level) {
  
  # Calculate degrees of freedom
  df <- n - 1
  
  # Calculate critical t-value
  t_crit <- qt(1 - (1 - conf_level) / 2, df)
  
  # Calculate standard error
  se <- se_paired_means(sd_diff, n)
  
  # Calculate confidence interval
  lower_bound <- xbar_diff - t_crit * se
  upper_bound <- xbar_diff + t_crit * se
  
  return(c(lower_bound, upper_bound))
}


# function to calculate the standard error for two means with unknown population std deviations
#' @param x1 Mean of sample 1
#' @param x2 Mean of sample 2
#' @param s1 Sample standard deviation for sample 1
#' @param s2 Sample standard deviation for sample 2
#' @param n1 Sample size for sample 1
#' @param n2 Sample size for sample 2
#' @return Standard error for the difference between two means
#' 
se_two_means_unknown_pop_sd <- function (x1,x2,s1,s2,n1,n2){
    return (sqrt((s1^2/n1) + (s2^2/n2)))
}


# Function to calculate the t-statistic for two independent means with unknown population std deviations
#' @param x1 Mean of sample 1
#' @param x2 Mean of sample 2
#' @param s1 Sample standard deviation for sample 1
#' @param s2 Sample standard deviation for sample 2
#' @param n1 Sample size for sample 1
#' @param n2 Sample size for sample 2
#' @return t-statistic for the difference in means
#' 
t_stat_two_means_unknown_pop_sd <- function(x1, x2, s1, s2, n1, n2) {
  
  # Calculate t-statistic
  t <- (x1 - x2) / se_two_means_unknown_pop_sd(x1, x2, s1, s2, n1, n2)
  
  # Return result
  return(t)
}


# Function to calculate the confidence interval for the difference between two means with unknown population std deviations
#' @param x1 Mean of sample 1
#' @param x2 Mean of sample 2
#' @param s1 Sample standard deviation for sample 1
#' @param s2 Sample standard deviation for sample 2
#' @param n1 Sample size for sample 1
#' @param n2 Sample size for sample 2
#' @param conf_level Confidence level (e.g., 0.95 for 95%)
#' @return A numeric vector with lower and upper bounds of the confidence interval
#' 
ci_two_means_unknown_pop_sd <- function(x1, x2, s1, s2, n1, n2, conf_level) {
  # Point estimate
  point_estimate <- x1 - x2
  
  # Degrees of freedom
  df <- min(n1 - 1, n2 - 1)
  
  # Critical t value
  alpha <- 1 - conf_level
  t_crit <- qt(1 - alpha / 2, df)
  
  # Standard error
  se <- se_two_means_unknown_pop_sd(x1, x2, s1, s2, n1, n2)
  
  # Confidence interval
  lower <- point_estimate - t_crit * se
  upper <- point_estimate + t_crit * se
  
  return(c(lower, upper))
}