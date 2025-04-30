# gwas_stats_utils.R
# This script contains function that calculate one of beta/se/pval from the other 2
# This script provides two methods for calculation, chisq (default) or normal
# `complete_stats` is the wrapper function for the calculation

# Function to complete gwas stats
complete_stats = function(){
  # input:
  #   - beta, se, pval should be vectors
  #   - method should be either "chisq" or "normal"
  function(beta = NA, se = NA, p_value = NA, method = "chisq") {
    
    # function(beta = NA, se = NA, p_value = NA, method = "chisq") {
    # Ensure exactly one value is missing
    missing_count <- sum(is.na(c(beta, se, p_value)))
    if (missing_count != 1) {
      stop("Exactly one of beta, se, or p_value must be NA.")
    }
    
    if (method == "normal") {
      if (is.na(p_value)) {
        return(c(p_value = calculate_pvalue_normal(beta, se)))
      } else if (is.na(beta)) {
        return(c(beta = calculate_beta_normal(se, p_value)))
      } else if (is.na(se)) {
        return(c(se = calculate_se_normal(beta, p_value)))
      }
    } else {  # Default to Chi-square
      if (is.na(p_value)) {
        return(c(p_value = calculate_pvalue_chisq(beta, se)))
      } else if (is.na(beta)) {
        return(c(beta = calculate_beta_chisq(se, p_value)))
      } else if (is.na(se)) {
        return(c(se = calculate_se_chisq(beta, p_value)))
      }
    }
  }
}

# Function to calculate p-value using the Chi-square distribution (df = 1)
calculate_pvalue_chisq = function(beta, se) {
  z_score = beta / se
  chi_square_stat = z_score^2
  p_value = pchisq(chi_square_stat, df = 1, lower.tail = FALSE)
  return(p_value)
}

# Function to calculate p-value using the normal distribution (two-tailed)
calculate_pvalue_normal = function(beta, se) {
  z_score = beta / se
  p_value = 2 * pnorm(-abs(z_score))  # Two-tailed p-value
  return(p_value)
}

# Function to calculate beta from standard error and p-value using Normal distribution
calculate_beta_normal = function(se, p_value) {
  z_score = qnorm(p_value / 2, lower.tail = FALSE)  # Inverse normal CDF
  beta = se * z_score
  return(beta)
}

# Function to calculate beta from standard error and p-value using Chi-square distribution
calculate_beta_chisq = function(se, p_value) {
  chi_square_stat = qchisq(p_value, df = 1, lower.tail = FALSE)  # Inverse chi-square CDF
  beta = se * sqrt(chi_square_stat) * sign(beta)  # Recover sign (if needed)
  return(beta)
}

# Function to calculate standard error from beta and p-value using Normal distribution
calculate_se_normal = function(beta, p_value) {
  z_score = qnorm(p_value / 2, lower.tail = FALSE)  # Inverse normal CDF
  se = abs(beta / z_score)
  return(se)
}

# Function to calculate standard error from beta and p-value using Chi-square distribution
calculate_se_chisq = function(beta, p_value) {
  chi_square_stat = qchisq(p_value, df = 1, lower.tail = FALSE)  # Inverse chi-square CDF
  se = abs(beta) / sqrt(chi_square_stat)
  return(se)
}


