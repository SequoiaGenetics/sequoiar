#' Complete missing GWAS summary statistics
#'
#' This function returns another function that completes one missing GWAS statistic (`beta`, `se`, or `p_value`)
#' from the other two using either the chi-squared or normal distribution. Exactly one of the three inputs must be `NA`.
#'
#' @return A function that accepts the arguments:
#' \describe{
#'   \item{beta}{Numeric value for effect size.}
#'   \item{se}{Numeric value for standard error.}
#'   \item{p_value}{Numeric value for p-value.}
#'   \item{method}{One of `"chisq"` (default) or `"normal"`, specifying the statistical method to use.}
#' }
#' and returns a named numeric vector of the computed missing value.
#'
#' @examples
#' fn <- sg.complete_stats()
#' fn(beta = 0.2, se = 0.05, p_value = NA)     # Calculate p-value
#' fn(beta = NA, se = 0.05, p_value = 1e-6)    # Calculate beta
#' fn(beta = 0.1, se = NA, p_value = 0.01)     # Calculate se
#'
#' @export
sg.complete_stats = function(beta = NA, se = NA, p_value = NA, method = "chisq") {
  
  # function(beta = NA, se = NA, p_value = NA, method = "chisq") {
  # Ensure exactly one value is missing
  missing_count = sum(is.na(c(beta, se, p_value)))
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


#'  Function to calculate p-value using the Chi-square distribution (df = 1)
#' @keywords internal
calculate_pvalue_chisq = function(beta, se) {
  z_score = beta / se
  chi_square_stat = z_score^2
  p_value = pchisq(chi_square_stat, df = 1, lower.tail = FALSE)
  return(p_value)
}

#'  Function to calculate p-value using the normal distribution (two-tailed)
#' @keywords internal
calculate_pvalue_normal = function(beta, se) {
  z_score = beta / se
  p_value = 2 * pnorm(-abs(z_score))  # Two-tailed p-value
  return(p_value)
}

#'  Function to calculate beta from standard error and p-value using Normal distribution
#' @keywords internal
calculate_beta_normal = function(se, p_value) {
  z_score = qnorm(p_value / 2, lower.tail = FALSE)  # Inverse normal CDF
  beta = se * z_score
  return(beta)
}

#'  Function to calculate beta from standard error and p-value using Chi-square distribution
#' @keywords internal
calculate_beta_chisq = function(se, p_value) {
  chi_square_stat = qchisq(p_value, df = 1, lower.tail = FALSE)  # Inverse chi-square CDF
  beta = se * sqrt(chi_square_stat) * sign(beta)  # Recover sign (if needed)
  return(beta)
}

#'  Function to calculate standard error from beta and p-value using Normal distribution
#' @keywords internal
calculate_se_normal = function(beta, p_value) {
  z_score = qnorm(p_value / 2, lower.tail = FALSE)  # Inverse normal CDF
  se = abs(beta / z_score)
  return(se)
}

# Function to calculate standard error from beta and p-value using Chi-square distribution
#' @keywords internal
calculate_se_chisq = function(beta, p_value) {
  chi_square_stat = qchisq(p_value, df = 1, lower.tail = FALSE)  # Inverse chi-square CDF
  se = abs(beta) / sqrt(chi_square_stat)
  return(se)
}


