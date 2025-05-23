% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwas_stats_utils.R
\name{sg.complete_stats}
\alias{sg.complete_stats}
\title{Complete missing GWAS summary statistics}
\usage{
sg.complete_stats()
}
\value{
A function that accepts the arguments:
\describe{
  \item{beta}{Numeric value for effect size.}
  \item{se}{Numeric value for standard error.}
  \item{p_value}{Numeric value for p-value.}
  \item{method}{One of `"chisq"` (default) or `"normal"`, specifying the statistical method to use.}
}
and returns a named numeric vector of the computed missing value.
}
\description{
This function returns another function that completes one missing GWAS statistic (`beta`, `se`, or `p_value`)
from the other two using either the chi-squared or normal distribution. Exactly one of the three inputs must be `NA`.
}
\examples{
fn <- sg.complete_stats()
fn(beta = 0.2, se = 0.05, p_value = NA)     # Calculate p-value
fn(beta = NA, se = 0.05, p_value = 1e-6)    # Calculate beta
fn(beta = 0.1, se = NA, p_value = 0.01)     # Calculate se

}
