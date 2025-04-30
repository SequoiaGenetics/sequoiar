# -----------------------------------------------------------------------------
# File Name: meta_analysis.R
# Author: Lin
# Date: 17/03/2025
# Description: 
#
# Input: 
#   - abcdefg
#
# Output: 
#   - abcdefg
#
# Key Steps:
# 1. abcdefg
#
# Required Libraries:
# - meta
#
# Notes:
# -----------------------------------------------------------------------------

library(meta)

# Example vectors (unnamed, with values being c(beta, se))
s1 = c(0.5, 0.1)
s2 = c(0.8, 0.2)
s3 = c(-0.3, 0.15)

# Combine into a list
s_list = list(s1, s2, s3)

# Extract beta (first element) and se (second element) into separate vectors
beta = sapply(s_list, function(x) x[1])
se = sapply(s_list, function(x) x[2])

print(beta)  # c(0.5, 0.8, -0.3)
print(se)    # c(0.1, 0.2, 0.15)

# Meta analysis
meta_res = metagen(beta, se)

# Print results
print(meta_res)

# Example for getting estimates
meta_res$TE.common # common effect beta
meta_res$seTE.common # common effect se
meta_res$zval.fixed # fixed effect Z score (the same as common (?))
meta_res$lower.random # random effect confidence interval lower bound
meta_res$pval.fixed # fixed effect p value of beta

meta_res$pval.Q # heterogeneity
meta_res$df.Q # Q statistics degree of freedom
meta_res$Q # Q statistics
