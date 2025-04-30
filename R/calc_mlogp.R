# calc_mlogp.R
# This script contains function to calculate mlogp for data.frame

library(Rmpfr)

cal_mlogp = function(df){
  # input:
  #   - df: a data.frame, need to have `pval` column
  #         `beta` and `se` column is required for Inf mlogp calculated
  # output:
  #   - df: the same df with mlogp calculated
  
  
  if (!'pval' %in% colnames(df)){
    stop("Error, `pval` not in df's columns")
  }
  
  df$mlogp = -log10(df$pva)
  
  # get rows with Inf mlogp
  inf_mlogp_df = df[is.infinite(df$mlogp), ]
  
  # use high precision calculation for Inf mlogp
  if (nrow(inf_mlogp_df) > 0){
    
    # check if beta and se is provided in df
    if (!'beta' %in% colnames(df)){
      stop("Error, Inf mlogp detected and `beta` is not supplied to calculate")
    }
    if (!'se' %in% colnames(df)){
      stop("Error, Inf mlogp detected and `se` is not supplied to calculate")
    }
    
    # set precision
    prec = 1000
    
    #  gradually increase prec
    while (any(is.infinite(inf_mlogp_df$mlogp))){
      calculated_mlogp = asNumeric(
        round(-log10(
          2 * (1 - pnorm(abs(mpfr(inf_mlogp_df$beta, prec) / mpfr(inf_mlogp_df$se, prec))))
        ), 5)
      )
      inf_mlogp_df$mlogp = calculated_mlogp
      prec = prec + 500
    }
    
    # update Inf mlogp
    df$mlogp[is.infinite(df$mlogp)] = inf_mlogp_df$mlogp
  }
  
  return(df)
}

