#' Calculate -log10(p-value) with high precision fallback
#'
#' Computes the `mlogp` column as `-log10(pval)` for a given data frame.
#' If `mlogp` is infinite due to extremely small p-values, it uses high-precision
#' computation based on `beta` and `se` columns to recalculate `mlogp` values.
#'
#' @param df A data frame containing a `pval` column. If any resulting `mlogp` values are `Inf`,
#' the `beta` and `se` columns must also be present for high-precision correction.
#'
#' @return A data frame identical to the input `df` but with an additional `mlogp` column.
#'
#' @importFrom Rmpfr mpfr asNumeric
#' @export
sg.cal_mlogp = function(df){

  if (!'pval' %in% colnames(df)){
    stop("Error, `pval` not in df's columns")
  }
  
  df$mlogp = -log10(df$pval)
  
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
      calculated_mlogp = Rmpfr::asNumeric(
        round(-log10(
          2 * (1 - pnorm(abs(Rmpfr::mpfr(inf_mlogp_df$beta, prec) / Rmpfr::mpfr(inf_mlogp_df$se, prec))))
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

