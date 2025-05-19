#' Harmonise exposure data with one or more outcome datasets
#'
#' This function standardises and aligns multiple outcome datasets to a common exposure dataset.
#' It handles INDEL filtering, allele ordering, and sign correction of `beta` values based on
#' matching effect alleles. All input data frames must contain `SNP`, `effect_allele`, and `other_allele` columns.
#'
#' @param exp A data frame representing the exposure GWAS data.
#' @param ... One or more outcome data frames to harmonise to the exposure data.
#'
#' @return A harmonised data frame containing exposure and outcome data with suffixes
#' corresponding to their variable names. Columns are aligned by `SNP`, alleles, and optional `chr`, `pos`, or `pos19`.
#'
#' @importFrom dplyr filter mutate rename_with left_join select
#' @export
sg.harmonise = function(exp, ...){
  # Convert additional arguments to a list of dataframes
  dataframes = list(...)
  
  # Get the name of the exp dataframe
  exp_name = deparse(substitute(exp))
  # get the name of other dataframe
  out_names = sapply(substitute(list(...))[-1], deparse)
  
  # Make sure there is no INDEL SNPs  
  exp = dplyr::filter(exp, nchar(exp$effect_allele) == 1 & nchar(exp$other_allele) == 1)
  exp$effect_allele = toupper(exp$effect_allele)
  exp$other_allele = toupper(exp$other_allele)
  exp = dplyr::filter(exp, !exp$effect_allele %in% c("I", "D") & !exp$other_allele %in% c("I", "D"))
  exp$allele_concat = paste(pmin(exp$effect_allele, exp$other_allele),
                             pmax(exp$effect_allele, exp$other_allele), sep = "")
  exp = dplyr::rename_with(exp, ~ paste0(., "_", exp_name), -c("SNP", "allele_concat", "effect_allele", "other_allele"))
  
  # Harmonise each dataframe to the exp dataframe
  harmonised_df = exp
  for (i in seq_along(dataframes)) {
    out = dataframes[[i]]
    out_name = out_names[i]
    suffix_out = paste0('_', out_name)
    
    # Make sure there is no INDEL SNPs
    out = dplyr::filter(out, nchar(out$effect_allele) == 1 & nchar(out$other_allele) == 1)
    out$effect_allele = toupper(out$effect_allele)
    out$other_allele = toupper(out$other_allele)
    out = dplyr::filter(out, !out$effect_allele %in% c("I", "D") & !out$other_allele %in% c("I", "D"))
    out$allele_concat = paste(pmin(out$effect_allele, out$other_allele),
                               pmax(out$effect_allele, out$other_allele), sep = "")
    out = dplyr::rename_with(out, ~ paste0(., suffix_out), -c("SNP", "allele_concat"))
    
    # Perform the join on SNP and allele_concat
    by_cols = c('SNP', 'allele_concat')
    
    # if chr in both exp and out, join using chr as well
    if ('chr' %in% colnames(exp) && 'chr' %in% colnames(out)){
      by_cols = c(by_cols, 'chr')
    }
    # if pos in both exp and out, join using pos as well
    if ('pos' %in% colnames(exp) && 'pos' %in% colnames(out)){
      by_cols = c(by_cols, 'pos')
    } else if ('pos19' %in% colnames(exp) && 'pos19' %in% colnames(out)){
      by_cols = c(by_cols, 'pos19')
    }
    
    
    # Perform the join on SNP and allele_concat, and add suffixes to distinguish columns
    harmonised_df = dplyr::left_join(harmonised_df, out, by = by_cols)

    harmonised_df$effect_allele_out_orig = harmonised_df[[paste0("effect_allele", suffix_out)]]
    harmonised_df$other_allele_out_orig = harmonised_df[[paste0("other_allele", suffix_out)]]
    
    harmonised_df[[paste0("beta", suffix_out)]] =
      ifelse(harmonised_df[[paste0("effect_allele_", exp_name)]] == harmonised_df$effect_allele_out_orig,
             harmonised_df[[paste0("beta", suffix_out)]],
             -1 * harmonised_df[[paste0("beta", suffix_out)]])
    
    harmonised_df[[paste0("effect_allele", suffix_out)]] =
      ifelse(harmonised_df[[paste0("effect_allele_", exp_name)]] == harmonised_df$effect_allele_out_orig,
             harmonised_df$effect_allele_out_orig,
             harmonised_df$other_allele_out_orig)
    
    harmonised_df[[paste0("other_allele", suffix_out)]] =
      ifelse(harmonised_df[[paste0("effect_allele_", exp_name)]] == harmonised_df$effect_allele_out_orig,
             harmonised_df$other_allele_out_orig,
             harmonised_df$effect_allele_out_orig)
    
    harmonised_df$effect_allele_out_orig = NULL
    harmonised_df$other_allele_out_orig = NULL
  }
  
  # Rename effect_allele and other_allele with suffix of exposure
  harmonised_df = dplyr::rename_with(harmonised_df, ~ paste0(., "_", exp_name),
                                      c("effect_allele", "other_allele"))
  harmonised_df$allele_concat = NULL
  
  return(harmonised_df)
}
