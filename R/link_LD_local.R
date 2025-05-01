#' @importFrom magrittr %>%
NULL

link_LD_local = function(loc_data){
  ld_mat = ld_matrix_local(loc_data$data$SNP,with_alleles=F,bfile='util/plink/EUR',plink_bin='util/plink/plink')
  ld_df = data.frame(
    SNP = colnames(ld_mat),
    ld = as.vector(ld_mat[, lead_SNP])# Ensure values are numeric
  )
  loc_data$data = loc_data$data %>%
    left_join(ld_df, by = 'SNP')
  return(loc_data)
}