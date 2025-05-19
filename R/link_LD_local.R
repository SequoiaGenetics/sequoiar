#' Append LD values to SNPs using local PLINK via `ieugwasr`
#'
#' This function uses `ieugwasr::ld_matrix_local()` to compute LD (`r²`) between a specified lead SNP
#' and other SNPs in a given locus. It appends the LD values to the `data` frame inside the `loc_data` list.
#' It assumes a PLINK binary and EUR reference panel are available in the installed package directory.
#'
#' @param loc_data A locus object when making locus plot containing a `data` data frame with a `SNP` column.
#' @param lead_SNP A character string specifying the SNP to be used as the LD reference.
#'                 Must exist in `loc_data$data$SNP`.
#'
#' @return A modified version of `loc_data` with an added `ld` column in `loc_data$data`.
#'
#' @importFrom dplyr left_join
#' @export
sg.link_LD_local = function(loc_data, lead_SNP) {
  # Validate input
  if (!"SNP" %in% colnames(loc_data$data)) {
    stop("loc_data$data must contain a 'SNP' column.")
  }
  if (!(lead_SNP %in% loc_data$data$SNP)) {
    stop("lead_SNP must be present in loc_data$data$SNP.")
  }
  
  # Resolve internal PLINK paths
  plink_dir = file.path(system.file(package = "sequoiar"), "plink")
  plink_bin = file.path(plink_dir, ifelse(.Platform$OS.type == "windows", "plink.exe", "plink"))
  bfile = file.path(plink_dir, "EUR")
  
  if (!file.exists(plink_bin)) {
    stop("PLINK binary not found in package plink directory. Please run `sg._copy_plink_win()` first.")
  }
  if (!all(file.exists(paste0(bfile, c(".bed", ".bim", ".fam"))))) {
    stop("EUR reference panel files not found in package plink directory.")
  }
  
  # Compute LD matrix using ieugwasr
  ld_mat = ieugwasr::ld_matrix_local(
    snps = loc_data$data$SNP,
    with_alleles = FALSE,
    bfile = bfile,
    plink_bin = plink_bin
  )
  
  # Create LD vector relative to lead SNP
  ld_df = data.frame(
    SNP = colnames(ld_mat),
    ld = as.vector(ld_mat[, lead_SNP])
  )
  
  # Merge LD info back to the original data
  loc_data$data = dplyr::left_join(loc_data$data, ld_df, by = "SNP")
  
  return(loc_data)
}