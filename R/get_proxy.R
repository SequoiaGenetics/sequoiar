#' Find Proxy SNPs Using PLINK and EUR Reference Panel
#'
#' This function identifies proxy SNPs in high linkage disequilibrium (LD) with a target SNP
#' using PLINK and a EUR reference panel. It relies on PLINK's `--r2` command with a
#' user-defined r² threshold. The reference panel and PLINK binary are expected to be in the
#' `plink/` subdirectory of the `sequoiar` package.
#'
#' @param rsid Character. The rsID of the target SNP for which to find proxies.
#' @param r2_threshold Numeric. Minimum r² value to report SNPs as proxies. Default is 0.5.
#' @param bfile Optional. Path to PLINK binary file prefix (without extension). If not provided, 
#'   the function uses the default EUR reference panel included in the package.
#'
#' @return A data frame of proxy SNPs with columns: `SNP` (proxy SNP ID), `chr`, `pos`, and `r2`.
#'   Returns an informative error if the target SNP is not found in the reference panel.
#'
#' @examples
#' \dontrun{
#' # Find proxies for rs146597587 with default EUR reference and r² ≥ 0.5
#' proxies <- sg.get_proxy("rs146597587")
#' 
#' # Use a custom reference panel and a higher LD threshold
#' proxies <- sg.get_proxy("rs123456", r2_threshold = 0.8, bfile = "path/to/custom_panel")
#' }
#'
#' @export
sg.get_proxy = function(rsid, r2_threshold = 0.5, bfile = NULL){
    # path to plink and reference panel
  plink_dir = file.path(system.file(package = "sequoiar"), "plink")
  plink_bin = file.path(plink_dir, ifelse(.Platform$OS.type == "windows", "plink.exe", "plink"))
  
  if (is.null(bfile)){
    message("bfile not provided, using default 1000G EUR as reference...")
    bfile = file.path(plink_dir, "EUR")  # Exclude file extension to let PLINK find .bed/.bim/.fam
  }
  
  message(paste0('Finding proxy SNPs with '))
  
  if (!file.exists(plink_bin)) {
    stop("PLINK binary not found in package plink directory. Please run `sg._copy_plink_win()` first.")
  }
  if (!all(file.exists(paste0(bfile, c(".bed", ".bim", ".fam"))))) {
    stop("EUR reference panel files not found in package plink directory.")
  }
  
  # Generate random file prefix in tempdir
  fn_prefix = file.path(tempdir(), paste0(sample(LETTERS, 10, replace = TRUE), collapse = ""))

  # Construct PLINK command
  shell_type = ifelse(.Platform$OS.type == "windows", "cmd", "sh")
  plink_cmd = paste0(
    shQuote(plink_bin, type = shell_type),
    " --bfile ", shQuote(bfile, type = shell_type),
    " --ld-snp ", rsid,
    " --ld-window ", 500*1000,
    " --ld-window-r2 ", r2_threshold,
    " --r2",
    " --out ", shQuote(fn_prefix, type = shell_type)
  )
  
  
  message(paste0("⏳  Getting Proxy SNPs above R2 threshold ", r2_threshold, ' ...'))

    # Run PLINK and capture output
  plink_output <- tryCatch(
    system(plink_cmd, intern = TRUE),
    error = function(e) return(character(0))
  )
  
  # Check for known PLINK error
  if (any(grepl("Error: No valid variants specified", plink_output))) {
    stop(paste0("❌ SNP '", rsid, "' not found in reference panel."))
  }
  
  # Check output file exists
  proxy_file <- paste0(fn_prefix, ".ld")
  if (!file.exists(proxy_file)) {
    stop("❌ PLINK failed: `.ld` file not found.")
  }
  
  # Read and process results
  res <- read.table(proxy_file, header = TRUE)
  unlink(paste0(fn_prefix, "*"))  # Clean temp files
  
  res <- dplyr::select(res, SNP = SNP_B, chr = CHR_B, pos = BP_B, r2 = R2)
  res <- dplyr::arrange(res, desc(r2))
  return(res)
}