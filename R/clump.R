#' Local LD clumping using PLINK and EUR reference panel
#'
#' This function performs LD clumping on a data frame using PLINK, based on SNP-level
#' `-log10(p-value)` significance and a provided LD threshold (`r2`). It expects PLINK
#' binary and a EUR reference panel to be placed in the package's `plink/` directory.
#'
#' @param df A data frame containing at least `SNP` and `mlogp` columns. `mlogp` must not contain Inf.
#' @param clump_r2 Numeric. The r² threshold for LD clumping (e.g., 0.1).
#'
#' @return A data frame of SNPs that remain after clumping, with the same structure as input `df`.
#' 
#' @export
sg.clump = function(df, clump_r2){

  # init constant
  clump_kb = 10000
  clump_p1 = 1
  clump_p2 = 1
  pval_column = "mlogp"
  
  # path to plink and reference panel
  plink_dir = file.path(system.file(package = "sequoiar"), "plink")
  plink_bin = file.path(plink_dir, ifelse(.Platform$OS.type == "windows", "plink.exe", "plink"))
  bfile = file.path(plink_dir, "EUR")  # Exclude file extension to let PLINK find .bed/.bim/.fam
  
  if (!file.exists(plink_bin)) {
    stop("PLINK binary not found in package plink directory. Please run `sg._copy_plink_win()` first.")
  }
  if (!all(file.exists(paste0(bfile, c(".bed", ".bim", ".fam"))))) {
    stop("EUR reference panel files not found in package plink directory.")
  }
  
  if (!pval_column %in% colnames(df)) {
    stop("Error: dataframe does not contain `mlogp` column for -log10(p-value).")
  }
  if (any(is.infinite(df$mlogp))) {
    stop("Error: `mlogp` column contains Inf values.")
  }
  
  # filter out rows with pval > 0.05 in case haven't
  sig_df = df[df$mlogp > -log10(0.05), ]
  
  # filter out rows without rsID and create input df
  input_df = sig_df[grepl("^rs", sig_df$SNP), c("SNP", "mlogp")]
  
  if (nrow(input_df) == 0) {
    stop("Clumping stopped: no significant SNPs with rsIDs.")
  }
  
  if (nrow(input_df) == 1) {
    return(subset(sig_df, sig_df$SNP == input_df$SNP))
  }
  

  # Generate random file prefix in tempdir
  fn_prefix = file.path(tempdir(), paste0(sample(LETTERS, 10, replace = TRUE), collapse = ""))
  clump_input_file = fn_prefix
  clump_output_prefix = fn_prefix

  # Write clump input file (PLINK expects p-values, so use 1/mlogp as proxy ranking)
  write.table(
    data.frame(SNP = input_df$SNP, P = 1 / input_df$mlogp),
    file = clump_input_file,
    row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  
  # Construct PLINK command
  shell_type = ifelse(.Platform$OS.type == "windows", "cmd", "sh")
  plink_cmd = paste0(
    shQuote(plink_bin, type = shell_type),
    " --bfile ", shQuote(bfile, type = shell_type),
    " --clump ", shQuote(clump_input_file, type = shell_type),
    " --clump-p1 ", clump_p1,
    " --clump-p2 ", clump_p2,
    " --clump-r2 ", clump_r2,
    " --clump-kb ", clump_kb,
    " --out ", shQuote(clump_output_prefix, type = shell_type)
  )
  

  message("⏳ Running PLINK clumping...")
  system(plink_cmd, intern = TRUE)
  
  clumped_file = paste0(clump_output_prefix, ".clumped")
  if (!file.exists(clumped_file)) {
    stop("PLINK clumping failed: `.clumped` file not found.")
  }
  
  # Read clumped results
  res = read.table(clumped_file, header = TRUE)
  unlink(paste0(fn_prefix, "*"))  # clean up temp files
  
  # Report removal
  removed_snps = setdiff(input_df$SNP, res$SNP)
  if (length(removed_snps) > 0) {
    message("Removing ", length(removed_snps), " of ", nrow(input_df),
            " SNPs due to LD or absence from reference panel.")
  }
  
  # Return clumped df
  output_df = subset(sig_df, SNP %in% res$SNP)
  return(output_df)
}



#' Copy PLINK binary and reference panel into package directory (Windows)
#'
#' This function copies all files (including `plink.exe`, `.bed`, `.bim`, `.fam`, etc.)
#' from a specified local directory into the internal `plink/` folder inside the installed
#' package directory (i.e., where `system.file(package = "sequoiar")` points to).
#'
#' @param local_plink_dir A string. The full path to the local folder containing `plink.exe`
#' and any reference panel files (e.g., `.bed`, `.bim`, `.fam`) you wish to use.
#'
#' @return Invisibly returns the path to the destination folder where files were copied.
#' 
#' @examples
#' \dontrun{
#' sg._copy_plink_win("V:/resource/tools/plink")
#' 
#' @export
sg._copy_plink_win = function(local_plink_dir) {
  # Destination: plink folder inside installed package directory
  plink_dir = file.path(system.file(package = "sequoiar"), "plink")
  
  if (!dir.exists(plink_dir)) dir.create(plink_dir, recursive = TRUE)
  
  message("Copying PLINK and reference files from: ", local_plink_dir)
  message("To: ", plink_dir)
  
  # Get list of all files (recursively) from source
  files_to_copy = list.files(local_plink_dir, recursive = TRUE, full.names = TRUE)
  
  # Destination paths preserving directory structure
  rel_paths = list.files(local_plink_dir, recursive = TRUE)
  dest_paths = file.path(plink_dir, rel_paths)
  
  # Create necessary subdirectories
  unique_dirs = unique(dirname(dest_paths))
  for (d in unique_dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  
  # Copy files
  success = file.copy(
    from = files_to_copy,
    to = dest_paths,
    overwrite = TRUE
  )
  
  if (all(success)) {
    message("✔ All files successfully copied.")
  } else {
    warning("⚠ Some files could not be copied.")
  }
  
  invisible(plink_dir)
}