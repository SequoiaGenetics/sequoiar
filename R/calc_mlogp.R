#' Calculate -log10(p-value) with high precision fallback
#'
#' Computes the `mlogp` column as `-log10(pval)` for a given data frame.
#' If any calculated `mlogp` values are `Inf` (often due to extremely small p-values
#' underflowing to 0 in standard double precision), this function attempts to
#' recalculate these `mlogp` values using high-precision arithmetic via the `Rmpfr`
#' package, provided that `beta` and `se` columns are available in the input data frame.
#'
#' @param df A data frame containing a `pval` column. If any resulting `mlogp`
#'   values are `Inf`, the `beta` and `se` columns must also be present in `df`
#'   for high-precision correction.
#' @param start_prec Numeric. Initial precision in bits for `Rmpfr` calculations.
#'   Default is 256.
#' @param max_prec Numeric. Maximum precision in bits for `Rmpfr` calculations.
#'   The function will not attempt precision beyond this value. Default is 4096.
#' @param prec_step Numeric. Increment in bits for precision in each iteration if
#'   `Inf` values persist. Default is 256.
#' @param max_iter Integer. Maximum number of iterations for increasing precision.
#'   Default is 10.
#' @param verbose Logical. If `TRUE`, prints messages about the high-precision
#'   calculation process. Default is `TRUE`.
#'
#' @return A data frame identical to the input `df` but with an additional
#'   `mlogp` column. If `mlogp` could not be resolved from `Inf` even with
#'   high precision (e.g., `beta` and `se` not provided, `Rmpfr` not installed,
#'   or p-value is truly zero even at `max_prec`), it will remain `Inf`.
#'   Values that were `NA` in `pval` will result in `NA` in `mlogp`.
#'
#' @importFrom Rmpfr mpfr asNumeric pnorm
#' @export
#'
#' @examples
#' # Create a sample data frame
#' df_test <- data.frame(
#'   id = 1:5,
#'   # pval for id=3 is small but non-zero, id=4,5 will be 0 in double precision
#'   pval = c(0.05, 1e-50, 1e-300, 1e-310, 1e-400),
#'   beta = c(0.5, 1.0, 10, 15, 20),
#'   se = c(0.1, 0.1, 0.1, 0.1, 0.1) # For Z=100, Z=150, Z=200
#' )
#' # For 1e-310, -log10() is 310. For 1e-400, pval becomes 0, -log10(0) is Inf.
#'
#' # Add a case with SE = 0 (should remain Inf or handle as per interpretation)
#' df_test_se0 <- data.frame(
#'   id = 6:7,
#'   pval = c(0, 0), # pval might be pre-set to 0 if SE was 0
#'   beta = c(1, 0),
#'   se = c(0, 0.005) # SE=0 for id=6, normal for id=7 (Z=0)
#' )
#' df_test <- rbind(df_test, df_test_se0)
#'
#' # Add a case where pval is NA
#' df_test_na <- data.frame(id=8, pval=NA, beta=1, se=0.1)
#' df_test <- rbind(df_test, df_test_na)
#'
#' print("Original df_test:")
#' print(df_test)
#'
#' # Calculate mlogp
#' # Ensure Rmpfr is installed for this example to run fully for small p-values
#' if (requireNamespace("Rmpfr", quietly = TRUE)) {
#'   df_with_mlogp <- sg.cal_mlogp(df_test, verbose = TRUE, max_iter = 3) # Limit iter for example speed
#'   print("df_test with mlogp:")
#'   print(df_with_mlogp)
#'
#'   # Example with only pval, no beta/se, and one Inf
#'   df_simple_inf <- data.frame(pval = c(0.01, 0))
#'   # This will attempt high-precision, then error if beta/se are missing.
#'   message("Example with Inf pval but no beta/se (expect error):")
#'   tryCatch({
#'     sg.cal_mlogp(df_simple_inf)
#'   }, error = function(e) {
#'     message("Caught expected error: ", e$message)
#'     # Show what happens without Rmpfr (manually)
#'     df_simple_inf$mlogp_manual <- -log10(df_simple_inf$pval)
#'     print(df_simple_inf)
#'   })
#'
#' } else {
#'   message("Rmpfr not installed. High-precision fallback will not be fully demonstrated.")
#'   # Calculate mlogp without Rmpfr fallback (will show Inf for pval=0)
#'   df_test$mlogp_no_rmpfr <- -log10(df_test$pval)
#'   print("df_test with mlogp (Rmpfr not available):")
#'   print(df_test)
#' }
#'
#' # Example of a p-value that doesn't need Rmpfr at all
#' df_no_inf <- data.frame(pval = c(0.1, 0.01, 0.0000000000000000000000000000000000000001))
#' # The last pval is 1e-40, -log10(1e-40) = 40. No Inf.
#' print("Example with no Inf needing correction:")
#' print(sg.cal_mlogp(df_no_inf, verbose = FALSE))
#'
sg.cal_mlogp = function(df, start_prec = 256, max_prec = 4096, prec_step = 256, max_iter = 10, verbose = TRUE){
  
  if (!'pval' %in% colnames(df)){
    stop("Error, `pval` not in df's columns")
  }
  
  df$mlogp = -log10(df$pval)
  
  # Identify rows with Inf mlogp
  inf_mlogp_indices = which(is.infinite(df$mlogp))
  
  # Use high precision calculation for Inf mlogp
  if (length(inf_mlogp_indices) > 0){
    if (verbose) {
      message(paste(length(inf_mlogp_indices), "rows with Inf mlogp found. Attempting high-precision calculation."))
    }
    
    # Check if Rmpfr is installed
    if (!requireNamespace("Rmpfr", quietly = TRUE)) {
      stop("Package 'Rmpfr' is needed for high-precision calculations but is not installed. Please install it via install.packages('Rmpfr').")
    }
    
    # Check if beta and se are provided in df
    if (!'beta' %in% colnames(df)){
      stop("Error, Inf mlogp detected and `beta` is not supplied to calculate.")
    }
    if (!'se' %in% colnames(df)){
      stop("Error, Inf mlogp detected and `se` is not supplied to calculate.")
    }
    
    subset_to_fix = df[inf_mlogp_indices, ]
    
    current_prec_bits = start_prec
    iter_count = 0
    
    while (any(is.infinite(subset_to_fix$mlogp)) && current_prec_bits <= max_prec && iter_count < max_iter) {
      iter_count <- iter_count + 1
      if (verbose) {
        message(paste("Iteration", iter_count, ": Attempting high-precision with precBits =", current_prec_bits,
                      "for", sum(is.infinite(subset_to_fix$mlogp)), "remaining Inf values."))
      }
      
      rows_still_inf_in_subset_logical = is.infinite(subset_to_fix$mlogp)
      
      beta_val = subset_to_fix$beta[rows_still_inf_in_subset_logical]
      se_val = subset_to_fix$se[rows_still_inf_in_subset_logical]
      
      valid_se_indices = which(se_val > 0 & !is.na(se_val) & !is.na(beta_val)) # Also ensure beta/se are not NA
      
      if (length(valid_se_indices) > 0) {
        # Process only those with valid SE and non-NA beta/se
        beta_to_process = beta_val[valid_se_indices]
        se_to_process = se_val[valid_se_indices]
        
        beta_mpfr = Rmpfr::mpfr(beta_to_process, precBits = current_prec_bits)
        se_mpfr   = Rmpfr::mpfr(se_to_process, precBits = current_prec_bits)
        
        z_mpfr = abs(beta_mpfr / se_mpfr)
        
        # Calculate two-tailed p-value using Rmpfr::pnorm
        p_val_mpfr = 2 * Rmpfr::pnorm(z_mpfr, lower.tail = FALSE) 
        
        # Calculate -log10(p), Convert to numeric, then round
        numeric_mlogp = Rmpfr::asNumeric(-log10(p_val_mpfr))
        numeric_mlogp = round(numeric_mlogp, 5) 
        
        # Update only the rows that were Inf and had valid SEs and non-NA beta/se
        temp_mlogp_updates <- subset_to_fix$mlogp[rows_still_inf_in_subset_logical]
        temp_mlogp_updates[valid_se_indices] <- numeric_mlogp
        subset_to_fix$mlogp[rows_still_inf_in_subset_logical] <- temp_mlogp_updates
      }
      
      if (!any(is.infinite(subset_to_fix$mlogp))) {
        if (verbose) message("All resolvable Inf mlogp values processed with current precision.")
        # Check if all originally Inf values are resolved before breaking entirely
        # Some might remain Inf due to SE<=0 or NA beta/se
        # The loop should continue if precision can still be increased for any remaining Inf values
        # that *could* be resolved (i.e. had SE > 0 initially but maybe Rmpfr still gave Inf)
        # The current logic of any(is.infinite(subset_to_fix$mlogp)) is okay, as it will try
        # higher precision for any remaining Infs. If they are Inf due to SE<=0, they won't change.
      }
      
      current_prec_bits = current_prec_bits + prec_step
    }
    
    df$mlogp[inf_mlogp_indices] = subset_to_fix$mlogp
    
    # Recalculate remaining Infs in the subset (subset_to_fix might have changed)
    remaining_inf_count_in_subset = sum(is.infinite(subset_to_fix$mlogp))
    if (remaining_inf_count_in_subset > 0) {
      last_attempted_prec = current_prec_bits - prec_step # Precision used in the last completed iteration
      if (iter_count == max_iter && last_attempted_prec <= max_prec) {
        # Hit max_iter
        warning(paste(remaining_inf_count_in_subset,
                      "mlogp values remain Inf after attempting high-precision calculation for",
                      max_iter, "iterations (last precision tried: precBits =", last_attempted_prec, "). Consider increasing max_iter or prec_step/max_prec."))
      } else if (last_attempted_prec >= max_prec && any(is.infinite(subset_to_fix$mlogp))) {
        # Hit max_prec
        warning(paste(remaining_inf_count_in_subset,
                      "mlogp values remain Inf after attempting high-precision calculation up to max_precBits =",
                      max_prec, ". Consider increasing max_prec."))
      } else if (all(subset_to_fix$se[is.infinite(subset_to_fix$mlogp)] <= 0 | 
                     is.na(subset_to_fix$se[is.infinite(subset_to_fix$mlogp)]) | 
                     is.na(subset_to_fix$beta[is.infinite(subset_to_fix$mlogp)])) &&
                 length(inf_mlogp_indices) > 0 && verbose) {
        # All remaining Infs are due to invalid SE/beta
        message(paste(remaining_inf_count_in_subset,
                      "mlogp values remain Inf due to non-positive SE or NA beta/se values."))
      }
      # General warning if specific conditions above not met, but still Infs
      else if (verbose) {
        warning(paste(remaining_inf_count_in_subset,
                      "mlogp values may remain Inf. Last precision tried:", last_attempted_prec,
                      "bits after", iter_count, "iterations."))
      }
    } else if (length(inf_mlogp_indices) > 0 && verbose) {
      message("All initially Inf mlogp values were successfully resolved or handled.")
    }
  }
  
  return(df)
}