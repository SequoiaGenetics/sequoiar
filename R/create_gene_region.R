#' Create a gene region list from RACER annotations
#'
#' Given a gene name and genome build, this function retrieves the chromosome,
#' start, and end coordinates of the specified protein-coding gene using
#' annotations from the `RACER` package.
#'
#' @param gene_name_chr A character string. The gene name to look up (must match RACER annotation).
#' @param build An integer, either `38` (default) or `37`, specifying the genome build to use.
#'
#' @return A named list containing `gene_name`, `chr`, `lower` (start), and `upper` (end) if the gene is found uniquely;
#' otherwise returns `NULL`.
#'
#' @importFrom dplyr filter pull
#' @export

sg.create_gene_region = function(gene_name_chr, build = 38){
    # Load gene annotation based on build
    if (build == 38) {
      gene_df = dplyr::filter(RACER::hg38, gene_type == "protein_coding")
    } else if (build == 37) {
      gene_df = dplyr::filter(RACER::hg19, gene_type == "protein_coding")
    } else {
      stop("Unsupported genome build. Use 37 or 38.")
    }
    
    # Filter by gene name
    gene_df = dplyr::filter(gene_df, gene_name == gene_name_chr)
    
    # Must find exactly one match
    if (nrow(gene_df) != 1) {
      message("Gene match not unique or not found.")
      return(NULL)
    }
    
    return(list(
      gene_name = gene_name_chr,
      chr = dplyr::pull(gene_df, chromosome),
      lower = dplyr::pull(gene_df, gene_start),
      upper = dplyr::pull(gene_df, gene_end)
    ))
  }