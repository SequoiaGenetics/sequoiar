#' Map a SNP to its gene region or closest gene
#'
#' This function checks whether a single SNP (row) falls within any gene regions from a provided
#' gene annotation data frame. It returns whether the SNP is located within a gene, in multiple
#' overlapping genes, or outside all genes (in which case the closest gene is reported).
#'
#' @param snp_row A one-row data frame representing a SNP. Must contain `SNP`, `chr`, and `pos` columns.
#' @param gene_df A data frame of gene annotations with columns: `gene`, `chr`, `gene_start_pos`, `gene_end_pos`.
#'
#' @return A data frame with one or more rows containing the SNP ID, chromosome, position,
#' region status (`"in_one_region"`, `"in_multiple_region"`, `"not_in_any_region"`, or `"no_genes_on_chr"`),
#' and the closest or overlapping gene(s) in a column called `closest_region`.
#'
#' @importFrom dplyr filter arrange slice mutate select rename

#' @examples
#' # Example SNP and gene data
#' df_genes = RACER::hg19 %>%
#'   filter(gene_type == 'protein_coding') %>%
#'   select(gene = gene_name, gene_id, chr = chromosome, gene_start_pos = gene_start, gene_end_pos = gene_end)
#' 
#' # get closest gene
#' mapped_snps = cluster_result %>%
#'   rowwise() %>%
#'   do(map_snp_to_gene(., df_genes)) %>%
#'   ungroup()
#'
#' @export
sg.map_snps_to_genes <- function(snp_row, gene_df) {
  # Subset to same chromosome
  chr_genes <- dplyr::filter(gene_df, chr == snp_row$chr)
  
  if (nrow(chr_genes) == 0) {
    return(data.frame(
      SNP = snp_row$SNP,
      chr = snp_row$chr,
      pos = snp_row$pos,
      region_status = "no_genes_on_chr",
      closest_region = NA
    ))
  }
  
  # Check if SNP is within any gene regions
  in_gene <- dplyr::filter(chr_genes,
                           gene_start_pos <= snp_row$pos & gene_end_pos >= snp_row$pos)
  
  if (nrow(in_gene) == 0) {
    # SNP is outside all genes → find closest gene
    chr_genes$distance <- pmin(abs(chr_genes$gene_start_pos - snp_row$pos),
                               abs(chr_genes$gene_end_pos - snp_row$pos))
    chr_genes <- chr_genes[order(chr_genes$distance), ]
    closest_gene <- chr_genes[1, ]
    
    return(data.frame(
      SNP = snp_row$SNP,
      chr = snp_row$chr,
      pos = snp_row$pos,
      region_status = "not_in_any_region",
      closest_region = closest_gene$gene
    ))
    
  } else if (nrow(in_gene) == 1) {
    return(data.frame(
      SNP = snp_row$SNP,
      chr = snp_row$chr,
      pos = snp_row$pos,
      region_status = "in_one_region",
      closest_region = in_gene$gene
    ))
    
  } else {
    # Multiple overlapping gene regions
    in_gene$SNP <- snp_row$SNP
    in_gene$chr <- snp_row$chr
    in_gene$pos <- snp_row$pos
    in_gene$region_status <- "in_multiple_region"
    
    result <- in_gene[, c("SNP", "chr", "pos", "region_status", "gene")]
    colnames(result)[colnames(result) == "gene"] <- "closest_region"
    
    return(result)
  }
}

