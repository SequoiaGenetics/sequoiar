#' @importFrom magrittr %>%
NULL

# Create a function to check if an SNP is within a gene region
map_snp_to_gene = function(snp_row, gene_df) {
  # Subset genes on the same chromosome
  chr_genes = gene_df %>% dplyr::filter(chr == snp_row$chr)
  
  if (nrow(chr_genes) == 0) {
    return(data.frame(
      SNP = snp_row$SNP,
      chr = snp_row$chr,
      pos = snp_row$pos,
      region_status = "no_genes_on_chr",
      closest_region = NA
    ))
  }
  
  # Check if the SNP position falls within any gene region
  in_gene = chr_genes %>% 
    dplyr::filter(gene_start_pos <= snp_row$pos & gene_end_pos >= snp_row$pos)
  
  # Determine which case: not_in_any_region, in_one_region, in_multiple_region
  if (nrow(in_gene) == 0) {
    # SNP is not in any gene region, find the closest gene
    chr_genes$distance = pmin(abs(chr_genes$gene_start_pos - snp_row$pos), abs(chr_genes$gene_end_pos - snp_row$pos))
    closest_gene = chr_genes %>% dplyr::arrange(distance) %>% dplyr::slice(1)
    return(data.frame(
      SNP = snp_row$SNP,
      chr = snp_row$chr,
      pos = snp_row$pos,
      region_status = "not_in_any_region",
      closest_region = closest_gene$gene
    ))
  } else if (nrow(in_gene) == 1) {
    # SNP is in exactly one gene region
    return(data.frame(
      SNP = snp_row$SNP,
      chr = snp_row$chr,
      pos = snp_row$pos,
      region_status = "in_one_region",
      closest_region = in_gene$gene
    ))
  } else {
    # SNP is in multiple gene regions
    result = in_gene %>%
      dplyr::mutate(SNP = snp_row$SNP, chr = snp_row$chr, pos = snp_row$pos, region_status = "in_multiple_region") %>%
      dplyr::select(SNP, chr, pos, region_status, gene) %>%
      dplyr::rename(closest_region = gene)
    
    return(result)
  }
}

# example
# df_genes = RACER::hg19 %>%
#   filter(gene_type == 'protein_coding') %>%
#   select(gene = gene_name, gene_id, chr = chromosome, gene_start_pos = gene_start, gene_end_pos = gene_end)
# 
# # get closest gene
# mapped_snps = cluster_result %>%
#   rowwise() %>%
#   do(map_snp_to_gene(., df_genes)) %>%
#   ungroup()