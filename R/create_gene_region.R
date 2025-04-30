library(RACER)

create_gene_region = function(gene_name_chr, build=38){
  if (build == 38){
    gene_df = hg38 %>% dplyr::filter(gene_type=='protein_coding')
  } else if (build == 37) {
    gene_df = hg19 %>% dplyr::filter(gene_type=='protein_coding')
  }
  gene_df = gene_df %>% dplyr::filter(gene_name == gene_name_chr)
  if (dim(gene_df)[1] != 1){
    print('length is not 1')
    return(NULL)
  }
  return(
    list(gene_name=gene_name_chr,
         chr=gene_df %>% pull(chromosome),
         lower=gene_df %>% pull(gene_start),
         upper=gene_df %>% pull(gene_end))
  )
}