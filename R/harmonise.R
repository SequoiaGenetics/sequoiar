# harmonise data

# Harmonise function for an arbitrary number of dataframes
# multi_harmonise = function(exp, ...) {
# Changing multi_harmonising to default harmonising
harmonise = function(exp, ...){
  # Convert additional arguments to a list of dataframes
  dataframes = list(...)
  
  # Get the name of the exp dataframe
  exp_name = deparse(substitute(exp))
  # get the name of other dataframe
  out_names = sapply(substitute(list(...))[-1], deparse)
  
  # Make sure there is no INDEL SNPs
  exp = exp %>%
    filter(nchar(effect_allele) == 1, nchar(other_allele) == 1) %>%
    mutate(effect_allele = toupper(effect_allele), other_allele = toupper(other_allele)) %>%
    filter(effect_allele != "I", 
           effect_allele != "D",
           other_allele != "I",
           other_allele != "D")
  
  # Add allele_concat column to exp dataframe
  exp = exp %>%
    mutate(allele_concat = paste(pmin(effect_allele, other_allele), 
                                 pmax(effect_allele, other_allele), sep = "")) %>%
    rename_with(~ paste0(., '_', exp_name), -c(SNP, allele_concat, effect_allele, other_allele))
  
  # Harmonise each dataframe to the exp dataframe
  harmonised_df = exp
  for (i in seq_along(dataframes)) {
    out = dataframes[[i]]
    out_name = out_names[i]
    suffix_out = paste0('_', out_name)
    
    # Make sure there is no INDEL SNPs
    out = out %>%
      filter(nchar(effect_allele) == 1, nchar(other_allele) == 1) %>%
      mutate(effect_allele = toupper(effect_allele), other_allele = toupper(other_allele)) %>%
      filter(effect_allele != "I", 
             effect_allele != "D",
             other_allele != "I",
             other_allele != "D")
    
    # Add allele_concat column to out dataframe
    out = out %>%
      mutate(allele_concat = paste(pmin(effect_allele, other_allele), 
                                   pmax(effect_allele, other_allele), sep = "")) %>%
      rename_with(~ paste0(., suffix_out), -c(SNP, allele_concat))
    
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
    harmonised_df = harmonised_df %>%
      left_join(out, by = by_cols, suffix = c('', suffix_out))
    
    harmonised_df = harmonised_df %>%
      # create a copy of original effect_allele/other_allele of outcome
      mutate(effect_allele_out_orig = get(paste0('effect_allele_', out_name)),
             other_allele_out_orig = get(paste0('other_allele_', out_name))) %>%
      # Adjust effect_allele/other_allele/beta of outcome based on exposure
      mutate(!!paste0('beta', suffix_out) := ifelse(effect_allele == effect_allele_out_orig,
                                                    get(paste0('beta', suffix_out)),
                                                    -1 * get(paste0('beta', suffix_out))),
             !!paste0('effect_allele', suffix_out) := ifelse(effect_allele == effect_allele_out_orig,
                                                             effect_allele_out_orig,
                                                             other_allele_out_orig),
             !!paste0('other_allele', suffix_out) := ifelse(effect_allele == effect_allele_out_orig,
                                                            other_allele_out_orig,
                                                            effect_allele_out_orig))  %>%
      select(-effect_allele_out_orig, -other_allele_out_orig)
  }
  
  # Rename effect_allele and other_allele with suffix of exposure
  harmonised_df = harmonised_df %>%
    rename_with(~ paste0(., '_', exp_name), c(effect_allele, other_allele)) %>%
    select(-allele_concat)
  
  return(harmonised_df)
}

harmonise_archive = function(exp, out){
  
  # make sure there are no INDEL SNPs
  exp = exp %>%
    filter(nchar(effect_allele) == 1, nchar(other_allele) == 1) %>%
    mutate(effect_allele = toupper(effect_allele), other_allele = toupper(other_allele)) %>%
    filter(effect_allele != "I", 
           effect_allele != "D",
           other_allele != "I",
           other_allele != "D")
  out = out %>%
    filter(nchar(effect_allele) == 1, nchar(other_allele) == 1) %>%
    mutate(effect_allele = toupper(effect_allele), other_allele = toupper(other_allele)) %>%
    filter(effect_allele != "I", 
           effect_allele != "D",
           other_allele != "I",
           other_allele != "D")
  
  
  # Add allele_concat column to exp and out dataframes
  exp = exp %>%
    mutate(allele_concat = paste(pmin(effect_allele, other_allele), 
                                 pmax(effect_allele, other_allele), sep = ""))
  
  out = out %>%
    mutate(allele_concat = paste(pmin(effect_allele, other_allele), 
                                 pmax(effect_allele, other_allele), sep = ""))
  
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
  
  harmonised_df = exp %>%
    left_join(out, by = by_cols, suffix = c('_exposure', '_outcome')) %>% # add suffixes to distinguish columns
    # create a copy of original effect_allele/other_allele of outcome
    mutate(effect_allele_out_orig = effect_allele_outcome,
           other_allele_out_orig = other_allele_outcome) %>%
    # Adjust effect_allele/other_allele/beta of outcome based on exposure
    mutate(beta_outcome = ifelse(effect_allele_exposure == effect_allele_out_orig,
                                 beta_outcome,
                                 -1*beta_outcome),
           effect_allele_outcome = ifelse(effect_allele_exposure == effect_allele_out_orig,
                                          effect_allele_out_orig,
                                          other_allele_out_orig),
           other_allele_outcome = ifelse(effect_allele_exposure == effect_allele_out_orig,
                                         other_allele_out_orig,
                                         effect_allele_out_orig)) %>%
    select(-effect_allele_out_orig, -other_allele_out_orig)
  
  return(harmonised_df)
}