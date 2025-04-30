# clump.R
# This script contains clump_local function to clump locally using plink and reference panel

# local clumping
clump_local = function(df, clump_r2){
  # input:
  #   - df: a data.frame, need to have at least SNP and mlogp column
  #         mlogp column shoud NOT contain Inf
  #   - clump_r2: clump r2 threshold
  # output:
  #   - output_df, a data.frame containing exact same columns as df with rows clumped
  
  # init constant
  clump_kb = 10000
  clump_p1 = 1
  clump_p2 = 1
  
  # path to plink and reference panel
  bfile = paste0("util/plink/EUR") # EUR reference panel
  plink_bin = paste0("util/plink/plink")
  
  pval_column = "mlogp"
  
  # check if mlogp col is in df
  if (!pval_column %in% colnames(df)){
    stop("Error: dataframe does not contain `mlogp` column for -log10(pval)")
  }
  
  # check if mlogp col contains Inf
  if (any(is.infinite(df$mlogp))) {
    stop("Error: dataframe contains infinate mlogp column")
  }
  
  # filter out rows with pval > 0.05 in case haven't
  sig_df = df[df$mlogp > -log10(0.05), ]
  
  # filter out rows without rsID and create input df
  input_df = sig_df[grepl("^rs", sig_df$SNP), c('SNP', 'mlogp')]
  
  if (nrow(input_df) == 0){
    stop("clumping stopped, no significant SNPs")
  }
  
  if (nrow(input_df)!=1){
    
    # select default shell by operating system
    shell=ifelse(Sys.info()["sysname"] == "Windows", "cmd", "sh")
    
    # generate temporary random filename 
    fn = paste0(
      getwd(),
      "/", 
      do.call(
        paste0, replicate(10, sample(LETTERS, 1, TRUE), FALSE)
      )
    )
    
    # save as a random file
    write.table(
      data.frame(
        SNP = input_df$SNP, 
        P = 1/input_df$mlogp # as a proxy for pval ranking
      ), 
      file = fn,
      row.names = F, 
      col.names = T, 
      quote = F
    )
    
    # command to execute
    fun2=paste0(
      shQuote(plink_bin, type = shell),
      " --bfile ",shQuote(bfile, type = shell),
      " --clump ",shQuote(fn, type = shell),
      " --clump-p1 ",clump_p1, 
      " --clump-r2 ",clump_r2,
      " --clump-kb ",clump_kb,
      " --out ",shQuote(fn,type = shell)
    )
    
    # execute command in shell
    print("Start clumping using PLINK")
    system(fun2,intern = T)
    
    # check if clump success
    if (file.exists(paste(fn, ".clumped", sep = ""))){
      
      # read output table
      res = read.table(paste(fn, ".clumped", sep = ""), header = T)
      
      # delete all temporary files
      unlink(paste(fn, "*", sep = ""))
      
      # check number of SNP removed
      removed_SNPs = subset(input_df, !input_df[["SNP"]] %in% res[["SNP"]])
      
      if (nrow(removed_SNPs) > 0){
        message("Removing ", 
                length(removed_SNPs[["SNP"]]), 
                " of ", 
                nrow(input_df), 
                " variants due to LD with other variants or absence from LD reference panel")
      }
      
      output_df = subset(sig_df, sig_df[["SNP"]] %in% res[["SNP"]])
      
    } else {
      
      unlink(paste(fn, "*", sep = ""))
      stop("PLINK clumping failed, cannot find `.clumped` output file")
      
    }
    
  } else if (nrow(input_df)!=1){
    
    output_df = subset(sig_df, sig_df[["SNP"]] %in% input_df[["SNP"]])
    
  }
  return(output_df)
}