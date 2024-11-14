args <- commandArgs(trailingOnly = TRUE)
print(args)
# Argument Descriptions

#1. Argument one is the directory. Example: `SampleData1`
#2. Argument two is the file name. Example: `SampleData1\\Fold_0`
#3. Argument three is the output file name. Example: `train_data`
#4. Argument four is the specific function to be called. Example: `train_data.QC.clumped.pruned`

#5. Argument five is LDpred-2 option. Example: `LDpred-2_full` or `LDpred-2_hapmap`
#6. Argument six is the size parameter. Example: `200`
#7. Argument seven is the alpha parameter. Example: `1`
#8. Argument eight is the thr_r2 parameter. Example: `0.1`
#9. Argument nine is the number of PCA. Example: `6`



if (args[5]=="1"){
  cran_mirror_url <- "https://cran.r-project.org"
  install.packages("remotes", repos = cran_mirror_url)
  library(remotes)
  #remotes::install_github("https://github.com/privefl/bigsnpr.git")
  library(bigsnpr)
  options(bigstatsr.check.parallel.blas = FALSE)
  options(default.nproc.blas = NULL)
  library(data.table)
  library(magrittr)
  info <- readRDS(runonce::download_file(
    "https://ndownloader.figshare.com/files/25503788",
    fname = "map_hm3_ldpred2.rds"))
  
  library(bigsnpr)
  options(bigstatsr.check.parallel.blas = FALSE)
  options(default.nproc.blas = NULL)
  library(data.table)
  library(magrittr)
  help(snp_cor)
}

if (args[5]=="2"){
  
  library(bigsnpr)
  options(bigstatsr.check.parallel.blas = FALSE)
  options(default.nproc.blas = NULL)
  library(data.table)
  library(magrittr)
  result <-paste(".",args[2],paste(args[3],toString(".PHENO"), sep = ""),sep="//")
  phenotype <- fread(result)
  result <-paste(".",args[2],paste(args[3],toString(".cov"), sep = ""),sep="//")
  covariate <- fread(result)
  result <-paste(".",args[2],paste(args[3],toString(".eigenvec"), sep = ""),sep="//")
  pcs <- fread(result)
  # rename columns
  colnames(pcs) <- c("FID","IID", paste0("PC",1:as.numeric(args[9])))
  # generate required table
  pheno <- merge(phenotype, covariate) %>%
    merge(., pcs)
  info <- readRDS(runonce::download_file(
    "https://ndownloader.figshare.com/files/25503788",
    fname = "map_hm3_ldpred2.rds"))
  # Read in the summary statistic file
  result <-paste(".",args[1],paste(args[1],toString(".txt"), sep = ""),sep="//")
  
  sumstats <- bigreadr::fread2(result) 
  # LDpred 2 require the header to follow the exact naming
  names(sumstats) <-
    c("chr",
      "pos",
      "rsid",
      "a1",
      "a0",
      "n_eff",
      "beta_se",
      "p",
      "BETA",
      "INFO",
      "MAF")
  # Transform the OR into log(OR)
  sumstats$beta <-  sumstats$BETA
  # Filter out hapmap SNPs
  sumstats <- sumstats[sumstats$rsid%in% info$rsid,]
  
  # Get maximum amount of cores
  NCORES <- nb_cores()
  # Open a temporary file
  
  result <-paste(".",args[2],"tmp-data",sep="//")

  if (dir.exists(result)) {
    # Delete the directory and its contents
    
    system(paste("rm -r", shQuote(result)))
    print(paste("Directory", result, "deleted."))
  }
  tmp <- tempfile(tmpdir = result)
  on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
    
  # Initialize variables for storing the LD score and LD matrix
  corr <- NULL
  ld <- NULL
  # We want to know the ordering of samples in the bed file 
  fam.order <- NULL
  # preprocess the bed file (only need to do once for each data set)
  result <-paste(".",args[2],paste(args[4],toString(".rds"), sep = ""),sep="//")
  if (file.exists(result)) {
    file.remove(result)
    print(paste("File", result, "deleted."))
  }
  result <-paste(".",args[2],paste(args[4],toString(".bk"), sep = ""),sep="//")
  if (file.exists(result)) {
    file.remove(result)
    print(paste("File", result, "deleted."))
  }
  
  
  
  result <-paste(".",args[2],paste(args[4],toString(".bed"), sep = ""),sep="//")
  
  snp_readBed(result)
  # now attach the genotype object
  result <-paste(".",args[2],paste(args[4],toString(".rds"), sep = ""),sep="//")
  
  obj.bigSNP <- snp_attach(result)
  
  # extract the SNP information from the genotype
  map <- obj.bigSNP$map[-3]
  
  names(map) <- c("chr", "rsid", "pos", "a1", "a0")
  
  # perform SNP matching
  info_snp <- snp_match(sumstats, map)
  help(snp_match)
  info_snp
  # Assign the genotype to a variable for easier downstream analysis
  genotype <- obj.bigSNP$genotypes
  # Rename the data structures
  CHR <- map$chr
  POS <- map$pos
  # get the CM information from 1000 Genome
  # will download the 1000G file to the current directory (".")
  #help(snp_asGeneticPos)
  POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
  
 
 check <-TRUE
  for (chr in 1:22) {
    # Extract SNPs that are included in the chromosome
    
    ind.chr <- which(info_snp$chr == chr)
    print(length(ind.chr))
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    ind.chr2
    
    if (length(ind.chr2) == 0) {
       
      next  
    }
    else{
       
      corr0 <- snp_cor(
      genotype,
      ind.col = ind.chr2,
      ncores = NCORES,
      infos.pos = POS2[ind.chr2],
      #size = 200,
      #thr_r2=0.1,
      #alpha = 1
      
      size = as.numeric(args[6]),
      alpha = as.numeric(args[7]),
      
      thr_r2=as.numeric(args[8]),
    )
    if (check==TRUE) {
      check <-FALSE
   
      ld <- Matrix::colSums(corr0^2)
      corr <- as_SFBM(corr0, tmp)
    } else {
      ld <- c(ld, Matrix::colSums(corr0^2))
      corr$add_columns(corr0, nrow(corr))
    }

    }
   
    }
  
  
  # We assume the fam order is the same across different chromosomes
  fam.order <- as.data.table(obj.bigSNP$fam)
  # Rename fam order
  setnames(fam.order,
           c("family.ID", "sample.ID"),
           c("FID", "IID"))
  
  df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
  
  length(df_beta$beta) 
  length(ld)
  help(snp_ldsc)
  ldsc <- snp_ldsc(ld, 
                   length(ld), 
                   chi2 = (df_beta$beta / df_beta$beta_se)^2,
                   sample_size = df_beta$n_eff, 
                   blocks = NULL)
  h2_est <- ldsc[["h2"]]
  h2_est
  
  result <-paste(".",args[2],"ldpred_h2_hapmap.txt",sep="//")
  if (file.exists(result)) {
    file.remove(result)
    print(paste("File", result, "deleted."))
  }
  write.table(h2_est, file = result, col.names = FALSE)
  
  result <-paste(".",args[2],"ldpred_h2_variants.txt",sep="//")
    
  if (file.exists(result)) {
    file.remove(result)
    print(paste("File", result, "deleted."))
  }
  write.table(length(ld), file = result, col.names = FALSE)
  
  
  delta <- c(0.001, 0.01, 0.1, 1)
  nlambda <- 10
   
  
  beta_lassosum2 <- snp_lassosum2(corr, df_beta,delta,nlambda )
  (params2 <- attr(beta_lassosum2, "grid_param"))
  params2
  gridparamters <- params2[,c("lambda", "delta", "sparsity")]
  result <-paste(".",args[2],paste(args[3],toString(".ldpred_lassosum_parameters"), sep = ""),sep="//")
  write.table(gridparamters, file = result, row.names = FALSE,sep=",", quote = FALSE)
  
  
  result <-paste(".",args[2],paste(args[3],toString(".ldpred_lassosum_betas"), sep = ""),sep="//")
  write.table(beta_lassosum2, file = result, row.names = FALSE,sep=",", quote = FALSE)
  newgwas <- info_snp[,c("rsid.ss", "a0", "beta")]
  
  result <-paste(".",args[2],paste(args[3],toString(".ldpred_lassosum_gwas"), sep = ""),sep="//")
  write.table(newgwas, file = result, row.names = FALSE,sep=",", quote = FALSE)
  
  
  
  
}
if (args[5]=="3"){
  
  library(bigsnpr)
  options(bigstatsr.check.parallel.blas = FALSE)
  options(default.nproc.blas = NULL)
  library(data.table)
  library(magrittr)
  result <-paste(".",args[2],paste(args[3],toString(".PHENO"), sep = ""),sep="//")
  phenotype <- fread(result)
  result <-paste(".",args[2],paste(args[3],toString(".cov"), sep = ""),sep="//")
  covariate <- fread(result)
  result <-paste(".",args[2],paste(args[3],toString(".eigenvec"), sep = ""),sep="//")
  pcs <- fread(result)
  # rename columns
  colnames(pcs) <- c("FID","IID", paste0("PC",1:as.numeric(args[9])))
  # generate required table
  pheno <- merge(phenotype, covariate) %>%
    merge(., pcs)
  info <- readRDS(runonce::download_file(
    "https://ndownloader.figshare.com/files/25503788",
    fname = "map_hm3_ldpred2.rds"))
  # Read in the summary statistic file
  result <-paste(".",args[1],paste(args[1],toString(".txt"), sep = ""),sep="//")
  
  sumstats <- bigreadr::fread2(result) 
  # LDpred 2 require the header to follow the exact naming
  names(sumstats) <-
    c("chr",
      "pos",
      "rsid",
      "a1",
      "a0",
      "n_eff",
      "beta_se",
      "p",
      "BETA",
      "INFO",
      "MAF")
  # Transform the OR into log(OR)
  sumstats$beta <-  sumstats$BETA 
  # Filter out hapmap SNPs
  # Turn off this line to ensure that all the SNPs from
  # the sumstats are included.
  #sumstats <- sumstats[sumstats$rsid%in% info$rsid,]
  
  # Get maximum amount of cores
  NCORES <- nb_cores()
  # Open a temporary file
  
  if (dir.exists("tmp-data")) {
    # Delete the directory and its contents
    
    system(paste("rm -r", shQuote("tmp-data")))
    print(paste("Directory", "tmp-data", "deleted."))
  }
  tmp <- tempfile(tmpdir = "tmp-data")
  on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
  # Initialize variables for storing the LD score and LD matrix
  corr <- NULL
  ld <- NULL
  # We want to know the ordering of samples in the bed file 
  fam.order <- NULL
  # preprocess the bed file (only need to do once for each data set)
  result <-paste(".",args[2],paste(args[4],toString(".rds"), sep = ""),sep="//")
  if (file.exists(result)) {
    file.remove(result)
    print(paste("File", result, "deleted."))
  }
  result <-paste(".",args[2],paste(args[4],toString(".bk"), sep = ""),sep="//")
  if (file.exists(result)) {
    file.remove(result)
    print(paste("File", result, "deleted."))
  }
  
  
  result <-paste(".",args[2],paste(args[4],toString(".bed"), sep = ""),sep="//")
  
  snp_readBed(result)
  # now attach the genotype object
  result <-paste(".",args[2],paste(args[4],toString(".rds"), sep = ""),sep="//")
  
  obj.bigSNP <- snp_attach(result)
  
  # extract the SNP information from the genotype
  map <- obj.bigSNP$map[-3]
  
  names(map) <- c("chr", "rsid", "pos", "a1", "a0")
  
  # perform SNP matching
  info_snp <- snp_match(sumstats, map)
  help(snp_match)
  info_snp
  # Assign the genotype to a variable for easier downstream analysis
  genotype <- obj.bigSNP$genotypes
  # Rename the data structures
  CHR <- map$chr
  POS <- map$pos
  # get the CM information from 1000 Genome
  # will download the 1000G file to the current directory (".")
  help(snp_asGeneticPos)
  POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
  
   

    check <-TRUE
  for (chr in 1:22) {
    # Extract SNPs that are included in the chromosome
    
    ind.chr <- which(info_snp$chr == chr)
    print(length(ind.chr))
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    
    print(length(ind.chr2))
    
    if (length(ind.chr2) == 0) {
       
      next  
    }
    else{
       
      corr0 <- snp_cor(
      genotype,
      ind.col = ind.chr,
      ncores = NCORES,
      infos.pos = POS2[ind.chr2],
      #size = 200,
      #thr_r2=0.1,
      #alpha = 1
      
      size = as.numeric(args[6]),
      alpha = as.numeric(args[7]),
      thr_r2=as.numeric(args[8]),
    )
    if (check==TRUE) {
      check <-FALSE
    
      ld <- Matrix::colSums(corr0^2)
      help(as_SFBM)
      corr <- as_SFBM(corr0, tmp)
    } else {
      ld <- c(ld, Matrix::colSums(corr0^2))
      corr$add_columns(corr0, nrow(corr))
    }

    }
   
    }
  
  # We assume the fam order is the same across different chromosomes
  fam.order <- as.data.table(obj.bigSNP$fam)
  # Rename fam order
  setnames(fam.order,
           c("family.ID", "sample.ID"),
           c("FID", "IID"))
  
  df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
  
  length(df_beta$beta) 
  length(ld)
  help(snp_ldsc)
  ldsc <- snp_ldsc(ld, 
                   length(ld), 
                   chi2 = (df_beta$beta / df_beta$beta_se)^2,
                   sample_size = df_beta$n_eff, 
                   blocks = NULL)
  h2_est <- ldsc[["h2"]]
 
  delta <- c(0.001, 0.01, 0.1, 1)
  nlambda <- 10
 
  
  beta_lassosum2 <- snp_lassosum2(corr, df_beta,delta,nlambda )
  (params2 <- attr(beta_lassosum2, "grid_param"))
  params2
  gridparamters <- params2[,c("lambda", "delta", "sparsity")]
  
  result <-paste(".",args[2],paste(args[3],toString(".ldpred_lassosum_parameters"), sep = ""),sep="//")
  write.table(gridparamters, file = result, row.names = FALSE,sep=",", quote = FALSE)
   
  
  result <-paste(".",args[2],paste(args[3],toString(".ldpred_lassosum_betas"), sep = ""),sep="//")
  write.table(beta_lassosum2, file = result, row.names = FALSE,sep=",", quote = FALSE)
  newgwas <- info_snp[,c("rsid.ss", "a0", "beta")]
  
  result <-paste(".",args[2],paste(args[3],toString(".ldpred_lassosum_gwas"), sep = ""),sep="//")
  write.table(newgwas, file = result, row.names = FALSE,sep=",", quote = FALSE)
  
   
    result <-paste(".",args[2],"ldpred_h2_full.txt",sep="//")
  if (file.exists(result)) {
    file.remove(result)
    print(paste("File", result, "deleted."))
  }
  write.table(h2_est, file = result, col.names = FALSE)
  
  result <-paste(".",args[2],"ldpred_h2_variants.txt",sep="//")
    
  if (file.exists(result)) {
    file.remove(result)
    print(paste("File", result, "deleted."))
  }
  write.table(length(ld), file = result, col.names = FALSE)
  
   
  
  
}

