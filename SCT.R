args <- commandArgs(trailingOnly = TRUE)
print(args)
# Argument Descriptions

#1. Argument one is the directory. Example: `SampleData1`
#2. Argument two is the file name. Example: `SampleData1\\Fold_0`
#3. Argument three is the output file name. Example: `train_data`
#4. Argument four is the specific function to be called. Example: `train_data.QC.clumped.pruned`

#5. Argument five is invokes specific section of the code.

 



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
  str(obj.bigSNP, max.level = 2, strict.width = "cut") 
   G   <- obj.bigSNP$genotypes
    CHR <- obj.bigSNP$map$chromosome
    POS <- obj.bigSNP$map$physical.pos
    y   <- obj.bigSNP$fam$affection - 1
    NCORES <- nb_cores()
    # Check some counts for the 10 first variants
    big_counts(G, ind.col = 1:10) 
   
  result <-paste(".",args[1],paste(args[1],toString("SCT.txt"), sep = ""),sep="//")
  
  sumstats <- bigreadr::fread2(result) 
  set.seed(1)
 
 ind.train <- sample(nrow(G) )
 
    

    names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "p")
    #sumstats$pos <- as.character(sumstats$pos)
    #map$pos <- as.character(map$pos)
    
    
    
map <- obj.bigSNP$map[,-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
    
info_snp <- snp_match(sumstats, map)
   
    info_snp <- snp_match(sumstats, map, strand_flip = FALSE)
    # beta and lpval need to have the same length as ncol(G), CHR and POS
    # -> one solution is to use missing values and use the 'exclude' parameter
    beta <- rep(NA, ncol(G))
    beta[info_snp$`_NUM_ID_`] <- info_snp$beta
    
    lpval <- rep(NA, ncol(G))
    lpval[info_snp$`_NUM_ID_`] <- -log10(info_snp$p)
    lpval[is.na(lpval)] <- 0
    beta[is.na(beta)] <- 0 
    
 
    
    # The clumping step might take some time to complete
    all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train,
                                  lpS = lpval, exclude = which(is.na(lpval)),
                                  ncores = NCORES)
        attr(all_keep, "grid")
         
    
      result <-paste(".",args[2],"public-data-scores.bk",sep="//")
      print(result)
      if (file.exists(result)) {
        file.remove(result)
        print(paste("File", result, "deleted."))
      }
    result <-paste(".",args[2],"public-data-scores",sep="//")
      
 
     
    multi_PRS <- snp_grid_PRS(G, all_keep, beta, lpval, ind.row = ind.train,
                          backingfile = result, 
                          n_thr_lpS = 1, ncores = NCORES)
      
dim(multi_PRS)  ## 4200 C+T scores for 400 individuals''
    
    final_mod <- snp_grid_stacking(multi_PRS, y[ind.train], ncores = NCORES, K = 4)
 
    
    new_beta <- final_mod$beta.G
    ind <- which(new_beta != 0)
    info_snp$newbeta <-new_beta
    print(head(new_beta))
  result <-paste(".",args[2],paste(args[3],toString(".sct"), sep = ""),sep="//")
  write.table(info_snp, file = result, row.names = FALSE,sep="\t", quote = FALSE)
   
  
  
}
 

