#!/usr/bin/env python
# coding: utf-8

# # LDpred-2-Grid
# 
# In this notebook, we will use LDpred-2-Grid to calculate the Polygenic Risk Score (PRS).
# 
# **Note:** LDpred-2 needs to be installed using R.
# 
# We will use the same flow of the code, and when we have to invoke LDpred-2 functions, we will call the script using Python. When LDpred-2 finishes executions, results will be stored in the specific directories for each fold and will be retrieved in Python for further calculations.
# 
# ## Install LDpred-2
# 
# Install LDpred-2 using the following commands:
# 
# 1. Activate the conda Environment.
# 2. Type R
# 3. Run the following commands:
# 
# ```R
# install.packages("remotes")
# library(remotes)
# remotes::install_github("https://github.com/privefl/bigsnpr.git")
# ```
# 
# The content in this notebook has been taken from the [PRS Tutorial](https://choishingwan.github.io/PRS-Tutorial/ldpred/), [LDpred-2 Documentation](https://privefl.github.io/bigsnpr/articles/LDpred2.html), and [Polygenic Scores (PGS)](https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html).
# 
# Thanks to Florian Privé for creating one of the best documentations explaining each aspect of the calculation. We will be using LDpred-2 for cross-validation design and hyper-parameter optimization, and this tutorial complements the existing one.
# 
# 

# ## LDpred-2  Hyperparameters
# 
# ### Hyperparameters for LDpred-2  performed using PLINK
# 
# #### Pruning Parameters
# 
# Informs Plink that we wish to perform pruning with a window size of 200 variants, sliding across the genome with a step size of 50 variants at a time, and filter out any SNPs with LD \( r^2 \) higher than 0.25.
# 
# 
# ```python
# 1. p_window_size = [200]
# 2. p_slide_size = [50]
# 3. p_LD_threshold = [0.25]
# ```
# 
# #### Clumping Parameters
# 
# The P-value threshold for an SNP to be included. 1 means to include all SNPs for clumping. SNPs having \( r^2 \) higher than 0.1 with the index SNPs will be removed. SNPs within 200k of the index SNP are considered for clumping.
# 
# ```python
# 1. clump_p1 = [1]
# 2. clump_r2 = [0.1]
# 3. clump_kb = [200]
# ```
# #### PCA
# Pca also affects the results evident from the initial analysis; however, including more PCA overfits the model. 
# 
# ### Hyperparameters for LDpred-2 
#  
# 
# 
# #### Heritability
# 
# To calculate h2 using LDpred-2, we followed the [LDpred-2 tutorial](https://choishingwan.github.io/PRS-Tutorial/ldpred/).
# The code for this part is in R, as LDpred-2 is in R.
# 
# 1. Using SNPs from the HapMap as preferred by the authors. The name in the code is `LDpred-2_hapmap`.
# 2. Using all the SNPs. The name in the code is `LDpred-2_full`.
# 
# This approach is computationally expensive, as the correlation between all the SNPs in a specific chromosome is being calculated.
# 
# ## LDpred-2 Models
# 
# Each model has its own set of hyperparameters, so first, we will generate those hyperparameters for those models and call them when invoking LDpred-2.
# 
# Use the following R code to know more about the hyperparameters:
# - `help(snp_ldpred2_inf)`
# - `help(snp_ldpred2_grid)`
# - `help(snp_ldpred2_auto)`
# - `help(snp_lassosum2)`
# 
# ### 1. LDpred2-inf: Infinitesimal Model
# 
# No Hyperparameters
# 
# ### 2. LDpred2(-grid): Grid of Models
# 
# Each model has its own set of parameters:
# - `burn_in = 100`
# - `num_iter = 10`
# - `p = c(1e-4, 1e-5)`
# - `h2 = c(0.1, 0.2)`
# - `sparse = c(FALSE, TRUE)`
# 
# ### 3. LDpred2-auto: Automatic Model
# 
# - `burn_in = 100`
# - `num_iter = 100`
# - `sparse = c(FALSE, TRUE)`
# - `alpha_bounds = c(1, 0.5)`
# - `use_MLE = c(FALSE, TRUE)`
# - `p = c(1e-4, 1e-5)`
# - `shrink_corr = 0.7`
# - `allow_jump_sign = c(FALSE, TRUE)`
# 
# ### 4. lassosum2: Grid of Models
# 
# - `lambda = c(0.001, 0.01, 0.1, 1)`
# - `delta = 30`
# 
# We will execute the code for each mkodel seperatly.
# 

# ### GWAS file processing for Plink.
# When the effect size relates to disease risk and is thus given as an odds ratio (OR) rather than BETA (for continuous traits), the PRS is computed as a product of ORs. To simplify this calculation, take the natural logarithm of the OR so that the PRS can be computed using summation instead.

# In[1]:


import os
import pandas as pd
import numpy as np
import sys
import os
import pandas as pd
import numpy as np

def check_phenotype_is_binary_or_continous(filedirec):
    # Read the processed quality controlled file for a phenotype
    df = pd.read_csv(filedirec+os.sep+filedirec+'_QC.fam',sep="\s+",header=None)
    column_values = df[5].unique()
 
    if len(set(column_values)) == 2:
        return "Binary"
    else:
        return "Continous"

    
filedirec = sys.argv[1]

#filedirec = "SampleData1"
#filedirec = "asthma_19"
#filedirec = "migraine_0"




# Read the GWAS file.
GWAS = filedirec + os.sep + filedirec+".gz"
df = pd.read_csv(GWAS,compression= "gzip",sep="\s+")


# We will save the new GWAS file containing only the betas, # LDpred-2 requires betas for calculation
# which is further used to calculate the new Betas.


if "BETA" in df.columns.to_list():
    # For Continous Phenotype.
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]

else:
    df["BETA"] = np.log(df["OR"])
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]

df.to_csv(filedirec + os.sep +filedirec+".txt",sep="\t",index=False)
print(df.head().to_markdown())
print("Length of DataFrame!",len(df))

 


# ### Define Hyperparameters
# 
# Define hyperparameters to be optimized and set initial values.
# 
# ### Extract Valid SNPs from Clumped File
# 
# For Windows, download `gwak`, and for Linux, the `awk` command is sufficient. For Windows, `GWAK` is required. You can download it from [here](https://sourceforge.net/projects/gnuwin32/). Get it and place it in the same directory.
# 
# 
# ### Execution Path
# 
# At this stage, we have the genotype training data `newtrainfilename = "train_data.QC"` and genotype test data `newtestfilename = "test_data.QC"`.
# 
# We modified the following variables:
# 
# 1. `filedirec = "SampleData1"` or `filedirec = sys.argv[1]`
# 2. `foldnumber = "0"` or `foldnumber = sys.argv[2]` for HPC.
# 
# Only these two variables can be modified to execute the code for specific data and specific folds. Though the code can be executed separately for each fold on HPC and separately for each dataset, it is recommended to execute it for multiple diseases and one fold at a time.
# Here’s the corrected text in Markdown format:
# 
#  
# ### P-values
# 
# PRS calculation relies on P-values. SNPs with low P-values, indicating a high degree of association with a specific trait, are considered for calculation.
# 
# You can modify the code below to consider a specific set of P-values and save the file in the same format.
# 
# We considered the following parameters:
# 
# - **Minimum P-value**: `1e-10`
# - **Maximum P-value**: `1.0`
# - **Minimum exponent**: `10`  (Minimum P-value in exponent)
# - **Number of intervals**: `100`  (Number of intervals to be considered)
# 
# The code generates an array of logarithmically spaced P-values:
# 
# ```python
# import numpy as np
# import os
# 
# minimumpvalue = 10  # Minimum exponent for P-values
# numberofintervals = 100  # Number of intervals to be considered
# 
# allpvalues = np.logspace(-minimumpvalue, 0, numberofintervals, endpoint=True)  # Generating an array of logarithmically spaced P-values
# 
# print("Minimum P-value:", allpvalues[0])
# print("Maximum P-value:", allpvalues[-1])
# 
# count = 1
# with open(os.path.join(folddirec, 'range_list'), 'w') as file:
#     for value in allpvalues:
#         file.write(f'pv_{value} 0 {value}\n')  # Writing range information to the 'range_list' file
#         count += 1
# 
# pvaluefile = os.path.join(folddirec, 'range_list')
# ```
# 
# In this code:
# - `minimumpvalue` defines the minimum exponent for P-values.
# - `numberofintervals` specifies how many intervals to consider.
# - `allpvalues` generates an array of P-values spaced logarithmically.
# - The script writes these P-values to a file named `range_list` in the specified directory.
#  

# ### 2. LDpred2(-grid): Grid of Models
# 
# The parameters for the grid model can be divided into two classes: the first one should be specified here in the python code and the other one should be specified in the corresponding R file.
# 
# ### Should be specified here.
# - `burn_in = 100`
# - `num_iter = 10`
# 
# ### Should be specified in R file - LDpred2-grid.R
# .
# - `p = c(1e-4, 1e-5)`
# - `h2 = c(0.1, 0.2)`
# - `sparse = c(FALSE, TRUE)`
# 
# For each combination of hyper parameters, LDpred-2 grid model will generate a set of betas, and we will be using those betas for each SNP to compute the polygenic risks using `plink –score`.
# 
# 
# 
# The following file contains the R code required to execute the code in this notebook.
# 
# ```r
# args <- commandArgs(trailingOnly = TRUE)
# print(args)
# # Argument Descriptions
# 
# #1. Argument one is the directory. Example: `SampleData1`
# #2. Argument two is the file name. Example: `SampleData1\\Fold_0`
# #3. Argument three is the output file name. Example: `train_data`
# #4. Argument four is the specific function to be called. Example: `train_data.QC.clumped.pruned`
# 
# #5. Argument five is LDpred-2 option. Example: `LDpred-2_full` or `LDpred-2_hapmap`
# #6. Argument six is the size parameter. Example: `200`
# #7. Argument seven is the alpha parameter. Example: `1`
# #8. Argument eight is the thr_r2 parameter. Example: `0.1`
# #9. Argument nine is burn_in. Example: `100`
# #10 Argument nine is num_iter. Example: `10`
# 
# 
# if (args[5]=="1"){
#   cran_mirror_url <- "https://cran.r-project.org"
#   install.packages("remotes", repos = cran_mirror_url)
#   library(remotes)
#   #remotes::install_github("https://github.com/privefl/bigsnpr.git")
#   library(bigsnpr)
#   options(bigstatsr.check.parallel.blas = FALSE)
#   options(default.nproc.blas = NULL)
#   library(data.table)
#   library(magrittr)
#   info <- readRDS(runonce::download_file(
#     "https://ndownloader.figshare.com/files/25503788",
#     fname = "map_hm3_ldpred2.rds"))
#   
#   library(bigsnpr)
#   options(bigstatsr.check.parallel.blas = FALSE)
#   options(default.nproc.blas = NULL)
#   library(data.table)
#   library(magrittr)
#   help(snp_cor)
# }
# 
# if (args[5]=="2"){
#   
#   library(bigsnpr)
#   options(bigstatsr.check.parallel.blas = FALSE)
#   options(default.nproc.blas = NULL)
#   library(data.table)
#   library(magrittr)
#   result <-paste(".",args[2],paste(args[3],toString(".PHENO"), sep = ""),sep="//")
#   phenotype <- fread(result)
#   result <-paste(".",args[2],paste(args[3],toString(".cov"), sep = ""),sep="//")
#   covariate <- fread(result)
#   result <-paste(".",args[2],paste(args[3],toString(".eigenvec"), sep = ""),sep="//")
#   pcs <- fread(result)
#   # rename columns
#   colnames(pcs) <- c("FID","IID", paste0("PC",1:as.numeric(args[9])))
#   # generate required table
#   pheno <- merge(phenotype, covariate) %>%
#     merge(., pcs)
#   info <- readRDS(runonce::download_file(
#     "https://ndownloader.figshare.com/files/25503788",
#     fname = "map_hm3_ldpred2.rds"))
#   # Read in the summary statistic file
#   result <-paste(".",args[1],paste(args[1],toString(".txt"), sep = ""),sep="//")
#   
#   sumstats <- bigreadr::fread2(result) 
#   # LDpred 2 require the header to follow the exact naming
#   names(sumstats) <-
#     c("chr",
#       "pos",
#       "rsid",
#       "a1",
#       "a0",
#       "n_eff",
#       "beta_se",
#       "p",
#       "BETA",
#       "INFO",
#       "MAF")
#   # Transform the OR into log(OR)
#   sumstats$beta <- sumstats$BETA
#   # Filter out hapmap SNPs
#   sumstats <- sumstats[sumstats$rsid%in% info$rsid,]
#   
#   # Get maximum amount of cores
#   NCORES <- nb_cores()
#   # Open a temporary file
#   
#   result <-paste(".",args[2],"tmp-data",sep="//")
# 
#   if (dir.exists(result)) {
#     # Delete the directory and its contents
#     
#     system(paste("rm -r", shQuote(result)))
#     print(paste("Directory", result, "deleted."))
#   }
#   tmp <- tempfile(tmpdir = result)
#   on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
#     
#     
#   # Initialize variables for storing the LD score and LD matrix
#   corr <- NULL
#   ld <- NULL
#   # We want to know the ordering of samples in the bed file 
#   fam.order <- NULL
#   # preprocess the bed file (only need to do once for each data set)
#   result <-paste(".",args[2],paste(args[4],toString(".rds"), sep = ""),sep="//")
#   if (file.exists(result)) {
#     file.remove(result)
#     print(paste("File", result, "deleted."))
#   }
#   result <-paste(".",args[2],paste(args[4],toString(".bk"), sep = ""),sep="//")
#   if (file.exists(result)) {
#     file.remove(result)
#     print(paste("File", result, "deleted."))
#   }
#   
#   
#   result <-paste(".",args[2],paste(args[4],toString(".bed"), sep = ""),sep="//")
#   
#   snp_readBed(result)
#   # now attach the genotype object
#   result <-paste(".",args[2],paste(args[4],toString(".rds"), sep = ""),sep="//")
#   
#   obj.bigSNP <- snp_attach(result)
#   
#   # extract the SNP information from the genotype
#   map <- obj.bigSNP$map[-3]
#   
#   names(map) <- c("chr", "rsid", "pos", "a1", "a0")
#   
#   # perform SNP matching
#   info_snp <- snp_match(sumstats, map)
#   help(snp_match)
#   info_snp
#   # Assign the genotype to a variable for easier downstream analysis
#   genotype <- obj.bigSNP$genotypes
#   # Rename the data structures
#   CHR <- map$chr
#   POS <- map$pos
#   # get the CM information from 1000 Genome
#   # will download the 1000G file to the current directory (".")
#   #help(snp_asGeneticPos)
#   POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
#   
#   
#   
#  
#  check <-TRUE
#   for (chr in 1:22) {
#     # Extract SNPs that are included in the chromosome
#     
#     ind.chr <- which(info_snp$chr == chr)
#     print(length(ind.chr))
#     ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
#     ind.chr2
#     
#     if (length(ind.chr2) == 0) {
#        
#       next  
#     }
#     else{
#        
#       corr0 <- snp_cor(
#       genotype,
#       ind.col = ind.chr2,
#       ncores = NCORES,
#       infos.pos = POS2[ind.chr2],
#       #size = 200,
#       #thr_r2=0.1,
#       #alpha = 1
#       
#       size = as.numeric(args[6]),
#       alpha = as.numeric(args[7]),
#       
#       thr_r2=as.numeric(args[8]),
#     )
#     if (check==TRUE) {
#       check <-FALSE
#       #print("FUCK")
#       ld <- Matrix::colSums(corr0^2)
#       corr <- as_SFBM(corr0, tmp)
#     } else {
#       ld <- c(ld, Matrix::colSums(corr0^2))
#       corr$add_columns(corr0, nrow(corr))
#     }
# 
#     }
#    
#     }
#   
#   
#   # We assume the fam order is the same across different chromosomes
#   fam.order <- as.data.table(obj.bigSNP$fam)
#   # Rename fam order
#   setnames(fam.order,
#            c("family.ID", "sample.ID"),
#            c("FID", "IID"))
#   
#   df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
#   
#   length(df_beta$beta) 
#   length(ld)
#   help(snp_ldsc)
#   ldsc <- snp_ldsc(ld, 
#                    length(ld), 
#                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
#                    sample_size = df_beta$n_eff, 
#                    blocks = NULL)
#   h2_est <- ldsc[["h2"]]
#   h2_est
#   
#   result <-paste(".",args[2],"ldpred_h2_hapmap.txt",sep="//")
#   if (file.exists(result)) {
#     file.remove(result)
#     print(paste("File", result, "deleted."))
#   }
#   write.table(h2_est, file = result, col.names = FALSE)
#   
#   result <-paste(".",args[2],"ldpred_h2_variants.txt",sep="//")
#     
#   if (file.exists(result)) {
#     file.remove(result)
#     print(paste("File", result, "deleted."))
#   }
#   write.table(length(ld), file = result, col.names = FALSE)
#   
#     
#     
#   help(snp_ldsc)
#   
#   # Here we have to specify the p and h2_seq
#   
#   p_seq <- signif(seq_log(1e-4, 1, length.out = 5), 2)
#   h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
#   
#   grid.param <-
#     expand.grid(p = p_seq,
#                 h2 = h2_seq,
#                 sparse = c(FALSE, TRUE))
#   
#   # Save the grid paramters.
#   gridparamters <- grid.param[,c("p", "h2", "sparse")]
#      result <-paste(".",args[2],paste(args[3],toString(".ldpred_grid_parameters"), sep = ""),sep="//")
#   write.table(gridparamters, file = result, row.names = FALSE,sep=",", quote = FALSE)
#    
#   
#   
#   # Get adjusted beta from grid model
#   beta_grid <-
#     snp_ldpred2_grid(corr, df_beta, grid.param,as.numeric(args[9]),as.numeric(args[10]))
#   
#    
#   result <-paste(".",args[2],paste(args[3],toString(".ldpred_grid_betas"), sep = ""),sep="//")
#   write.table(beta_grid, file = result, row.names = FALSE,sep=",", quote = FALSE)
#   newgwas <- info_snp[,c("rsid.ss", "a0", "beta")]
#  
#   result <-paste(".",args[2],paste(args[3],toString(".ldpred_grid_gwas"), sep = ""),sep="//")
#   write.table(newgwas, file = result, row.names = FALSE,sep=",", quote = FALSE)
#   
#   #print(ldsc)
#   #exit(0)
#   
#   
# }
# if (args[5]=="3"){
#   
#   library(bigsnpr)
#   options(bigstatsr.check.parallel.blas = FALSE)
#   options(default.nproc.blas = NULL)
#   library(data.table)
#   library(magrittr)
#   result <-paste(".",args[2],paste(args[3],toString(".PHENO"), sep = ""),sep="//")
#   phenotype <- fread(result)
#   result <-paste(".",args[2],paste(args[3],toString(".cov"), sep = ""),sep="//")
#   covariate <- fread(result)
#   result <-paste(".",args[2],paste(args[3],toString(".eigenvec"), sep = ""),sep="//")
#   pcs <- fread(result)
#   # rename columns
#   colnames(pcs) <- c("FID","IID", paste0("PC",1:as.numeric(args[9])))
#   # generate required table
#   pheno <- merge(phenotype, covariate) %>%
#     merge(., pcs)
#   info <- readRDS(runonce::download_file(
#     "https://ndownloader.figshare.com/files/25503788",
#     fname = "map_hm3_ldpred2.rds"))
#   # Read in the summary statistic file
#   result <-paste(".",args[1],paste(args[1],toString(".txt"), sep = ""),sep="//")
#   
#   sumstats <- bigreadr::fread2(result) 
#   # LDpred 2 require the header to follow the exact naming
#   names(sumstats) <-
#     c("chr",
#       "pos",
#       "rsid",
#       "a1",
#       "a0",
#       "n_eff",
#       "beta_se",
#       "p",
#       "BETA",
#       "INFO",
#       "MAF")
#   # Transform the OR into log(OR)
#   sumstats$beta <- sumstats$BETA
#   # Filter out hapmap SNPs
#   # Turn off this line to ensure that all the SNPs from
#   # the sumstats are included.
#   #sumstats <- sumstats[sumstats$rsid%in% info$rsid,]
#   
#   # Get maximum amount of cores
#   NCORES <- nb_cores()
#   # Open a temporary file
#   
#   result <-paste(".",args[2],"tmp-data",sep="//")
# 
#   if (dir.exists(result)) {
#     # Delete the directory and its contents
#     
#     system(paste("rm -r", shQuote(result)))
#     print(paste("Directory", result, "deleted."))
#   }
#   tmp <- tempfile(tmpdir = result)
#   on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
#     
#     
#   # Initialize variables for storing the LD score and LD matrix
#   corr <- NULL
#   ld <- NULL
#   # We want to know the ordering of samples in the bed file 
#   fam.order <- NULL
#   # preprocess the bed file (only need to do once for each data set)
#   result <-paste(".",args[2],paste(args[4],toString(".rds"), sep = ""),sep="//")
#   if (file.exists(result)) {
#     file.remove(result)
#     print(paste("File", result, "deleted."))
#   }
#   result <-paste(".",args[2],paste(args[4],toString(".bk"), sep = ""),sep="//")
#   if (file.exists(result)) {
#     file.remove(result)
#     print(paste("File", result, "deleted."))
#   }
#   
#   
#   result <-paste(".",args[2],paste(args[4],toString(".bed"), sep = ""),sep="//")
#   
#   snp_readBed(result)
#   # now attach the genotype object
#   result <-paste(".",args[2],paste(args[4],toString(".rds"), sep = ""),sep="//")
#   
#   obj.bigSNP <- snp_attach(result)
#   
#   # extract the SNP information from the genotype
#   map <- obj.bigSNP$map[-3]
#   
#   names(map) <- c("chr", "rsid", "pos", "a1", "a0")
#   
#   # perform SNP matching
#   info_snp <- snp_match(sumstats, map)
#   help(snp_match)
#   info_snp
#   # Assign the genotype to a variable for easier downstream analysis
#   genotype <- obj.bigSNP$genotypes
#   # Rename the data structures
#   CHR <- map$chr
#   POS <- map$pos
#   # get the CM information from 1000 Genome
#   # will download the 1000G file to the current directory (".")
#   help(snp_asGeneticPos)
#   POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
#   
#   
# 
#     check <-TRUE
#   for (chr in 1:22) {
#     # Extract SNPs that are included in the chromosome
#     
#     ind.chr <- which(info_snp$chr == chr)
#     print(length(ind.chr))
#     ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
#     
#     print(length(ind.chr2))
#     
#     if (length(ind.chr2) == 0) {
#        
#       next  
#     }
#     else{
#        
#       corr0 <- snp_cor(
#       genotype,
#       ind.col = ind.chr,
#       ncores = NCORES,
#       infos.pos = POS2[ind.chr2],
#       #size = 200,
#       #thr_r2=0.1,
#       #alpha = 1
#       
#       size = as.numeric(args[6]),
#       alpha = as.numeric(args[7]),
#       thr_r2=as.numeric(args[8]),
#     )
#     if (check==TRUE) {
#       check <-FALSE
#      
#       ld <- Matrix::colSums(corr0^2)
#       help(as_SFBM)
#       corr <- as_SFBM(corr0, tmp)
#     } else {
#       ld <- c(ld, Matrix::colSums(corr0^2))
#       corr$add_columns(corr0, nrow(corr))
#     }
# 
#     }
#    
#     }
#   
#   # We assume the fam order is the same across different chromosomes
#   fam.order <- as.data.table(obj.bigSNP$fam)
#   # Rename fam order
#   setnames(fam.order,
#            c("family.ID", "sample.ID"),
#            c("FID", "IID"))
#   
#   df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
#   
#   length(df_beta$beta) 
#   length(ld)
#   help(snp_ldsc)
#   ldsc <- snp_ldsc(ld, 
#                    length(ld), 
#                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
#                    sample_size = df_beta$n_eff, 
#                    blocks = NULL)
#   h2_est <- ldsc[["h2"]]
#   # Here we have to specify the p and h2_seq
#   
#   p_seq <- signif(seq_log(1e-4, 1, length.out = 5), 2)
#   h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
#   
#   grid.param <-
#     expand.grid(p = p_seq,
#                 h2 = h2_seq,
#                 sparse = c(FALSE, TRUE))
#   
#   # Save the grid paramters.
#   gridparamters <- grid.param[,c("p", "h2", "sparse")]
#   
#     
#   result <-paste(".",args[2],paste(args[3],toString(".ldpred_grid_parameters"), sep = ""),sep="//")
#   write.table(gridparamters, file = result, row.names = FALSE,sep=",", quote = FALSE)
#   
#   
#   
#   # Get adjusted beta from grid model
#   beta_grid <-
#     snp_ldpred2_grid(corr, df_beta, grid.param,as.numeric(args[9]),as.numeric(args[10]))
#   
#   
#   result <-paste(".",args[2],paste(args[3],toString(".ldpred_grid_betas"), sep = ""),sep="//")
#   write.table(beta_grid, file = result, row.names = FALSE,sep=",", quote = FALSE)
#   newgwas <- info_snp[,c("rsid.ss", "a0", "beta")]
#   
#   result <-paste(".",args[2],paste(args[3],toString(".ldpred_grid_gwas"), sep = ""),sep="//")
#   write.table(newgwas, file = result, row.names = FALSE,sep=",", quote = FALSE)
#   
#   
#   
#   
#   
#   
#   
#   result <-paste(".",args[2],"ldpred_h2_full.txt",sep="//")
#   if (file.exists(result)) {
#     file.remove(result)
#     print(paste("File", result, "deleted."))
#   }
#   write.table(h2_est, file = result, col.names = FALSE)
#   
#   result <-paste(".",args[2],"ldpred_h2_variants.txt",sep="//")
#     
#   if (file.exists(result)) {
#     file.remove(result)
#     print(paste("File", result, "deleted."))
#   }
#   write.table(length(ld), file = result, col.names = FALSE)
#   
# }
# 
# 
# 
# ```

# In[2]:


from operator import index
import pandas as pd
import numpy as np
import os
import subprocess
import sys
import pandas as pd
import statsmodels.api as sm
import pandas as pd
from sklearn.metrics import roc_auc_score, confusion_matrix
from statsmodels.stats.contingency_tables import mcnemar

def create_directory(directory):
    """Function to create a directory if it doesn't exist."""
    if not os.path.exists(directory):  # Checking if the directory doesn't exist
        os.makedirs(directory)  # Creating the directory if it doesn't exist
    return directory  # Returning the created or existing directory

 
foldnumber = sys.argv[2]
#foldnumber = "0"  # Setting 'foldnumber' to "0"

folddirec = filedirec + os.sep + "Fold_" + foldnumber  # Creating a directory path for the specific fold
trainfilename = "train_data"  # Setting the name of the training data file
newtrainfilename = "train_data.QC"  # Setting the name of the new training data file

testfilename = "test_data"  # Setting the name of the test data file
newtestfilename = "test_data.QC"  # Setting the name of the new test data file

# Number of PCA to be included as a covariate.
numberofpca = ["6"]  # Setting the number of PCA components to be included

# Clumping parameters.
clump_p1 = [1]  # List containing clump parameter 'p1'
clump_r2 = [0.1]  # List containing clump parameter 'r2'
clump_kb = [200]  # List containing clump parameter 'kb'

# Pruning parameters.
p_window_size = [200]  # List containing pruning parameter 'window_size'
p_slide_size = [50]  # List containing pruning parameter 'slide_size'
p_LD_threshold = [0.25]  # List containing pruning parameter 'LD_threshold'

# Kindly note that the number of p-values to be considered varies, and the actual p-value depends on the dataset as well.
# We will specify the range list here.


minimumpvalue = 10  # Minimum p-value in exponent
numberofintervals = 20  # Number of intervals to be considered
allpvalues = np.logspace(-minimumpvalue, 0, numberofintervals, endpoint=True)  # Generating an array of logarithmically spaced p-values
count = 1
with open(folddirec + os.sep + 'range_list', 'w') as file:
    for value in allpvalues:
        file.write(f'pv_{value} 0 {value}\n')  # Writing range information to the 'range_list' file
        count = count + 1

pvaluefile = folddirec + os.sep + 'range_list'

# Initializing an empty DataFrame with specified column names
# Initializing an empty DataFrame with specified column names
prs_result = pd.DataFrame(columns=["clump_p1", "clump_r2", "clump_kb", "p_window_size", "p_slide_size", "p_LD_threshold",
                                   "pvalue", "numberofpca","burn_in","num_iter","p","h2","sparse","h2model","numberofvariants","Train_pure_prs", "Train_null_model", "Train_best_model",
                                   "Test_pure_prs", "Test_null_model", "Test_best_model"])


# ### Define Helper Functions
# 
# 1. **Perform Clumping and Pruning**
# 2. **Calculate PCA Using Plink**
# 3. **Fit Binary Phenotype and Save Results**
# 4. **Fit Continuous Phenotype and Save Results**
# 

# In[3]:


import os
import subprocess
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import explained_variance_score

prs_result = pd.DataFrame()
def perform_clumping_and_pruning_on_individual_data(traindirec, newtrainfilename,numberofpca, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
    
    command = [
    "./plink",
    "--bfile", traindirec+os.sep+newtrainfilename,
    "--indep-pairwise", p1_val, p2_val, p3_val,
    "--out", traindirec+os.sep+trainfilename
    ]
    subprocess.run(command)
    # First perform pruning and then clumping and the pruning.

    command = [
    "./plink",
    "--bfile", traindirec+os.sep+newtrainfilename,
    "--clump-p1", c1_val,
    "--extract", traindirec+os.sep+trainfilename+".prune.in",
    "--clump-r2", c2_val,
    "--clump-kb", c3_val,
    "--clump", filedirec+os.sep+filedirec+".txt",
    "--clump-snp-field", "SNP",
    "--clump-field", "P",
    "--out", traindirec+os.sep+trainfilename
    ]    
    subprocess.run(command)

    # Extract the valid SNPs from th clumped file.
    # For windows download gwak for linux awk commmand is sufficient.
    ### For windows require GWAK.
    ### https://sourceforge.net/projects/gnuwin32/
    ##3 Get it and place it in the same direc.
    #os.system("gawk "+"\""+"NR!=1{print $3}"+"\"  "+ traindirec+os.sep+trainfilename+".clumped >  "+traindirec+os.sep+trainfilename+".valid.snp")
    #print("gawk "+"\""+"NR!=1{print $3}"+"\"  "+ traindirec+os.sep+trainfilename+".clumped >  "+traindirec+os.sep+trainfilename+".valid.snp")

    #Linux:
    command = f"awk 'NR!=1{{print $3}}' {traindirec}{os.sep}{trainfilename}.clumped > {traindirec}{os.sep}{trainfilename}.valid.snp"
    os.system(command)
    
    
    command = [
    "./plink",
    "--make-bed",
    "--bfile", traindirec+os.sep+newtrainfilename,
    "--indep-pairwise", p1_val, p2_val, p3_val,
    "--extract", traindirec+os.sep+trainfilename+".valid.snp",
    "--out", traindirec+os.sep+newtrainfilename+".clumped.pruned"
    ]
    subprocess.run(command)
    
    command = [
    "./plink",
    "--make-bed",
    "--bfile", traindirec+os.sep+testfilename,
    "--indep-pairwise", p1_val, p2_val, p3_val,
    "--extract", traindirec+os.sep+trainfilename+".valid.snp",
    "--out", traindirec+os.sep+testfilename+".clumped.pruned"
    ]
    subprocess.run(command)    
    
    
 
def calculate_pca_for_traindata_testdata_for_clumped_pruned_snps(traindirec, newtrainfilename,p):
    
    # Calculate the PRS for the test data using the same set of SNPs and also calculate the PCA.


    # Also extract the PCA at this point.
    # PCA are calculated afer clumping and pruining.
    command = [
        "./plink",
        "--bfile", folddirec+os.sep+testfilename+".clumped.pruned",
        # Select the final variants after clumping and pruning.
        "--extract", traindirec+os.sep+trainfilename+".valid.snp",
        "--pca", p,
        "--out", folddirec+os.sep+testfilename
    ]
    subprocess.run(command)


    command = [
    "./plink",
        "--bfile", traindirec+os.sep+newtrainfilename+".clumped.pruned",
        # Select the final variants after clumping and pruning.        
        "--extract", traindirec+os.sep+trainfilename+".valid.snp",
        "--pca", p,
        "--out", traindirec+os.sep+trainfilename
    ]
    subprocess.run(command)


# This function fit the binary model on the PRS.
def fit_binary_phenotype_on_PRS(traindirec, newtrainfilename,h2model, p,b,iterr, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val, Name, pvaluefile,heritability,numberofvariants,tempp,temph2,tempsparse,count_grid):
    threshold_values = allpvalues

    # Merge the covariates, pca and phenotypes.
    tempphenotype_train = pd.read_table(traindirec+os.sep+newtrainfilename+".clumped.pruned"+".fam", sep="\s+",header=None)
    phenotype_train = pd.DataFrame()
    phenotype_train["Phenotype"] = tempphenotype_train[5].values
    pcs_train = pd.read_table(traindirec+os.sep+trainfilename+".eigenvec", sep="\s+",header=None, names=["FID", "IID"] + [f"PC{str(i)}" for i in range(1, int(p)+1)])
    covariate_train = pd.read_table(traindirec+os.sep+trainfilename+".cov",sep="\s+")
    covariate_train.fillna(0, inplace=True)
    covariate_train = covariate_train[covariate_train["FID"].isin(pcs_train["FID"].values) & covariate_train["IID"].isin(pcs_train["IID"].values)]
    covariate_train['FID'] = covariate_train['FID'].astype(str)
    pcs_train['FID'] = pcs_train['FID'].astype(str)
    covariate_train['IID'] = covariate_train['IID'].astype(str)
    pcs_train['IID'] = pcs_train['IID'].astype(str)
    covandpcs_train = pd.merge(covariate_train, pcs_train, on=["FID","IID"])
    covandpcs_train.fillna(0, inplace=True)


    ## Scale the covariates!
    from sklearn.preprocessing import MinMaxScaler
    from sklearn.metrics import explained_variance_score
    scaler = MinMaxScaler()
    normalized_values_train = scaler.fit_transform(covandpcs_train.iloc[:, 2:])
    #covandpcs_train.iloc[:, 2:] = normalized_values_test 
    
    
    tempphenotype_test = pd.read_table(traindirec+os.sep+testfilename+".clumped.pruned"+".fam", sep="\s+",header=None)
    phenotype_test= pd.DataFrame()
    phenotype_test["Phenotype"] = tempphenotype_test[5].values
    pcs_test = pd.read_table(traindirec+os.sep+testfilename+".eigenvec", sep="\s+",header=None, names=["FID", "IID"] + [f"PC{str(i)}" for i in range(1, int(p)+1)])
    covariate_test = pd.read_table(traindirec+os.sep+testfilename+".cov",sep="\s+")
    covariate_test.fillna(0, inplace=True)
    covariate_test = covariate_test[covariate_test["FID"].isin(pcs_test["FID"].values) & covariate_test["IID"].isin(pcs_test["IID"].values)]
    covariate_test['FID'] = covariate_test['FID'].astype(str)
    pcs_test['FID'] = pcs_test['FID'].astype(str)
    covariate_test['IID'] = covariate_test['IID'].astype(str)
    pcs_test['IID'] = pcs_test['IID'].astype(str)
    covandpcs_test = pd.merge(covariate_test, pcs_test, on=["FID","IID"])
    covandpcs_test.fillna(0, inplace=True)
    normalized_values_test  = scaler.transform(covandpcs_test.iloc[:, 2:])
    #covandpcs_test.iloc[:, 2:] = normalized_values_test     
    
    
    
    
    tempalphas = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    l1weights = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

    tempalphas = [0.1]
    l1weights = [0.1]

    phenotype_train["Phenotype"] = phenotype_train["Phenotype"].replace({1: 0, 2: 1}) 
    phenotype_test["Phenotype"] = phenotype_test["Phenotype"].replace({1: 0, 2: 1})
      
    for tempalpha in tempalphas:
        for l1weight in l1weights:

            
            try:
                null_model =  sm.Logit(phenotype_train["Phenotype"], sm.add_constant(covandpcs_train.iloc[:, 2:])).fit_regularized(alpha=tempalpha, L1_wt=l1weight)
                #null_model =  sm.Logit(phenotype_train["Phenotype"], sm.add_constant(covandpcs_train.iloc[:, 2:])).fit()
            
            except:
                print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
                continue

            train_null_predicted = null_model.predict(sm.add_constant(covandpcs_train.iloc[:, 2:]))
            
            from sklearn.metrics import roc_auc_score, confusion_matrix
            from sklearn.metrics import r2_score
            
            test_null_predicted = null_model.predict(sm.add_constant(covandpcs_test.iloc[:, 2:]))
            
           
            
            global prs_result 
            for i in threshold_values:
                try:
                    prs_train = pd.read_table(traindirec+os.sep+Name+os.sep+"train_data.pv_"+f"{i}.profile", sep="\s+", usecols=["FID", "IID", "SCORE"])
                except:
                    
                    continue

                prs_train['FID'] = prs_train['FID'].astype(str)
                prs_train['IID'] = prs_train['IID'].astype(str)
                try:
                    prs_test = pd.read_table(traindirec+os.sep+Name+os.sep+"test_data.pv_"+f"{i}.profile", sep="\s+", usecols=["FID", "IID", "SCORE"])
                except:
                    continue
                prs_test['FID'] = prs_test['FID'].astype(str)
                prs_test['IID'] = prs_test['IID'].astype(str)
                pheno_prs_train = pd.merge(covandpcs_train, prs_train, on=["FID", "IID"])
                pheno_prs_test = pd.merge(covandpcs_test, prs_test, on=["FID", "IID"])
        
                try:
                    model = sm.Logit(phenotype_train["Phenotype"], sm.add_constant(pheno_prs_train.iloc[:, 2:])).fit_regularized(alpha=tempalpha, L1_wt=l1weight)
                    #model = sm.Logit(phenotype_train["Phenotype"], sm.add_constant(pheno_prs_train.iloc[:, 2:])).fit()
                
                except:
                    print("Did not work!")
                    continue


                
                train_best_predicted = model.predict(sm.add_constant(pheno_prs_train.iloc[:, 2:]))    
 

                test_best_predicted = model.predict(sm.add_constant(pheno_prs_test.iloc[:, 2:])) 
 
        
                from sklearn.metrics import roc_auc_score, confusion_matrix

                prs_result = prs_result._append({
                    "clump_p1": c1_val,
                    "clump_r2": c2_val,
                    "clump_kb": c3_val,
                    "p_window_size": p1_val,
                    "p_slide_size": p2_val,
                    "p_LD_threshold": p3_val,
                    "pvalue": i,
                    "numberofpca":p, 

                    "tempalpha":str(tempalpha),
                    "l1weight":str(l1weight),
                      
                     
                    "numberofvariants": numberofvariants,
       
                    "heritability_model":h2model,
                    "gridh2":heritability,
                    "burn_in":b,
                    "num_iter":iterr,
                    "sparse":tempsparse,
                    "unique_h2":count_grid,
                    "grid_pvalue":tempp,                    
                    
                    "Train_pure_prs":roc_auc_score(phenotype_train["Phenotype"].values,prs_train['SCORE'].values),
                    "Train_null_model":roc_auc_score(phenotype_train["Phenotype"].values,train_null_predicted.values),
                    "Train_best_model":roc_auc_score(phenotype_train["Phenotype"].values,train_best_predicted.values),
                    
                    "Test_pure_prs":roc_auc_score(phenotype_test["Phenotype"].values,prs_test['SCORE'].values),
                    "Test_null_model":roc_auc_score(phenotype_test["Phenotype"].values,test_null_predicted.values),
                    "Test_best_model":roc_auc_score(phenotype_test["Phenotype"].values,test_best_predicted.values),
                    
                }, ignore_index=True)

          
                prs_result.to_csv(traindirec+os.sep+Name+os.sep+"Results.csv",index=False)
     
    return

# This function fit the binary model on the PRS.
def fit_continous_phenotype_on_PRS(traindirec, newtrainfilename,h2model, p,b,iterr, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val, Name, pvaluefile,heritability,numberofvariants,tempp,temph2,tempsparse,count_grid):
    threshold_values = allpvalues

    # Merge the covariates, pca and phenotypes.
    tempphenotype_train = pd.read_table(traindirec+os.sep+newtrainfilename+".clumped.pruned"+".fam", sep="\s+",header=None)
    phenotype_train = pd.DataFrame()
    phenotype_train["Phenotype"] = tempphenotype_train[5].values
    pcs_train = pd.read_table(traindirec+os.sep+trainfilename+".eigenvec", sep="\s+",header=None, names=["FID", "IID"] + [f"PC{str(i)}" for i in range(1, int(p)+1)])
    covariate_train = pd.read_table(traindirec+os.sep+trainfilename+".cov",sep="\s+")
    covariate_train.fillna(0, inplace=True)
    covariate_train = covariate_train[covariate_train["FID"].isin(pcs_train["FID"].values) & covariate_train["IID"].isin(pcs_train["IID"].values)]
    covariate_train['FID'] = covariate_train['FID'].astype(str)
    pcs_train['FID'] = pcs_train['FID'].astype(str)
    covariate_train['IID'] = covariate_train['IID'].astype(str)
    pcs_train['IID'] = pcs_train['IID'].astype(str)
    covandpcs_train = pd.merge(covariate_train, pcs_train, on=["FID","IID"])
    covandpcs_train.fillna(0, inplace=True)


    ## Scale the covariates!
    from sklearn.preprocessing import MinMaxScaler
    from sklearn.metrics import explained_variance_score
    scaler = MinMaxScaler()
    normalized_values_train = scaler.fit_transform(covandpcs_train.iloc[:, 2:])
    #covandpcs_train.iloc[:, 2:] = normalized_values_test 
    
    tempphenotype_test = pd.read_table(traindirec+os.sep+testfilename+".clumped.pruned"+".fam", sep="\s+",header=None)
    phenotype_test= pd.DataFrame()
    phenotype_test["Phenotype"] = tempphenotype_test[5].values
    pcs_test = pd.read_table(traindirec+os.sep+testfilename+".eigenvec", sep="\s+",header=None, names=["FID", "IID"] + [f"PC{str(i)}" for i in range(1, int(p)+1)])
    covariate_test = pd.read_table(traindirec+os.sep+testfilename+".cov",sep="\s+")
    covariate_test.fillna(0, inplace=True)
    covariate_test = covariate_test[covariate_test["FID"].isin(pcs_test["FID"].values) & covariate_test["IID"].isin(pcs_test["IID"].values)]
    covariate_test['FID'] = covariate_test['FID'].astype(str)
    pcs_test['FID'] = pcs_test['FID'].astype(str)
    covariate_test['IID'] = covariate_test['IID'].astype(str)
    pcs_test['IID'] = pcs_test['IID'].astype(str)
    covandpcs_test = pd.merge(covariate_test, pcs_test, on=["FID","IID"])
    covandpcs_test.fillna(0, inplace=True)
    normalized_values_test  = scaler.transform(covandpcs_test.iloc[:, 2:])
    #covandpcs_test.iloc[:, 2:] = normalized_values_test     
    
    
    
    
    tempalphas = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    l1weights = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

    tempalphas = [0.1]
    l1weights = [0.1]

    #phenotype_train["Phenotype"] = phenotype_train["Phenotype"].replace({1: 0, 2: 1}) 
    #phenotype_test["Phenotype"] = phenotype_test["Phenotype"].replace({1: 0, 2: 1})
      
    for tempalpha in tempalphas:
        for l1weight in l1weights:

            
            try:
                #null_model =  sm.OLS(phenotype_train["Phenotype"], sm.add_constant(covandpcs_train.iloc[:, 2:])).fit_regularized(alpha=tempalpha, L1_wt=l1weight)
                null_model =  sm.OLS(phenotype_train["Phenotype"], sm.add_constant(covandpcs_train.iloc[:, 2:])).fit()
                #null_model =  sm.OLS(phenotype_train["Phenotype"], sm.add_constant(covandpcs_train.iloc[:, 2:])).fit()
            except:
                print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
                continue

            train_null_predicted = null_model.predict(sm.add_constant(covandpcs_train.iloc[:, 2:]))
            
            from sklearn.metrics import roc_auc_score, confusion_matrix
            from sklearn.metrics import r2_score
            
            test_null_predicted = null_model.predict(sm.add_constant(covandpcs_test.iloc[:, 2:]))
            
            
            
            global prs_result 
            for i in threshold_values:
                try:
                    prs_train = pd.read_table(traindirec+os.sep+Name+os.sep+"train_data.pv_"+f"{i}.profile", sep="\s+", usecols=["FID", "IID", "SCORE"])
                except:
                    continue

                prs_train['FID'] = prs_train['FID'].astype(str)
                prs_train['IID'] = prs_train['IID'].astype(str)
                try:
                    prs_test = pd.read_table(traindirec+os.sep+Name+os.sep+"test_data.pv_"+f"{i}.profile", sep="\s+", usecols=["FID", "IID", "SCORE"])
                except:
                    continue
                prs_test['FID'] = prs_test['FID'].astype(str)
                prs_test['IID'] = prs_test['IID'].astype(str)
                pheno_prs_train = pd.merge(covandpcs_train, prs_train, on=["FID", "IID"])
                pheno_prs_test = pd.merge(covandpcs_test, prs_test, on=["FID", "IID"])
        
                try:
                    #model = sm.OLS(phenotype_train["Phenotype"], sm.add_constant(pheno_prs_train.iloc[:, 2:])).fit_regularized(alpha=tempalpha, L1_wt=l1weight)
                    model = sm.OLS(phenotype_train["Phenotype"], sm.add_constant(pheno_prs_train.iloc[:, 2:])).fit()
                
                except:
                    continue


                
                train_best_predicted = model.predict(sm.add_constant(pheno_prs_train.iloc[:, 2:]))    
                test_best_predicted = model.predict(sm.add_constant(pheno_prs_test.iloc[:, 2:])) 
 
        
                from sklearn.metrics import roc_auc_score, confusion_matrix

                prs_result = prs_result._append({
                    "clump_p1": c1_val,
                    "clump_r2": c2_val,
                    "clump_kb": c3_val,
                    "p_window_size": p1_val,
                    "p_slide_size": p2_val,
                    "p_LD_threshold": p3_val,
                    "pvalue": i,
                    "numberofpca":p, 

                    "tempalpha":str(tempalpha),
                    "l1weight":str(l1weight),
                    "numberofvariants": numberofvariants,
       
                    "heritability_model":h2model,
                    "h2":heritability,
                    "burn_in":b,
                    "num_iter":iterr,
                    "sparse":tempsparse,
                    "unique_h2":count_grid,
                    "grid_pvalue":tempp,

                    "Train_pure_prs":explained_variance_score(phenotype_train["Phenotype"],prs_train['SCORE'].values),
                    "Train_null_model":explained_variance_score(phenotype_train["Phenotype"],train_null_predicted),
                    "Train_best_model":explained_variance_score(phenotype_train["Phenotype"],train_best_predicted),
                    
                    "Test_pure_prs":explained_variance_score(phenotype_test["Phenotype"],prs_test['SCORE'].values),
                    "Test_null_model":explained_variance_score(phenotype_test["Phenotype"],test_null_predicted),
                    "Test_best_model":explained_variance_score(phenotype_test["Phenotype"],test_best_predicted),
                    
                }, ignore_index=True)

          
                prs_result.to_csv(traindirec+os.sep+Name+os.sep+"Results.csv",index=False)
     
    return


# ## Execute LDpred-2-Grid

# In[10]:


def transform_ldpred2_grid_data(traindirec, newtrainfilename,models,p,b,iterr, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
    ### First perform clumping on the file and save the clumpled file.
    #perform_clumping_and_pruning_on_individual_data(traindirec, newtrainfilename,p, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile)
   
        
    # Also extract the PCA at this point for both test and training data.
    #calculate_pca_for_traindata_testdata_for_clumped_pruned_snps(traindirec, newtrainfilename,p)

    #Extract p-values from the GWAS file.
    # Command for Linux.
    os.system("awk "+"\'"+"{print $3,$8}"+"\'"+" ./"+filedirec+os.sep+filedirec+".txt >  ./"+traindirec+os.sep+"SNP.pvalue")

    # Define the paths to the files
    files_to_remove = [
        traindirec+os.sep+"train_data.ldpred_grid_parameters",
        traindirec+os.sep+"train_data.ldpred_grid_betas",
        traindirec+os.sep+"train_data.ldpred_grid_gwas",
        traindirec+os.sep+"train_data.ldpred_grid_gwas_final",
        traindirec+os.sep+"ldpred_h2_hapmap.txt",
        traindirec+os.sep+"ldpred_h2_full.txt"
        

    ]

    # Loop through the files and remove them if they exist
    for file_path in files_to_remove:
        if os.path.exists(file_path):
            os.remove(file_path)
            print(f"Removed: {file_path}")
        else:
            print(f"File does not exist: {file_path}")
                  
  
    if models=="LDpred-2_full":
        try:
            os.system("Rscript LDpred2-grid.R "+os.path.join(filedirec)+"  "+traindirec+" "+trainfilename+" "+newtrainfilename+".clumped.pruned"+ " "+"3"+" "+c3_val+" "+c1_val+" "+c2_val+" "+p+" "+b+" "+i)
            heritability = pd.read_csv(traindirec+os.sep+"ldpred_h2_full.txt",sep="\s+",header=None)[1].values[0]
            numberofvariants =  pd.read_csv(traindirec+os.sep+"ldpred_h2_variants.txt",sep="\s+",header=None)[1].values[0] 
        except:
            print("For model ",models,"it did not work! may be h2 is negative.")
            return

    if models=="LDpred-2_hapmap":
        try:
            os.system("Rscript LDpred2-grid.R "+os.path.join(filedirec)+"  "+traindirec+" "+trainfilename+" "+newtrainfilename+".clumped.pruned"+ " "+"2"+" "+c3_val+" "+c1_val+" "+c2_val+" "+p+" "+b+" "+i)
            heritability = pd.read_csv(traindirec+os.sep+"ldpred_h2_hapmap.txt",sep="\s+",header=None)[1].values[0]
            numberofvariants =  pd.read_csv(traindirec+os.sep+"ldpred_h2_variants.txt",sep="\s+",header=None)[1].values[0] 
        except:
            print("For model ",models,"it did not work! may be h2 is negative.")
            return
        

    

            
    # After obtaining betas from the LDpred-2
    # Append the new betas.
    # Read all the betas from the file.
    # First read the paramters.
    gridparameters = pd.read_csv(traindirec+os.sep+"train_data.ldpred_grid_parameters",sep=",")
    allbetas = pd.read_csv(traindirec+os.sep+"train_data.ldpred_grid_betas",sep=",")
    allbetas = allbetas.fillna(0)
    
    
    count = 1
    for index, row in gridparameters.iterrows():
        gwas = pd.read_csv(traindirec+os.sep+"train_data.ldpred_grid_gwas",sep=",")
        # Accessing individual elements in the row
        tempp = row['p']
        temph2 = row['h2']
        tempsparse = row['sparse']
        gwas["newbetas"] = allbetas.iloc[:, index].values
        
        if check_phenotype_is_binary_or_continous(filedirec)=="Binary":
            gwas["newbetas"] = np.exp(gwas["newbetas"])
        else:
            pass

        print(gwas.head())
        gwas.rename(columns={'rsid.ss': 'SNP', 'a0': 'A1','newbetas':'BETA'}, inplace=True)
        gwas[["SNP","A1","BETA"]].to_csv(traindirec+os.sep+"train_data.ldpred_grid_gwas_final",sep="\t",index=False)
         
        
        count=count+1
        
        print(gwas.head())
       




        command = [
            "./plink",
             "--bfile", traindirec+os.sep+newtrainfilename+".clumped.pruned",
            ### SNP column = 1, Effect allele column 2 = 4, Effect column=4
            "--score", traindirec+os.sep+"train_data.ldpred_grid_gwas_final", "1", "2", "3", "header",
            # This q-score range is not required as we already used the p values in the LDpred2-grid
            "--q-score-range", traindirec+os.sep+"range_list",traindirec+os.sep+"SNP.pvalue",
            "--extract", traindirec+os.sep+trainfilename+".valid.snp",
            "--out", traindirec+os.sep+Name+os.sep+trainfilename
        ]
        #exit(0)
        subprocess.run(command)

 

        command = [
            "./plink",
            "--bfile", folddirec+os.sep+testfilename,
            ### SNP column = 3, Effect allele column 1 = 4, Beta column=12
            "--score", traindirec+os.sep+"train_data.ldpred_grid_gwas_final", "1", "2", "3", "header",
            "--q-score-range", traindirec+os.sep+"range_list",traindirec+os.sep+"SNP.pvalue",
            "--extract", traindirec+os.sep+trainfilename+".valid.snp",
            "--out", folddirec+os.sep+Name+os.sep+testfilename
        ]
        subprocess.run(command)
        # Simulate a KeyboardInterrupt in Jupyter Notebook
        
    
        if check_phenotype_is_binary_or_continous(filedirec)=="Binary":
            print("Binary Phenotype!")
            fit_binary_phenotype_on_PRS(traindirec, newtrainfilename,model, p,b,i, str(p1_val), str(p2_val), str(p3_val), str(c1_val), str(c2_val), str(c3_val), "ldpred2_grid", pvaluefile,heritability,numberofvariants,tempp,temph2,tempsparse,count)
        else:
            print("Continous Phenotype!")
            fit_continous_phenotype_on_PRS(traindirec, newtrainfilename,model, p,b,i, str(p1_val), str(p2_val), str(p3_val), str(c1_val), str(c2_val), str(c3_val), "ldpred2_grid", pvaluefile,heritability,numberofvariants,tempp,temph2,tempsparse,count)
    return
     
     
h2models = ["LDpred-2_full","LDpred-2_hapmap"]
#h2models = ["LDpred-2_hapmap"]


burn_in = ['200']
num_iter = ['10']
result_directory = "ldpred2_grid"

# Nested loops to iterate over different parameter values
create_directory(folddirec+os.sep+result_directory)
for p1_val in p_window_size:
 for p2_val in p_slide_size: 
  for p3_val in p_LD_threshold:
   for c1_val in clump_p1:
    for c2_val in clump_r2:
     for c3_val in clump_kb:
      for p in numberofpca:
       for model in  h2models:
        for b in burn_in:
         for i in num_iter:
          transform_ldpred2_grid_data(folddirec, newtrainfilename,model, p,b,i, str(p1_val), str(p2_val), str(p3_val), str(c1_val), str(c2_val), str(c3_val), result_directory, pvaluefile)



# ### Repeat the process for each fold.
# 
# Change the `foldnumber` variable.
# 
# ```python
# #foldnumber = sys.argv[1]
# foldnumber = "0"  # Setting 'foldnumber' to "0"
# ```
# 
# Or uncomment the following line:
# ```python
# # foldnumber = sys.argv[1]
# python LDpred-2-Grid.py 0
# python LDpred-2-Grid.py 1
# python LDpred-2-Grid.py 2
# python LDpred-2-Grid.py 3
# python LDpred-2-Grid.py 4
# ```
# 
# The following files should exist after the execution:
# 
# 1. `SampleData1/Fold_0/ldpred2_grid/Results.csv`
# 2. `SampleData1/Fold_1/ldpred2_grid/Results.csv`
# 3. `SampleData1/Fold_2/ldpred2_grid/Results.csv`
# 4. `SampleData1/Fold_3/ldpred2_grid/Results.csv`
# 5. `SampleData1/Fold_4/ldpred2_grid/Results.csv`
# 

# ### Check the results file for each fold.

# In[4]:


import os


# List of file names to check for existence
f = [
    "./"+filedirec+"/Fold_0"+os.sep+result_directory+"Results.csv",
    "./"+filedirec+"/Fold_1"+os.sep+result_directory+"Results.csv",
    "./"+filedirec+"/Fold_2"+os.sep+result_directory+"Results.csv",
    "./"+filedirec+"/Fold_3"+os.sep+result_directory+"Results.csv",
    "./"+filedirec+"/Fold_4"+os.sep+result_directory+"Results.csv",
]

 

# Loop through each file name in the list
for loop in range(0,5):
    # Check if the file exists in the specified directory for the given fold
    if os.path.exists(filedirec+os.sep+"Fold_"+str(loop)+os.sep+result_directory+os.sep+"Results.csv"):
        temp = pd.read_csv(filedirec+os.sep+"Fold_"+str(loop)+os.sep+result_directory+os.sep+"Results.csv")
        print("Fold_",loop, "Yes, the file exists.")
        #print(temp.head())
        print("Number of P-values processed: ",len(temp))
        # Print a message indicating that the file exists
    
    else:
        # Print a message indicating that the file does not exist
        print("Fold_",loop, "No, the file does not exist.")



# ### Sum the results for each fold.

# In[5]:


print("We have to ensure when we sum the entries across all Folds, the same rows are merged!")

def sum_and_average_columns(data_frames):
    """Sum and average numerical columns across multiple DataFrames, and keep non-numerical columns unchanged."""
    # Initialize DataFrame to store the summed results for numerical columns
    summed_df = pd.DataFrame()
    non_numerical_df = pd.DataFrame()
    
    for df in data_frames:
        # Identify numerical and non-numerical columns
        numerical_cols = df.select_dtypes(include=[np.number]).columns
        non_numerical_cols = df.select_dtypes(exclude=[np.number]).columns
        
        # Sum numerical columns
        if summed_df.empty:
            summed_df = pd.DataFrame(0, index=range(len(df)), columns=numerical_cols)
        
        summed_df[numerical_cols] = summed_df[numerical_cols].add(df[numerical_cols], fill_value=0)
        
        # Keep non-numerical columns (take the first non-numerical entry for each column)
        if non_numerical_df.empty:
            non_numerical_df = df[non_numerical_cols]
        else:
            non_numerical_df[non_numerical_cols] = non_numerical_df[non_numerical_cols].combine_first(df[non_numerical_cols])
    
    # Divide the summed values by the number of dataframes to get the average
    averaged_df = summed_df / len(data_frames)
    
    # Combine numerical and non-numerical DataFrames
    result_df = pd.concat([averaged_df, non_numerical_df], axis=1)
    
    return result_df

from functools import reduce

import os
import pandas as pd
from functools import reduce

def find_common_rows(allfoldsframe):
    # Define the performance columns that need to be excluded
    performance_columns = [
        'Train_null_model', 'Train_pure_prs', 'Train_best_model',
        'Test_pure_prs', 'Test_null_model', 'Test_best_model'
    ]
    
    # Based on these columns a unique combination of hyperparameters can be defined for all PRS Tools. 
    
    important_columns = [
        'clump_p1',
        'clump_r2',
        'clump_kb',
        'p_window_size',
        'p_slide_size',
        'p_LD_threshold',
        'pvalue',
        'referencepanel',
        'PRSice-2_Model',
        'effectsizes',
        'h2model',
        'lambda',
        'delta',
        'model',
        'numberofpca',
        'tempalpha',
        'l1weight',
        'LDpred-funct-bins',
        'heritability_model',
         
       
    ]
    # Function to remove performance columns from a DataFrame
    def drop_performance_columns(df):
        return df.drop(columns=performance_columns, errors='ignore')
    
    def get_important_columns(df ):
        existing_columns = [col for col in important_columns if col in df.columns]
        if existing_columns:
            return df[existing_columns].copy()
        else:
            return pd.DataFrame()

    # Drop performance columns from all DataFrames in the list
    allfoldsframe_dropped = [drop_performance_columns(df) for df in allfoldsframe]
    
    # Get the important columns.
    allfoldsframe_dropped = [get_important_columns(df) for df in allfoldsframe_dropped]    
    
    # Iteratively find common rows and track unique and common rows
    common_rows = allfoldsframe_dropped[0]
    for i in range(1, len(allfoldsframe_dropped)):
        # Get the next DataFrame
        next_df = allfoldsframe_dropped[i]

        # Count unique rows in the current DataFrame and the next DataFrame
        unique_in_common = common_rows.shape[0]
        unique_in_next = next_df.shape[0]

        # Find common rows between the current common_rows and the next DataFrame
        common_rows = pd.merge(common_rows, next_df, how='inner')
    
        # Count the common rows after merging
        common_count = common_rows.shape[0]

        # Print the unique and common row counts
        print(f"Iteration {i}:")
        print(f"Unique rows in current common DataFrame: {unique_in_common}")
        print(f"Unique rows in next DataFrame: {unique_in_next}")
        print(f"Common rows after merge: {common_count}\n")
    # Now that we have the common rows, extract these from the original DataFrames
 
    extracted_common_rows_frames = []
    for original_df in allfoldsframe:
        # Merge the common rows with the original DataFrame, keeping only the rows that match the common rows
        extracted_common_rows = pd.merge(common_rows, original_df, how='inner', on=common_rows.columns.tolist())
        
        # Add the DataFrame with the extracted common rows to the list
        extracted_common_rows_frames.append(extracted_common_rows)

    # Print the number of rows in the common DataFrames
    for i, df in enumerate(extracted_common_rows_frames):
        print(f"DataFrame {i + 1} with extracted common rows has {df.shape[0]} rows.")

    # Return the list of DataFrames with extracted common rows
    return extracted_common_rows_frames


# Example usage (assuming allfoldsframe is populated as shown earlier):
allfoldsframe = []

# Loop through each file name in the list
for loop in range(0, 5):
    # Check if the file exists in the specified directory for the given fold
    file_path = os.path.join(filedirec, "Fold_" + str(loop), result_directory, "Results.csv")
    if os.path.exists(file_path):
        allfoldsframe.append(pd.read_csv(file_path))
        # Print a message indicating that the file exists
        print("Fold_", loop, "Yes, the file exists.")
    else:
        # Print a message indicating that the file does not exist
        print("Fold_", loop, "No, the file does not exist.")

# Find the common rows across all folds and return the list of extracted common rows
extracted_common_rows_list = find_common_rows(allfoldsframe)
 
# Sum the values column-wise
# For string values, do not sum it the values are going to be the same for each fold.
# Only sum the numeric values.

divided_result = sum_and_average_columns(extracted_common_rows_list)
  
print(divided_result)

 


# ## Results
# 
# ### 1. **Reporting Based on Best Training Performance:**
#    - One can report the results based on the best performance of the training data. For example, if for a specific combination of hyperparameters, the training performance is high, report the corresponding test performance.
#    - Example code:
#      ```python
#      df = divided_result.sort_values(by='Train_best_model', ascending=False)
#      print(df.iloc[0].to_markdown())
#      ```
#  
# #### Binary Phenotypes Result Analysis
# 
# You can find the performance quality for binary phenotype using the following template:
# 
# ![PerformanceBinary](PerformanceBinary.PNG)
#  
# 
# This figure shows the 8 different scenarios that can exist in the results, and the following table explains each scenario.
# 
#  
# 
# We classified performance based on the following table:
# 
# | Performance Level      | Range    |
# |------------------------|----------|
# | **Low Performance**    | 0 to 0.5 |
# | **Moderate Performance** | 0.6 to 0.7 |
# | **High Performance**   | 0.8 to 1 |
# 
#  
# You can match the performance based on the following scenarios:
# 
# | Scenario                     | What's Happening                                                                                                   | Implication                                                                                             |
# |------------------------------|--------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
# | **High Test, High Train**     | The model performs well on both training and test datasets, effectively learning the underlying patterns.          | The model is well-tuned, generalizes well, and makes accurate predictions on both datasets.              |
# | **High Test, Moderate Train** | The model generalizes well but may not be fully optimized on training data, missing some underlying patterns.      | The model is fairly robust but may benefit from further tuning or more training to improve its learning. |
# | **High Test, Low Train**      | An unusual scenario, potentially indicating data leakage or overestimation of test performance.                    | The model’s performance is likely unreliable; investigate potential data issues or random noise.        |
# | **Moderate Test, High Train** | The model fits the training data well but doesn’t generalize as effectively, capturing only some test patterns.     | The model is slightly overfitting; adjustments may be needed to improve generalization on unseen data.   |
# | **Moderate Test, Moderate Train** | The model shows balanced but moderate performance on both datasets, capturing some patterns but missing others. | The model is moderately fitting; further improvements could be made in both training and generalization. |
# | **Moderate Test, Low Train**  | The model underperforms on training data and doesn’t generalize well, leading to moderate test performance.         | The model may need more complexity, additional features, or better training to improve on both datasets. |
# | **Low Test, High Train**      | The model overfits the training data, performing poorly on the test set.                                           | The model doesn’t generalize well; simplifying the model or using regularization may help reduce overfitting. |
# | **Low Test, Low Train**       | The model performs poorly on both training and test datasets, failing to learn the data patterns effectively.      | The model is underfitting; it may need more complexity, additional features, or more data to improve performance. |
# 
# ##### Recommendations for Publishing Results
# 
# When publishing results, scenarios with moderate train and moderate test performance can be used for complex phenotypes or diseases. However, results showing high train and moderate test, high train and high test, and moderate train and high test are recommended.
# 
# For most phenotypes, results typically fall in the moderate train and moderate test performance category.
# 
#  
# #### Continuous Phenotypes Result Analysis
# 
# You can find the performance quality for continuous phenotypes using the following template:
# 
# ![PerformanceContinous](PerformanceContinous.PNG)
# 
# This figure shows the 8 different scenarios that can exist in the results, and the following table explains each scenario.
# 
#  
# 
# We classified performance based on the following table:
# 
# | Performance Level      | Range        |
# |------------------------|--------------|
# | **Low Performance**    | 0 to 0.2     |
# | **Moderate Performance** | 0.3 to 0.7 |
# | **High Performance**   | 0.8 to 1     |
# 
#  
# 
# You can match the performance based on the following scenarios:
# 
# | Scenario                     | What's Happening                                                                                                   | Implication                                                                                             |
# |------------------------------|--------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
# | **High Test, High Train**     | The model performs well on both training and test datasets, effectively learning the underlying patterns.          | The model is well-tuned, generalizes well, and makes accurate predictions on both datasets.              |
# | **High Test, Moderate Train** | The model generalizes well but may not be fully optimized on training data, missing some underlying patterns.      | The model is fairly robust but may benefit from further tuning or more training to improve its learning. |
# | **High Test, Low Train**      | An unusual scenario, potentially indicating data leakage or overestimation of test performance.                    | The model’s performance is likely unreliable; investigate potential data issues or random noise.        |
# | **Moderate Test, High Train** | The model fits the training data well but doesn’t generalize as effectively, capturing only some test patterns.     | The model is slightly overfitting; adjustments may be needed to improve generalization on unseen data.   |
# | **Moderate Test, Moderate Train** | The model shows balanced but moderate performance on both datasets, capturing some patterns but missing others. | The model is moderately fitting; further improvements could be made in both training and generalization. |
# | **Moderate Test, Low Train**  | The model underperforms on training data and doesn’t generalize well, leading to moderate test performance.         | The model may need more complexity, additional features, or better training to improve on both datasets. |
# | **Low Test, High Train**      | The model overfits the training data, performing poorly on the test set.                                           | The model doesn’t generalize well; simplifying the model or using regularization may help reduce overfitting. |
# | **Low Test, Low Train**       | The model performs poorly on both training and test datasets, failing to learn the data patterns effectively.      | The model is underfitting; it may need more complexity, additional features, or more data to improve performance. |
# 
# ##### Recommendations for Publishing Results
# 
# When publishing results, scenarios with moderate train and moderate test performance can be used for complex phenotypes or diseases. However, results showing high train and moderate test, high train and high test, and moderate train and high test are recommended.
# 
# For most continuous phenotypes, results typically fall in the moderate train and moderate test performance category.
# 
# 
# 
# 
# 
# 
#  
# 
# ### 2. **Reporting Generalized Performance:**
#    - One can also report the generalized performance by calculating the difference between the training and test performance, and the sum of the test and training performance. Report the result or hyperparameter combination for which the sum is high and the difference is minimal.
#    - Example code:
#      ```python
#      df = divided_result.copy()
#      df['Difference'] = abs(df['Train_best_model'] - df['Test_best_model'])
#      df['Sum'] = df['Train_best_model'] + df['Test_best_model']
# 
#      sorted_df = df.sort_values(by=['Sum', 'Difference'], ascending=[False, True])
#      print(sorted_df.iloc[0].to_markdown())
#      ```
# 
# 
# ### 3. **Reporting Hyperparameters Affecting Test and Train Performance:**
#    - Find the hyperparameters that have more than one unique value and calculate their correlation with the following columns to understand how they are affecting the performance of train and test sets:
#      - `Train_null_model`
#      - `Train_pure_prs`
#      - `Train_best_model`
#      - `Test_pure_prs`
#      - `Test_null_model`
#      - `Test_best_model`
# 
# 
# 
# ### 4. Other Analysis
# 1. Once you have the results, you can find how hyperparameters affect the model performance.
# 2. Analysis, like overfitting and underfitting, can be performed as well.
# 3. The way you are going to report the results can vary.
# 4. Results can be visualized, and other patterns in the data can be explored.
# 

# In[7]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = divided_result.sort_values(by='Train_best_model', ascending=False)
print("1. Reporting Based on Best Training Performance:\n")
print(df.iloc[0].to_markdown())


 
df = divided_result.copy()

# Plot Train and Test best models against p-values
plt.figure(figsize=(10, 6))
plt.plot(df['pvalue'], df['Train_best_model'], label='Train_best_model', marker='o', color='royalblue')
plt.plot(df['pvalue'], df['Test_best_model'], label='Test_best_model', marker='o', color='darkorange')

# Highlight the p-value where both train and test are high
best_index = df[['Train_best_model']].sum(axis=1).idxmax()
best_pvalue = df.loc[best_index, 'pvalue']
best_train = df.loc[best_index, 'Train_best_model']
best_test = df.loc[best_index, 'Test_best_model']

# Use dark colors for the circles
plt.scatter(best_pvalue, best_train, color='darkred', s=100, label=f'Best Performance (Train)', edgecolor='black', zorder=5)
plt.scatter(best_pvalue, best_test, color='darkblue', s=100, label=f'Best Performance (Test)', edgecolor='black', zorder=5)

# Annotate the best performance with p-value, train, and test values
plt.text(best_pvalue, best_train, f'p={best_pvalue:.4g}\nTrain={best_train:.4g}', ha='right', va='bottom', fontsize=9, color='darkred')
plt.text(best_pvalue, best_test, f'p={best_pvalue:.4g}\nTest={best_test:.4g}', ha='right', va='top', fontsize=9, color='darkblue')

# Calculate Difference and Sum
df['Difference'] = abs(df['Train_best_model'] - df['Test_best_model'])
df['Sum'] = df['Train_best_model'] + df['Test_best_model']

# Sort the DataFrame
sorted_df = df.sort_values(by=['Sum', 'Difference'], ascending=[False, True])
#sorted_df = df.sort_values(by=[ 'Difference','Sum'], ascending=[  True,False])

# Highlight the general performance
general_index = sorted_df.index[0]
general_pvalue = sorted_df.loc[general_index, 'pvalue']
general_train = sorted_df.loc[general_index, 'Train_best_model']
general_test = sorted_df.loc[general_index, 'Test_best_model']

plt.scatter(general_pvalue, general_train, color='darkgreen', s=150, label='General Performance (Train)', edgecolor='black', zorder=6)
plt.scatter(general_pvalue, general_test, color='darkorange', s=150, label='General Performance (Test)', edgecolor='black', zorder=6)

# Annotate the general performance with p-value, train, and test values
plt.text(general_pvalue, general_train, f'p={general_pvalue:.4g}\nTrain={general_train:.4g}', ha='left', va='bottom', fontsize=9, color='darkgreen')
plt.text(general_pvalue, general_test, f'p={general_pvalue:.4g}\nTest={general_test:.4g}', ha='left', va='top', fontsize=9, color='darkorange')

# Add labels and legend
plt.xlabel('p-value')
plt.ylabel('Model Performance')
plt.title('Train vs Test Best Models')
plt.legend()
plt.show()
 




print("2. Reporting Generalized Performance:\n")
df = divided_result.copy()
df['Difference'] = abs(df['Train_best_model'] - df['Test_best_model'])
df['Sum'] = df['Train_best_model'] + df['Test_best_model']
sorted_df = df.sort_values(by=['Sum', 'Difference'], ascending=[False, True])
print(sorted_df.iloc[0].to_markdown())


print("3. Reporting the correlation of hyperparameters and the performance of 'Train_null_model', 'Train_pure_prs', 'Train_best_model', 'Test_pure_prs', 'Test_null_model', and 'Test_best_model':\n")

print("3. For string hyperparameters, we used one-hot encoding to find the correlation between string hyperparameters and 'Train_null_model', 'Train_pure_prs', 'Train_best_model', 'Test_pure_prs', 'Test_null_model', and 'Test_best_model'.")

print("3. We performed this analysis for those hyperparameters that have more than one unique value.")

correlation_columns = [
 'Train_null_model', 'Train_pure_prs', 'Train_best_model',
 'Test_pure_prs', 'Test_null_model', 'Test_best_model'
]

hyperparams = [col for col in divided_result.columns if len(divided_result[col].unique()) > 1]
hyperparams = list(set(hyperparams+correlation_columns))
 
# Separate numeric and string columns
numeric_hyperparams = [col for col in hyperparams if pd.api.types.is_numeric_dtype(divided_result[col])]
string_hyperparams = [col for col in hyperparams if pd.api.types.is_string_dtype(divided_result[col])]


# Encode string columns using one-hot encoding
divided_result_encoded = pd.get_dummies(divided_result, columns=string_hyperparams)

# Combine numeric hyperparams with the new one-hot encoded columns
encoded_columns = [col for col in divided_result_encoded.columns if col.startswith(tuple(string_hyperparams))]
hyperparams = numeric_hyperparams + encoded_columns
 

# Calculate correlations
correlations = divided_result_encoded[hyperparams].corr()
 
# Display correlation of hyperparameters with train/test performance columns
hyperparam_correlations = correlations.loc[hyperparams, correlation_columns]
 
hyperparam_correlations = hyperparam_correlations.fillna(0)

# Plotting the correlation heatmap
plt.figure(figsize=(12, 8))
ax = sns.heatmap(hyperparam_correlations, annot=True, cmap='viridis', fmt='.2f', cbar=True)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='right')

# Rotate y-axis labels to horizontal
#ax.set_yticklabels(ax.get_yticklabels(), rotation=0, va='center')

plt.title('Correlation of Hyperparameters with Train/Test Performance')
plt.show() 

sns.set_theme(style="whitegrid")  # Choose your preferred style
pairplot = sns.pairplot(divided_result_encoded[hyperparams],hue = 'Test_best_model', palette='viridis')

# Adjust the figure size
pairplot.fig.set_size_inches(15, 15)  # You can adjust the size as needed

for ax in pairplot.axes.flatten():
    ax.set_xlabel(ax.get_xlabel(), rotation=90, ha='right')  # X-axis labels vertical
    #ax.set_ylabel(ax.get_ylabel(), rotation=0, va='bottom')  # Y-axis labels horizontal

# Show the plot
plt.show()




# In[ ]:




