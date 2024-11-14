LDpred-2-Auto
=============

In this notebook, we will use LDpred-2-Auto to calculate the Polygenic
Risk Score (PRS).

**Note:** LDpred-2 needs to be installed using R. **Note:**
LDpred-2-Auto calculates multiple Betas and then uses a formula to
derive a single Beta value for all SNPs. However, we used each set of
Betas estimated by the LDpred-2-Auto model.

We will use the same flow of the code, and when we have to invoke
LDpred-2 functions, we will call the script using Python. When LDpred-2
finishes executions, results will be stored in the specific directories
for each fold and will be retrieved in Python for further calculations.

Install LDpred-2
----------------

Install LDpred-2 using the following commands:

1. Activate the conda Environment.
2. Type R
3. Run the following commands:

.. code:: r

   install.packages("remotes")
   library(remotes)
   remotes::install_github("https://github.com/privefl/bigsnpr.git")

The content in this notebook has been taken from the `PRS
Tutorial <https://choishingwan.github.io/PRS-Tutorial/ldpred/>`__,
`LDpred-2
Documentation <https://privefl.github.io/bigsnpr/articles/LDpred2.html>`__,
and `Polygenic Scores
(PGS) <https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html>`__.

Thanks to Florian Privé for creating one of the best documentations
explaining each aspect of the calculation. We will be using LDpred-2 for
cross-validation design and hyper-parameter optimization, and this
tutorial complements the existing one.

LDpred-2 Hyperparameters
------------------------

Hyperparameters for LDpred-2 performed using PLINK
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pruning Parameters
^^^^^^^^^^^^^^^^^^

Informs Plink that we wish to perform pruning with a window size of 200
variants, sliding across the genome with a step size of 50 variants at a
time, and filter out any SNPs with LD ( r^2 ) higher than 0.25.

.. code:: python

   1. p_window_size = [200]
   2. p_slide_size = [50]
   3. p_LD_threshold = [0.25]

Clumping Parameters
^^^^^^^^^^^^^^^^^^^

The P-value threshold for an SNP to be included. 1 means to include all
SNPs for clumping. SNPs having ( r^2 ) higher than 0.1 with the index
SNPs will be removed. SNPs within 200k of the index SNP are considered
for clumping.

.. code:: python

   1. clump_p1 = [1]
   2. clump_r2 = [0.1]
   3. clump_kb = [200]

PCA
^^^

Pca also affects the results evident from the initial analysis; however,
including more PCA overfits the model.

Hyperparameters for LDpred-2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Heritability
^^^^^^^^^^^^

To calculate h2 using LDpred-2, we followed the `LDpred-2
tutorial <https://choishingwan.github.io/PRS-Tutorial/ldpred/>`__. The
code for this part is in R, as LDpred-2 is in R.

1. Using SNPs from the HapMap as preferred by the authors. The name in
   the code is ``LDpred-2_hapmap``.
2. Using all the SNPs. The name in the code is ``LDpred-2_full``.

This approach is computationally expensive, as the correlation between
all the SNPs in a specific chromosome is being calculated.

LDpred-2 Models
---------------

Each model has its own set of hyperparameters, so first, we will
generate those hyperparameters for those models and call them when
invoking LDpred-2.

Use the following R code to know more about the hyperparameters: -
``help(snp_ldpred2_inf)`` - ``help(snp_ldpred2_grid)`` -
``help(snp_ldpred2_auto)`` - ``help(snp_lassosum2)``

1. LDpred2-inf: Infinitesimal Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No Hyperparameters

2. LDpred2(-grid): Grid of Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each model has its own set of parameters: - ``burn_in = 100`` -
``num_iter = 10`` - ``p = c(1e-4, 1e-5)`` - ``h2 = c(0.1, 0.2)`` -
``sparse = c(FALSE, TRUE)``

3. LDpred2-auto: Automatic Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``burn_in = 100``
-  ``num_iter = 100``
-  ``sparse = c(FALSE, TRUE)``
-  ``alpha_bounds = c(1, 0.5)``
-  ``use_MLE = c(FALSE, TRUE)``
-  ``p = c(1e-4, 1e-5)``
-  ``shrink_corr = 0.7``
-  ``allow_jump_sign = c(FALSE, TRUE)``

4. lassosum2: Grid of Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``lambda = c(0.001, 0.01, 0.1, 1)``
-  ``delta = 30``

We will execute the code for each mkodel seperatly.

GWAS file processing for LDpred-2.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LDpred-2 accepts betas so we converted OR to betas.

.. code:: ipython3

    
    import os
    import pandas as pd
    import numpy as np
    
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
    
        
    #filedirec = sys.argv[1]
    
    filedirec = "SampleData1"
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
    
     


.. parsed-literal::

    |    |   CHR |     BP | SNP        | A1   | A2   |      N |         SE |        P |        BETA |     INFO |      MAF |
    |---:|------:|-------:|:-----------|:-----|:-----|-------:|-----------:|---------:|------------:|---------:|---------:|
    |  0 |     1 | 756604 | rs3131962  | A    | G    | 388028 | 0.00301666 | 0.483171 | -0.00211532 | 0.890558 | 0.36939  |
    |  1 |     1 | 768448 | rs12562034 | A    | G    | 388028 | 0.00329472 | 0.834808 |  0.00068708 | 0.895894 | 0.336846 |
    |  2 |     1 | 779322 | rs4040617  | G    | A    | 388028 | 0.00303344 | 0.42897  | -0.00239932 | 0.897508 | 0.377368 |
    |  3 |     1 | 801536 | rs79373928 | G    | T    | 388028 | 0.00841324 | 0.808999 |  0.00203363 | 0.908963 | 0.483212 |
    |  4 |     1 | 808631 | rs11240779 | G    | A    | 388028 | 0.00242821 | 0.590265 |  0.00130747 | 0.893213 | 0.45041  |
    Length of DataFrame! 499617
    

Define Hyperparameters
~~~~~~~~~~~~~~~~~~~~~~

Define hyperparameters to be optimized and set initial values.

Extract Valid SNPs from Clumped File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For Windows, download ``gwak``, and for Linux, the ``awk`` command is
sufficient. For Windows, ``GWAK`` is required. You can download it from
`here <https://sourceforge.net/projects/gnuwin32/>`__. Get it and place
it in the same directory.

Execution Path
~~~~~~~~~~~~~~

At this stage, we have the genotype training data
``newtrainfilename = "train_data.QC"`` and genotype test data
``newtestfilename = "test_data.QC"``.

We modified the following variables:

1. ``filedirec = "SampleData1"`` or ``filedirec = sys.argv[1]``
2. ``foldnumber = "0"`` or ``foldnumber = sys.argv[2]`` for HPC.

Only these two variables can be modified to execute the code for
specific data and specific folds. Though the code can be executed
separately for each fold on HPC and separately for each dataset, it is
recommended to execute it for multiple diseases and one fold at a time.
Here’s the corrected text in Markdown format:

P-values
~~~~~~~~

PRS calculation relies on P-values. SNPs with low P-values, indicating a
high degree of association with a specific trait, are considered for
calculation.

You can modify the code below to consider a specific set of P-values and
save the file in the same format.

We considered the following parameters:

-  **Minimum P-value**: ``1e-10``
-  **Maximum P-value**: ``1.0``
-  **Minimum exponent**: ``10`` (Minimum P-value in exponent)
-  **Number of intervals**: ``100`` (Number of intervals to be
   considered)

The code generates an array of logarithmically spaced P-values:

.. code:: python

   import numpy as np
   import os

   minimumpvalue = 10  # Minimum exponent for P-values
   numberofintervals = 100  # Number of intervals to be considered

   allpvalues = np.logspace(-minimumpvalue, 0, numberofintervals, endpoint=True)  # Generating an array of logarithmically spaced P-values

   print("Minimum P-value:", allpvalues[0])
   print("Maximum P-value:", allpvalues[-1])

   count = 1
   with open(os.path.join(folddirec, 'range_list'), 'w') as file:
       for value in allpvalues:
           file.write(f'pv_{value} 0 {value}\n')  # Writing range information to the 'range_list' file
           count += 1

   pvaluefile = os.path.join(folddirec, 'range_list')

In this code: - ``minimumpvalue`` defines the minimum exponent for
P-values. - ``numberofintervals`` specifies how many intervals to
consider. - ``allpvalues`` generates an array of P-values spaced
logarithmically. - The script writes these P-values to a file named
``range_list`` in the specified directory.

2. LDpred2-auto: Automatic Model
--------------------------------

The parameters for the grid model can be divided into two classes: the
first one should be specified here in the python code and the other one
should be specified in the corresponding R file.

Should be specified here.
~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``burn_in = 100``
-  ``num_iter = 10``
-  ``shrink_corr = 0.7``
-  ``use_MLE = c(FALSE, TRUE)``
-  ``sparse = c(FALSE, TRUE)``
-  ``allow_jump_sign = c(FALSE, TRUE)``

Should be specified in R file.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``alpha_bounds = c(1, 0.5)``
-  ``p = c(1e-4, 1e-5)``

For each combination of hyper parameters, LDpred-2 grid model will
generate a set of betas, and we will be using those betas for each SNP
to compute the polygenic risks using ``plink –score``.

The following file contains the R code required to execute the code in
this notebook.

.. code:: r

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
   #9. Argument nine is burn_in. Example: `100`
   #10 Argument ten is num_iter. Example: `10`
   #11 Argument 11 is shrink_corr = 0.7`
   #12 Argument 12 is use_MLE = c(FALSE, TRUE)`
   #13 Argument 13 is sparse = c(FALSE, TRUE)`
   #14 Argument 14 is allow_jump_sign = c(FALSE, TRUE)`


    







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
     sumstats$beta <- sumstats$BETA
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
         #print("FUCK")
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
     
     help(snp_ldsc)
     
     # Here we have to specify the p and h2_seq
     multi_auto <- snp_ldpred2_auto(
       corr,
       df_beta,
       h2_init = h2_est,
       vec_p_init = seq_log(1e-4, 1, length.out = 10),
       burn_in = as.numeric(args[9]),
       num_iter = as.numeric(args[10]),
       
       sparse = tolower(args[11]) == "true",
       allow_jump_sign = tolower(args[12]) == "true",
       #shrink_corr = as.numeric(args[13]),
       use_MLE = tolower(args[14]) == "true",
       # Here we have specify the alpha bounds.
       alpha_bounds = c(-1.5, 0.5),
       ind.corr = cols_along(corr)
     )
     print("Sex")
    
     beta_auto <- sapply(multi_auto, function(auto)  auto$beta_est)
     beta_auto
     print("Sex")
     p_seq <- signif(seq_log(1e-4, 1, length.out = 10), 2)
     
     grid.param <-
       expand.grid(p = p_seq
                )
     
     # Save the grid paramters.
     gridparamters <- grid.param[,c("p" )]
     result <-paste(".",args[2],paste(args[3],toString(".ldpred_auto_parameters"), sep = ""),sep="//")
     write.table(gridparamters, file = result, row.names = FALSE,sep=",", quote = FALSE)
     
     
      
     
     result <-paste(".",args[2],paste(args[3],toString(".ldpred_auto_betas"), sep = ""),sep="//")
     write.table(beta_auto, file = result, row.names = FALSE,sep=",", quote = FALSE)
     newgwas <- info_snp[,c("rsid.ss", "a0", "beta")]
     
     result <-paste(".",args[2],paste(args[3],toString(".ldpred_auto_gwas"), sep = ""),sep="//")
     write.table(newgwas, file = result, row.names = FALSE,sep=",", quote = FALSE)
     print("Sex")
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
         #print("FUCK")
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
     # Here we have to specify the p and h2_seq
     
     # Here we have to specify the p and h2_seq
     multi_auto <- snp_ldpred2_auto(
       corr,
       df_beta,
       h2_init = h2_est,
       vec_p_init = seq_log(1e-4, 1, length.out = 10),
       burn_in = as.numeric(args[9]),
       num_iter = as.numeric(args[10]),
       
       sparse = tolower(args[11]) == "true",
       allow_jump_sign = tolower(args[12]) == "true",
       #shrink_corr = as.numeric(args[13]),
       use_MLE = tolower(args[14]) == "true",
       # Here we have specify the alpha bounds.
       alpha_bounds = c(-1.5, 0.5),
       ind.corr = cols_along(corr)
     )
     print("Sex")
     
     beta_auto <- sapply(multi_auto, function(auto)  auto$beta_est)
     beta_auto
     print("Sex")
     p_seq <- signif(seq_log(1e-4, 1, length.out = 10), 2)
     
     grid.param <-
       expand.grid(p = p_seq
       )
     
     # Save the grid paramters.
     gridparamters <- grid.param[,c("p" )]
     result <-paste(".",args[2],paste(args[3],toString(".ldpred_auto_parameters"), sep = ""),sep="//")
     write.table(gridparamters, file = result, row.names = FALSE,sep=",", quote = FALSE)
     
     
     
     
     result <-paste(".",args[2],paste(args[3],toString(".ldpred_auto_betas"), sep = ""),sep="//")
     write.table(beta_auto, file = result, row.names = FALSE,sep=",", quote = FALSE)
     newgwas <- info_snp[,c("rsid.ss", "a0", "beta")]
     
     result <-paste(".",args[2],paste(args[3],toString(".ldpred_auto_gwas"), sep = ""),sep="//")
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
     
     #print(ldsc)
     #exit(0)
     
     
   }

.. code:: ipython3

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
    
     
    #foldnumber = sys.argv[1]
    foldnumber = "0"  # Setting 'foldnumber' to "0"
    
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
    #prs_result = pd.DataFrame(columns=["clump_p1", "clump_r2", "clump_kb", "p_window_size", "p_slide_size", "p_LD_threshold",
    #                                   "pvalue", "numberofpca","h2model","numberofvariants","Train_pure_prs", "Train_null_model", "Train_best_model",
    #                                   "Test_pure_prs", "Test_null_model", "Test_best_model"])

Define Helper Functions
~~~~~~~~~~~~~~~~~~~~~~~

1. **Perform Clumping and Pruning**
2. **Calculate PCA Using Plink**
3. **Fit Binary Phenotype and Save Results**
4. **Fit Continuous Phenotype and Save Results**

.. code:: ipython3

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
    def fit_binary_phenotype_on_PRS(traindirec, newtrainfilename,h2model,tempp, p,b,iterr,shrink_corr,use_MLE,sparse,allow_jump_sign, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val, Name, pvaluefile,heritability,numberofvariants):
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
                  
                     
    
    
                        "burn_in":b,
                        "num_iter":iterr,
                        "sparse":sparse,
                        "allow_jump_sign":allow_jump_sign,
                        "shrink_corr":shrink_corr,
                        "use_MLE":use_MLE,
     
                        "temp_pvalue":tempp,                 
                         
                        "numberofvariants": numberofvariants,
           
                        "heritability_model":h2model,
                        "h2":heritability,
    
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
    def fit_continous_phenotype_on_PRS(traindirec, newtrainfilename,h2model,tempp, p,b,iterr,shrink_corr,use_MLE,sparse,allow_jump_sign, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val, Name, pvaluefile,heritability,numberofvariants):
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
                      
                        "temp_pvalue":tempp,
                        "burn_in":b,
                        "num_iter":iterr,
                        "sparse":sparse,
                        "allow_jump_sign":allow_jump_sign,
                        "shrink_corr":shrink_corr,
                        "use_MLE":use_MLE,
     
                        "Train_pure_prs":explained_variance_score(phenotype_train["Phenotype"],prs_train['SCORE'].values),
                        "Train_null_model":explained_variance_score(phenotype_train["Phenotype"],train_null_predicted),
                        "Train_best_model":explained_variance_score(phenotype_train["Phenotype"],train_best_predicted),
                        
                        "Test_pure_prs":explained_variance_score(phenotype_test["Phenotype"],prs_test['SCORE'].values),
                        "Test_null_model":explained_variance_score(phenotype_test["Phenotype"],test_null_predicted),
                        "Test_best_model":explained_variance_score(phenotype_test["Phenotype"],test_best_predicted),
                        
                    }, ignore_index=True)
    
              
                    prs_result.to_csv(traindirec+os.sep+Name+os.sep+"Results.csv",index=False)
         
        return
    

Execute LDpred-2-Auto
---------------------

.. code:: ipython3

    def transform_ldpred2_auto_data(traindirec, newtrainfilename,model,p,b,iterr,shrink_corr,use_MLE,sparse,allow_jump_sign, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
        ### First perform clumping on the file and save the clumpled file.
        #perform_clumping_and_pruning_on_individual_data(traindirec, newtrainfilename,p, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile)
       
            
        # Also extract the PCA at this point for both test and training data.
        #calculate_pca_for_traindata_testdata_for_clumped_pruned_snps(traindirec, newtrainfilename,p)
    
        #Extract p-values from the GWAS file.
        # Command for Linux.
        os.system("awk "+"\'"+"{print $3,$8}"+"\'"+" ./"+filedirec+os.sep+filedirec+".txt >  ./"+traindirec+os.sep+"SNP.pvalue")
     
        
        # Define the paths to the files
        files_to_remove = [
            traindirec+os.sep+"ldpred_h2_full.txt",
            traindirec+os.sep+"ldpred_h2_variants.txt",
            traindirec+os.sep+"train_data.ldpred_auto_parameters",
            traindirec+os.sep+"train_data.ldpred_auto_betas",
            traindirec+os.sep+"train_data.ldpred_auto_gwas",
            traindirec+os.sep+"ldpred_h2_full.txt",
            traindirec+os.sep+"train_data.ldpred_auto_gwas_final",
            
    
        ]
    
        # Loop through the files and remove them if they exist
        for file_path in files_to_remove:
            if os.path.exists(file_path):
                os.remove(file_path)
                print(f"Removed: {file_path}")
            else:
                print(f"File does not exist: {file_path}")
    
     
        if model=="LDpred-2_full":
            command = [
                "Rscript",
                "LDpred2-auto.R",
                os.path.join(filedirec),
                traindirec,
                trainfilename,
                newtrainfilename + ".clumped.pruned",
                "3",
                c3_val,
                c1_val,
                c2_val,
                p,
                b,
                iterr,
                sparse,
                allow_jump_sign,
                shrink_corr,
                use_MLE
                
                
            ]
            try:
                subprocess.run(command)
                heritability = pd.read_csv(traindirec+os.sep+"ldpred_h2_full.txt",sep="\s+",header=None)[1].values[0]
                numberofvariants =  pd.read_csv(traindirec+os.sep+"ldpred_h2_variants.txt",sep="\s+",header=None)[1].values[0] 
            except:
                print("For model ",models,"it did not work! may be h2 is negative.")
                return            
    
        if model=="LDpred-2_hapmap":
            command = [
                "Rscript",
                "LDpred2-auto.R",
                os.path.join(filedirec),
                traindirec,
                trainfilename,
                newtrainfilename + ".clumped.pruned",
                "2",
                c3_val,
                c1_val,
                c2_val,
                p,
                b,
                iterr,
                sparse,
                allow_jump_sign,
                shrink_corr,
                use_MLE
            ]
            try:
                subprocess.run(command)
                heritability = pd.read_csv(traindirec+os.sep+"ldpred_h2_hapmap.txt",sep="\s+",header=None)[1].values[0]
                numberofvariants =  pd.read_csv(traindirec+os.sep+"ldpred_h2_variants.txt",sep="\s+",header=None)[1].values[0] 
            except:
                print("For model ",models,"it did not work! may be h2 is negative.")
                return            
    
            
        
        # After obtaining betas from the LDpred-2
        # Append the new betas.
        # Read all the betas from the file.
        # First read the paramters.
        gridparameters = pd.read_csv(traindirec+os.sep+"train_data.ldpred_auto_parameters",sep=",")
        allbetas = pd.read_csv(traindirec+os.sep+"train_data.ldpred_auto_betas",sep=",")
        allbetas = allbetas.fillna(0)
        
    
        for index, row in gridparameters.iterrows():
            gwas = pd.read_csv(traindirec+os.sep+"train_data.ldpred_auto_gwas",sep=",")
            # Accessing individual elements in the row
            tempp = row['x']
        
            
            gwas["newbeta"] = allbetas.iloc[:, index].values
            
            if check_phenotype_is_binary_or_continous(filedirec)=="Binary":
                gwas["newbeta"] = np.exp(gwas["newbeta"])
            else:
                pass
            
      
            gwas.rename(columns={'rsid.ss': 'SNP', 'a0': 'A1','newbeta':'BETA'}, inplace=True)
            gwas[["SNP","A1","BETA"]].to_csv(traindirec+os.sep+"train_data.ldpred_auto_gwas_final",sep="\t",index=False)
                 
             
     
    
            command = [
                "./plink",
                 "--bfile", traindirec+os.sep+newtrainfilename+".clumped.pruned",
                ### SNP column = 1, Effect allele column 2 = 4, Effect column=4
                "--score", traindirec+os.sep+"train_data.ldpred_auto_gwas_final", "1", "2", "3", "header",
                "--q-score-range", traindirec+os.sep+"range_list",traindirec+os.sep+"SNP.pvalue",
                #"--extract", traindirec+os.sep+trainfilename+".valid.snp",
                "--out", traindirec+os.sep+Name+os.sep+trainfilename
            ]
    
            subprocess.run(command)
    
    
            command = [
                "./plink",
                "--bfile", folddirec+os.sep+testfilename,
                ### SNP column = 3, Effect allele column 1 = 4, Beta column=12
                "--score", traindirec+os.sep+"train_data.ldpred_auto_gwas_final", "1", "2", "3", "header",
                "--q-score-range", traindirec+os.sep+"range_list",traindirec+os.sep+"SNP.pvalue",
                "--out", folddirec+os.sep+Name+os.sep+testfilename
            ]
            subprocess.run(command)
    
            if check_phenotype_is_binary_or_continous(filedirec)=="Binary":
                print("Binary Phenotype!")
                fit_binary_phenotype_on_PRS(folddirec, newtrainfilename,model,tempp, p,b,iterr,shrink_corr,use_MLE,sparse,allow_jump_sign, str(p1_val), str(p2_val), str(p3_val), str(c1_val), str(c2_val), str(c3_val), Name, pvaluefile,heritability,numberofvariants)
            else:
                print("Continous Phenotype!")
                fit_continous_phenotype_on_PRS(folddirec, newtrainfilename,model,tempp, p,b,iterr,shrink_corr,use_MLE,sparse,allow_jump_sign, str(p1_val), str(p2_val), str(p3_val), str(c1_val), str(c2_val), str(c3_val), Name, pvaluefile,heritability,numberofvariants)
    
    
     
     
     
    #h2models = ["LDpred-2_full","LDpred-2_hapmap"]
    h2models = ["LDpred-2_hapmap","LDpred-2_full"]
    
    burn_in = ['200']
    num_iter = ['10']
    shrink_corrs = ['0.1','0.2','0.3','0.4','0.5','0.6']
    shrink_corrs = ['0.1']
    
    use_MLEs = ["FALSE"]
    sparses = ["FALSE"]
    allow_jump_signs = ["TRUE"]
    
    
    result_directory = "ldpred2_auto"
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
              for shrink_corr in shrink_corrs:
               for use_MLE in use_MLEs:
                for sparse in sparses:
                 for allow_jump_sign in allow_jump_signs:
                  transform_ldpred2_auto_data(folddirec, newtrainfilename,model, p,b,i,shrink_corr,use_MLE,sparse,allow_jump_sign, str(p1_val), str(p2_val), str(p3_val), str(c1_val), str(c2_val), str(c3_val), result_directory, pvaluefile)
    
    


.. parsed-literal::

    File does not exist: SampleData1/Fold_0/ldpred_h2_full.txt
    File does not exist: SampleData1/Fold_0/ldpred_h2_variants.txt
    File does not exist: SampleData1/Fold_0/train_data.ldpred_auto_parameters
    File does not exist: SampleData1/Fold_0/train_data.ldpred_auto_betas
    File does not exist: SampleData1/Fold_0/train_data.ldpred_auto_gwas
    File does not exist: SampleData1/Fold_0/ldpred_h2_full.txt
    File does not exist: SampleData1/Fold_0/train_data.ldpred_auto_gwas_final
     [1] "SampleData1"                  "SampleData1/Fold_0"          
     [3] "train_data"                   "train_data.QC.clumped.pruned"
     [5] "2"                            "200"                         
     [7] "1"                            "0.1"                         
     [9] "6"                            "200"                         
    [11] "10"                           "FALSE"                       
    [13] "TRUE"                         "0.1"                         
    [15] "FALSE"                       
    

.. parsed-literal::

    Loading required package: bigstatsr
    trying URL 'https://ndownloader.figshare.com/files/25503788'
    Content type 'application/octet-stream' length 35191166 bytes (33.6 MB)
    ==================================================
    downloaded 33.6 MB
    
    

.. parsed-literal::

    [1] "Directory .//SampleData1/Fold_0//tmp-data deleted."
    [1] "File .//SampleData1/Fold_0//train_data.QC.clumped.pruned.rds deleted."
    [1] "File .//SampleData1/Fold_0//train_data.QC.clumped.pruned.bk deleted."
    

.. parsed-literal::

    132,375 variants to be matched.
    0 ambiguous SNPs have been removed.
    35,527 variants have been matched; 0 were flipped and 847 were reversed.
    

.. parsed-literal::

    [1] 2975
    

.. parsed-literal::

    Creating directory ".//SampleData1/Fold_0//tmp-data" which didn't exist..
    

.. parsed-literal::

    [1] 2599
    [1] 2230
    [1] 1924
    [1] 2121
    [1] 2192
    [1] 1794
    [1] 1715
    [1] 1702
    [1] 1917
    [1] 1675
    [1] 1802
    [1] 1352
    [1] 1230
    [1] 1178
    [1] 1222
    [1] 1253
    [1] 1164
    [1] 1066
    [1] 1100
    [1] 601
    [1] 715
    [1] "File .//SampleData1/Fold_0//ldpred_h2_hapmap.txt deleted."
    

.. parsed-literal::

    Loading required package: foreach
    Loading required package: rngtools
    

.. parsed-literal::

    [1] "Sex"
    [1] "Sex"
    [1] "Sex"
    

.. parsed-literal::

    Warning message:
    package ‘magrittr’ was built under R version 4.1.3 
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 35527 valid predictors loaded.
    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    

.. parsed-literal::

    Warning: 464091 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    File does not exist: SampleData1/Fold_0/ldpred_h2_full.txt
    Removed: SampleData1/Fold_0/ldpred_h2_variants.txt
    Removed: SampleData1/Fold_0/train_data.ldpred_auto_parameters
    Removed: SampleData1/Fold_0/train_data.ldpred_auto_betas
    Removed: SampleData1/Fold_0/train_data.ldpred_auto_gwas
    File does not exist: SampleData1/Fold_0/ldpred_h2_full.txt
    Removed: SampleData1/Fold_0/train_data.ldpred_auto_gwas_final
     [1] "SampleData1"                  "SampleData1/Fold_0"          
     [3] "train_data"                   "train_data.QC.clumped.pruned"
     [5] "3"                            "200"                         
     [7] "1"                            "0.1"                         
     [9] "6"                            "200"                         
    [11] "10"                           "FALSE"                       
    [13] "TRUE"                         "0.1"                         
    [15] "FALSE"                       
    

.. parsed-literal::

    Loading required package: bigstatsr
    trying URL 'https://ndownloader.figshare.com/files/25503788'
    Content type 'application/octet-stream' length 35191166 bytes (33.6 MB)
    ==================================================
    downloaded 33.6 MB
    
    

.. parsed-literal::

    [1] "Directory .//SampleData1/Fold_0//tmp-data deleted."
    [1] "File .//SampleData1/Fold_0//train_data.QC.clumped.pruned.rds deleted."
    [1] "File .//SampleData1/Fold_0//train_data.QC.clumped.pruned.bk deleted."
    

.. parsed-literal::

    499,617 variants to be matched.
    0 ambiguous SNPs have been removed.
    172,878 variants have been matched; 0 were flipped and 1,662 were reversed.
    

.. parsed-literal::

    [1] 14013
    [1] 14013
    

.. parsed-literal::

    Creating directory ".//SampleData1/Fold_0//tmp-data" which didn't exist..
    

.. parsed-literal::

    [1] 13813
    [1] 13813
    [1] 11785
    [1] 11785
    [1] 11041
    [1] 11041
    [1] 10632
    [1] 10632
    [1] 10068
    [1] 10068
    [1] 9496
    [1] 9496
    [1] 8867
    [1] 8867
    [1] 7768
    [1] 7768
    [1] 8824
    [1] 8824
    [1] 8420
    [1] 8420
    [1] 8198
    [1] 8198
    [1] 6350
    [1] 6350
    [1] 5742
    [1] 5742
    [1] 5569
    [1] 5569
    [1] 6069
    [1] 6069
    [1] 5723
    [1] 5723
    [1] 5578
    [1] 5578
    [1] 4364
    [1] 4364
    [1] 4916
    [1] 4916
    [1] 2811
    [1] 2811
    [1] 2831
    [1] 2831
    

.. parsed-literal::

    Loading required package: foreach
    Loading required package: rngtools
    

.. parsed-literal::

    [1] "Sex"
    [1] "Sex"
    

.. parsed-literal::

    Warning message:
    package ‘magrittr’ was built under R version 4.1.3 
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    Continous Phenotype!
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --out SampleData1/Fold_0/ldpred2_auto/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999891.
    172878 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/ldpred2_auto/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --out SampleData1/Fold_0/ldpred2_auto/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/train_data.ldpred_auto_gwas_final 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/ldpred2_auto/test_data.*.profile.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    Continous Phenotype!
    

Repeat the process for each fold.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Change the ``foldnumber`` variable.

.. code:: python

   #foldnumber = sys.argv[1]
   foldnumber = "0"  # Setting 'foldnumber' to "0"

Or uncomment the following line:

.. code:: python

   # foldnumber = sys.argv[1]
   python LDpred-2-Auto.py 0
   python LDpred-2-Auto.py 1
   python LDpred-2-Auto.py 2
   python LDpred-2-Auto.py 3
   python LDpred-2-Auto.py 4

The following files should exist after the execution:

1. ``SampleData1/Fold_0/ldpred2_auto/Results.csv``
2. ``SampleData1/Fold_1/ldpred2_auto/Results.csv``
3. ``SampleData1/Fold_2/ldpred2_auto/Results.csv``
4. ``SampleData1/Fold_3/ldpred2_auto/Results.csv``
5. ``SampleData1/Fold_4/ldpred2_auto/Results.csv``

Check the results file for each fold.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    import os
     
    import pandas as pd
    import numpy as np
    
    
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
    
    


.. parsed-literal::

    Fold_ 0 Yes, the file exists.
    Number of P-values processed:  200
    Fold_ 1 Yes, the file exists.
    Number of P-values processed:  200
    Fold_ 2 Yes, the file exists.
    Number of P-values processed:  200
    Fold_ 3 Yes, the file exists.
    Number of P-values processed:  200
    Fold_ 4 Yes, the file exists.
    Number of P-values processed:  200
    

Sum the results for each fold.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

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
            "heritability_model",
            "unique_h2",
            "grid_pvalue",
            "burn_in", 
            "num_iter",
            "sparse",
            "temp_pvalue",              
            "allow_jump_sign" ,
            "shrink_corr" ,
            "use_MLE" ,
            #"sparsity",
            "lasso_parameters_count",]
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
    
     


.. parsed-literal::

    We have to ensure when we sum the entries across all Folds, the same rows are merged!
    Fold_ 0 Yes, the file exists.
    Fold_ 1 Yes, the file exists.
    Fold_ 2 Yes, the file exists.
    Fold_ 3 Yes, the file exists.
    Fold_ 4 Yes, the file exists.
    Iteration 1:
    Unique rows in current common DataFrame: 200
    Unique rows in next DataFrame: 200
    Common rows after merge: 200
    
    Iteration 2:
    Unique rows in current common DataFrame: 200
    Unique rows in next DataFrame: 200
    Common rows after merge: 200
    
    Iteration 3:
    Unique rows in current common DataFrame: 200
    Unique rows in next DataFrame: 200
    Common rows after merge: 200
    
    Iteration 4:
    Unique rows in current common DataFrame: 200
    Unique rows in next DataFrame: 200
    Common rows after merge: 200
    
    DataFrame 1 with extracted common rows has 200 rows.
    DataFrame 2 with extracted common rows has 200 rows.
    DataFrame 3 with extracted common rows has 200 rows.
    DataFrame 4 with extracted common rows has 200 rows.
    DataFrame 5 with extracted common rows has 200 rows.
         clump_p1  clump_r2  clump_kb  p_window_size  p_slide_size  \
    0         1.0       0.1     200.0          200.0          50.0   
    1         1.0       0.1     200.0          200.0          50.0   
    2         1.0       0.1     200.0          200.0          50.0   
    3         1.0       0.1     200.0          200.0          50.0   
    4         1.0       0.1     200.0          200.0          50.0   
    ..        ...       ...       ...            ...           ...   
    195       1.0       0.1     200.0          200.0          50.0   
    196       1.0       0.1     200.0          200.0          50.0   
    197       1.0       0.1     200.0          200.0          50.0   
    198       1.0       0.1     200.0          200.0          50.0   
    199       1.0       0.1     200.0          200.0          50.0   
    
         p_LD_threshold        pvalue  numberofpca  tempalpha  l1weight  ...  \
    0              0.25  1.000000e-10          6.0        0.1       0.1  ...   
    1              0.25  3.359818e-10          6.0        0.1       0.1  ...   
    2              0.25  1.128838e-09          6.0        0.1       0.1  ...   
    3              0.25  3.792690e-09          6.0        0.1       0.1  ...   
    4              0.25  1.274275e-08          6.0        0.1       0.1  ...   
    ..              ...           ...          ...        ...       ...  ...   
    195            0.25  7.847600e-03          6.0        0.1       0.1  ...   
    196            0.25  2.636651e-02          6.0        0.1       0.1  ...   
    197            0.25  8.858668e-02          6.0        0.1       0.1  ...   
    198            0.25  2.976351e-01          6.0        0.1       0.1  ...   
    199            0.25  1.000000e+00          6.0        0.1       0.1  ...   
    
         Train_pure_prs  Train_null_model  Train_best_model  Test_pure_prs  \
    0         -0.000043           0.23001          0.236510      -0.000049   
    1         -0.000042           0.23001          0.237036      -0.000046   
    2         -0.000046           0.23001          0.240716      -0.000052   
    3         -0.000046           0.23001          0.242414      -0.000046   
    4         -0.000044           0.23001          0.245855      -0.000043   
    ..              ...               ...               ...            ...   
    195       -0.000012           0.23001          0.327007      -0.000013   
    196       -0.000009           0.23001          0.335955      -0.000010   
    197       -0.000006           0.23001          0.339013      -0.000006   
    198       -0.000004           0.23001          0.341792      -0.000004   
    199       -0.000002           0.23001          0.341717      -0.000002   
    
         Test_null_model  Test_best_model  heritability_model  sparse  \
    0           0.118692         0.124425       LDpred-2_full   False   
    1           0.118692         0.126365       LDpred-2_full   False   
    2           0.118692         0.132428       LDpred-2_full   False   
    3           0.118692         0.129256       LDpred-2_full   False   
    4           0.118692         0.135784       LDpred-2_full   False   
    ..               ...              ...                 ...     ...   
    195         0.118692         0.253170       LDpred-2_full   False   
    196         0.118692         0.280130       LDpred-2_full   False   
    197         0.118692         0.287047       LDpred-2_full   False   
    198         0.118692         0.297977       LDpred-2_full   False   
    199         0.118692         0.293601       LDpred-2_full   False   
    
         allow_jump_sign  use_MLE  
    0               True    False  
    1               True    False  
    2               True    False  
    3               True    False  
    4               True    False  
    ..               ...      ...  
    195             True    False  
    196             True    False  
    197             True    False  
    198             True    False  
    199             True    False  
    
    [200 rows x 26 columns]
    

.. parsed-literal::

    /tmp/ipykernel_781499/3471353363.py:24: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      non_numerical_df[non_numerical_cols] = non_numerical_df[non_numerical_cols].combine_first(df[non_numerical_cols])
    /tmp/ipykernel_781499/3471353363.py:24: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      non_numerical_df[non_numerical_cols] = non_numerical_df[non_numerical_cols].combine_first(df[non_numerical_cols])
    /tmp/ipykernel_781499/3471353363.py:24: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      non_numerical_df[non_numerical_cols] = non_numerical_df[non_numerical_cols].combine_first(df[non_numerical_cols])
    /tmp/ipykernel_781499/3471353363.py:24: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      non_numerical_df[non_numerical_cols] = non_numerical_df[non_numerical_cols].combine_first(df[non_numerical_cols])
    

Results
-------

1. **Reporting Based on Best Training Performance:**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  One can report the results based on the best performance of the
   training data. For example, if for a specific combination of
   hyperparameters, the training performance is high, report the
   corresponding test performance.
-  Example code:
   ``python      df = divided_result.sort_values(by='Train_best_model', ascending=False)      print(df.iloc[0].to_markdown())``

Binary Phenotypes Result Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can find the performance quality for binary phenotype using the
following template:

.. figure:: PerformanceBinary.PNG
   :alt: PerformanceBinary

   PerformanceBinary

This figure shows the 8 different scenarios that can exist in the
results, and the following table explains each scenario.

We classified performance based on the following table:

======================== ==========
Performance Level        Range
======================== ==========
**Low Performance**      0 to 0.5
**Moderate Performance** 0.6 to 0.7
**High Performance**     0.8 to 1
======================== ==========

You can match the performance based on the following scenarios:

+-------+--------------------------------+-----------------------------+
| Sce   | What’s Happening               | Implication                 |
| nario |                                |                             |
+=======+================================+=============================+
| *     | The model performs well on     | The model is well-tuned,    |
| *High | both training and test         | generalizes well, and makes |
| Test, | datasets, effectively learning | accurate predictions on     |
| High  | the underlying patterns.       | both datasets.              |
| Tr    |                                |                             |
| ain** |                                |                             |
+-------+--------------------------------+-----------------------------+
| *     | The model generalizes well but | The model is fairly robust  |
| *High | may not be fully optimized on  | but may benefit from        |
| Test, | training data, missing some    | further tuning or more      |
| Mod   | underlying patterns.           | training to improve its     |
| erate |                                | learning.                   |
| Tr    |                                |                             |
| ain** |                                |                             |
+-------+--------------------------------+-----------------------------+
| *     | An unusual scenario,           | The model’s performance is  |
| *High | potentially indicating data    | likely unreliable;          |
| Test, | leakage or overestimation of   | investigate potential data  |
| Low   | test performance.              | issues or random noise.     |
| Tr    |                                |                             |
| ain** |                                |                             |
+-------+--------------------------------+-----------------------------+
| **Mod | The model fits the training    | The model is slightly       |
| erate | data well but doesn’t          | overfitting; adjustments    |
| Test, | generalize as effectively,     | may be needed to improve    |
| High  | capturing only some test       | generalization on unseen    |
| Tr    | patterns.                      | data.                       |
| ain** |                                |                             |
+-------+--------------------------------+-----------------------------+
| **Mod | The model shows balanced but   | The model is moderately     |
| erate | moderate performance on both   | fitting; further            |
| Test, | datasets, capturing some       | improvements could be made  |
| Mod   | patterns but missing others.   | in both training and        |
| erate |                                | generalization.             |
| Tr    |                                |                             |
| ain** |                                |                             |
+-------+--------------------------------+-----------------------------+
| **Mod | The model underperforms on     | The model may need more     |
| erate | training data and doesn’t      | complexity, additional      |
| Test, | generalize well, leading to    | features, or better         |
| Low   | moderate test performance.     | training to improve on both |
| Tr    |                                | datasets.                   |
| ain** |                                |                             |
+-------+--------------------------------+-----------------------------+
| **Low | The model overfits the         | The model doesn’t           |
| Test, | training data, performing      | generalize well;            |
| High  | poorly on the test set.        | simplifying the model or    |
| Tr    |                                | using regularization may    |
| ain** |                                | help reduce overfitting.    |
+-------+--------------------------------+-----------------------------+
| **Low | The model performs poorly on   | The model is underfitting;  |
| Test, | both training and test         | it may need more            |
| Low   | datasets, failing to learn the | complexity, additional      |
| Tr    | data patterns effectively.     | features, or more data to   |
| ain** |                                | improve performance.        |
+-------+--------------------------------+-----------------------------+

Recommendations for Publishing Results
''''''''''''''''''''''''''''''''''''''

When publishing results, scenarios with moderate train and moderate test
performance can be used for complex phenotypes or diseases. However,
results showing high train and moderate test, high train and high test,
and moderate train and high test are recommended.

For most phenotypes, results typically fall in the moderate train and
moderate test performance category.

Continuous Phenotypes Result Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can find the performance quality for continuous phenotypes using the
following template:

.. figure:: PerformanceContinous.PNG
   :alt: PerformanceContinous

   PerformanceContinous

This figure shows the 8 different scenarios that can exist in the
results, and the following table explains each scenario.

We classified performance based on the following table:

======================== ==========
Performance Level        Range
======================== ==========
**Low Performance**      0 to 0.2
**Moderate Performance** 0.3 to 0.7
**High Performance**     0.8 to 1
======================== ==========

You can match the performance based on the following scenarios:

+-------+--------------------------------+-----------------------------+
| Sce   | What’s Happening               | Implication                 |
| nario |                                |                             |
+=======+================================+=============================+
| *     | The model performs well on     | The model is well-tuned,    |
| *High | both training and test         | generalizes well, and makes |
| Test, | datasets, effectively learning | accurate predictions on     |
| High  | the underlying patterns.       | both datasets.              |
| Tr    |                                |                             |
| ain** |                                |                             |
+-------+--------------------------------+-----------------------------+
| *     | The model generalizes well but | The model is fairly robust  |
| *High | may not be fully optimized on  | but may benefit from        |
| Test, | training data, missing some    | further tuning or more      |
| Mod   | underlying patterns.           | training to improve its     |
| erate |                                | learning.                   |
| Tr    |                                |                             |
| ain** |                                |                             |
+-------+--------------------------------+-----------------------------+
| *     | An unusual scenario,           | The model’s performance is  |
| *High | potentially indicating data    | likely unreliable;          |
| Test, | leakage or overestimation of   | investigate potential data  |
| Low   | test performance.              | issues or random noise.     |
| Tr    |                                |                             |
| ain** |                                |                             |
+-------+--------------------------------+-----------------------------+
| **Mod | The model fits the training    | The model is slightly       |
| erate | data well but doesn’t          | overfitting; adjustments    |
| Test, | generalize as effectively,     | may be needed to improve    |
| High  | capturing only some test       | generalization on unseen    |
| Tr    | patterns.                      | data.                       |
| ain** |                                |                             |
+-------+--------------------------------+-----------------------------+
| **Mod | The model shows balanced but   | The model is moderately     |
| erate | moderate performance on both   | fitting; further            |
| Test, | datasets, capturing some       | improvements could be made  |
| Mod   | patterns but missing others.   | in both training and        |
| erate |                                | generalization.             |
| Tr    |                                |                             |
| ain** |                                |                             |
+-------+--------------------------------+-----------------------------+
| **Mod | The model underperforms on     | The model may need more     |
| erate | training data and doesn’t      | complexity, additional      |
| Test, | generalize well, leading to    | features, or better         |
| Low   | moderate test performance.     | training to improve on both |
| Tr    |                                | datasets.                   |
| ain** |                                |                             |
+-------+--------------------------------+-----------------------------+
| **Low | The model overfits the         | The model doesn’t           |
| Test, | training data, performing      | generalize well;            |
| High  | poorly on the test set.        | simplifying the model or    |
| Tr    |                                | using regularization may    |
| ain** |                                | help reduce overfitting.    |
+-------+--------------------------------+-----------------------------+
| **Low | The model performs poorly on   | The model is underfitting;  |
| Test, | both training and test         | it may need more            |
| Low   | datasets, failing to learn the | complexity, additional      |
| Tr    | data patterns effectively.     | features, or more data to   |
| ain** |                                | improve performance.        |
+-------+--------------------------------+-----------------------------+

.. _recommendations-for-publishing-results-1:

Recommendations for Publishing Results
''''''''''''''''''''''''''''''''''''''

When publishing results, scenarios with moderate train and moderate test
performance can be used for complex phenotypes or diseases. However,
results showing high train and moderate test, high train and high test,
and moderate train and high test are recommended.

For most continuous phenotypes, results typically fall in the moderate
train and moderate test performance category.

2. **Reporting Generalized Performance:**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  One can also report the generalized performance by calculating the
   difference between the training and test performance, and the sum of
   the test and training performance. Report the result or
   hyperparameter combination for which the sum is high and the
   difference is minimal.

-  Example code:

   .. code:: python

      df = divided_result.copy()
      df['Difference'] = abs(df['Train_best_model'] - df['Test_best_model'])
      df['Sum'] = df['Train_best_model'] + df['Test_best_model']

      sorted_df = df.sort_values(by=['Sum', 'Difference'], ascending=[False, True])
      print(sorted_df.iloc[0].to_markdown())

3. **Reporting Hyperparameters Affecting Test and Train Performance:**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Find the hyperparameters that have more than one unique value and
   calculate their correlation with the following columns to understand
   how they are affecting the performance of train and test sets:

   -  ``Train_null_model``
   -  ``Train_pure_prs``
   -  ``Train_best_model``
   -  ``Test_pure_prs``
   -  ``Test_null_model``
   -  ``Test_best_model``

4. Other Analysis
~~~~~~~~~~~~~~~~~

1. Once you have the results, you can find how hyperparameters affect
   the model performance.
2. Analysis, like overfitting and underfitting, can be performed as
   well.
3. The way you are going to report the results can vary.
4. Results can be visualized, and other patterns in the data can be
   explored.

.. code:: ipython3

    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    %matplotlib notebook
    
    import matplotlib
    import numpy as np
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
    
    
    


.. parsed-literal::

    1. Reporting Based on Best Training Performance:
    
    |                    | 98                      |
    |:-------------------|:------------------------|
    | clump_p1           | 1.0                     |
    | clump_r2           | 0.1                     |
    | clump_kb           | 200.0                   |
    | p_window_size      | 200.0                   |
    | p_slide_size       | 50.0                    |
    | p_LD_threshold     | 0.25                    |
    | pvalue             | 0.2976351441631313      |
    | numberofpca        | 6.0                     |
    | tempalpha          | 0.1                     |
    | l1weight           | 0.1                     |
    | burn_in            | 200.0                   |
    | num_iter           | 10.0                    |
    | temp_pvalue        | 0.006                   |
    | shrink_corr        | 0.1                     |
    | numberofvariants   | 173108.8                |
    | h2                 | 0.2032220489181616      |
    | Train_pure_prs     | -3.3108670615789036e-06 |
    | Train_null_model   | 0.2300103041419889      |
    | Train_best_model   | 0.342791259851683       |
    | Test_pure_prs      | -3.6300773897490755e-06 |
    | Test_null_model    | 0.11869244971792384     |
    | Test_best_model    | 0.29928645894907757     |
    | heritability_model | LDpred-2_full           |
    | sparse             | False                   |
    | allow_jump_sign    | True                    |
    | use_MLE            | False                   |
    


.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAyAAAAHgCAYAAABdBwn1AAAAAXNSR0IArs4c6QAAIABJREFUeF7snQd4VMXXxt/0QiCNTgJI711Aeu8gXUCQIoJKlaIIWPgLKsVCU0ARlI4UASlSQ+899BYSamghJNn0fN+ZS3rbu9mdvZuc+zw8gezcmTO/M/cy7545M1ZxcXFx4IsJMAEmwASYABNgAkyACTABJiCBgBULEAmUuQkmwASYABNgAkyACTABJsAEBAEWIDwQmAATYAJMgAkwASbABJgAE5BGgAWINNTcEBNgAkyACTABJsAEmAATYAIsQHgMMAEmwASYABNgAkyACTABJiCNAAsQaai5ISbABJgAE2ACTIAJMAEmwARYgPAYYAJMgAkwASbABJgAE2ACTEAaARYg0lBzQ0yACTABJsAEmAATYAJMgAmwAOExwASYABNgAkyACTABJsAEmIA0AixApKHmhpgAE2ACTIAJMAEmwASYABNgAcJjgAkwASbABJgAE2ACTIAJMAFpBFiASEPNDTEBJsAEmAATYAJMgAkwASbAAoTHABNgAkyACTABJsAEmAATYALSCLAAkYaaG2ICTIAJMAEmwASYABNgAkyABQiPASbABJgAE2ACTIAJMAEmwASkEWABIg01N8QEmAATYAJMgAkwASbABJgACxAeA0yACTABJsAEmAATYAJMgAlII8ACRBpqbogJMAEmwASYABNgAkyACTABFiA8BpgAE2ACTIAJMAEmwASYABOQRoAFiDTU3BATYAJMgAkwASbABJgAE2ACLEB4DDABJsAEmAATYAJMgAkwASYgjQALEGmouSEmwASYABNgAkyACTABJsAEWIDwGGACTIAJMAEmwASYABNgAkxAGgEWINJQc0NMgAkwASbABJgAE2ACTIAJsADhMcAEmAATYAJMgAkwASbABJiANAIsQKSh5oaYABNgAkyACTABJsAEmAATYAHCY4AJMAEmwASYABNgAkyACTABaQRYgEhDzQ0xASbABJgAE2ACTIAJMAEmwAKExwATYAJMgAkwASbABJgAE2AC0giwAJGGmhtiAkyACTABJsAEmAATYAJMgAUIjwEmwASYABNgAkyACTABJsAEpBFgASINNTfEBJgAE2ACTIAJMAEmwASYAAsQHgNMgAkwASbABJgAE2ACTIAJSCPAAkQaam6ICTABJsAEmAATYAJMgAkwARYgPAaYABNgAkyACTABJsAEmAATkEaABYg01NwQE2ACTIAJMAEmwASYABNgAixAeAwwASbABJgAE2ACTIAJMAEmII0ACxBpqLkhJsAEmAATYAJMgAkwASbABFiA8BhgAkyACTABJsAEmAATYAJMQBoBFiDSUHNDTIAJMAEmwASYABNgAkyACbAA4THABJgAE2ACTIAJMAEmwASYgDQCLECkoeaGmAATyMkEvvzyS8TExGDatGk5GUOO6Xt0dDQqVqyIv/76C3Xq1NGr32XLlsWSJUtQr149vcpzISbABJiApRJgAWKpnmO7mQATMBmB6tWrJ9QdFRUlhIOjo2PC73777TfUqlXLZO2bsuLNmzfjq6++SmhCp9PBzs4Otra24neFCxfG1q1bDTJBnwn0hg0bMHHiRDg5OYk2qO3y5cvj008/FRP2rF737t1D8+bNsXPnThQrVizd6ubOnYt58+ahffv2+PHHH5OV69u3L06ePImpU6eiR48eBpnEAsQgbHwTE2ACOYQAC5Ac4mjuJhNgAoYR+Omnn3DmzBksW7YszQoiIyNhb29vWOUauKtRo0YYPXo0unbtmmVr9BUgP//8Mw4cOCDai4iIEAKARM/BgwdhZWWVJTvUCJAtW7bgyZMn2LNnDzw8PES7N2/exLvvviuE0ahRo1iAZMkbfDMTYAJMIG0CLEB4ZDABJsAEMiCQUoDQN+dHjhwRy2rWrVuHPHnyYNu2bZg9e7b4GRgYKH7XokULjBs3LuGb/gkTJoC+FZ81a5ZorVmzZmLSf/78eZw6dQqenp747LPP0LJlyzSt6dWrFxo2bIhhw4YlfE7f8tPSLprMP336VEQ2qD6K2BQqVAhff/11ppGalALk1q1bmDFjBi5evAgbGxsRTaDohLOzM+Li4jBnzhysX78er169Qu7cudG5c2eMGTNGRBJo8u7g4CDuq1mzJn7//fdUfaEISFIBQgWuXr2Kt99+G0ePHk0QAtQPYnX9+nXRNn0+fPhwEakh0fftt99i165dCAsLg7u7OwYOHIh+/fqhWrVqoKgORVhIzHTs2BH/+9//UtkR78f8+fOLyMuQIUNEmW+++UaID2L70UcfJQgQ4vLdd98JLtTHJk2aYPz48YIBXc+ePRP8jx07BldXVyFe6POkS7Ay6hPVkVTAPXjwwCB/8sPMBJgAE7AEAixALMFLbCMTYAJmI5CWAPn111/x8ccf44MPPkBsbKyY7P7zzz+oXbu2mPjTRJw+b9OmDcaOHStsT0uA0L3z588XS5D+/PNP0KSYxISLi0uq/pLYoXZ3796dECUYPHgwSpYsic8//1y0Q3Z88cUXYgLt5+cnfnp7e2fILqkAef78Odq1a4cPP/wQffr0QUhIiBAXXl5eYjnS4cOHRVurV68WS7WCgoJw584dxC9ZMyQCQgKCGFPdFAUh0XD79m1069ZNtNm6dWs8fvxYiIG2bduKn2vXrsXKlSvxxx9/CMFCUQwSfiQk1ERASEiSUJg8ebIQMyRcGjdujL///huDBg1KECDEgdomkUXlSXwRFxKav/zyi+BL5cn2H374QfybxKSPj0+CAMmsTykFiKH+NNuDwg0zASbABFQQYAGiAhYXZQJMIOcRSEuA0ASYhEJGy4WWLl0Kyregb/zTEyAUAaFv9emiiThN5KnuqlWrpgJNnzdo0EDkLVCS8sOHD0UUhZYRlSpVSgiDFy9eCCFC/9Z3KVNSAUIJ0Dt27MCaNWsS2j99+jT69+8vIiv09xEjRogICUWAkubFpJxApzdS4nNA4kUWTe4pYkBRkbfeekvcRsKDxFDS3AxiSQKNhMLGjRuFGKOEfop4kNCKv9QKkFWrVgmhOGnSJMGUoljkO2IbHwH5999/hU2HDh1KyJW5fPkyunTpIn5HQpI4ko0kwui6du0aOnXqlCBAMutTSn6G+jPnPaHcYybABCyRAAsQS/Qa28wEmIA0AmkJEPq2nqIASS/6N03caekMLbWiP25ubti/f3+6AiTpEh99JvAU3QgNDRUTcxIiNPmNt4PEB03Kqb3g4GCxRIjESN68eTNklVSA0BIiWl6VVFjQsita8kSRlwIFCohlZxTtoQk4RW4o0lO/fn3RhiEREEry3759uxAAVDfVQZGd48ePi6VO8RdN8smWs2fPiiVmJBJILFFkgYQb5bFUqlRJdQSEBAjVRUnn5LuhQ4cKQZJUgNCmA2RjvJgkm0g40TIzipaQXT179hS5Qrly5RImk59q1KiRIEAy61NKfob6U9qDwQ0xASbABLJAgAVIFuDxrUyACWR/AunlgNDENf6iSTHtnLR48WIxKaVv5GlSS0uE4pOt01qCpVaAXLhwQSRIk8jo3r27mPzTz5QXLVmi/BPKb4hfEpSep5IKEBI1J06cEJPmzC4SJStWrBBiiHI3KKJRrlw50eeMtpFNKweE2qKICk3+aSkTfftPF+VcZHbRRJ/yUkggEGsSEU2bNtVrFyxagkV+fPnypVh6RX2gZVOUZ5IyAkLRFkqSj98t7MqVKyL/JWkEhKJRZcqUESZT7grln8TngOjTp/QEnBp/ZsaLP2cCTIAJaIEACxAteIFtYAJMQLME9BEgNPGl5HCKHtAE9NKlSyJXgCbpxhQgBIkmtUWKFBFCgSa/lKBNF+VPVK5cWeRr0HItap+iH9OnT8+QbVIBQpN3SvaO3xWLIiGPHj0S/aGkehJAtGsVtUPRCfr2n5YWUbSC8k9oiRiJKhJJ6V0pBQhFMyiSQTkVJNpoGRa1M2DAALHEipLgKamdllbdvXtXLHWKFzw0Ybe2thZ5NBSV2bdvH8LDw0VEhKJBFAVK74pPQo8XkpRcTruZxS+hSipAKNpBURFaUjVy5MiEHBCKdixYsEA0QfaS8IzfZIByQMieeAGSWZ+ojqQCxFB/avZBYsOYABNgAkkIsADh4cAEmAATyICAPgKElgfRt/WUA0BLrygvIT6fw9gChCa0NDGnJT+0Y1P8RZEOylWgpTskHGgiT0u24reXTa+Lae2CRVENyvmgpGxadkXJ1ySwaOI/c+ZMkXhOE//ixYsLoUN10EW5GbQbWPzypIULF6ZqNuU5IJSrQoKK8kySRnOofaqLlnoRUyrTu3dv0G5gNDkngXH//n0RkaCJO+3UVaVKFdHeokWLxMSfxEiHDh3EbmApr5QCJOXnSQUIfXbjxg18//338PX1FUKFIia0yxXlr9BFifC0hI3EGP2ORFxau2Cl16eUAsRQf/LDzASYABOwBAIsQCzBS2wjE2ACTIAJMAEmwASYABPIJgRYgGQTR3I3mAATYAJMgAkwASbABJiAJRBgAWIJXmIbmYBEArSjz95Ro3Bn2zbYOjuj7dKlKFCjRjILosLCsLlHD7y8dQtWNjYo2bEjGn3/vSjz8u5d/DdoEMKePIGjhwfaL1+O3F5e8N+3D/s++SShnudXr6LD6tUo3bmz2EXo0OTJuP7336K+ah99hBojR2ap11TnqNWjsO3iNjjbO2PpwKWoUSx5P8IiwtBjYQ/cenILNlY26Fi1I77vpvTj7rO7GLR0EJ68egKPXB5Y/v5yeHl4Yd/VffhkTWI/rj66itVDVqNzdaUfk/+ZjL9P/Q0baxt81OQjjGyetX5kCQLfzASYABNgAkxAgwRYgGjQKWwSEzAngdvbtuHM3Lnotm0bHh4/LsRI3+PHUwkQ+qxo06aIiYzE2ubNUWfiRJRo21YIkxIdOqBS//7w37sXvkuWoN2yZcnu1z1/jsWlSmHovXuwc3bGxSVLELBvnxA7VtbWCA0MRK78+bOEgYTH3L1zsW3kNhy/fRyj1ozC8YnJ+0EC5Pid42harikioyPR/IfmmNhuItpWboseC3qgQ5UO6F+vP/Ze2YslR5Zg2fvJ+/E89DlKTSyFezPuwdnBGUsOLxEChcQO5UgEBgcif56s9SNLEPhmJsAEmAATYAIaJMACRINOYZOYgDEIvPTzw/o2bVCgZk08PnMGnhUrot3/JzDThD+ja+fQofBu0gTle/cWxRaXLYt3fHzgUqhQureRSMlbqRKqfPABllSsiG47diCPt7eICMx1dcXI4OBk955ftAj39u9H+xUrxO+X166N9itXwr1UqVRt+D31Q5vZbVCzaE2c8T+DioUr4q9Bf4kJf0bX0GVD0aRME/Suo/Sj7OSy8Bnng0Ju6feDIiaVClfCB40+QMUvK2LH6B3w9lD64TrSFcFzk/dj0YFF2H9tP1Z8oPSj9rTaWPnBSpTKn7ofxvAp18EEmAATYAJMIDsQYAGSHbzIfWACaRAgAfLbG2+g96FDKFK/PnYMGgTPChUQcv++WA6V8irXqxfqTJiADR06oPaECfBq0EAUoehGo+nTUbBWrTQ5hwcFYVmNGuixezfcSpTAv336oFCdOqg5ahSub9iAzd26YdjTp3Dy9Ey4f02zZqg1ZgxK/v8ORXTN8/QU/76xcSOc8uVD8zlz4F66tPiMBMgbn7+BQ58dQv1S9cWyqAqFKuB+0H0RbUh59ardCxPaTkCHOR3EzwallX5QdGN6t+moVTztfgSFBaHGNzWwe8xulMhXAn1+64M6b9TBqBajsOHMBnT7tRue/vQUni6J/Wg2qxnGtByDDlWVfniO9hT/3nh2I/Llzoc5veagdAGlH3wxASbABJgAE2ACCgEWIBoaCbTVJB2IRfvr0/INvphAVggE372Lf1q1wnvXrolq7vn44OKvv6LtmjUZVru1WzfUGDsWherVE+U2tWuHt6ZORf4UeSD0WWx0NLZ17w7vFi1QdfhwUT704UMcGDMGr/z8UKh+fdzetAm9Tp6Eg5tbwudr6tZF/5s3YWNnJ363KH9+1J40CdVGjcKtTZtwYd48dNm1S3xGuRit5rTCtSlKP3yu++DXA79izeCM+9FtYTeMbTEW9Uoq/Wg3rx2mdpqKGkWT54HQZ9Ex0ej+W3e0KNcCw5so/Xj48iHGrBsDv2d+qF+yPjad34STE07CzVnpB31ed3pd3PzmJuxslH7kH58fk9pOwqhmo0T5eT7zsGuU0g++mAATYAJMIHsQoK3X6Uwk2nI7/nDS7NEzeb1gASKPdaYtPXv2DH5+fpmW4wJMQB8Cuv8/VO7skCGo9++/oviLkydxb80aOBYqhKDTp1NVkb9VKxQbMADXpk2DW82aKNCmjShzrGtXVF+0CA5586a658qUKbBxdkaZ8ePTNCk6LAwnundHvW3bEj4PWLUKobdvo9ykSQm/O96tG6rMmQOnIkXEcqeDTZqg0f794vMHwQ8wZPMQ/NtX6cfJ+yex5uIaFMpdCKcfpO5Hq1KtMKD6AEzbPw01C9dEm9JKP7qu6opFnRYhb67U/Ziybwqc7ZwxvkHa/QiLCkP31d2xrV9iP1ZdWIXbL25jUuPEfnRb3Q1z2s1BkTxKP5osaYL9g5R+8MUEmAATYALZiwCdheSZJLqfvXpn2t7keAESHBwsDuuiw8LoVNvBgweLE21TXgEBAeKkXjqJl5RvqVKlMG7cONR6vSyFDp+ig7ToNOD4a+jQofjwww/19iAd3nXt2jVxuFfSevSuQGVBOoH4+vXr4uRmOmmYL8sjkJEPKQKyvEIFdN27FwXr1MG+jz+Ge9myIsqQ0eW3Ywd8FyxA+40b8fjkSRwaNw7dDxxIdcvxKVPw4to1tF6+XCSOx1+6p0/F7lf0u2Nffw1rGxvU/uKLhM/XN2mCulOmoEjjxgm/O/rFF3ArVQrl+/fH/QMHcGTSJPQ4eFB8ThGQCv+rgL2j94olUR+v+hhlC5QVUYaMrh2XdmDBwQXYOHQjTt49iXHrx+HA2NT9mLJ1Cq49uoblA5cnizw+DXkKD2cP8buv//1a7Gr1RbvEfjT5sQmmdJyCxqUT+/HF5i9E/kf/uv1x4MYBTNo0CQfHKf1I7+Ln0PKeu6QWs/8s239kPfuQfaiWAB3SSl8Y0yGoLi4uam/n8rwEC0JEhIaGitN96VRdEh902i2dcpv0InFAEQpvb2/Qyb27du3CpEmTcPjwYXEqLgkQEij0b0OvsLAwXLlyBeXLl4dzJonChraR8j/Oc+fOiVObWYAYg6j8Oug/zvR8mJCEXqsWHp8+LfI/aDeqzJLQ6Zv7PcOH486OHaJsmyVLEvI//qxWDf3PncOre/ew0NsbHuXKwcbBQXS8+vDhqDJ4MK6tW4eDn38unhOvRo3QfP582L4uQzatql8fQwMCkokWyiPZ+u67eOXvDzsXF7RcsAD5q1YV9cYnodcqVgun755GhcIVsGzQskyT0Kkfw1cOBwkR2oZ3yYAlCfkf1aZUw7mvzuHe83vw/swb5QqWg4Ot0o/hzYZjcMPBWHd6HT7f8DmsYIVGZRphfp/5cLBTypBN9afXR8D0gGSihfJI3v39Xfg/94eLgwsW9F2Aqt5KPzISIPwcyn92jNViRs+gsdrgekxLgH1oWr4yapftQ9nzNRkMZbeRoyMgNIBq166NDRs2iCgAXT/99BPu3LmDOXPmpOsLioDs3bsXw4YNE5GTAgUKsACRPXK5PUEgMwFCCeUDfX0tmhZN9jvM7QDfKZbdDxYgFj0MWUBmT/dl+h7Nxt3OVl1jAWJ57szRAuTy5cvo0aMHLl26lOC57du3C/FBP9O6mjZtisDAQFDCeNeuXfHdd9+JYhQBGTRoENzc3EREpGHDhiIiQv/W94pX1CSGZEVALl68iMqVK3MERF8naawcvXTT8yFFGza9/TbeO39eY1arM4cEyNu/vI3zX1p2PzISIPwcqhsTWiqd0TOoJTvZlvQJsA8tf3TI9iHN12gJu6wVK5bvodQ9yNEC5NSpUyKKQeIh/qIlVJ9//rmIbKR30c4HW7duFUtMunTpIoo9efIEQf+/HWnJkiXx+PFjfPXVV2JpxoIFC/QeN/ECRO8buCATYAJMgAkwASbABJiAWQiwADEce44WIBQB6dmzJ3yTLFHZsWMHZs+enW4EJCnqVq1aiWhJuXLlUnmAktbp8zNnzuidUM4REMMHcna6U/fsGda3aiW6FProEaxsbOCcL5/4d++jR2Fjb5/Q3ZTf+jw6dQpXli9H059/zhKS8OfPsbV3b1Aie55ixdB+9Wo4ursnq5M+29K9O+JiYxETFYVqw4ah6tChycps6twZL+/cSYjCHPj0U9zeulX0wbVECbRavBiObm7i/l1DhiDw7FnERUejfN++4iySZyHP0OpnhcWj4EewsbIR52vQdXTCUdjbJrJI2eFTd09h+bHl+PmdrLGg0857/9ZbJMMX8yyG1R+shnuu5Czos+4LuiM2LhZRMVEY1nQYhjZKzqLzL51x5+mdhEjOp+s/xdYLW8UWvnnt82LNsDXwzO0p7h+ybAjO+p9FdGw0+tbtiwltJmTJn3yz6QjI/ubVdD3JuTWzDy3f97J9yBGQrI+ZHC1A4nNANm7ciNKvDz3TJwckHjstx5o4cSJatmyZyhOU0N68eXMhQPRdTiU7qUn2msmsD9ecV8Phr7+GvYsL3hw3LqHzdPaGta2t+LepfLj/00/FTlZ0MOHx779H+IsXaDx9ejIHxERGiq1mKcE8MiQESytVQp8jR+BSuLAoR4cQXl+3Dk8uXEjIQ/HbuRNFmzUT9u//7DNRjuq9snIlbm7ejI6rVyMqLAxLKlQQp6+7Fi+e0ObXm78Wid3jWieyoPM7bG0UFqa6Pl33KTxyeYhDDb/f/j1ehL7A9O7JWURGKywoST0kPASVvq6EIxOOoLCbwoIOMqSk9gv3LiTksuy8tBPNyjUTSe6DFgxC/vz5MbPHTKw8vhKbz2/G6iGrERYRhgpfVRAnuBfPm8jCVH3letUTMNUzqN4SvsNQAuxDQ8lp5z7ZPpQ9X9MOaeNZkqMFCGEcO3YsaDu1GTNm4MGDBxg4cCC+/fbbVLtgHT16VGzTW6FCBURFRWHx4sXiD0VMKAn92LFj8PLyQpEiRfD06VN8+eWXiIyMFGX0vWQPaNkPrL4cuFwigXgB8tTXF7aOjnh89qw41ZxOLd87ahSiw8MRERuLLqtWIV+FCvD38cGpWbPQ9d9/QffSrlJBt2+LnzVHj0aNkSP1wru4bFkhAFwKFULIw4dY06QJ3n99oGFaFVDU5q/q1fHusWNCgJAgWd+mDVouWoQtPXummQhPp56TQGm/YgWurFolREjnjRsR8fIlVr71FvocOwYnD49UAsT3vi8c7RxxNuCsOCCQTj4ftXoUwqPC4WTvJHa7KluwLHyu+WDWf7Pw78h/QeKFdqa6/eS2+Dm6xWiMbK4fi7KTywoBUMitEB4GPUSTWU1wbapyKGJaF0Vtqn9THcc+PyYECAmSNrPbYFG/Rei5sGeqZHp6Dn/a+BPOPD+DlUNWYtXxVVh5YiU2frwRL3Uv8dZ3b+HYxGNCBPGlPQL8HtWeT9RaxD5US0x75WX7UPZ8TXvEs25RjhcgdA7I5MmTcfDgwVTngFSvXh2//fabOOtj9+7dYocsEimUZE57P48cOTLhHJAlS5Zg6dKlIg8kT548Igmdtvj1SDKBysxdsge07Ac2s/7z56kJJBUgdL5G502bxLkaEcHBYovcOCsr7F24EKG7d6Pzhg2pBMjdnTvRc98+RL56hT/KlsVHjx6J08dXNWwofpfyajJrFoq1aIG5bm4YERQkPqZv9ue5uyf8O+k9wQEB2NC+PYJu3kTjmTNRfdgw8fG+Tz4RW/Dmr14d6e3EtaFjR5R75x1U6NtXLMHa1q8f/PfswYOQENxtUA+nHt1HyKsQuOR2QduWbRFTNgZFixcFCRA6o2PT8E3ibI5gXbDYZpciIbsv78av+3/F+o/WpxIgOy/vxL6x+/Aq/BXKflEWj2Y9gp2tHRpObyh+l/Ka1WMWWlRoAbeRbgiak8jCfZR7wr+T3hPwPADt57THzSc3MbP7TLEMi65P1nyCRqUboXrR6mnu5kXPYdPvmmJws8F4r957iIqOQr8/+mHPlT0IiwzDT+/8hCGNhvDjoVEC/B7VqGNUmMU+VAFLo0Vl+1D2fE2j2LNkVo4XIFmiZ+SbZQ9o2Q+skXHliOqSChDvpk1RqX9/0W+a+O8dORIvbtxAeHg47G1t8f7Vq6kECImNuq9PHP+jfHn02LULub28MmWXVIBQ4bkkQF68SPe+kAcP8E/nzuiyZQtCHz7E4S+/RJfNm0E7caUlQI5NmwbKV3l7wwaxmcP9w4dx6OefsVT3Crv+2wUUiUVMSQCU4hEJ2AXYIfp2NMrUKoMafWugdY3W6F9PYUET/5GrRuJG4A1RF+VQXP3maioBQrkWk9orp5aX/6I8dn2yC14embNIKkDoXhIgL2anz+JB0AN0nt8ZW0ZsERGTLzd/ic3DN4uzQ9LaTvibLd9g74W92PXZLtja2uLwzcP4Zd8vWDpwKV6EvUDDGQ2xfdR2lMhXIlO/cQH5BPg9Kp+5sVtkHxqbqPz6ZPtQ9nxNPlHTt8gCxPSM9W5B9oCW/cDqDYILJhBIKkBKdOiAst27i8+2DxiAAjVqoOqwYTi6bRsujxiBIX5+qQRI0vyRJZUqiaVZlFeRWQRE7RIssmnHoEF4o1076J48wdFvvhGJ5pSvEhYYiML16qGXj4+w3XfpUpxfuBA99+xJOBTxn8GDMXrHVjywe4aoelGAcxqDIAywOWSDXNG5MGfFHPRvogiQAX8MQI1iNcSSKprk0xIpv+/9UgmQpPkjlb6qhH9H/CvyKjKLgKhdgkU2DVo6CO0qt8OTV0/wzb/fiGR5ylcJfBWIeiXrwWe8wmLp4aVYsH8gRu5cAAAgAElEQVQBfmj2A+q+WVdshz1sxTDULVEX/d7qJ8pQXW0qtkHPN3vyk6FBAvwe1aBTVJrEPlQJTIPFZftQ9nxNg8izbBILkCwjNF4Fsge07AfWeKRyTk3pCZB/unQRS5dKdu6MjR99hOc7d6oSIJkR9Bk/Hk6enolJ6M+fo/GMGcluo9PQHT09YefkJJLUV9Spg07r1yNf5coJ5VJGQOh09X1jxqDX/v0JO3tR4XoVy+Jk8E1Et4oFrDOwjj7ebY2qxarizIEzomCX+V3ETlHdanYTuR5LjyxVJUAyYzH+7/HwdPFMSEKnXbFmdE/Ogk5UpzKUg0JJ6nW+rSOWgVX2SmSRMgKyw3cHxqwdg71j9uL+rfuoVq2aECDTt0/H1UdXsWTgEoRGhOLNaW+KhPQqXlUyM5U/NwMBfo+aAbqRm2QfGhmo5OpomfCFGzocPX0bb9UsgSqlnUQ03JSX7PmaKftirrpZgJiLfBrtyh7Q/NLVkPPTMSU9AfLg6FFs798ftrlywblGDbz4/9wJNRGQzHpOSeWUPB7s7y+24e24dq1ICKdlU+cXLEDr33+H365d8Bk7Vrzo6T+A6sOHo+qQ5LkKKQXI76VKISYiQggXugrXrYs3xo1DufLlENMzJu3IR0pjwwDrtda4fvW6OHfn6K2j6P9Hf+RyyIX2ldtj+fHlRhUglFROyeOUvE7b8K4dulYkhJ/yOyWiF7/3/x27Lu/C2LWJLIY3G54qbyOlACk1sRQioiPgmctTbITRpEITLHxvoUhaH7h0IC4/uIw4xGFg/YEY33p8Zi7jz81EgN+jZgJvxGbZh0aEKbmqR8+i8dm8QDx4Ek0ZiwCsUDifLaYPz4+CnqbbIVH2fE0yVinNsQCRglm/RmQPaH7p6ucXLZfKDj78ZMwnmP/vfEQ1jtIbtZ2PHYZ1HIaffvxJ73u0WjA7+FCrbGXYxf6TQdm0bbAPTcvXVLXTF199v3qIh09JfCS/Cue1xbIphUwWCZE9XzMVQ3PWywLEnPRTtC17QPNLV0PON9AUS/NhbGwcIqLiEBGZ+LNJk0q4U+Q6oCbH+jZQ2K805vx2ShU5U4flVRnzunBcXCzu3bsntvG2ts5o/Zn62o25CsHoCxqMXKGRq4O+7GhMB9wLgLeXN6yt07dCi2MvfkTp21d9R6C5fKGvfUhhYGxsLAL8/eFdtKjqZ1DrfTW2fXo/GHo6Iyv2+T+KxIr/Uu9gGN/0z5/kR5XSjnpaoq6Y7PmaOussozQLEA35SfaAtrTJq4ZcpRlTjOFD+hYpOgYJoiA8MjaZQAiPjENkVBzoZwR99lpAJP99vKCIFeUSyycKDao3KvUXVTi8sg4i33oEZL4hVSL3e4D90YKo3+e4ZnzBhjABJsAEmIB2CLzTPBeGdlOW+xr7kj1fM7b9WqiPBYgWvPDaBtkD2hiTV3Pjo8mz760I3H8SjSL5bFGppIPJQq6y+xofLUhvMk9CQBcegxu3/FGgoBcio62SCIdEEaEIhySRhxQiguqJjZXdu8T2jv3dFLoqt1VHQJwulEDdHvvMZzi3zASYABNgApol4JXPCn9N8TaJfbLnaybphJkrZQFiZgckbV72gLZ0AULJZ+NnP8L9p4mz5yJ5rTFzVEGTJZ8lRAsSogCx6XzbHx8FMDyaQMIjJ1w3jkzB/eBliGuufw6I1R47FMnzHkrX+zInIOI+MgEmwASYgEoCjnbAttlFVd6lX3HZ8zX9rLKsUixANOQv2QPakgUICYEunwYgODS1A12cgbnjCiIyynjLiuIjCBRNMGe0QEPD1WimhL30w/HVTYFesXrvgoXV1qjbywdOrsWMZgdXxASYABNgAtmHgLMD8O9PLEC06lEWIBryDAuQ9J1BgiM4NFZstUfLrU5eCsGukxEa8h6bog8BO1sgt7M1cueyQR7x01r8e9ncXrgauB8xLaIyPQfEbo8d6pVuhsVLNqdqMi6NoBH9Kq3fJ72ZxlfKK+WvEkqkE5jKqJ302qd2KQH2xo0bKF26tEiAzdTeOGWzSePYm3Zn0rc3bS9naPPrJlK2lCZzqj6DwB/ZpU/f461Mtx+inYz7ntretPseExOLu3fvolixYqkSmNOzN2Xb8W1lbG867NPlm1759AGnb2/yujK1N0MfJn6YtJjaMRfvwjSbSvHLRHvTNoyWuz548BCFChVKtZFAUruS3W3Qe0C/5y0zvvQQpIc4fXvTb1ufupK/L9N4/2T0nn1tVFrtvAiOxYnLOoSFZzBo9PnPJY0yNcrYYdboQgbenfFtsudrJumEmStlAWJmByRtXvaA1loEhCYkz17GCIHx4Gk0HgQqP+nft+9HikRpvoxLgCb/Ls7WiWIgV5K/0+9zJf1cEQ25nKxgZ2uV4Y4/aq0MCgpCnXp14KfzQ2S9yHRPQrc/Yo/iTsVx4ugJuLq6qm1Gk+W19hxqEpKGjWL/adg5eprGPtQTlBGKkdi75h+J4746HPcNF3839Gpc3QH7z6b/ReRPo/OiahlnQ6vP8D7Z8zWTdMLMlbIAMbMDLFGA0De2mw+G4EZAFEp726FTQxe9ty6MiY1D4PMYRVgERiWIjYBHUfB/nMYWSRryj6lNcbS3eh0dSC4ClIgBiQGbVJ87O8Th6uULqFFDOUXbki8SIX369cF/23cirqg14t6IAuwBRAJWd+xg5R+L1m1bYeWylXBzc7PkriaznSc/lu1K9p9l+4+sZx+a1oevwmJx8jIJDh2OXNAh1MBoR4nCdviouzuql3EQX4DRl5Y9J97Hs5epd1HJ62aNNdOKmGxTGhYgWR8zLECyztBoNcge0Ia8dC/eCsfoHwKThX9pH++fx+ZH5ZLKftuUPP3ouRLBuP8kKiGKcfl2BEJ0xg+zGs0BSSqi7fzjlwdRFCBBBLxePiT+nSAMlL+LaIGTNWxssrKzubreGOJDdS3IL33z5k3MnTcXK9dvQ1hYCJydXdCnWzuMHDFSnHye3a7s6MPs5qOM+sP+s3xvsw+N60MSBrfvRylRjkvhuHjL8OXSvVvlwTstc4sv4NK6xEnocwPFl5lxiIMVrFAkvy1mjMiPAh58ErpxPWvc2liAGJdnlmrTugChyEeL4ffS7WOZora47q+NKEbt8naYMrQAHOyNe7Bblhxsgpv5P04TQJVcJftQMnAjN8f+MzJQM1THPsw69LDwWJy5Go7jl5Qox4tXhu3tXsrLDh91c0e1MvpvqU+C58INHY6evo23apZAldJOJot8xJOSPV/Luoe0VwMLEA35RPaAVvvS/XvXc/y6MUQqMY881qhe1hHVyziKnwU9bcSLJTo6Gq1GPkjXlp1zCsPW1nTffkiFkEFjan2oFbvZjkQC7EPLHg3sP8v2H1nPPlTvQ5r0BwRGv87l0OHMNcOjHO+2pihHHpGPaOgl24ey52uGctHyfSxANOQd2QNa7QPbaaw/QnTGBWZvZyXWcwqRUdYRtMZT3yVMmw4EY/bqoFQGfdLHDR0b5DGuoRqtTa0PNdqNHG0W+9Cy3c/+s2z/sQDR338RkbE4dyNCyeW4qBP5nIZcZYra46NubqhSSv8oR2btyH4OZc/XMuu/JX7OAkRDXpM9oNU+sC2G+SPWgBSOiiXsE6IYZYvZw9nR8G85UrqLIiHfLn2Bm/ciUcrLHhMHuOeIyEc8B7U+1NBwZ1NeE2AfWvZQYP9Ztv9YgGTsP8qxUHI5dDjmG26ws99rlwfdm+cReZKmuGQ/h7Lna6ZgZu46WYCY2wNJ2pc9oNU+sB3G+CMsg/ePoz3w19eF4emqLJPiy/QE1PrQ9BZxC2oJsA/VEtNWefaftvxhiDXsw0Rq0TFxuHgzQgiOoxd1CDBwd8ryxSnK4Q76AlLGfEC2D2XP1wwZ11q/hwWIhjwke0CreWBpZ6s2owIypPVxt9zo3txdQ0SzvylqfJj9aVhmD9mHlum3eKvZf5btP46AQJy/deKSEuU4cNbwddYDOriie7PcRl3loO/okv0cyp6v6cvBksqxANGQt2QPaDUP7EffP8rwwCCKd+ya56X3eSAawm7RpqjxoUV3NBsbzz60bOey/yzbfzlRgNB5XFf9IpVlVRd1uHkvyiAnUnSDohwU7ZAR5cjISNnPoez5mkEO0vhNLEA05CDZA1rfB3bPyVBMW/IsXVI21sDPY/KjYgnlHBC+5BHQ14fyLOKW1BJgH6olpq3y7D9t+cMQa3KCD1+GxODUlXAc89Vh3+kwxBq2Sy7e7+SKrk1yw8mIuZyG+CzlPbJ9KHu+ZgxGWquDBYiGPCJ7QOvzwFICWp8vkm93u3uel8EnoWsId7YwRR8fZouOZuNOsA8t27nsP8v2X3aNgNA2uTcCokSUg5LIL9+JNMhRtFPVh13dQBvImDvKwREQg1yo2ZtYgGjINVoTIFHRcWg9Mnnex9YfvTT3zYeGXCjdFJ78SEdu9AbZh0ZHKrVC9p9U3CZpLLv4MET3+jBAXx18zoRBF2HAtpUAPujshs6NXeDkYJodq0zhRNk+lD1fMwUzc9fJAsTcHkjSvuwBndkDO3jaQ9y+n7g2dPaY/KhcipdZaWjI8AFaWnKGgbZk9hwaWC3fJokA+08SaBM2Y6k+pCjH3UeJhwHSGR2GXHTq+Idd3UHnc1jqJduHsudrluqXjOxmAaIhr8oe0Bk9sCnzPlrWdsbnA/JqiBabQgRkv3SZuvEJsA+Nz1Rmjew/mbRN05Yl+TCcDgO8pmyTu/9MGIJCDEvmGNrFDW83doGjveVEOTLyvmwfyp6vmWbkm7dWFiDm5Z+sddkDOuUDS/9euPElLt6KwLW7yXfF2PtLUQ2RYlPiCch+6TJ54xNgHxqfqcwa2X8yaZumLa378MHTxCjHicuGHQZYs5wjSHSU8rbcKAcLENOMf3PVygLEXOTTaNecAuTQhXBM+S3tna7+m+MNO1s+WFBDQyXBFK3/x6lFZlqziX2oNY+os4f9p46XFktrzYd07hZ9EUjJ4wfPheHx8xiDsH3c3Q0dG7jAIZtEOViAGDQMNHsTCxANucZcAqRy5cpoM/pRuiR2zS0CGxsbDZFiUzgCkn3GgNYmP9mHrJyesP/kcDZlK1rw4ZMX0Th+KVwsrTp83rDDAGtXcMSQLm4oUSR7RjlYgJjyKZBfNwsQ+czTbdFcAuTA1cLYcjj9g4g6NbTH6N4FNUSKTWEBkn3GgBYmP9mHpvyesP/kMzd2i+bwYUxMHC7foVyOcBw+HyaSyQ25hvdwR4cGLrC3y9mrFGT7UPZ8zZCxofV7WIBoyEOyB3T8Azvpz3yIzODdZ28H7JjNOSAaGioJpsh+6WqRgaXbxD60bA+y/yzbf2S9LB8GvYoB5XDQ0io6DNCQq05FJcrxRuGcF+XgCIghI0a797AA0ZBvzCVAPvsjH2Iz2C7c2grYPZ8FiIaGCgsQLTrDQJtkTX4MNI9vy4QA+8/yh4ipfBgbS4cBRiZEOehgQLWXjTUwvKc72r7FUQ4WIGpHj7bLswDRkH/MJUBmrs+PwKD0FUh+N2us/tZLQ6TYlHgCpvqPkwnLI8A+lMfaFC2x/0xBVW6dxvRhSFgsTl3RCdHx37FQgzpSr4oThnR2Q9GCdgbdnxNvMqYP9eEne76mj02WVoYFiIY8JntAxz+wsY6l8Nn8l+mSmDXCHTXK59YQKTaFBUj2GQOy/+PMPuS00RP2nzb8kBUrsuJDOgzwzoMoITiOXtTB95b6wwDtbIERPT3Qum4u3nHSQEdmxYeGNCl7vmaIjVq/hwWIhjwke0DHP7BVq1ZF188eICSNjTdyOwP/zPSGlVXOTnDT0DBJZorsl65WOViyXexDS/aevPwBy6akbevVPoO68FicuR6OE77h2Hk8FBFRGaxhTqfrDao64YPObvAuwFEOY4wOtT7Mapuy52tZtVeL97MA0ZBXZA/opA9s4ItYvPvlw2Q0iuSzxqxRBVHAw1ZDlNiUpARkv3SZvvEJsA+Nz1Rmjew/mbRN05Y+PrwXGCWSx4/5huP0VfWHATo6WGF4d3e0qpsLtjb8hZ6xPamPD43Zpuz5mjFt10pdLEC04gkAsgd00gfW/3EM3p+qnAXyaT8PFMlni0olHTjyoaHxkZYpsl+6GsdhkeaxDy3SbQlGs/8s239kfVo+pMMAz99QdqzafTIMwaGxqjvauIYzBndyRZH8HOVQDU/lDbKfQ9nzNZU4LKI4CxANuUn2gE76wHYe/wCh4UoYee8vvOOVhoZFhqbIfulaChdLspN9aEneSm0r+8+y/ZdUgBQuWgknr0Ti2OtIh9qe5XK0wrAe7mhRm6Mcatlltbzs51D2fC2rfLR4PwsQDXlF9oBO+sC2HHE/gQQLEA0NikxMkf3StRwylmMp+9ByfJWWpew/y/VfdEwcLt2KwFHfMOw8+hJBoTaqO9O0FkU53FAoLy9VVg3PiDfIfg5lz9eMiEozVbEA0YwrzLsEiwWIhgaCClNkv3RVmMZF9STAPtQTlEaLsf806ph0zHr+kg4DpFwOHQ6cTWPnlUy6k9vZWkQ5mtdyhg3ncmjG+bKfQxYgWXc9C5CsMzRaDbIHNEdAjOY6s1Uk+6Vrto5m44bZh5btXPaftv0XExuH63fpMEAdfE6Hwf9xtGqDW7zpjEGd3FDQk6McquFJukH2cyh7viYJo9RmWIBIxZ1xY7IHdFoCxCOPNdZ9z4cOamhYZGiK7JeupXCxJDvZh5bkrdS2sv+057/g0BicuhIuohy7T4SpNtDNRYlyNKnBUQ7V8Mx0g+znUPZ8zUxYTdosCxCT4lVXuewBnZYAWfZ1Id6xQ53bzFpa9kvXrJ3Npo2zDy3bsew/8/uPDgO8dY8OA6RlVWG4ERCl2qiapcMxpl8xFMrroPpevsH8BGQ/h7Lna+YnbHwLWIAYn6nBNcoe0PEPbOUqVdFm1ANhNyegG+w+s9wo+6Vrlk5m80bZh5btYPafefwXFh4rzuMg0bHtcKhqIyjaP6y7OxrVcAbiYnHu3DlUq1YNNjbqE9FVN843GJ2A7OdQ9nzN6MA0UCELEA04Id4E2QM6/oF1dKuAET88YQGiobGgrymyX7r62sXl9CfAPtSflRZLsv/keIWiHAGPo4XgOHROh4u3IlQ33PatXBjQ0RX53JLncrAPVaPU3A2yfSh7vqY54EYwiAWIESAaqwrZAzr+gR2/OF9CFzgCYixvyqlH9ktXTq9yVivsQ8v2N/vPdP6LiIzFuRsR4jDAfw+FIDpGXVt53WwwvIc7GlR1grV1+qePsw/VcdViadk+lD1f0yLzrNrEAiSrBI14v+wBzQLEiM4zU1WyX7pm6ma2bpZ9aNnuZf8Z13+PnkULwXHovE4ssVJ7dWjggv7tXeHpqv9SKvahWsraKy/bh7Lna9ojnnWLWIBknaHRapA9oFmAGM11ZqtI9kvXbB3Nxg2zDy3buey/rPkvKjpOLKci0bH9SAhCdHGqKizgYSNyOepVyTjKkVGl7ENVyDVZWLYPZc/XNAk9i0axAMkiQGPeLntAswAxpvfMU5fsl655epm9W2UfWrZ/2X/q/fc0KBonLoXjyEUdjlxQfxhgp4YueK+9Kzzy6B/lYAGi3k+WdIfs51D2fM2SfKGvrSxA9CUloZzsAZ1SgJT0ssNvEwtJ6Ck3YSwCsl+6xrKb60kkwD607NHA/svcf3QY4JU7ymGA/x0LxdMgdckchTxtxLkcdSsZHuVgAZK5nyy5hOznUPZ8zZJ9k57tLEA05FXZAzo6Ohr//HcJv2x1FxRWTy2E/B52GiLCpmRGQPZLNzN7+HP1BNiH6plp6Q72X9reeBkSg5OXw3H0og77Tqs/DLBLExf0a+sKt9zGiXKwANHSU2N8W2Q/h7Lna8YnZv4aWYCY3wcJFsgc0JToN272Qzx4GgtA2R2kSF5rzBxVEAU9k29RqCFEbEoKArJfuuwA4xNgHxqfqcwa2X8K7djYONykwwB9ddh1IhT3AqNVuaFIPlsR5ahT0RFWVunvWKWqUj0Lsw/1BKXhYrJ9KHO+pmHsWTKNBUiW8Bn3ZlkDmvZT7/JpAILTOLspTy5g4wxv6f8BGJdkzqlN9ks355CV11P2oTzWpmgpJ/svRBeL01fCccxXWVql9urWLDf6tskDVxfTRzkysi0n+1Ctz7RaXrYPZc3XtMrbGHaxADEGRSPVIWtAn70agrFznqdr9Q8jPVC9nIuResXVmJKA7JeuKfuSU+tmH1q253OS/+jLK7+HUTh+KRx7T4XiZkCUKud5F7AV53LUKi8/ysECRJWrLK6w7OdQ1nzN4hyhwmAWICpgmbqorAHd70t/3H+afm+K5AWW/a+oqbvL9RuBgOyXrhFM5ipSEGAfWvaQyO7+00XE4tz1CBHl2HIwRLWzejTPjXfb5EGeXOaNcrAAUe06i7pB9nMoa75mUU5QaSwLEJXATFlc1oBuOdwfMZT6kc5lYw3smscCxJS+Nlbdsl+6xrKb60kkwD607NGQHf13/wnlcoTD50wYfG9FqHJQ8UJ2IspRvayDxSzlzY4+VOW0bFBYtg9lzdeygWvS7QILEA15V9aAbjPSH5EZ5Afa2wI75rAA0dDQSNcU2S9dS2BiaTayDy3NY8ntzQ7+i4xSDgOkKMf6va9UO6RXqzzo0yoPXJytVd+rhRuygw+1wNGcNsj2oaz5mjmZmrptFiCmJqyiflkD+suFj3DofGS6ljWoao//DS2ownIuai4Csl+65upndm6XfWjZ3rVU/z15ES1yOfafCcPpq+GqnFCiiBLlqFracqIcGXXQUn2oymnZvLBsH8qar2Vnt7EA0ZB3ZQ1oOv+j1cgH6fZ855zCsLXlrXg1NDQ4AmIJzjDQRtn/cRpoJt+WDgFL8V9MTBwu3YkQS6s27HuFiKg4VT7t0zoPKNLh4mSZUQ4WIKrcbXGFZT+HsuZrFucIFQazAFEBy9RFZQ7oTQeCMXt1UKoufdLHDR0b5DF1V7l+IxGQ/dI1ktlcTRIC7EPLHg5a9t+LVzE4eUmHg+d1OHxepwp0aW87fNzdHVVKZY8oBwsQVe63uMKyn0OZ8zWLc4aeBrMA0ROUjGKyB3RERATafvI4oWsc+ZDhZeO2Ifula1zruTYiwD607HGgJf/RYYDX/SPF0qrNB17hxasMdhtJA3u/tnnwTss8cHbMflEOFiCW/ZxlZr3s51D2fC2z/lvi5yxANOQ12QOaHtiWI+4nENj7Cyeea2g46GWK7JeuXkZxIVUE2IeqcGmusLn99yqMDgPU4dB5HfaeClPFp2wxewzr7o6KJewtZscqVR3Us7C5fainmVwsAwKyfSh7vpYdnc8CRENelT2gkwoQR3srbPvZW0M02BR9CMh+6epjE5dRR4B9qI6X1krL9h8dBnj7Ph0GqMPWw6F4+DSDLQ3TgNW/vSt6Ns8NpxwW5eAIiNaeHOPaI/s5lD1fMy4tbdTGAkQbfhBWyB7QSQVIqzq5MKG/p4ZosCn6EJD90tXHJi6jjgD7UB0vrZWW4T9deCzOXAvHkQs6bD8aqgpBhTfsRS5HhTccVN2XkwrL8GFO4mmOvsr2oez5mjmYmrpNFiCmJqyiftkDOqkAmTIkLxpWc1ZhLRfVAgHZL10t9Dm72cA+tGyPmsp/AY+VKMeOo6Ei4qHmGtTRFd2a5YaTQ87K5VDDKGlZU/nQUHv4PvUEZPtQ9nxNPRHt38ECREM+kj2gkwqQjTOKwNXFRkM02BR9CMh+6epjE5dRR4B9qI6X1koby390GOD5G+E4clGHTftDVHWzckkHfNTNDeWKc5RDFbjXhY3lQ0Pa5nuMQ0C2D2XP14xDSVu1sADRkD9kD+ikAmTPfO8cnYSooWGgyhTZL11VxhlY+ObNm1jwy3z47N6OVyEhyO3igiYt2uLDj4ehVKlSBtaq3duyow+1S9v4lmXFf4+eRePEJR12Hg/F5TvpHw6bltWD33ZFt6a54WDPUY6sejUrPsxq23y/cQjI9qHs+ZpxKGmrFhYgGvKH7AGdVIDwDlgaGggqTJH90lVhmuqiQUFBGPReH2zdsRNvV7JGj0pRcHMCgnTA37522OQbi/ZtWuGPv1bCzc1Ndf1avSE7+VCrjE1plxr/RcfEwfdWBI5e1OHvPa9UmUWnjn/UzR1litqruo8LZ05AjQ8zr41LmIOAbB/Knq+Zg6mp28zxAiQ4OBhffPEFDhw4gFy5cmHw4MEYMGBAKu4BAQEYM2YM7t69i9jYWPFN7Lhx41CrVq2EssuXL8fChQsREhKCBg0aYOrUqXB1ddXbh7IHNAsQvV2j2YKyX7qmAkHio3H9OiiCO1jcPQqF0jgL82Ew8P46e9xHcew/fDzbiJDs4kNTjQ2t15uZ/56/jMHxyzrsPRmG01fDVXVnSBc3dG2SG/Z2Vqru48LqCGTmQ3W1cWlzEJDtQ9nzNXMwNXWbOV6AkIgIDQ3FzJkzcf/+fSE+vv/+ezRu3DgZexIVz549g7e3slRp165dmDRpEg4fPgx7e3vxkwTKH3/8gWLFionPqNzPP/+stw9lD2gWIHq7RrMFZb90TQWia6d2CL+1B5v7R8I2g1Sk6Big05/2cCzZHBs2bzOVOVLrzS4+lApNQ42l9F9MbByu3Y3EMV8dlm8PVmVpjbJKlKOkF0c5VIHLYmF+BrMIUAO3y/ah7PmaBhAb3YQcLUBoANWuXRsbNmxAmTJlBNyffvoJd+7cwZw5c9KFTRGQvXv3YtiwYSJyUqBAAYwdOxb58+fHZ599Ju7z8/ND+/btcezYMeTOnVsvx8ke0CxA9HKLpgvJfumaAgblfFSsUA5+n8ekGflI2SZFQop/Z4PLV66hZMmSpjBJap3ZwYdSgWmsMfLf4WPnEWFbGvvPKlvlqrkoefztRhzlUMPM2GX5GTQ2Ufn1yfah7PmafKKmbzFHC5DLly+jR48euHTpUgLp7du3C/FBP9O6mkq99YMAACAASURBVDZtisDAQERHR6Nr16747rvvRLFOnTrh/fffx9tvv51wW7Vq1fDnn3+iatWqenkyfkCTGHJ2Nv2WuPTAthn9SNi2a24RvWzkQtoiQD68ePEiKleuDBsby9zFbPy4Mbi3/1esfVf/rUZ7rLCDd+OPMXPWD9pyiAHWZAcfGtBti76FDgO8JQ4DjMDqna8QHhmnd39qlnPA0C6ueKOwnd73cEHTEuBn0LR8ZdQu24c0X7t+/TrKly8vZb4mg6HsNnK0ADl16pSIYhw/fjyBOy2l+vzzz0VkI70rIiICW7duFUusunTpIoq1aNFCLLsigRJ/NWzYENOnT0e9evX08mu8ANGrsJEKjV+cT9Q08/0nRqqRq2EC6gj079MNX9W5ix766XRR+dpzwLQjnvhnzijEWdkgztpW+WlFP5U/SPh3ys/TKpdYhu7jiwmkJBAeaYUbD+xw4Y4Dzt12VAXo7bohqFtOl+HyQlUVcmEmwAQ0QYAFiOFuyNEChCIgPXv2hK+vbwLBHTt2YPbs2elGQJKibtWqlYiWlCtXTkRAKIGdfsZf1atXx9KlSzkCYvj45DszIWCKb33o291PPvHB9u134Oxsh8WLW6FGjQLJLAkLi8I772zF7dtBsLGxQvv2JfDddw1Fmbt3gzF48E48faqDu7sj/vqrDby8cmPfvgCMG7cfiIsFYiJx9UYIVk6Nw+c/fY15HUPRUlkFqde18xow4h/gmrLi0ahXHKwAa0XAIPZ1VMbKGnD0BGwdgZgoIPwpQOWsrIDocMC9DJCroHJPsB8Q9lh8HudWCvCsAJAgsrFL8tMGsLYT7ZBYehT4FAUKe8Oayrz+fcqfcQm/t1XsS6dcmr9P2r64l0VWZoOGngP/x7RNbjjW7Q3B8+DYzG5J+PzN8kqUo1ghjnLoDc2MBU3xHjVjd3Jk07J9yBGQrA+zHC1A4nNANm7ciNKlSwua+uSAxGOnaMfEiRPRsmXLVDkgtFtWu3btNJ0DQsvIWo18ILrD2/Bm/WEyRw2mWPe6bdttzJ17Btu2dcPx4w8xatReHD/eN5UAoc+aNi2KyMgYNG++FhMn1kHbtiXQo8dmdGhXHP07O2PvtrNYstIfy8bcAV7cAIJuAMF38TzUAaW+m4B7X0xFo1+i8FlTqI6AzPABTo02B/Xs0OZrkZWeiEkmllKKp3T+rVYUmaJ8vN0GiqzwyFicux6BQ+fCsO1IqCpHj+jhio6N8sDWhnesUgVOA4VN8R7VQLdylAmyfcg5IFkfXjlagBA+Sh7X6XSYMWMGHjx4gIEDB+Lbb79NtQvW0aNHxTa9FSpUQFRUFBYvXiz+UMSEktDjd8FasmSJ2AVr8uTJoG/QtLwLVkhYFDqNe8gCJOvPkdlqyOil6+f3Em3arEfNmgVw5sxjVKzoib/+aieiGhldQ4fuRJMm3ujdu7woVrbsYvj4vINChVyS3xYbA7zyF8Ji1MRLqFToKT6ofw4VR1TFjvcXwdv1GeLiANfJ3yB42hfJ7l10rA723yqBFe+uwrgtwK1n1tg4QP9vmHv8BRRzB2Z1NBt6bljzBPQTWVFxtgiNsMHTV1bQRdogGraIiXv9B7aIpr/H/3z9d/o8r4cDKpV2gUsuOzx+8gwFCnnB2tZeXWTKEBGWRZGlebeZwUDZk1czdDHbNynbhyxAsj6kcrwAoXNASCwcPHgw1TkgtITqt99+E2d97N69W0RHSKTQtrtly5bFyJEjU50DsmDBArGtb/369TFt2jRNnwNy6nIYPp1HS0k4ApL1R8k8NWQmQN544zccOtQb9esXwaBBO1Chgifu3w/Bvn3+qQzu1ascJkyogw4dNmDChNpo0MBLLJdq3mQ5po9zR62i95UoRnwk4+VtsZQqSOeIGj+Nxu6hi1DC8zn6rOiDOkX9MarhIWy4WAnd/uyPp1O+gmeusIQ2m/06FGMaH0CHCldw8ylQboYVAibH6b8L1jTg8nigZF7zcOdWmYA2COgnspIv/8skomWIKFKzHFCf+iWLLNmTV22MnWxkRVwcYgIOIODCPnhXaQob70bK8lgTXixAsg43xwuQrCM0Xg2yB/TPq55h80FlmQEvwTKeH2XWlJkAadRoNfz9hyo+3uuPOXPO4J9/Oic3kcIUoQ9fi4vr6DA0EBM6XUODQmeAl7fQfN57mN5+G2p530vVtegYa3T8YyBal72G0Y0Oic8fvMyD4Rs7485zDzQqcRvrL1aG77gf4OakHML2MDg3qvwwBg++/AZ2NkrUI9fEMmhQ/Aa2vh+X+TkgSwBHW2BD6vNCZaLntpgAE5BCwPQiKxbWCHz6AvkLFoa1rYP6HCt9RFVKkSZZZElxlTkaCb4LrGuNuJe3ERdnBSurOFi5lgC6/wfkKWYyi2TP10zWETNWzALEjPBTNi17QHefcC8hsZIFiIYGggpTMhMgjRuvxt27Q0FrofZuO4+5c0+jeD4d9h0LAWIiXv+JFJGOXtXOYUKzfRi6rhualLyF3tXPCUvKTh8Pn48WoFCeV6ksG7SmB1wcIjGn86Y0rQ6JsEe5GeNx74tpCZ/PPtgAlx4VwKIe6xN+V/q7EbCzmY/iHrFY3BPpn4S+1gr3Iz1w4DNvuLrnBZ5fUZLaxR/aCjUWoP98PMoC9w8BDq6AowcQpQMCzwC5vZXE8mhd4p/YaBXEuSgTYAJMwNgETC+yVG9cYYioki2y6J2/uBRA0fiUl2tJ4P0bJouEyJ6vGXvEaaE+FiBa8MJrG2QP6GYfJy7D6d7MRezaYqlnSWjIjVJNSSZArK0B3TMl0fvFDfhduYk3euTBkS924q18xzB4eWuUzx+IsU3S32KajN96uRzmHa6PbYMX47h/UYz8522cGDU3Vb8mb2+NK4EF8He/ZbC2TjwH4WmoMzycdOJ3k7a3gY1VLP7XZmfC/XXnDMd3791D0/Y1geKtgXxVMWHSMXh52WDvf1OxdfsOtC1njT7VY+HmBATpgL8vAJsvAe3btsEfA/LCrVg1oNbYjFnf3gqcnQd03QY8PA7sGwm8eyL1PYcmK0Km498A7XYVf4U9BZw8lN8dmqTscFX/f4mfr6wLNPgOKJq49TYOTFB2xKo8CAjwAfaPB/qeTLyH8mZI8JAIop8xEYiJCMEV33MoX7oYbGIjk4sjEkox4crvSEQlFU5J/x5fJunvUpaH/mdVSB3E3BgTYAI5gIAJRJbuKeD3X/rs3jkAeCm7Mxr7kj1fM7b9WqiPBYgWvGAGAbL/bCim/PYsVe+/+sATjavn0hAVNiUVgfAXCXkYsc+u4cWdk/CwegqroJtARFBCcb/n7mjz22CxdOr0vSKoUCAQy3qvgrN9xgf+0ZdKwzd2wY5rZeFsF4kl76xNWH5V7cdPcG7MT7gX5ArvqZNRLv9jONjGiDaH1z+MwXVOYN35yvh8eztY2TmjUW0XzP+lNRwKlhNbv1JifP36qxAQMBTW1olrdIOCwvHuu1vh7/8KtjG3UcX5B1x6HIZXEUBuB6BJSeCjLm+iJEVhaFvbtssAu0wO66SO7BkO+O1QyrZeAhSspfD5qxrw3jng1T1gkTfgUQ6wcVA+qzYcqDIYuL4OOPi58g1akUZA8/kALc+g66UfsLo+MCQguWgJDwK2vask59u5AC0WAPkzPuAk260/J+4UkUoQWZGvRdRrIZWmiFL/WVy0DlHhYYgO1yE6UgcHq3A4WEfwC4MJMAEmoBCo9RnQ+HuT0GABknWsLECyztBoNcga0DThaTnifrp206noHAkxmlsNqyjyVfKEb5H8fV35XXhq4ZhWIyRAOiweBN/xJjotnJY2FW+jRDGKNgNcihgv3B0TAxwYDzw+BbiVBh4dBwYkntdjGFRt3pXtBIgJMYfoYnH6Sjj+OxaCY75KTpE+V+MazhjS2Q2F8trqU1y/MkJoxSAmKgIXzp1GlYrlYBOXRgQrQXCpF1lKxCuN++LPp9HPUi7FBHImgXy1gPeSRKCNSEHWfM2IJmuuKhYgGnKJrAH986pH2HwwMt2ed2poj9G9C2qITDY1JSoUeHEzYclUwu5SJDLEQXZZu4wiQCifggQG/SnSUFmSJPuiaMPGDixAZHPXQHu0lbnfwygcvajD75teqrJoQn9PNH/TGTZJIm2qKtCzcLYXkK+FljiAk4RPyiWBaS3/ixdd8csAMyqTnsii31MkjS8mYCgBew9ghH5f2KltQtZ8Ta1dllSeBYiGvCVrQLcZ5Y/IDFbh2NsBO2YX1RAZCzaF/gN+eSt1NIPyNEKUQyDNfuWvnigyCr4J2PESPJk+yfYTWJUwdRGxOHstHLtPhMHnTOLWzZlV07SWEuUo4GHEKEdmjQJg/+kBSeNFhA/PnkW1KpVgYxWnbM6RNP8qrRwsffKu4oVYhgJM/0iexjFqzzz7PMAIdV9c6NsJWfM1fe2xxHIsQDTkNVkDusUwf8RmkI9KXxjuns8CRO+hERMJBN1OjGSIJPDrwPNrQEj6S930rt8YBb2bAMUoktESyFsZsLE3Rq1chxEI8AQWuB8YhWO+Oize8hLhEfony08a6ImmNZ2T5RMZwSWqqmD/qcKlycI5yociP4s2w6DlgkmX96W1RDCDZYNqN8eg/ydlX+4VgEGXTNKqrPmaSYzXSKUsQDTiCDJD1oDuNekeAl+kH9rO726N1dO8NERGA6bQ8oNgv8RD+OIP5Ht2yeQi41moM5ovVM7yePQqt9hVKp+Lcn7LiZFzYG8bgzgbB1gl5GM0B9xKiqRvunx8/GFvb4N69YroDfLPP30xdeoxUX7y5Lro379SqnvHj/fBli23YW9vjZIl3bBkSRu4uTkiMjIGdJr6qVOPxcRw9uymaNJEEbRr1lzFtGnHEBMThw4dSmD69Mbi9wsWnMP8+edgY2MFFxd7LFrUEhUq5FV29fq7udJ26COlT075lH/TjlaZCSnaicraHihST+++49KfwLGpSvm6k4GK/VPfS7tb3dqitE+sKcHd0U0czIhdQ5XcFdo9q+lsgMQfXVfXAMenibwBlOgANJoufh179hdEHP8Rjk4usLJ3AVotUhLts/EVGRWHCzfDse9UGLYfVcayPleL2s744G035HOXG+XIyLYcNXnVx0kWWIZ9qCGnRbwEzi8EzvysnE+VlavpHKDGiKzUkO69suZrJjFeI5WyANGII8gMWQP6zJUQjJv7PN2ezxrhgRrlXTREJgNTaH3y/cMA7QDlVgooUt/wRGj6VogONXq9jW1CTkbguay/CI1B0zk/vj70HlwKl8a4L9sDLoURExuLc+fOoVq1auluHPD114fFpH7cuDf1suL5cx1q1VqOU6f6wsrKCjVrLsPp0/3g7u6Y7P6dO/3QrFlR2Npa47PP9ovPSFDMn38Wp049wpIlbREYGIq2bTfg5Mm+ePEiHNWr/yXqypfPGf37b8N771VE8+bFEBwcgTx5lB2mNm++iV9+OYcdO7ont/fI18rOUm+O06sfopDae3TPgRW1gHdJQFgBy2sCfU8Dju7J2/TbqSTe0175Bz5TPiNBcXa+Ij7aLAHCAoH1bZVteGnnsmXVlbqc8wHb+wMV3gOKNUdM2Aucu3Jb8eGdrcD5X4BuO/Tvo4WUDHwejeOXdPhrWzCevVR2TtPn+mKQJyiJPOmuafrcJ6sMT15lkTZdO+xD07HVu+aQh8CZ2cD5X4HIYL1vE9uj05c6KS8rW2B0BEDb05vgkjVfM4HpmqmSBYhmXCFPgFBiZ+fx9/AqLPVSh9zOVvhnppeYeGr+en0CqhAfZC+JERIhGZ2ASuFn2nqVIhgJQuO6ck6E7on5u+xe+vVSqdZA4Xqpkr7jxUTTpt4YM8YHISGRsLePwt9/d4eXl6s46ZyiCSQKKlTwxPffN0LduitgY2ONfPmcMHduczRsmHF0a9WqK/DxCcDCha0ED4pmNGnijd69y6fLZ+PGG1i37jpWrGiPYcN2o27dQujXr6Io37z5Wnz3XUPhogkTDmLPnp7i98uWXcLRow/wyy8tk9VL7f/11yVs356OAKFzN3zGAJEhgFNeoM1SwKUQcGYOcH6BIgoogtDwe4DO6oiPmjSbm/me8FdWAfd8gJYLFZsomuHVBCjfO/2xcWOjsmVv+xXA7mFA4bpAhX5KeYre0Fkh1PmDE4Aee5TfX14GPDgKtPgleQ7B9bXA5b+AbtvNPxazaAFFuS7diYDP6TD8sz9E79pa182F9zu5Iq+bdqIcGRnPk1e9XavZguxDM7qG/i8+ORPw/SNtIZGWabTjYp/jQO4iyheQa5siLsnOcFYU9X7HByj8lsk6xgIk62hZgGSdodFqkDmgHz2LRp8vkidBe+WzxsxRBaUncRoEkMTGb8WV8xZSXrmLAr0OAy9vJiZ/U07G/QPKN9HmvgrUSkz6pgRwWnaj50UCJFcuO2zceBObNnWGh4cDZsz4D9evW4mIQ+HCv+LOnQ/g4GALOluDlkSljICsWHEZM2em3pqwVCk3rFv3NmbNOonw8GhMnqy8vL/55iicnGwzjKB07LgB77xTDn37VsCiReexa9ddrFrVAQEBwahefRkWL24toiWVKy/FoUO94eWVG++8s0Us19qypatohyInP/54CpGRsdi7tydKl04RdRDRjFwATfjp5HWKJNCyJjqIqs0fwILCwOA7ylkddB4HLYlKGQG5skL5zy7lRcK10/+LiJOzlF1+aOkVXUe/AWydMo66bOwIlH0HqNAXuLAIuLsLaL8KeBWgRD1aLVaiJX9WBnodAnJ7Af++oyzX6rJFCJD7WyfCO3AdrOh3PfcCJEQt8HoeHIMTl3RYtTMYAY/1P2H+y8F50aiak2ajHCxALHAwqjCZBYgKWMYq+ugkcGI6cGO9/jXSboydNqQ+Ayo2FjHn5uPFlT1wL98cNtWGmSzyEW+szPma/oAsqyQLEA35S/aATnoS+uwx+VGppINlRD7IZwEHgLVK/oAmr6LNFZFRrCXgUT7xALssGktigqIbM2acRIkSrqK2kJAwvPGGJ3bt6ok2bdbBxcUOnTuXRufOpcTSK7VLsNQKEMrpoCVXGza8LcZPdHQsxo/fj337/FGsWB5ERcViyJAqwqYtW25h6tSjYqJZr15h3Lr1Ev/80zkZlZUrr+C//+7gzz/bJadFYoKiGydnALQ9MF0Ues9VCOi+E1jfRlmiVaqz8oeEndolWGoFyLFpypIr+k+Rohx0ujnlhwTsA/IUU7YtrTwEKN1ZyRmh3BLKDaHoVtAtoPM/KSIgrwVV2z+zOFLk3B4bG4dr/pE4eDYMq3e90rvRdvVyYVBHN3i4KnlKlnzx5NWSvafYzj6U5EP64vDuTkV40DtS36vmJ0DjWckPfU1xr2wfyp6v6YvKksqxANGQt2QP6KQCZO8vFrTrVdhTYH17IPCE+bxH38THn4/h3QxwfSMh6duURpGYoOVU27bdxtGj76b6jzMmJhYHDtwTE/3t2+/g4sUBYsKfNAckswiImiVYS5f6YuHC82JZlbOzXZpdr1dvJX7/vZWSVJ7kokjJzZtBmDEjuZCkSa27+1y8fDkytQCh9b53tgF9jqZui3J47h0Abm8B7mwH+l9UJvxJ80Yyi4CoWYLluxS4sFBZVpXeqewr6wGtf0+dVE6REjoDpvGM5D6kLejmuZts60hjjM1XYbE4dVmHtXte4dpd/Xe2mTIkLxpUdbKcLzn0hCV74qOnWVxMBQH2oQpYhhSlL2au/a18efTknP41tPwNqDJYr/KyfSh7vqYXBAsrxAJEQw6TPaAtSoDQNre3NgHn5ivfHMu4XAon5mN4NVS+aTdzbgwJEJroL1p0AcuWtUPt2gVw6tQZODkVRaVK+eHvH4zixV0RFRWDYsUW4fLlgVi8+CKCgyMxZUp9vahREjolnp85854oX6OGkjju4eGU7P4dO+5gzJh92L+/l0gqj7/CwqJAeUa5ctlj1y4/fPPNMRw40Et8TEnp+fPnEgnpTZuuwdq1HVGmjAdu3HiRsOSKxNOUKUdw6tTrPIr4iimaYesMXFwEtF2mrO+l3cloeZ1neSDYH3Atrvzut2LAgMuA72IgIhioP0WvvoOS0CnxvN8ZpfyyGkrieMoDGO/sUPJQ3tmvLAWLv6Lo3Io4ZamY3y7g2DdArwPKp5SU/v8bCYhlgGubAh3WAh5lEPP0Ks7dDVWS0P22AUenAH1P6WevhFLky9v3o3DovA5/btV/T/329XNhUCc3uOe2/ChHRphlT3wkuDzHNcE+NJHL6X3ouwQ4/QPw8o7+jfTYC1Cun4pLtg9lz9dUoLCYoixANOQq2QNa0wKEQrWBZ4Gb/yiTOFNdtDwqPpJRsLZ5TvpW0bf45VQtWhTDyJF78fJluFiC9emn9TBoUGU0bboWL19GCAFA+RgTJtTB9evP0b37ZrHsSZ8kdDLnjz8u4ttvjwvLJk2qg4EDK4u/Dx78Hz78sCpq1SqIUqV+R0REDDw9ld2x6tYtjAULWsLP7yVat14n2itSxEXkfxQrpiwX6937X5w/Hyj+/uWX9dCrVznx91Gj9mL37ruws7MWu23Nm9ccFSsmj5gkLKcq1gLYOxKIfKkseaoxGqg4QJnU0+9o7JTvC9SZADy/DmzproTu9UlCJ2Mu/gGc+FbxSp1JQKWByt//GwxU/RAoWAtYXAqIjgCcPJXPCtUFWi4A6NT29a2V9ihRsvViZSkWXf/2Bp6cV/7+1pdAOUWUxe4ZiYjrW+HonBtWTu5As3lAXiWB31yXLjwWZ66FY6PPK5y5FqG3Gd98mBf1Kme/KEdGAGRPfPR2BhfUmwD7UG9U+hWkL3Loy8KzcwDdU/3uoS3N+/sanP8m24ey52v6QbSsUixANOQvmQOaJqjNhwUk9H7PfG/zL42gb67vH1R2FKKt+Ix10eQwPh8jXxXAPrexajZ7PbJfumbvcDY0QAs+pPdBQGA0jpzXYdE/QXpT7tTQBQM7usLVJXtHOViA6D0kLLKgFp5BiwSX0ujgAOD0j8CFBcohh/pc+WsAPXan3upcn3uTlJHtQ5nzNZUoLKY4CxANuUrWgKYdsMb9/BAPniVuw1vY0wqzRhdCQU/JW1/SVqq0i9GV5Uq0w1hX2T5AhxXGqk2z9ch+6WoWhAUbZi4fRkTG4vyNCGw5FILD53V6EaQUlakf5UOdio7m/8JCL4tNX8hc/jN9z3JOC+zDLPr66SUlv4O2ENf3onOQWv0O2KSdO6hvNfHlZPtQ1nxNLQdLKs8CREPekjGg6ZvOt8cFICSN+YaLE7Dp/9j7DvCoiu79N41A6C10pIn4Q6RIkyqCoiCIUgQLVRCkl0+sICq9FxVEpYhIkyJFQOmdhN6kCYReQgkJCSmb/3/umpCQ3b1379577szuzPP4fJ/mzJlz3nfuZN5MG0ewEsL2wrMbgdgqx/V95jDQZgtQrJ45vjnySj3ocpS614RCySH748OuI7GYukj7ddQt6mdDh6a+vcrhqrNR8uc1nZ6zRCSHOgm5tB0IGw38u0q7A/YuUvXBhp+npOaQYr6mHVQxLaUA4Yg3ig69/8R9DJrqfPIxrnduVHnahC1K7MafM8uAPcOBh9oPsqajp2A1oMFk+357VthBY/bWwuOFvQPS9bzhAxxHXSU1FOpBl0cMRI/JTA4Tk5Jx9OxDrNkRjb/D2AF59RLgD4z4MD+qPi1XOdTRkle4asGIdxszv0Hec3c7PvaY79lVduFxZaf26s2X2a8jN6lQc0gxXzMJKm7cSgHCDRU0L6G3+SQCt1zM//PlBBaNVLmSlx3yZa+PshfI2QNuRWpnnOwzG/Y+wumlwN5R+lGu0heo8TkQ8tiBZOaRvYS++GXgHrsVi73cnvIS+nogh0DXCutHR95f7wF2vFQ1+hdn5L0k7D4ai4nzb8P2aJely3TfbJBdWeXIHuLPCyzCxGE0f8Ik7kWBSg41kMkeST0xHwgfC0Qe11DhP5N39wMFKmu312lJzaEUIDqJSlNNChDPMTTMA0WHbtQzwuWkhO3x/vtbF5N3Nulf0tguPtiVtExoMBHSap39mtqLW4B/5gPHZuvDJXMeoOF3QNmW9kfn1IoWMaTmQ+CfUw+6AkPFbeiecphkS8Y/5+Oxfk8MVm6L1pRnUCAw8sNQVH5KoMdHNWVGb+Qpf/QRyxYfR0By6KJPxN8HDs8E9k0Eoi9p6zxsF8Lbu4FshbTZG2BFzSHFfM0AWLh2IQUIR/RQdOgm/SIQ5+LtsMyZgDWTnAgQNtmfWQK4H5ERtYBgIEn7dZ3pHJRqBtQfC+R5iiM2xAiFetAVAxWxotTD4b3oJIQdj8O0xXcQFWPTlHDrhtnx3qs5kU2ucmjCS6uRHv60+pZ2NAhIDh3gzM5q7p8CHJymfdt0yVeBZkucP8xqIp3UHFLM10yEiwvXUoBwQYM9CIoOPWXhDSzf4vx6vBb1M6PPW6GOUbm4FViU/tVq3fCxg2hVetsfbJNFNwLUg67uQGVFpwho4ZBdHnHmUgI2hMVg0d/3NaEZHOSHkT3zo1JZ+zstspiDgBb+zGlZejUKAclhGiTZo7/h44BD0+3bmrWUqoOAeqPt7x9ZVKg5pJivWQQlWbNSgJBBrd4QRYdmH+lLvS87DeavqUUQEODkTv8V7YAzC9QTcWSR9//sD8EVa+ATh8P1geR+LepB1/0IZQ01BJxxGBNrw75/4jBj2V1cvZWo5kb5eZtG9lWOrFmsmwhoCtSLjOQ3KD6ZkkMA1w/YD5afXKid0Jd/Aip01m5voiU1hxTzNRPh4sK1FCBc0GAPgqpDf//7bSzekHGveOuG2dCjZR7niEzNDcRrf6QMFboCtb8CshbkCGXvCoV60PUu9PjIJoXDihUr4vLNZGzaF4O5a6I0BZc52A8jP8yPik/KVQ5NgJlgJL9BE0AldumzHLJt1Rc3AXtHAxfWa0e9zWagmEG7IbS36tKSmkOqJgkLeQAAIABJREFU+ZpB8HDpRgoQjmih6NAerYBMzAzYVM55NP4ZePpdwx4X4ogeLkOhHnS5BEHgoOLibdh34gGmL7mOy5HaHuRq93IOvPNKDoRklqscPFAvv0EeWPAsBp/j0JZkv6GSPR7IbqvUUtg5zw5HgdxltFiT21BzSDFfIweRuEEpQIgBd9UcRYeeOP8qVm5PcBpGszpB6P+2k5srvi8MPLjqPIWQQkCPKxwh6v2hUA+63o+o+RleuZWIrfsf4Ifl2lYTs2ZmZzlC8UzpYPODky24jYD8Bt2GjLsKPsNhYpz9tfKwsfabLLUU9v5Wy/VA5lxarC2zoeaQYr5mGZhEDUsBQgS0lmYoOnTjPhFIcLGdnF3PuW6Kk1uwwqcAW/o6T6X+ZKBqHy2pShuDEKAedA0K26fcJCQm48iZh5i39h4OntJ2Uxxb4Xi7cQ5kCZarHLx3FvkN8s6Qenxez2HcXfuh8v2TgAfX1QFhFuU7Ai/P1HYdvjaPplpRc0gxXzMVMA6cSwHCAQkpIVB06IY9I5SnO5wV9rTHBmfvgNhswET2NocjB/5A/wTAX06YKLsU9aBLmZvIbd28m4jtB2MxddEdTWlkyWTDqF4FUKFMFk320ogfBOQ3yA8XeiPxWg6jrwD7JgEHpmi/Jr/eGIDdasUmAwIVag4p5msCwa8rVClAdMFmTiWKDv3u0Mu4cjPJaQKF8wdg3rAizhOcHAIkxqb/uV8Q0HYLUPh5c4CRXp0iQD3oSiocI5CUlIwT5+Ox8K8o7Dj82PfhBLT2TXKg7cs5EBSQjIMHD6JSpUrOb6CTwHOLgPwGuaVGc2Bex+Htk/ZtVkd/0owBXl8BlGmu3Z4zS2oOKeZrnEFseDhSgBgOqX6HFB364cOHeLW/8yXYPycWQHCwi73m49P8VeSZLkBoFaBid7nyoZ92j2pSD7oeBetlle/eT8LOw7EY9+ttTZnlyOqP0b3y46kn0n9fkkNN8HFrJPnjlhrNgXkNh1f32G+0OrNMc+547yAQWlG7PaeW1BxSzNc4hdqwsKQAMQxKzx1Rdehx825hzc4HGQJ+rU4IBrydz3kibO/WhDRbrAZqfKTIc2ikBycIUA+6vkyEzZaM0xfjsXTTffy1N+P34wibjq/lxFuNsiM4k/OtiZJDsXuV5E9s/lj0QnPIfi+fX2sXHpe2aCMjxxPA27u96op8ag6p5mvaCBXTyisECHsl+ObNmwgNdfKCtyDcUHbotCsh2UOARcNVVj4YhreOAnMq2NEs2wZo5saDRYJwIFqY1IOuaPh4Gm/0Axv2HIvF8FmRmlzlyuaPUb1CUbZ4Jk32wk9+NGfpvYbyGxSfWyE5tCXaHw1kV+nePKyNhFJNgdcWA0Hed9aMmkPK+Zo2csWzElqAxMbGYuTIkVi2bJmyd5rto/77779x+vRp9OjRQzg2qDv0ix9GKBh1aZYD77yq4Yq9iZkA239X+Pa8w/21fMJ1AB0BUw+6OkIUqgr7Y8b5qwn4Y1s0VmzJ+Fino2Q6N8uJNo1yIFOQvkObkkOhukiGYCV/YvMn3B8BEmKAIz8D+8YDURe0gV9tMFB3pHAHy7UlZ7ei/g6p52vuYCGKrdACZNiwYbhw4QJ69uyJDz74AOHh4bh69Sq6du2KVatWicJBapyUHZp9rC/1vqy0/UypIEzsH6p+ADbt+Q+5/YqL/kU96HKRtMFBxD5kjwHGYcgPtzR5zp3DH6N7hqJMMe2rHK4cSw41wc6tkeSPW2o0ByYEh7GRwIFpwP6JwMN72nJrPAt4pqM2W8GtqDmknK8JTo3T8IUWIPXr18eKFSuQK1cuVK9eHXv37lUSrVatGsLCwoTjjKpDbzkQg2EzM24pGdo1L+pXzuoYN3n+g8v+RD3ocgmCjqAu3UjAmp0xWLA+SlPt91/PiVYv6l/lkAJEE8xCGslvUEja0gXNNYdslSN8gv0qXa3lra1A0bparb3CjppDqvmaV5DjJAmhBUjt2rWxZcsWBAYGpgqQuLg4NGrUCNu3bxeON4oOnXblwxFAf00t4ngl5MYh4JdK9ir/9x7w6lzh8PXGgKkHXVExjE9IxsFTcfjqp1t4EKd+eULenAHKjVWlihizyiEFiKg9Rz1u+Q2qY8S7BZcc3jxiP99xYp42+AIzAx2OAblKabP3MitqDinma15GUYZ0hBYg3bt3R61atdC+fftUATJ//nzs3LkT06ZNE447ig49cf5VrNz+3zkOBwg1qxOE/m8XyviTtNuvet0DgnMIh683Bkw96IqE4fXbiVi/JwazVmrbrtDtjVxo2SA7ggL1neXQi43kUC9yfNST/PHBgydRcMMh22lweZv9Rqtza7SlVKgG0HIdEJxTm72XWlFzSDFf81KqUtMSWoCcPXsW7777Lp544gkcPXoUVatWxYkTJ7BgwQKULFlSOO4oOnTjPhFISHQOTVAgsG5KcdcCRJ7/4KZvUQ+63CTuIJDEpGQcO/sQo+ZG4vpt549tplTNlysAY3qHokShIEvTkhxaCr/HjUv+PIbQcgeWc5hsA878AYSNBq7u1obHM52Bl2YA/oHa7L3cippDivmal1MGoQUII+fOnTvKOZDz588jX758aNmyJQoVcvAXfAGYpOjQDXtGgP2RxVnx8wM2fPuYAJHnP7jtPdSDLm9A3I5KwsbwGHy35K6m0Hq0zIU3XsiOwADaVQ5Xwfk6h5qI49hI8scxORpDs4zDxIfAiV+B8LHA7X+0RVt/HFB1oDZbH7Ki5pBivubt9AkvQLyJIIoO/e7Qy7hy0/lfhwvnD8C8YUXSw3p9PzDvOft/q9AVePkHb4Jd6FyoB12rwWKPAZ68EI8Jv93G2UvOtxKmxBmaOwBj+oSieAFrVzmkALG655jXvq99g+YhaZ1ncg4fRgGHfwDCxwEPrmtLvMVKoPRr2mx90IqaQ4r5mrfTKLQAGTFiBF5++WVl61VKYbdfbdiwAR9//LFw3FF06LQPEDoC6M+JDh4kTHv+o3cUkCm7cNh6a8DUg64VOEbFJGHrgVhMmH9bU/M9W+VCi/rZEcDRKocUIJqoE9LIF75BIYlxI2gyDmOuA/sn21c82EOCWkr7Q0D+Z7VY+rQNGYf/oUwxX/N2QoUWIHXq1MH69esREhKSylNMTAwaN24sb8Fy0XPHzbuFNTsfZLB4rU4IBrydL2NN+f4Ht+MA9aBLAQR7DPDfywn4dskdHDz1ULXJAnkCMLZPKIqG8rvKIQWIKo3CGnjjNygsGToDN53DO2fsqx2HZ2iLMGdJoN0uIGsBbfbSSj5EKGAfEFqAsJUP9vaHv79/KvQ2m015B2Tfvn3C0UGpqNOuhOTK5offvg5FcHBwRszk+Q+u+5HpvziJsn8QZ8OuI7EYPivj+zSOQujdJjea18uGAH9+znLohcpbONSbv+j1JH+iM2jiK9rX99lvtDq1WBtIpZsDry0E2JW6sriFAPV3SDlfcwsIgYyFFiAtWrRQtlrVrFkzFfLdu3eDbc36448/BKLBHip1h37xwwil3f7tcqFZXSfX6l4LB36tZg+wUi+g4VThcPXmgKkHXaOwZKscF28k4sfld7H9UKyq20J57Wc5iuQXc5XDVYKicqhKmo8YSP7EJ9pQDtkf7S78bb/RKmKDNnBqfArU/gZgt8DIogsBQznUEAH1fE1DSMKZCC1Ali1bhrFjx6Jbt24oUaKEchPWzJkzMWDAAOU2LNEKdYfWJEDSbr/qEw0EOXkpXTSwvSRe6kHXE9gextsQdjwOQ364pclN37a58Vod71jlkAJEE+VCGon0DQoJMEHQhnBoSwJOLbE/Hnhjv7ao2YO+7GFfWTxGwBAO3YiCer7mRmjCmAotQBjKixcvxty5c3Hp0iUUKVJEeZSwTZs2whCQNlDqDu22AJHvf3DXr6gHXXcBuBaZiDmr72Hd7hjVqoXyBWJ831AUzOtb99rzzqEqcT5uIPkTvwN4xGFCLHBstv1g+b1z2sBoux0oUlubrbTShIBHHGpqIb0R9XxNR4jcVxFegHCPsBsBUndoVQHCHkeaEPAoAylA3GCTxpR60FXLKiExGQdPxWHwtJtqpsrP+7fLjaa1s8HfC85yaErYgRFvHOrNw1frSf7EZ14Xh3F3gIPfAXtHAQnR6iAEZgE6HgPYAXNZDEdAF4ceREE9X/MgVG6reoUAuX37NtjtV2lLsWLFuAXdWWDUHTpFgAx6Jxea1HZwBuTqHmD+f+drnhsIvDBOOEy9PWDqQdcRnpH3krDgryj8vvG+KtxF8gdiXN9QFMjjW6scroDhgUNV4qSBUwQkf+J3Drc4vH8J2DcR2DdBW+KFngdargWCnZyz1OZFWqkg4BaHBqBJPV8zIGTuXAgtQNhNV//73/9w9epVsEOtfn5+qf974sQJ7sBWC4i6Q6cIkI/b58bLNR287ZHu/EcMEPToumO1XOTPaRCgHnRZVkm2ZBz/9yH6TrihKcmB7+TBq89n9elVDilANHUVIY2s+AaFBIrjoDVxGHkCCBsLHJulLRP2aG+j7wH/NLsItNWUVjoQ0MShDr/OqlDP1wwMnRtXQguQ5s2bKzdgsTMfWbJkSQcqOw8iWqHu0CkC5JMOufFSDRUBIrdfcdmdqAbde9FJWLb5PuauiVLFoWio/SxH/txylUMVLCbokpJw8OBBVKpUCQEBcrKiBTOebCR/PLGhLxaXHF7ZZb9K9+wKbc5fmAg810+brbQyDAHq75B6vmYYUBw5ElqAVK5cWXnvI+07IBxh63Yo1B3apQCR5z/c5s+KCmYNujZbMk5fjEeP0dc1pfW/9/LglZpZlVVIWdxDwCwO3YtCWutFQPKnFzl+6mXgkF2le26NXXhc3qYt0DdWA6WaaLOVVoYjQP0dUs/XDAeMA4dCC5B27dphzJgxEPG8hyPuqTu0SwFyeSew4L9bOtgd5XWGc9BdZQiPI2DkoBsda8OfO6Px/e93VYFmqxwT+oUiXy65yqEKloqBkRx6Gous7z4Ckj/3MeOtRiqHFcoj4Mxiu/CIPKYtzA5HgHzPaLOVVqYhQP0dUs/XTAPOQsdCC5B58+ZhyZIl6Ny5M/Lnz58Oxueff95CWPU1Td2hXQqQdOc/HgBB6be46ctQ1jIaAU8GXXZu6vzVBPQcex1xD5NVQ/u4fR68VEOucqgC5aaBJxy62ZQ0NwEByZ8JoBK7TIqLwpV1X6HolV/g90DD2bacpYC3dwEhocSRyuacIUD9HVLP17yReaEFSLly5RxywraByEPorrtrbGwsmg60X5WaJROwZGT+9Odo0goQef6D22/f3UE3Lt6Gv/c+wIT5t1VzKl6ArXIUQJ6c8lyCKlgeGLjLoQdNyaomICD5MwFUKpcPbgEHpgK7v9LWYpkWQNMFQGCwNntpRYYA9XcoBYjn1AotQDxPny8PVB36s2+vY9exhxmSf758MIb3LADI8x98dQwX0WgZdC/fTMDAyTdw43aSal6fdsyLhtVC5FkOVaSMM9DCoXGtSU9GIyD5MxpRAn/3zgPh44GD07Q1VvNzoNZXgDzjpg0vC6yov0Oq+ZoFUJI1KQUIGdTqDVF06LQrH44iWj0+P7JEhgML69l//PyXQK2h6sFLC0sQcDToxickY9vBBxg+K1I1JrbKMXFAAeTOLlc5VMEyyYD6F6dJafisW8mfQNTfOASEjQH+ma8t6CbzgKff0WYrrSxFgPo7pJivWQooQePCC5CzZ89i9+7diIyMVN4ASSl9+/YlgM/YJig6dIv/RSAq/ZuN6ZLIkRVYXvCJR/+tb5xcbjaWZkO9pQy6RZ6ogCEzI3HmYoKq/88750WD5+QqhypQRAbUvziJ0vKZZiR/nFPN5gWXttgPlp9fqy3YdjuBwuKdI9WWnHdaUX+HFPM172TqUVZCC5C1a9di0KBBKF26NJgQYf975swZVKlSBb/88otw3FF06IYfRsDVcWN2ieqG0mkEiDz/wWU/SkpKxt7jcfjse/s5HlfliYKBmDSgAHJmk6scalhZ8XPqX5xW5OjNbUr+OGXXlmR/u4MJj2t7VYO0+QcjueNxBOQupWorDfhDgPo7pJiv8YeysREJLUBef/11tG/fHi1btkS1atUQFhamCI/bt29DroA47ihN+0UgNt55J8qaKQkri6UZgKUAMfaL88DbnftJGP7zLew/mfH8zuNuh7yfD/UrZ5FnOTzAm6oq9S9Oqrx8pR3JH2dMJz4Ejv8C7BkORJ1XD65wLSS1WI2Dx8/Kx0DV0eLWgvo7lALE864gtABhKx1MdLDXg6tWrYrw8HDEx8ejUaNG2Lp1q+foEHug6NDzVl3Bz2sSnWb2yQt78dLF1vaf1xkB1PiEGAXZXAoC7DHAQ6cfKgfI1UrxgoGYLFc51GDi8ufUvzi5BEHgoCR/nJD38B5waAawayiQGKce1LMfAA2/BfwDIDlUh4t3C2oOKeZrvGPuaXxCC5DatWtjw4YNyJw5Mxo2bIj58+cjR44cqFWrFg4cOOApNuT1KTp0k34RiHO6AmLDxtIlH+XdJxYIykyOgy83GP3AhvHzb2PL/geqMHzZNR9qPxuMgwcPyr/cqaLFrwH1L05+kRAzMsmfxbxFXwX2TwbCRmsLpMEUoErvdLaSQ23Q8WxFzSHFfI1nvI2ITWgB0r17d7Rq1UpZ8fjiiy9w7tw5RYzExcWBPVIoWqHo0M7OgJQPDseEIq0R5GdLA1sA0HYLUOS/F9FFA1SAeNnFCf9ciEfPMddVoy1WIBDT/lcQ2UP8U22pB13VIKWB2whIDt2GjKsKkj+L6LhzGggbCxyZqS2AN/8ESr7i0FZyqA1Cnq2oOaSYr/GMtxGxCS1Abt68CZvNhgIFCiAqKgrjxo1DdHQ0+vXrh+LFixuBD6kPig79xkcXcS/68WPoNqwvVQoBSHZwzbk/0D8B8H806SUFxQsbi42zYfrSu1i5PVo1u6+65UOdSiFO7agHXdWApYHbCEgO3YaMqwqSP2I6roXZD5af/l1bwx2OAvnKu7SVHGqDkmcrag4p5ms8421EbEILECMA4MkHRYd29A5I8xw/o2++Yc7fWKo3Cagm3rXGPHF77ko8unxzTTWkoqGB+O6jgsiWZpXDVSXqQVc1AWngNgKSQ7ch46qC5I+ADnaV7oX1duFxcZN6g7lKA+12ASH51W0BeQZEE0p8G1F/hxTzNb4R9zw64QXI9evXcezYMcTEpH/colmzZp6jQ+yBqkM//hL6wuJVkC8w0rkAyRwK9FTfIkQMF9fNsccAZ6+6iwV/3VeN85vu+VDrWeerHFKAqEIotAH1L06hweIweMmfiaTYEoGTi4G9I4BbR9UbevJNoMl8t9+ukhyqQ8u7BTWHVPM13nH3JD6hBcjChQvx9ddfIzg4GFmyZEmHw/bt2z3BxZK6lB067UrI2hKlERSQCPYGiMPiFwQMcHF3ryVo8dfolVuJeHfIFdXAiuQPxPSPCyJrFs+3tVEPuqrJSQO3EZAcug0ZVxUkfybQkfAAODoL2DkEiLut3kDNIUCtL+H8r2iuXUgO1SHm3YKaQ8r5Gu/Y641PaAFSt25dfP7552jcuLHe/JWzI+wAO7u2N2vWrHj//ffRsWPHDP7YTUNTp07F0aP2v8JUrFgRn376KUqUKKH8+549e9ChQ4d0QuiDDz4AOyivtVB36Bc/jFBCW1e2AoKSopyHmSkn0Puu1jR8xi4xKRmL/o7CjyvuqeY8okd+1KyQXiSrVtJgQD3oaghJmriJgOTQTcA4M5f8GUhI7G3g4Ld24aGlsNWOp9tpsXRpIzn0GELLHVBzSD1fsxxgEwIQWoDUqFEDu3fv9uixNfaSOtu+NXbsWFy+fFkRH6NGjUL9+vXTwb1lyxbFjoketuIyefJkbNy4EX/++WeqABkwYAB27NihmybKDp2QkIDGfa8qsQ4uPgwvB/3sfAWkYh+g0WTdeXlTxVt3E9HmU/VVjkJ5AzDzs0IIyez5Kocr/KgHXW/ikpdcJIe8MKEvDsmfPtzS1Yq6COybAOyfpM0ZO99RuKY2Ww1WkkMNIHFuQs0h5XyNc+h1hye0APnyyy9Rp04d5RpePYV1oOrVq2Pp0qUoW7as4mLixInKdb5Tpkxx6TIyMlJ5b4QJoNy5cysrIKIIkF/+vItZK9OueCRhQ6lSzlev+yUCAQF6IBa+TpItGau2R2PygjuquYzqmR/Vyxu/yiEFiCr0QhtQ/+IUGiwOg5f8eUDKrWNA2Bjg+Fx1J4EhQKcTQA7jb7iUHKrDz7sFNYdSgHjeI4QWIGxF4q233kLBggWRP3/62y5Gjhypis7x48fRunVr5RB7SmErGkx8pKxsOHPCfj58+HCknDVhAqRz587IlSsXMmXKpKyUMEHC/l1rSenQTAyFhOg7lKzWFlv5aDIg48vadULW4KtCPZTqaS/ptTVdAjzZQs2tV/08KiYJLT9Wv7EqNHcAfvwsFFmCzV3lUBMgR44cQYUKFRDgoyJR9M7HfnFKDsVlUfKng7vL2+EfPhZ+51arVk4uXAu2FquBTNlVbfUaSA71IsdPPWoO2Xzt1KlTePrpp02br/GDrjmRCC1APvnkE2UbVLVq1TIcQmdbqtRKeHg4evbsqaxepBS2hYr5ZWdCnJWLFy8qwoedP2nSpIlixt4kuXv3LkqXLg12M9fQoUPh7++P6dOnq4WR+vMUAaK5gg7DaSty4MKtTICDDVcbSz+heIzOXh7ROSrhcqk+PrHywW54PPhvMOZvzqGK6PuN7+KpogmqdtJAIiARkAhIBNIgkGxDzshtKHhhFrLdV7/R6kbh1rj45CDAzzdX32XfEQMBKUD08yS0AKlcuTJWrlyJokWL6kKArYC0adMm9WA5c7J27VrlfIezFZCrV6/i3XffVf7p1KmTS5Hy8ssvY//+/RnEkbNKFCsgr/S9jKS0j52nBpOMjaXtB+qT+ibovk1EFxEWVIqJtaHFR/YzMK5K3pz+mD2kADJnsm6Vw1V81H/1UcNL/tx9BCSH7mPGUw3JnwobSfHwO7kAfru+hN99+8UnroqtwRQkV/xQzczQn0sODYXTEmfUHMoVEM9pFlqAvPjii4pgYFue9JSUMyDLli3Dk08+qbhwdQbk2rVraN++PVq1aoVu3bq5bJIdaG/YsKEiQLRup6LYU/hSrwiHAiRXwC0sLfGcPaeBj7+UrgddvuokJydj99E4fPb9TdXAxvYJxXPlMqva8WBAve+Vh5y9LQbJodiMSv6c8BcfDRyZCWwe+NjGXif2LdcBJV62pDNIDi2B3dBGqTmkmK8ZChCHzoQWIPPnz8elS5cwcOBA3fvfWV32JsaYMWNw5coVZVVjxIgRGW7BYtuq3nvvPTRv3hy9evXKQCU7jM5WYooUKYJbt25hyJAhiI+Px08//aSZdooO3XvMRRw7n1FgvJ5jDvrm/+/qQy8RIA/jbWj18WXExLkWVLmz+2P+14URzOkqh6sORD3oau7M0lAzApJDzVBxaSj5e4yWBzeA/VOAPcO18dXxOJD3aW22JllJDk0CltAtNYcU8zVC+CxpSmgBwq7KZZP9wMBA5SaqtGXz5s2aAGXvgLCzHNu2bcvwDgjb4jVz5kxUrVoV06ZNU94BeXw1Y/Xq1ShcuDBmzZqF2bNnK+dAcuTIoRxCZ1f85smTR1MczIiiQzfpG4E4B0cYUs5/KMEKLEAOnY5D/4kZD9k/TsK4vqGo8pQYqxxSgGj+hIQ0pP7FKSRIHAct+fuPnLv/AuHjgUPfqbOVqzTQbjcQkk/dlsBCckgAsslNUHNIMV8zGTLL3QstQNjWKWfljTfesBxcdwOg6NApjw8+HpuoAiQhMRkdhl3Btcgkl3BnD/HH4pFFkCnI6Xvv7tLFhT31oMtF0l4WhORQbEJ9nr/rB+xX6Z5coE5k2VZAk1+BAH3bptUb0Gfh8xzqg42rWtQcUszXuALYhGCEFSDsOll2DS67sYo9DOgNhaJDe4MAORURj+6j1K/J9ZZVDld9m3rQ9YbvjLccJIe8MeJePD7JH7s68OImYM9IIOJvdcCe/xJ4fgi3l5v4JIfqrAllQc0hxXxNKAJ0BCusAGG5spfQ016hqyN/rqpQdGgRBQh7DLDnmOtgwsNVyRzsh+VjinrdKocUIFx9poYHQ/2L0/AEfNyhT/FnSwLOLAN2DgUij6sz33QBUO4tdTuLLXyKQ4uxNqt5ag4p5mtmYcWLX6EFCHvoj221YuctvKFQdOhmAyIQE5cRLd62YEVcS0DHr9SvyR3fNxSVveAsh97+Sz3o6o1T1nOOgORQ7N7hE/wlxtlfK2c3WiVEqxP29h6gUHV1O04sfIJDTrA2KwxqDinma2ZhxYtfoQXI119/jaVLl4Jdx8tuoGIP/6WUvn378oKx5jgoOnTTfhGIdbCQkCJAbMl+8B/k8KEQzXnoMbTZkjF42k3s+8eBOkrjkFH856RiCAr0rrMcejBjdagHXb1xynpSgHhrH/DqbzDuLnBoOrD9E3X6ArMAnU4COYqp23Jm4dUccoa1WeFQc0gxXzMLK178Ci1A2LW4joqfnx/mzp3LC8aa46Do0GpbsIZe+x7DxnbXHLMnhtdvJ6Ld51dUXfjCWQ5VEJwYUA+6euOU9aQA8dY+4JXfYPQVYN8kIHysOm1F6gBv/glkyqZuy6mFV3LIKdZmhUXNIcV8zSysePErtADhBUSj4qDo0GoCpMm/x7HmW3PuZGePAX79UyQ273+gCtn6qcUQGCBXOdSAoh501eKRP3cfAcmh+5jxVMOr+Lt9EggbCxzV8H5VpZ7Ai1MAv0c7D3jixZ1YvIpDdxL3IltqDinma15Ej8NUvEKAREZGKo8Isvc48ubNKyxnFB3akQAJRDzWl7a/BP/i2QvY+F1xwzC8HZWkPAaoVsb1CUUVQV4fV8uF8ufUgy5lbr7SluRQbKa9gr+re4BW0fMwAAAgAElEQVS9o4Azy9XJaPgdUKmHup1AFl7BoUB4mxEqNYcU8zUzcOLJp9ACJDo6GoMHD8aGDRsUTNnWK3YeZNSoUciePTtPOGuKhaJDOxIgZTIdxQ/FmhomQCbMv41V29UPKv41tRgC5CqHpr7hzIh60PUoWFnZIQKSQ7E7hrD8sat0z68Fdg0DmABRKy3XAyVeUrMS8ufCcigk2uYETc0hxXzNHKT48Sq0ABkyZAjOnTuHTz/9FMWLF0dERARGjhyJEiVK4KuvvuIHZY2RUHRoRwJkeMHOeD6rXcTpWQGJfmBD80GXVLMc2ycUz8lVDlWc3DGgHnTdiU3aakNAcqgNJ16thOPPlgicXAhsGQTEqL+nhI4ngLzleIXfkLiE49CQrL3LCTWHFPM172IoYzZCC5AXXngBS5YsQb58+VIzu3nzJlq1aoUtW7YIxx1Fh36lTwTiE9NDk/YK3lcuXMDaKepbsGYsvYOFf99XxViucqhC5JEB9aDrUbCyskMEJIdidwxh+Et4ABz5CdjURx3wXKUBdpVuFnG3NKsn+chCGA7dScrHbKk5pJiveTuFQguQmjVrYuvWrciUKVMqTw8fPkS9evWEfKCQokNPWXADy7emv+o2rQCZUuQ6+rQNzdDvYx/a0LS/+irH6F75Ue3/snj7d8NNftSDLjeJe1EgkkOxyeSev9hI4MA0YNeX6kCXbQ00mQcEPPqdql5JfAvuORQfYtMzoOaQYr5mOmgWNyC0AOnatSuefPJJDBo0SHkDxGazYcKECThx4gR++knDLR4Wg/948xQdmn2kL/VOfyg8rQBJ6peIgIAAJbR5f97DzyvvqaL017RiCPCXN1apAmWCAfWga0IKPu9Scih2F+CWv6gLQPgE4MAUdYBrfw3U+IwdpFS39UILbjn0QqzNSomaQ4r5mllY8eJXaAFy9uxZdOrUCYmJicoNWFevXlUmz7NmzULp0qV5wVhzHBQdmrXx2qBb6WJKK0DYGRC1MrxHfjxfQa5yqOFE8XPqQZciJ19rQ3IoNuPc8XfzCBA2GjjxqzqwzRYDZVup23m5BXccejneZqRHzSHFfM0MnHjyKbQAYUDGxMRg8+bNivgoVKgQ6tevj2zZxHwQiaJDN+kXgbjHXkLXIkDkKgdPn+2jWKgHXT5REDsqyaHkz2ME2I1Wl7cBu74CIuwXirgs74QBBauqWfnMz+U3KD7V1BxSzNfEZ8V1BsIJEHbwnAkOVr788kvlH28pFB3a0S1YKQJkW3RjDL3+gx3brvlQr3KIt0DrtXlQD7peC6SFiUkOLQTfgKYt5S/ZBpz5A9g2GLhzynU2gVmAzqeA7EUNyNq7XFjKoXdBaVk21BxSzNcsA5OoYeEESLVq1bBr1y4EBgaiSpUq2L9/PxFU5jdD0aFdCZAuF9fhXHw5Qx8iNB81326BetD1bbTNyV5yaA6uVF4t4S8pHjg+D/irG5Cc5DrVInWBln8CQVmpIBGuHUs4FA4lvgOm5pBivsY34p5HJ5wAYQfP2cvnZcqUwZo1a9CkSROHKIwZM8ZzdIg9UHToJn0jEJeQklgynskchilFWiv/odHZs8gUFIg1k9Wv4SWGRjbnBAHqQVcSYTwCkkPjMaX0SMrfwyjg8A/A1v+pp1i5D9BgIuDnr27r4xakHPo41malT80hxXzNLKx48SucAImKisKCBQuURweXL1+OZs2aOcSSPUgoWqHo0CmH0AsEXsLoQu+hcFAEAv3sD4NExJdC/s6rkaWAdz86JVq/cBUv9aDrTdjxkovkkBcm9MVBwl/MdWD/ZGCvht9rL80Anu2mLxkfrUXCoY9iS5U2NYcU8zUq7KxqRzgBkhaoDh06YM6cOVZhZ3i7VB160ISr6B9fE4UDI9LdupgMwC9naaDLaZ+9jtFwUk12SD3ompyOT7qXHIpNu6n83TkDhI8DDs9QB6n1BqD4i+p20iIDAqZyKPEmQYCaQ6r5Ggl4FjUirABJSEjAG2+8gd9//x3BwcEWwWdss2Qd+tI2JC+sB6c3vr+1FSha19jkpDdTEKAedE1JwsedSg7F7gCm8Hd9H7D7G+DMcnVwOp0E8pRVt5MWThEwhUOJNykC1BySzddIUaRtTFgBwmCqW7cuNmzYkO4ldFr4jG2NrENv+gjYP9Z58FX+BzQQ7wyNsWyI4Y160BUDFbGilByKxdfj0RrGH7tKl12hu+1jgAkQVyVXaeDtvUCWPGKDx0n0hnHIST6+GAY1h2TzNS8mU2gBMnnyZEV89OjRwysoIuvQ04sDMRedY5a1GNA9wisw9fYkqAddb8fTivwkh1agblybHvNnSwJO/w783R2Iu+M6sLKtgSbzgIBMxiUgPcFjDiWGliNAzSHZfM1yZM0LQGgB8vbbb+Pw4cPImzev8hK6v/+j2z5+/VXDK7Dm4arLM1mHnhQCJMU6jzEgC9Dvga4cZCVaBKgHXdrsfKM1yaHYPOvmLyEWODYb2PChOgC1vwFqfCrP5qkjpctCN4e6WpOVzECAmkOy+ZoZYHHiU2gBMm3aNKcw9urVixOItYdB1qGn5QMeRjoPLDgv0OuW9sClpWUIUA+6liXqxQ1LDsUm123+2CrHoe+B7Z+pJ978d+DJN9XtpIVHCLjNoUetycpmIEDNIdl8zQywOPEptADhBEPDwiDr0CvfAU7Ndx532beBZuKtIBlGhECOqAddgaARJlTJoTBUOQxUM3/3LwP7JgL7xqsn/G44UOA5dTtpYQgCmjk0pDXpxAwEqDkkm6+ZARYnPoUXINHR0di0aROuXbsG9kjhrVu3kJycjPz583MCsfYwyDr0xa3AovrOA2uzBShWT3vg0tIyBKgHXcsS9eKGJYdik6vKX+QJ+/sdx39xnWhAMNDlLJC9iNiACBi9KocC5uRrIVNzSDZf82IihRYgJ06cQJcuXZA9e3bcuHEDBw4cwNatW7FkyRJMmTJFONrIOjS7bWVmCSTfj8h4FW+OJ4D3z8m9xoL0HupBVxBYhApTcigUXRmCdcrflV32bVYXN7lOsEhdoOVaIChEbCAEjl5+gwKT91/o1BySzdfEp8ZpBkILkPfeew+vvPIK3nnnHVSrVg1hYWFgKyJNmjRRhIhohbRDR11A8pKX4Xfn1COYcpcDWq0DchQXDTqfjZd60PVZoE1MXHJoIrgErtPxxy5CObcG2NATiLrguvXKfYAGEwG/R5enEIQrm3CAgPwGxe8W1BySztfEp8dhBkILkBo1amDXrl3K7VfVq1fH3r17lSSfe+457Nunco86h4RSd+ikxEQETA56hMQAm1z54LBfuAqJetAVDB4hwpUcCkGT0yAV/vaHo3KWf+C/rqN6Mi/9ADzbVd1OWpAhIL9BMqhNa4iaQ+r5mmnAWehYaAHy0ksvgV23GxoamipArly5gk6dOmHdunUWwqqvaeoOzT7YgEmBj4IdmKwvcFnLMgSoB13LEvXihiWHApObEAPboR/gv2WAehJtNgHFXlC3kxbkCMhvkBxywxuk5pB6vmY4YBw4FFqATJgwAUePHsWQIUPQpk0b/PHHHxg2bBjKly8PeQ2veu+SAkQdI94tqAdd3vEQMT7JoYCsPbgFHJgK7P5KPfjOp4DcT6rbSQvLEJDfoGXQG9YwNYdSgHhOndACJD4+XhEfy5cvV5Dw8/NDo0aNMH78eOWFdNEKdYeWAkS0HpIxXupBV3zE+MtAcsgfJ04junceCBtjf8fDVclZCmBX6WbOLVByvhuq/AbF556aQ+r5mvgMZcxAaAGSks7du3cRERGBfPnyKS+ii1qoO7QUIKL2lEdxUw+64iPGXwaSQ/44yRDRjUPAzqHA2RWugy3bBmgyDwhIc7ZOgPR8PUT5DYrfA6g5pJ6vic+QlwgQ9ubHRx99hGPHjqFixYoYNWqUcg5E9ELdoaUAEb3HANSDrviI8ZeB5JA/TpSI2HXll7YAm/oCNw+7DPJyyZ4o2HwSAgLTnKnjNC0ZVkYE5Dcofq+g5pB6viY+Q14iQPr374/r16+jadOmWLVqFYoXL47Ro0cLzw91h04nQPJVADq4/iUrPMBemAD1oOuFEFqekuTQcgrSB5BsA84sB1a3A5LiXQfXfCmSSjXHwYMHUalSJQQEBHCWjAxHCwLyG9SCEt821BxSz9f4Rl9fdEJuwapTp45y7oNtuWKrIW+//TY2btyoDwGOalF36HQCpHwH4JXZHKEhQ9GCAPWgqyUmaeMeApJD9/AyzTrxof218r80XJH77n6gQGUlFMmfaYyQOZYckkFtWkPUHFLP10wDzkLHQgqQKlWqYP/+/amwsfdA9uzZYyGMxjRN3aHTCZAGU4AqvY1JRHohQ4B60CVLzIcakhxaTPbDe8ChGcC2wa4DCcgEvH8OyJb+nKHkz2L+DGhecmgAiBa7oOaQer5mMbymNC+sAGEPDSYnJyv/PP/884oAYf8/pbDHCUUr1B06nQB5axtQtI5okPl8vNSDrs8DbgIAkkMTQNXiMvoqsG8iED7WtXWRukDLtUBQiEM7yZ8WsPm2kRzyzY+W6Kg5pJ6vacFANBshBUi5cuWUK3dTChMeaf+d/fcTJ06IxgWoO3Q6AdLjBhCSXzjMfD1g6kHX1/E2I3/JoRmouvB55zSwaxhw4lfXDVfpC7wwAfBz/ccsyR8xfyY0Jzk0AVRil9QcUs/XiOEkaU5IAbJ3715VcKpXr65qw5sBdYdOJ0D6JwD+8gYX3vqEWjzUg65aPPLn7iMgOXQfM101roUBmwcAl7e7rv7yj0CFLpqbkPxphopbQ8kht9RoDoyaQ+r5mmYgBDIUUoAIhK9boVJ36HQCZOCj7WtuBS2NLUWAetC1NFkvbVxyaCKxbFvuhfXAmveA2JuuG2qzGShW3+1gJH9uQ8ZdBckhd5S4HRA1h9TzNbcBEaCCFCAckUTdoaUA4Yh8naFQD7o6w5TVXCAgOTShe9gSgZOLgTVvqzvvfBrIXUbdzomF5E83dNxUlBxyQ4XuQKg5pJ6v6QaG44pSgHBEDnWHlgKEI/J1hkI96OoMU1aTAoSmDyQ8AI7OAjb2ct1ezpLAu/uAzLk9jkt+gx5DaLkDyaHlFHgcADWH1PM1jwHi0IEUIByRQt2hpQDhiHydoVAPujrDlNWkADG3D8TeBg5OA3YOdd1O2TZAk3lAQJBh8chv0DAoLXMkObQMesMapuaQer5mGFAcOZIChCMyqDu0FCAcka8zFOpBV2eYspoUIOb0gaiLwN4RwKHprv3XGQlUHwykuT3RqIDkN2gUktb5kRxah71RLVNzSD1fMwonnvwIJ0BsNpsm/OQ7IOowSQGijhHvFtSDLu94iBif5FAHa7eOAVs/As6tcV359eVAmdd1NKC9iuRPO1a8WkoOeWVGe1zUHEoBop0bZ5bCCZDH3wBxlph8B0S9c0gBoo4R7xbUgy7veIgYn+TQDdYubQfWdwbYWx6uynsHgNBKbjjWbyr5048dLzUlh7wwoT8Oag6lANHPVUpN4QSIljdAWHLyHRD1ziEFiDpGvFtQD7q84yFifJJDFdaSbcDZVcAKlZUM9o5R1wggWyHSbiD5I4XblMYkh6bASuqUmkMpQDynVzgB4nnK/Hqg7tBSgPDbF7RGRj3oao1L2mlHQHLoBKukeOCf34C1HV2DWaQO0HI9EJRFO+gGWkr+DATTIleSQ4uAN7BZag6p52sGQsWNK+EFyP79+7Fs2TLcvHkT06dPx7FjxxAbG4uqVatyA7LWQKg7tBQgWpnh14560OUXCXEjkxw+xl18NHD4B2DLQNekVukLvDAB8PO3lHzJn6XwG9K45NAQGC11Qs0h9XzNUnBNalxoAbJmzRp88cUXaNKkCVavXg0mRo4cOYKxY8di7ty5JkFmnlvqDi0FiHlcUnmmHnSp8vKldiSH/7H94AYQPh4IG+Oa/pd/Aip05qaLSP64oUJ3IJJD3dBxU5GaQ+r5GjdAGxiI0AKkWbNmGDZsGKpUqYJq1aohLCwM8fHxqF+/Pnbt2mUgTDSuqDu0FCA0vJrZCvWga2Yuvurb5zm8+y+w/TPg5ALXXeCtrUDRutx1E5/njztG3A9Icug+ZrzVoOaQer7GG95GxCO0AGHbrMLDwxUc2KFzdkA9OTkZNWrUUP6/aIW6Q0sBIloPyRgv9aArPmL8ZeCzHF4/APzVDbhuH8Odli5ngFyl+SPuv4h8lj9uGXE/MMmh+5jxVoOaQ+r5Gm94GxGP0ALk9ddfx8iRI/F///d/qQKEbcEaOnQoli5dagQ+pD6oO3RSYiICJv/3IvDAZNJcZWPGIEA96BoTtfSSFgGf4jA5Gbi4CVjaBEh66Lwj5HgCeO8gkDkX953Fp/jjng19AUoO9eHGUy1qDqnnazxhbVQsQguQ5cuXY8qUKejevTtGjRqlCI/vvvsOvXv3xmuvvWYURmR+qDt0UnwsAqaG2POTAoSMZyMboh50jYxd+rIj4BMc2pKAM8uAla1d0162DdD0V4BdqStI8Qn+BOFCb5iSQ73I8VOPmkPq+Ro/SBsXidAChMHAVjpmz56NCxcuIH/+/Gjfvr3yj4iFukMnRV1CwMxiUoCI2Fn+i5l60BUYKm5D92oOE+OA43OBvz5wjX/dUUC1jwA/P255chaYV/MnHBv6ApYc6sONp1rUHFLP13jC2qhYhBcgRgHBgx/qDp10/m8E/P6SFCA8kK8zBupBV2eYspoLBLySw4f3gP1TgJ1DXHPf4g+gdDOh+4dX8ic0I+4HLzl0HzPealBzSD1f4w1vI+KRAsQIFA3yQd2hbXvHwn/bR1KAGMSfFW6oB10rcvT2Nr2Kw+grwM4vgSMzXdPW/hCQ/1mvoNar+PMKRtxPQnLoPma81aDmkHq+xhveRsQjnAApV64c/DQs0584ccIIfEh9UHdo26q34X/yNylASFk2tjHqQdfY6KU3hoBXcHj7JLDhQyBio3NS2YOBH1wGshb0KuK9gj+vYsT9ZCSH7mPGWw1qDqnna7zhbUQ8wgmQtO97nDx5EvPnz1fOfBQtWhSXLl3CL7/8gnbt2qFjx45G4EPqg7pDJ//8FPzunJIChJRlYxujHnSNjV56E16AXN0DLH8deHDdOZmFawOt/gKCsngl4fIbFJ9WyaHk0F0EqOdr7sYngr1wAiQtqK1atcLo0aNRuvSjO+LPnj2LwYMHY8mSJSLgny5G8g49Ps2BT3kLlnD9RfjJq5CIGx+0cJMfdpXu+bX2q3RdlSp9gRcmAGzlw4uLcPx5MRd6U5Mc6kWOn3rUHJLP1/iB2rBIhBYg7AX0PXv2ICjov7csAOUl9Jo1a2L//v2GgUTliLxDSwFCRa1p7VAPuqYl4sOOheHQlgicXAisedc1W41/Bp7p5DOMCsOfzzDifqKSQ/cx460GNYfk8zXeADcgHqEFSJs2bVC3bl306tVLORfCXkH/9ttvsXnzZrkCoqVzSAGiBSWubagHXa7BEDQ47jlMeAAcngFsHuAa4be2AUXrCMqC/rC5509/aj5TU3IoPtXUHEoB4nmfEVqAHD16FF27dkVAQAAKFiyIa9euKQc6f/jhB1SoUMFzdIg9kHfotALkxW+Bit0Bf+/eLkFMqenNUQ+6pifkgw1wy2FsJLBnJLBvvGtW3v8XyFnSB5mzp8wtfz7LiPuJSw7dx4y3GtQcks/XeAPcgHiEFiAs/+joaGzcuBHXr19XREiDBg2QLVs2A6Chd0HaoS/vABY8/tfKAKDtFqBIbfrkZYu6EKAedHUFKSu5RIA7DqMuAJv6218ud1ayFwc6HAaCc/o8u9zx5/OMuA+A5NB9zHirQc0h6XyNN7ANikd4AWIQDly4IevQNhswkZ2bsTnI2x/onyBXQrjoEepBUA+66hFJC3cR4IbDm0eAla2BOyedp1C2DdD0V8A/0N00vdaeG/68FmHzE5Mcmo+x2S1Qc0g2XzMbOAv9Cy9A5s6diwULFuDq1asoVKgQ3nrrLeVaXi1vhViIu8OmyTp0+BRgS1/n6defDFTtwxs8Mh4HCFAPupIE4xGwlEN2o9XlbcDC+q4TqzsKqPYRoOENJuMR4tujpfzxDY0w0UkOhaHKaaDUHJLN18SnxmkGQguQn3/+GbNmzUKXLl3wxBNP4MKFC2D/rUOHDsp/01KioqLwxRdfYOvWrciaNSvef/99h2+IHDx4EFOnTgU7d8JKxYoV8emnn6JEiRKpzcybNw8zZsxQtoXVqVMH33zzDXLm1L5FgaxDTy8OxFx0Dk/WYkD3CC3wSRuLEaAedC1O1yubt4TDZBtw5g/gjzdcY9piJVD6Na/E3aikLOHPqOClHwUByaH4HYGaQ7L5mvjUeKcAeeWVVzB+/HiUL18+NcHjx4+jf//+WLdunSbaBg0ahJiYGIwdOxaXL19WxMeoUaNQv376vwhu2bJFsWO3bgUHB2Py5MnK2ZM///xTaWfHjh0YMGCAIoCYGPrss8+UVZhJkyZpioMZkXXoicGALd55XP6ZgP4PNcctDa1DgHrQtS5T722ZlMOkeODYHOCvbq4BbX8YyC/eRR5W9BJS/qxI0AfalByKTzI1h2TzNfGp8U4BUrVqVezduxf+aW5ustlsqF69OsLDw1VpYx2I2S5duhRly5ZV7CdOnIhz585hypQpLutHRkaiVq1a2L17N3Lnzo2BAwciNDRUeQSRlfPnz6Np06bKz7Nnz64aC6kAmZwdSIx2HlNgNqDvfU0xSyNrEaAedK3N1jtbJ+Ew/j4QPgHY9aVrELtfA7IW8E6gTcqKhD+TYpdu7QhIDsXvCdQcSgHieZ8RegvWm2++iQ8++ACNGzdOReKvv/7C999/r4gKtcJWS1q3bo1jx46lmrIVDSY+UlY2nPlgPx8+fDi2b9+umDRv3lzZ9vX666+nVqlUqRLmzJmjbNfSUlI6NBNDISEhWqros1n5FvzP/o4076Cn+klmR9NLtwSaLdTnW9YiRYANukeOHFGunWbXUcsiHgKmchhzHX47PoX/8TlOgUkuVBO2ln8DgZnFA4+DiE3lj4P8fCEEyaH4LFNzyOZrp06dwtNPP23ufE18apxmILQAYQ8OskcI2dW7xYoVw8WLF5VHCJmAYP9NrbBVkp49eyqvqacUtpXqk08+Uc6EOCusHXbY/fPPP0eTJk0Us0aNGinbrtK2y7ZrjR49Wlkp0VJSBIgWW49sEhNRZXvNDAKEiQ9W9tfZDQTKW248wlhWlghYiEDwg4soeeJzZL3/6I8rj4dzvUg7XCrTH/CTb/9YSJVsWiIgERAYASlA9JMntABhaR8+fFh59Zw9QsjeAWnVqhWeffZZTYiwFRD2mnrKwXJWae3atcr5DmcrIOy2rXfffVf5p1OnTqntsBUQdoCd/W9KqVy5MmbPns3fCggL8NAMBGzqmRpriviwvfg98GxXTfhJI+sRoP6rj/UZe18EhnJ4Yz/8F9SCny3RKVC2l35CcvkO3gekRRkZyp9FOfh6s5JD8XsANYdyBcTzPiO8APEEgpQzIMuWLcOTTz6puHJ1BoSJHHbFLxM53bqlP8T5+BkQdiMXWx3h8gxICmhpX0Jn9/u/+qtc+fCkQ1lQl3rfqwUpen2THnPIrtKN2AAseck1Vm23y0dGTehNHvNnQkzSpXsISA7dw4tHa2oO5RkQz3uBkAKEbYFSK2xLlpbChENsbCzGjBmDK1euKKsaI0aMyHALFntp/b333lNWONi2r8dLyi1Y7FpgdgsW256VnJzM5y1YjgTIwJQ1EC2oSRteEKAedHnJ25vi0M2hLQk4uQhY87ZrON4/B+R8dF24N2HHQy66+eMheBmDgoDkUPyOQM2hFCCe9xkhBQjbc5dS2CQ/pbBrb9m/s/89ceKEJnTYOyBMLGzbti3DOyBsC9XMmTPBbtuaNm2a8g7I44fDV69ejcKFCyttsXdApk+frlzXW7t2beWQOpfvgEgBoqlviGBEPeiKgIloMbrNYUIscHg6sHmA81SzFQU6HgOCc4gGh3Dxus2fcBl6f8CSQ/E5puZQChDP+4yQAoQJg6JFi6Jt27aoV69eumt4UyApUqSI5+gQeyDv0Gm3YMkVEGK2jWmOetA1JmrpJS0CmjmMuwPs/BI44OKK8LKtgKa/Af7yEgmqXqaZP6qAZDtuIyA5dBsy7ipQc0g+X+MOcc8DElKAsBWGFStW4LffflNWG9iNVC1btkS+fPk8R8RCD+QdWgoQC9k2pmnqQdeYqKUXtwTI/cvAus7AhfXOgas7Cqj2EeDn6HJtibeZCMhv0Ex0aXxLDmlwNrMVag7J52tmgmeRbyEFSFqs2FW6CxYsUF4lZ6+Xf/zxxyhQQMyHtMg7tBQgFn12xjVLPegaF7n0lIKAUw4jTwAL6wOxN52D1WIlUPo1CaaFCMhv0ELwDWpacmgQkBa6oeaQfL5mIbZmNS28AGHAREdH48cff8SMGTPADoHXrFnTLLxM9UveoaUAMZVPCufUgy5FTr7WRgYOr+wCflN5O6jDESDfM74GFZf58vINsvOPKf9wCRTHQTEO5YOuHBOkITSjOWRniVP+cdQ8+XxNAwaimQgtQNgbIGwbFnv9vEaNGsqZEPb4n6iFvENLASJqV0mNm5fJj/BAWpiAwuGBA6iU8xIC/njDdSQ9rgMhoRZGK5t+HAGrv0GbzYYbN27g7t27igCRxX0EGG4JCQkICgpSJp2yiIeAGRyyvpArVy6EhoZmOGtMPl8TjxLViIUUIIsWLVK2XUVGRqJ169bKP6Juu0rLEHmHlgJE9QPh3cDqyQ/v+HAfX1ICbEd+hv+G7s5DLVQTaLMJCMzMfTq+GKDV3+C5c+eUyRH7Hcgm0LK4jwCbvLLr+LNkySIFiPvwcVHDDA6ZKGVPMDCRX7JkyXR5ks/XuEDZ2CCEFCDlypVTHg5s0KABAgICHCLSt29fY5Ei8EbeoaUAIWDV3ELF/E4AACAASURBVCasnvyYm50Xe0+IAfaOAnZ/4zzJyn2ABhMBP38vBkL81Kz8BtnE6OTJk8rvw8BAefOZ3t7EJq/s9y+7Zl+ugOhF0dp6ZnGYmJiI06dP46mnnkq3CkI+X7MWXlNaF1KAsAcBXRU2gMydO9cUwMx0St6hpQAxk04S31ZOfkgS9LZGHtwCNnwInFrsPLPGs4BnOnpb5l6bj5XfIGv71KlTKFu2rNM/xnkt8AYmZtbk1cAQpSsVBMzi0Nk3Rj5f88IeIKQA8UIelJTIO7QUIMJ3JSsnP8KDR5nAvfPA7y8Dd047bTWpzVYEFBP3DBslnDy1ZeU3KAWIMT3BrMmrMdFJL1oQMItDKUC0oK/PRgoQfbiZUksKEFNg9WqnVk5+vBpYo5K7cQj4pZJLb0mdz+Lgv3dQqVIl+Vdso3An9GPlNygFiDFEmzV5NSY66UULAmZxKAWIFvT12UgBog83U2pJAWIKrF7t1MrJj1cD60ly7Caii5uBxS8695KtMNDxBBCcA5JDT8C2vq6V/HkqQNik7ejZh7h8MxFF8gfimdLB3J6BmD59Os6fP49Ro0bpJn3Pnj0YMGAAduzYkc6HWZNX3YF6ULF///4oVaoUevfureqFbWdv0qQJ2rVrp2rLu4FZHEoBYh7zUoCYh63bnkkFiM0GTExzgL9/EuAvD7u6TZrFFayc/FicOn/NJ9uAfxYCa952HtuTLYHXFgD+jw4MSw75o9KdiKzkzxMBci0yEYOn3cDVW4kIDPRDYmIyCuULxOheoSiY15gD7ZUrV06FMi4uTjkon3JY/oMPPkD37i5uf3OHBI22ZgmQqVOn4t9//8XEiRM1RmKemRQgxl4kIAWIeX1VChDzsHXbM5kAubwDWFAPgC1NjP5A261Akdpuxy0rWIeAlZMf67LmrOXEh8D+ycC2wc4DqzMSqD4YcPDGgOSQMz7dDMdK/vQKEPbX4o5fXcWVm4lISvNrIMAfykrIrCGFDF8JadOmjfJW15tvvpkBYXbTEMUtXlKApIderoCof+xSgKhjpNdCChC9yJlQj0SAPL7ykS4PP6B/olwJMYFbs1xaOfkxKydh/D6MArYMBI786DzkFiuB0q+5TElyKAzjDgO1kj9Hk6O5a+5h074HLkGNi7fhemSSU5sCeQOQOZPzFfEGz4WgfZOcbhGXVoBcunQJDRs2xMiRIzFt2jRkzpwZa9asUf593bp1uHfvHkqUKIFPPvkE1atXV9pJu8qQUn/06NGYMmUK7t+/jzfeeAOffvqpy5hSBEinTp3w008/Ke+mdOnSBe3bt1cugWHvgMyZM0d5Z+zOnTuoWLEivv76axQsWFB55HHMmDFYsWIF2GoOe5xu2LBhSj223Yn9PFOmTMidOzc2btzoNA6WB7u5LFu2bEqu+fLlw7hx45QVlMmTJyv+evTogY4d7TfhxcfHKysrq1evBhNq9erVw2effYbs2bMrP9+9e7cS45UrV/DSSy8p9dmVsSlbsLZt26bUj4iIQLFixZS6VatWVepKAaLehaUAUcdIr4UUIHqRM6EeiQDZOwnY1t959HUnAtX7mZCddGkGAlZOfszIRwif0VeBFS2Aa3udh9vhKJCvvKZ0JIeaYOLWyEr+9AqQ6Ac23I5KgqOH09kiXZ4cAcgWYr4AYecP2OSZrX4wEfLHH3+gTp06yJEjB+bNm4cZM2Yok3kmDBwJkBYtWmDo0KHKo8RMgDCb559/3mlfYQKEiQ8mhpi4Ye87dO7cWZn4P/vss/j999+xfPlyRdSwhx2ZONq7dy9+++03sIn8559/jsWLFyvi4+LFi0o7bFLvzhYsZsvyYm2+8MILivhYu3Yt6tatqwgo9q7LO++8g7/++guFChVSfDMM2BkY9k7J//73P2TNmhXjx4/H3bt38eKLL2LIkCF47bXXsHLlSkVgMAHDBMg///yDDh064Ntvv0WVKlWwZcsWfPzxx/jzzz+RJ08eKUA0jCpSgGgASaeJFCA6gTOjGokAmZobiL/rPPxMuYDed8xIT/o0AQErJz8mpMO3S3aF7s9lXcfY4wYQkt+tPCSHbsHFnbGV/OndgnXkTBwGTr6BRAeLIEEBwLi+oahQJrOhWDtaAWGrHqVLl3baTrVq1TBr1iw888wzDgXIhg0bULRoUaU+O0/C7NmKhrPCBAibkIeFhaWuIIwYMQIxMTHK5L9169YYPHgw6tevr7hgKw7sdrr169crKwj9+vXDhAkTlHbSvjrvrgBhouaXX35R2jh27JiyLW379u3In98+djRu3FgRCuyxZbaqwWJq1KiR8rOzZ8+iWbNmOHjwoLJqxN48W7p0aWrKTJSx1SUmQL788ktFrDDRklJY/kysMTu5AqLexaUAUcdIr4UUIHqRM6EeiQCZEAgkO196h18AMCDRhOykSzMQsHLyY0Y+XPq8uheYX8N5aAWrA29tAQL1Tdgkh1yyrjkoK/nTK0CsPgOSsoXq0KFDyspHSmHbopYsWYIbN24oZ1Cio6Mxc+ZMZXXA0QrI4cOHERwcrFTXcviaCZBevXopAiSlzJ49W5n8sxWJWrVqKVdh+6e5kIVtgWI2bAXh119/xbJly5TbuFhMTCSwlRJ3BUjaA+tMULCVILbykVJef/11dOvWDU2bNlVWZhYtWoRy5copP3748KHy3zZt2oRVq1aBYchWOFLKhx9+iKeffloRIF27dlVWcNKKJSaqmA3zLwWI+mcuBYg6RnotpADRi5wJ9UgEyLS8wMPbzqMPzgP0ijQhO+nSDASsnPyYkQ83PtnelH9XA8ubOQ+pcm+gwSTAz7Pb4ySH3LCuKxAr+dMrQFiiKbdgXbuViIBAPySl3ILVOxQF8hhzC1ZaQB2tgKQVEOHh4ejZs6dyBoO97M5EAFtpYFuN2LkHowQIWwFgbbEzGKywcydM6LAVkJYtWyrnOmrUcPEHh/8vdtgZlS+++EIRP2PHjlW2ajEhoeUWrMfFipoAcbQCwrZbMeHhaAWErW6wbVlMgLCtWWxVxdmVvFKAqH/yUoCoY6TXQgoQvciZUI9EgIRPAbb0dR59/clA1T4mZCddmoGAlZMfM/Kx3Kct0X6o/O8ezkNp/DPwTCfDQpUcGgalJY6s5M8TAcLAonwHRE2AsPMJ7FwGW2FgB7nZagg7i8HOSxgpQNgZEHYbF1u9OHPmjHImhAkHduCcrb6w8xHssHnx4sUVocHeDGErFEwsMbzLl7ef7WLnQZhIYu+SsDMiLG52eD3t6omjDumuAGEYbN68OfUMyEcffaSciWHCjB2UZ2Ljq6++wquvvqqsiLAzIGw7GhMdbHsXu+540qRJygoOW81hW7fYAX92sF4KEPUhQwoQdYz0WkgBohc5E+qRCJCkJGCSi79u9UsEAtK8D2JCntKlcQhYOfkxLgsOPCXEAju+APaNdx5M2x1AkVqGBys5NBxSUodW8uepAKEESk2AsFzYpJ7dDMUOW7OVCjaxZ+cYjBQg7CHClFuw2OF3dmaEtcV+/7LtYPPnz1f+uX79unIYnm3LYudEdu3apayWsMPn7LYrNqFnE3+2wsCEAFu9YYfac+bMib///tsptO4KECYa2LmTlFuw2NYvhhOLjZWdO3fim2++wdWrVx3egsUEFBMxbNsXy5dt32IH9wsXLiwFiIYPQAoQDSDpNJECRCdwZlQjESDnNwG/u3ihueVGoEQDM9KTPk1AwMrJjwnp0LuMvQ2sbgtc+Mt5210vADmKmxab5NA0aEkcW8mfSAKEhAydjZj1irbOcGQ1HQiYxaEUIDrI0FhFChCNQFGYkQiQmWWBqNPO08nxJND1FEW6sg0DELBy8mNA+Na5iLoIzCkPxN93HEOW/MD7Z4FM9rv2zSySQzPRNd+3lfxJAWIMv2ZNXo2JTnrRgoBZHEoBogV9fTZSgOjDzZRaJAJkQiYgOcF5/H5BwIB4U/KTTo1HwMrJj/HZEHi8dQyY84zzhsq8ATRbBPgbfwjXWaOSQwLeTWzCSv6kAMlILDt4zd7DeLw899xz+PFHx4+GGj15ZbdXsYcBHy/sPAY7nyGL8QgYzWFKhFKAGM9VikcpQMzD1m3PJAJkUhYgKc55bAGZgX6xbscuK1iDgJWTH2sy1tnqxc3AIhdbC+uMAKp/DLBX2IiL5JAYcIObs5I/KUCMIdOsyasx0UkvWhAwi0MpQLSgr89GChB9uJlSi0SALG8FnP3defylWwItlpiSn3RqPAJWTn6Mz8Zgj8k24J8FwJp3nDtusRIo/ZrBDbvnTnLoHl68WVvJnxQgxvQGsyavxkQnvWhBwCwOpQDRgr4+GylA9OFmSi0SAZKYCEwOch5/3wQgkG77iSlA+pBTKyc/3MKcFA/sHQ3sHOI8xI7HgLz/x0UKkkMuaNAdhJX8SQGim7Z0Fc2avBoTnfSiBQGzOJQCRAv6+mykANGHmym1SAQIi/zgdGDDo3cOkgEoG08azQAqdjMlN+nUHASsnPyYk5EHXuOjgXVdgFOLnDvpcRMIyedBI8ZXlRwajymlRyv5kwLEGKbNmrwaE530ogUBsziUAkQL+vpspADRh5sptcgECIs+zUqILe8z8H/3gFz5MIVVc51aOfkxNzM3vD+4AcyrBtyPcFwptArQbicQGOyGUzpTySEd1ma0ZCV/UoAYw6hZk1djopNetCBgFodSgGhBX5+NFCD6cDOlFqkAYRmMtx+4tdUeDv+an5qSk3RqLgJWTn7MzUyD97v/Aj+Vdm5YqRfw4mTAz1+DM+tMfJpD62A3rGUr+ZMCxBgazZq8GhOd9KIFAbM4lAJEC/r6bKQA0YebKbWkADEFVq92auXkxzJgr4UDv1Zz3nzjn4FnOlkWnrsN+ySH7oLEsb2V/HksQJKTgcs7gLtngFxlgCK1LbkJzgx6H39x3FUbZk1ezcgrrc/atWsrr6TXqFFDtamnnnoKa9asQenSLv5oo+qFXwOzOJQCxDzOpQAxD1u3PUsB4jZkPl/ByskPKfhsonR2JbDidefNtt0BFKlFGpYRjfkMh0aAxaEPK/nzSIBEXQCWNAbunQMCMgHs8oacJYFW64AcTxiCdOXKlVP9xMXFITAwUPmHFT1vYrgjKtyxdWfy+t5776FJkyZo166dIRh54kQKkEfoucOhO5hLAeIOWu7ZSgHiHl6mWlsmQOqMhH+Nj03NTTo3BwErJz/mZPSYV1sScHAasKmf8+a6XgByFCcJx4xGvJ5DM0DjyKeV/OkWIEzQz3oauHsWSE58hKZfIJC7DNDxuOErIW3atEHbtm3x5ptv6mbPHVHhjq07k1cpQHTTZ2pFdzh0JxApQNxByz1bKUDcw8tUa8sESL2x8K82yNTcpHNzELBy8mNORv95TYwDNvYBjsx03ExwTqDbRSBTdlPDoHDutRxSgMdBG1by53BytOsr4ORC18gkPADYCgjYHYiPFz/7CkhQiHMfT70FPO/immsHNdMKEDZZnD17NhYsWIA7d+6gYsWK+Prrr1GwYEGwn40ZMwYrVqwAWzUJDQ3FsGHDwH4/9u7dW/l5pkyZkDt3bmzcuNFpjEyAnDx5EpkzZ8aGDRtQuHBhDB06FNWrV1fqREdHK+1s3rwZDMcXX3wRn332mWLPYvr0008RHh6u2JYoUQIzZszATz/9hJ9//jl1JadRo0YYO3as0xiYWKlSpQr27duHo0ePokKFCpg0aRJmzpyJpUuXInv27BgxYkTqFqqbN2/iq6++wt69e5E1a1YwzLp16wZ/f/s5NoYZiyExMRFdu3ZV/n/KFixXmLK6cguWvsFCChB9uGmpJQWIFpSIbCwTIPUnwL9qf6IsZTNGImDl5MfIPFJ9PbwHLKgL3Dri2H3p14HmSwB/73mrxus4NKVj8OvUSv50C5C4u0DMVecCJGshIHMu0wTIL7/8gmXLlmHKlCkoUKAApk2bpky6f/vtN2zbtg2ff/45Fi9erIiPixcvKnEUK1YM7qxqMNvvv/9emeC/9tprWLlyJYYPH66IkZw5cypihk3yWVsMxz59+ihCaMCAARg/fjxOnz6NiRMnKmLnxIkTKFmypGLvzgoIs718+bIiOJgA6tKlC65fv44ePXrgjTfeUATFokWLsG7dOiVHZl+0aFEMGTIEN27cwPvvv68IDSZEduzYocTGBFCZMmWUXJYsWYJZs2YpAsYVplKA6B8/pADRj51aTSlA1BAi/DmpALHZgIkBSna2Mi3h32wR8N9fWQhTlk15iICVkx8PQ09fPfoKMKOIc5d1hgPVPzF8W4ihOeh05jUc6sxf9GpW8qd7C9al7cDihoAtPiP87DxIqw1A0TqGUpN2BYSdoRg8eDDq16+vtMH+ol+pUiWsX78eERER6Nevn/KX/WrVqiEo6NHDue4KkE2bNikrDSmlRYsW6NSpE+rUqYN69ephz549yJYtm7Kqsn37dnz55ZeKQGHCaOfOncq/lytXLh0O7goQtuLCxA4rTCwsXLgQa9euVf6diREWx/79+3H//n288MILSkxMILHCBBkTTvPnz8cnn3yi/PePP7Zvl46KilJWc+bMmaMIEFeYMvEjV0D0dWcpQPThpqWWFCBaUCKyIRMg7NaTBWzgT1Iysz9EGAC03WK/BUUWYRCwcvJjCEiRx4HZ5Z27avEHULqZIU3x6kR4DnkFliguK/nTLUAsPgPCVhoCAgJStxYxquLj45UVAbZl6ddff1VWSM6fP4+6desqk262UuKuAPnnn3/w7bffpvaEDz/8UBE6NWvWVFYVmPhIKUyE2Gw2HDhwADExMUq9v//+W9n61bx5c/Tv318RQ+4KkLQH1pmgYDdRsdWKFBHBRBYTP1euXFFWPMLCwlJj2rp1q7L9jIkitnrCxEqHDh1Sf161alUlTiZA1DCVAkTfgCAFiD7ctNSSAkQLSkQ2JAJEWflg21cc7f31B/onyJUQIr6NaMbKyY9H8UdsAha/6NxFx2NA3v/zqAlRKgvLoSgAmxynlfzpFiDK7Pe/W7CizgH+meyrITlLAS3ZLVjGX+qQdgXklVdeUSbWatfH3rt3D1988QWCg4OVsxZsq9bZs2eVrVFqhYmVx1dA2Lanjh07olatWmjQoIGy8sC2WLk6wMxWZNg2KCYOWrdubZoAYVyyFRC2FS1HjhxKeu6sgKhhKgWIWo9x/HMpQPThpqWWFCBaUCKyIREgYZOBrS5uFKo3CajWlyhj2YynCFg5+XE7dvZX12OzgXWdnVftcRMIyee2a5ErCMWhyECbFLuV/HkkQBgehO+ApBUgbNvQn3/+qRwCL168OJjQYGcc2GrB4cOHlTMZ5cvbV0bZGQ12CHvUqFHKhJytjLDD6ykHs53RmnIGZPTo0Xj11VexatUq5aA7W03IlSsXevbsiXz58innKthh8HPnzinnNdgqAxMu7OD5E088gbt37yqig61AsBu82EoIOyzPtpCplcdXS1ytgOTPnx/vvvuu0iYTXSlnQFi7b731lnI2ZtCgQcoqUalSpTBy5Ejl/EjKGRBXmLI4pQBRY0sKEH0I6a8lBYh+7AyvSSJAvi0AxN1wHnvmUKDndcNzkw7NQcDKyY/mjGyJwNaPgX3jHVfJVwF4JwwIDNbs0psMheDQmwA3OBcr+fNYgBiMhSt3aQUI2+rEtlmxsw3sHAT7iz9blWAHxnft2qVMrtnhc7Y6wbZksZuh2ASd3U7FhAM7IM7OQ7AtUq4ESNpbsAoVKqRM7J9//nmlCrsFi91IxQQJE0DswDt724NtcWKTfLZNKjIyUtmmlXK+gm0bY1u02Jaw27dvKzdnMYHjrLgrQBgWLFd2+1ZISAhatWqlHFhPEVvsADoTHI5uwXKFqRQg+ju6XAHRj51aTSlA1BAi/DmJAJmQCUhOcJ6VXxAwwMHBREIcZFPaEbBy8qMaJbvq8/fGwOXtjk0rfgg0nAr42a+Y9NXCNYe+SoobeVvJn0gCxA1IyU3NekOCPBEfbtAsDqUAMa9TSQFiHrZueyYRIFNzAfH3nMeWKSfQ+67bscsK1iBg5eTHacaxkcB3LrZRvfwjUKGLNYBx2CqXHHKIE68hWcmfFCDG9AqzJq/GRCe9aEHALA6lANGCvj4bKUD04WZKLRIB8ndf4NAU5/FX7AM0mmxKftKp8QhYOfnJkM29c8CPpZwn2Xa7vGXNATpccWh8F/V6j1by5+sCpGnTpsrtUY+XDz74AN27d9fc9/ROXlnbLAZHhb39wW6pkoUGAb0cqkUnBYgaQvp/LgWIfuwMr0kiQJKSgEkuHnHrlwgE2N8HkYV/BKyc/KSicy0M+NX+urDD0jUCyFGMfzAtipALDi3K3RuatZI/XxcgRvUfsyavRsUn/agjYBaHUoCoY6/XQgoQvciZUI9EgLC4T/0OrGyVmoH9HRAAzZcCT75hQmbSpVkIWDb5YbfnnF6arh+ly5G9VN7zDpDp0T37ZmEgul/LOBQdOE7it5I/KUCM6QRmTV6NiU560YKAWRxKAaIFfX02UoDow82UWmQChEWfZiXEVqwR/FuulSsfprBqrlPyyU+yDdgzAtjxhePESr4KsMcDmQCRRRMC5BxqikoaaUXASv6kANHKkms7syavxkQnvWhBwCwOpQDRgr4+GylA9OFmSi1SAcIyGK+se8DWaAb8K3YzJSfp1FwEyCY/SfHAH28C/652nFDtr4EanwF+9j4li3YEyDjUHpK0dAMBK/mTAsQNolyYmjV5NSY66UULAmZxKAWIFvT12UgBog83U2pJAWIKrF7t1PTJT/x9YFpuIDnJMY6vrwDKNPdqjM1OznQOzU7Ax/1byZ8UIMZ0PrMmr8ZEJ71oQcAsDqUA0YK+PhspQPThZkotKUBMgdWrnZo2+Ym+Cswo7By7DkeBfPaXimXxDAHTOPQsLFlbIwJW8icFiEaSVMy0Tl7Pnz+vvJzOXk1///33lUcRfaVMmzYN7CFG9jCikYXh2LhxY7Ru3VrVLXukkb3ozl6Gf7xo5VC1kccMpABxFzHt9lKAaMfKdEspQEyH2OsaMHzyc+sYMOcZ5zh9eAvIktfrcLQyIcM5tDIZH2zbSv48ESBswhYWFoalS5ci8uZN5M2fH2+++SaqV3dxo50Oftlr4AcPHkRgoP1cWMmSJZWXxD1th8W9YMECLFq0yGlUU6dOxfTp05UX1f38/PDEE0+gX79+qF+/fro6Wievn3/+OYKCgjB06FAdSIhb5e7du8p1w+vWrcM///yDrl27pibD5i1ZsmRR8GXFzOuHr127hjZt2iiv1zMe0hatHLrLghQg7iKm3V4KEO1YmW4pBYjpEHtdA4ZNfs6vt79a7qjkKgOwFY/AYK/Dj4eEDOOQh2R8MAYr+dMrQI4ePYr32rbFqVOn8ExyMrInJuJ+YCCO+vmhbNmymLdwIcqXN2aFkwmQJk2aoF27drDZbFiyZAnGjRuHnTt3pooSPd1GqwD5999/MXHiRDCsfvnlF0yaNAlbtmxBzpw5U5tVm7wmJiYqsXbs2FH5az3Lxd2S4sPdejzYz549G6zPMN7SlocPH+LZZ59VBEHRokUzhGpGzoyDtm3b4pVXXpEChIfO4UEMUoB4AJ7RVaUAMRpR7/fn0eSHXaV7YCqwqa9joCp0BV6aDvj5ez+QFmboEYcWxi2btiNgJX96BAibSNauWRPVYmPxgs2GzGmIjAOw2d8fYVmyYOeePYaIkLQChDUVGxuLSpUqYdOmTShc2L7Nc8WKFcpfztlfuJkAGjZsGJ588knlZz/99BPmzp2LqKgo5MmTR9kCVaZMGbRq1Qpsgps5sz2DHTt2ICQkJF23ZCsgKQKE/YD9jq1cubIigthWHrataNWqVYiJiUHNmjWVdnPlyoVLly6hYcOGGDlypGLD2sidOzf279+vCBH2z5w5c1CqVCmMGDFCETRse9Krr76KgQMHKisujnywlRMWf7f/1953gElVNF2fXRZWMhLEAAqoJAmSDCBKjhIVUFRQokRBUJSkRJUoIJIUkSThFVQkKwICviQlS3gFRBBBEAkCC7vw/6f57jo7zOzce+eGubvVz7MP7E53dfWpOzN1uqq627VT67169Sp69uyJIkWKoE+fPmpMlSpV1Lyc48KFC+r17du3q7USN+qo4UZsy5Qpg61btyqCwDWNHDkSd911l8KBa6esHTt2qN8Zxejfv7/6//fff6+I2ZEjR5A3b141f7CLE+n0P/nkkwpz3+ZPQDRSWK5cORVZq169Ol5//fWQa9AIqjb+0UcfVdEt4vjaa6+hfv1/6wwnTJgApsK99957SXQJRSLNfl5JBMQscqHHCQEJjZFjPVwjINU/RnSJVo6tUyayDgFTzs+1BGBpC2Dv7MCKVJ8MlPg3xG6dtiIpEAKmbChQRgwCbtrPKAGhk1a6eHHk+Pln1Lp2LSiGy6KjcbpIEfy0a1fYOPsSEOo7d+5c5XyvXLlSOdmrVq3C4MGDQceSxGL+/Pnq9aVLlyqHvFGjRli4cKFy9k+ePKmICPsZjYDQ2Z85cyZISkgYxo8fjwMHDmDYsGGKuPTr1w/Ehw68Rh7oGA8aNEjpSRLiT6aYSvbnn38qR57OeMeOHVG+fHl07949oAwSiZdeeknVj3Tu3Bnr1q3DK6+8ggoVKoD1DZyHTj7H16tXT62VkSKmjDF6xBQw+gmTJk1SdqE+1JW/58uXTzn61JNrIqmi/owWcE42khSSDKZRtWzZUmFQunRphQfXQsxJ8vwbCYHWNxQBoY4kWZRPe1+5ciXkGnwJCO1AXRhl4rPRq1cvRZYyZbpxp9SKFSsUKfzqq6+EgIT97nRXgBAQd/FPMrtrBKTmp4gu1iKCkBBV9CJgyPm5egmYnAe4/Fdg8c+sA+6qoHdq6WcRAoZsaNGcIsY6BNy0n1ECsmnTJlR+7DH0vnoVySVUMhIyNG1arFm/HtzNDqfRSeYOPHP26aSz0lsX+wAAIABJREFUMWpAB5uN9QSVK1dG8+bNE6epVq2a2rnPnTu32nkfPny4csK1aAc76iUgdM45js49nXQ6/Ix2MBLy+eef495771XEg9GCBg0aqGjD8ePHVQRkyZIl6nWt+ZOpkiVLqmhK4cKFVZe1a9cqIkOHXiMxvjI2btyI1q1b46effkqsYSABIMlhdIKN4+ls0/H2byQOdMw5no36sJamS5cu6nc652PHjlVRncWLFytSx//7t7fffhsZM2ZU0QWtkTCQ7DVs2PCm/kzHIwlkdMq3BYqAjBo1ShEGrSZEzxp8CQjJBYmH1ogPo03FixdXf2Kki9Ga1atXJxEtEZBw3qXujBUC4g7uAWd1jYDUmoHoB56PICREFb0I6HJ+Lv0FfJhM4XjbX4Esd+udUvpZjIAuG1o8p4izDgE37WeUgHBn+buRI9EkPj4kAPNjYlC5Rw+8++67Ifsm18HXaaeTuG/fPhUBIMGoWLGi2qWnw88UJq0xWjFkyBBFPujAz549G3v27FG793TMSQr0EhDfFCxN/unTp1WkInPmzIlzUjfu1n/zzTcqNYoEhGTEl/T4ruXUqVMqcrFly5ZEOb/88osiVrt378axY8dukkECwugAnWitUQad9ocfflj9ieumY8+ICNPViBMd+rNnz6rXGdnYuXOnSk/yj8iQAJFc0IFnFIlE5cMPP7zJPCR9JKO+hdxM8WIEh+lh/s1IBOSzzz5TUSytGVlDIJv64yMRkLDejhE1WAhIBJlDCEgEGcMjqiTr/Px9EPj43927m5bU5TyQ7kZYW5p7CLjpwLq36pQzs5v2M0pA2rZujUNTpyJp+W5gWywFcG/r1pj80UdhGcvfSaawrl27IleuXGq3nxEB7v7zBK7kGh1ZpkeRiJCQcEeezm6oU7ACERCmMzECwigBi6f9d8+16AUjN7Gx/8aKzERAfGUYJSBMeyJZYeH8bbfdplKnGKXRZCZHQLg2ngC2aNGim2BlHQjx1yInoQxstAbE1yZG1qCHgEgNSChreed1ISARZCshIBFkDI+oEtD5ObYBmBMkleruKsBTy4HoG0diSnMfATcdWPdX730N3LSfUQLidgSE1v7f//6n6gPo/LI+gREHnq5EJ5tF1PwepKPO1CLWfLAwnYXWTKEaM2aM2tXnaVZauhNrSRgNCNT8i9B9+zANjFEKFobTGSfpYHSG6V96CAhlMRrz119/KWLE6AnvBWEkg1GOQDKMEhDWcpB0MIrBqAgJG4/C1UNAWMDOoniSlBYtbqRYazUgjNC0b99eYc4UJ+rOo5KZonb77bffBCVPwSLxoz6+LVgRui8BMbIGPQSE9Sy8M4SRM98mKVje+ywVAhJBNhMCEkHG8IgqSZyfn6cDy4McJvDoWwB//u+sdo8sL1Wo6aYDmyoAtnmRbtrPKAFxqwbE9x4QnibFU40YBYmOvnHCHusUmDJEp513SpBwMBWJv9PpJmlhilbRokUVYWAKFp1myuDJVIxokJCEOgXL91HgeNaHsJiZ6VQsvmYkJhh54Fj/iMP58+cTT8HiWng0LOsqGDWxgoCcOHFCFZaTODACQueb69dDQDSyR6LF8dSPKW0sEmdjZIX1IowQkdzxOF3K1k7Y8sXK9x4QrRicr+shIEbWEIqAUBbJB0mrP+kUAmLzB50N4oWA2ACqWZFCQMwil3rHJcTH4695zZDr+ILAIDT4ArivQeoFyAMrd9OB9QA8Ea+im/YzSkDcOAUr4g0I3JSC5QWdndaR0SQSFatvQjeyDhbr83hmRs78mxAQI0hGRl8hIJFhB6WFEJAIMkakq5JwFZhaEDh3OLCmLXcCOZO50TzS15eK9HPTgU1FMNu2VDftZ5SAEASm35R/+OFk7wHZkj491lt0D4htwFso2C7n1UIVRVQIBOyyYbD3mOP+Wgp8AoSARJBRHX+gR0ap1V+TU7Ai6CkIocqV88C4LME7dTwFpE/mxCvvrDTVaOqmA5tqQLZxoW7azwwB0UjI882aJd6Enik+HhdiYrAzKgqFLL4J3UboLRNtl/NqmYIiKCQCdtlQCEhI6E13EAJiGjrrBwoBsR7TFCPxwu/ApBu32/q3q2mzIbr9MaSJTXoLcIpZewpfiJsObAqH1pHluWk/swREA2bz5s3qLoy/WAORMyeeeuqpsO/9cAR0iyexy3m1WE0RlwwCdtlQCIh9j50QEPuwNSzZLQKSUGcO0hRpZlhfGeAAAie3ATNKBZ7ogZZIqDYF27bvwIMPPpjkHH0HNJMpLELATQfWoiWkajFu2i9cApKqDeezeLucV8HXOQTssqEQEPtsKATEPmwNS3aNgDz5OdIUSv4MdsOLkQHhIbD/P8CiJoFlVJsIlGyvXnPT+QlvgTJaQ0Bs6O1nwU37CQGx5tmxy3m1RjuRogcBu2woBEQP+ub6CAExh5sto4SA2AKrd4Revw6sfR3YMiKwzs3WAnkqJnnNTefHO8BGtqZiw8i2Tyjt3LSfEJBQ1tH3ul3Oq77ZpZcVCNhlQyEgVlgnsAwhIPZha1iyEBDDkKWMAdevAdMfBE7tDLyetoeBLPcEfM1N5ydlgO/+KsSG7tsgHA3ctJ8QkHAs9+9Yu5xXa7QTKXoQsMuGQkD0oG+ujxAQc7jZMspRAsLd9lE3LoFKKD8YaR7pLZfU2WLVZITGXwbGpA/eoct5IF2mZLVy0/lxGq6UOp/Y0NuWddN+QkCseXbscl6t0U6k6EHALhsKAdGDvrk+QkDM4WbLKMcIyLlfgXnVgbMH1DquA4jKej/QdGXQnXZbFpxahV46DXyYM/Dqc5cBmv8XiI7RhY6bzo8uBaVTSATEhiEhiugObtpPCIg1j4Ze5/Xw4cPqpvRDhw6hTZs26NSpkzUKeEDKBx98oA46ceIiQtqjcePGGDFihLr1Xk/Ta0M9snz7CAExipj+/kJA9GNle09HCAgjH+NzAHFnbl5P7K1Ap9MSCbHL0n/tAz4pHFj6w72BCoMNY++m82MXTKlNrtjQ2xZ3037hEBA6bDyGd8GCBfjzz9PIlSuHcvoeeughSw3ywgsvYNu2beoWbbb8+fPjjTfeCHse6j1nzhzMmzcvqL68vXvixIlIly4doqKicM8996Bbt2544oknkozR67z27dsXadOmxVtvvWUpRpEu7O+//0bdunWxfPly7N27F23btk1UmX5L+vTpFb5sU6ZMQdmyZQ0tic9InTp18OyzzyaOW7RoEVauXImxY8fqkqXXhrqE+XQSAmIUMf39Uz0BOXfuHPr164e1a9ciY8aMalfjxRdfvAnBK1euoGfPnti1axeOHTum3mSPP/54Yr+NGzeiZcuW6o2otfbt2+Pll1/WbQ1HCMihVcCCqsF1avwtkL+Kbp2low4EDi0FFtQJ3LH+58D95k8gc9P50bFy6aIDAbGhDpAiuIub9jNLQPg99swzL6iLCK9fL4b4+MyIiTmPqKhdKFiwIObOnYkHHnjAEtR9nctr167hP//5j9rZ3rBhQyIpMTORXgJy8OBBjB49Wp0YOGPGDLz//vtYs2YNsmbNmjhtKOc1Pj5e6UrfoGbNmkkcZb26azL09o+kftOmTVO+D+3m2+Li4lCiRAl8++23yJMnj2mVAxGQy5cvKx9r8eLFyJUrV0jZoWwYUkCQDkJAzCIXelyqJyAkFf/88w+GDx+uiAU/YN59992bdkhIQGbPno1ixYqhR48eGDRo0E0EhKHZ9evXh0Y9SA9HCMiEvMDFo8F1zJAH6PCb6TXIQB8EfhgIbAiyU9ZiB5CreNhwuen8hK28CFAIiA29/SC4aT8zBISO5COPVMClS+Vw7VolALf4GOAyoqNXI336zdi4cYMlJMTfubx06ZK6t+i7777DnXfeqeb+8ssv1abeH3/8oQjQgAEDcP/996vXPv74Y0yfPh3cLMyePbtKgbrvvvvw9NNPg079Lbfc0J/fvRkyJL2MlREQjYCwD79jS5UqpUhQoUKFwLSir7/+WvkAjzzyiJo3W7ZsOHr0KKpWrYp33nlH9eEct956K3788UdFRPjz6aefokCBAhg6dKgiNExPql27tvIPGHEJJIORE+rfrl07td6rV6+qjc0iRYqgT58+akyVKlXUvJzjwoUL6vXt27ertRI36qjhRmzLlCmDrVu3KoLANY0cORJ33XXj0lqunbJ27NihfmcUo3///ur/33//vSJmR44cQd68edX8wSIX9IuefPJJhblv8ycg1HfYsGFYvXq1+lwjWevVqxdiY2Nx5swZ9O7dG1u2bFEi8uXLh0mTJin7Tp06NRHXatWqKX+M7aWXXlLz8oLMUE0ISCiEIu/1VE1A+GHEcDN3Uvihx8Y3JPM7kwv78QPi7bff9iYBGRUDXE8I/iRGpQFejY+8J9UrGjHFbc5jwO8bAmvc4U8gQ5D6DxNrdNP5MaGuDAmAgNjQ24+Fm/YzSkDopBUvXho//5wD167VCgp8dPQyFClyGrt2/RS2cXwJCPWdO3eucr6ZXkMne9WqVRg8eDAmTJigiMX8+fPV60uXLlUOeaNGjbBw4ULl7J88eVIREfYzGgGhsz9z5kyQlJAwjB8/HgcOHFAOM4kLMyGIDx14jTwwLYibjdSTJMSfTDGV7M8//1R+A53xjh07onz58ujevXtAGSQSdKqZadG5c2esW7cOr7zyCipUqICBAweqeejkc3y9evXUWhkpYsoYo0dMAaPfQsedjfpQV/5Oh55khXpyTSRV1P+ZZ55Rc7KRpJBkMI2KGRvEoHTp0goProWYk+T5t0cffTSxb3IEpEuXLiqThHpSX66Dm7b8l7gSb2JFgvbzzz+rdDz2DxQB4Tx8LtgoL1QTAhIKoch7PVUTkD179qBJkybYvXt3omX4BiT54L/BWjAC0qpVK7V7wjdXxYoV1U4Hf9fbtAgIyZD/To5eGaH6RY3LgqiEi7iRrZm0sRj9epoMuN7lXCgx8ro/AtfikWas705i0g4Jnf8BYmItx41f6Dt37kTx4sXlJnTL0XVGoNjQGZztmsVN+3Hu//3vf2oDjTvwodqmTZtQsWLl/7/z3gdAcp9Hl5E27VCsW7ca5cqVCyU22dfpXPIzirUTdNLZhgwZohxsNkYDKleunCStqXr16iqykDt3btWPDjWdcC3awXEkICQz/AnWSDYmT56sxtG5p5NOh//hhx9WjjcjISxypvPKjccGDRqoepXjx4+DO/FLlixRxEdrXAujCHTqiT0jEiRMhQvfqO1jVIFEhhEAEgN/GUzVJvlgJIV4sDGCQfJBuWwcnylTJhU58G8kDs2bN1fj2agPN1Hp+LOtWLFCESzWTzB1ifUv/L9/YxSFjj8Ji9YY5WjYsKH68W8kEb4btdrrtGfJkiXxzTffqPRz2ui///2v0p+N0Q4SG75Ov4pkinMzUuPbfHH1/TvJCqNi7733XlAbay/Qhoyu+dajhByko4NG8kl6fd9j9NeYwsjolV3+mg71PN0lVRMQvjl4igU/FLTGMO6bb76pakKCtUAEhLsgLNTih9mJEydUkVp0dLT6ANDbNAKit7+ZfnftG47cx+cGJSAn7miGY4VeMyM6VY6Jjr+AUuuYxnBzO5OzEg4+MAyIunHcsTRBQBAQBKxGgI41nWR+34RqTL8ZO3Yd4uObhOqKmJj56Nr1MeUch9NYsMxUHO7s00nkLjh3/5lFwGgB02voZPo6d0w3oq61atVSTjULzfft26fSp7ibzp3zr776ShEIpmcFa/z+5clVTKv2bX/99ZciB5qjrL3GVGvK5fxM/aHD7Et6fNdy+vRpkCgxepA5c2YlgiSmadOmINEjifGXQZ+D/gWjP1qjDKZJaelPTD+iHkyJokPNyAH1OH/+vBrCyAadfG50+urD1+i/kLiRfLBug6lXo0aNugkeEhambWkHA7AD10x5WrTEdxDT0SiHZMO3kYAwOsI0NmLKqAqJjW+jA0+96N+QDDL1jutidIb+F4mY/zq08SSeLG5/7TX3fBJGcpjKRnwCNSEg5j8dUjUBYQSEHxYMS2pt2bJlGDNmjOEIiL8JfvvtN9SoUUPtVPgWpidnKiciIEhIQPS4GztfvlEQRj/YrnWJA3TspJl/5FLIyHO/Is3UwMcDXqs0BtcfdOZ4Rjd3X1OIJV1fhtjQdROEpYCb9jMaAWndui0++eQwgODpV/+CsRStWxfAlCmTw8In0O42oxA5c+ZUu/2MCNAZ5QlcyTXNGWfqzqxZs1RaFk/BChUBISnwd8LpVDICQseZxdP+u+da9IIpU6xf0JqZCIivDG52skaEqVdae+yxxxTJYFSGjQSCjj0jBR9++KFy3hkJuO2221TqFCMUmkx/bBmB4eYn09qSi4CwD/HXIiehDBysFsM3AkKcuDlLYkNylFxj3QlJB39ITFu0aJGYLuY7jlkljAzprQGRCEgoS0bW66magGg1IPwg0wrezNaA+JuVBe3cNSAB0Ruec6QInYru/xxY9G8xmboHhH+vvwC4v1FkPaGRps1va4B5gSMeaLoayJv0eEe71Xcz/9zutaUW+WJDb1vaTfsZrQFhOszIkd/pjoD06FH5puiBUWv55/czZYw75XR+mcrE9ByersTTqZiaw+9BOupMLWLNB6MjTFPibj03B3/66Sd1mhWzFEhgGE0I5vD6F6H76k5Hn9/TdMZ5yhJJB6MsjIxoNSCMIPgTEN/jYpkmxZ1/EghGLbijTyLB9OtAMrgu/8NqWP9BgqQREKan0bFn5IkRAJIOEhH+jevlUbiaXv7YEhNGlkhAWBDOonj2oYPPptWAMO2cp3QScxIx6s7UM6ao3X777TeZmNEUbthSH9/mX4TO9ZPYcI1ZsmRRtmPEi6dZMfJB+TwKmdki1Kt169aKeDKqxXl9084om6nsTCFjKl6oJjUgoRCKvNdTNQGhObgbQdbMN9bvv/+uwo/8YPI/J5x9+SblQ86wMD8I+MHB8CFD3wyJcieFp0+cOnVKhY/Znyc86G2OERAqlJAAvH/jXPaEexshTb35EvlIzlBbRgJr/s2XTdK1zSEgaz69Zra0n5vOj6ULScXCxIbeNr6b9jNKQJga9NhjrAHprasGZP36NZbUgPjeA8LTpOrXr4+uXbsmpo0xEsHCczrtzBgg4aAjzt/5XUvSwhStokWLKsLAVGd+v1IGN/kY0aDzHeoULN8njeNZvM2UK35ns/iau+3ByAPH+jv8TIvSTsGiH0DfgOlCJC1WEBCmc2vH/zMCQv+E69dDQKgvcaN+JB7UjylhWkE3Iyusy2B6Eckdj9OlbO2ELV+sfO8B8U1bC3QKFkkNj+U9e/asIhUkmSRAJDEkjkxdowwSORIO2pWkkuSYZI5RFNZ8MILDjBSSSD1NCIgelCKrT6onIDxlgm9Ihi797wFhvqnvpTp8Y3DHxLcx/5Q7F5988ol6g/GNSuZP5s4PjkAnSgR7BBwlIFRi5I0krIRGS5CmQO3IejIjRZv/1AR+XRFYmy7ngXQ3iu3cam46P26tOaXNKzb0tkXdtJ9RAuLGKVhesK5dzqsX1q5XRxIBEhUnb0JnPQyLv/U0u2wY7D3muL+mBwSP9Un1BCSS7OX4Ay0EJLD5r18DRgU5USZDbqD9USD6RvTI7eam8+P22lPK/GJDb1vSTfsZJSBEmuk3Dz9cPsQ9IFuwceN6S+4B8YJ17XJevbD2lKKjXTYUAmLfEyIExD5sDUt2lIBICtbN9rl6CRib9CKrxE4sKq8yDogKdICxYVNbNsBN58eyRaRyQWJDbz8AbtrPDAHRSEizZs/73ISeCTExFxAVtRMFCxay9CZ0L1jXLufVC2tPKTraZUMhIPY9IUJA7MPWsGTHCEiwIvR6/wEKhr5x1PDCIn3AhePApBs38t7U6s4BCjeL2BW46fxELCgeU0xs6DGD+anrpv3MEhBtCZs3b8bnn3+OU6f+Qs6c2dVpQ+He++FFa9rlvHoRC6/qbJcNhYDY90QIAbEPW8OSHSEgPpGPgAp2i089xejHNwKzHwlspxbbgVwlDNvQ6QFuOj9OrzWlzic29LZl3bRfuATE28hbp71dzqt1GoqkUAjYZUMhIKGQN/+6EBDz2Fk+0hECsrwLsOuD4LoX6wzU1HfqhOUAOCVw+0Tgmw6BZ+twEsiQyylNwp7HTecnbOVFgEJAbOjtB8FN+wkBsebZsct5tUY7kaIHAbtsKARED/rm+ggBMYebLaMcISDvZwASLgXXP016oNtFW9bnutAvGgK/fBlYjVcuAzH/Xjjluq46FXDT+dGponQLgYDY0NuPiJv2EwJizbNjl/NqjXYiRQ8CdtlQCIge9M31EQJiDjdbRjlCQEalBa7HB9c/KgZ49aot63NF6PXrwKjowFPnLgM8twmICvK6Kwobm9RN58eYptI7GAJiQ28/G27aTwiINc+OXc6rNdqJFD0I2GVDISB60DfXRwiIOdxsGeUIAZl0N3Dht+D6Z8oLtD9iy/ocFZpwBXg/SESjwiDgkb6OqmPXZG46P3atKbXJFRt62+Ju2k8IiDXPjl3OqzXaiRQ9CNhlQyEgetA310cIiDncbBnlCAE5tApYUDW4/o2/BfJXsWV9jgi9eAqYEKSG4+mVwD3VHFHDqUncdH6cWmNKn0ds6G0Lu2m/cAgIb8ke/+F4LF25FBfOX0CmzJlQu3ptdOrYSfflb9623L/a+zqvvGy4atWqibeN271G3ho+ePBgdYnxxIkT1cXGqaXxhvRXXnlF3XwfbvO1YenSpbFgwQLkz58/WbEc07hxY4wYMQL33ntvwL5CQMK1TPDxQkDsw9awZEcICFOSxucA4s7crF9sdqDTqYi760IXkCd+AmaWDty1zUEga/IfRLrmiMBObjo/EQiHJ1USG3rSbIlKu2k/MwSEjm7zF5pjxbIViCkQg7g8cUA6AFeA2KOxiD8Yjxq1amD2jNnIli2bJcZZvnw5PvnkE+zbtw/p0qXD7bffjjp16oAOaPr06S2ZIxwhRgjICy+8gG3btqlbwbmWYsWKoU+fPihQoIApFWrUqIHu3bujdu3apsZ7ddC6desU4Zo5c6b6d9KkSWopfKavXLmS5Ln46aefQi7TbARk0aJFWLlyJcaOHSsEJCTK1nYQAmItnmFJc4SAUMNzvwLzqgNnDyh9rwOIyno/0PQbIMvdYa3B8cG7PgGWtwo8bZdzQLrMjqvk5IRuOj9OrjMlzyU29LZ13bSfUQJC8vFw+Yfx6+VfEfdoHBDo3tWLQOwPsbjnlnuwccPGsEnItGnTMGHCBPTt2xeVKlVC5syZ8csvv+Czzz5D06ZNUbBgQcceADqp165dQ5o0aZLMaZSAkDw9++yz+Oeff9CvXz8wajJ37lxD64iPj1ckpmjRoqATHGwHPjmhmgxDE0dI544dO6Jy5cpo0qRJEo3Wrl2Lt99+G6tWrQqoabA1myUgly9fxuOPP47FixcjV66bsyckAmLfAyMExD5sDUt2jIBQs/h4YEzaGzsOeSohzVMrgZgYwzq7NmBRU2D//MDTd78KRHtoLWGA6KbzE4baMtQHAbGhtx8HN+1nlIDUqVcHq3avQlzlOCC5szeuAbHfxaLKA1WwZNES0wY6f/68cu6GDBmiIh7BGp1HEpU5c+bgzJkzKFmyJAYNGqQiJWyFChXCwIEDVRTlzz//VI7r0KFDVQSC7fvvv8fo0aNx5MgR5M2bV0UkypYtq15jxKJUqVLgLvqOHTvUPOfOncP777+PX3/9VRGiRo0aoW3btsiQIYMiE8mlYFGeRkAof/Xq1SqCQfnUjWvdtGkT0qZNi6effhqdOnVCdHS0Sgni+njRI//POej08nufUaDY2Fhs3LgRhw4dUmvdtWsXsmfPjjZt2iQ66f4yqlevrpzm/fv3I1OmTGCkKWfOnCql6ODBgxgzZoyS36FDB7z44osKj507dyodmYLHOSmjd+/eiViGwnrNmjVKLrHjnEyhYhpTKBv62v7q1avKPl9//bWyl2/zJyBvvPGG0u306dPYsGGDei7uueeem9bw5ptvguSENixcuDCWLFmiSB3H33LLLTh16hTWr1+PPHnyYNiwYShSpEjitC+99BKefPJJdRGnfxMCYvrtH3KgEJCQEDnXwTEC4tWb0JM70arAk0DDr7yZPhbGI+am8xOG2jJUCEiKeQbcfA8aISB0OAsXKYyEpgmBIx/+FrkIpJmXBvv27jO1O09xdCZffvllbN++XTnkwdqMGTOwcOFClQaTO3dufPDBB8qJZ5REIyAVKlTAyJEjVQSjWbNmSi4d/L1796Jly5YYP348mPtPB5lO59KlS5UDT8JAp37KlCmKyNBJpT4kHoy+HDhwAK1atcLrr7+O+vXrGyIgFy5cUBGQkydPgmtgRId6cnef0ab27dujefPm6u8kD4wCvfrqq0pf2o6OMXXSnGU65nSESXBIGpiyRgJCclW+fPmAMrgupi+RFDDCRPKxbNkyVKxYURELynjuuedUmtEdd9yBPXv2gLv+JUqUwIkTJ9CuXTtFIFq3bh0Sa5Ii4kl9KP/s2bP4448/VBQnlA19bc9nsWHDhopk+bdABITEilE01sfExcUpcuW/BpJIRqUCEZBvvvlGpXmRiL777rsKg1mzZiVOzRocNtpHCIhzH81CQJzDOuRMjhAQL96Efi0eGB3ky6v6JKBEu5DYptQObjo/KRVTp9clNnQacWvnc9N+RghI91e7Y8LiCYh7PE43ALFrY9GhbgeMHjVa9xjfjl9++aXabebOs9bonJNc0Nmmw8cUHDrcvXr1whNPPKG6kSQ8+OCDWLFiBe68807lpE+fPj2xQJsOI4lI//79VbpOxowZ8dprryXOQQefDimdXDrMjKj07Nkz6BoYTaFDO2DAAF0EhJEU7sozgkDZJDyM3HBtjMZERUWpub744gtFGqg7/x01alSS1zVypRGQLVv8c2XzAAAgAElEQVS2oHPnzgovLU2MhIJEYfjw4QFljBs3TuFJAsC2e/duRShYY6GlFNWsWVPpyMiRf/v000/xww8/KAdd0ycY1sSbaWP817+FsqFv/61btyqsGPHxb4EICEkHSU+wxjUwOkJ8AxEQjiPx0PB5/vnnVcRKa5R9/Phx9az6N4mAmHrr6xokBEQXTM50coSALOsM7B4ffEEPdAJqJXNTujNQ3Jjl8t/A+FsDz/jcZuD2GyH21NzcdH5SM+5Wrl1saCWazsty035GCEjh4oWx77Z9gJFa6YNAoZOFsHfnXlPAJhcB8U1lohNPh5upSlpjITLTpRjV8I0S8HU65UypoVPJ1Ckt5UkbSwJDB5e7+5yHBd6MRGiNERDKYPSDRIhz0TlnFMFoCpYmkySCJIcOsNZIkhh1YKoVCQgjOvPnJ00d9l0bZUyePFkRF61xDKM5Gonxl0ECwoiA5qCzvoZkgJEPrTVo0EBhUbduXRUNIm6MPly6dElFYpiypNWwhMKaER4tncv3oQhlQ9++RiMgOXLkSEIwA62BejNFLxABYVqaRkAD4SMREFNv77AHCQEJG0LrBDhCQEYyZza5iwbTAj2uWLcoM5JO7QY+LRZ4ZIcTQIbbzEhNkWPcdH5SJKAuLEps6ALoFk7ppv2MEJA8+fLgWJFjQB4Diz8K3PXzXTh6+KiBQf92ZQ0IU3Xeeeedm0558iUgtWrVUtGHYEfQJucUczeeO/1dunQJqKN/zQY7VatWTaXrMDWJaVCMgHAHnClgZgkIT8ZielWw4mmtfmPevHlJ9PRdm54ICOtIfGUYJSCMDjH1jLUbrOFg9IAESZMZCutgEZBQNvRdNAkfa2H01oD4EgjKCbYGswREakBMvb3DHiQEJGwIrRPgDAG5ERpOtvXguVgutD0zgKUtAk/8yiUg5hYXlIrsKd10fiIbGe9oJzb0jq0Caeqm/YwQEDciIMRr6tSpqv6C6VZMsWK6FAuYWcNBR5JEgE4wd/mZAnP33Xer2gKmIWmF68k5xUw5Yq0Fi8oZLaFzSzKQL18+VcQeiIA8+uijiiww/YtF2RzPuyjCISC0BWs9qlSpAjq0JDYsimd9yEMPPZRYhJ4cAWE0hlGKevXqKZ0YoWFtBqM1jz32WEAZRgkI62ZICrt27YrDhw+rSBHrYfQQEGJFmzHaQn1oJ6aHsaA7lA393zuscWEhPvXxbYFSsPwJSLA1mCEgTO8iHjyJjPVH/k1SsOz7fBYCYh+2hiWnWgLy9TPAviBHGL6aAEQld1yLYZhT1AA3nZ8UBaSLixEbugi+BVO7aT8jBMSNGhANXqYWMYWIBeOsmyAxoKPNtCjuwjNViUXBs2fPVg5tlixZVNE1IxNsyREQvk6yQvLAVCTu0LPA+q233lL1I4EICIu033vvPVUoTnLAfkzpCoeAUA+egsVaDdYjML2JJzwxRYxr1RMBoQymCPEULBIrFtGTgLDoni2QDKMEZPPmzapwXiMOjDoRPz0EhDowwkOcSCKzZs2qIimstwllQ/+3GmtlmG6m1a5or+shIMHWYIaAMPrD54E4Bttg4CljjBr5Ht/siL9mwedTJIsQAhJB1nHkgR6dGbh2IfiqozMB3c/bj0pyJ1oVawXU/Nh+HVLADG46PykAvohYgtgwIsxgWgk37WeEgLhxCpZpUB0eaPYOCYfVTHHTkRx269bN8pvQtUMAQgGm3YRO0njfffcJAQkFmMWvCwGxGNBwxDlCQC5eBCZkDK5mh38AnyK6cNYTcOy1BGB0kDs6npwHFEp6KZHl86cwgW46PykMSteWIzZ0DXpLJnbTfkYICBfr9D0glgDsgBAhIA6AbPMUdtlQUrDsM5wQEPuwNSzZEQJy/jwwOUtw3dqdAzLbcHt43Fngg2yB531pH5DdudtwDRsmgge46fxEMCyeUk1s6Clz3aSsm/YzSkCM3oS+6YdNKs0mpTe7nNeUjlskrc8uGwoBsc/KQkDsw9awZEcIyMg0AK4lo1s00CPBsO5BB5z+GZhWNPDLnc8CscmQIeu0SLGS3HR+UiyoDi9MbOgw4BZP56b9jBIQLp0kpPkLzbFi2QrEFIhBXJ44gIcjXgFij8Yi/mA8atSqgdkzZiNbtiCbRhZj6LY4u5xXt9eVmua3y4ZCQOx7ioSA2IetYcnOEBCHTsH6eRaw5PnAGHS/CkQHScMyjFrqHuCm85O6kbdu9WJD67B0Q5Kb9jNDQDSMWBMy/sPxWLpyKS6cv4BMmTOhdvXa6Nyps+mbz93A34o57XJerdBNZOhDwC4bCgHRh7+ZXkJAzKBm05gUQUAWNQP2Jz3nXMEVmw3o9BfwfzfE2gRhqhPrpvOT6sC2acFiQ5uAdUism/YLh4A4BI8nprHLefXE4lOIknbZUAiIfQ+IEBD7sDUs2RECMqkAcOFQcN0y5QfaHzSsO0YGiayU6QFUGmFcnozQhYCbzo8uBaVTSATEhiEhiugObtpPCIg1j4Zdzqs12okUPQjYZUMhIHrQN9dHCIg53GwZ5QgBOXMGmJo9uP6t/gJuvVXf+q7FA6PTBu771HIgXw19cqSXaQTcdH5MKy0DkyAgNvT2A+Gm/YSAWPPs2OW8WqOdSNGDgF02FAKiB31zfYSAmMPNllGOEJBx2YArZ4Prny4r0OXv5NcXdw74IMjJKG2PAFny2oKPCL0ZATedH7GHNQiIDa3B0S0pbtpPCIg1VrfLebVGO5GiBwG7bCgERA/65voIATGHmy2jHCEgwVKlfFfU43rg9Z3eA0x7IPBrXS8CadPbgosIDY6Am86P2MUaBMSG1uDolhQ37RcOAWER+sQPx2P1N0tx/sIFZM6UCZWq1cbLHTsFvZTNLYztntfXeT127BiqVq2KHTt2qFvb7W7ffvstBg8erE4nmzhxIngzeWppLVq0UDeplylTJuwlhyIgvPn9xRdfxJdffol06Xjsm74mBEQfTmZ6CQExg5pNYyKWgOyZASxtEXjVryYAUdE2ISJiQyHgpvMTSjd5XR8CYkN9OEVqLzftZ4aA0NFt1aI5Fi9bgUYlYvBU0ThkSw/8fQn4fE8sFu6IR91aNTB1unXH8C5fvhyffPIJ9u3bp5y/22+/HXXq1AEd0PTp3d+4MkJAeHv3tm3bEBMTo9ZSrFgx9OnTBwUKFDD1iNaoUQPdu3dH7dq1TY336qB169YpwjVz5kz176RJk9RS+ExfuXIlyXPx008/hVymPwEpVKgQlixZkuREt379+oF/f/75ICd0BphFCEhI6E13EAJiGjrrB0YcAfmyMfC/hTcvNGt+oI2JQnXrIUv1Et10flI9+BYBIDa0CEiXxLhpP6MEhOTjiQoPI0/0r/iocRzuCHAN0/FzQJsFsTh67R6sWb8x7LtApk2bhgkTJqBv376oVKkSMmfOjF9++QWfffYZmjZtioIFnbuElk7qtWvXkCYN78P6txklICRPzz77LP755x/QqWXUZO7cuYaewPj4eEViihYtikWLFpk6+liTYWjiCOncsWNHVK5cGU2aNEmi0dq1a/H2229j1apVhjTVQ0C2bNmiZH/99de6ZQsB0Q2V4Y5CQAxDZt+AiCEgwZZYfgDwaH/7ABDJhhFw0/kxrKwMCIiA2NDbD4ab9jNKQBrXr4O4Q6vw5QtxiEnqgycxQnwC0GBGLGLzV8GCr5aYNtD58+fx+OOPY8iQISriEazReSRRmTNnDs6cOYOSJUti0KBBKlLCxl3rgQMHqijKn3/+qRzXoUOHJqbSfP/99xg9ejSOHDmCvHnzqohE2bJl1VhGLEqVKgXuojOtivOcO3cO77//Pn799VdFiBo1aoS2bdsiQ4YMikwkl4JFeRoBofzVq1erCAblUzeuddOmTUibNi2efvppdOrUCdHR0ViwYIFaX7ly5dT/OcfixYvB731GgZjutXHjRhw6dEitddeuXciePTvatGmT6KT7y6hevTpy5cqF/fv3I1OmTGCkKWfOnBgxYgQOHjyIMWPGKPkdOnRQ6UdsO3fuVDoyBY9zUkbv3r0TsQyF9Zo1a5RcYsc5mULVuHFjhLKhr+2vXr2q7EMiQHv5Nn8CcvjwYZWiRr1pq5deegnPPfdc4lqIFQkto1GPPPKIeg5IDmkP4hoVFYU33ngDzZo1Awkb070CzRvs2RQCYvrtH3KgEJCQEDnXIWIJSLO1QJ6KzgEhM+lGwE3nR7eS0jFZBMSG3n5A3LSfEQJCh/OBooVx+M2EgJEPfyswEpLvnTTY8/M+U7vzlEdn8uWXX8b27duVQx6szZgxAwsXLsTYsWORO3dufPDBB8qJZ5REIyAVKlTAyJEjVQSDziTl0sHfu3cvWrZsifHjx6N06dKgg0yHc+nSpcqBJ2GgUz9lyhRFZOiEUh86s4y+HDhwAK1atcLrr7+O+vXrGyIgFy5cUBGQkydPgmtgRId6cnef0ab27dujefPm6u8kD4wCvfrqq0pf2u6WW25ROmmpQnTMn3zySUVwSBqYskYCQqe6fPnyAWVwXUxfIilghInkY9myZahYsaIiFpRBh33lypW44447sGfPHly+fBklSpQA6yLatWunCETr1q1DYk1SRDypD+WfPXsWf/zxh4rihLKhr+35LDZs2FCRLP/mS0AuXbqkUtOoIyMlv/32m7IVCRRx5nPANRNnpm1t3bpV4UTSESgFi3PVq1cPXbt2VcRLTxMCogclc32EgJjDzZZRjhCQ97MACef16d/hBJDhNn19pZcrCLjp/Liy4BQ4qdjQ20Z1035GCEjPV7vj6NoJmPNsnG7Am82ORd4nOmDEqNG6x/h2ZMHvsGHDsH79+sQ/0zknuaCzTYecjiUd7l69euGJJ55Q/UgSHnzwQaxYsQJ33nmncianT5+eWKDNHXESkf79+6uUmowZM+K1115LnIMOPqMadHLpMDOi0rNnz6BrYDSFTvmAAQN0ERBGUrjjzggCZZPwMHLDtTEaQweY7YsvvlCkgbrz31GjRiV5XSNXGgFhilDnzp0VXlqaGAkFicLw4cMDyhg3bpzCkwSAbffu3YpQsMaC0RG2mjVrKh0ZOfJvn376KX744QdVh6HpEwxr4s20Mf7r30LZ0Lc/iQKxYsTHv/kSEOLCqNf8+fMTu5FoHj16FO+8846q5cifP7+KMpG40odiFCs5AvLMM88o4sofPU0IiB6UzPURAmION1tGOUJARmcGrl1IXv9uV4A0wXerbFm8CDWFgJvOjymFZdBNCIgNvf1QuGk/IwSkbInC6PXgPjQpqR/veduAYTsKYcv2vfoH+fRMLgLim8pEJ54ON1OVtMYdbaZLMarhv5tNp/zUqVN49913VeqUlvKkjSWBoYPLnXPOw110RiK0xggIZTD6QSLEueicM4pgNAVLk0lnmSSHDrDWSJIYdWCqFQkIIzq+zrTm8GsEhP9OnjxZERetcQyjORqJ8ZdBAsJ0K0Yl2JiORDLAyIfWGjRooLCoW7euigYRN0YfGGHgM1S4cOHEGpZQWDPyoKVz+T4UoWzo21dvBITRHdqEkSKtUV+mb/E1ptwxakbCxmgXbcwfiYCYers6PkgIiOOQB5/QEQIykh/wQY7ZVapFAT2uRRAqokpyCLjp/IhlrEFAbGgNjm5JcdN+RghIoQJ58EH1Y6huoOZ7xT6gyzd3Yd/Bo6bgZQ0IU3W4W+1/ypMvAalVq5aKPgQ7gjY5p5i78dzp79KlS0Ad/Ws22KlatWqqToCpSXRuGQE5fvy4cmbNEhCejMX0qmDF01r9xrx585Lo6bs2PREQ1pH4yjBKQBgdYuoZazdYw8EICAmSJjMU1sEiIKFs6LtoEj7WwoSqAaFeXK8W3Qn2ELL+hNEUppFxTL58+QKmYEkNiKm3sW2DhIDYBq1xwY4QkHG3AleSuWgwXTagyxnjyssIVxBw0/lxZcEpcFKxobeN6qb9jBAQNyIgtOzUqVPVbjXTrZhixXQpFjCzhoPOMIkAnWDu8jNd6+6771a1BdzV1grXk3OKmXLEGgAWlTNaQueWZIBOKIvYAxGQRx99VJEFpn+xuJnjWZwcDgGhLVjrUaVKFVUoTWLDHXrWhzz00EOJRejJERBGYxilYJ0CdWKEhk41ozWPPfZYQBlGCQhTj0gKWQfBAm9GilgPo4eAECvajNEW6kM7MT2sSJEiIW3o/y5njQsL8f1ToXxTsHjKGLFgX0ZxGCFjhCcuLk7VsDBSxLXkyJFD1bYQfz5HLGxnpIa1IqwR0RpTv1izw0iT3iYpWHqRMt5PCIhxzGwb4QgBuXgRmJAx+Bo6/AP4hJBtW6wItgQBN50fSxYgQlQKBB0m5rz7Hw8q8EQ+Am7azwgBcaMGRLMeHT6mELFgnHUTJAZ0tJkuw114pirNmjULs2fPVg5tlixZVDExIxNsyREQvk6yQvLAVCTu0NM5feutt1T9SCACwiLt9957TxWKkxywH1O6wiEg1IOnYLFWY8OGDSq9iY4wU8S4Vj0REMqgg82TnUismFZEAsJia7ZAMowSkM2bNysnXCMOjDoRPz0EhDowwkOcSCKzZs2qIimstwllQ/93MmtlmG7mH90IdAoWbcXPSEYweN9Kt27dQBLJuh/qTqwZBSM50lKwGDnhYQas7eEBAyQnjJbdf//96pnQ24SA6EXKeD8hIMYxs22EIwSE2s+tBhz9NkkiliqZu7sm0GSZbesTwdYj4KbzY/1qUqdEsaG37e6m/YwQEDdOwfKKZUPdou2VdXhNTxIBkgm5Cd1rlrNGXyEg1uBoiRTHCAi1vXgRCR/fiagrZ3E9XVakaf27RD4ssaKzQtx0fpxdacqdTWzobdu6aT8jBIQoO30PiFcsKwTEK5YKrqddNpQIiH3PhhAQ+7A1LNlRAgJI6odhC0XeADedn8hDw5saiQ29aTdNazftZ5SAGL0Jfe2GTSrNJqU3u5zXlI5bJK3PLhsKAbHPykJA7MPWsGQhIIYhS/UD3HR+Uj34FgEgNrQISJfEuGk/owSEEJGEtGrRHIuXrUDD4jF4qmgcsqUH/r4EfL4nFl/sjEfdWjUwdfpsZMuWzSVUnZ3WLufV2VWk7tnssqEQEPueKyEg9mFrWLIQEMOQpfoBbjo/qR58iwAQG1oEpEti3LSfGQKiwcSakIkfjsfqb5fi/PkLyJw5EypVrY0OnTqbvvncJROEPa1dzmvYiokA3QjYZUMhILpNYLijEBDDkNk3QAiIfdimVMluOj8pFVOn1yU2dBpxa+dz037hEBBrUfC2NLucV2+j4i3t7bKhEBD7ngMhIPZha1iyEBDDkKX6AW46P6kefIsAEBtaBKRLYty0H48+5Y3XPFqUx89KM4eAXc6rOW1klBkE7LIhj/7lfSw8Cpr3kGjNaX/NDCaRPkYISARZyOkH2s0vzgiC3dOqiA09bT6lvNjQ2zZ0236HDh1SjlHu3LmRNm1ab4PpkvZ0XnmXRPr06REVpQ6ll+YxBOywIS+G5H0pJPr58+dPgojT/prHzKFLXSEgumByppPTD7TbX5zOoJqyZxEbet++YkNv29Bt+9E54m3bLC6nEybNOALEjc4mCZwQEOP4RcIIO2zIZ4EHMdx2221Joh9cr9P+WiRgbLUOQkCsRjQMeU4/0G5/cYYBlQz9PwTEht5/FMSG3rZhpNiPDpj2421EndeeNty5cyeKFy+ONGnSOK+AzBg2AlbbkORD+wmknNP+WtgARaAAISARZBSnH+hI+eKMIBN4ThWxoedMdpPCYkNv21Ds5237UXuxodjQKAJO+2tG9fNCfyEgEWQlpx9o+dCNIOObVEVsaBK4CBomNowgY5hQRexnArQIGyI2jDCDmFDHaRs67a+ZgCTihwgBiSATOf1AO/2GjSCoU4wqYkPvm1Js6G0biv28bT+JgHjffm7Y0Gl/LWVYKekqhIBEkFUvXLigjlTMly+fOo3D7sYvzv3796NgwYKS92o32DbJFxvaBKyDYsWGDoJtw1RiPxtAdVik2NBhwG2Yzmkb8tS0w4cPq+N5M2XKZMOKUr5IISARZOPTp0+rB1qaICAICAKCgCAgCAgCgkBkI8AN4xw5ckS2khGqnRCQCDIML7w5e/YsYmNjbzryLYLUFFUEAUFAEBAEBAFBQBBItQjw+Ou4uDhkzZpVLgE1+RQIATEJnAwTBAQBQUAQEAQEAUFAEBAEBAHjCAgBMY6ZjBAEBAFBQBAQBAQBQUAQEAQEAZMICAExCZwMEwQEAUFAEBAEBAFBQBAQBAQB4wgIATGOmYwQBAQBQUAQEAQEAUFAEBAEBAGTCAgBMQmcDBMEBAFBQBAQBAQBQUAQEAQEAeMICAExjpmMEAQEAUFAEBAEBAFBQBAQBAQBkwgIATEJnAwTBAQBQUAQEAQEAUFAEBAEBAHjCAgBMY6ZjBAEBAFBQBAQBAQBQUAQEAQEAZMICAExCZwMEwQEAUFAEBAEBAFBQBAQBAQB4wgIATGOmYwQBAQBQUAQEAQEAUFAEBAEBAGTCAgBMQmcV4adO3cO/fr1w9q1a5ExY0a0adMGL774YkD1N23ahIEDB+K3337DfffdhyFDhqBw4cJeWWqK1VOvDbdt24Zx48Zh165dCouSJUuid+/eyJcvX4rFxisL02tD3/UsWLAAb775Jt5++208++yzXllqitTTiP3i4uIwfPhwLF68GPw/33/Tp09HpkyZUiQ2XlmUERsuWbIEH3zwAY4fP45cuXKhbdu2aNKkiVeWmiL1nDlzJviZuH//flSvXh2jR48Ouk7xZbzxCAgB8YadTGvZs2dP/PPPP+oL8dixY4p8vPvuu3jiiSeSyDxz5ox6U/ft2xd16tTBrFmz1Jfm8uXLkS5dOtPzy8DwEdBrwzVr1ihbV6xYEbGxsRgzZgxWrVqFpUuXhq+ESAgLAb021Cbh+/GZZ55B2rRp8dxzzwkBCQv98Acbsd8bb7yBixcvon///siePbtymAoUKCCfo+GbISwJem34+++/o0aNGhg7diwqV66M7du346WXXlLfiUWLFg1LBxlsHoEVK1YgOjoaGzZsAD8fgxEQ8WXMY+z0SCEgTiPu4Hz8EnzooYfUrkHBggXVzHzTHjp0SH24+rZ58+Zhzpw5qi/b9evXUalSJQwYMED9K80dBIzY0F/D06dPo3z58vjvf/+LW2+91Z0FyKzKGdX7PtTgYuSDESzuonNDQCIg7j1IRux38OBBPPXUU+BmQJYsWdxTWmZOgoARG27duhWdO3fGDz/8kCiDNn3hhRfQsGFDQdZlBBjl5/ssGAERX8ZlAxmYXgiIAbC81nXPnj0qbLx79+5E1bkbTvLhvys+ePBgXLp0SaVdaa1du3YoW7Ys+K80dxAwYkN/DWlj2nPdunXuKC+zKgSM2pDpAyNGjFAbAi1bthQC4vJzZMR+X3zxBT7++GNUqFAB/D+Jf6tWrSR9x0M2jI+PR4sWLVTUo2rVqvjxxx/RqVMntTl31113ubwSmT4UARFfxjvPiBAQ79jKsKZbtmxRH5wbN25MHLt+/XqVV86aEN/GWoGsWbOiV69eiX/u0aMH8uTJg+7duxueWwZYg4ARG/rOyDqeZs2aJabUWaONSDGDgBEbXrlyBY0bN8Z7772HBx54QO26SgTEDOrWjTFiv4kTJ6qd2Zdffll99u7du1cRkA8//FBFwaS5g4ARG1LD+fPnY+jQoaqGJyoqStVhSQ2IO7bznzUUARFfJjLspEcLISB6UPJoH+7cNW3aNLEomctYtmyZqg0IFAG5fPkyuHugtfbt26NMmTISAXHR/kZsqKnJwsnnn39e/XAXT5q7CBixIR3VU6dOqfoBNiEg7tqOsxux37Rp0zBs2DDwQAitdo51dZkzZ06yueP+qlKXBkZsyIgxN90mT56s0iB/+eUX8LuQ70lJR3b/uQlFQOjDiC/jvp30aCAERA9KHu2j5b0uXLgQ999/v1pFcjUgc+fOxeeff676sQaEBXjc+ZEPXfceACM2pJZ//PGHSh94+umnhTi6Z7YkMxuxIQkHi5ZjYmKUjLNnz6oDBVgU+84770TIilKXGkbsx7qB1q1bCwGJsEfEiA2ZQsc0yEmTJiWugrWQjIRoGwMRtrxUpU4oAsIaEPFlvPFICAHxhp1Ma8k0KtZ2cFeOp3twR5yh5WCnYPHI3tq1a2P27Nngbh5PnpBTsEzDb8lAvTY8ceKE2jGvX7++KqKUFjkI6LXh33//jatXryYq3qVLF5WHznQ6KWp2z5567cf6gbp166qfjh07KjLJOh5GtsqVK+feAmRm6LXh5s2bVfrcRx99hBIlSqiCZx5fz7Q6ZhRIcwcBvrcSEhLUe+nw4cPKp+GpWDwp0Ldpp2CJL+OOnYzMKgTECFoe7Muzz5kC8P333990D0ipUqUwZcoUVWjOxlqRQYMG4ciRIypiwlBmkSJFPLjqlKWyXhvy3HruDmXIkCEJADxJ6c4770xZoHhsNXpt6L8sScGKDEMbsR9Tduj8MO3ntttuU+k7PEVJmrsIGLHhZ599pjbgTp48qWojuanTrVs35fBKcwcBfrfxO863NWrUSF0rIL6MOzYJd1YhIOEiKOMFAUFAEBAEBAFBQBAQBAQBQUA3AkJAdEMlHQUBQUAQEAQEAUFAEBAEBAFBIFwEhICEi6CMFwQEAUFAEBAEBAFBQBAQBAQB3QgIAdENlXQUBAQBQUAQEAQEAUFAEBAEBIFwERACEi6CMl4QEAQEAUFAEBAEBAFBQBAQBHQjIAREN1TSURAQBAQBQUAQEAQEAUFAEBAEwkVACEi4CMp4QUAQEAQEAUFAEBAEBAFBQBDQjYAQEN1QSUdBQBAQBAQBqxDgHSelS5dG9+7drRIpcgQBQUAQEAQ8goAQEI8YStQUBAQBQSAlISAEJCVZU9YiCAgCgoAxBISAGMNLegsCgoAgIAhYgIAQEAtAFBGCgCAgCHgUARbhuowAAAUJSURBVCEgHjWcqC0ICAKCgJUIkBDcf//9OHnyJNavX4+cOXOiU6dOaNiw4U3TjBo1Ctu3b8enn36a+NqpU6dQqVIlzJ07Fw888AD69eun5Jw5cwa33nqrktO5c2dER0erMb4E5OjRo6hatSpWrFiBe+65R72+ceNGtGjRArt370ZMTIz62xdffIGpU6fi2LFjuP3229GxY0fUrVvXShhEliAgCAgCgoADCAgBcQBkmUIQEAQEgUhHgISApOL999/H448/rsgDCQhJRpkyZZKof+TIEdSsWVMRhrx586rXJk+ejKVLl2LhwoXq9/nz56Ny5crIkSOHktu+fXtV7/HMM8+YIiALFizA2LFjMW7cOEVwfvzxR7Rr107NW7Zs2UiHV/QTBAQBQUAQ8EFACIg8DoKAICAICAIqIpEtWzbl4GutW7duyJAhA4YOHXoTQoxOlCpVKrGInISEf3vuuecCojl48GCcOHEiUb7RCEi9evWUbI3AcJK+ffvi+vXrGDJkiFhQEBAEBAFBwEMICAHxkLFEVUFAEBAE7EKAhKBYsWLo1atX4hQjRozAnj17FHH4/fff1d9JBAYOHIhFixZh+PDh+O6777B161a0bdsW69atQ+bMmRUpmDhxIr7++muV0sXf4+LiULx4ccyePVvJMUpASpQogTRp0qgfrSUkJKjox5QpU+yCReQKAoKAICAI2ICAEBAbQBWRgoAgIAh4DQESAtZqMM1Ja0yZuuWWW/DOO+/ctBwSiooVK2LYsGFYvHgxoqKi1P/ZSDwY8fjoo49QpEgRRRr4O+s5Pvvss5sICOtEHnnkEXz55ZcoXLiwep0Ep2fPnok1IFWqVEHXrl0D1qR4DWvRVxAQBASB1I6AEJDU/gTI+gUBQUAQ+L+IxI4dOzBmzBhFLFgDwiLvadOmBa2xIKk4cOAAtm3bpshGuXLlFJYkGePHj1d1ILlz58amTZvwyiuvoECBAgEJCMewCJ1F7G+++SaOHz+uyAajL1oROmtRZs2aBUZlGKmJj4/H3r17VVE7f5cmCAgCgoAg4B0EhIB4x1aiqSAgCAgCtiHgfwoWi8c7dOiAp556KuicJAANGjRAvnz5sHz58sR+jI706dMHq1evVgShfPny6lStYBEQDtyyZQsGDBgAnohVqFChxFQv31OwvvrqK0WIWATPqErBggUVUdGIj23giGBBQBAQBAQBSxEQAmIpnCJMEBAEBAFvIiD3cnjTbqK1ICAICAJeREAIiBetJjoLAoKAIGAxAkJALAZUxAkCgoAgIAgERUAIiDwcgoAgIAgIAklOpRI4BAFBQBAQBAQBOxEQAmInuiJbEBAEBAFBQBAQBAQBQUAQEASSICAERB4IQUAQEAQEAUFAEBAEBAFBQBBwDAEhII5BLRMJAoKAICAICAKCgCAgCAgCgoAQEHkGBAFBQBAQBAQBQUAQEAQEAUHAMQSEgDgGtUwkCAgCgoAgIAgIAoKAICAICAJCQOQZEAQEAUFAEBAEBAFBQBAQBAQBxxAQAuIY1DKRICAICAKCgCAgCAgCgoAgIAgIAZFnQBAQBAQBQUAQEAQEAUFAEBAEHENACIhjUMtEgoAgIAgIAoKAICAICAKCgCAgBESeAUFAEBAEBAFBQBAQBAQBQUAQcAwBISCOQS0TCQKCgCAgCAgCgoAgIAgIAoKAEBB5BgQBQUAQEAQEAUFAEBAEBAFBwDEE/h8DMiDBNXhD9gAAAABJRU5ErkJggg==" width="999.999985098839">


.. parsed-literal::

    2. Reporting Generalized Performance:
    
    |                    | 98                      |
    |:-------------------|:------------------------|
    | clump_p1           | 1.0                     |
    | clump_r2           | 0.1                     |
    | clump_kb           | 200.0                   |
    | p_window_size      | 200.0                   |
    | p_slide_size       | 50.0                    |
    | p_LD_threshold     | 0.25                    |
    | pvalue             | 0.2976351441631313      |
    | numberofpca        | 6.0                     |
    | tempalpha          | 0.1                     |
    | l1weight           | 0.1                     |
    | burn_in            | 200.0                   |
    | num_iter           | 10.0                    |
    | temp_pvalue        | 0.006                   |
    | shrink_corr        | 0.1                     |
    | numberofvariants   | 173108.8                |
    | h2                 | 0.2032220489181616      |
    | Train_pure_prs     | -3.3108670615789036e-06 |
    | Train_null_model   | 0.2300103041419889      |
    | Train_best_model   | 0.342791259851683       |
    | Test_pure_prs      | -3.6300773897490755e-06 |
    | Test_null_model    | 0.11869244971792384     |
    | Test_best_model    | 0.29928645894907757     |
    | heritability_model | LDpred-2_full           |
    | sparse             | False                   |
    | allow_jump_sign    | True                    |
    | use_MLE            | False                   |
    | Difference         | 0.043504800902605445    |
    | Sum                | 0.6420777188007606      |
    3. Reporting the correlation of hyperparameters and the performance of 'Train_null_model', 'Train_pure_prs', 'Train_best_model', 'Test_pure_prs', 'Test_null_model', and 'Test_best_model':
    
    3. For string hyperparameters, we used one-hot encoding to find the correlation between string hyperparameters and 'Train_null_model', 'Train_pure_prs', 'Train_best_model', 'Test_pure_prs', 'Test_null_model', and 'Test_best_model'.
    3. We performed this analysis for those hyperparameters that have more than one unique value.
    


.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAA8AAAAKACAYAAABT+RFDAAAAAXNSR0IArs4c6QAAIABJREFUeF7snQd0FFUbht80SqgqoDRFkPJbKNKkidKbGHqTJh1BuhQBkV6ld6RXQVB6EQHpTQEBpSNSpIi00CH/+W7ckLIhu3Fn9wv7zjme/w+ZnXvnee9s5pnvzoxPSEhICLiQAAmQAAmQAAmQAAmQAAmQAAmQwDNOwIcC/IwnzN0jARIgARIgARIgARIgARIgARIwBCjAHAgkQAIkQAIkQAIkQAIkQAIkQAJeQYAC7BUxcydJgARIgARIgARIgARIgARIgAQowBwDJEACJEACJEACJEACJEACJEACXkGAAuwVMXMnSYAESIAESIAESIAESIAESIAEKMAcAyRAAiRAAiRAAiRAAiRAAiRAAl5BgALsFTFzJ0mABEiABEiABEiABEiABEiABCjAHAMkQAIkQAIkQAIkQAIkQAIkQAJeQYAC7BUxcydJgARIgARIgARIgARIgARIgAQowBwDJEACJEACJEACJEACJEACJEACXkGAAuwVMXMnSYAESIAESIAESIAESIAESIAEKMAcAyRAAiRAAiRAAiRAAiRAAiRAAl5BgALsFTFzJ0mABEiABEiABEiABEiABEiABCjAHAMkQAIkQAIkQAIkQAIkQAIkQAJeQYAC7BUxcydJgARIgARIgARIgARIgARIgAQowBwDJEACJEACJEACJEACJEACJEACXkGAAuwVMXMnSYAESIAESIAESIAESIAESIAEKMAcAyRAAiRAAiRAAiRAAiRAAiRAAl5BgALsFTFzJ0mABEiABEiABEiABEiABEiABCjAHAMkQAIkQAIkQAIkQAIkQAIkQAJeQYAC7BUxcydJgARIgARIgARIgARIgARIgAQowBwDJEACJEACJEACJEACJEACJEACXkGAAuwVMXMnSYAESIAESIAESIAESIAESIAEKMAcAyRAAiRAAiRAAiRAAiRAAiRAAl5BgALsFTFzJ0mABEiABEiABEiABEiABEiABCjAHAMkQAIkQAIkQAIkQAIkQAIkQAJeQYAC7BUxcydJgARIgARIgARIgARIgARIgAQowBwDJEACJEACJEACJEACJEACJEACXkGAAuwVMXMnSYAESIAESIAESIAESIAESIAEKMAcAyRAAiRAAiRAAiRAAiRAAiRAAl5BgALsFTFzJ0mABEiABEiABEiABEiABEiABCjAHAMkQAIkQAIkQAIkQAIkQAIkQAJeQYAC7BUxcydJgARIgARIgARIgARIgARIgAQowBwDJEACJEACJEACJEACJEACJEACXkGAAuwVMXMnSYAESIAESIAESIAESIAESIAEKMAcAyRAAiRAAiRAAiRAAiRAAiRAAl5BgALsFTFzJ0mABEiABEiABEiABEiABEiABCjAHAMkQAIkQAIkQAIkQAIkQAIkQAJeQYAC7BUxcydJgARIgARIgARIgARIgARIgAQowBwDJEACJEACJEACJEACJEACJEACXkGAAuwVMXMnSYAESIAESIAESIAESIAESIAEKMAcAyRAAiRAAiRAAiRAAiRAAiRAAl5BgALsFTFzJ0mABEiABEiABEiABEiABEiABCjAHAMkQAIkQAIkQAIkQAIkQAIkQAJeQYAC7BUxcydJgARIgARIgARIgARIgARIgAQowBwDJEACJEACJEACJEACJEACJEACXkGAAuwVMXMnSYAESIAESIAESIAESIAESIAEKMAcAyRAAiRAAiRAAiRAAiRAAiRAAl5BgALsFTFzJ0mABEiABEiABEiABEiABEiABCjAHAMkQAIkQAIkQAIkQAIkQAIkQAJeQYAC7BUxcydJgARIgARIgARIgARIgARIgAQowBwDJEACJEACJEACJEACJEACJEACXkGAAuwVMXMnSYAESIAESIAESIAESIAESIAEKMAcAyRAAiRAAiRAAiRAAiRAAiRAAl5BgALsFTFzJ0mABEiABEiABEiABEiABEiABCjAHAMkQAIkQAIkQAIkQAIkQAIkQAJeQYAC7BUxcydJgARIgARIgARIgARIgARIgAQowBwDJEACJEACJEACJEACJEACJEACXkGAAuwVMXMnSYAESIAESIAESIAESIAESIAEKMAcAyTgRQTOnj2L4sWLY+3atXjllVdiveeNGzdGrly58Mknn8R6G//1gw8ePECXLl2wefNmyP//6aefkCRJkv+6WX7eCwjI2J0wYQLy588f7d4WK1YMLVq0QLVq1TxGxJF+eqxzHmr48OHD6NatG/744w+8//77+OqrrzzUEzZLAiRAAiQQVwlQgONqcux3nCFw5MgRc7K9e/duBAcH47nnnkPOnDnRqFEjvPHGG27dD2cF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" width="1199.9999821186068">



.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <img src="data:image/png;base64,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PDEWdBJgPOmkybqcxBPpJZYABazE8ufVSYAESIAESIAESIAESIAESIAESIAESIAE4hCggMUQIQESIAESIAESIAESIAESIAESIAESIAES8DUBCli+bh4aRwIkQAIkQAIkQAIkQAIkQAIkQAIkQAIkQAGLMUACJEACJEACJEACJEACJEACJEACJEACJOBrAhSwfN08NI4ESIAESIAESIAESIAESIAESIAESIAESIACFmOABEiABEiABEiABEiABEiABEiABEiABEjA1wQoYPm6eWgcCZAACZAACZAACZAACZAACZAACZAACZAABSzGAAmQAAmQAAmQAAmQAAmQAAmQAAmQAAmQgK8JUMDydfPQOBIgARIgARIgARIgARIgARIgARIgARIgAQpYv8XAiRMn5Nlnn5X//Oc/UqBAAXn88celTZs2WUbI2rVr5YUXXpC9e/dK2bJlpX///lKhQoVARZMVf6+77jrJly+fXHLJJcrHqlWryuuvvx4Yf6dMmSKzZ8+W7du3S+PGjWX48OHZ2q6rba3w1XXNRDaIFX8ZT8HqK2LjKjU1VRCzGzdulJSUlESGna1ro+969913Zd++ffK73/1O7rnnHvn73/8uOXPmtFWfl4Ws3Gde2uX1tYIQg0GIs0TFk9nrfv755zJy5EjZvHmzCrGbbrpJ+vTpI6VLl1b/v2bNGnnkkUfU+MRI7du3lw4dOngdkhmuZ9Y/jCG7d+8uu3fvll9++UWNJ3v27CnVqlWL1Ifxyz//+U85efKk1KlTR1566SXVbzFlT8As/59++knxRnzt379fxo8fL7fddlukYowZBw0apH7/7rvvMj3zEJtjx46VPHnyRMqgjuj2c9JOuvyYM2eOII527dql7pWGDRtKr1691LsOEjggrubPn6+egy1btlRxaYz5nfiAsl75YaU9zNoUrw+Kx9YpO7PldfmzaNEiGTZsmBw5ckTFQvXq1dW7cdGiRc2aoiWfLn+ijXnmmWcE7YU4L1OmjBY7WUniCFDA+o09HmKnTp2SIUOGqAcZxKuBAwdKvXr1MrTO8ePHlQjSt29fadq0qUydOlUmTZokCxcuzPAQS1yTmruyWX9RGwSHIN/wH374oeTIkUNWrlwpaL/sBCydbWuWr85rmmt5d3KZ9ZfxFLy+Ijpi8PCfNWuWrFu3LrAC1rhx46RWrVrqowMGaX/7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" width="1499.9999776482584">


