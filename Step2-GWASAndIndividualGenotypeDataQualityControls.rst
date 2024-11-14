Step2-GWASAndIndividualGenotypeDataQualityControls
==================================================

Note on GWAS File
-----------------

1. It is important to note that different polygenic risk tools accept
   Genome-Wide Association summary statistic files in various formats.
2. The information we have encompasses all the columns required by PRS
   tools (where some fields are missing, we will highlight it).
3. For each individual tool, we will process this GWAS file so that it
   can be consumed by a specific PRS tool for prediction. The specific
   number of fields and the names of those fields are highlighted in
   each tool’s documentation.

**Fields in the GWAS File:**

These fields are also explained in `this
link <https://choishingwan.github.io/PRS-Tutorial/base/>`__.

+---+---+------+---+---+---+-----+-----+----------+---------+----------+
| C | B | SNP  | A | A | N | SE  | P   | OR       | INFO    | MAF      |
| H | P |      | 1 | 2 |   |     |     |          |         |          |
| R |   |      |   |   |   |     |     |          |         |          |
+===+===+======+===+===+===+=====+=====+==========+=========+==========+
| 1 | 7 | r    | A | G | 3 | 0   | 0.  | 0        | 0.8     | 0        |
|   | 5 | s313 |   |   | 8 | .00 | 483 | .9978869 | 9055794 | .3693895 |
|   | 6 | 1962 |   |   | 8 | 301 | 171 | 15712657 | 1364774 | 92764921 |
|   | 6 |      |   |   | 0 | 666 |     |          |         |          |
|   | 0 |      |   |   | 2 |     |     |          |         |          |
|   | 4 |      |   |   | 8 |     |     |          |         |          |
+---+---+------+---+---+---+-----+-----+----------+---------+----------+
| 1 | 7 | rs   | A | G | 3 | 0   | 0.  | 1.000687 | 0.8     | 0        |
|   | 6 | 1256 |   |   | 8 | .00 | 834 | 31609353 | 9589351 | .3368457 |
|   | 8 | 2034 |   |   | 8 | 329 | 808 |          | 1351165 | 54096289 |
|   | 4 |      |   |   | 0 | 472 |     |          |         |          |
|   | 4 |      |   |   | 2 |     |     |          |         |          |
|   | 8 |      |   |   | 8 |     |     |          |         |          |
+---+---+------+---+---+---+-----+-----+----------+---------+----------+
| 1 | 7 | r    | G | A | 3 | 0   | 0   | 0        | 0.8     | 0        |
|   | 7 | s404 |   |   | 8 | .00 | .42 | .9976035 | 9750829 | .3773680 |
|   | 9 | 0617 |   |   | 8 | 303 | 897 | 56067569 | 0615237 | 10940814 |
|   | 3 |      |   |   | 0 | 344 |     |          |         |          |
|   | 2 |      |   |   | 2 |     |     |          |         |          |
|   | 2 |      |   |   | 8 |     |     |          |         |          |
+---+---+------+---+---+---+-----+-----+----------+---------+----------+

Genome-Wide Association Study (GWAS) (Base Data) - Quality Controls
-------------------------------------------------------------------

Please note that some tools have built-in quality control procedures for
both GWAS and genotype data.

Set the directory where the files are located
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   filedirec = "SampleData1"

In case you are running the jobs on HPC and want to use parallel
computing, you can replace it with ``sys.argv[1]``. All the files for a
specific phenotype will be produced in that specific directory, which in
this case is ``SampleData1``.

.. code:: ipython3

    import os
    import pandas as pd
    import subprocess
    
    # Set the directory where the files are located
    filedirec = "SampleData1"
    #filedirec = "asthma_19"
    #filedirec = "migraine_0"
    
    # Define file paths for different data files
    BED = filedirec + os.sep + filedirec
    BIM = filedirec + os.sep + filedirec+".bim"
    FAM = filedirec + os.sep + filedirec+".fam"
    COV = filedirec + os.sep + filedirec+".cov"
    Height = filedirec + os.sep + filedirec+".height"
    GWAS = filedirec + os.sep + filedirec+".gz"
    
    # Read GWAS data from a compressed file using pandas
    df = pd.read_csv(GWAS, compression="gzip", sep="\s+")
    
    # Display the initial number of rows in the dataframe
    print("Initial number of SNPs:", len(df))
    
    # Apply quality control steps: Filter SNPs based on Minor Allele Frequency (MAF) and Imputation Information Score (INFO)
    df = df.loc[(df['MAF'] > 0.01) & (df['INFO'] > 0.8)]
    
    # Display the number of rows after applying the filters
    print("Number of SNPs after quality control:", len(df))
    
     
    # Display the number of rows after removing duplicate SNPs
    print("SNPs in GWAS after removing duplicate SNPs:", len(df))
    
    # Remove ambiguous SNPs with complementary alleles (C/G or A/T) to avoid potential errors
    df = df[~((df['A1'] == 'A') & (df['A2'] == 'T') |
              (df['A1'] == 'T') & (df['A2'] == 'A') |
              (df['A1'] == 'G') & (df['A2'] == 'C') |
              (df['A1'] == 'C') & (df['A2'] == 'G'))]
    
    # Display the final number of SNPs after removing ambiguous SNPs
    print("Final number of SNPs after removing ambiguous SNPs:", len(df))
    
    # Save the data.
    df.to_csv(GWAS,compression="gzip",sep="\t",index=None)
    
    df = pd.read_csv(GWAS,compression= "gzip",sep="\s+")
    print(len(df))
    print(df.head().to_markdown())
    


.. parsed-literal::

    Initial number of SNPs: 499617
    Number of SNPs after quality control: 499617
    SNPs in GWAS after removing duplicate SNPs: 499617
    Final number of SNPs after removing ambiguous SNPs: 499617
    499617
    |    |   CHR |     BP | SNP        | A1   | A2   |      N |         SE |        P |       OR |     INFO |      MAF |
    |---:|------:|-------:|:-----------|:-----|:-----|-------:|-----------:|---------:|---------:|---------:|---------:|
    |  0 |     1 | 756604 | rs3131962  | A    | G    | 388028 | 0.00301666 | 0.483171 | 0.997887 | 0.890558 | 0.36939  |
    |  1 |     1 | 768448 | rs12562034 | A    | G    | 388028 | 0.00329472 | 0.834808 | 1.00069  | 0.895894 | 0.336846 |
    |  2 |     1 | 779322 | rs4040617  | G    | A    | 388028 | 0.00303344 | 0.42897  | 0.997604 | 0.897508 | 0.377368 |
    |  3 |     1 | 801536 | rs79373928 | G    | T    | 388028 | 0.00841324 | 0.808999 | 1.00204  | 0.908963 | 0.483212 |
    |  4 |     1 | 808631 | rs11240779 | G    | A    | 388028 | 0.00242821 | 0.590265 | 1.00131  | 0.893213 | 0.45041  |
    

Match Variants Between GWAS and Individual Genotype Data
--------------------------------------------------------

If RSID is present in the GWAS, the following step can be skipped.

Steps for Handling RSIDs in GWAS and Genotype Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. If RSIDs are not present for SNPs, put ``X`` in the SNP column in the
   GWAS file.

2. Read the ``genotype.bim`` file and extract the RSIDs from the
   genotype data.

3. If RSIDs are not present in the genotype data, use HapMap3 or another
   reference panel to obtain the RSIDs.

Some PRS tools use different criteria to create unique variants and
match them between GWAS and individual genotype data:

-  **CHR:BP:A1:A2**: Some PRS tools use this format to define a unique
   variant.
-  **RSID**: Some PRS tools use RSID/SNP to define a unique variant.
-  **CHR:BP**: Some PRS tools use this format to define a unique
   variant.

We have highlighted which criteria are necessary for each tool.

.. code:: ipython3

    
    
    bimfile = pd.read_csv(BIM, sep="\s+", header=None)
    print("Columns of BIM file:")
    print(bimfile.columns)
    print("First 10 rows of BIM file:")
    
    
    print("Removing SNPs for which even a single row does not contain the required value:", len(df))
    
    
    # If RSID's are not present for SNPs, put X in the SNP column in the GWAS file.
    # Read the genotype.bim file, and extract the RSID from the genotype data.
    # If RSID are not present in the genotype data, use HapMap3 or other reference panel to get the RSIDs.
    
    if (df['SNP'] == 'X').all():
        print("RSIDs are missing!")
        bimfile = pd.read_csv(filedirec+os.sep+filedirec+".bim", sep="\s+", header=None)
        
        # create a unique variant using CHR:BP:A1:A2.
        
        bimfile["match"] = bimfile[0].astype(str)+"_"+bimfile[3].astype(str)+"_"+bimfile[4].astype(str)+"_"+bimfile[5].astype(str)
        df["match"] = df["CHR"].astype(str)+"_"+df["BP"].astype(str)+"_"+df["A1"].astype(str)+"_"+df["A2"].astype(str)
        
    
      
    
        df.drop_duplicates(subset='match', inplace=True)
        bimfile.drop_duplicates(subset='match', inplace=True)
    
        df = df[df['match'].isin(bimfile['match'].values)]
        bimfile = bimfile[bimfile['match'].isin(df['match'].values)]
        df = df[df['match'].isin(bimfile['match'].values)]
        bimfile = bimfile[bimfile['match'].isin(df['match'].values)]
     
        
        df = df.sort_values(by='BP')
        bimfile = bimfile.sort_values(by=3)
        
        print(df.head())
        print(bimfile.head())
    
        df["SNP"] = bimfile[1].values
        print("match",len(df))
    
    
        df.drop_duplicates(subset='match', inplace=True)
        bimfile.drop_duplicates(subset='match', inplace=True)  
    
        print(len(df))
        print(len(bimfile))
        print(df.head())
        print(bimfile.head())
        
        del df["match"]
        # Just save the modified GWAS file.
        # If bim, file is modified, the genotype data will be considered as corupt by Plink.
        df.to_csv(GWAS,compression="gzip",sep="\t",index=None)   
        print("Total SNPs", len(df))
    
        pass
    else:
        df.drop_duplicates(subset='SNP', inplace=True)
        df.to_csv(GWAS,compression="gzip",sep="\t",index=None)
        print("RSID is present!")
        print("Total SNPs",len(df))
        pass
     


.. parsed-literal::

    Columns of BIM file:
    Index([0, 1, 2, 3, 4, 5], dtype='int64')
    First 10 rows of BIM file:
    Removing SNPs for which even a single row does not contain the required value: 499617
    RSID is present!
    Total SNPs 499617
    

Individual genotype data (Target Data) Processing
-------------------------------------------------

Ensure that the phenotype file, FAM file, and covariate file contain an
identical number of samples. Remove any missing samples based on your
data. Note that the extent of missingness in phenotypes and covariates
may vary.

**Note:** Plink needs to be installed or placed in the same directory as
this notebook.

`Download Plink <https://www.cog-genomics.org/plink/>`__

We recommend using Linux. In cases where Windows is required due to
package installation issues on Linux, we provide the following guidance:

1. For Windows, use ``plink``.
2. For Linux, use ``./plink``.

Remove people with missing Phenotype
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Modify the fam file, make bed file, and modify the covariates files as
well.

.. code:: ipython3

    # New files to be saved with QC suffix
    newfilename = filedirec + "_QC"
    
    # Read information from FAM file
    f = pd.read_csv(FAM, header=None, sep="\s+", names=["FID", "IID", "Father", "Mother", "Sex", "Phenotype"])
    print("FAM file contents:")
    print(f.head())
    print("Total number of people in FAM file:", len(f))
    
    # Append the Height phenotype values to FAM file
    # Height file is basically the phenotype file.
    h = pd.read_csv(Height, sep="\t")
    print("Phenotype information is available for:", len(h), "people")
    print(len(h))
    result = pd.merge(f, h, on=['FID', 'IID'])
    
    # Replace 'Phenotype' column with 'Height' and save to a new PeopleWithPhenotype.txt file
    # Ensure that the input Phenotype file has teh header Height.
    result["Phenotype"] = result["Height"].values
    del result["Height"]
    
    # Remove NA or missing in the phenotype column
    result = result.dropna(subset=["Phenotype"])
    
    
    print(result)
    result.to_csv(filedirec + os.sep + "PeopleWithPhenotype.txt", index=False, header=False, sep="\t")
    
    # Use plink to keep only the people with phenotype present
    plink_command = [
        './plink',
        '--bfile', filedirec + os.sep + filedirec,
        '--keep', filedirec + os.sep + "PeopleWithPhenotype.txt",
        '--make-bed',
        '--out', filedirec + os.sep + newfilename
    ]
    subprocess.run(plink_command)
    
    # Update the phenotype information in the new FAM file
    f = pd.read_csv(filedirec + os.sep + newfilename + ".fam", header=None, sep="\s+",
                    names=["FID", "IID", "Father", "Mother", "Sex", "Phenotype"])
    f["Phenotype"] = result["Phenotype"].values
    f.to_csv(filedirec + os.sep + newfilename + ".fam", index=False, header=False, sep="\t")
    
    # Update the covariate file as well
    covfile = filedirec + os.sep + filedirec + '.cov'
    covfile = pd.read_csv(covfile, sep="\s+")
    
    print("Covariate file contents:")
    print(covfile.head())
    print("Total number of people in Covariate file:", len(covfile))
    
    # Match the FID and IID from covariate and height file
    covfile = covfile[covfile['FID'].isin(f["FID"].values) & covfile['IID'].isin(f["IID"].values)]
    print("Covariate file contents after matching with FAM file:")
    print(covfile.head())
    print("Total number of people in Covariate file after matching:", len(covfile))
    covfile.to_csv(filedirec + os.sep + newfilename + ".cov", index=None, sep="\t")
    


.. parsed-literal::

    FAM file contents:
           FID      IID  Father  Mother  Sex  Phenotype
    0  HG00096  HG00096       0       0    1         -9
    1  HG00097  HG00097       0       0    2         -9
    2  HG00099  HG00099       0       0    2         -9
    3  HG00100  HG00100       0       0    2         -9
    4  HG00101  HG00101       0       0    1         -9
    Total number of people in FAM file: 503
    Phenotype information is available for: 475 people
    475
             FID      IID  Father  Mother  Sex   Phenotype
    0    HG00096  HG00096       0       0    1  169.132169
    1    HG00097  HG00097       0       0    2  171.256259
    2    HG00099  HG00099       0       0    2  171.534380
    3    HG00101  HG00101       0       0    1  169.850176
    4    HG00102  HG00102       0       0    2  172.788361
    ..       ...      ...     ...     ...  ...         ...
    470  NA20822  NA20822       0       0    2  170.405056
    471  NA20826  NA20826       0       0    2  168.523029
    472  NA20827  NA20827       0       0    1  170.975735
    473  NA20828  NA20828       0       0    2  170.222028
    474  NA20832  NA20832       0       0    2  169.431705
    
    [475 rows x 6 columns]
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/SampleData1_QC.log.
    Options in effect:
      --bfile SampleData1/SampleData1
      --keep SampleData1/PeopleWithPhenotype.txt
      --make-bed
      --out SampleData1/SampleData1_QC
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    503 people (240 males, 263 females) loaded from .fam.
    --keep: 475 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 475 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate in remaining samples is 0.999896.
    551892 variants and 475 people pass filters and QC.
    Note: No phenotypes present.
    --make-bed to SampleData1/SampleData1_QC.bed + SampleData1/SampleData1_QC.bim +
    SampleData1/SampleData1_QC.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    Covariate file contents:
           FID      IID  Sex
    0  HG00096  HG00096    1
    1  HG00097  HG00097    2
    2  HG00099  HG00099    2
    3  HG00100  HG00100    2
    4  HG00101  HG00101    1
    Total number of people in Covariate file: 503
    Covariate file contents after matching with FAM file:
           FID      IID  Sex
    0  HG00096  HG00096    1
    1  HG00097  HG00097    2
    2  HG00099  HG00099    2
    4  HG00101  HG00101    1
    5  HG00102  HG00102    2
    Total number of people in Covariate file after matching: 475
    

Split Data into Test and Train and Perform Quality Controls on Training Data
----------------------------------------------------------------------------

1. We adopt a cross-validation design to evaluate polygenic risk scores.
2. The base data is divided into two sets: training and test sets. The
   training set is used to find the best combination of parameters or
   hyperparameters offered by each tool, along with the summary
   statistic file. The data is split into 5 folds, and further
   processing is performed on the first fold. ``Fold_0``
3. Data is divided into training and test sets, and hyperparameter
   optimization is performed on the training data. Operations such as
   clumping or pruning are applied to the training data, and the same
   remaining SNPs should be extracted from the test set rather than
   separately using pruning or clumping on test data.
4. Regarding hyperparameters from individual tools, all of those should
   be applied to the training set, and the performance should be
   measured for both the null model and the complete model (including
   covariates and polygenic risk scores). Since it is a continuous
   phenotype, Explained Variance is considered to assess the performance
   of the PRS model. For binary phenotypes, we considered AUC to
   evaluate the PRS performance.

**Note:** We will divide the data into training and test sets, perform
quality controls on the training data. The data for each fold will be
saved separately for further processing. .

Pruning
~~~~~~~

Pruning is an integral part of the analysis, but it has been skipped as
a quality control on the training set for the following reasons. In our
initial analysis, we observed that pruning and clumping can affect the
performance of the polygenic risk score model. Rather than performing it
at this stage, we consider it as one of the hyperparameters and it will
be performed at a later stage. If we perform pruning at this stage, the
process that passed the pruning step would limit us to use other values
for pruning at the latest stage. However, even at the latest stage,
pruning is only performed on the training set.

Cross-validation Designs
~~~~~~~~~~~~~~~~~~~~~~~~

There are multiple ways of doing cross-validation design, and one of
them is to use all the cohorts except one for the training and then use
the last cohort as the validation set. We performed quality controls and
hyperparameter optimization on the training set and found the best
combination across all folds, reporting the performance of the best
hyperparameter combination on the test set.

A Simple Analysis
~~~~~~~~~~~~~~~~~

If you have a separate GWAS file and our training data that you will be
using to optimize the hyperparameter on a subset without
cross-validation, then follow the original tutorial presented in the
first cell `Shing Wan Choi’s PRS
Tutorial <https://choishingwan.github.io/PRS-Tutorial/base/>`__.

**Quality Controls considered for Training Data:** 1. GWAS studies,
e.g., removing SNPs with low genotyping rate, low minor allele
frequency, out of Hardy-Weinberg Equilibrium, removing individuals with
low genotyping rate. 2. Pruning was skipped. 3. Heterozygosity check. 4.
Sex chromosomes. 5. Relatedness.

**R Script - Module1.R**

This file contains the code presented in this
`tutorial <https://choishingwan.github.io/PRS-Tutorial/base/>`__ to
assist in performing quality controls on the training data. Kindly,
follow their instructions for better understanding, as quality controls
are not the main focus of this research.

.. code:: ipython3

    from IPython.display import FileLink
    
    # R file used in to execute the following code.
    FileLink('Module1.R')




.. raw:: html

    <a href='Module1.R' target='_blank'>Module1.R</a><br>



**R Script - Module1.R**

 
This file contains the code presented in this [tutorial](https://choishingwan.github.io/PRS-Tutorial/base/) to assist in performing quality controls on the training data. Kindly, follow their instructions for better understanding, as quality controls are not the main focus of this research.

args <- commandArgs(trailingOnly = TRUE)
print(args)
# Argument one is going to be the directory.
# Argument two is going to be the file name.
# Argument three is going to be the output file name.
# Argument four is going to be the specific function to be called.

if (args[4]=="1"){

  result <-paste(".",args[1],paste(args[3],toString(".het"), sep = ""),sep="//")
  dat <- read.table(result, header=T) # Read in the EUR.het file, specify it has header
  m <- mean(dat$F) # Calculate the mean  
  s <- sd(dat$F) # Calculate the SD
  valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
  result <-paste("./",args[1],paste(args[2],toString(".valid.sample"), sep = ""),sep="//")
  write.table(valid[,c(1,2)], result, quote=F, row.names=F) # print FID and IID for valid samples
  result <-paste("./",args[1],paste(args[2],toString(".bim"), sep = ""),sep="//")
  bim <- read.table(result)
  
  colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
  
  # Read in QCed SNPs
  result <-paste("./",args[1],paste(args[3],toString(".snplist"), sep = ""),sep="//")
  
  qc <- read.table(result, header = F, stringsAsFactors = F)
  # Read in the GWAS data
  path_parts <- strsplit(args[1], "/|\\\\")[[1]]
  print(path_parts)
  print(path_parts[1])
  result <-paste("./",path_parts[1],paste(path_parts[1],toString(".gz"), sep = ""),sep="//")
  #result <-paste("./",args[1],"train",toString("train.assoc.fisher"),sep="//")
  
  height <-read.table(gzfile(result),
               header = T,
               stringsAsFactors = F, 
               sep="\t")
  # Change all alleles to upper case for easy comparison
  

  height$A1 <- toupper(height$A1)
  height$A2 <- toupper(height$A2)
  bim$B.A1 <- toupper(bim$B.A1)
  bim$B.A2 <- toupper(bim$B.A2)
  info <- merge(bim, height, by = c("SNP", "CHR", "BP"))
  # Filter QCed SNPs
  print(length(info))

  info <- info[info$SNP %in% qc$V1,]
  print(length(info))
  # Function for finding the complementary allele
  
  complement <- function(x) {
    switch (
      x,
      "A" = "T",
      "C" = "G",
      "T" = "A",
      "G" = "C",
      return(NA)
    )
  }
  
  # Get SNPs that have the same alleles across base and target
  info.match <- subset(info, A1 == B.A1 & A2 == B.A2)
  print(length(info.match))
   
  # Identify SNPs that are complementary between base and target
  info$C.A1 <- sapply(info$B.A1, complement)
  info$C.A2 <- sapply(info$B.A2, complement)
  info.complement <- subset(info, A1 == C.A1 & A2 == C.A2)
  # Update the complementary alleles in the bim file
  # This allow us to match the allele in subsequent analysis
  
  complement.snps <- bim$SNP %in% info.complement$SNP
  bim[complement.snps,]$B.A1 <-
    sapply(bim[complement.snps,]$B.A1, complement)
  bim[complement.snps,]$B.A2 <-
    sapply(bim[complement.snps,]$B.A2, complement)
  
  # identify SNPs that need recoding
  info.recode <- subset(info, A1 == B.A2 & A2 == B.A1)
  # Update the recode SNPs
  recode.snps <- bim$SNP %in% info.recode$SNP
  tmp <- bim[recode.snps,]$B.A1
  bim[recode.snps,]$B.A1 <- bim[recode.snps,]$B.A2
  bim[recode.snps,]$B.A2 <- tmp
  
  # identify SNPs that need recoding & complement
  info.crecode <- subset(info, A1 == C.A2 & A2 == C.A1)
  # Update the recode + strand flip SNPs
  com.snps <- bim$SNP %in% info.crecode$SNP
  tmp <- bim[com.snps,]$B.A1
  bim[com.snps,]$B.A1 <- as.character(sapply(bim[com.snps,]$B.A2, complement))
  bim[com.snps,]$B.A2 <- as.character(sapply(tmp, complement))
  result <-paste("./",args[1],paste(args[2],toString(".a1"), sep = ""),sep="//")
  
  # Output updated bim file
  write.table(
    bim[,c("SNP", "B.A1")],
    result,
    quote = F,
    row.names = F,
    col.names = F,
    sep="\t"
  )
  mismatch <-
    bim$SNP[!(bim$SNP %in% info.match$SNP |
                bim$SNP %in% info.complement$SNP | 
                bim$SNP %in% info.recode$SNP |
                bim$SNP %in% info.crecode$SNP)]
  result <-paste("./",args[1],paste(args[2],toString(".mismatch"), sep = ""),sep="//")
  
  write.table(
    mismatch,
    result,
    quote = F,
    row.names = F,
    col.names = F
  )
  
  
}
if (args[4]=="2"){
  result <-paste("./",args[1],paste(args[2],toString(".valid.sample"), sep = ""),sep="//")
  valid <- read.table(result, header=T)
  result <-paste("./",args[1],paste(args[3],toString(".sexcheck"), sep = ""),sep="//")
  dat <- read.table(result, header=T)
  valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
  result <-paste("./",args[1],paste(args[3],toString(".valid"), sep = ""),sep="//")
  print(result)
  write.table(valid[,c("FID", "IID")], result, row.names=F, col.names=F, sep="\t", quote=F) 
}





.. code:: ipython3

    
    import os
    import pandas as pd
    from sklearn.model_selection import StratifiedKFold
    import subprocess
    from sklearn.model_selection import KFold, cross_val_score
    
    # Step 1: Read the Fam file
    
    input_file_path = filedirec+os.sep+newfilename+'.fam'
    df = pd.read_csv(input_file_path,sep="\s+",header=None)
    
    # Step 2: Create 5 directories for storing fold information
    output_directory_base = filedirec
    os.makedirs(output_directory_base, exist_ok=True)
    
    # Step 3: Split the data into 5 folds using cross validation 
    fold_column = 5  # fifth column contains phenotypes
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    phenotype_col = 5
    
    # The following code is for binary phenotype.
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    column_values = df[phenotype_col].unique()
    
    if set(column_values) == {1, 2}:
        print("The column contains only 0 and 1.")
        #exit(0)
        for fold_id, (train_index, test_index) in enumerate(skf.split(df, df[phenotype_col])):
            fold_directory = os.path.join(output_directory_base, f'Fold_{fold_id}')
    
            train_file_name = "train_data"
            test_file_name = "test_data"
    
            new_train_file_name = "train_data.QC"
            new_test_file_name = "test_data.QC"
    
            os.makedirs(fold_directory, exist_ok=True)
            #"""
            # Save train and test data to separate CSV files
            train_data = df.iloc[train_index]
            test_data = df.iloc[test_index]
            
            train_data.to_csv(os.path.join(fold_directory, 'train_data.fam'),sep="\t",header=False,index=False)
            test_data.to_csv(os.path.join(fold_directory, 'test_data.fam'),sep="\t",header=False,index=False)
            #print(train_data)
            #exit(0)
            
            #exit(0)
            # Step 4: Use PLINK to extract test and train samples for each fold
            plink_train_command = [
                './plink',
                '--bfile', filedirec+os.sep+newfilename,
                '--keep', os.path.join(fold_directory, train_file_name+'.fam'),
                '--make-bed',
                '--out', os.path.join(fold_directory, train_file_name)
            ]
    
            plink_test_command = [
                './plink',
                '--bfile', filedirec+os.sep+newfilename,
                '--keep', os.path.join(fold_directory, test_file_name+'.fam'),
                '--make-bed',
                '--out', os.path.join(fold_directory, test_file_name)
            ]
    
            subprocess.run(plink_train_command)
            subprocess.run(plink_test_command)
            
            covfile = filedirec+os.sep+newfilename+'.cov'
            covfile = pd.read_csv(covfile)
    
            cov_train_data = covfile.iloc[train_index]
            cov_test_data = covfile.iloc[test_index]
    
            cov_train_data.to_csv(os.path.join(fold_directory, train_file_name+'.cov'),sep=",",index=False)
            cov_test_data.to_csv(os.path.join(fold_directory, test_file_name+'.cov'),sep=",",index=False)
    
            #exit(0)
            ### perform Quality controls on the training data only.
            plink_command_1 = [
                './plink',
                '--bfile', os.path.join(fold_directory, train_file_name),
                '--maf', '0.01',
                '--hwe', '1e-6',
                '--geno', '0.1',
                '--mind', '0.1',
                '--write-snplist',
                '--make-just-fam',
                '--out', os.path.join(fold_directory, new_train_file_name)
            ]
        
            subprocess.run(plink_command_1)
            
            # Command 2
            # Perform pruning. Skip it.
            plink_command_2 = [
                './plink',
                '--bfile', os.path.join(fold_directory, train_file_name),
                '--keep', os.path.join(fold_directory, new_train_file_name+'.fam'),
                '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
                '--indep-pairwise', '200', '50', '0.25',
                '--out', os.path.join(fold_directory, new_train_file_name)
            ]
        
            #subprocess.run(plink_command_2)
            
            # Command 3
            plink_command_3 = [
                './plink',
                '--bfile', os.path.join(fold_directory, train_file_name),
                #'--extract', os.path.join(fold_directory, new_train_file_name+'.prune.in'),
                '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
                '--keep', os.path.join(fold_directory, new_train_file_name+'.fam'),
                '--het',
                '--out', os.path.join(fold_directory, new_train_file_name)
            ]
        
            subprocess.run(plink_command_3)
            # Invoked R functions.
            
            os.system("Rscript Module1.R "+os.path.join(fold_directory)+"  "+train_file_name+" "+new_train_file_name+ " "+"1")
            print("Rscript Module1.R "+os.path.join(fold_directory)+"  "+train_file_name+" "+new_train_file_name+ " "+"1")
             
            
            # Code for sex check: The sample data have 22 chromosomes, so this operation is skipped.
            # If the Chromosome X is not available, then do not execute the following two commands. 
            plink_command = [
            './plink',
            '--bfile', os.path.join(fold_directory, 'train_data'),
            #'--extract', os.path.join(fold_directory, 'train_data.QC.prune.in'),
            '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
            '--keep', os.path.join(fold_directory, 'train_data.valid.sample'),
            '--check-sex',
            '--out', os.path.join(fold_directory, 'train_data.QC')
            ]
            
            # Invoke the PLINK command using subprocess for sex check
            #subprocess.run(plink_command)
            #os.system("Rscript Module1.R "+os.path.join(fold_directory)+"  "+train_file_name+" "+new_train_file_name+ " "+"2")
            
            #"""
        
            plink_command_1 = [
            './plink',
            '--bfile', os.path.join(fold_directory, train_file_name),
            #'--extract', os.path.join(fold_directory, new_train_file_name+'.prune.in'),
            # If you use sex check, then use the following line. Otherwise, it is not required.
            '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
    
            # Uncomment the following line if the sex Chromosome is available.
            # Kindly not, if you perform Sex check, and some people are removed, you have to
            # remove people from covariate file and Phenotype file.
                
            #'--keep', os.path.join(fold_directory, 'train_data.QC.valid'),
            '--rel-cutoff', '0.125',
            '--out', os.path.join(fold_directory, new_train_file_name)
            ]
            
            subprocess.run(plink_command_1)
             
             
            plink_command_2 = [
                './plink',
                '--bfile', os.path.join(fold_directory, train_file_name),
                '--make-bed',
                '--keep', os.path.join(fold_directory, new_train_file_name+'.rel.id'),
                '--out', os.path.join(fold_directory, new_train_file_name),
                '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
                '--exclude', os.path.join(fold_directory, train_file_name+'.mismatch'),
                '--a1-allele', os.path.join(fold_directory, train_file_name+'.a1')
            ]
        
            subprocess.run(plink_command_2)
             
    else:
        print("The column does not contain only 0 and 1.")
        for fold_id, (train_index, test_index) in enumerate(kf.split(df, df[phenotype_col])):
            fold_directory = os.path.join(output_directory_base, f'Fold_{fold_id}')
    
            train_file_name = "train_data"
            test_file_name = "test_data"
    
            new_train_file_name = "train_data.QC"
            new_test_file_name = "test_data.QC"
    
            os.makedirs(fold_directory, exist_ok=True)
            #"""
            # Save train and test data to separate CSV files
            train_data = df.iloc[train_index]
            test_data = df.iloc[test_index]
            
            train_data.to_csv(os.path.join(fold_directory, 'train_data.fam'),sep="\t",header=False,index=False)
            test_data.to_csv(os.path.join(fold_directory, 'test_data.fam'),sep="\t",header=False,index=False)
            #print(train_data)
            #exit(0)
            
            #exit(0)
            # Step 4: Use PLINK to extract test and train samples for each fold
            plink_train_command = [
                './plink',
                '--bfile', filedirec+os.sep+newfilename,
                '--keep', os.path.join(fold_directory, train_file_name+'.fam'),
                '--make-bed',
                '--out', os.path.join(fold_directory, train_file_name)
            ]
    
            plink_test_command = [
                './plink',
                '--bfile', filedirec+os.sep+newfilename,
                '--keep', os.path.join(fold_directory, test_file_name+'.fam'),
                '--make-bed',
                '--out', os.path.join(fold_directory, test_file_name)
            ]
    
            subprocess.run(plink_train_command)
            subprocess.run(plink_test_command)
            
            covfile = filedirec+os.sep+newfilename+'.cov'
            covfile = pd.read_csv(covfile)
    
            cov_train_data = covfile.iloc[train_index]
            cov_test_data = covfile.iloc[test_index]
    
            cov_train_data.to_csv(os.path.join(fold_directory, train_file_name+'.cov'),sep=",",index=False)
            cov_test_data.to_csv(os.path.join(fold_directory, test_file_name+'.cov'),sep=",",index=False)
    
            #exit(0)
            ### perform Quality controls on the training data only.
            plink_command_1 = [
                './plink',
                '--bfile', os.path.join(fold_directory, train_file_name),
                '--maf', '0.01',
                '--hwe', '1e-6',
                '--geno', '0.1',
                '--mind', '0.1',
                '--write-snplist',
                '--make-just-fam',
                '--out', os.path.join(fold_directory, new_train_file_name)
            ]
        
            subprocess.run(plink_command_1)
            
            # Command 2
            plink_command_2 = [
                './plink',
                '--bfile', os.path.join(fold_directory, train_file_name),
                '--keep', os.path.join(fold_directory, new_train_file_name+'.fam'),
                '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
                '--indep-pairwise', '200', '50', '0.25',
                '--out', os.path.join(fold_directory, new_train_file_name)
            ]
        
            #subprocess.run(plink_command_2)
            
            # Command 3
            plink_command_3 = [
                './plink',
                '--bfile', os.path.join(fold_directory, train_file_name),
                #'--extract', os.path.join(fold_directory, new_train_file_name+'.prune.in'),
                '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
                '--keep', os.path.join(fold_directory, new_train_file_name+'.fam'),
                '--het',
                '--out', os.path.join(fold_directory, new_train_file_name)
            ]
        
            subprocess.run(plink_command_3)
            # Invoked R functions.
            
            os.system("Rscript Module1.R "+os.path.join(fold_directory)+"  "+train_file_name+" "+new_train_file_name+ " "+"1")
            print("Rscript Module1.R "+os.path.join(fold_directory)+"  "+train_file_name+" "+new_train_file_name+ " "+"1")
            
            
            
            # Code for sex check: The sample data have 22 chromosomes, so this operation is skipped.
            # If the Chromosome X is not available, then do not execute the following two commands. 
            #"""
            plink_command = [
            './plink',
            '--bfile', os.path.join(fold_directory, 'train_data'),
            #'--extract', os.path.join(fold_directory, 'train_data.QC.prune.in'),
            '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
            '--keep', os.path.join(fold_directory, 'train_data.valid.sample'),
            '--check-sex',
            '--out', os.path.join(fold_directory, 'train_data.QC')
            ]
            
            # Invoke the PLINK command using subprocess for sex check
             
            #subprocess.run(plink_command)
            #os.system("Rscript Module1.R "+os.path.join(fold_directory)+"  "+train_file_name+" "+new_train_file_name+ " "+"2")
            
            #"""
        
            plink_command_1 = [
            './plink',
            '--bfile', os.path.join(fold_directory, train_file_name),
            #'--extract', os.path.join(fold_directory, new_train_file_name+'.prune.in'),
            # If you use sex check then use the following line. Otherwise it is not required.
            '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
            #'--keep', os.path.join(fold_directory, 'train_data.QC.valid'),
            '--rel-cutoff', '0.125',
            '--out', os.path.join(fold_directory, new_train_file_name)
            ]
            
            subprocess.run(plink_command_1)
            #exit(0)
            
            #exit(0)
            # Command 2
            plink_command_2 = [
                './plink',
                '--bfile', os.path.join(fold_directory, train_file_name),
                '--make-bed',
                '--keep', os.path.join(fold_directory, new_train_file_name+'.rel.id'),
                '--out', os.path.join(fold_directory, new_train_file_name),
                '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
                '--exclude', os.path.join(fold_directory, train_file_name+'.mismatch'),
                '--a1-allele', os.path.join(fold_directory, train_file_name+'.a1')
            ]
        
            subprocess.run(plink_command_2)
            #exit(0)
     
    
    
    
     
    
    


.. parsed-literal::

    The column does not contain only 0 and 1.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.log.
    Options in effect:
      --bfile SampleData1/SampleData1_QC
      --keep SampleData1/Fold_0/train_data.fam
      --make-bed
      --out SampleData1/Fold_0/train_data
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    475 people (227 males, 248 females) loaded from .fam.
    475 phenotype values loaded from .fam.
    --keep: 380 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate in remaining samples is 0.999896.
    551892 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.bed +
    SampleData1/Fold_0/train_data.bim + SampleData1/Fold_0/train_data.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.log.
    Options in effect:
      --bfile SampleData1/SampleData1_QC
      --keep SampleData1/Fold_0/test_data.fam
      --make-bed
      --out SampleData1/Fold_0/test_data
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    475 people (227 males, 248 females) loaded from .fam.
    475 phenotype values loaded from .fam.
    --keep: 95 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate in remaining samples is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.bed +
    SampleData1/Fold_0/test_data.bim + SampleData1/Fold_0/test_data.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data
      --geno 0.1
      --hwe 1e-6
      --maf 0.01
      --make-just-fam
      --mind 0.1
      --out SampleData1/Fold_0/train_data.QC
      --write-snplist
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    0 people removed due to missing genotype data (--mind).
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    1 variant removed due to missing genotype data (--geno).
    --hwe: 822 variants removed due to Hardy-Weinberg exact test.
    8415 variants removed due to minor allele threshold(s)
    (--maf/--max-maf/--mac/--max-mac).
    542654 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    List of variant IDs written to SampleData1/Fold_0/train_data.QC.snplist .
    --make-just-fam to SampleData1/Fold_0/train_data.QC.fam ... done.
    

.. parsed-literal::

    Warning: --hwe observation counts vary by more than 10%, due to the X
    chromosome.  You may want to use a less stringent --hwe p-value threshold for X
    chromosome variants.
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data
      --extract SampleData1/Fold_0/train_data.QC.snplist
      --het
      --keep SampleData1/Fold_0/train_data.QC.fam
      --out SampleData1/Fold_0/train_data.QC
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 542654 variants remaining.
    --keep: 380 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999898.
    542654 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --het: 521371 variants scanned, report written to
    SampleData1/Fold_0/train_data.QC.het .
    [1] "SampleData1/Fold_0" "train_data"         "train_data.QC"     
    [4] "1"                 
    [1] "SampleData1" "Fold_0"     
    [1] "SampleData1"
    [1] 14
    [1] 14
    [1] 14
    Rscript Module1.R SampleData1/Fold_0  train_data train_data.QC 1
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data
      --extract SampleData1/Fold_0/train_data.QC.snplist
      --out SampleData1/Fold_0/train_data.QC
      --rel-cutoff 0.125
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 542654 variants remaining.
    Using up to 8 threads (change this with --threads).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999898.
    542654 variants and 380 people pass filters and QC (before --rel-cutoff).
    Phenotype data is quantitative.
    Excluding 21283 variants on non-autosomes from relationship matrix calc.
    Relationship matrix calculation complete.
    0 people excluded by --rel-cutoff.
    Remaining sample IDs written to SampleData1/Fold_0/train_data.QC.rel.id .
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.log.
    Options in effect:
      --a1-allele SampleData1/Fold_0/train_data.a1
      --bfile SampleData1/Fold_0/train_data
      --exclude SampleData1/Fold_0/train_data.mismatch
      --extract SampleData1/Fold_0/train_data.QC.snplist
      --keep SampleData1/Fold_0/train_data.QC.rel.id
      --make-bed
      --out SampleData1/Fold_0/train_data.QC
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 542654 variants remaining.
    --exclude: 491952 variants remaining.
    --keep: 380 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999894.
    --a1-allele: 491952 assignments made.
    491952 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.bed +
    SampleData1/Fold_0/train_data.QC.bim + SampleData1/Fold_0/train_data.QC.fam ...
    101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_1/train_data.log.
    Options in effect:
      --bfile SampleData1/SampleData1_QC
      --keep SampleData1/Fold_1/train_data.fam
      --make-bed
      --out SampleData1/Fold_1/train_data
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    475 people (227 males, 248 females) loaded from .fam.
    475 phenotype values loaded from .fam.
    --keep: 380 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate in remaining samples is 0.999922.
    551892 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_1/train_data.bed +
    SampleData1/Fold_1/train_data.bim + SampleData1/Fold_1/train_data.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_1/test_data.log.
    Options in effect:
      --bfile SampleData1/SampleData1_QC
      --keep SampleData1/Fold_1/test_data.fam
      --make-bed
      --out SampleData1/Fold_1/test_data
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    475 people (227 males, 248 females) loaded from .fam.
    475 phenotype values loaded from .fam.
    --keep: 95 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate in remaining samples is 0.999794.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_1/test_data.bed +
    SampleData1/Fold_1/test_data.bim + SampleData1/Fold_1/test_data.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_1/train_data.QC.log.
    Options in effect:
      --bfile SampleData1/Fold_1/train_data
      --geno 0.1
      --hwe 1e-6
      --maf 0.01
      --make-just-fam
      --mind 0.1
      --out SampleData1/Fold_1/train_data.QC
      --write-snplist
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (178 males, 202 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    0 people removed due to missing genotype data (--mind).
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999922.
    1 variant removed due to missing genotype data (--geno).
    --hwe: 855 variants removed due to Hardy-Weinberg exact test.
    7923 variants removed due to minor allele threshold(s)
    (--maf/--max-maf/--mac/--max-mac).
    543113 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    List of variant IDs written to SampleData1/Fold_1/train_data.QC.snplist .
    --make-just-fam to SampleData1/Fold_1/train_data.QC.fam ... done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_1/train_data.QC.log.
    Options in effect:
      --bfile SampleData1/Fold_1/train_data
      --extract SampleData1/Fold_1/train_data.QC.snplist
      --het
      --keep SampleData1/Fold_1/train_data.QC.fam
      --out SampleData1/Fold_1/train_data.QC
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (178 males, 202 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 543113 variants remaining.
    --keep: 380 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697%

.. parsed-literal::

    Warning: --hwe observation counts vary by more than 10%, due to the X
    chromosome.  You may want to use a less stringent --hwe p-value threshold for X
    chromosome variants.
    

.. parsed-literal::

    989 done.
    Total genotyping rate is 0.999923.
    543113 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --het: 521800 variants scanned, report written to
    SampleData1/Fold_1/train_data.QC.het .
    [1] "SampleData1/Fold_1" "train_data"         "train_data.QC"     
    [4] "1"                 
    [1] "SampleData1" "Fold_1"     
    [1] "SampleData1"
    [1] 14
    [1] 14
    [1] 14
    Rscript Module1.R SampleData1/Fold_1  train_data train_data.QC 1
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_1/train_data.QC.log.
    Options in effect:
      --bfile SampleData1/Fold_1/train_data
      --extract SampleData1/Fold_1/train_data.QC.snplist
      --out SampleData1/Fold_1/train_data.QC
      --rel-cutoff 0.125
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (178 males, 202 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 543113 variants remaining.
    Using up to 8 threads (change this with --threads).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999923.
    543113 variants and 380 people pass filters and QC (before --rel-cutoff).
    Phenotype data is quantitative.
    Excluding 21313 variants on non-autosomes from relationship matrix calc.
    Relationship matrix calculation complete.
    0 people excluded by --rel-cutoff.
    Remaining sample IDs written to SampleData1/Fold_1/train_data.QC.rel.id .
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_1/train_data.QC.log.
    Options in effect:
      --a1-allele SampleData1/Fold_1/train_data.a1
      --bfile SampleData1/Fold_1/train_data
      --exclude SampleData1/Fold_1/train_data.mismatch
      --extract SampleData1/Fold_1/train_data.QC.snplist
      --keep SampleData1/Fold_1/train_data.QC.rel.id
      --make-bed
      --out SampleData1/Fold_1/train_data.QC
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (178 males, 202 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 543113 variants remaining.
    --exclude: 492382 variants remaining.
    --keep: 380 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.99992.
    --a1-allele: 492382 assignments made.
    492382 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_1/train_data.QC.bed +
    SampleData1/Fold_1/train_data.QC.bim + SampleData1/Fold_1/train_data.QC.fam ...
    101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_2/train_data.log.
    Options in effect:
      --bfile SampleData1/SampleData1_QC
      --keep SampleData1/Fold_2/train_data.fam
      --make-bed
      --out SampleData1/Fold_2/train_data
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    475 people (227 males, 248 females) loaded from .fam.
    475 phenotype values loaded from .fam.
    --keep: 380 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate in remaining samples is 0.999896.
    551892 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_2/train_data.bed +
    SampleData1/Fold_2/train_data.bim + SampleData1/Fold_2/train_data.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_2/test_data.log.
    Options in effect:
      --bfile SampleData1/SampleData1_QC
      --keep SampleData1/Fold_2/test_data.fam
      --make-bed
      --out SampleData1/Fold_2/test_data
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    475 people (227 males, 248 females) loaded from .fam.
    475 phenotype values loaded from .fam.
    --keep: 95 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate in remaining samples is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_2/test_data.bed +
    SampleData1/Fold_2/test_data.bim + SampleData1/Fold_2/test_data.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_2/train_data.QC.log.
    Options in effect:
      --bfile SampleData1/Fold_2/train_data
      --geno 0.1
      --hwe 1e-6
      --maf 0.01
      --make-just-fam
      --mind 0.1
      --out SampleData1/Fold_2/train_data.QC
      --write-snplist
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (180 males, 200 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    0 people removed due to missing genotype data (--mind).
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    1 variant removed due to missing genotype data (--geno).
    --hwe: 847 variants removed due to Hardy-Weinberg exact test.
    8536 variants removed due to minor allele threshold(s)
    (--maf/--max-maf/--mac/--max-mac).
    542508 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    List of variant IDs written to SampleData1/Fold_2/train_data.QC.snplist .
    --make-just-fam to SampleData1/Fold_2/train_data.QC.fam ... done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_2/train_data.QC.log.
    Options in effect:
      --bfile SampleData1/Fold_2/train_data
      --extract SampleData1/Fold_2/train_data.QC.snplist
      --het
      --keep SampleData1/Fold_2/train_data.QC.fam
      --out SampleData1/Fold_2/train_data.QC
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (180 males, 200 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 542508 variants remaining.
    --keep: 380 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 1011121314151617181920212223242526272829303132333435363738394041424344454647484950515253545556575859606162636465666768697071727374757677787980818283848586%

.. parsed-literal::

    Warning: --hwe observation counts vary by more than 10%, due to the X
    chromosome.  You may want to use a less stringent --hwe p-value threshold for X
    chromosome variants.
    

.. parsed-literal::

    8788899091929394959697989 done.
    Total genotyping rate is 0.999898.
    542508 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --het: 521405 variants scanned, report written to
    SampleData1/Fold_2/train_data.QC.het .
    [1] "SampleData1/Fold_2" "train_data"         "train_data.QC"     
    [4] "1"                 
    [1] "SampleData1" "Fold_2"     
    [1] "SampleData1"
    [1] 14
    [1] 14
    [1] 14
    Rscript Module1.R SampleData1/Fold_2  train_data train_data.QC 1
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_2/train_data.QC.log.
    Options in effect:
      --bfile SampleData1/Fold_2/train_data
      --extract SampleData1/Fold_2/train_data.QC.snplist
      --out SampleData1/Fold_2/train_data.QC
      --rel-cutoff 0.125
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (180 males, 200 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 542508 variants remaining.
    Using up to 8 threads (change this with --threads).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999898.
    542508 variants and 380 people pass filters and QC (before --rel-cutoff).
    Phenotype data is quantitative.
    Excluding 21103 variants on non-autosomes from relationship matrix calc.
    Relationship matrix calculation complete.
    0 people excluded by --rel-cutoff.
    Remaining sample IDs written to SampleData1/Fold_2/train_data.QC.rel.id .
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_2/train_data.QC.log.
    Options in effect:
      --a1-allele SampleData1/Fold_2/train_data.a1
      --bfile SampleData1/Fold_2/train_data
      --exclude SampleData1/Fold_2/train_data.mismatch
      --extract SampleData1/Fold_2/train_data.QC.snplist
      --keep SampleData1/Fold_2/train_data.QC.rel.id
      --make-bed
      --out SampleData1/Fold_2/train_data.QC
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (180 males, 200 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 542508 variants remaining.
    --exclude: 491965 variants remaining.
    --keep: 380 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999894.
    --a1-allele: 491965 assignments made.
    491965 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_2/train_data.QC.bed +
    SampleData1/Fold_2/train_data.QC.bim + SampleData1/Fold_2/train_data.QC.fam ...
    101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_3/train_data.log.
    Options in effect:
      --bfile SampleData1/SampleData1_QC
      --keep SampleData1/Fold_3/train_data.fam
      --make-bed
      --out SampleData1/Fold_3/train_data
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    475 people (227 males, 248 females) loaded from .fam.
    475 phenotype values loaded from .fam.
    --keep: 380 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate in remaining samples is 0.999871.
    551892 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_3/train_data.bed +
    SampleData1/Fold_3/train_data.bim + SampleData1/Fold_3/train_data.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_3/test_data.log.
    Options in effect:
      --bfile SampleData1/SampleData1_QC
      --keep SampleData1/Fold_3/test_data.fam
      --make-bed
      --out SampleData1/Fold_3/test_data
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    475 people (227 males, 248 females) loaded from .fam.
    475 phenotype values loaded from .fam.
    --keep: 95 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate in remaining samples is 0.999998.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_3/test_data.bed +
    SampleData1/Fold_3/test_data.bim + SampleData1/Fold_3/test_data.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_3/train_data.QC.log.
    Options in effect:
      --bfile SampleData1/Fold_3/train_data
      --geno 0.1
      --hwe 1e-6
      --maf 0.01
      --make-just-fam
      --mind 0.1
      --out SampleData1/Fold_3/train_data.QC
      --write-snplist
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (181 males, 199 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    0 people removed due to missing genotype data (--mind).
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999871.
    1 variant removed due to missing genotype data (--geno).
    --hwe: 841 variants removed due to Hardy-Weinberg exact test.
    8339 variants removed due to minor allele threshold(s)
    (--maf/--max-maf/--mac/--max-mac).
    542711 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    List of variant IDs written to SampleData1/Fold_3/train_data.QC.snplist .
    --make-just-fam to SampleData1/Fold_3/train_data.QC.fam ... done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_3/train_data.QC.log.
    Options in effect:
      --bfile SampleData1/Fold_3/train_data
      --extract SampleData1/Fold_3/train_data.QC.snplist
      --het
      --keep SampleData1/Fold_3/train_data.QC.fam
      --out SampleData1/Fold_3/train_data.QC
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (181 males, 199 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 542711 variants remaining.
    --keep: 380 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182%

.. parsed-literal::

    Warning: --hwe observation counts vary by more than 10%, due to the X
    chromosome.  You may want to use a less stringent --hwe p-value threshold for X
    chromosome variants.
    

.. parsed-literal::

    838485868788899091929394959697989 done.
    Total genotyping rate is 0.999872.
    542711 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --het: 521511 variants scanned, report written to
    SampleData1/Fold_3/train_data.QC.het .
    [1] "SampleData1/Fold_3" "train_data"         "train_data.QC"     
    [4] "1"                 
    [1] "SampleData1" "Fold_3"     
    [1] "SampleData1"
    [1] 14
    [1] 14
    [1] 14
    Rscript Module1.R SampleData1/Fold_3  train_data train_data.QC 1
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_3/train_data.QC.log.
    Options in effect:
      --bfile SampleData1/Fold_3/train_data
      --extract SampleData1/Fold_3/train_data.QC.snplist
      --out SampleData1/Fold_3/train_data.QC
      --rel-cutoff 0.125
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (181 males, 199 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 542711 variants remaining.
    Using up to 8 threads (change this with --threads).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999872.
    542711 variants and 380 people pass filters and QC (before --rel-cutoff).
    Phenotype data is quantitative.
    Excluding 21200 variants on non-autosomes from relationship matrix calc.
    Relationship matrix calculation complete.
    0 people excluded by --rel-cutoff.
    Remaining sample IDs written to SampleData1/Fold_3/train_data.QC.rel.id .
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_3/train_data.QC.log.
    Options in effect:
      --a1-allele SampleData1/Fold_3/train_data.a1
      --bfile SampleData1/Fold_3/train_data
      --exclude SampleData1/Fold_3/train_data.mismatch
      --extract SampleData1/Fold_3/train_data.QC.snplist
      --keep SampleData1/Fold_3/train_data.QC.rel.id
      --make-bed
      --out SampleData1/Fold_3/train_data.QC
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (181 males, 199 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 542711 variants remaining.
    --exclude: 492072 variants remaining.
    --keep: 380 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999868.
    --a1-allele: 492072 assignments made.
    492072 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_3/train_data.QC.bed +
    SampleData1/Fold_3/train_data.QC.bim + SampleData1/Fold_3/train_data.QC.fam ...
    101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_4/train_data.log.
    Options in effect:
      --bfile SampleData1/SampleData1_QC
      --keep SampleData1/Fold_4/train_data.fam
      --make-bed
      --out SampleData1/Fold_4/train_data
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    475 people (227 males, 248 females) loaded from .fam.
    475 phenotype values loaded from .fam.
    --keep: 380 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate in remaining samples is 0.999896.
    551892 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_4/train_data.bed +
    SampleData1/Fold_4/train_data.bim + SampleData1/Fold_4/train_data.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_4/test_data.log.
    Options in effect:
      --bfile SampleData1/SampleData1_QC
      --keep SampleData1/Fold_4/test_data.fam
      --make-bed
      --out SampleData1/Fold_4/test_data
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    475 people (227 males, 248 females) loaded from .fam.
    475 phenotype values loaded from .fam.
    --keep: 95 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate in remaining samples is 0.999896.
    551892 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_4/test_data.bed +
    SampleData1/Fold_4/test_data.bim + SampleData1/Fold_4/test_data.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_4/train_data.QC.log.
    Options in effect:
      --bfile SampleData1/Fold_4/train_data
      --geno 0.1
      --hwe 1e-6
      --maf 0.01
      --make-just-fam
      --mind 0.1
      --out SampleData1/Fold_4/train_data.QC
      --write-snplist
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (186 males, 194 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    0 people removed due to missing genotype data (--mind).
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    1 variant removed due to missing genotype data (--geno).
    --hwe: 828 variants removed due to Hardy-Weinberg exact test.
    8533 variants removed due to minor allele threshold(s)
    (--maf/--max-maf/--mac/--max-mac).
    542530 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    List of variant IDs written to SampleData1/Fold_4/train_data.QC.snplist .
    --make-just-fam to SampleData1/Fold_4/train_data.QC.fam ... done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_4/train_data.QC.log.
    Options in effect:
      --bfile SampleData1/Fold_4/train_data
      --extract SampleData1/Fold_4/train_data.QC.snplist
      --het
      --keep SampleData1/Fold_4/train_data.QC.fam
      --out SampleData1/Fold_4/train_data.QC
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (186 males, 194 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 542530 variants remaining.
    --keep: 380 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596%

.. parsed-literal::

    Warning: --hwe observation counts vary by more than 10%, due to the X
    chromosome.  You may want to use a less stringent --hwe p-value threshold for X
    chromosome variants.
    

.. parsed-literal::

    97989 done.
    Total genotyping rate is 0.999898.
    542530 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --het: 521557 variants scanned, report written to
    SampleData1/Fold_4/train_data.QC.het .
    [1] "SampleData1/Fold_4" "train_data"         "train_data.QC"     
    [4] "1"                 
    [1] "SampleData1" "Fold_4"     
    [1] "SampleData1"
    [1] 14
    [1] 14
    [1] 14
    Rscript Module1.R SampleData1/Fold_4  train_data train_data.QC 1
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_4/train_data.QC.log.
    Options in effect:
      --bfile SampleData1/Fold_4/train_data
      --extract SampleData1/Fold_4/train_data.QC.snplist
      --out SampleData1/Fold_4/train_data.QC
      --rel-cutoff 0.125
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (186 males, 194 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 542530 variants remaining.
    Using up to 8 threads (change this with --threads).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999898.
    542530 variants and 380 people pass filters and QC (before --rel-cutoff).
    Phenotype data is quantitative.
    Excluding 20973 variants on non-autosomes from relationship matrix calc.
    Relationship matrix calculation complete.
    0 people excluded by --rel-cutoff.
    Remaining sample IDs written to SampleData1/Fold_4/train_data.QC.rel.id .
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_4/train_data.QC.log.
    Options in effect:
      --a1-allele SampleData1/Fold_4/train_data.a1
      --bfile SampleData1/Fold_4/train_data
      --exclude SampleData1/Fold_4/train_data.mismatch
      --extract SampleData1/Fold_4/train_data.QC.snplist
      --keep SampleData1/Fold_4/train_data.QC.rel.id
      --make-bed
      --out SampleData1/Fold_4/train_data.QC
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    380 people (186 males, 194 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 542530 variants remaining.
    --exclude: 492135 variants remaining.
    --keep: 380 people remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999894.
    --a1-allele: 492135 assignments made.
    492135 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_4/train_data.QC.bed +
    SampleData1/Fold_4/train_data.QC.bim + SampleData1/Fold_4/train_data.QC.fam ...
    101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    

Check final files.
------------------

.. code:: ipython3

    import os
    
    # List of file names to check for existence
    files = [
        "train_data.QC.bed",
        "train_data.QC.bim",
        "train_data.QC.fam",
        "train_data.cov",
        "test_data.bed",
        "test_data.bim",
        "test_data.fam",
        "test_data.cov"
    ]
    
    # Directory where the files are expected to be found
    
    # Print the table header
    print("{:<20} {:<5} {:<5} {:<5} {:<5} {:<5}".format("File Name", "Fold 0", "Fold 1", "Fold 2", "Fold 3", "Fold 4"))
    
    # Loop through each file name in the list
    for file in files:
        # Create a list to store the existence status for each fold
        status = []
        # Check for each fold from 0 to 4
        for fold_number in range(5):
            # Check if the file exists in the specified directory for the given fold
            #print(os.path.join("./",filedirec, f"Fold_{fold_number}", file))
            if os.path.exists(filedirec + os.sep + "Fold_" + str(fold_number) + os.sep + file):
                status.append("yes")
            else:
                status.append("no")
        
        # Print the file name and its status for each fold
        print("{:<20} {:<5} {:<5} {:<5} {:<5} {:<5}".format(file, *status))
    


.. parsed-literal::

    File Name            Fold 0 Fold 1 Fold 2 Fold 3 Fold 4
    train_data.QC.bed    yes   yes   yes   yes   yes  
    train_data.QC.bim    yes   yes   yes   yes   yes  
    train_data.QC.fam    yes   yes   yes   yes   yes  
    train_data.cov       yes   yes   yes   yes   yes  
    test_data.bed        yes   yes   yes   yes   yes  
    test_data.bim        yes   yes   yes   yes   yes  
    test_data.fam        yes   yes   yes   yes   yes  
    test_data.cov        yes   yes   yes   yes   yes  
    

We will have the following directories if everything works fine.

::

   ├── Fold_0
   │   ├── test_data.bed
   │   ├── test_data.bim
   │   ├── test_data.cov
   │   ├── test_data.fam
   │   ├── test_data.log
   │   ├── train_data.a1
   │   ├── train_data.bed
   │   ├── train_data.bim
   │   ├── train_data.cov
   │   ├── train_data.fam
   │   ├── train_data.log
   │   ├── train_data.mismatch
   │   ├── train_data.QC.bed
   │   ├── train_data.QC.bim
   │   ├── train_data.QC.fam
   │   ├── train_data.QC.het
   │   ├── train_data.QC.log
   │   ├── train_data.QC.rel.id
   │   ├── train_data.QC.sexcheck
   │   ├── train_data.QC.snplist
   │   ├── train_data.QC.valid
   │   └── train_data.valid.sample
   ├── Fold_1
    
   ├── Fold_2
    
   ├── Fold_3
    
   ├── Fold_4

   ├── PeopleWithPhenotype.txt
   ├── SampleData1.bed
   ├── SampleData1.bim
   ├── SampleData1.cov
   ├── SampleData1.fam
   ├── SampleData1.gz
   ├── SampleData1.height
   ├── SampleData1_QC.bed
   ├── SampleData1_QC.bim
   ├── SampleData1_QC.cov
   ├── SampleData1_QC.fam
   └── SampleData1_QC.log

Important Note
--------------

1. Kindly ensure you have all the files required for the next step after
   the completion of this step.
2. It is better to pass the dataset on which quality controls have
   already been performed.
3. We considered genotype files for which chromosomes 1 to 22 are
   available, and sex information is present. If you want to use the
   function in the script to check the sex information, ensure the
   genotype file includes other chromosomes as well.
4. Go through the logs generated by the code if an error occurs. Even if
   no error occurs, it is always good to ensure that the log file does
   not produce any errors.

