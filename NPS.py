#!/usr/bin/env python
# coding: utf-8

#  
# # NPS
# 
# NPS implements a non-parametric polygenic risk prediction algorithm. NPS transforms genetic data into an orthogonal domain called "eigenlocus space". Then, it re-weights GWAS effect sizes by partitioning genetic variations into trenches and measuring the predictive power of each trench in an independent training cohort.
# 
# ## Installation
# 
# ```bash
# https://github.com/sgchun/nps
# wget https://github.com/sgchun/nps/archive/1.1.1.tar.gz
# tar -zxvf  1.1.1.tar.gz
# cd nps-1.1.1/
# make .
# ```
# 
# Copy the files in `nps-1.1.1` into the current working directory. It ensures smooth execution, as some code files rely on bash files from NPS, and placing the files in the current directory helps with execution.
# 
# ## Genotype and GWAS Data Processing
# 
# NPS requires both GWAS and genotype data in a specific format. Genotype data should be in a `dosage.gz` file, and you may need `qctool` ([documentation](https://www.chg.ox.ac.uk/~gav/qctool_v2/documentation/download.html)) to process the information.
# 
# 1. Genotype data was processed using the following quality controls:
#    ```python
#    command = [
#        "./plink",
#        "--maf", str(0.2),
#        "--geno", str(0.001),
#        "--hwe", str(0.00001),
#        "--bfile", traindirec + os.sep + newtrainfilename,
#        "--indep-pairwise", p1_val, p2_val, p3_val,
#        "--out", traindirec + os.sep + trainfilename
#    ]
#    subprocess.run(command)
#    ```
#    
# 2. The number of SNPs in the GWAS and genotype data should be the same, with no missing values or `NA` in the genotype data.
# 
# 3. The name of the genotype data (train/test) should follow specific prefixes; otherwise, it may not work as it is a bash-based tool and can be difficult to troubleshoot.
# 
# 4. Genotype data for all chromosomes should be present, or the GWAS tool may crash as it runs for all chromosomes.
# 
# ## Hyperparameter  
# 
# NPS considers multiple hyperparameters. One of the hyperparameters we considered is the window size, and new betas are generated for each window size.
# 
# 
# ## Train/Testing sets
# 
# NPS generates multiple sets of betas. We used these betas and calculated the PRS using Plink, but one can also use the default method provided by NPS to calculate the performance for the test set.
# 
# 
# 
# ## Input files for NPS
# 
# **Taken from NPS documentation.**
# 
# 
# 
# To run NPS, you need the following set of files: 
# 
# 1. **GWAS summary statistics.** This is a *tab-delimited* text file format with the following seven required columns: 
#      - **chr**: chromosome name starting with "chr." Currently, NPS expects only chromosomes 1-22. Chromosome names should be designated by "chr1", ..., "chr22".
#      - **pos**: base position of SNP.
#      - **ref** and **alt**: reference and alternative alleles of SNP, respectively.
#      - **reffreq**: allele frequency of reference allele in the discovery GWAS cohort. 
#      - **pval**: p-value of association. 
#      - **effalt**: estimated *per-allele* effect size of *the alternative allele*. For case/control GWAS, log(Odds Ratio) should be used. NPS will convert **effalt** to effect sizes relative to *the standardized genotype* internally using **reffreq**.  
#      ```
#      chr        pos     ref     alt     reffreq pval    effalt
#      chr1       11008   C       G       0.9041  0.1126  -0.0251
#      chr1       11012   C       G       0.9041  0.1126  -0.0251
#      chr1       13116   T       G       0.8307  0.615   0.0071
#      chr1       13118   A       G       0.8307  0.615   0.0071
#      chr1       14464   A       T       0.8386  0.476   -0.0105
#      ...
#      ```
# 
# 2. **Training genotype files.** Training genotype files should be in the qctool dosage format and named as "chrom*N*.*DatasetID*.dosage.gz" for each chromosome. Genotype files in bgen format can be converted to the dosage files by running [qctool](https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/examples/converting.html) with the `-ofiletype dosage` option. NPS allows only *biallelic* variants. 
#    
# 3. **Training sample file (.fam).** Sample information of training cohort should be provided in [PLINK FAM format](https://www.cog-genomics.org/plink2/formats#fam). The samples in the .fam file should appear in the exactly same order as in the genotype dosage files. The sex of sample (5-th column) is optional ("0" or "-9" for missing; "1" for male; "2" for female). If the sex is provided, NPS will incorporate the sex covariate in the PRS model. The 6-th column is for phenotype data and can be specified here or in a separeate phenotype file. 
#    ```
#    trainF2  trainI2  0  0  1 -9
#    trainF3  trainI3  0  0  2 -9
#    trainF39 trainI39 0  0  1 -9
#    trainF41 trainI41 0  0  2 -9
#    trainF58 trainI58 0  0  1 -9
#    ```
# 4. **Training phenotype file (.phen).** Phenotypes of the .fam file can be overridden by a .phen file (use `nps_init.R --train-phen` option). This is a tab-delimited file with three columns: "FID", "IID", and "Outcome". FID and IID correspond to the family and individual IDs in the .fam file. The name of phenotype should be "Outcome". Binary phenotypes (case/control) are specified by "1" and "2", respectively; "0" and "-9" denote missing phenotype. For quantitative phenotypes, "-9" represents a missing phenotype value. 
#    ```
#    FID   IID    Outcome
#    trainF2  trainI2  1
#    trainF39 trainI39 1
#    trainF3  trainI3  2
#    trainF41 trainI41 2
#    trainF58 trainI58 1
#    ```
# 5. **Validation genotype files.** Validation genotypes can be in the dosage or .bgen format. If they are in .bgen format, the files should be named as "chrom*N*.*DatasetID*.bgen". 
# 
# 6. **Validation sample file (.fam).** Similar to the training .fam file. 
# 
# 7. **Training phenotype file (.phen).** Similar to training .phen file. 

# ## GWAS File Processing for NPS
# 
# When the effect size relates to disease risk and is given as an odds ratio (OR) rather than BETA (for continuous traits), the polygenic risk score (PRS) is computed as a product of ORs. To simplify this calculation, take the natural logarithm of the OR so that the PRS can be computed using summation instead.
# 
# We will modify the GWAS accordingly. NPS requires modifying the GWAS in a particular format, and we modified it when calculating PRS using NPS.
# 
# 
# 

# In[9]:


import os
import pandas as pd
import numpy as np
import sys

filedirec = sys.argv[1]

#filedirec = "SampleData1"
 
#filedirec = "migraine_0"

def check_phenotype_is_binary_or_continous(filedirec):
    # Read the processed quality controlled file for a phenotype
    df = pd.read_csv(filedirec+os.sep+filedirec+'_QC.fam',sep="\s+",header=None)
    column_values = df[5].unique()
 
    if len(set(column_values)) == 2:
        return "Binary"
    else:
        return "Continous"

# Read the GWAS file.
GWAS = "./"+filedirec + os.sep + filedirec+".gz"
df = pd.read_csv(GWAS,compression= "gzip",sep="\s+")

if check_phenotype_is_binary_or_continous(filedirec)=="Binary":
    if "BETA" in df.columns.to_list():
        # For Binary Phenotypes.
        df["OR"] = np.exp(df["BETA"])
        df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'OR', 'INFO', 'MAF']]
 
    else:
        # For Binary Phenotype.
        df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'OR', 'INFO', 'MAF']]
 
elif check_phenotype_is_binary_or_continous(filedirec)=="Continous":
    if "BETA" in df.columns.to_list():
        # For Continous Phenotype.
        df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]

    else:
        df["BETA"] = np.log(df["OR"])
        df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]


df.to_csv(filedirec + os.sep +filedirec+"_NPS.txt",sep="\t",index=False)
print(df.head().to_markdown())
print("Length of DataFrame!",len(df))

df.to_csv(filedirec + os.sep +filedirec+"_NPS.gz", index=False,sep="\t",compression='gzip')


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

# In[10]:


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
prs_result = pd.DataFrame(columns=["clump_p1", "clump_r2", "clump_kb", "p_window_size", "p_slide_size", "p_LD_threshold",
                                   "pvalue", "numberofpca","numberofvariants","Train_pure_prs", "Train_null_model", "Train_best_model",
                                   "Test_pure_prs", "Test_null_model", "Test_best_model"])


# ### Define Helper Functions
# 
# 1. **Perform Clumping and Pruning**
# 2. **Calculate PCA Using Plink**
# 3. **Fit Binary Phenotype and Save Results**
# 4. **Fit Continuous Phenotype and Save Results**
# 

# In[11]:


import os
import subprocess
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import explained_variance_score


def perform_clumping_and_pruning_on_individual_data(traindirec, newtrainfilename,numberofpca, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
    
    # Perform additional quality controls for NPS.
    #    "--maf",str(0.2),
    #    "--geno",str(0.001),
    #    "--hwe",str(0.00001),
    
    command = [
    "./plink",
    "--maf",str(0.2),
    "--geno",str(0.001),
    "--hwe",str(0.00001),
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
def fit_binary_phenotype_on_PRS(traindirec, newtrainfilename,p,window_shift, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
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
          
                    "window_shift":window_shift,
                    
                     

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
def fit_continous_phenotype_on_PRS(traindirec, newtrainfilename,p,window_shift, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
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
                 
                    "window_shift":window_shift,                     

                    "Train_pure_prs":explained_variance_score(phenotype_train["Phenotype"],prs_train['SCORE'].values),
                    "Train_null_model":explained_variance_score(phenotype_train["Phenotype"],train_null_predicted),
                    "Train_best_model":explained_variance_score(phenotype_train["Phenotype"],train_best_predicted),
                    
                    "Test_pure_prs":explained_variance_score(phenotype_test["Phenotype"],prs_test['SCORE'].values),
                    "Test_null_model":explained_variance_score(phenotype_test["Phenotype"],test_null_predicted),
                    "Test_best_model":explained_variance_score(phenotype_test["Phenotype"],test_best_predicted),
                    
                }, ignore_index=True)

          
                prs_result.to_csv(traindirec+os.sep+Name+os.sep+"Results.csv",index=False)
     
    return


# ## Execute NPS

# In[ ]:


import shutil
# Define a global variable to store results
prs_result = pd.DataFrame()

          
def transform_nps_data(traindirec, newtrainfilename,p,windowsize, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
 
    ### First perform clumping on the file and save the clumpled file.
    perform_clumping_and_pruning_on_individual_data(traindirec, newtrainfilename,p, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile)
    gwas_file = filedirec + os.sep + filedirec + '_NPS.gz'
    bim_file = traindirec + os.sep + newtrainfilename+".clumped.pruned.bim"
 
    # Read the files
    df = pd.read_csv(gwas_file, sep="\s+",compression='gzip')
    bim = pd.read_csv(bim_file, delim_whitespace=True, header=None)

    print(len(df))
    print(len(bim))

    # Create a 'match' column to find common SNPs
    bim['match'] = bim[0].astype(str) + "_" + bim[3].astype(str) + "_" + bim[4].astype(str) + "_" + bim[5].astype(str)
    df['match'] = df['CHR'].astype(str) + "_" + df['BP'].astype(str) + "_" + df['A1'].astype(str) + "_" + df['A2'].astype(str)

    # Drop duplicates based on the 'match' column
    df.drop_duplicates(subset='match', inplace=True)
    bim.drop_duplicates(subset='match', inplace=True)

    # Filter dataframes to keep only matching SNPs
    df = df[df['match'].isin(bim['match'].values)]
    bim = bim[bim['match'].isin(df['match'].values)]

    print(len(bim))
    print(len(df))
    
    print(bim.head())
    print(df.head())
    
    del df["match"]
    del bim["match"]
    df.to_csv(traindirec+os.sep+filedirec+"_NPS.gz",compression="gzip",sep="\t",index=None)   
    bim.to_csv(traindirec + os.sep +  "commonsnps.txt",sep="\t",index=None)
    
    
    command = [
    './plink', 
    '--bfile', traindirec+os.sep+newtrainfilename,
    '--extract', traindirec + os.sep +  "commonsnps.txt", 
    '--make-bed', 
    '--out', traindirec+os.sep+newtrainfilename+".clumped.pruned.NPS"
    ]
    subprocess.run(command)
 

        
    # Also extract the PCA at this point for both test and training data.
    calculate_pca_for_traindata_testdata_for_clumped_pruned_snps(traindirec, newtrainfilename,p)

     
    # Extract the exact list of the SNPs common between the GWAS and the cleaned genotype Data.
    df1  = pd.read_csv(traindirec+os.sep+newtrainfilename+".clumped.pruned.NPS.bim",sep="\t",names=["CHR","SNP","X","POS","A1","A2"])
    # Read the original GWAS 
    df2 = pd.read_csv(traindirec+os.sep+filedirec+"_NPS.gz",compression= "gzip",sep="\s+")
    columns_to_compare = ['SNP', 'A1', 'A2']
    merged_df = pd.merge(df1, df2, on=columns_to_compare)
    merged_df["SNP"].to_csv(traindirec+os.sep+"exactmatch.txt",header=False,index=False)
     
    
    # Delete files generated in the previous iteration.
    files_to_remove = [
        traindirec + os.sep + filedirec + ".NPS",
        traindirec + os.sep + filedirec + ".NPSfinal",
        traindirec + os.sep + "npsdat",
    ]

    # Loop through the files and directories and remove them if they exist
    for file_path in files_to_remove:
        if os.path.exists(file_path):
            if os.path.isfile(file_path):
                os.remove(file_path)
                print(f"Removed file: {file_path}")
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
                print(f"Removed directory: {file_path}")
        else:
            print(f"File or directory does not exist: {file_path}")   
    
    
    for chromosome in range(1, 23):
         
        # Recode training data.
        plink_command = [
            "./plink",
            "--bfile", traindirec+os.sep+newtrainfilename+".clumped.pruned.NPS",
            "--chr", str(chromosome),
            '--extract', traindirec + os.sep +  "commonsnps.txt", 
            "--recode", "oxford",
        
            "--out", traindirec+os.sep+"chrom"+str(chromosome)+".train",
        ]
        try:
            subprocess.run(plink_command)
        except:
            pass

        # Extract the Same SNPs from the test Data.
        plink_command = [
            "./plink",
            "--bfile", folddirec+os.sep+testfilename,
            "--chr", str(chromosome),
     
            '--extract', traindirec + os.sep +  "commonsnps.txt", 
            "--recode", "oxford",
 
            "--out", traindirec+os.sep+"chrom"+str(chromosome)+".test",
        ]
        try:
            subprocess.run(plink_command)
        except:
            pass
        
        # Convert data dosage information. Ensure the prefix and the sufix are the same.
        
        plink_command = [
            "./qctools",
            "-g", traindirec+os.sep+"chrom"+str(chromosome)+".train.gen",
            "-og", traindirec+os.sep+"chrom"+str(chromosome)+".train.dosage.gz"
        ]
        try:
            subprocess.run(plink_command)
        except:
            pass
        
        plink_command = [
            "./qctools",
            "-g", traindirec+os.sep+"chrom"+str(chromosome)+".test.gen",
            "-og", traindirec+os.sep+"chrom"+str(chromosome)+".test.dosage.gz"
        ]
        try:
            subprocess.run(plink_command)
        except:
            pass
        
        # Once we have the data in dosage format. 
        # We read the file for chromosomes.
        # Convert the chromosome to 1 to 01. or 2 to 02
        # Converted position to string.
        # Converted RSID  to chromosome_BP_A1_A2
        # Rename the files.
        def change_the_dosage_files():
            def format_two_digits(number):
                return f"{number:02d}"
            df = pd.read_csv(traindirec+os.sep+"chrom"+str(chromosome)+".train.dosage.gz", compression='gzip',sep="\s+")
            df['chromosome'] = str(format_two_digits(chromosome))
            #print(df.iloc[0:6].head())
            df['position'] = df['position'].astype(str)
            df['SNPID'] = str(chromosome)+"_"+df['position']+":"+df['position']+":"+df['alleleA']+":"+df['alleleB']
            #print(df.iloc[0:6].head())
            #exit(0)
            df['rsid'] = str(chromosome)+"_"+df['position']+":"+df['position']+":"+df['alleleA']+":"+df['alleleB']
            df['position'] = df['position'].astype(str)
            df.iloc[:, 6:] = df.iloc[:, 6:].apply(pd.to_numeric, errors='coerce')
            df = df.astype(str)
            df.to_csv(traindirec+os.sep+"chrom"+str(chromosome)+".change.train.dosage.gz", compression='gzip', index=False,sep=" ")

            df = pd.read_csv(traindirec+os.sep+"chrom"+str(chromosome)+".test.dosage.gz", compression='gzip',sep="\s+")
            df['chromosome'] = str(format_two_digits(chromosome))
            #print(df.iloc[0:6].head())
            df['position'] = df['position'].astype(str)
            df['SNPID'] = str(chromosome)+"_"+df['position']+":"+df['position']+":"+df['alleleA']+":"+df['alleleB']
            #print(df.iloc[0:6].head())
            #exit(0)
            df['rsid'] = str(chromosome)+"_"+df['position']+":"+df['position']+":"+df['alleleA']+":"+df['alleleB']
            df['position'] = df['position'].astype(str)
            df.iloc[:, 6:] = df.iloc[:, 6:].apply(pd.to_numeric, errors='coerce')
            df = df.astype(str)
            df.to_csv(traindirec+os.sep+"chrom"+str(chromosome)+".change.test.dosage.gz", compression='gzip', index=False,sep=" ")
            #raise
        
        try:
            change_the_dosage_files()
            
        except:
            pass     
        
 
    
    
    # Execute NPS. Standardize genotypes 
    command = "./nps-1.1.1/run_all_chroms.sh ./nps-1.1.1/sge/nps_stdgt.job "+traindirec+" "+"change.train"
    #try:
    subprocess.run(command, shell=True)
    #except:
    #    print("First step did not work!")
    #    return
    
    
    # Modify the fam file for train. As required by the NPS. The name should be the same as the dosage file.
    tempphenotype_train = pd.read_table(traindirec+os.sep+newtrainfilename+".fam", sep="\s+",header=None)
    x = tempphenotype_train[[0, 1, 5]] 
    x.columns = ['FID', 'IID', 'Outcome']
    file_path = os.path.join(traindirec, "change.train.pheno")
    x.to_csv(file_path, index=False, sep='\t')
    tempphenotype_train.to_csv(traindirec+os.sep+"change.train.fam", sep="\t",index=False,header=None)
    
    
    tempphenotype_train = pd.read_table(traindirec+os.sep+testfilename+".fam", sep="\s+",header=None)
    x = tempphenotype_train[[0, 1, 5]] 
    x.columns = ['FID', 'IID', 'Outcome']
    file_path = os.path.join(traindirec, "change.test.pheno")
    x.to_csv(file_path, index=False, sep='\t')
    tempphenotype_train.to_csv(traindirec+os.sep+"change.test.fam", sep="\t",index=False,header=None)
    
    
    
    # Read the original GWAS and find the common SNPs between the GWAS and the Genotype data.
    
    GWAS = traindirec+os.sep+filedirec+"_NPS.gz"
    original_df = pd.read_csv(GWAS,compression= "gzip",sep="\s+")
    bim = pd.read_csv(traindirec+os.sep+newtrainfilename+".clumped.pruned.NPS"+".bim",header=None,sep="\s+")[1].values
    original_df = original_df[original_df["SNP"].isin(bim)]
    
    
       
    if "BETA" in original_df.columns.to_list():
        # For Continous Phenotype.
        #original_df = original_df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]
        transformed_df = pd.DataFrame({
            'chr': 'chr' + original_df['CHR'].astype(str),
            'pos': original_df['BP'],
            'ref': original_df['A1'],
            'alt': original_df['A2'],
            'reffreq': original_df['MAF'],
            #'reffreq': 1 - original_df['MAF'],    
            'SNP': original_df['SNP'],
            'info': original_df['INFO'],
            'rs': original_df['CHR'].astype(str)+'_' + original_df['BP'].astype(str) + ':' + original_df['BP'].astype(str) + ':' + original_df['A1'] + ':' + original_df['A2'],
            'pval': original_df['P'],
            'effalt': original_df['BETA']
        })
        
    else:
        transformed_df = pd.DataFrame({
            'chr': 'chr' + original_df['CHR'].astype(str),
            'pos': original_df['BP'],
            'ref': original_df['A1'],
            'alt': original_df['A2'],
            'reffreq': original_df['MAF'],
            #'reffreq':1 - original_df['MAF'],    
            'SNP': original_df['SNP'],
            'info': original_df['INFO'],
            'rs': original_df['CHR'].astype(str)+'_' + original_df['BP'].astype(str) + ':' + original_df['BP'].astype(str) + ':' + original_df['A1'] + ':' + original_df['A2'],
            'pval': original_df['P'],
            'effalt': original_df['OR']
        })

    transformed_df.to_csv(traindirec + os.sep +filedirec+".NPS",sep="\t",index=False)
 
    hell = 0
    for chromosome in range(1, 23):
        try:
            hell = hell+len(pd.read_csv(traindirec+os.sep+"chrom"+str(chromosome)+".change.train.snpinfo"))
        except:
            pass
    print("ss",hell)
    #Ensure the number of entries in the chrom.ALL.change.train.snpinfo is the same as the number entries in GWAS file
        
    #raise
    
    windowsize = windowsize
    # 2 Configure an NPS run
    command = [
        'Rscript',
        'npsR/nps_init.R',
        '--gwas', traindirec + os.sep +filedirec+".NPS",
        '--train-dir', traindirec,
        '--train-fam',traindirec+os.sep+"change.train.fam",
        '--train-phen',traindirec+os.sep+"change.train.pheno",
        '--train-dataset', "change.train",
        '--window-size', str(windowsize),
        '--out', traindirec+os.sep+'npsdat'
    ]
   
    subprocess.run(command)
    print(" ".join(command))
    #raise
    
    #3. Set up a special partition for GWAS-significant SNPs
    
    command = [
        './run_all_chroms.sh',
        'sge/nps_gwassig.job',
        traindirec+os.sep+'npsdat'
    ]
    subprocess.run(command)
    print(" ".join(command)) 
    window_shifts = [0,windowsize*0.25,windowsize*0.5,windowsize*0.75]
    
    # 4 Set up the decorrelated "eigenlocus" space.
    for window_shift in window_shifts:
        command = [
            './run_all_chroms.sh',
            'sge/nps_decor_prune.job',
            traindirec+os.sep+'npsdat',
            str(window_shift)
        ]
        print(" ".join(command))
        subprocess.run(command)
        
    #5. Partition the rest of genome.
    command = [
        'Rscript',
        'npsR/nps_prep_part.R',
        traindirec+os.sep+'npsdat',
        '10',  # Number of partitions on eigenvalues
        '10'  # Number of partitions on estimated effects
    ]
    print(" ".join(command))
    subprocess.run(command)
    # Calculate partitioned polygenic risk scores in the training cohort
 
   
    for window_shift in window_shifts:
        command = [
            './run_all_chroms.sh',
            'sge/nps_part.job',
            traindirec+os.sep+'npsdat',
            str(window_shift)
        ]
        print(" ".join(command))
        subprocess.run(command)
    
    
    
    
    #6. Estimate shrinkage weights for each partition
    command = [
        'Rscript',
        'npsR/nps_reweight.R',
        traindirec+os.sep+'npsdat'
    ]
    print(" ".join(command))
   
    subprocess.run(command)
    
    #raise
    # 7. Evaluate the accuracy of trained prediction model in a validation cohort. 
    for window_shift in window_shifts:
        command = [
            './run_all_chroms.sh',
            'sge/nps_score.dosage.job',
              traindirec+os.sep+'npsdat'+os.sep,
            traindirec+os.sep,
               "change.test", 
            
            str(int(window_shift))
        ]
        print(" ".join(command))
        subprocess.run(command)
    
    command = [
        'Rscript',
        'npsR/nps_val.R',
        '--out',traindirec+os.sep+'npsdat',
        "--val-dataset", "change.test"
        
    ]
    
    subprocess.run(command)
    
    
    
    #Rscript npsR/nps_val.R --out testdata/Test1/npsdat --val-dataset Test1.val 
    import glob
    for window_shift in [0,windowsize*0.25,windowsize*0.5,windowsize*0.75]:
        file_pattern = 'change.train.win_'+str(int(window_shift))+'.adjbetahat_pg.chrom*.txt'
        print(file_pattern)
        file_paths = glob.glob(os.path.join(traindirec+os.sep+"npsdat", file_pattern))
        print(file_paths)
        for x in file_paths:
            print(x)
        if len(file_paths)==0:
            print("files for specific window :",window_shift," does not exist")
            continue
        dataframes_list = []
        
        for file_path in file_paths:
            try:
                df = pd.read_csv(file_path, delimiter='\t', header=None)  # Adjust delimiter if needed
                dataframes_list.append(df)
            except:
                pass
            
        # Concatenate all DataFrames vertically
        result_df = pd.concat(dataframes_list, ignore_index=True)
        print(result_df)
        tempgwas = pd.read_csv(traindirec+os.sep+filedirec+".NPS",sep="\t")
        tempgwas["NEWBETAS"] = result_df[0].values 
        
        if check_phenotype_is_binary_or_continous(filedirec)=="Binary":
            tempgwas["NEWBETAS"] = np.exp(tempgwas["NEWBETAS"])
        else:
            pass
        

        tempgwas.iloc[:,[5,2,10]].to_csv(traindirec+os.sep+filedirec+".NPSfinal",index=False,sep="\t")
        print(tempgwas.iloc[:,[5,2,10]].head())
         
    


        # Caluclate Plink Score.
        command = [
            "./plink",
             "--bfile", traindirec+os.sep+newtrainfilename+".clumped.pruned",
            ### SNP column = 3, Effect allele column 1 = 4, OR column=9
            "--score", traindirec+os.sep+filedirec+".NPSfinal", "1", "2", "3", "header",
            "--q-score-range", traindirec+os.sep+"range_list",traindirec+os.sep+"SNP.pvalue",
            "--extract", traindirec+os.sep+trainfilename+".valid.snp",
            "--out", traindirec+os.sep+Name+os.sep+trainfilename
        ]

        subprocess.run(command)



        command = [
            "./plink",
            "--bfile", folddirec+os.sep+testfilename+".clumped.pruned",
            ### SNP column = 3, Effect allele column 1 = 4, Beta column=12
            "--score", traindirec+os.sep+filedirec+".NPSfinal", "1", "2", "3", "header",
            "--q-score-range", traindirec+os.sep+"range_list",traindirec+os.sep+"SNP.pvalue",
            "--extract", traindirec+os.sep+trainfilename+".valid.snp",
            "--out", folddirec+os.sep+Name+os.sep+testfilename
        ]
        subprocess.run(command)

 

        if check_phenotype_is_binary_or_continous(filedirec)=="Binary":
            print("Binary Phenotype!")
            fit_binary_phenotype_on_PRS(traindirec, newtrainfilename,p,window_shift, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile)
        else:
            print("Continous Phenotype!")
            fit_continous_phenotype_on_PRS(traindirec, newtrainfilename,p,window_shift, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile)



result_directory = "NPS"
# Nested loops to iterate over different parameter values
windowsizes = [80]
create_directory(folddirec+os.sep+result_directory)
for p1_val in p_window_size:
 for p2_val in p_slide_size: 
  for p3_val in p_LD_threshold:
   for c1_val in clump_p1:
    for c2_val in clump_r2:
     for c3_val in clump_kb:
      for p in numberofpca:
       for windowsize in windowsizes: 
        transform_nps_data(folddirec, newtrainfilename, p,windowsize, str(p1_val), str(p2_val), str(p3_val), str(c1_val), str(c2_val), str(c3_val), result_directory, pvaluefile)


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
# python NPS.py 0
# python NPS.py 1
# python NPS.py 2
# python NPS.py 3
# python NPS.py 4
# ```
# 
# The following files should exist after the execution:
# 
# 1. `SampleData1/Fold_0/NPS/Results.csv`
# 2. `SampleData1/Fold_1/NPS/Results.csv`
# 3. `SampleData1/Fold_2/NPS/Results.csv`
# 4. `SampleData1/Fold_3/NPS/Results.csv`
# 5. `SampleData1/Fold_4/NPS/Results.csv`
# 

# ### Check the results file for each fold.

# In[ ]:


import os
 
result_directory = "NPS"
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

# In[ ]:


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
    important_columns = [
        'clump_p1',
        'clump_r2',
        'clump_kb',
        'p_window_size',
        'p_slide_size',
        'p_LD_threshold',
        'pvalue',
        'referencepanel',
   
        
        'window_shift',
        'numberofpca',
        'tempalpha',
        'l1weight',
         
       
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

# In[ ]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'notebook')

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




# In[ ]:




