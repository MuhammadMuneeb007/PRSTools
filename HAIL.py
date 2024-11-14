#!/usr/bin/env python
# coding: utf-8

#  
# # HAIL
# 
# In this notebook, we used [Hail](https://github.com/hail-is/hail) to calculate Polygenic Risk Scores (PRS). Hail does not calculate new betas but instead uses the existing weights from the GWAS and applies a custom formula for the calculation.
# 
# We followed this tutorial to calculate the PRS:  
# [https://nbviewer.org/github/ddbj/imputation-server-wf/blob/main/Notebooks/hail-prs-tutorial.ipynb](https://nbviewer.org/github/ddbj/imputation-server-wf/blob/main/Notebooks/hail-prs-tutorial.ipynb)
# 
# ## Basic Process
# 
# 1. Read the genotype data.
# 2. Convert it to VCF.
# 3. Use Beagle to convert the data to Beagle format.
# 4. Convert the Beagle format to Hail format.
# 5. Pass the data to Hail, GWAS, and genotype.
# 6. Calculate PRS using Hail.
#      
# ### Genotype Data Processing
# 
# 1. **Convert genotype data to VCF format.**  
#    Hail requires data in `beagle.vcf.gz` format. The first step is to convert the `bed`, `bim`, and `fam` files to VCF format. A simple approach is to extract genotype data for each chromosome and then convert the data to VCF.
# 
#     ```bash
#     plink --bfile traindirec/newtrainfilename.clumped.pruned --chr 22 --make-bed --out traindirec/newtrainfilename.clumped.pruned.22
#     plink --bfile traindirec/newtrainfilename.clumped.pruned.22 --chr 22 --recode vcf --out traindirec/hail.train.22
#     ```
# 
# 2. **Download the phased reference panel.**  
#    To run the following commands, download the phased reference panel from 1000 Genomes and place it in the current working directory:  
#    [1000 Genomes Reference Panel](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)
# 
#     ```bash
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#     ```
# 
# 3. **Alternatively, download the reference panel from Hail.**  
#    - Reference panels:  
#      [Hail Reference Panel](https://sc.ddbj.nig.ac.jp/en/advanced_guides/imputation_server_tutorial)
# 
#    - Files (2.6 GB):  
#      - [test-data.GRCh37.vcf.gz](https://zenodo.org/records/6650681/files/test-data.GRCh37.vcf.gz?download=1) (1.3 GB, `md5:aff8bca4689cc70f6dbc1c3296590458`)  
#      - [test-data.GRCh38.vcf.gz](https://zenodo.org/records/6650681/files/test-data.GRCh38.vcf.gz?download=1) (1.3 GB, `md5:d28a741e820444ca926f7b0d5ac2e196`)
# 
# 4. **Download genetic distances from the Beagle website.**  
#    [Beagle Genetic Maps](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)
# 
# 5. **Download Beagle.**  
#    - Beagle documentation:  
#      [Beagle Documentation](https://faculty.washington.edu/browning/beagle/beagle_5.4_18Mar22.pdf)
#    - Beagle download page:  
#      [Download Beagle](https://faculty.washington.edu/browning/beagle/beagle.html)
# 
# 6. **Run the Beagle command.**  
#    Run the following command to perform phasing and imputation:
# 
#     ```bash
#     java -Xmx50g -jar beagle gt=traindirec/hail.train.22.vcf ref=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr22.GRCh37.map out=traindirec/beagle.hail.train.22
#     ```
# 
# 7. **Follow the remaining process.**  
#    The rest of the process is straightforward and can be followed from this GitHub repository:  
#    [PRS on Hail GitHub Repository](https://github.com/hacchy1983/prs-on-hail-public)
#  
#  

# ## GWAS file processing for HAIL for Binary Phenotypes.
# When the effect size relates to disease risk and is thus given as an odds ratio (OR) rather than BETA (for continuous traits), the PRS is computed as a product of ORs. To simplify this calculation, take the natural logarithm of the OR so that the PRS can be computed using summation instead.

# In[24]:


import os
import pandas as pd
import numpy as np
import sys

filedirec = sys.argv[1]

#filedirec = "SampleData1"
#filedirec = "asthma_19"
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
GWAS = filedirec + os.sep + filedirec+".gz"
df = pd.read_csv(GWAS,compression= "gzip",sep="\s+")

 
if "BETA" in df.columns.to_list():
    # For Continous Phenotype.
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]

else:
    df["BETA"] = np.log(df["OR"])
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]


    
df = df.rename(columns={
    'CHR': 'chr_name',
    'BP': 'chr_position',
    'A1': 'effect_allele',
    'A2': 'other_allele',
    'BETA': 'effect_weight'
})

# Selecting the relevant columns
df = df[['chr_name', 'chr_position', 'effect_allele', 'other_allele', 'effect_weight']]
# Remove duplicates based on 'chr_name' and 'chr_position'
print("Length of DataFrame!",len(df))
df = df.drop_duplicates(subset=['chr_name', 'chr_position', 'effect_allele', 'other_allele'])
print("Length of DataFrame!",len(df))

df.to_csv(filedirec + os.sep +"Hail.txt",sep="\t",index=False)
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

# In[20]:


import os
import subprocess
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import explained_variance_score


def perform_clumping_and_pruning_on_individual_data(traindirec, newtrainfilename,numberofpca, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
    
    command = [
    "./plink",
    "--bfile", traindirec+os.sep+newtrainfilename,
    "--maf",str(0.2),
    "--geno",str(0.001),
    "--hwe",str(0.00001),
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
def fit_binary_phenotype_on_PRS(traindirec, newtrainfilename,p, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
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

            prs_train = pd.read_table(traindirec+os.sep+Name+os.sep+"train_hail_prs.csv",sep=",",index_col=0)
            print(prs_train)
            prs_train[['FID', 'IID']] = prs_train['subjectID'].str.split('_', expand=True)

            # Drop the original 'Combined' column if it's no longer needed
            prs_train = prs_train.drop(columns=['subjectID'])
            prs_train.rename(columns={'PGS002724': 'SCORE'}, inplace=True)
   
            
            

            prs_train['FID'] = prs_train['FID'].astype(str)
            prs_train['IID'] = prs_train['IID'].astype(str)
         
             
            prs_test = pd.read_table(traindirec+os.sep+Name+os.sep+"test_hail_prs.csv",sep=",",index_col=0)
            print(prs_test)
            prs_test[['FID', 'IID']] = prs_test['subjectID'].str.split('_', expand=True)

            # Drop the original 'Combined' column if it's no longer needed
            prs_test = prs_test.drop(columns=['subjectID'])
            prs_test.rename(columns={'PGS002724': 'SCORE'}, inplace=True)
   
        
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
                #"pvalue": i,
                "numberofpca":p, 

                "tempalpha":str(tempalpha),
                "l1weight":str(l1weight),
              


                "Train_pure_prs":roc_auc_score(phenotype_train["Phenotype"].values,prs_train['SCORE'].values),
                "Train_null_model":roc_auc_score(phenotype_train["Phenotype"].values,train_null_predicted.values),
                "Train_best_model":roc_auc_score(phenotype_train["Phenotype"].values,train_best_predicted.values),

                "Test_pure_prs":roc_auc_score(phenotype_test["Phenotype"].values,prs_test['SCORE'].values),
                "Test_null_model":roc_auc_score(phenotype_test["Phenotype"].values,test_null_predicted.values),
                "Test_best_model":roc_auc_score(phenotype_test["Phenotype"].values,test_best_predicted.values),

            }, ignore_index=True)


            prs_result.to_csv(traindirec+os.sep+Name+os.sep+"Results.csv",index=False)
 

# This function fit the binary model on the PRS.
def fit_continous_phenotype_on_PRS(traindirec, newtrainfilename,p, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
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

            prs_train = pd.read_table(traindirec+os.sep+Name+os.sep+"train_hail_prs.csv",sep=",",index_col=0)
            print(prs_train)
            prs_train[['FID', 'IID']] = prs_train['subjectID'].str.split('_', expand=True)

            # Drop the original 'Combined' column if it's no longer needed
            prs_train = prs_train.drop(columns=['subjectID'])
            prs_train.rename(columns={'PGS002724': 'SCORE'}, inplace=True)

            

            prs_train['FID'] = prs_train['FID'].astype(str)
            prs_train['IID'] = prs_train['IID'].astype(str)
         
             
            prs_test = pd.read_table(traindirec+os.sep+Name+os.sep+"test_hail_prs.csv",sep=",",index_col=0)
            
            prs_test[['FID', 'IID']] = prs_test['subjectID'].str.split('_', expand=True)

            # Drop the original 'Combined' column if it's no longer needed
            prs_test = prs_test.drop(columns=['subjectID'])
            prs_test.rename(columns={'PGS002724': 'SCORE'}, inplace=True)
   
            print(prs_test)
        
            prs_test['FID'] = prs_test['FID'].astype(str)
            prs_test['IID'] = prs_test['IID'].astype(str)
            pheno_prs_train = pd.merge(covandpcs_train, prs_train, on=["FID", "IID"])
            pheno_prs_test = pd.merge(covandpcs_test, prs_test, on=["FID", "IID"])
            print(pheno_prs_train)
            
            

            try:
                #model = sm.OLS(phenotype_train["Phenotype"], sm.add_constant(pheno_prs_train.iloc[:, 2:])).fit_regularized(alpha=tempalpha, L1_wt=l1weight)
                model = sm.OLS(phenotype_train["Phenotype"], sm.add_constant(pheno_prs_train.iloc[:, 2:])).fit()

            except:
                print("Model did not fit")
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
             
                "numberofpca":p, 

                "tempalpha":str(tempalpha),
                "l1weight":str(l1weight),
                

                "Train_pure_prs":explained_variance_score(phenotype_train["Phenotype"],prs_train['SCORE'].values),
                "Train_null_model":explained_variance_score(phenotype_train["Phenotype"],train_null_predicted),
                "Train_best_model":explained_variance_score(phenotype_train["Phenotype"],train_best_predicted),

                "Test_pure_prs":explained_variance_score(phenotype_test["Phenotype"],prs_test['SCORE'].values),
                "Test_null_model":explained_variance_score(phenotype_test["Phenotype"],test_null_predicted),
                "Test_best_model":explained_variance_score(phenotype_test["Phenotype"],test_best_predicted),

            }, ignore_index=True)

            print(prs_result)
            prs_result.to_csv(traindirec+os.sep+Name+os.sep+"Results.csv",index=False)

    return


# ## Execute HAIL

# In[23]:


import hail as hl
import shutil 
import os
# Define a global variable to store results
prs_result = pd.DataFrame()
def transform_hail_data(traindirec, newtrainfilename,p, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
    ### First perform clumping on the file and save the clumpled file.
    perform_clumping_and_pruning_on_individual_data(traindirec, newtrainfilename,p, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile)
    
    #newtrainfilename = newtrainfilename+".clumped.pruned"
    #testfilename = testfilename+".clumped.pruned"
    
    
    #clupmedfile = traindirec+os.sep+newtrainfilename+".clump"
    #prunedfile = traindirec+os.sep+newtrainfilename+".clumped.pruned"

        
    # Also extract the PCA at this point for both test and training data.
    calculate_pca_for_traindata_testdata_for_clumped_pruned_snps(traindirec, newtrainfilename,p)

 
    # Delete the files generated in the previous iteration
 

    # Loop through each chromosome and delete the corresponding files
    import time
    
    def remove_path(path):
        try:
            if os.path.isdir(path):
                shutil.rmtree(path)  # Remove directory
                print(f'Deleted directory: {path}')
            else:
                os.remove(path)  # Remove file
                print(f'Deleted file: {path}')
        except OSError as e:
            print(f'Error deleting {path}: {e}')
            print('Retrying...')
          
        

    # Loop through each chromosome and delete the corresponding files and directories
    for chromosome in range(1, 23):
        train_path = os.path.join(traindirec, f'hail.train.{chromosome}.mt')
        test_path = os.path.join(traindirec, f'hail.test.{chromosome}.mt')

        # Delete train path if it exists
        if os.path.exists(train_path):
            remove_path(train_path)
            pass
        
        # Delete test path if it exists
        if os.path.exists(test_path):
            remove_path(test_path)
            pass
        
    

    #"""
    # Here we will perform processing on the genotype data and convert it to specific format
    # as required by HAIL.
    # Read the genotype data
    # Convert it to VCF
    # Use beagle to convert the data to  beagle format
    # Convert the beagle format to Hail format.
    
    
    
    
    for chromosome in range(1, 23):
        # Plink command to split by chromosome
        plink_command = [
            "./plink",
            "--bfile", traindirec+os.sep+newtrainfilename+".clumped.pruned",
            
            "--chr", str(chromosome),
            "--make-bed",
            "--out", traindirec+os.sep+newtrainfilename+".clumped.pruned."+str(chromosome),
        ]
        try:
            subprocess.run(plink_command, check=True)
        except:
            pass
        
        # Convert chromosomes to VCF format.
        plink_command = [
            "./plink",
            "--bfile", traindirec+os.sep+newtrainfilename+".clumped.pruned."+str(chromosome),
            "--chr", str(chromosome),
            "--recode","vcf",
            "--out", traindirec+os.sep+"hail.train."+str(chromosome),
        ]
        try:
            subprocess.run(plink_command)
        except:
            pass
        
        # Run beagle and convert the data to beagle format.
        #A “DR2” subfield with the estimated squared correlation between the estimated allele
        #dose and the true allele dose
        #• An “AF” subfield with the estimated alternate allele frequencies in the target samples
        #• The “IMP” flag if the marker is imputed
        
        beagle_command = [
            "java", 
            "-Xmx50g", 
            "-jar", "beagle", 
            #"gp=true",
            #"impute=false",
            #"burnin=1",
            #"window=100",
            "gt="+traindirec+os.sep+"hail.train."+str(chromosome)+".vcf",
            "ref="+"ALL.chr"+str(chromosome)+".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
            "map="+"plink.chr"+str(chromosome)+".GRCh37.map",
            "out="+traindirec+os.sep+"beagle.hail.train."+str(chromosome)
        ]
        
        # Execute the command
        try:
            subprocess.run(beagle_command)
            print(" ".join(beagle_command))
        except:
            pass
        
        
        # Convert beagle format to HAIL format.
        try:
            hl.import_vcf(traindirec+os.sep+"beagle.hail.train."+str(chromosome)+".vcf.gz", force_bgz=True).write(traindirec+os.sep+"hail.train."+str(chromosome)+".mt", overwrite=True)
        except:
            pass
        
    # Read the genotype data       
 
    mt = ""
     
    filecount = 0
    for chromosome in range(1, 23):
        file = os.path.join(traindirec, f"hail.train.{chromosome}.mt")

        if os.path.exists(file):  # Check if the file exists
            if filecount == 0:
                mt = hl.read_matrix_table(file)
                print(mt.rows().show(5))
                print(chromosome, file, mt.count(), mt.count()) 
                filecount += 1  # Increment filecount after the first file is processed
            else:
                tmpmt = hl.read_matrix_table(file)
                mt = mt.union_rows(tmpmt)
                print(chromosome, file, tmpmt.count(), mt.count())
        else:
            print("File does not exist:", file)
            
 
    
    print(mt.count())

    mt.write(traindirec+os.sep+'HAILTRAIN.mt', overwrite=True)
    mt = hl.read_matrix_table(traindirec+os.sep+'HAILTRAIN.mt')
    print(mt.rows().show(5))
    
    mt = mt.annotate_rows(variantID = (hl.str(mt.locus.contig) + ":" + hl.str(mt.locus.position)) )    
    
    # We skipped the cleaning process as the data is already clean.
    #print(mt.rows().show(5))
    #mt1 = mt.filter_rows(hl.len(mt.info.DR2) == 1)  
    #mtnot1 = mt.filter_rows(hl.len(mt.info.DR2) > 1)
    #mt1_filt = mt1.filter_rows(mt1.info.DR2.first()>=0.3)
    #print(mt1_filt.rows().show(5))   
    
    model_PGS002724  = hl.import_table(filedirec+os.sep+"Hail.txt", impute=True, force=True, comment='#')
    model_PGS002724 = model_PGS002724.annotate(
    variantID = hl.str(model_PGS002724.chr_name) + ":" + hl.str(model_PGS002724.chr_position) 
    )
    model_PGS002724 = model_PGS002724.key_by('variantID')
    mt_match_PGS002724 = mt.annotate_rows(**model_PGS002724[mt.variantID])
    mt_match_PGS002724 = mt_match_PGS002724.filter_rows(hl.is_defined(mt_match_PGS002724.effect_weight))
    flip_PGS002724 = hl.case().when( 
        (mt_match_PGS002724.effect_allele == mt_match_PGS002724.alleles[0])
        & (mt_match_PGS002724.other_allele == mt_match_PGS002724.alleles[1]), True ).when( 
        (mt_match_PGS002724.effect_allele == mt_match_PGS002724.alleles[1])
        & (mt_match_PGS002724.other_allele == mt_match_PGS002724.alleles[0]), False ).or_missing()
    
    mt_match_PGS002724 = mt_match_PGS002724.annotate_rows(flip=flip_PGS002724)
    prs_PGS002724 = hl.agg.sum(hl.float64(mt_match_PGS002724.effect_weight) * 
                    hl.if_else( mt_match_PGS002724.flip, 
                                2 - mt_match_PGS002724.DS.first(),
                                mt_match_PGS002724.DS.first()))
     
    mt_match_PGS002724 = mt_match_PGS002724.annotate_cols(prs=prs_PGS002724)
    mt_match_PGS002724.cols().export(traindirec+os.sep+Name+os.sep+'train_PRS.txt')
    prs_PGS002724 = hl.import_table(traindirec+os.sep+Name+os.sep+'train_PRS.txt', impute=True, force=True)
   
    prs_PGS002724 = prs_PGS002724.key_by('s')
    prs_merge = prs_PGS002724.rename({'s':'subjectID', 'prs':'PGS002724'}) 
    prs_merge_pandas = prs_merge.to_pandas()
    print(prs_merge_pandas.head())
    prs_merge_pandas.to_csv(traindirec+os.sep+Name+os.sep+"train_hail_prs.csv")
     
    # Save the data.
    prs_merge_pandas = pd.read_csv(traindirec+os.sep+Name+os.sep+"train_hail_prs.csv",index_col=0)
    print(prs_merge_pandas.head())
    #"""
    
    #"""
    
    # Repeat the process for test dataset.
    # testfilename
    
    for chromosome in range(1, 23):
        # Plink command to split by chromosome
        plink_command = [
            "./plink",
            "--bfile", traindirec+os.sep+testfilename+".clumped.pruned",
 
            "--chr", str(chromosome),
            "--make-bed",
            "--out", traindirec+os.sep+testfilename+".clumped.pruned."+str(chromosome),
        ]
        try:
            subprocess.run(plink_command)
        except:
            pass
        
        plink_command = [
            "./plink",
            "--bfile",traindirec+os.sep+testfilename+".clumped.pruned."+str(chromosome),
            "--chr", str(chromosome),
            "--recode","vcf",
            "--out", traindirec+os.sep+"hail.test."+str(chromosome),
        ]
        try:
            subprocess.run(plink_command, check=True)
        except:
            pass
        
        beagle_command = [
            "java", 
            "-Xmx50g", 
            "-jar", "beagle", 
            #"gp=true",
            #"burnin=1",
            #"window=100",
            "gt="+traindirec+os.sep+"hail.test."+str(chromosome)+".vcf",
            "ref="+"ALL.chr"+str(chromosome)+".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
            "map="+"plink.chr"+str(chromosome)+".GRCh37.map",
            "out="+traindirec+os.sep+"beagle.hail.test."+str(chromosome)
        ]

        # Execute the command
        try:
            subprocess.run(beagle_command)
            print(" ".join(beagle_command))
        except:
            pass
        
        #raise
        try:
            hl.import_vcf(traindirec+os.sep+"beagle.hail.test."+str(chromosome)+".vcf.gz", force_bgz=True).write(traindirec+os.sep+"hail.test."+str(chromosome)+".mt", overwrite=True)
            #raise
        except:
            pass
 
    mt = ""
     
    filecount = 0
    for chromosome in range(1, 23):
        file = os.path.join(traindirec, f"hail.test.{chromosome}.mt")

        if os.path.exists(file):  # Check if the file exists
            if filecount == 0:
                mt = hl.read_matrix_table(file)
                print(mt.rows().show(5))
                print(chromosome, file, mt.count(), mt.count()) 
                filecount += 1  # Increment filecount after the first file is processed
            else:
                tmpmt = hl.read_matrix_table(file)
                mt = mt.union_rows(tmpmt)
                print(chromosome, file, tmpmt.count(), mt.count())
        else:
            print("File does not exist:", file)    
    
    print(mt.count())
    mt.write(traindirec+os.sep+'HAILTEST.mt', overwrite=True)
    
    #"""
    #convert_data_in_hail_format()
    #"""
    mt = hl.read_matrix_table(traindirec+os.sep+'HAILTEST.mt')
    print(mt.rows().show(5))
    
    mt = mt.annotate_rows(variantID = (hl.str(mt.locus.contig) + ":" + hl.str(mt.locus.position)) )    
    #print(mt.rows().show(5))
    #mt1 = mt.filter_rows(hl.len(mt.info.DR2) == 1)  
    #mtnot1 = mt.filter_rows(hl.len(mt.info.DR2) > 1)
    #mt1_filt = mt1.filter_rows(mt1.info.DR2.first()>=0.3)
    #print(mt1_filt.rows().show(5))   
    
    model_PGS002724  = hl.import_table(filedirec+os.sep+"Hail.txt", impute=True, force=True, comment='#')
    model_PGS002724 = model_PGS002724.annotate(
    variantID = hl.str(model_PGS002724.chr_name) + ":" + hl.str(model_PGS002724.chr_position) 
    )
    model_PGS002724 = model_PGS002724.key_by('variantID')
    mt_match_PGS002724 = mt.annotate_rows(**model_PGS002724[mt.variantID])
    mt_match_PGS002724 = mt_match_PGS002724.filter_rows(hl.is_defined(mt_match_PGS002724.effect_weight))
    flip_PGS002724 = hl.case().when( 
        (mt_match_PGS002724.effect_allele == mt_match_PGS002724.alleles[0])
        & (mt_match_PGS002724.other_allele == mt_match_PGS002724.alleles[1]), True ).when( 
        (mt_match_PGS002724.effect_allele == mt_match_PGS002724.alleles[1])
        & (mt_match_PGS002724.other_allele == mt_match_PGS002724.alleles[0]), False ).or_missing()
    
    mt_match_PGS002724 = mt_match_PGS002724.annotate_rows(flip=flip_PGS002724)
    prs_PGS002724 = hl.agg.sum(hl.float64(mt_match_PGS002724.effect_weight) * 
                    hl.if_else( mt_match_PGS002724.flip, 
                                2 - mt_match_PGS002724.DS.first(),
                                mt_match_PGS002724.DS.first()))
     
    mt_match_PGS002724 = mt_match_PGS002724.annotate_cols(prs=prs_PGS002724)
    mt_match_PGS002724.cols().export(traindirec+os.sep+Name+os.sep+'test_PRS.txt')
    prs_PGS002724 = hl.import_table(traindirec+os.sep+Name+os.sep+'test_PRS.txt', impute=True, force=True)
   
    prs_PGS002724 = prs_PGS002724.key_by('s')
    prs_merge = prs_PGS002724.rename({'s':'subjectID', 'prs':'PGS002724'}) 
    prs_merge_pandas = prs_merge.to_pandas()
    print(prs_merge_pandas.head())
    prs_merge_pandas.to_csv(traindirec+os.sep+Name+os.sep+"test_hail_prs.csv")
    
    prs_merge_pandas = pd.read_csv(traindirec+os.sep+Name+os.sep+"test_hail_prs.csv",index_col=0)
    print(prs_merge_pandas.head())
    #"""
        
    
     
 

    # Load the PCA and Load the Covariates for trainingdatafirst.
    
    if check_phenotype_is_binary_or_continous(filedirec)=="Binary":
        print("Binary Phenotype!")
        fit_binary_phenotype_on_PRS(traindirec, newtrainfilename,p, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile)
    else:
        print("Continous Phenotype!")
        fit_continous_phenotype_on_PRS(traindirec, newtrainfilename,p, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile)
            
 


 
result_directory = "HAIL"
# Nested loops to iterate over different parameter values
create_directory(folddirec+os.sep+result_directory)
for p1_val in p_window_size:
 for p2_val in p_slide_size: 
  for p3_val in p_LD_threshold:
   for c1_val in clump_p1:
    for c2_val in clump_r2:
     for c3_val in clump_kb:
      for p in numberofpca:
        
        transform_hail_data(folddirec, newtrainfilename, p, str(p1_val), str(p2_val), str(p3_val), str(c1_val), str(c2_val), str(c3_val), result_directory, pvaluefile)


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
# python HAIL.py 0
# python HAIL.py 1
# python HAIL.py 2
# python HAIL.py 3
# python HAIL.py 4
# ```
# 
# The following files should exist after the execution:
# 
# 1. `SampleData1/Fold_0/HAIL/Results.csv`
# 2. `SampleData1/Fold_1/HAIL/Results.csv`
# 3. `SampleData1/Fold_2/HAIL/Results.csv`
# 4. `SampleData1/Fold_3/HAIL/Results.csv`
# 5. `SampleData1/Fold_4/HAIL/Results.csv`
# 

# ### Check the results file for each fold.

# In[8]:


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

# In[9]:


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
        'PRSice-2_Model',
        'effectsizes',
        'gemmamodel',
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

# In[10]:


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




