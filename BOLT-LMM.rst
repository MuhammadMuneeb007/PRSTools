BOLT-LMM
========

In this notebook, we will use BOLT to calculate the Polygenic Risk Score
(PRS).

Installation
------------

**Note:** BOLT needs to be installed or placed in the same directory as
this notebook.

1. Download BOLT and extract the files:

   .. code:: bash

      wget https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.4.1.tar.gz
      tar -xvf BOLT-LMM_v2.4.1.tar.gz

2. Copy all BOLT files to the current working directory:

   .. code:: bash

      cd BOLT-LMM_v2.4.1/
      cp -r * ../

**Documentation:** Available at `BOLT-LMM
Manual <https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.4.1_manual.pdf>`__

GWAS for BOLT
-------------

BOLT does not accept the GWAS file directly. Instead, it computes the
GWAS file using the genotype data.

**Note:** When using BOLT, make sure the covariate files use a specific
kind of prefix. If you do not specify a covariate file, the code still
works; you can comment out the covariate line. We renamed the covariate
and made them consistent with the prefix, as required by BOLT.

Another important point is that BOLT does not use the GWAS summary file.
Instead, it uses the individual genotype data to estimate effect sizes.
These effect sizes are then used to calculate PRS using Plink for both
training and test data.

GWAS file processing for BOLT for Binary Phenotypes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the effect size relates to disease risk and is thus given as an
odds ratio (OR) rather than BETA (for continuous traits), the PRS is
computed as a product of ORs. To simplify this calculation, take the
natural logarithm of the OR so that the PRS can be computed using
summation instead.

Sample
------

When executing BOLT, specify the path to the conda environment:

.. code:: python

   plink_command = [
       './bolt',
       '--bfile=' + traindirec + os.sep + newtrainfilename + ".clumped.pruned",
       '--phenoFile=' + traindirec + os.sep + trainfilename + ".PHENO",
       '--phenoCol=PHENO',
       mm,
       '--LDscoresFile=tables/LDSCORE.1000G_EUR.tab.gz',
       '--covarFile=' + traindirec + os.sep + trainfilename + ".COV_PCA",
       #'--covarCol=' + "Sex",
       '--qCovarCol=COV_{1:' + str(len(columns) - 2) + '}',
       #'--statsFile=' + filedirec + os.sep + filedirec2 + "." + mm.replace("-", ""),
       '--statsFile=' + traindirec + os.sep + filedirec + "." + mm.replace("-", "") + "_stat",
       #'--predBetasFile=' + filedirec + os.sep + filedirec + "." + mm.replace("-", "") + "_pred"
   ]

   print(" ".join(plink_command))

   os.system("LD_LIBRARY_PATH=/data/ascher01/uqmmune1/miniconda3/envs/genetics/lib/ " + " ".join(plink_command))

Possible error
--------------

::

   ERROR: Heritability estimate is close to 0; LMM may not correct confounding
          Instead, use PC-corrected linear/logistic regression on unrelateds

Possible solution
~~~~~~~~~~~~~~~~~

Include more samples and remove related individuals using Plink.

.. code:: ipython3

    print("GWAS file is not required when using BOLt! It generates it's own GWAS file include BETA estimates!")


.. parsed-literal::

    GWAS file is not required when using BOLt! It generates it's own GWAS file include BETA estimates!
    

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
    def check_phenotype_is_binary_or_continous(filedirec):
        # Read the processed quality controlled file for a phenotype
        df = pd.read_csv(filedirec+os.sep+filedirec+'_QC.fam',sep="\s+",header=None)
        column_values = df[5].unique()
     
        if len(set(column_values)) == 2:
            return "Binary"
        else:
            return "Continous"
    def create_directory(directory):
        """Function to create a directory if it doesn't exist."""
        if not os.path.exists(directory):  # Checking if the directory doesn't exist
            os.makedirs(directory)  # Creating the directory if it doesn't exist
        return directory  # Returning the created or existing directory
    
    #filedirec = sys.argv[1]
    filedirec = "SampleData1"
    
    #foldnumber = sys.argv[2]
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
    prs_result = pd.DataFrame(columns=["clump_p1", "clump_r2", "clump_kb", "p_window_size", "p_slide_size", "p_LD_threshold",
                                       "pvalue", "numberofpca","numberofvariants","Train_pure_prs", "Train_null_model", "Train_best_model",
                                       "Test_pure_prs", "Test_null_model", "Test_best_model"])

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
    
    
    def perform_clumping_and_pruning_on_individual_data(traindirec, newtrainfilename,numberofpca, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
        
        command = [
        "./plink",
        "--maf", "0.0001",
    
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
    def fit_binary_phenotype_on_PRS(traindirec, newtrainfilename,p,mm, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
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
                        "numberofvariants": len(pd.read_csv(traindirec+os.sep+newtrainfilename+".clumped.pruned.bim")),
                        
                         "BOLTmodel":mm,
    
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
    def fit_continous_phenotype_on_PRS(traindirec, newtrainfilename,p,mm, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
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
                         
                        
                        "BOLTmodel":mm,
                        "Train_pure_prs":explained_variance_score(phenotype_train["Phenotype"],prs_train['SCORE'].values),
                        "Train_null_model":explained_variance_score(phenotype_train["Phenotype"],train_null_predicted),
                        "Train_best_model":explained_variance_score(phenotype_train["Phenotype"],train_best_predicted),
                        
                        "Test_pure_prs":explained_variance_score(phenotype_test["Phenotype"],prs_test['SCORE'].values),
                        "Test_null_model":explained_variance_score(phenotype_test["Phenotype"],test_null_predicted),
                        "Test_best_model":explained_variance_score(phenotype_test["Phenotype"],test_best_predicted),
                        
                    }, ignore_index=True)
    
              
                    prs_result.to_csv(traindirec+os.sep+Name+os.sep+"Results.csv",index=False)
         
        return
    

Execute BOLT-LMM
----------------

.. code:: ipython3

    
    def transform_bolt_lmm_data(traindirec, newtrainfilename,mm,numberofpca, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
        ### First perform clumping on the file and save the clumpled file.
        #perform_clumping_and_pruning_on_individual_data(traindirec, newtrainfilename,p, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile)
        
        #newtrainfilename = newtrainfilename+".clumped.pruned"
        #testfilename = testfilename+".clumped.pruned"
        
        
        #clupmedfile = traindirec+os.sep+newtrainfilename+".clump"
        #prunedfile = traindirec+os.sep+newtrainfilename+".clumped.pruned"
    
            
        # Also extract the PCA at this point for both test and training data.
        #calculate_pca_for_traindata_testdata_for_clumped_pruned_snps(traindirec, newtrainfilename,p)
    
        #Extract p-values from the GWAS file.
        # Command for Linux.
        #os.system("awk "+"\'"+"{print $3,$8}"+"\'"+" ./"+filedirec+os.sep+filedirec+".txt >  ./"+traindirec+os.sep+"SNP.pvalue")
    
        # Command for windows.
        ### For windows get GWAK.
        ### https://sourceforge.net/projects/gnuwin32/
        ##3 Get it and place it in the same direc.
        #os.system("gawk "+"\""+"{print $3,$8}"+"\""+" ./"+filedirec+os.sep+filedirec+".txt >  ./"+traindirec+os.sep+"SNP.pvalue")
        #print("gawk "+"\""+"{print $3,$8}"+"\""+" ./"+filedirec+os.sep+filedirec+".txt >  ./"+traindirec+os.sep+"SNP.pvalue")
    
        #exit(0)
        # Delete files generated in the previous iteration.
        files_to_remove = [
            traindirec+os.sep+filedirec+"."+mm.replace("-","")+"_stat",
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
                
        
        
        
        tempphenotype_train = pd.read_table(traindirec+os.sep+newtrainfilename+".clumped.pruned"+".fam", sep="\s+",header=None)
        phenotype = pd.DataFrame()
        phenotype = tempphenotype_train[[0,1,5]]
        phenotype.to_csv(traindirec+os.sep+trainfilename+".PHENO",sep="\t",header=['FID', 'IID', 'PHENO'],index=False)
     
        pcs_train = pd.read_table(traindirec+os.sep+trainfilename+".eigenvec", sep="\s+",header=None, names=["FID", "IID"] + [f"PC{str(i)}" for i in range(1, int(p)+1)])
        covariate_train = pd.read_table(traindirec+os.sep+trainfilename+".cov",sep="\s+")
        covariate_train.fillna(0, inplace=True)
        print(covariate_train.head())
        print(len(covariate_train))
        covariate_train = covariate_train[covariate_train["FID"].isin(pcs_train["FID"].values) & covariate_train["IID"].isin(pcs_train["IID"].values)]
        print(len(covariate_train))
     
        covariate_train['FID'] = covariate_train['FID'].astype(str)
        pcs_train['FID'] = pcs_train['FID'].astype(str)
        covariate_train['IID'] = covariate_train['IID'].astype(str)
        pcs_train['IID'] = pcs_train['IID'].astype(str)
        covandpcs_train = pd.merge(covariate_train, pcs_train, on=["FID","IID"])
        covandpcs_train.fillna(0, inplace=True)    
    
        print(covandpcs_train)
        
        # Here ensure that the first and second column is FID and IID. For Bolt we need the covariates to start from the 
        # specific prefix like COV_X
        original_columns = covandpcs_train.columns
    
        # Rename columns other than 'FID' and 'IID'
        new_columns = ['FID', 'IID'] + [f'COV_{i+1}' for i in range(len(original_columns) - 2)]
    
        # Create a mapping of old column names to new names
        rename_mapping = dict(zip(original_columns, new_columns))
    
        # Rename the columns
        covandpcs_train.rename(columns=rename_mapping, inplace=True)
    
        # Reorder columns to ensure 'FID' and 'IID' are first
        columns = ['FID', 'IID'] + [f'COV_{i+1}' for i in range(len(original_columns) - 2)]
        covandpcs_train = covandpcs_train[columns]
        covandpcs_train.to_csv(traindirec+os.sep+trainfilename+".COV_PCA",sep="\t",index=False)    
        
        
        command = [
            
        './bolt',
        '--bfile='+traindirec+os.sep+newtrainfilename+".clumped.pruned",
        '--phenoFile='+traindirec+os.sep+trainfilename+".PHENO" ,
        '--phenoCol=PHENO',
        mm,
        '--LDscoresFile=tables/LDSCORE.1000G_EUR.tab.gz',
        '--covarFile='+traindirec+os.sep+trainfilename+".COV_PCA",
        #'--covarCol='+"COV_1",
        # TO include the first covariate which is sex use the following code.
        # Here i assumed that the first covariate is the sex. For our data the first covariate is sex.
        
        #
        #ERROR: Heritability estimate is close to 0; LMM may not correct confounding
        #   Instead, use PC-corrected linear/logistic regression on unrelateds
        #ERROR: Heritability estimate is close to 1; LMM may not correct confounding
        #   Instead, use PC-corrected linear/logistic regression on unrelateds
        
        #'--qCovarCol=COV_{1:'+str(len(columns)-len(columns)+1)+'}',
        
        # To include all the covariate use the following code, but not that it may crash the code as the heritability
        # from the geneotype data may reach to 0 and the BOLT-LMM may not work.
        # If heriability is close 0 or close to 1 the BOLT-LMM may not work.
        '--qCovarCol=COV_{1:'+str(len(columns)-4)+'}',
            
        
        #'--statsFile='+filedirec+os.sep+filedirec2+"."+mm.replace("-","")
        '--statsFile='+traindirec+os.sep+filedirec+"."+mm.replace("-","")+"_stat",
        #'--predBetasFile='+filedirec+os.sep+filedirec+"."+mm.replace("-","")+"_pred"
        ]
        print(" ".join(command))
        
    
        os.system("LD_LIBRARY_PATH=/data/ascher01/uqmmune1/miniconda3/envs/genetics/lib/ "+" ".join(command))
        
        #return
        gwas = pd.read_csv(traindirec+os.sep+filedirec+"."+mm.replace("-","")+"_stat",sep="\s+")
        print(gwas.head())
        
        if check_phenotype_is_binary_or_continous(filedirec)=="Binary":
            gwas["BETA"] = np.exp(gwas["BETA"])
        else:
            pass
    
        gwas.iloc[:,[0,4,8]].to_csv(traindirec+os.sep+filedirec+"."+mm.replace("-","")+"_stat",sep="\t",index=False)
         
        
        command = [
            "./plink",
            "--bfile", traindirec+os.sep+newtrainfilename+".clumped.pruned",
            "--score", traindirec+os.sep+filedirec+"."+mm.replace("-","")+"_stat", "1", "2", "3", "header",
            "--q-score-range", traindirec+os.sep+"range_list",traindirec+os.sep+"SNP.pvalue",
            "--extract", traindirec+os.sep+trainfilename+".valid.snp",
            "--out", traindirec+os.sep+Name+os.sep+trainfilename
        ]
        subprocess.run(command)
    
    
        command = [
            "./plink",
            "--bfile", traindirec+os.sep+testfilename+".clumped.pruned",
            ### SNP column = 3, Effect allele column 1 = 4, Beta column=12
            "--score", traindirec+os.sep+filedirec+"."+mm.replace("-","")+"_stat", "1", "2", "3", "header",
            "--q-score-range", traindirec+os.sep+"range_list",traindirec+os.sep+"SNP.pvalue",
            "--extract", traindirec+os.sep+trainfilename+".valid.snp",
            "--out", folddirec+os.sep+Name+os.sep+testfilename
        ]
        subprocess.run(command)
        
        
        if check_phenotype_is_binary_or_continous(filedirec)=="Binary":
            print("Binary Phenotype!")
            fit_binary_phenotype_on_PRS(traindirec, newtrainfilename,p,mm.replace("-",""), p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile)
        else:
            print("Continous Phenotype!")
            fit_continous_phenotype_on_PRS(traindirec, newtrainfilename,p,mm.replace("-",""), p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile)
                
     
    
     
     
     
    result_directory = "BOLT-LMM"
    models = ["--lmm","--lmmInfOnly","--lmmForceNonInf"]
    # Nested loops to iterate over different parameter values
    create_directory(folddirec+os.sep+result_directory)
    for p1_val in p_window_size:
     for p2_val in p_slide_size: 
      for p3_val in p_LD_threshold:
       for c1_val in clump_p1:
        for c2_val in clump_r2:
         for c3_val in clump_kb:
          for p in numberofpca:
           for mm in models: 
            transform_bolt_lmm_data(folddirec, newtrainfilename,mm, p, str(p1_val), str(p2_val), str(p3_val), str(c1_val), str(c2_val), str(c3_val), result_directory, pvaluefile)
    
     


.. parsed-literal::

    Removed file: SampleData1/Fold_0/SampleData1.lmm_stat
           FID      IID  Sex
    0  HG00097  HG00097    2
    1  HG00099  HG00099    2
    2  HG00101  HG00101    1
    3  HG00102  HG00102    2
    4  HG00103  HG00103    1
    380
    380
             FID      IID  Sex       PC1       PC2       PC3       PC4       PC5  \
    0    HG00097  HG00097    2 -0.001453  0.084820  0.006792  0.013653  0.027149   
    1    HG00099  HG00099    2 -0.002017  0.089514 -0.022355  0.001888 -0.000037   
    2    HG00101  HG00101    1 -0.000380  0.096056 -0.018231 -0.016026  0.012093   
    3    HG00102  HG00102    2  0.000292  0.071832  0.018087 -0.045180  0.028123   
    4    HG00103  HG00103    1 -0.008372  0.065005 -0.009089 -0.026468 -0.009184   
    ..       ...      ...  ...       ...       ...       ...       ...       ...   
    375  NA20818  NA20818    2 -0.047156 -0.040644 -0.052693  0.021050 -0.013389   
    376  NA20826  NA20826    2 -0.042629 -0.059404 -0.066130  0.006495 -0.009525   
    377  NA20827  NA20827    1 -0.044060 -0.053125 -0.065463  0.015030 -0.004314   
    378  NA20828  NA20828    2 -0.047621 -0.050577 -0.043164  0.003004 -0.016823   
    379  NA20832  NA20832    2 -0.041535 -0.049826 -0.047877  0.005951 -0.003770   
    
              PC6  
    0    0.032581  
    1    0.009107  
    2    0.019296  
    3   -0.003620  
    4   -0.030565  
    ..        ...  
    375 -0.047403  
    376  0.010779  
    377  0.003873  
    378  0.015832  
    379 -0.023086  
    
    [380 rows x 9 columns]
    ./bolt --bfile=SampleData1/Fold_0/train_data.QC.clumped.pruned --phenoFile=SampleData1/Fold_0/train_data.PHENO --phenoCol=PHENO --lmm --LDscoresFile=tables/LDSCORE.1000G_EUR.tab.gz --covarFile=SampleData1/Fold_0/train_data.COV_PCA --qCovarCol=COV_{1:2} --statsFile=SampleData1/Fold_0/SampleData1.lmm_stat
                          +-----------------------------+
                          |                       ___   |
                          |   BOLT-LMM, v2.4.1   /_ /   |
                          |   November 16, 2022   /_/   |
                          |   Po-Ru Loh            //   |
                          |                        /    |
                          +-----------------------------+
    
    Copyright (C) 2014-2022 Harvard University.
    Distributed under the GNU GPLv3 open source license.
    
    Compiled with USE_SSE: fast aligned memory access
    Compiled with USE_MKL: Intel Math Kernel Library linear algebra
    Boost version: 1_58
    
    Command line options:
    
    ./bolt \
        --bfile=SampleData1/Fold_0/train_data.QC.clumped.pruned \
        --phenoFile=SampleData1/Fold_0/train_data.PHENO \
        --phenoCol=PHENO \
        --lmm \
        --LDscoresFile=tables/LDSCORE.1000G_EUR.tab.gz \
        --covarFile=SampleData1/Fold_0/train_data.COV_PCA \
        --qCovarCol=COV_{1:2} \
        --statsFile=SampleData1/Fold_0/SampleData1.lmm_stat 
    
    Setting number of threads to 1
    fam: SampleData1/Fold_0/train_data.QC.clumped.pruned.fam
    bim(s): SampleData1/Fold_0/train_data.QC.clumped.pruned.bim
    bed(s): SampleData1/Fold_0/train_data.QC.clumped.pruned.bed
    
    === Reading genotype data ===
    
    Total indivs in PLINK data: Nbed = 380
    Total indivs stored in memory: N = 380
    Reading bim file #1: SampleData1/Fold_0/train_data.QC.clumped.pruned.bim
        Read 172878 snps
    Total snps in PLINK data: Mbed = 172878
    
    Breakdown of SNP pre-filtering results:
      172878 SNPs to include in model (i.e., GRM)
      0 additional non-GRM SNPs loaded
      0 excluded SNPs
    Allocating 172878 x 380/4 bytes to store genotypes
    Reading genotypes and performing QC filtering on snps and indivs...
    Reading bed file #1: SampleData1/Fold_0/train_data.QC.clumped.pruned.bed
        Expecting 16423410 (+3) bytes for 380 indivs, 172878 snps
    

.. parsed-literal::

    WARNING: Genetic map appears to be in cM units; rescaling by 0.01
    

.. parsed-literal::

    Total indivs after QC: 380
    Total post-QC SNPs: M = 172878
      Variance component 1: 172878 post-QC SNPs (name: 'modelSnps')
    Time for SnpData setup = 1.51395 sec
    
    === Reading phenotype and covariate data ===
    
    Read data for 380 indivs (ignored 0 without genotypes) from:
      SampleData1/Fold_0/train_data.COV_PCA
    Read data for 380 indivs (ignored 0 without genotypes) from:
      SampleData1/Fold_0/train_data.PHENO
    Number of indivs with no missing phenotype(s) to use: 380
        Using quantitative covariate: COV_1
        Using quantitative covariate: COV_2
        Using quantitative covariate: CONST_ALL_ONES
    Number of individuals used in analysis: Nused = 380
    Singular values of covariate matrix:
        S[0] = 36.3836
        S[1] = 5.21943
        S[2] = 0.995369
    Total covariate vectors: C = 3
    Total independent covariate vectors: Cindep = 3
    
    === Initializing Bolt object: projecting and normalizing SNPs ===
    
    

.. parsed-literal::

    NOTE: Using all-1s vector (constant term) in addition to specified covariates
    

.. parsed-literal::

    Number of chroms with >= 1 good SNP: 22
    Average norm of projected SNPs:           375.206815
    Dimension of all-1s proj space (Nused-1): 379
    Time for covariate data setup + Bolt initialization = 0.778457 sec
    
    Phenotype 1:   N = 380   mean = 170.14   std = 0.945865
    
    === Computing linear regression (LINREG) stats ===
    
    Time for computing LINREG stats = 0.201925 sec
    
    === Estimating variance parameters ===
    
    Using CGtol of 0.005 for this step
    Using default number of random trials: 15 (for Nused = 380)
    
    Estimating MC scaling f_REML at log(delta) = 1.09384, h2 = 0.25...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.43  rNorms/orig: (0.02,0.02)  res2s: 882.719..149.989
      iter 2:  time=0.42  rNorms/orig: (0.0005,0.001)  res2s: 883.769..150.228
      Converged at iter 2: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 38.6%, memory/overhead = 61.4%
      MCscaling: logDelta = 1.09, h2 = 0.250, f = 0.00255218
    
    Estimating MC scaling f_REML at log(delta) = -0.00476781, h2 = 0.5...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.42  rNorms/orig: (0.04,0.05)  res2s: 206.958..66.4433
      iter 2:  time=0.42  rNorms/orig: (0.002,0.005)  res2s: 207.859..66.8373
      iter 3:  time=0.42  rNorms/orig: (9e-05,0.0003)  res2s: 207.871..66.8446
      Converged at iter 3: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 36.0%, memory/overhead = 64.0%
      MCscaling: logDelta = -0.00, h2 = 0.500, f = 0.000475493
    
    Estimating MC scaling f_REML at log(delta) = -0.256314, h2 = 0.562557...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.41  rNorms/orig: (0.04,0.05)  res2s: 142.182..50.8157
      iter 2:  time=0.37  rNorms/orig: (0.002,0.006)  res2s: 142.949..51.1908
      iter 3:  time=0.37  rNorms/orig: (0.0001,0.0005)  res2s: 142.962..51.1995
      Converged at iter 3: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 36.1%, memory/overhead = 63.9%
      MCscaling: logDelta = -0.26, h2 = 0.563, f = 5.80701e-05
    
    Estimating MC scaling f_REML at log(delta) = -0.291308, h2 = 0.571149...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.37  rNorms/orig: (0.04,0.05)  res2s: 134.763..48.8336
      iter 2:  time=0.38  rNorms/orig: (0.002,0.007)  res2s: 135.51..49.2044
      iter 3:  time=0.38  rNorms/orig: (0.0001,0.0005)  res2s: 135.523..49.2132
      Converged at iter 3: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 36.2%, memory/overhead = 63.8%
      MCscaling: logDelta = -0.29, h2 = 0.571, f = 3.23365e-06
    
    Secant iteration for h2 estimation converged in 2 steps
    Estimated (pseudo-)heritability: h2g = 0.571
    To more precisely estimate variance parameters and estimate s.e., use --reml
    Variance params: sigma^2_K = 0.406791, logDelta = -0.291308, f = 3.23365e-06
    
    Time for fitting variance components = 5.53623 sec
    
    === Computing mixed model assoc stats (inf. model) ===
    
    Selected 30 SNPs for computation of prospective stat
    Tried 30; threw out 0 with GRAMMAR chisq > 5
    Assigning SNPs to 22 chunks for leave-out analysis
    Each chunk is excluded when testing SNPs belonging to the chunk
      Batch-solving 52 systems of equations using conjugate gradient iteration
      iter 1:  time=0.70  rNorms/orig: (0.04,0.07)  res2s: 48.7902..70.329
      iter 2:  time=0.69  rNorms/orig: (0.002,0.007)  res2s: 49.1781..70.6337
      iter 3:  time=0.68  rNorms/orig: (9e-05,0.0005)  res2s: 49.1875..70.6346
      Converged at iter 3: rNorms/orig all < CGtol=0.0005
      Time breakdown: dgemm = 64.0%, memory/overhead = 36.0%
    
    AvgPro: 0.844   AvgRetro: 0.839   Calibration: 1.006 (0.003)   (30 SNPs)
    Ratio of medians: 1.013   Median of ratios: 1.002
    
    Time for computing infinitesimal model assoc stats = 2.26405 sec
    
    === Estimating chip LD Scores using 400 indivs ===
    
    Reducing sample size to 376 for memory alignment
    

.. parsed-literal::

    WARNING: Only 380 indivs available; using all
    

.. parsed-literal::

    
    Time for estimating chip LD Scores = 0.636621 sec
    
    === Reading LD Scores for calibration of Bayesian assoc stats ===
    
    Looking up LD Scores...
      Looking for column header 'SNP': column number = 1
      Looking for column header 'LDSCORE': column number = 5
    Found LD Scores for 171753/172878 SNPs
    
    Estimating inflation of LINREG chisq stats using MLMe as reference...
    Filtering to SNPs with chisq stats, LD Scores, and MAF > 0.01
    # of SNPs passing filters before outlier removal: 171753/172878
    Masking windows around outlier snps (chisq > 20.0)
    # of SNPs remaining after outlier window removal: 171753/171753
    Intercept of LD Score regression for ref stats:   1.003 (0.005)
    Estimated attenuation: 0.409 (0.534)
    Intercept of LD Score regression for cur stats: 1.004 (0.005)
    Calibration factor (ref/cur) to multiply by:      0.999 (0.000)
    LINREG intercept inflation = 1.00087
    
    === Estimating mixture parameters by cross-validation ===
    
    Setting maximum number of iterations to 250 for this step
    Max CV folds to compute = 5 (to have > 10000 samples)
    
    ====> Starting CV fold 1 <====
    
        Using quantitative covariate: COV_1
        Using quantitative covariate: COV_2
        Using quantitative covariate: CONST_ALL_ONES
    Number of individuals used in analysis: Nused = 304
    Singular values of covariate matrix:
        S[0] = 32.5614
        S[1] = 4.66512
        S[2] = 0.894787
    Total covariate vectors: C = 3
    Total independent covariate vectors: Cindep = 3
    
    === Initializing Bolt object: projecting and normalizing SNPs ===
    
    

.. parsed-literal::

    NOTE: Using all-1s vector (constant term) in addition to specified covariates
    

.. parsed-literal::

    Number of chroms with >= 1 good SNP: 22
    Average norm of projected SNPs:           299.575111
    Dimension of all-1s proj space (Nused-1): 303
      Beginning variational Bayes
      iter 1:  time=1.09 for 18 active reps
      iter 2:  time=0.80 for 18 active reps  approxLL diffs: (37.61,49.27)
      iter 3:  time=0.80 for 18 active reps  approxLL diffs: (1.23,1.36)
      iter 4:  time=0.81 for 18 active reps  approxLL diffs: (0.14,0.21)
      iter 5:  time=0.80 for 18 active reps  approxLL diffs: (0.01,0.01)
      iter 6:  time=0.19 for  1 active reps  approxLL diffs: (0.00,0.00)
      Converged at iter 6: approxLL diffs each have been < LLtol=0.01
      Time breakdown: dgemm = 25.7%, memory/overhead = 74.3%
    Computing predictions on left-out cross-validation fold
    Time for computing predictions = 0.452455 sec
    
    Average PVEs obtained by param pairs tested (high to low):
      f2=0.5, p=0.5: 0.006720
      f2=0.3, p=0.5: 0.006720
      f2=0.5, p=0.2: 0.006720
                ...
     f2=0.3, p=0.01: 0.006710
    
    Detailed CV fold results:
      Absolute prediction MSE baseline (covariates only): 0.586308
      Absolute prediction MSE using standard LMM:         0.582368
      Absolute prediction MSE, fold-best  f2=0.5, p=0.5:  0.582368
        Absolute pred MSE using   f2=0.5, p=0.5: 0.582368
        Absolute pred MSE using   f2=0.5, p=0.2: 0.582368
        Absolute pred MSE using   f2=0.5, p=0.1: 0.582368
        Absolute pred MSE using  f2=0.5, p=0.05: 0.582368
        Absolute pred MSE using  f2=0.5, p=0.02: 0.582369
        Absolute pred MSE using  f2=0.5, p=0.01: 0.582372
        Absolute pred MSE using   f2=0.3, p=0.5: 0.582368
        Absolute pred MSE using   f2=0.3, p=0.2: 0.582368
        Absolute pred MSE using   f2=0.3, p=0.1: 0.582368
        Absolute pred MSE using  f2=0.3, p=0.05: 0.582369
        Absolute pred MSE using  f2=0.3, p=0.02: 0.582370
        Absolute pred MSE using  f2=0.3, p=0.01: 0.582374
        Absolute pred MSE using   f2=0.1, p=0.5: 0.582368
        Absolute pred MSE using   f2=0.1, p=0.2: 0.582368
        Absolute pred MSE using   f2=0.1, p=0.1: 0.582368
        Absolute pred MSE using  f2=0.1, p=0.05: 0.582369
        Absolute pred MSE using  f2=0.1, p=0.02: 0.582371
        Absolute pred MSE using  f2=0.1, p=0.01: 0.582370
    
    ====> End CV fold 1: 18 remaining param pair(s) <====
    
    Estimated proportion of variance explained using inf model: 0.007
    Relative improvement in prediction MSE using non-inf model: 0.000
    
    Exiting CV: non-inf model does not substantially improve prediction
    Optimal mixture parameters according to CV: f2 = 0.5, p = 0.5
    Bayesian non-infinitesimal model does not fit substantially better
    => Not computing non-inf assoc stats (to override, use --lmmForceNonInf)
    
    Time for estimating mixture parameters = 32.8537 sec
    
    Calibration stats: mean and lambdaGC (over SNPs used in GRM)
      (note that both should be >1 because of polygenicity)
    Mean BOLT_LMM_INF: 1.0074 (172878 good SNPs)   lambdaGC: 1.01336
    
    === Streaming genotypes to compute and write assoc stats at all SNPs ===
    
    Time for streaming genotypes and writing output = 1.68225 sec
    
    Total elapsed time for analysis = 45.4673 sec
               SNP  CHR      BP    GENPOS ALLELE1 ALLELE0    A1FREQ  F_MISS  \
    0   rs79373928    1  801536  0.587220       G       T  0.014474     0.0   
    1    rs4970382    1  840753  0.620827       C       T  0.406579     0.0   
    2   rs13303222    1  849998  0.620827       A       G  0.196053     0.0   
    3   rs72631889    1  851390  0.620827       T       G  0.034210     0.0   
    4  rs192998324    1  862772  0.620827       G       A  0.027632     0.0   
    
           BETA        SE  P_BOLT_LMM_INF  
    0  0.015560  0.258427           0.950  
    1 -0.060199  0.059829           0.310  
    2 -0.006768  0.078287           0.930  
    3  0.315642  0.172246           0.067  
    4 -0.227920  0.190562           0.230  
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/BOLT-LMM/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --extract SampleData1/Fold_0/train_data.valid.snp
      --out SampleData1/Fold_0/BOLT-LMM/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/SampleData1.lmm_stat 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 172878 variants remaining.
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
    Results written to SampleData1/Fold_0/BOLT-LMM/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/BOLT-LMM/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --extract SampleData1/Fold_0/train_data.valid.snp
      --out SampleData1/Fold_0/BOLT-LMM/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/SampleData1.lmm_stat 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    --extract: 172878 variants remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    172878 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/BOLT-LMM/test_data.*.profile.
    Continous Phenotype!
    

.. parsed-literal::

    /tmp/ipykernel_1257867/115566430.py:347: FutureWarning: The behavior of DataFrame concatenation with empty or all-NA entries is deprecated. In a future version, this will no longer exclude empty or all-NA columns when determining the result dtypes. To retain the old behavior, exclude the relevant entries before the concat operation.
      prs_result = prs_result._append({
    

.. parsed-literal::

    Removed file: SampleData1/Fold_0/SampleData1.lmmInfOnly_stat
           FID      IID  Sex
    0  HG00097  HG00097    2
    1  HG00099  HG00099    2
    2  HG00101  HG00101    1
    3  HG00102  HG00102    2
    4  HG00103  HG00103    1
    380
    380
             FID      IID  Sex       PC1       PC2       PC3       PC4       PC5  \
    0    HG00097  HG00097    2 -0.001453  0.084820  0.006792  0.013653  0.027149   
    1    HG00099  HG00099    2 -0.002017  0.089514 -0.022355  0.001888 -0.000037   
    2    HG00101  HG00101    1 -0.000380  0.096056 -0.018231 -0.016026  0.012093   
    3    HG00102  HG00102    2  0.000292  0.071832  0.018087 -0.045180  0.028123   
    4    HG00103  HG00103    1 -0.008372  0.065005 -0.009089 -0.026468 -0.009184   
    ..       ...      ...  ...       ...       ...       ...       ...       ...   
    375  NA20818  NA20818    2 -0.047156 -0.040644 -0.052693  0.021050 -0.013389   
    376  NA20826  NA20826    2 -0.042629 -0.059404 -0.066130  0.006495 -0.009525   
    377  NA20827  NA20827    1 -0.044060 -0.053125 -0.065463  0.015030 -0.004314   
    378  NA20828  NA20828    2 -0.047621 -0.050577 -0.043164  0.003004 -0.016823   
    379  NA20832  NA20832    2 -0.041535 -0.049826 -0.047877  0.005951 -0.003770   
    
              PC6  
    0    0.032581  
    1    0.009107  
    2    0.019296  
    3   -0.003620  
    4   -0.030565  
    ..        ...  
    375 -0.047403  
    376  0.010779  
    377  0.003873  
    378  0.015832  
    379 -0.023086  
    
    [380 rows x 9 columns]
    ./bolt --bfile=SampleData1/Fold_0/train_data.QC.clumped.pruned --phenoFile=SampleData1/Fold_0/train_data.PHENO --phenoCol=PHENO --lmmInfOnly --LDscoresFile=tables/LDSCORE.1000G_EUR.tab.gz --covarFile=SampleData1/Fold_0/train_data.COV_PCA --qCovarCol=COV_{1:2} --statsFile=SampleData1/Fold_0/SampleData1.lmmInfOnly_stat
                          +-----------------------------+
                          |                       ___   |
                          |   BOLT-LMM, v2.4.1   /_ /   |
                          |   November 16, 2022   /_/   |
                          |   Po-Ru Loh            //   |
                          |                        /    |
                          +-----------------------------+
    
    Copyright (C) 2014-2022 Harvard University.
    Distributed under the GNU GPLv3 open source license.
    
    Compiled with USE_SSE: fast aligned memory access
    Compiled with USE_MKL: Intel Math Kernel Library linear algebra
    Boost version: 1_58
    
    Command line options:
    
    ./bolt \
        --bfile=SampleData1/Fold_0/train_data.QC.clumped.pruned \
        --phenoFile=SampleData1/Fold_0/train_data.PHENO \
        --phenoCol=PHENO \
        --lmmInfOnly \
        --LDscoresFile=tables/LDSCORE.1000G_EUR.tab.gz \
        --covarFile=SampleData1/Fold_0/train_data.COV_PCA \
        --qCovarCol=COV_{1:2} \
        --statsFile=SampleData1/Fold_0/SampleData1.lmmInfOnly_stat 
    
    Setting number of threads to 1
    fam: SampleData1/Fold_0/train_data.QC.clumped.pruned.fam
    bim(s): SampleData1/Fold_0/train_data.QC.clumped.pruned.bim
    bed(s): SampleData1/Fold_0/train_data.QC.clumped.pruned.bed
    
    === Reading genotype data ===
    
    Total indivs in PLINK data: Nbed = 380
    Total indivs stored in memory: N = 380
    Reading bim file #1: SampleData1/Fold_0/train_data.QC.clumped.pruned.bim
        Read 172878 snps
    Total snps in PLINK data: Mbed = 172878
    
    Breakdown of SNP pre-filtering results:
      172878 SNPs to include in model (i.e., GRM)
      0 additional non-GRM SNPs loaded
      0 excluded SNPs
    Allocating 172878 x 380/4 bytes to store genotypes
    Reading genotypes and performing QC filtering on snps and indivs...
    Reading bed file #1: SampleData1/Fold_0/train_data.QC.clumped.pruned.bed
        Expecting 16423410 (+3) bytes for 380 indivs, 172878 snps
    

.. parsed-literal::

    WARNING: Genetic map appears to be in cM units; rescaling by 0.01
    

.. parsed-literal::

    Total indivs after QC: 380
    Total post-QC SNPs: M = 172878
      Variance component 1: 172878 post-QC SNPs (name: 'modelSnps')
    Time for SnpData setup = 1.19219 sec
    
    === Reading phenotype and covariate data ===
    
    Read data for 380 indivs (ignored 0 without genotypes) from:
      SampleData1/Fold_0/train_data.COV_PCA
    Read data for 380 indivs (ignored 0 without genotypes) from:
      SampleData1/Fold_0/train_data.PHENO
    Number of indivs with no missing phenotype(s) to use: 380
        Using quantitative covariate: COV_1
        Using quantitative covariate: COV_2
        Using quantitative covariate: CONST_ALL_ONES
    Number of individuals used in analysis: Nused = 380
    Singular values of covariate matrix:
        S[0] = 36.3836
        S[1] = 5.21943
        S[2] = 0.995369
    Total covariate vectors: C = 3
    Total independent covariate vectors: Cindep = 3
    
    === Initializing Bolt object: projecting and normalizing SNPs ===
    
    

.. parsed-literal::

    NOTE: Using all-1s vector (constant term) in addition to specified covariates
    

.. parsed-literal::

    Number of chroms with >= 1 good SNP: 22
    Average norm of projected SNPs:           375.206815
    Dimension of all-1s proj space (Nused-1): 379
    Time for covariate data setup + Bolt initialization = 0.375053 sec
    
    Phenotype 1:   N = 380   mean = 170.14   std = 0.945865
    
    === Computing linear regression (LINREG) stats ===
    
    Time for computing LINREG stats = 0.145117 sec
    
    === Estimating variance parameters ===
    
    Using CGtol of 0.005 for this step
    Using default number of random trials: 15 (for Nused = 380)
    
    Estimating MC scaling f_REML at log(delta) = 1.09384, h2 = 0.25...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.42  rNorms/orig: (0.02,0.02)  res2s: 882.719..149.989
      iter 2:  time=0.40  rNorms/orig: (0.0005,0.001)  res2s: 883.769..150.228
      Converged at iter 2: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 36.3%, memory/overhead = 63.7%
      MCscaling: logDelta = 1.09, h2 = 0.250, f = 0.00255218
    
    Estimating MC scaling f_REML at log(delta) = -0.00476781, h2 = 0.5...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.37  rNorms/orig: (0.04,0.05)  res2s: 206.958..66.4433
      iter 2:  time=0.37  rNorms/orig: (0.002,0.005)  res2s: 207.859..66.8373
      iter 3:  time=0.38  rNorms/orig: (9e-05,0.0003)  res2s: 207.871..66.8446
      Converged at iter 3: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 36.0%, memory/overhead = 64.0%
      MCscaling: logDelta = -0.00, h2 = 0.500, f = 0.000475493
    
    Estimating MC scaling f_REML at log(delta) = -0.256314, h2 = 0.562557...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.37  rNorms/orig: (0.04,0.05)  res2s: 142.182..50.8157
      iter 2:  time=0.37  rNorms/orig: (0.002,0.006)  res2s: 142.949..51.1908
      iter 3:  time=0.37  rNorms/orig: (0.0001,0.0005)  res2s: 142.962..51.1995
      Converged at iter 3: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 36.1%, memory/overhead = 63.9%
      MCscaling: logDelta = -0.26, h2 = 0.563, f = 5.80701e-05
    
    Estimating MC scaling f_REML at log(delta) = -0.291308, h2 = 0.571149...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.37  rNorms/orig: (0.04,0.05)  res2s: 134.763..48.8336
      iter 2:  time=0.37  rNorms/orig: (0.002,0.007)  res2s: 135.51..49.2044
      iter 3:  time=0.37  rNorms/orig: (0.0001,0.0005)  res2s: 135.523..49.2132
      Converged at iter 3: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 36.1%, memory/overhead = 63.9%
      MCscaling: logDelta = -0.29, h2 = 0.571, f = 3.23365e-06
    
    Secant iteration for h2 estimation converged in 2 steps
    Estimated (pseudo-)heritability: h2g = 0.571
    To more precisely estimate variance parameters and estimate s.e., use --reml
    Variance params: sigma^2_K = 0.406791, logDelta = -0.291308, f = 3.23365e-06
    
    Time for fitting variance components = 5.17613 sec
    
    === Computing mixed model assoc stats (inf. model) ===
    
    Selected 30 SNPs for computation of prospective stat
    Tried 30; threw out 0 with GRAMMAR chisq > 5
    Assigning SNPs to 22 chunks for leave-out analysis
    Each chunk is excluded when testing SNPs belonging to the chunk
      Batch-solving 52 systems of equations using conjugate gradient iteration
      iter 1:  time=0.68  rNorms/orig: (0.04,0.07)  res2s: 48.7902..70.329
      iter 2:  time=0.68  rNorms/orig: (0.002,0.007)  res2s: 49.1781..70.6337
      iter 3:  time=0.68  rNorms/orig: (9e-05,0.0005)  res2s: 49.1875..70.6346
      Converged at iter 3: rNorms/orig all < CGtol=0.0005
      Time breakdown: dgemm = 63.7%, memory/overhead = 36.3%
    
    AvgPro: 0.844   AvgRetro: 0.839   Calibration: 1.006 (0.003)   (30 SNPs)
    Ratio of medians: 1.013   Median of ratios: 1.002
    
    Time for computing infinitesimal model assoc stats = 2.25897 sec
    
    === Estimating chip LD Scores using 400 indivs ===
    
    Reducing sample size to 376 for memory alignment
    

.. parsed-literal::

    WARNING: Only 380 indivs available; using all
    

.. parsed-literal::

    
    Time for estimating chip LD Scores = 0.682689 sec
    
    === Reading LD Scores for calibration of Bayesian assoc stats ===
    
    Looking up LD Scores...
      Looking for column header 'SNP': column number = 1
      Looking for column header 'LDSCORE': column number = 5
    Found LD Scores for 171753/172878 SNPs
    
    Estimating inflation of LINREG chisq stats using MLMe as reference...
    Filtering to SNPs with chisq stats, LD Scores, and MAF > 0.01
    # of SNPs passing filters before outlier removal: 171753/172878
    Masking windows around outlier snps (chisq > 20.0)
    # of SNPs remaining after outlier window removal: 171753/171753
    Intercept of LD Score regression for ref stats:   1.003 (0.005)
    Estimated attenuation: 0.409 (0.534)
    Intercept of LD Score regression for cur stats: 1.004 (0.005)
    Calibration factor (ref/cur) to multiply by:      0.999 (0.000)
    LINREG intercept inflation = 1.00087
    
    Calibration stats: mean and lambdaGC (over SNPs used in GRM)
      (note that both should be >1 because of polygenicity)
    Mean BOLT_LMM_INF: 1.0074 (172878 good SNPs)   lambdaGC: 1.01336
    
    === Streaming genotypes to compute and write assoc stats at all SNPs ===
    
    Time for streaming genotypes and writing output = 27.6488 sec
    
    Total elapsed time for analysis = 37.4789 sec
               SNP  CHR      BP    GENPOS ALLELE1 ALLELE0    A1FREQ  F_MISS  \
    0   rs79373928    1  801536  0.587220       G       T  0.014474     0.0   
    1    rs4970382    1  840753  0.620827       C       T  0.406579     0.0   
    2   rs13303222    1  849998  0.620827       A       G  0.196053     0.0   
    3   rs72631889    1  851390  0.620827       T       G  0.034210     0.0   
    4  rs192998324    1  862772  0.620827       G       A  0.027632     0.0   
    
           BETA        SE  P_BOLT_LMM_INF  
    0  0.015560  0.258427           0.950  
    1 -0.060199  0.059829           0.310  
    2 -0.006768  0.078287           0.930  
    3  0.315642  0.172246           0.067  
    4 -0.227920  0.190562           0.230  
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/BOLT-LMM/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --extract SampleData1/Fold_0/train_data.valid.snp
      --out SampleData1/Fold_0/BOLT-LMM/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/SampleData1.lmmInfOnly_stat 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 172878 variants remaining.
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
    Results written to SampleData1/Fold_0/BOLT-LMM/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/BOLT-LMM/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --extract SampleData1/Fold_0/train_data.valid.snp
      --out SampleData1/Fold_0/BOLT-LMM/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/SampleData1.lmmInfOnly_stat 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    --extract: 172878 variants remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

     done.
    Total genotyping rate is 0.999891.
    172878 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/BOLT-LMM/test_data.*.profile.
    Continous Phenotype!
    File or directory does not exist: SampleData1/Fold_0/SampleData1.lmmForceNonInf_stat
           FID      IID  Sex
    0  HG00097  HG00097    2
    1  HG00099  HG00099    2
    2  HG00101  HG00101    1
    3  HG00102  HG00102    2
    4  HG00103  HG00103    1
    380
    380
             FID      IID  Sex       PC1       PC2       PC3       PC4       PC5  \
    0    HG00097  HG00097    2 -0.001453  0.084820  0.006792  0.013653  0.027149   
    1    HG00099  HG00099    2 -0.002017  0.089514 -0.022355  0.001888 -0.000037   
    2    HG00101  HG00101    1 -0.000380  0.096056 -0.018231 -0.016026  0.012093   
    3    HG00102  HG00102    2  0.000292  0.071832  0.018087 -0.045180  0.028123   
    4    HG00103  HG00103    1 -0.008372  0.065005 -0.009089 -0.026468 -0.009184   
    ..       ...      ...  ...       ...       ...       ...       ...       ...   
    375  NA20818  NA20818    2 -0.047156 -0.040644 -0.052693  0.021050 -0.013389   
    376  NA20826  NA20826    2 -0.042629 -0.059404 -0.066130  0.006495 -0.009525   
    377  NA20827  NA20827    1 -0.044060 -0.053125 -0.065463  0.015030 -0.004314   
    378  NA20828  NA20828    2 -0.047621 -0.050577 -0.043164  0.003004 -0.016823   
    379  NA20832  NA20832    2 -0.041535 -0.049826 -0.047877  0.005951 -0.003770   
    
              PC6  
    0    0.032581  
    1    0.009107  
    2    0.019296  
    3   -0.003620  
    4   -0.030565  
    ..        ...  
    375 -0.047403  
    376  0.010779  
    377  0.003873  
    378  0.015832  
    379 -0.023086  
    
    [380 rows x 9 columns]
    ./bolt --bfile=SampleData1/Fold_0/train_data.QC.clumped.pruned --phenoFile=SampleData1/Fold_0/train_data.PHENO --phenoCol=PHENO --lmmForceNonInf --LDscoresFile=tables/LDSCORE.1000G_EUR.tab.gz --covarFile=SampleData1/Fold_0/train_data.COV_PCA --qCovarCol=COV_{1:2} --statsFile=SampleData1/Fold_0/SampleData1.lmmForceNonInf_stat
                          +-----------------------------+
                          |                       ___   |
                          |   BOLT-LMM, v2.4.1   /_ /   |
                          |   November 16, 2022   /_/   |
                          |   Po-Ru Loh            //   |
                          |                        /    |
                          +-----------------------------+
    
    Copyright (C) 2014-2022 Harvard University.
    Distributed under the GNU GPLv3 open source license.
    
    Compiled with USE_SSE: fast aligned memory access
    Compiled with USE_MKL: Intel Math Kernel Library linear algebra
    Boost version: 1_58
    
    Command line options:
    
    ./bolt \
        --bfile=SampleData1/Fold_0/train_data.QC.clumped.pruned \
        --phenoFile=SampleData1/Fold_0/train_data.PHENO \
        --phenoCol=PHENO \
        --lmmForceNonInf \
        --LDscoresFile=tables/LDSCORE.1000G_EUR.tab.gz \
        --covarFile=SampleData1/Fold_0/train_data.COV_PCA \
        --qCovarCol=COV_{1:2} \
        --statsFile=SampleData1/Fold_0/SampleData1.lmmForceNonInf_stat 
    
    Setting number of threads to 1
    fam: SampleData1/Fold_0/train_data.QC.clumped.pruned.fam
    bim(s): SampleData1/Fold_0/train_data.QC.clumped.pruned.bim
    bed(s): SampleData1/Fold_0/train_data.QC.clumped.pruned.bed
    
    === Reading genotype data ===
    
    Total indivs in PLINK data: Nbed = 380
    Total indivs stored in memory: N = 380
    Reading bim file #1: SampleData1/Fold_0/train_data.QC.clumped.pruned.bim
        Read 172878 snps
    Total snps in PLINK data: Mbed = 172878
    
    Breakdown of SNP pre-filtering results:
      172878 SNPs to include in model (i.e., GRM)
      0 additional non-GRM SNPs loaded
      0 excluded SNPs
    Allocating 172878 x 380/4 bytes to store genotypes
    Reading genotypes and performing QC filtering on snps and indivs...
    Reading bed file #1: SampleData1/Fold_0/train_data.QC.clumped.pruned.bed
        Expecting 16423410 (+3) bytes for 380 indivs, 172878 snps
    

.. parsed-literal::

    WARNING: Genetic map appears to be in cM units; rescaling by 0.01
    

.. parsed-literal::

    Total indivs after QC: 380
    Total post-QC SNPs: M = 172878
      Variance component 1: 172878 post-QC SNPs (name: 'modelSnps')
    Time for SnpData setup = 1.15584 sec
    
    === Reading phenotype and covariate data ===
    
    Read data for 380 indivs (ignored 0 without genotypes) from:
      SampleData1/Fold_0/train_data.COV_PCA
    Read data for 380 indivs (ignored 0 without genotypes) from:
      SampleData1/Fold_0/train_data.PHENO
    Number of indivs with no missing phenotype(s) to use: 380
        Using quantitative covariate: COV_1
        Using quantitative covariate: COV_2
        Using quantitative covariate: CONST_ALL_ONES
    Number of individuals used in analysis: Nused = 380
    Singular values of covariate matrix:
        S[0] = 36.3836
        S[1] = 5.21943
        S[2] = 0.995369
    Total covariate vectors: C = 3
    Total independent covariate vectors: Cindep = 3
    
    === Initializing Bolt object: projecting and normalizing SNPs ===
    
    

.. parsed-literal::

    NOTE: Using all-1s vector (constant term) in addition to specified covariates
    

.. parsed-literal::

    Number of chroms with >= 1 good SNP: 22
    Average norm of projected SNPs:           375.206815
    Dimension of all-1s proj space (Nused-1): 379
    Time for covariate data setup + Bolt initialization = 0.391653 sec
    
    Phenotype 1:   N = 380   mean = 170.14   std = 0.945865
    
    === Computing linear regression (LINREG) stats ===
    
    Time for computing LINREG stats = 0.152288 sec
    
    === Estimating variance parameters ===
    
    Using CGtol of 0.005 for this step
    Using default number of random trials: 15 (for Nused = 380)
    
    Estimating MC scaling f_REML at log(delta) = 1.09384, h2 = 0.25...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.42  rNorms/orig: (0.02,0.02)  res2s: 882.719..149.989
      iter 2:  time=0.42  rNorms/orig: (0.0005,0.001)  res2s: 883.769..150.228
      Converged at iter 2: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 36.4%, memory/overhead = 63.6%
      MCscaling: logDelta = 1.09, h2 = 0.250, f = 0.00255218
    
    Estimating MC scaling f_REML at log(delta) = -0.00476781, h2 = 0.5...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.42  rNorms/orig: (0.04,0.05)  res2s: 206.958..66.4433
      iter 2:  time=0.41  rNorms/orig: (0.002,0.005)  res2s: 207.859..66.8373
      iter 3:  time=0.41  rNorms/orig: (9e-05,0.0003)  res2s: 207.871..66.8446
      Converged at iter 3: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 35.4%, memory/overhead = 64.6%
      MCscaling: logDelta = -0.00, h2 = 0.500, f = 0.000475493
    
    Estimating MC scaling f_REML at log(delta) = -0.256314, h2 = 0.562557...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.41  rNorms/orig: (0.04,0.05)  res2s: 142.182..50.8157
      iter 2:  time=0.38  rNorms/orig: (0.002,0.006)  res2s: 142.949..51.1908
      iter 3:  time=0.37  rNorms/orig: (0.0001,0.0005)  res2s: 142.962..51.1995
      Converged at iter 3: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 35.9%, memory/overhead = 64.1%
      MCscaling: logDelta = -0.26, h2 = 0.563, f = 5.80701e-05
    
    Estimating MC scaling f_REML at log(delta) = -0.291308, h2 = 0.571149...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.38  rNorms/orig: (0.04,0.05)  res2s: 134.763..48.8336
      iter 2:  time=0.37  rNorms/orig: (0.002,0.007)  res2s: 135.51..49.2044
      iter 3:  time=0.36  rNorms/orig: (0.0001,0.0005)  res2s: 135.523..49.2132
      Converged at iter 3: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 35.9%, memory/overhead = 64.1%
      MCscaling: logDelta = -0.29, h2 = 0.571, f = 3.23365e-06
    
    Secant iteration for h2 estimation converged in 2 steps
    Estimated (pseudo-)heritability: h2g = 0.571
    To more precisely estimate variance parameters and estimate s.e., use --reml
    Variance params: sigma^2_K = 0.406791, logDelta = -0.291308, f = 3.23365e-06
    
    Time for fitting variance components = 5.40124 sec
    
    === Computing mixed model assoc stats (inf. model) ===
    
    Selected 30 SNPs for computation of prospective stat
    Tried 30; threw out 0 with GRAMMAR chisq > 5
    Assigning SNPs to 22 chunks for leave-out analysis
    Each chunk is excluded when testing SNPs belonging to the chunk
      Batch-solving 52 systems of equations using conjugate gradient iteration
      iter 1:  time=0.68  rNorms/orig: (0.04,0.07)  res2s: 48.7902..70.329
      iter 2:  time=0.68  rNorms/orig: (0.002,0.007)  res2s: 49.1781..70.6337
      iter 3:  time=0.68  rNorms/orig: (9e-05,0.0005)  res2s: 49.1875..70.6346
      Converged at iter 3: rNorms/orig all < CGtol=0.0005
      Time breakdown: dgemm = 63.8%, memory/overhead = 36.2%
    
    AvgPro: 0.844   AvgRetro: 0.839   Calibration: 1.006 (0.003)   (30 SNPs)
    Ratio of medians: 1.013   Median of ratios: 1.002
    
    Time for computing infinitesimal model assoc stats = 2.23798 sec
    
    === Estimating chip LD Scores using 400 indivs ===
    
    Reducing sample size to 376 for memory alignment
    

.. parsed-literal::

    WARNING: Only 380 indivs available; using all
    

.. parsed-literal::

    
    Time for estimating chip LD Scores = 0.614924 sec
    
    === Reading LD Scores for calibration of Bayesian assoc stats ===
    
    Looking up LD Scores...
      Looking for column header 'SNP': column number = 1
      Looking for column header 'LDSCORE': column number = 5
    Found LD Scores for 171753/172878 SNPs
    
    Estimating inflation of LINREG chisq stats using MLMe as reference...
    Filtering to SNPs with chisq stats, LD Scores, and MAF > 0.01
    # of SNPs passing filters before outlier removal: 171753/172878
    Masking windows around outlier snps (chisq > 20.0)
    # of SNPs remaining after outlier window removal: 171753/171753
    Intercept of LD Score regression for ref stats:   1.003 (0.005)
    Estimated attenuation: 0.409 (0.534)
    Intercept of LD Score regression for cur stats: 1.004 (0.005)
    Calibration factor (ref/cur) to multiply by:      0.999 (0.000)
    LINREG intercept inflation = 1.00087
    
    === Estimating mixture parameters by cross-validation ===
    
    Setting maximum number of iterations to 250 for this step
    Max CV folds to compute = 5 (to have > 10000 samples)
    
    ====> Starting CV fold 1 <====
    
        Using quantitative covariate: COV_1
        Using quantitative covariate: COV_2
        Using quantitative covariate: CONST_ALL_ONES
    Number of individuals used in analysis: Nused = 304
    Singular values of covariate matrix:
        S[0] = 32.5614
        S[1] = 4.66512
        S[2] = 0.894787
    Total covariate vectors: C = 3
    Total independent covariate vectors: Cindep = 3
    
    === Initializing Bolt object: projecting and normalizing SNPs ===
    
    

.. parsed-literal::

    NOTE: Using all-1s vector (constant term) in addition to specified covariates
    

.. parsed-literal::

    Number of chroms with >= 1 good SNP: 22
    Average norm of projected SNPs:           299.575111
    Dimension of all-1s proj space (Nused-1): 303
      Beginning variational Bayes
      iter 1:  time=1.10 for 18 active reps
      iter 2:  time=0.80 for 18 active reps  approxLL diffs: (37.61,49.27)
      iter 3:  time=0.80 for 18 active reps  approxLL diffs: (1.23,1.36)
      iter 4:  time=0.83 for 18 active reps  approxLL diffs: (0.14,0.21)
      iter 5:  time=0.80 for 18 active reps  approxLL diffs: (0.01,0.01)
      iter 6:  time=0.19 for  1 active reps  approxLL diffs: (0.00,0.00)
      Converged at iter 6: approxLL diffs each have been < LLtol=0.01
      Time breakdown: dgemm = 26.0%, memory/overhead = 74.0%
    Computing predictions on left-out cross-validation fold
    Time for computing predictions = 0.452956 sec
    
    Average PVEs obtained by param pairs tested (high to low):
      f2=0.5, p=0.5: 0.006720
      f2=0.3, p=0.5: 0.006720
      f2=0.5, p=0.2: 0.006720
                ...
     f2=0.3, p=0.01: 0.006710
    
    Detailed CV fold results:
      Absolute prediction MSE baseline (covariates only): 0.586308
      Absolute prediction MSE using standard LMM:         0.582368
      Absolute prediction MSE, fold-best  f2=0.5, p=0.5:  0.582368
        Absolute pred MSE using   f2=0.5, p=0.5: 0.582368
        Absolute pred MSE using   f2=0.5, p=0.2: 0.582368
        Absolute pred MSE using   f2=0.5, p=0.1: 0.582368
        Absolute pred MSE using  f2=0.5, p=0.05: 0.582368
        Absolute pred MSE using  f2=0.5, p=0.02: 0.582369
        Absolute pred MSE using  f2=0.5, p=0.01: 0.582372
        Absolute pred MSE using   f2=0.3, p=0.5: 0.582368
        Absolute pred MSE using   f2=0.3, p=0.2: 0.582368
        Absolute pred MSE using   f2=0.3, p=0.1: 0.582368
        Absolute pred MSE using  f2=0.3, p=0.05: 0.582369
        Absolute pred MSE using  f2=0.3, p=0.02: 0.582370
        Absolute pred MSE using  f2=0.3, p=0.01: 0.582374
        Absolute pred MSE using   f2=0.1, p=0.5: 0.582368
        Absolute pred MSE using   f2=0.1, p=0.2: 0.582368
        Absolute pred MSE using   f2=0.1, p=0.1: 0.582368
        Absolute pred MSE using  f2=0.1, p=0.05: 0.582369
        Absolute pred MSE using  f2=0.1, p=0.02: 0.582371
        Absolute pred MSE using  f2=0.1, p=0.01: 0.582370
    
    ====> End CV fold 1: 18 remaining param pair(s) <====
    
    Estimated proportion of variance explained using inf model: 0.007
    Relative improvement in prediction MSE using non-inf model: 0.000
    
    ====> Starting CV fold 2 <====
    
        Using quantitative covariate: COV_1
        Using quantitative covariate: COV_2
        Using quantitative covariate: CONST_ALL_ONES
    Number of individuals used in analysis: Nused = 304
    Singular values of covariate matrix:
        S[0] = 32.3248
        S[1] = 4.70299
        S[2] = 0.886997
    Total covariate vectors: C = 3
    Total independent covariate vectors: Cindep = 3
    
    === Initializing Bolt object: projecting and normalizing SNPs ===
    
    

.. parsed-literal::

    NOTE: Using all-1s vector (constant term) in addition to specified covariates
    

.. parsed-literal::

    Number of chroms with >= 1 good SNP: 22
    Average norm of projected SNPs:           299.595029
    Dimension of all-1s proj space (Nused-1): 303
      Beginning variational Bayes
      iter 1:  time=1.10 for 18 active reps
      iter 2:  time=0.80 for 18 active reps  approxLL diffs: (34.03,41.62)
      iter 3:  time=0.80 for 18 active reps  approxLL diffs: (1.12,1.16)
      iter 4:  time=0.80 for 18 active reps  approxLL diffs: (0.13,0.17)
      iter 5:  time=0.80 for 18 active reps  approxLL diffs: (0.01,0.01)
      Converged at iter 5: approxLL diffs each have been < LLtol=0.01
      Time breakdown: dgemm = 26.3%, memory/overhead = 73.7%
    Computing predictions on left-out cross-validation fold
    Time for computing predictions = 0.45387 sec
    
    Average PVEs obtained by param pairs tested (high to low):
      f2=0.5, p=0.5: 0.003876
      f2=0.3, p=0.5: 0.003875
      f2=0.5, p=0.2: 0.003874
                ...
     f2=0.1, p=0.01: 0.003516
    
    Detailed CV fold results:
      Absolute prediction MSE baseline (covariates only): 0.884599
      Absolute prediction MSE using standard LMM:         0.883686
      Absolute prediction MSE, fold-best  f2=0.5, p=0.5:  0.883686
        Absolute pred MSE using   f2=0.5, p=0.5: 0.883686
        Absolute pred MSE using   f2=0.5, p=0.2: 0.883690
        Absolute pred MSE using   f2=0.5, p=0.1: 0.883699
        Absolute pred MSE using  f2=0.5, p=0.05: 0.883717
        Absolute pred MSE using  f2=0.5, p=0.02: 0.883771
        Absolute pred MSE using  f2=0.5, p=0.01: 0.883859
        Absolute pred MSE using   f2=0.3, p=0.5: 0.883687
        Absolute pred MSE using   f2=0.3, p=0.2: 0.883697
        Absolute pred MSE using   f2=0.3, p=0.1: 0.883715
        Absolute pred MSE using  f2=0.3, p=0.05: 0.883752
        Absolute pred MSE using  f2=0.3, p=0.02: 0.883861
        Absolute pred MSE using  f2=0.3, p=0.01: 0.884043
        Absolute pred MSE using   f2=0.1, p=0.5: 0.883691
        Absolute pred MSE using   f2=0.1, p=0.2: 0.883709
        Absolute pred MSE using   f2=0.1, p=0.1: 0.883739
        Absolute pred MSE using  f2=0.1, p=0.05: 0.883800
        Absolute pred MSE using  f2=0.1, p=0.02: 0.883989
        Absolute pred MSE using  f2=0.1, p=0.01: 0.884320
    
    ====> End CV fold 2: 18 remaining param pair(s) <====
    
    ====> Starting CV fold 3 <====
    
        Using quantitative covariate: COV_1
        Using quantitative covariate: COV_2
        Using quantitative covariate: CONST_ALL_ONES
    Number of individuals used in analysis: Nused = 304
    Singular values of covariate matrix:
        S[0] = 32.9844
        S[1] = 4.5875
        S[2] = 0.882915
    Total covariate vectors: C = 3
    Total independent covariate vectors: Cindep = 3
    
    === Initializing Bolt object: projecting and normalizing SNPs ===
    
    

.. parsed-literal::

    NOTE: Using all-1s vector (constant term) in addition to specified covariates
    

.. parsed-literal::

    Number of chroms with >= 1 good SNP: 22
    Average norm of projected SNPs:           299.590344
    Dimension of all-1s proj space (Nused-1): 303
      Beginning variational Bayes
      iter 1:  time=1.10 for 18 active reps
      iter 2:  time=0.80 for 18 active reps  approxLL diffs: (34.43,42.51)
      iter 3:  time=0.80 for 18 active reps  approxLL diffs: (1.13,1.20)
      iter 4:  time=0.80 for 18 active reps  approxLL diffs: (0.13,0.17)
      iter 5:  time=0.80 for 18 active reps  approxLL diffs: (0.01,0.01)
      Converged at iter 5: approxLL diffs each have been < LLtol=0.01
      Time breakdown: dgemm = 26.4%, memory/overhead = 73.6%
    Computing predictions on left-out cross-validation fold
    Time for computing predictions = 0.454183 sec
    
    Average PVEs obtained by param pairs tested (high to low):
      f2=0.5, p=0.5: 0.007021
      f2=0.3, p=0.5: 0.007020
      f2=0.5, p=0.2: 0.007016
                ...
     f2=0.1, p=0.01: 0.006257
    
    Detailed CV fold results:
      Absolute prediction MSE baseline (covariates only): 0.83708
      Absolute prediction MSE using standard LMM:         0.825938
      Absolute prediction MSE, fold-best  f2=0.5, p=0.5:  0.825938
        Absolute pred MSE using   f2=0.5, p=0.5: 0.825938
        Absolute pred MSE using   f2=0.5, p=0.2: 0.825947
        Absolute pred MSE using   f2=0.5, p=0.1: 0.825969
        Absolute pred MSE using  f2=0.5, p=0.05: 0.826011
        Absolute pred MSE using  f2=0.5, p=0.02: 0.826137
        Absolute pred MSE using  f2=0.5, p=0.01: 0.826335
        Absolute pred MSE using   f2=0.3, p=0.5: 0.825940
        Absolute pred MSE using   f2=0.3, p=0.2: 0.825965
        Absolute pred MSE using   f2=0.3, p=0.1: 0.826007
        Absolute pred MSE using  f2=0.3, p=0.05: 0.826091
        Absolute pred MSE using  f2=0.3, p=0.02: 0.826336
        Absolute pred MSE using  f2=0.3, p=0.01: 0.826722
        Absolute pred MSE using   f2=0.1, p=0.5: 0.825949
        Absolute pred MSE using   f2=0.1, p=0.2: 0.825991
        Absolute pred MSE using   f2=0.1, p=0.1: 0.826062
        Absolute pred MSE using  f2=0.1, p=0.05: 0.826201
        Absolute pred MSE using  f2=0.1, p=0.02: 0.826608
        Absolute pred MSE using  f2=0.1, p=0.01: 0.827253
    
    ====> End CV fold 3: 18 remaining param pair(s) <====
    
    ====> Starting CV fold 4 <====
    
        Using quantitative covariate: COV_1
        Using quantitative covariate: COV_2
        Using quantitative covariate: CONST_ALL_ONES
    Number of individuals used in analysis: Nused = 304
    Singular values of covariate matrix:
        S[0] = 32.5614
        S[1] = 4.66474
        S[2] = 0.892081
    Total covariate vectors: C = 3
    Total independent covariate vectors: Cindep = 3
    
    === Initializing Bolt object: projecting and normalizing SNPs ===
    
    

.. parsed-literal::

    NOTE: Using all-1s vector (constant term) in addition to specified covariates
    

.. parsed-literal::

    Number of chroms with >= 1 good SNP: 22
    Average norm of projected SNPs:           299.591407
    Dimension of all-1s proj space (Nused-1): 303
      Beginning variational Bayes
      iter 1:  time=1.09 for 18 active reps
      iter 2:  time=0.80 for 18 active reps  approxLL diffs: (35.87,44.89)
      iter 3:  time=0.80 for 18 active reps  approxLL diffs: (1.18,1.29)
      iter 4:  time=0.80 for 18 active reps  approxLL diffs: (0.14,0.19)
      iter 5:  time=0.80 for 18 active reps  approxLL diffs: (0.01,0.01)
      iter 6:  time=0.19 for  1 active reps  approxLL diffs: (0.00,0.00)
      Converged at iter 6: approxLL diffs each have been < LLtol=0.01
      Time breakdown: dgemm = 25.8%, memory/overhead = 74.2%
    Computing predictions on left-out cross-validation fold
    Time for computing predictions = 0.454353 sec
    
    Average PVEs obtained by param pairs tested (high to low):
      f2=0.5, p=0.5: 0.006447
      f2=0.3, p=0.5: 0.006445
      f2=0.5, p=0.2: 0.006440
                ...
     f2=0.1, p=0.01: 0.005512
    
    Detailed CV fold results:
      Absolute prediction MSE baseline (covariates only): 0.721953
      Absolute prediction MSE using standard LMM:         0.718542
      Absolute prediction MSE, fold-best  f2=0.5, p=0.5:  0.718542
        Absolute pred MSE using   f2=0.5, p=0.5: 0.718542
        Absolute pred MSE using   f2=0.5, p=0.2: 0.718550
        Absolute pred MSE using   f2=0.5, p=0.1: 0.718565
        Absolute pred MSE using  f2=0.5, p=0.05: 0.718597
        Absolute pred MSE using  f2=0.5, p=0.02: 0.718689
        Absolute pred MSE using  f2=0.5, p=0.01: 0.718837
        Absolute pred MSE using   f2=0.3, p=0.5: 0.718544
        Absolute pred MSE using   f2=0.3, p=0.2: 0.718562
        Absolute pred MSE using   f2=0.3, p=0.1: 0.718594
        Absolute pred MSE using  f2=0.3, p=0.05: 0.718656
        Absolute pred MSE using  f2=0.3, p=0.02: 0.718840
        Absolute pred MSE using  f2=0.3, p=0.01: 0.719137
        Absolute pred MSE using   f2=0.1, p=0.5: 0.718551
        Absolute pred MSE using   f2=0.1, p=0.2: 0.718582
        Absolute pred MSE using   f2=0.1, p=0.1: 0.718634
        Absolute pred MSE using  f2=0.1, p=0.05: 0.718738
        Absolute pred MSE using  f2=0.1, p=0.02: 0.719050
        Absolute pred MSE using  f2=0.1, p=0.01: 0.719587
    
    ====> End CV fold 4: 18 remaining param pair(s) <====
    
    ====> Starting CV fold 5 <====
    
        Using quantitative covariate: COV_1
        Using quantitative covariate: COV_2
        Using quantitative covariate: CONST_ALL_ONES
    Number of individuals used in analysis: Nused = 304
    Singular values of covariate matrix:
        S[0] = 32.2773
        S[1] = 4.70983
        S[2] = 0.892752
    Total covariate vectors: C = 3
    Total independent covariate vectors: Cindep = 3
    
    === Initializing Bolt object: projecting and normalizing SNPs ===
    
    

.. parsed-literal::

    NOTE: Using all-1s vector (constant term) in addition to specified covariates
    

.. parsed-literal::

    Number of chroms with >= 1 good SNP: 22
    Average norm of projected SNPs:           299.583395
    Dimension of all-1s proj space (Nused-1): 303
      Beginning variational Bayes
      iter 1:  time=1.12 for 18 active reps
      iter 2:  time=0.80 for 18 active reps  approxLL diffs: (37.87,48.73)
      iter 3:  time=0.91 for 18 active reps  approxLL diffs: (1.23,1.38)
      iter 4:  time=0.90 for 18 active reps  approxLL diffs: (0.14,0.21)
      iter 5:  time=0.82 for 18 active reps  approxLL diffs: (0.01,0.01)
      iter 6:  time=0.19 for  1 active reps  approxLL diffs: (0.00,0.00)
      Converged at iter 6: approxLL diffs each have been < LLtol=0.01
      Time breakdown: dgemm = 25.6%, memory/overhead = 74.4%
    Computing predictions on left-out cross-validation fold
    Time for computing predictions = 0.453788 sec
    
    Average PVEs obtained by param pairs tested (high to low):
      f2=0.5, p=0.5: 0.007272
      f2=0.3, p=0.5: 0.007270
      f2=0.5, p=0.2: 0.007265
                ...
     f2=0.1, p=0.01: 0.006266
    
    Detailed CV fold results:
      Absolute prediction MSE baseline (covariates only): 0.550938
      Absolute prediction MSE using standard LMM:         0.545112
      Absolute prediction MSE, fold-best  f2=0.5, p=0.5:  0.545112
        Absolute pred MSE using   f2=0.5, p=0.5: 0.545112
        Absolute pred MSE using   f2=0.5, p=0.2: 0.545117
        Absolute pred MSE using   f2=0.5, p=0.1: 0.545127
        Absolute pred MSE using  f2=0.5, p=0.05: 0.545148
        Absolute pred MSE using  f2=0.5, p=0.02: 0.545211
        Absolute pred MSE using  f2=0.5, p=0.01: 0.545314
        Absolute pred MSE using   f2=0.3, p=0.5: 0.545114
        Absolute pred MSE using   f2=0.3, p=0.2: 0.545125
        Absolute pred MSE using   f2=0.3, p=0.1: 0.545146
        Absolute pred MSE using  f2=0.3, p=0.05: 0.545188
        Absolute pred MSE using  f2=0.3, p=0.02: 0.545312
        Absolute pred MSE using  f2=0.3, p=0.01: 0.545519
        Absolute pred MSE using   f2=0.1, p=0.5: 0.545118
        Absolute pred MSE using   f2=0.1, p=0.2: 0.545138
        Absolute pred MSE using   f2=0.1, p=0.1: 0.545173
        Absolute pred MSE using  f2=0.1, p=0.05: 0.545242
        Absolute pred MSE using  f2=0.1, p=0.02: 0.545452
        Absolute pred MSE using  f2=0.1, p=0.01: 0.545825
    
    ====> End CV fold 5: 18 remaining param pair(s) <====
    
    Optimal mixture parameters according to CV: f2 = 0.5, p = 0.5
    
    Time for estimating mixture parameters = 53.0219 sec
    
    === Computing Bayesian mixed model assoc stats with mixture prior ===
    
    Assigning SNPs to 22 chunks for leave-out analysis
    Each chunk is excluded when testing SNPs belonging to the chunk
      Beginning variational Bayes
      iter 1:  time=1.20 for 22 active reps
      iter 2:  time=0.91 for 22 active reps  approxLL diffs: (44.97,45.17)
      iter 3:  time=0.98 for 22 active reps  approxLL diffs: (1.48,1.49)
      iter 4:  time=1.00 for 22 active reps  approxLL diffs: (0.17,0.17)
      iter 5:  time=0.90 for 22 active reps  approxLL diffs: (0.01,0.01)
      Converged at iter 5: approxLL diffs each have been < LLtol=0.01
      Time breakdown: dgemm = 26.3%, memory/overhead = 73.7%
    Filtering to SNPs with chisq stats, LD Scores, and MAF > 0.01
    # of SNPs passing filters before outlier removal: 171753/172878
    Masking windows around outlier snps (chisq > 20.0)
    # of SNPs remaining after outlier window removal: 171753/171753
    Intercept of LD Score regression for ref stats:   1.003 (0.005)
    Estimated attenuation: 0.409 (0.534)
    Intercept of LD Score regression for cur stats: 0.997 (0.005)
    Calibration factor (ref/cur) to multiply by:      1.006 (0.000)
    
    Time for computing Bayesian mixed model assoc stats = 5.7944 sec
    
    Calibration stats: mean and lambdaGC (over SNPs used in GRM)
      (note that both should be >1 because of polygenicity)
    Mean BOLT_LMM_INF: 1.0074 (172878 good SNPs)   lambdaGC: 1.01336
    Mean BOLT_LMM: 1.00721 (172878 good SNPs)   lambdaGC: 1.01326
    
    === Streaming genotypes to compute and write assoc stats at all SNPs ===
    
    Time for streaming genotypes and writing output = 1.77796 sec
    
    Total elapsed time for analysis = 70.5482 sec
               SNP  CHR      BP    GENPOS ALLELE1 ALLELE0    A1FREQ  F_MISS  \
    0   rs79373928    1  801536  0.587220       G       T  0.014474     0.0   
    1    rs4970382    1  840753  0.620827       C       T  0.406579     0.0   
    2   rs13303222    1  849998  0.620827       A       G  0.196053     0.0   
    3   rs72631889    1  851390  0.620827       T       G  0.034210     0.0   
    4  rs192998324    1  862772  0.620827       G       A  0.027632     0.0   
    
           BETA        SE  P_BOLT_LMM_INF  P_BOLT_LMM  
    0  0.015560  0.258427           0.950       0.950  
    1 -0.060199  0.059829           0.310       0.310  
    2 -0.006768  0.078287           0.930       0.930  
    3  0.315642  0.172246           0.067       0.067  
    4 -0.227920  0.190562           0.230       0.230  
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/BOLT-LMM/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --extract SampleData1/Fold_0/train_data.valid.snp
      --out SampleData1/Fold_0/BOLT-LMM/train_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/SampleData1.lmmForceNonInf_stat 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 172878 variants remaining.
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
    Results written to SampleData1/Fold_0/BOLT-LMM/train_data.*.profile.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/BOLT-LMM/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --extract SampleData1/Fold_0/train_data.valid.snp
      --out SampleData1/Fold_0/BOLT-LMM/test_data
      --q-score-range SampleData1/Fold_0/range_list SampleData1/Fold_0/SNP.pvalue
      --score SampleData1/Fold_0/SampleData1.lmmForceNonInf_stat 1 2 3 header
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    --extract: 172878 variants remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    172878 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172878 valid predictors loaded.
    

.. parsed-literal::

    Warning: 326740 lines skipped in --q-score-range data file.
    

.. parsed-literal::

    --score: 20 ranges processed.
    Results written to SampleData1/Fold_0/BOLT-LMM/test_data.*.profile.
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
   python BOLT-LMM.py 0
   python BOLT-LMM.py 1
   python BOLT-LMM.py 2
   python BOLT-LMM.py 3
   python BOLT-LMM.py 4

The following files should exist after the execution:

1. ``SampleData1/Fold_0/BOLT-LMM/Results.csv``
2. ``SampleData1/Fold_1/BOLT-LMM/Results.csv``
3. ``SampleData1/Fold_2/BOLT-LMM/Results.csv``
4. ``SampleData1/Fold_3/BOLT-LMM/Results.csv``
5. ``SampleData1/Fold_4/BOLT-LMM/Results.csv``

Check the results file for each fold.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

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
            #print(temp.tail())
            print("Number of P-values processed: ",len(temp))
            # Print a message indicating that the file exists
        
        else:
            # Print a message indicating that the file does not exist
            print("Fold_",loop, "No, the file does not exist.")
    
    


.. parsed-literal::

    Fold_ 0 Yes, the file exists.
    Number of P-values processed:  60
    Fold_ 1 No, the file does not exist.
    Fold_ 2 Yes, the file exists.
    Number of P-values processed:  60
    Fold_ 3 Yes, the file exists.
    Number of P-values processed:  60
    Fold_ 4 Yes, the file exists.
    Number of P-values processed:  60
    

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
            'BOLTmodel',
         
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
    
     


.. parsed-literal::

    We have to ensure when we sum the entries across all Folds, the same rows are merged!
    Fold_ 0 Yes, the file exists.
    Fold_ 1 No, the file does not exist.
    Fold_ 2 No, the file does not exist.
    Fold_ 3 No, the file does not exist.
    Fold_ 4 No, the file does not exist.
    DataFrame 1 with extracted common rows has 60 rows.
        clump_p1  clump_r2  clump_kb  p_window_size  p_slide_size  p_LD_threshold  \
    0        1.0       0.1     200.0          200.0          50.0            0.25   
    1        1.0       0.1     200.0          200.0          50.0            0.25   
    2        1.0       0.1     200.0          200.0          50.0            0.25   
    3        1.0       0.1     200.0          200.0          50.0            0.25   
    4        1.0       0.1     200.0          200.0          50.0            0.25   
    5        1.0       0.1     200.0          200.0          50.0            0.25   
    6        1.0       0.1     200.0          200.0          50.0            0.25   
    7        1.0       0.1     200.0          200.0          50.0            0.25   
    8        1.0       0.1     200.0          200.0          50.0            0.25   
    9        1.0       0.1     200.0          200.0          50.0            0.25   
    10       1.0       0.1     200.0          200.0          50.0            0.25   
    11       1.0       0.1     200.0          200.0          50.0            0.25   
    12       1.0       0.1     200.0          200.0          50.0            0.25   
    13       1.0       0.1     200.0          200.0          50.0            0.25   
    14       1.0       0.1     200.0          200.0          50.0            0.25   
    15       1.0       0.1     200.0          200.0          50.0            0.25   
    16       1.0       0.1     200.0          200.0          50.0            0.25   
    17       1.0       0.1     200.0          200.0          50.0            0.25   
    18       1.0       0.1     200.0          200.0          50.0            0.25   
    19       1.0       0.1     200.0          200.0          50.0            0.25   
    20       1.0       0.1     200.0          200.0          50.0            0.25   
    21       1.0       0.1     200.0          200.0          50.0            0.25   
    22       1.0       0.1     200.0          200.0          50.0            0.25   
    23       1.0       0.1     200.0          200.0          50.0            0.25   
    24       1.0       0.1     200.0          200.0          50.0            0.25   
    25       1.0       0.1     200.0          200.0          50.0            0.25   
    26       1.0       0.1     200.0          200.0          50.0            0.25   
    27       1.0       0.1     200.0          200.0          50.0            0.25   
    28       1.0       0.1     200.0          200.0          50.0            0.25   
    29       1.0       0.1     200.0          200.0          50.0            0.25   
    30       1.0       0.1     200.0          200.0          50.0            0.25   
    31       1.0       0.1     200.0          200.0          50.0            0.25   
    32       1.0       0.1     200.0          200.0          50.0            0.25   
    33       1.0       0.1     200.0          200.0          50.0            0.25   
    34       1.0       0.1     200.0          200.0          50.0            0.25   
    35       1.0       0.1     200.0          200.0          50.0            0.25   
    36       1.0       0.1     200.0          200.0          50.0            0.25   
    37       1.0       0.1     200.0          200.0          50.0            0.25   
    38       1.0       0.1     200.0          200.0          50.0            0.25   
    39       1.0       0.1     200.0          200.0          50.0            0.25   
    40       1.0       0.1     200.0          200.0          50.0            0.25   
    41       1.0       0.1     200.0          200.0          50.0            0.25   
    42       1.0       0.1     200.0          200.0          50.0            0.25   
    43       1.0       0.1     200.0          200.0          50.0            0.25   
    44       1.0       0.1     200.0          200.0          50.0            0.25   
    45       1.0       0.1     200.0          200.0          50.0            0.25   
    46       1.0       0.1     200.0          200.0          50.0            0.25   
    47       1.0       0.1     200.0          200.0          50.0            0.25   
    48       1.0       0.1     200.0          200.0          50.0            0.25   
    49       1.0       0.1     200.0          200.0          50.0            0.25   
    50       1.0       0.1     200.0          200.0          50.0            0.25   
    51       1.0       0.1     200.0          200.0          50.0            0.25   
    52       1.0       0.1     200.0          200.0          50.0            0.25   
    53       1.0       0.1     200.0          200.0          50.0            0.25   
    54       1.0       0.1     200.0          200.0          50.0            0.25   
    55       1.0       0.1     200.0          200.0          50.0            0.25   
    56       1.0       0.1     200.0          200.0          50.0            0.25   
    57       1.0       0.1     200.0          200.0          50.0            0.25   
    58       1.0       0.1     200.0          200.0          50.0            0.25   
    59       1.0       0.1     200.0          200.0          50.0            0.25   
    
              pvalue  numberofpca  tempalpha  l1weight  numberofvariants  \
    0   1.000000e-10          6.0        0.1       0.1               0.0   
    1   3.359818e-10          6.0        0.1       0.1               0.0   
    2   1.128838e-09          6.0        0.1       0.1               0.0   
    3   3.792690e-09          6.0        0.1       0.1               0.0   
    4   1.274275e-08          6.0        0.1       0.1               0.0   
    5   4.281332e-08          6.0        0.1       0.1               0.0   
    6   1.438450e-07          6.0        0.1       0.1               0.0   
    7   4.832930e-07          6.0        0.1       0.1               0.0   
    8   1.623777e-06          6.0        0.1       0.1               0.0   
    9   5.455595e-06          6.0        0.1       0.1               0.0   
    10  1.832981e-05          6.0        0.1       0.1               0.0   
    11  6.158482e-05          6.0        0.1       0.1               0.0   
    12  2.069138e-04          6.0        0.1       0.1               0.0   
    13  6.951928e-04          6.0        0.1       0.1               0.0   
    14  2.335721e-03          6.0        0.1       0.1               0.0   
    15  7.847600e-03          6.0        0.1       0.1               0.0   
    16  2.636651e-02          6.0        0.1       0.1               0.0   
    17  8.858668e-02          6.0        0.1       0.1               0.0   
    18  2.976351e-01          6.0        0.1       0.1               0.0   
    19  1.000000e+00          6.0        0.1       0.1               0.0   
    20  1.000000e-10          6.0        0.1       0.1               0.0   
    21  3.359818e-10          6.0        0.1       0.1               0.0   
    22  1.128838e-09          6.0        0.1       0.1               0.0   
    23  3.792690e-09          6.0        0.1       0.1               0.0   
    24  1.274275e-08          6.0        0.1       0.1               0.0   
    25  4.281332e-08          6.0        0.1       0.1               0.0   
    26  1.438450e-07          6.0        0.1       0.1               0.0   
    27  4.832930e-07          6.0        0.1       0.1               0.0   
    28  1.623777e-06          6.0        0.1       0.1               0.0   
    29  5.455595e-06          6.0        0.1       0.1               0.0   
    30  1.832981e-05          6.0        0.1       0.1               0.0   
    31  6.158482e-05          6.0        0.1       0.1               0.0   
    32  2.069138e-04          6.0        0.1       0.1               0.0   
    33  6.951928e-04          6.0        0.1       0.1               0.0   
    34  2.335721e-03          6.0        0.1       0.1               0.0   
    35  7.847600e-03          6.0        0.1       0.1               0.0   
    36  2.636651e-02          6.0        0.1       0.1               0.0   
    37  8.858668e-02          6.0        0.1       0.1               0.0   
    38  2.976351e-01          6.0        0.1       0.1               0.0   
    39  1.000000e+00          6.0        0.1       0.1               0.0   
    40  1.000000e-10          6.0        0.1       0.1               0.0   
    41  3.359818e-10          6.0        0.1       0.1               0.0   
    42  1.128838e-09          6.0        0.1       0.1               0.0   
    43  3.792690e-09          6.0        0.1       0.1               0.0   
    44  1.274275e-08          6.0        0.1       0.1               0.0   
    45  4.281332e-08          6.0        0.1       0.1               0.0   
    46  1.438450e-07          6.0        0.1       0.1               0.0   
    47  4.832930e-07          6.0        0.1       0.1               0.0   
    48  1.623777e-06          6.0        0.1       0.1               0.0   
    49  5.455595e-06          6.0        0.1       0.1               0.0   
    50  1.832981e-05          6.0        0.1       0.1               0.0   
    51  6.158482e-05          6.0        0.1       0.1               0.0   
    52  2.069138e-04          6.0        0.1       0.1               0.0   
    53  6.951928e-04          6.0        0.1       0.1               0.0   
    54  2.335721e-03          6.0        0.1       0.1               0.0   
    55  7.847600e-03          6.0        0.1       0.1               0.0   
    56  2.636651e-02          6.0        0.1       0.1               0.0   
    57  8.858668e-02          6.0        0.1       0.1               0.0   
    58  2.976351e-01          6.0        0.1       0.1               0.0   
    59  1.000000e+00          6.0        0.1       0.1               0.0   
    
        Train_pure_prs  Train_null_model  Train_best_model  Test_pure_prs  \
    0         0.002293          0.227477          0.602341   2.369136e-04   
    1         0.002266          0.227477          0.638998   2.858599e-04   
    2         0.002376          0.227477          0.697855   3.838787e-04   
    3         0.002409          0.227477          0.734995   2.228255e-04   
    4         0.002395          0.227477          0.776750   6.394029e-05   
    5         0.002321          0.227477          0.812544   6.596809e-05   
    6         0.002306          0.227477          0.835596   1.222748e-04   
    7         0.002205          0.227477          0.864166   1.204769e-04   
    8         0.002158          0.227477          0.883590   9.661765e-05   
    9         0.002167          0.227477          0.909570   5.687403e-05   
    10        0.002181          0.227477          0.933095   4.772806e-05   
    11        0.002209          0.227477          0.947256   7.291404e-05   
    12        0.002186          0.227477          0.961601   8.539601e-07   
    13        0.002202          0.227477          0.975360  -3.215678e-05   
    14        0.002203          0.227477          0.982642  -1.351810e-05   
    15        0.002183          0.227477          0.989179  -6.292480e-07   
    16        0.002147          0.227477          0.993029   1.580836e-05   
    17        0.002146          0.227477          0.995876   2.640523e-05   
    18        0.002128          0.227477          0.997875   1.470785e-05   
    19        0.002117          0.227477          0.999334   1.192114e-05   
    20        0.002293          0.227477          0.602341   2.369136e-04   
    21        0.002266          0.227477          0.638998   2.858599e-04   
    22        0.002376          0.227477          0.697855   3.838787e-04   
    23        0.002409          0.227477          0.734995   2.228255e-04   
    24        0.002395          0.227477          0.776750   6.394029e-05   
    25        0.002321          0.227477          0.812544   6.596809e-05   
    26        0.002306          0.227477          0.835596   1.222748e-04   
    27        0.002205          0.227477          0.864166   1.204769e-04   
    28        0.002158          0.227477          0.883590   9.661765e-05   
    29        0.002167          0.227477          0.909570   5.687403e-05   
    30        0.002181          0.227477          0.933095   4.772806e-05   
    31        0.002209          0.227477          0.947256   7.291404e-05   
    32        0.002186          0.227477          0.961601   8.539601e-07   
    33        0.002202          0.227477          0.975360  -3.215678e-05   
    34        0.002203          0.227477          0.982642  -1.351810e-05   
    35        0.002183          0.227477          0.989179  -6.292480e-07   
    36        0.002147          0.227477          0.993029   1.580836e-05   
    37        0.002146          0.227477          0.995876   2.640523e-05   
    38        0.002128          0.227477          0.997875   1.470785e-05   
    39        0.002117          0.227477          0.999334   1.192114e-05   
    40        0.002293          0.227477          0.602341   2.369136e-04   
    41        0.002266          0.227477          0.638998   2.858599e-04   
    42        0.002376          0.227477          0.697855   3.838787e-04   
    43        0.002409          0.227477          0.734995   2.228255e-04   
    44        0.002395          0.227477          0.776750   6.394029e-05   
    45        0.002321          0.227477          0.812544   6.596809e-05   
    46        0.002306          0.227477          0.835596   1.222748e-04   
    47        0.002205          0.227477          0.864166   1.204769e-04   
    48        0.002158          0.227477          0.883590   9.661765e-05   
    49        0.002167          0.227477          0.909570   5.687403e-05   
    50        0.002181          0.227477          0.933095   4.772806e-05   
    51        0.002209          0.227477          0.947256   7.291404e-05   
    52        0.002186          0.227477          0.961601   8.539601e-07   
    53        0.002202          0.227477          0.975360  -3.215678e-05   
    54        0.002203          0.227477          0.982642  -1.351810e-05   
    55        0.002183          0.227477          0.989179  -6.292480e-07   
    56        0.002147          0.227477          0.993029   1.580836e-05   
    57        0.002146          0.227477          0.995876   2.640523e-05   
    58        0.002128          0.227477          0.997875   1.470785e-05   
    59        0.002117          0.227477          0.999334   1.192114e-05   
    
        Test_null_model  Test_best_model       BOLTmodel  
    0           0.14297         0.181276             lmm  
    1           0.14297         0.199960             lmm  
    2           0.14297         0.209241             lmm  
    3           0.14297         0.158463             lmm  
    4           0.14297         0.103519             lmm  
    5           0.14297         0.126593             lmm  
    6           0.14297         0.142460             lmm  
    7           0.14297         0.166682             lmm  
    8           0.14297         0.217751             lmm  
    9           0.14297         0.185704             lmm  
    10          0.14297         0.182166             lmm  
    11          0.14297         0.184267             lmm  
    12          0.14297         0.162743             lmm  
    13          0.14297         0.149901             lmm  
    14          0.14297         0.171236             lmm  
    15          0.14297         0.191330             lmm  
    16          0.14297         0.199744             lmm  
    17          0.14297         0.216609             lmm  
    18          0.14297         0.207610             lmm  
    19          0.14297         0.210739             lmm  
    20          0.14297         0.181276      lmmInfOnly  
    21          0.14297         0.199960      lmmInfOnly  
    22          0.14297         0.209241      lmmInfOnly  
    23          0.14297         0.158463      lmmInfOnly  
    24          0.14297         0.103519      lmmInfOnly  
    25          0.14297         0.126593      lmmInfOnly  
    26          0.14297         0.142460      lmmInfOnly  
    27          0.14297         0.166682      lmmInfOnly  
    28          0.14297         0.217751      lmmInfOnly  
    29          0.14297         0.185704      lmmInfOnly  
    30          0.14297         0.182166      lmmInfOnly  
    31          0.14297         0.184267      lmmInfOnly  
    32          0.14297         0.162743      lmmInfOnly  
    33          0.14297         0.149901      lmmInfOnly  
    34          0.14297         0.171236      lmmInfOnly  
    35          0.14297         0.191330      lmmInfOnly  
    36          0.14297         0.199744      lmmInfOnly  
    37          0.14297         0.216609      lmmInfOnly  
    38          0.14297         0.207610      lmmInfOnly  
    39          0.14297         0.210739      lmmInfOnly  
    40          0.14297         0.181276  lmmForceNonInf  
    41          0.14297         0.199960  lmmForceNonInf  
    42          0.14297         0.209241  lmmForceNonInf  
    43          0.14297         0.158463  lmmForceNonInf  
    44          0.14297         0.103519  lmmForceNonInf  
    45          0.14297         0.126593  lmmForceNonInf  
    46          0.14297         0.142460  lmmForceNonInf  
    47          0.14297         0.166682  lmmForceNonInf  
    48          0.14297         0.217751  lmmForceNonInf  
    49          0.14297         0.185704  lmmForceNonInf  
    50          0.14297         0.182166  lmmForceNonInf  
    51          0.14297         0.184267  lmmForceNonInf  
    52          0.14297         0.162743  lmmForceNonInf  
    53          0.14297         0.149901  lmmForceNonInf  
    54          0.14297         0.171236  lmmForceNonInf  
    55          0.14297         0.191330  lmmForceNonInf  
    56          0.14297         0.199744  lmmForceNonInf  
    57          0.14297         0.216609  lmmForceNonInf  
    58          0.14297         0.207610  lmmForceNonInf  
    59          0.14297         0.210739  lmmForceNonInf  
    

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
    
    |                  | 59                    |
    |:-----------------|:----------------------|
    | clump_p1         | 1.0                   |
    | clump_r2         | 0.1                   |
    | clump_kb         | 200.0                 |
    | p_window_size    | 200.0                 |
    | p_slide_size     | 50.0                  |
    | p_LD_threshold   | 0.25                  |
    | pvalue           | 1.0                   |
    | numberofpca      | 6.0                   |
    | tempalpha        | 0.1                   |
    | l1weight         | 0.1                   |
    | numberofvariants | 0.0                   |
    | Train_pure_prs   | 0.0021171291324699    |
    | Train_null_model | 0.2274767654484456    |
    | Train_best_model | 0.9993342811094132    |
    | Test_pure_prs    | 1.192113731751654e-05 |
    | Test_null_model  | 0.1429700622403754    |
    | Test_best_model  | 0.2107386792427508    |
    | BOLTmodel        | lmmForceNonInf        |
    


.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAu4AAAHCCAYAAACwm0waAAAAAXNSR0IArs4c6QAAIABJREFUeF7snQd4VEX3xt/0XugISFGqVFEUBKR3pIkgCNJFwc8PC4iCn6igItgRQXq30HsHEZQiHaQjht7Te/n/z+CGJWyS3ezevVveeZ5oSOZO+Z1J8s65Z854ZGRkZICFBEiABEiABEiABEiABEjAoQl4ULg7tH04OBIgARIgARIgARIgARJQBCjcuRBIgARIgARIgARIgARIwAkIULg7gZE4RBIgARIgARIgARIgARKgcOcaIAESIAESIAESIAESIAEnIEDh7gRG4hBJgARIgARIgARIgARIgMKda4AESIAESIAESIAESIAEnIAAhbsTGIlDJAESIAESIAESIAESIAEKd64BEiABEiABEiABEiABEnACAhTuTmAkDpEESIAESIAESIAESIAEKNy5BkiABEiABEiABEiABEjACQhQuDuBkThEEiABEiABEiABEiABEqBw5xogARIgARIgARIgARIgAScgQOHuBEbiEEmABEiABEiABEiABEiAwp1rgARIgARIgARIgARIgAScgACFuxMYiUMkARIgARIgARIgARIgAQp3rgESIAESIAESIAESIAEScAICFO5OYCQOkQRIgARIgARIgARIgAQo3LkGSIAESIAESIAESIAESMAJCFC4O4GROEQSIAESIAESIAESIAESoHDnGiABEiABEiABEiABEiABJyBA4e4ERuIQSYAESIAESIAESIAESIDCnWuABEiABEiABEiABEiABJyAAIW7ExiJQyQBEiABEiABEiABEiABCneuARIgARIgARIgARIgARJwAgIU7k5gJA6RBEiABEiABEiABEiABCjcuQZIgARIgARIgARIgARIwAkIULg7gZE4RBIgARIgARIgARIgARKgcOcaIAESIAESIAESIAESIAEnIEDh7gRG4hBJgARIgARIgARIgARIgMKda4AESIAESIAESIAESIAEnIAAhbsTGIlDJAESIAESIAESIAESIAEKd64BEiABEiABEiABEiABEnACAhTuTmAkDpEESIAESIAESIAESIAEKNy5BkiABEiABEiABEiABEjACQhQuDuBkThEEiABEiABEiABEiABEqBw5xogARIgARIgARIgARIgAScgQOHuBEbiEEmA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j5IIpBoBnqjDeHj/hU/+/kkT2Vs6pYU7jz/iH0Cev0qxziuKY01cfuzRo4f5I8rlx6FDh4YNs7n11lvNEhCvW46WNEijEdL7IiACIpC+BOidVBKBVCMQbhVJmsje0ikt3BkjyAsFmAIHi/9UmePHj0ckxVisKlWqmPh4DqrWrVuD4pzXDwcnxtHzhjEejSThbj8AlVMEREAEREAERCC9CUi429s/pYW7PYbQORctWmRi4XlRAxMvGOCVvxTygWnPnj1o0qQJRo4cidtuu03CPafg9bwIiIAIiIAIiEDaEJBwtze1hHsEVhMnTsSkSZPMqTRMa9euRbt27cyVyYGJVw936tTJbLbgrWbyuNsPQOUUAREQAREQARFIbwIS7vb2l3CPwIqn0EyePDmicJ81axbGjBkDXnCwZcuWsMKdnnp++ROvdN6xY4e9pZRTBERABERABERABFKQgIS7vVEl3COwChUq8/rrr2P69OlZTz3++OPmqt9zzjkHJ0+exM6dO1GtWjWsXLkyohU0SO0HqXKKgAiIgAiIgAikLgFpInvbSrhHYHXq1ClUrlwZvNLXvzn1lltuQceOHUM+FcnjHvyABqn9IFVOERABERABERCB1CUgTWRvWwn3KKw+/fRT9OzZ0xwHWa9ePQwbNgwzZ840G1SHDx9+xtMS7vYDTzlFQAREQAREQAREgAQk3O3HgYS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNtNp0/5AAAgAElEQVSzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStHc2qQOopThYmACIiACIiACCQpAWkie8NJuNuzcjSnBqmjOFWYCIiACIiACIhAkhKQJrI3nIS7PStH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EAEREAERiD8BCff4M1eNIiACIiACIiACIiACIhAzAQn3mJHpAREQAREQAREQAREQARGIPwEJ9/gzV40iIAIikJAE/vOf/2D+/PlZR8gmZCPVKBEQARFIYwIS7mlsfHVdBERABAIJSLhrPIiACIhAYhOQcE9s+6h1IiACInAWAV588uijj2L69Onm4pq33377jIu7+MCIESPMTaFDhw41zy9ZssRcRvTll1/i+eefx+TJk3H8+HFzqyifL1asGAKFO2/vZT28BIypQYMG5ucbbrgB33//PXr06IGffvrJvPfyyy+bi9+UREAEREAE3CUg4e4uX5UuAiIgAo4TyJUrFwYPHozu3btj6tSp4O2TX3311Rn1HDhwANWrVwdvd6W453XylSpVMoJ77969RqgzvfLKK9i9e7e5SdlWuN94440YOHCgKX/r1q1GzG/atAlsl5IIiIAIiIB7BCTc3WOrkkVABETAFQIUyD///DMKFixoyi9atCh++OGHrN/9lbZo0QIdOnTAbbfdhosvvtiI+xIlShhPPYX64cOHcezYMfzhD3/AjBkzrIT7FVdcYUR/tWrVsvq2Z88eLF26FBdddJEr/VWhIiACIiACPgIS7hoJIiACIpBkBMIJ91atWmH79u34/e9/j2nTpmHixIkYP3487r33XhPOMnPmTGRkZKB06dImjOaSSy4xIn7AgAGYM2fOGcL92WefRWZmJp566ilDp27duuBrl19+OcqXL4/9+/cnGTU1VwREQASSn4CEe/LbUD0QARFIMwIU7m+++Sa6detmPOWMRw8OlSESinR62im6mzZtijvvvBMHDx5EhQoVTIhLgQIFQK88vffBwn306NFG9E+ZMsWEwdSqVQsffvihCYtheZ06dTJfTMuWLTOCXkkEREAERMBdAhLu7vJV6SIgAiLgOAFuGn3ssceMt5w/MzadwjpUuueeezBmzBjs2rUrK5SGMfEjR45E8eLFUb9+fSO8g4X7L7/8gubNm5sYeZbNUJznnnsua3Nqz549jfg/ceIE6tSpY8pTEgEREAERcJeAhLu7fFW6CIiACDhOgGL95MmTjperAkVABERABBKbgIR7YttHrRMBERCBswhIuGtQiIAIiEB6EpBwT0+7q9ciIAIiIAIiIAIiIAJJRkDCPckMpuaKgAiIgAiIgAiIgAikJwEJ9/S0u3otAiIgAiIgAiIgAiKQZAQk3JPMYGquCIiACIiACIiACIhAehKQcE9Pu6vXIiACIiACIiACIiACSUZAwj3JDKbmioAIiIAIiIAIiIAIpCcBCff0tLt6LQIiIAIiIAIiIAIikGQEJNyTzGBqrgiIgAiIgAiIgAiIQHoSkHBPT7ur1yIgAiIgAiIgAiIgAklGQMI9yQym5oqACIiACIiACIiACKQnAQn39LS7ei0CIiACIiACIiACIpBkBCTck8xgaq4IiIAIiIAIiIAIiEB6EpBwT0+7q9ciIAIiIAIiIAIiIAJJRkDCPckMpuaKgAiIgAiIgAiIgAikJwEJ9/S0u3otAiIgAiIgAiIgAiKQZAQk3JPMYGquCIiACIiACIiACIhAehKQcE9Pu6vXIiACIiACIiACIiACSUZAwj3JDKbmioAIiIAIiIAIiIAIpCeB/wcClL22hkcXMAAAAABJRU5ErkJggg==" width="1000">


.. parsed-literal::

    2. Reporting Generalized Performance:
    
    |                  | 17                    |
    |:-----------------|:----------------------|
    | clump_p1         | 1.0                   |
    | clump_r2         | 0.1                   |
    | clump_kb         | 200.0                 |
    | p_window_size    | 200.0                 |
    | p_slide_size     | 50.0                  |
    | p_LD_threshold   | 0.25                  |
    | pvalue           | 0.0885866790410083    |
    | numberofpca      | 6.0                   |
    | tempalpha        | 0.1                   |
    | l1weight         | 0.1                   |
    | numberofvariants | 0.0                   |
    | Train_pure_prs   | 0.0021456782406613    |
    | Train_null_model | 0.2274767654484456    |
    | Train_best_model | 0.9958762340331948    |
    | Test_pure_prs    | 2.640522907393361e-05 |
    | Test_null_model  | 0.1429700622403754    |
    | Test_best_model  | 0.2166086418607637    |
    | BOLTmodel        | lmm                   |
    | Difference       | 0.7792675921724311    |
    | Sum              | 1.2124848758939586    |
    3. Reporting the correlation of hyperparameters and the performance of 'Train_null_model', 'Train_pure_prs', 'Train_best_model', 'Test_pure_prs', 'Test_null_model', and 'Test_best_model':
    
    3. For string hyperparameters, we used one-hot encoding to find the correlation between string hyperparameters and 'Train_null_model', 'Train_pure_prs', 'Train_best_model', 'Test_pure_prs', 'Test_null_model', and 'Test_best_model'.
    3. We performed this analysis for those hyperparameters that have more than one unique value.
    


.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAA4QAAAJYCAYAAAA6xSjbAAAAAXNSR0IArs4c6QAAIABJREFUeF7snQd0VMXbh38hEKrSO6EECL33EjpIV5AiTUSa0kEQRUCayh8pSlVERbpSQocoUgJIk470JgEU6SC95Dvv5Nu4CZu2uXf2Zvc35+RIsve+M/PM7HqffefO9QoNDQ0FCwmQAAmQAAmQAAmQAAmQAAmQgMcR8KIQetyYs8MkQAIkQAIkQAIkQAIkQAIkoAhQCDkRSIAESIAESIAESIAESIAESMBDCVAIPXTg2W0SIAESIAESIAESIAESIAESoBByDpAACZAACZAACZAACZAACZCAhxKgEHrowLPbJEACJEACJEACJEACJEACJEAh5BwgARIgARIgARIgARIgARIgAQ8lQCH00IFnt0mABEiABEiABEiABEiABEiAQsg5QAIkQAIkQAIkQAIkQAIkQAIeSoBC6KEDz26TAAmQAAmQAAmQAAmQAAmQAIWQc4AESIAESIAESIAESIAESIAEPJQAhdBDB57dJgESIAESIAESIAESIAESIAEKIecACZAACZAACZAACZAACZAACXgoAQqhhw48u00CJEACJEACJEACJEACJEACFELOARIgARIgARIgARIgARIgARLwUAIUQg8deHabBEiABEiABEiABEiABEiABCiEnAMkQAIkQAIkQAIkQAIkQAIk4KEEKIQeOvDsNgmQAAmQAAmQAAmQAAmQAAlQCDkHSIAESIAESIAESIAESIAESMBDCVAIPXTg2W0SIAESIAESIAESIAESIAESoBByDpAACZAACZAACZAACZAACZCAhxKgEHrowLPbJEACJEACJEACJEACJEACJEAh5BwgARIgARIgARIgARIgARIgAQ8lQCH00IFnt0mABEiABEiABEiABEiABEiAQsg5QAIkQAIkQAIkQAIkQAIkQAIeSoBC6KEDz26TAAmQAAmQAAmQAAmQAAmQAIWQc4AESIAESIAESIAESIAESIAEPJQAhdBDB57dJgESIAESIAESIAESIAESIAEKIecACZAACZAACZAACZAACZAACXgoAQqhhw48u00CJEACJEACJEACJEACJEACFELOARIgARIgARIgARIgARIgARLwUAIUQg8deHabBEiABEiABEiABEiABEiABCiEnAMkQAIkQAIkQAIkQAIkQAIk4KEEKIQeOvDsNgmQAAmQAAmQAAmQAAmQAAlQCDkHSIAESIAESIAESIAESIAESMBDCVAIPXTg2W0SIAESIAESIAESIAESIAESoBByDpAACZAACZAACZAACZAACZCAhxKgEHrowLPbJEACJEACJEACJEACJEACJEAh5BwgARIgARIgARIgARIgARIgAQ8lQCH00IFnt0mABEiABEiABEiABEiABEiAQsg5QAIkQAIkQAIkQAIkQAIkQAIeSoBC6KEDz26TAAmQAAmQAAmQAAmQAAmQAIWQc4AESIAESIAESIAESIAESIAEPJQAhdBDB57dJgESIAESIAESIAESIAESIAEKIecACZAACZAACZAACZAACZAACXgoAQqhhw48u00CJEACJEACJEACJEACJEACFELOARIgARIgARIgARIgARIgARLwUAIUQg8deHabBEiABEiABEiABEiABEiABCiEnAMkQAIkQAIkQAIkQAIkQAIk4KEEKIQeOvDsNgmQAAmQAAmQAAmQAAmQAAlQCDkHSIAESIAESIAESIAESIAESMBDCVAIPXTg2W0SIAESIAESIAESIAESIAESoBByDpAACZAACZAACZAACZAACZCAhxKgEHrowLPbJEACJEACJEACJEACJEACJEAh5BwgARIgARIgARIgARIgARIgAQ8lQCH00IFnt0mABEiABEiABEiABEiABEiAQsg5QAIJkICXlxf279+PkiVLOtX6IkWK4H//+x8aN27s1PnxOWnVqlXo3bs3rl27hnnz5uG1116LTziem8AJzJ8/H9OmTcNvv/0WZU/iO9+dRRSbtjkbO6Gct3PnTrz99tsICQnBJ598gj59+iSUprOdJEACJEACsSRAIYwlKB5GAjYC27ZtUxdGcqEUGhqKXLlyoV27dujXrx98fHy0gIrLBXKNGjWUdEn7rFDy5cuHjz/+GB06dHDYnLfeegtp0qTBF198EeH13Llzq79RIJ0bxYTCz1E7Yzvft27digYNGoQDunfvHpInT45EiRKpvw0ZMkT9mFkuXbqEEiVKIGPGjEqipDx69AjSB9vng3xm/PHHH3FuRmzGUI65cuUKvL29kSxZMlSqVEm9b/LmzRvn+uSEOnXqoHr16hg2bJhT5/MkEiABEiAB6xOgEFp/jNhCCxFYvXo12rRpg9GjR6N9+/bI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" width="1200">



.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <img src="data:image/png;base64,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" width="1500">


