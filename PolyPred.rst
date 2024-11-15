PolyPred
========

In this notebook, we will use PolyPred to calculate the PRS.

The installation process for PolyPred is well explained on their `GitHub
Repository <https://github.com/omerwe/polyfun>`__. PolyPred exploits
fine-mapping to improve cross-population polygenic risk scores by
predicting with causal effect estimates instead of tagging effect
estimates.

One needs to create a separate conda environment to execute this tool.
If you have already installed conda, use the following steps:

.. code:: bash

   git clone https://github.com/omerwe/polyfun
   cd polyfun
   conda env create -f polyfun.yml
   conda activate polyfun

The installation process is also available
`here <https://github.com/omerwe/polyfun>`__.

After installation, we followed the tutorial on GitHub for PRS
calculation, which can be found
`here <https://github.com/omerwe/polyfun/wiki/6.-Trans-ethnic-polygenic-risk-prediction-with-PolyPred>`__.

You also need to download **BOLT** to calculate SNP tagging. Please
refer to the BOLT Notebook for the installation process of BOLT.

Important Notes regarding executing PolyPred
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For fine mapping, PolyPred downloads the LD files from Broad Institute.
The process must be executed separately for each fold and phenotype
(dataset), which is time-consuming. Once the batch file of jobs is
created for a specific dataset, it can be executed separately. This
process can take 2-3 days to complete. The tool has two primary issues:

1. It is time-consuming, and even on HPC, sending multiple requests can
   slow down execution. If jobs are run for all datasets, some LD panels
   might be missed.

2. Although this tool optimizes the PRS score, including it in the final
   prediction or classification can sometimes decrease performance. The
   PRS score calculated by this tool should be used directly for
   prediction or classification. Here’s the provided content in Markdown
   format:

3. ``polyfunoutput+"polypred_bash.sh"`` Execute the file at the
   following path separately to execute the fine mapping jobs required
   by PolyPred for calculations. The following code snippet shows the
   content of the ``polypred_bash.sh`` file:

.. code:: python

   job_file = "bash "+ polyfunoutput+'jobs.txt'

   bash_script_content = f"""#!/bin/bash
   #SBATCH --job-name=ttemp
   #SBATCH --nodes=1
   #SBATCH --partition=ascher
   #SBATCH --time=200:00:00
   #SBATCH --output=job.%A_%a.out
   #SBATCH --error=job.%A_%a.err
   #SBATCH --array=1-1

   {job_file}
   """

   # File name to write the content
   file_name = polyfunoutput+"polypred_bash.sh"

   # Write the content to the file
   with open(file_name, "w") as file:
       for line in bash_script_content.split('\n'):
           file.write(line.lstrip() + '\n')
   print(f"Content written to {file_name}")

   #os.system("sbatch "+ file_name)

4. Execute the file using the following command. The path of
   ``polyfunoutput`` will be different for each phenotype and fold. For
   example, in the case of ``SampleData1``, it was
   ``traindirec+os.sep+"polyfun/"`` or ``SampleData1/Fold_0/polyfun``.

.. code:: bash

   bash polyfunoutput+"polypred_bash.sh"

GWAS File Processing for PolyPred
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PolyPred accepts GWAS in multiple formats, and the GWAS file contains
the column PolyPred requires for PRS calculation.

.. code:: ipython2

    import os
    import pandas as pd
    import numpy as np
    import sys
    
    #filedirec = sys.argv[1]
    
    filedirec = "SampleData1"
    #filedirec = "asthma"
    #filedirec = "irritable_bowel_syndrome"
    
    
    # Check the status of the Phenotype. 
    # If Phenotype is Binary, Plink requires OR, and if it continous, it requires BETAS in the GWAS file.
    
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
    
    N = df["N"].mean()
    N = int(N)
    print(N)
    
    df.to_csv(filedirec + os.sep +filedirec+".txt",sep="\t",index=False)
    
    df.to_csv(filedirec + os.sep +filedirec+"_Plink.txt",sep="\t",index=False)
    
    print(df.head().to_markdown())
    print("Length of DataFrame!",len(df))
    


.. parsed-literal::

    388028
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

.. code:: ipython2

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
    
    print("Minimum P-value",allpvalues[0])
    print("Maximum P-value",allpvalues[-1])
    print("Number of P-value",len(allpvalues))
    
    
    
    
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


.. parsed-literal::

    Minimum P-value 1e-10
    Maximum P-value 1.0
    Number of P-value 20
    

Define Helper Functions
~~~~~~~~~~~~~~~~~~~~~~~

1. **Perform Clumping and Pruning**
2. **Calculate PCA Using Plink**
3. **Fit Binary Phenotype and Save Results**
4. **Fit Continuous Phenotype and Save Results**

.. code:: ipython2

    import os
    import subprocess
    import pandas as pd
    import statsmodels.api as sm
    from sklearn.metrics import explained_variance_score
    
    
    # The following function is used to perform clumping and pruing on the Genotype Data.
    # It is almost the same in for all the phenotypes.
    
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
        
        
    # This function calculats PCA for test and train genotype data, which is then Appended with covariates, PRS, for prediction.
    
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
    
    # This function used for prediction of binary phenotypes.
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
        
        
        
        # One can use multiple ranges of alphas and weights for Binary Phenotype.
        tempalphas = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
        l1weights = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    
        tempalphas = [0.1]
        l1weights = [0.1]
    
        # The following transformation is required by Logit.
        phenotype_train["Phenotype"] = phenotype_train["Phenotype"].replace({1: 0, 2: 1}) 
        phenotype_test["Phenotype"] = phenotype_test["Phenotype"].replace({1: 0, 2: 1})
          
        for tempalpha in tempalphas:
            for l1weight in l1weights:
    
                
                try:
                    # Fit the null model.
                    null_model =  sm.Logit(phenotype_train["Phenotype"], sm.add_constant(covandpcs_train.iloc[:, 2:])).fit_regularized(alpha=tempalpha, L1_wt=l1weight)
                    #null_model =  sm.Logit(phenotype_train["Phenotype"], sm.add_constant(covandpcs_train.iloc[:, 2:])).fit()
                
                except:
                    print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
                    continue
                
                # Prediction using null model on training data.
                train_null_predicted = null_model.predict(sm.add_constant(covandpcs_train.iloc[:, 2:]))
                
                from sklearn.metrics import roc_auc_score, confusion_matrix
                from sklearn.metrics import r2_score
                
                # Prediction using null model on test data.
                test_null_predicted = null_model.predict(sm.add_constant(covandpcs_test.iloc[:, 2:]))
            
           
            
          
                # Get the PRS for a specific range.
                prs_train = pd.read_table(traindirec+os.sep+'Train_PRS_polypred.prs', sep="\s+", usecols=["FID", "IID", "PRS"])
             
    
                prs_train['FID'] = prs_train['FID'].astype(str)
                prs_train['IID'] = prs_train['IID'].astype(str)
               
                prs_test = pd.read_table(traindirec+os.sep+'Test_PRS_polypred.prs', sep="\s+", usecols=["FID", "IID", "PRS"])
        
                prs_test['FID'] = prs_test['FID'].astype(str)
                prs_test['IID'] = prs_test['IID'].astype(str)
                # Append PRS with covaraites and PCA for both datasets.
                pheno_prs_train = pd.merge(covandpcs_train, prs_train, on=["FID", "IID"])
                pheno_prs_test = pd.merge(covandpcs_test, prs_test, on=["FID", "IID"])
        
                try:
                    # Fit the model.
                    model = sm.Logit(phenotype_train["Phenotype"], sm.add_constant(pheno_prs_train.iloc[:, 2:])).fit_regularized(alpha=tempalpha, L1_wt=l1weight)
                    #model = sm.Logit(phenotype_train["Phenotype"], sm.add_constant(pheno_prs_train.iloc[:, 2:])).fit()
                
                except:
                    continue
    
                global prs_result 
                
                train_best_predicted = model.predict(sm.add_constant(pheno_prs_train.iloc[:, 2:]))    
                test_best_predicted = model.predict(sm.add_constant(pheno_prs_test.iloc[:, 2:])) 
    
        
                from sklearn.metrics import roc_auc_score, confusion_matrix
                
                # Save the information for all hyperparameters.
                prs_result = prs_result._append({
                    "clump_p1": c1_val,
                    "clump_r2": c2_val,
                    "clump_kb": c3_val,
                    "p_window_size": p1_val,
                    "p_slide_size": p2_val,
                    "p_LD_threshold": p3_val,
                    "pvalue":str(1),
                    "numberofpca":p, 
    
                    "tempalpha":str(tempalpha),
                    "l1weight":str(l1weight),
                     
                     
                    # Here one can use the different evaluation metrix.
                    # Same goes for continous phenotype in the following function.
                    "Train_pure_prs":roc_auc_score(phenotype_train["Phenotype"].values,prs_train['PRS'].values),
                    "Train_null_model":roc_auc_score(phenotype_train["Phenotype"].values,train_null_predicted.values),
                    "Train_best_model":roc_auc_score(phenotype_train["Phenotype"].values,train_best_predicted.values),
                    
                    "Test_pure_prs":roc_auc_score(phenotype_test["Phenotype"].values,prs_test['PRS'].values),
                    "Test_null_model":roc_auc_score(phenotype_test["Phenotype"].values,test_null_predicted.values),
                    "Test_best_model":roc_auc_score(phenotype_test["Phenotype"].values,test_best_predicted.values),
                    
                }, ignore_index=True)
    
          
                prs_result.to_csv(traindirec+os.sep+Name+os.sep+"Results.csv",index=False)
         
        return
    
    # This function used for prediction of continous phenotypes.
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
    
                 
                # If explained variance is negative, use fit regualized rather than fit.
                #null_model =  sm.OLS(phenotype_train["Phenotype"], sm.add_constant(covandpcs_train.iloc[:, 2:])).fit_regularized(alpha=tempalpha, L1_wt=l1weight)
                null_model =  sm.OLS(phenotype_train["Phenotype"], sm.add_constant(covandpcs_train.iloc[:, 2:])).fit()
                #null_model =  sm.OLS(phenotype_train["Phenotype"], sm.add_constant(covandpcs_train.iloc[:, 2:])).fit()
          
                train_null_predicted = null_model.predict(sm.add_constant(covandpcs_train.iloc[:, 2:]))
                
                from sklearn.metrics import roc_auc_score, confusion_matrix
                from sklearn.metrics import r2_score
                
                test_null_predicted = null_model.predict(sm.add_constant(covandpcs_test.iloc[:, 2:]))
                
                
                
                global prs_result 
          
                prs_train = pd.read_table(traindirec+os.sep+'Train_PRS_polypred.prs', sep="\s+", usecols=["FID", "IID", "PRS"])
          
    
                prs_train['FID'] = prs_train['FID'].astype(str)
                prs_train['IID'] = prs_train['IID'].astype(str)
                
                
                prs_test = pd.read_table(traindirec+os.sep+'Test_PRS_polypred.prs', sep="\s+", usecols=["FID", "IID", "PRS"])
                prs_test['FID'] = prs_test['FID'].astype(str)
                prs_test['IID'] = prs_test['IID'].astype(str)
                pheno_prs_train = pd.merge(covandpcs_train, prs_train, on=["FID", "IID"])
                pheno_prs_test = pd.merge(covandpcs_test, prs_test, on=["FID", "IID"])
         
                #model = sm.OLS(phenotype_train["Phenotype"], sm.add_constant(pheno_prs_train.iloc[:, 2:])).fit_regularized(alpha=tempalpha, L1_wt=l1weight)
                model = sm.OLS(phenotype_train["Phenotype"], sm.add_constant(pheno_prs_train.iloc[:, 2:])).fit()
          
    
                
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
                    "pvalue":str(1),
                    "tempalpha":str(tempalpha),
                    "l1weight":str(l1weight),
                      
    
                    "Train_pure_prs":explained_variance_score(phenotype_train["Phenotype"],prs_train['PRS'].values),
                    "Train_null_model":explained_variance_score(phenotype_train["Phenotype"],train_null_predicted),
                    "Train_best_model":explained_variance_score(phenotype_train["Phenotype"],train_best_predicted),
                    
                    "Test_pure_prs":explained_variance_score(phenotype_test["Phenotype"],prs_test['PRS'].values),
                    "Test_null_model":explained_variance_score(phenotype_test["Phenotype"],test_null_predicted),
                    "Test_best_model":explained_variance_score(phenotype_test["Phenotype"],test_best_predicted),
                    
                }, ignore_index=True)
    
          
                prs_result.to_csv(traindirec+os.sep+Name+os.sep+"Results.csv",index=False)
     
        return
    

Execute PolyPred
----------------

.. code:: ipython2

    
    # Define a global variable to store results
    prs_result = pd.DataFrame()
    def transform_polypred_data(traindirec, newtrainfilename,p, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile):
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
     
        import subprocess
    
        # Define the command as a list
        command = [
            'python',
            'polyfun/munge_polyfun_sumstats.py',
            '--sumstats', filedirec+os.sep+filedirec+".gz",
            '--out', traindirec+os.sep+filedirec+".temppolyfun.parquet",
        ]
    
        subprocess.run(command)
        
        create_directory(traindirec+os.sep+"polyfun")
        
        polyfunoutput = traindirec+os.sep+"polyfun/"
        command = [
            'python',
            'polyfun/create_finemapper_jobs.py',
            '--sumstats', traindirec+os.sep+filedirec+".temppolyfun.parquet",
            '--method', 'susie',
            '--n',str(N),
            '--non-funct',
       
            '--max-num-causal','5',
            '--allow-missing',
            
            '--out-prefix', polyfunoutput+'polyfun_output',
            '--jobs-file', polyfunoutput+'jobs.txt'
        ]
        print(" ".join(command))
        subprocess.run(command)
        
        job_file = "bash "+ polyfunoutput+'jobs.txt'
        
        bash_script_content = f"""#!/bin/bash
        #SBATCH --job-name=ttemp
        #SBATCH --nodes=1
        #SBATCH --partition=ascher
        #SBATCH --time=200:00:00
        #SBATCH --output=job.%A_%a.out
        #SBATCH --error=job.%A_%a.err
        #SBATCH --array=1-1
     
        {job_file}
        """
    
        # File name to write the content
        file_name = polyfunoutput+"polypred_bash.sh"
    
        # Write the content to the file
        with open(file_name, "w") as file:
            for line in bash_script_content.split('\n'):
                file.write(line.lstrip() + '\n')
        print(f"Content written to {file_name}")
    
        #os.system("sbatch "+ file_name)
    
        
        
        
        # Kindly note that first execute the bash file polyfunoutput+"polypred_bash.sh", and after that execute the following code. 
        command = [
            'python',
            'polyfun/aggregate_finemapper_results.py',
            '--out-prefix', polyfunoutput+'polyfun_output',
            '--sumstats', traindirec+os.sep+filedirec+".temppolyfun.parquet",
            
            # If some of the jobs, failed, use the following tag.
            '--allow-missing-jobs',
            
            '--out', traindirec+os.sep+filedirec+".finalpolypred",
            '--adjust-beta-freq'
        ]
    
        # Run the command
        subprocess.run(command)
        
        
    
        # Calculate GWAS using BOlT
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
    
        mm = "--lmmForceNonInf"
        
        command = [
            
        './bolt',
        '--bfile='+traindirec+os.sep+newtrainfilename+".clumped.pruned",
        '--phenoFile='+traindirec+os.sep+trainfilename+".PHENO" ,
        '--phenoCol=PHENO',
        mm,
        '--LDscoresFile=tables/LDSCORE.1000G_EUR.tab.gz',
        '--covarFile='+traindirec+os.sep+trainfilename+".COV_PCA",
        '--LDscoresMatchBp',
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
    
        command = [
            'python', 'polyfun/polypred.py',
            '--estimate-mixweights',
            '--betas', traindirec+os.sep+filedirec+"."+mm.replace("-","")+"_stat"+","+traindirec+os.sep+filedirec+".finalpolypred",
            '--pheno', traindirec+os.sep+trainfilename+".PHENO",
            '--output-prefix', traindirec+os.sep+'GWAS_polypred',
            '--plink-exe', './plink',
            traindirec+os.sep+newtrainfilename+".clumped.pruned"+".bed"
            #traindirec+os.sep+newtrainfilename+".bed"
            
        ]
        
        # Run the command
        subprocess.run(command)
        #raise
        
        command = [
             'python', 'polyfun/polypred.py',
            '--predict',
            '--betas',  traindirec+os.sep+filedirec+"."+mm.replace("-","")+"_stat"+","+traindirec+os.sep+filedirec+".finalpolypred",
            '--mixweights-prefix',  traindirec+os.sep+'GWAS_polypred',
            '--output-prefix', traindirec+os.sep+'Train_PRS_polypred',
            '--plink-exe', './plink',
            traindirec+os.sep+newtrainfilename+".clumped.pruned"+".bed"
        ]
        
        # Run the command
        subprocess.run(command)    
        
        command = [
             'python', 'polyfun/polypred.py',
            '--predict',
            '--betas',  traindirec+os.sep+filedirec+"."+mm.replace("-","")+"_stat"+","+traindirec+os.sep+filedirec+".finalpolypred",
            '--mixweights-prefix',  traindirec+os.sep+'GWAS_polypred',
            '--output-prefix', traindirec+os.sep+'Test_PRS_polypred',
            '--plink-exe', './plink',
            traindirec+os.sep+testfilename+".bed"
        ]
        
        # Run the command
        subprocess.run(command)    
        
     
        # At this stage the scores are finalizied. 
        # The next step is to fit the model 
        
        if check_phenotype_is_binary_or_continous(filedirec)=="Binary":
            print("Binary Phenotype!")
            fit_binary_phenotype_on_PRS(traindirec, newtrainfilename,p, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile)
        else:
            print("Continous Phenotype!")
            fit_continous_phenotype_on_PRS(traindirec, newtrainfilename,p, p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name,pvaluefile)
                
     
    
     
    
    result_directory = "PolyPred"
    # Nested loops to iterate over different parameter values
    # Create a directory to save the results.
    create_directory(folddirec+os.sep+result_directory)
    for p1_val in p_window_size:
     for p2_val in p_slide_size: 
      for p3_val in p_LD_threshold:
       for c1_val in clump_p1:
        for c2_val in clump_r2:
         for c3_val in clump_kb:
          for p in numberofpca:
            #pass
            transform_polypred_data(folddirec, newtrainfilename, p, str(p1_val), str(p2_val), str(p3_val), str(c1_val), str(c2_val), str(c3_val), result_directory, pvaluefile)
    


.. parsed-literal::

    [INFO]  Reading sumstats file...
    [INFO]  Done in 0.93 seconds
    [INFO]  Converting OR column to log-odds
    [INFO]  499617 SNPs are in the sumstats file
    [INFO]  Removing 5737 HLA SNPs
    [INFO]  493880 SNPs with sumstats remained after all filtering stages
    [INFO]  Saving munged sumstats of 493880 SNPs to SampleData1/Fold_0/SampleData1.temppolyfun.parquet
    [INFO]  Done
    python polyfun/create_finemapper_jobs.py --sumstats SampleData1/Fold_0/SampleData1.temppolyfun.parquet --method susie --n 388028 --non-funct --max-num-causal 5 --allow-missing --out-prefix SampleData1/Fold_0/polyfun/polyfun_output --jobs-file SampleData1/Fold_0/polyfun/jobs.txt
    [INFO]  Wrote fine-mapping commands to SampleData1/Fold_0/polyfun/jobs.txt
    Content written to SampleData1/Fold_0/polyfun/polypred_bash.sh
    [INFO]  Aggregating results...
    

.. parsed-literal::

    632it [00:00, 2041.28it/s]

.. parsed-literal::

    [WARNING]  output file for chromosome 1 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 91000001-94000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 92000001-95000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 93000001-96000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 94000001-97000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 95000001-98000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 96000001-99000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 97000001-100000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 98000001-101000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 99000001-102000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 100000001-103000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 101000001-104000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 102000001-105000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 103000001-106000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 104000001-107000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 105000001-108000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 106000001-109000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 107000001-110000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 108000001-111000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 109000001-112000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 110000001-113000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 111000001-114000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 112000001-115000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 113000001-116000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 114000001-117000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 115000001-118000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 116000001-119000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 117000001-120000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 118000001-121000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 119000001-122000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 120000001-123000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 121000001-124000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 142000001-145000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 143000001-146000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 144000001-147000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 145000001-148000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 146000001-149000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 147000001-150000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 148000001-151000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 149000001-152000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 150000001-153000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 151000001-154000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 152000001-155000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 153000001-156000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 154000001-157000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 155000001-158000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 156000001-159000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 157000001-160000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 158000001-161000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 159000001-162000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 160000001-163000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 161000001-164000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 162000001-165000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 163000001-166000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 164000001-167000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 165000001-168000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 166000001-169000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 167000001-170000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 168000001-171000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 169000001-172000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 170000001-173000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 171000001-174000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 172000001-175000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 173000001-176000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 174000001-177000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 175000001-178000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 176000001-179000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 177000001-180000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 178000001-181000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 179000001-182000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 180000001-183000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 181000001-184000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 182000001-185000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 183000001-186000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 184000001-187000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 185000001-188000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 186000001-189000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 187000001-190000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 188000001-191000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 189000001-192000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 190000001-193000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 191000001-194000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 192000001-195000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 193000001-196000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 194000001-197000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 195000001-198000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 196000001-199000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 197000001-200000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 198000001-201000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 199000001-202000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 200000001-203000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 201000001-204000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 202000001-205000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 203000001-206000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 204000001-207000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 205000001-208000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 206000001-209000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 207000001-210000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 208000001-211000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 209000001-212000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 210000001-213000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 211000001-214000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 212000001-215000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 213000001-216000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 214000001-217000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 215000001-218000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 216000001-219000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 217000001-220000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 218000001-221000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 219000001-222000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 220000001-223000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 221000001-224000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 222000001-225000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 223000001-226000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 224000001-227000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 225000001-228000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 226000001-229000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 227000001-230000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 228000001-231000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 229000001-232000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 230000001-233000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 231000001-234000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 232000001-235000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 233000001-236000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 234000001-237000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 235000001-238000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 236000001-239000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 237000001-240000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 238000001-241000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 239000001-242000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 240000001-243000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 241000001-244000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 242000001-245000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 243000001-246000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 244000001-247000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 245000001-248000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 246000001-249000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 247000001-250000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 248000001-251000001 doesn't exist
    [WARNING]  output file for chromosome 1 bp 249000001-252000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 1000001-4000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 7000001-10000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 8000001-11000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 9000001-12000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 91000001-94000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 92000001-95000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 93000001-96000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 94000001-97000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 95000001-98000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 96000001-99000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 97000001-100000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 98000001-101000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 99000001-102000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 100000001-103000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 101000001-104000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 102000001-105000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 103000001-106000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 104000001-107000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 105000001-108000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 106000001-109000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 107000001-110000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 108000001-111000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 109000001-112000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 110000001-113000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 111000001-114000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 112000001-115000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 113000001-116000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 114000001-117000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 115000001-118000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 116000001-119000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 117000001-120000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 118000001-121000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 119000001-122000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 120000001-123000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 121000001-124000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 122000001-125000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 123000001-126000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 124000001-127000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 125000001-128000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 126000001-129000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 127000001-130000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 128000001-131000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 129000001-132000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 130000001-133000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 131000001-134000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 132000001-135000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 133000001-136000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 134000001-137000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 135000001-138000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 136000001-139000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 137000001-140000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 138000001-141000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 139000001-142000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 140000001-143000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 141000001-144000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 142000001-145000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 143000001-146000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 144000001-147000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 145000001-148000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 146000001-149000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 147000001-150000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 148000001-151000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 149000001-152000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 150000001-153000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 151000001-154000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 152000001-155000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 153000001-156000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 154000001-157000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 155000001-158000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 156000001-159000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 157000001-160000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 158000001-161000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 159000001-162000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 160000001-163000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 161000001-164000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 162000001-165000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 163000001-166000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 164000001-167000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 165000001-168000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 166000001-169000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 167000001-170000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 168000001-171000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 169000001-172000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 170000001-173000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 171000001-174000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 172000001-175000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 173000001-176000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 174000001-177000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 175000001-178000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 176000001-179000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 177000001-180000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 178000001-181000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 179000001-182000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 180000001-183000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 181000001-184000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 182000001-185000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 183000001-186000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 184000001-187000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 185000001-188000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 186000001-189000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 187000001-190000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 188000001-191000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 189000001-192000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 190000001-193000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 191000001-194000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 192000001-195000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 193000001-196000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 194000001-197000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 195000001-198000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 196000001-199000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 197000001-200000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 198000001-201000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 199000001-202000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 200000001-203000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 201000001-204000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 202000001-205000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 203000001-206000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 204000001-207000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 205000001-208000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 206000001-209000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 207000001-210000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 208000001-211000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 209000001-212000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 210000001-213000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 211000001-214000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 212000001-215000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 213000001-216000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 214000001-217000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 215000001-218000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 216000001-219000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 217000001-220000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 218000001-221000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 219000001-222000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 220000001-223000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 221000001-224000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 222000001-225000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 223000001-226000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 224000001-227000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 225000001-228000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 226000001-229000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 227000001-230000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 228000001-231000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 229000001-232000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 230000001-233000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 231000001-234000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 232000001-235000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 233000001-236000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 234000001-237000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 235000001-238000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 236000001-239000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 237000001-240000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 238000001-241000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 239000001-242000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 240000001-243000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 241000001-244000001 doesn't exist
    [WARNING]  output file for chromosome 2 bp 242000001-245000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 1000001-4000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 7000001-10000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 8000001-11000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 9000001-12000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 91000001-94000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 92000001-95000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 93000001-96000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 94000001-97000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 95000001-98000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 96000001-99000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 97000001-100000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 98000001-101000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 99000001-102000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 100000001-103000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 101000001-104000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 102000001-105000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 103000001-106000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 104000001-107000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 105000001-108000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 106000001-109000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 107000001-110000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 108000001-111000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 109000001-112000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 110000001-113000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 111000001-114000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 112000001-115000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 113000001-116000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 114000001-117000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 115000001-118000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 116000001-119000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 117000001-120000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 118000001-121000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 119000001-122000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 120000001-123000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 121000001-124000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 122000001-125000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 123000001-126000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 124000001-127000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 125000001-128000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 126000001-129000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 127000001-130000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 128000001-131000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 129000001-132000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 130000001-133000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 131000001-134000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 132000001-135000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 133000001-136000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 134000001-137000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 135000001-138000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 136000001-139000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 137000001-140000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 138000001-141000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 139000001-142000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 140000001-143000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 141000001-144000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 142000001-145000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 143000001-146000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 144000001-147000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 145000001-148000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 146000001-149000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 147000001-150000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 148000001-151000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 149000001-152000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 150000001-153000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 151000001-154000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 152000001-155000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 153000001-156000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 154000001-157000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 155000001-158000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 156000001-159000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 157000001-160000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 158000001-161000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 159000001-162000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 160000001-163000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 161000001-164000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 162000001-165000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 163000001-166000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 164000001-167000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 165000001-168000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 166000001-169000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 167000001-170000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 168000001-171000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 169000001-172000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 170000001-173000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 171000001-174000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 172000001-175000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 173000001-176000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 174000001-177000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 175000001-178000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 176000001-179000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 177000001-180000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 178000001-181000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 179000001-182000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 180000001-183000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 181000001-184000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 182000001-185000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 183000001-186000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 184000001-187000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 185000001-188000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 186000001-189000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 187000001-190000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 188000001-191000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 189000001-192000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 190000001-193000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 191000001-194000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 192000001-195000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 193000001-196000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 194000001-197000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 195000001-198000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 196000001-199000001 doesn't exist
    [WARNING]  output file for chromosome 3 bp 197000001-200000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 1000001-4000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 7000001-10000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 8000001-11000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 9000001-12000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 91000001-94000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 92000001-95000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 93000001-96000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 94000001-97000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 95000001-98000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 96000001-99000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 97000001-100000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 98000001-101000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 99000001-102000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 100000001-103000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 101000001-104000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 102000001-105000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 103000001-106000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 104000001-107000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 105000001-108000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 106000001-109000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 107000001-110000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 108000001-111000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 109000001-112000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 110000001-113000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 111000001-114000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 112000001-115000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 113000001-116000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 114000001-117000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 115000001-118000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 116000001-119000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 117000001-120000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 118000001-121000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 119000001-122000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 120000001-123000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 121000001-124000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 122000001-125000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 123000001-126000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 124000001-127000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 125000001-128000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 126000001-129000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 127000001-130000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 128000001-131000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 129000001-132000001 doesn't exist
    

.. parsed-literal::

    956it [00:00, 2470.11it/s]

.. parsed-literal::

    [WARNING]  output file for chromosome 4 bp 130000001-133000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 131000001-134000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 132000001-135000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 133000001-136000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 134000001-137000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 135000001-138000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 136000001-139000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 137000001-140000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 138000001-141000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 139000001-142000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 140000001-143000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 141000001-144000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 142000001-145000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 143000001-146000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 144000001-147000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 145000001-148000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 146000001-149000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 147000001-150000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 148000001-151000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 149000001-152000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 150000001-153000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 151000001-154000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 152000001-155000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 153000001-156000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 154000001-157000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 155000001-158000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 156000001-159000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 157000001-160000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 158000001-161000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 159000001-162000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 160000001-163000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 161000001-164000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 162000001-165000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 163000001-166000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 164000001-167000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 165000001-168000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 166000001-169000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 167000001-170000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 168000001-171000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 169000001-172000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 170000001-173000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 171000001-174000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 172000001-175000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 173000001-176000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 174000001-177000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 175000001-178000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 176000001-179000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 177000001-180000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 178000001-181000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 179000001-182000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 180000001-183000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 181000001-184000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 182000001-185000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 183000001-186000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 184000001-187000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 185000001-188000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 186000001-189000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 187000001-190000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 188000001-191000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 189000001-192000001 doesn't exist
    [WARNING]  output file for chromosome 4 bp 190000001-193000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 1000001-4000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 7000001-10000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 8000001-11000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 9000001-12000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 91000001-94000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 92000001-95000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 93000001-96000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 94000001-97000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 95000001-98000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 96000001-99000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 97000001-100000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 98000001-101000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 99000001-102000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 100000001-103000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 101000001-104000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 102000001-105000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 103000001-106000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 104000001-107000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 105000001-108000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 106000001-109000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 107000001-110000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 108000001-111000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 109000001-112000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 110000001-113000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 111000001-114000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 112000001-115000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 113000001-116000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 114000001-117000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 115000001-118000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 116000001-119000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 117000001-120000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 118000001-121000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 119000001-122000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 120000001-123000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 121000001-124000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 122000001-125000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 123000001-126000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 124000001-127000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 125000001-128000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 126000001-129000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 127000001-130000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 128000001-131000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 129000001-132000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 130000001-133000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 131000001-134000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 132000001-135000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 133000001-136000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 134000001-137000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 135000001-138000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 136000001-139000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 137000001-140000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 138000001-141000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 139000001-142000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 140000001-143000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 141000001-144000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 142000001-145000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 143000001-146000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 144000001-147000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 145000001-148000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 146000001-149000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 147000001-150000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 148000001-151000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 149000001-152000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 150000001-153000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 151000001-154000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 152000001-155000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 153000001-156000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 154000001-157000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 155000001-158000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 156000001-159000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 157000001-160000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 158000001-161000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 159000001-162000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 160000001-163000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 161000001-164000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 162000001-165000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 163000001-166000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 164000001-167000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 165000001-168000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 166000001-169000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 167000001-170000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 168000001-171000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 169000001-172000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 170000001-173000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 171000001-174000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 172000001-175000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 173000001-176000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 174000001-177000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 175000001-178000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 176000001-179000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 177000001-180000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 178000001-181000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 179000001-182000001 doesn't exist
    [WARNING]  output file for chromosome 5 bp 180000001-183000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 1000001-4000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 7000001-10000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 8000001-11000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 9000001-12000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 91000001-94000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 92000001-95000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 93000001-96000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 94000001-97000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 95000001-98000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 96000001-99000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 97000001-100000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 98000001-101000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 99000001-102000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 100000001-103000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 101000001-104000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 102000001-105000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 103000001-106000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 104000001-107000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 105000001-108000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 106000001-109000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 107000001-110000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 108000001-111000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 109000001-112000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 110000001-113000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 111000001-114000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 112000001-115000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 113000001-116000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 114000001-117000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 115000001-118000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 116000001-119000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 117000001-120000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 118000001-121000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 119000001-122000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 120000001-123000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 121000001-124000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 122000001-125000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 123000001-126000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 124000001-127000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 125000001-128000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 126000001-129000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 127000001-130000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 128000001-131000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 129000001-132000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 130000001-133000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 131000001-134000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 132000001-135000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 133000001-136000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 134000001-137000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 135000001-138000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 136000001-139000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 137000001-140000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 138000001-141000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 139000001-142000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 140000001-143000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 141000001-144000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 142000001-145000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 143000001-146000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 144000001-147000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 145000001-148000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 146000001-149000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 147000001-150000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 148000001-151000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 149000001-152000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 150000001-153000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 151000001-154000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 152000001-155000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 153000001-156000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 154000001-157000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 155000001-158000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 156000001-159000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 157000001-160000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 158000001-161000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 159000001-162000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 160000001-163000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 161000001-164000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 162000001-165000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 163000001-166000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 164000001-167000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 165000001-168000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 166000001-169000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 167000001-170000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 168000001-171000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 169000001-172000001 doesn't exist
    [WARNING]  output file for chromosome 6 bp 170000001-173000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 1000001-4000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 7000001-10000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 8000001-11000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 9000001-12000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 91000001-94000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 92000001-95000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 93000001-96000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 94000001-97000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 95000001-98000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 96000001-99000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 97000001-100000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 98000001-101000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 99000001-102000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 100000001-103000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 101000001-104000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 102000001-105000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 103000001-106000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 104000001-107000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 105000001-108000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 106000001-109000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 107000001-110000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 108000001-111000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 109000001-112000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 110000001-113000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 111000001-114000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 112000001-115000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 113000001-116000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 114000001-117000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 115000001-118000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 116000001-119000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 117000001-120000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 118000001-121000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 119000001-122000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 120000001-123000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 121000001-124000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 122000001-125000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 123000001-126000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 124000001-127000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 125000001-128000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 126000001-129000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 127000001-130000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 128000001-131000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 129000001-132000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 130000001-133000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 131000001-134000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 132000001-135000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 133000001-136000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 134000001-137000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 135000001-138000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 136000001-139000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 137000001-140000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 138000001-141000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 139000001-142000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 140000001-143000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 141000001-144000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 142000001-145000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 143000001-146000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 144000001-147000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 145000001-148000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 146000001-149000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 147000001-150000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 148000001-151000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 149000001-152000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 150000001-153000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 151000001-154000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 152000001-155000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 153000001-156000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 154000001-157000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 155000001-158000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 156000001-159000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 157000001-160000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 158000001-161000001 doesn't exist
    [WARNING]  output file for chromosome 7 bp 159000001-162000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 1000001-4000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 7000001-10000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 91000001-94000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 92000001-95000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 93000001-96000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 94000001-97000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 95000001-98000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 96000001-99000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 97000001-100000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 98000001-101000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 99000001-102000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 100000001-103000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 101000001-104000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 102000001-105000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 103000001-106000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 104000001-107000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 105000001-108000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 106000001-109000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 107000001-110000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 108000001-111000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 109000001-112000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 110000001-113000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 111000001-114000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 112000001-115000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 113000001-116000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 114000001-117000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 115000001-118000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 116000001-119000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 117000001-120000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 118000001-121000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 119000001-122000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 120000001-123000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 121000001-124000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 122000001-125000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 123000001-126000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 124000001-127000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 125000001-128000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 126000001-129000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 127000001-130000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 128000001-131000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 129000001-132000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 130000001-133000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 131000001-134000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 132000001-135000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 133000001-136000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 134000001-137000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 135000001-138000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 136000001-139000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 137000001-140000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 138000001-141000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 139000001-142000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 140000001-143000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 141000001-144000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 142000001-145000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 143000001-146000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 144000001-147000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 145000001-148000001 doesn't exist
    [WARNING]  output file for chromosome 8 bp 146000001-149000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 1000001-4000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 7000001-10000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 8000001-11000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 9000001-12000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 91000001-94000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 92000001-95000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 93000001-96000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 94000001-97000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 95000001-98000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 96000001-99000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 97000001-100000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 98000001-101000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 99000001-102000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 100000001-103000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 101000001-104000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 102000001-105000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 103000001-106000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 104000001-107000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 105000001-108000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 106000001-109000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 107000001-110000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 108000001-111000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 109000001-112000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 110000001-113000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 111000001-114000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 112000001-115000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 113000001-116000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 114000001-117000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 115000001-118000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 116000001-119000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 117000001-120000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 118000001-121000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 119000001-122000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 120000001-123000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 121000001-124000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 122000001-125000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 123000001-126000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 124000001-127000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 125000001-128000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 126000001-129000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 127000001-130000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 128000001-131000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 129000001-132000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 130000001-133000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 131000001-134000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 132000001-135000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 133000001-136000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 134000001-137000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 135000001-138000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 136000001-139000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 137000001-140000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 138000001-141000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 139000001-142000001 doesn't exist
    [WARNING]  output file for chromosome 9 bp 140000001-143000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 1000001-4000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 7000001-10000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 8000001-11000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 9000001-12000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 91000001-94000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 92000001-95000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 93000001-96000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 94000001-97000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 95000001-98000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 96000001-99000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 97000001-100000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 98000001-101000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 99000001-102000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 100000001-103000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 101000001-104000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 102000001-105000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 103000001-106000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 104000001-107000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 105000001-108000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 106000001-109000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 107000001-110000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 108000001-111000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 109000001-112000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 110000001-113000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 111000001-114000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 112000001-115000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 113000001-116000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 114000001-117000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 115000001-118000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 116000001-119000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 117000001-120000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 118000001-121000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 119000001-122000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 120000001-123000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 121000001-124000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 122000001-125000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 123000001-126000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 124000001-127000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 125000001-128000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 126000001-129000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 127000001-130000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 128000001-131000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 129000001-132000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 130000001-133000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 131000001-134000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 132000001-135000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 133000001-136000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 134000001-137000001 doesn't exist
    [WARNING]  output file for chromosome 10 bp 135000001-138000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 1000001-4000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 7000001-10000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 8000001-11000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 9000001-12000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 38000001-41000001 doesn't ex

.. parsed-literal::

    2                          t/s]

.. parsed-literal::

    ist
    [WARNING]  output file for chromosome 11 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 91000001-94000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 92000001-95000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 93000001-96000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 94000001-97000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 95000001-98000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 96000001-99000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 97000001-100000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 98000001-101000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 99000001-102000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 100000001-103000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 101000001-104000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 102000001-105000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 103000001-106000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 104000001-107000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 105000001-108000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 106000001-109000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 107000001-110000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 108000001-111000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 109000001-112000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 110000001-113000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 111000001-114000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 112000001-115000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 113000001-116000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 114000001-117000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 115000001-118000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 116000001-119000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 117000001-120000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 118000001-121000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 119000001-122000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 120000001-123000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 121000001-124000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 122000001-125000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 123000001-126000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 124000001-127000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 125000001-128000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 126000001-129000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 127000001-130000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 128000001-131000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 129000001-132000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 130000001-133000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 131000001-134000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 132000001-135000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 133000001-136000001 doesn't exist
    [WARNING]  output file for chromosome 11 bp 134000001-137000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 1000001-4000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 8000001-11000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 9000001-12000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 91000001-94000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 92000001-95000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 93000001-96000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 94000001-97000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 95000001-98000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 96000001-99000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 97000001-100000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 98000001-101000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 99000001-102000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 100000001-103000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 101000001-104000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 102000001-105000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 103000001-106000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 104000001-107000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 105000001-108000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 106000001-109000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 107000001-110000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 108000001-111000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 109000001-112000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 110000001-113000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 111000001-114000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 112000001-115000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 113000001-116000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 114000001-117000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 115000001-118000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 116000001-119000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 117000001-120000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 118000001-121000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 119000001-122000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 120000001-123000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 121000001-124000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 122000001-125000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 123000001-126000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 124000001-127000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 125000001-128000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 126000001-129000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 127000001-130000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 128000001-131000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 129000001-132000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 130000001-133000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 131000001-134000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 132000001-135000001 doesn't exist
    [WARNING]  output file for chromosome 12 bp 133000001-136000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 91000001-94000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 92000001-95000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 93000001-96000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 94000001-97000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 95000001-98000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 96000001-99000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 97000001-100000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 98000001-101000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 99000001-102000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 100000001-103000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 101000001-104000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 102000001-105000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 103000001-106000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 104000001-107000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 105000001-108000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 106000001-109000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 107000001-110000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 108000001-111000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 109000001-112000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 110000001-113000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 111000001-114000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 112000001-115000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 113000001-116000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 114000001-117000001 doesn't exist
    [WARNING]  output file for chromosome 13 bp 115000001-118000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 91000001-94000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 92000001-95000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 93000001-96000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 94000001-97000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 95000001-98000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 96000001-99000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 97000001-100000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 98000001-101000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 99000001-102000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 100000001-103000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 101000001-104000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 102000001-105000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 103000001-106000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 104000001-107000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 105000001-108000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 106000001-109000001 doesn't exist
    [WARNING]  output file for chromosome 14 bp 107000001-110000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 91000001-94000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 92000001-95000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 93000001-96000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 94000001-97000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 95000001-98000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 96000001-99000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 97000001-100000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 98000001-101000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 99000001-102000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 100000001-103000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 101000001-104000001 doesn't exist
    [WARNING]  output file for chromosome 15 bp 102000001-105000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 100       FID      IID  Sex
    0  HG00097  HG00097    2
    1  HG00099  HG00099    2
    2  HG00101  HG00101    1
    3  HG00102  HG00102    2
    4  HG00103  HG00103    1
    380
    380
             FID      IID  Sex       PC1       PC2       PC3       PC4       PC5  \
    0    HG00097  HG00097    2 -0.000939  0.083605  0.007354 -0.010632  0.017002   
    1    HG00099  HG00099    2 -0.002231  0.088623 -0.016674 -0.000493 -0.012635   
    2    HG00101  HG00101    1 -0.000295  0.095897 -0.016722  0.019230 -0.000330   
    3    HG00102  HG00102    2  0.000205  0.073966  0.015507  0.046601  0.005690   
    4    HG00103  HG00103    1 -0.008589  0.064142 -0.007911  0.027496 -0.001118   
    ..       ...      ...  ...       ...       ...       ...       ...       ...   
    375  NA20818  NA20818    2 -0.047710 -0.040114 -0.053486 -0.017914  0.009163   
    376  NA20826  NA20826    2 -0.042617 -0.059276 -0.067061 -0.015577 -0.007602   
    377  NA20827  NA20827    1 -0.043531 -0.052465 -0.066099 -0.013848  0.004271   
    378  NA20828  NA20828    2 -0.048143 -0.049573 -0.040417 -0.016462 -0.021911   
    379  NA20832  NA20832    2 -0.041337 -0.047921 -0.045979  0.001853  0.004596   
    
              PC6  
    0   -0.037442  
    1   -0.007065  
    2   -0.029482  
    3   -0.020718  
    4    0.023350  
    ..        ...  
    375  0.055850  
    376 -0.003985  
    377  0.000613  
    378 -0.009084  
    379  0.026654  
    
    [380 rows x 9 columns]
    ./bolt --bfile=SampleData1/Fold_0/train_data.QC.clumped.pruned --phenoFile=SampleData1/Fold_0/train_data.PHENO --phenoCol=PHENO --lmmForceNonInf --LDscoresFile=tables/LDSCORE.1000G_EUR.tab.gz --covarFile=SampleData1/Fold_0/train_data.COV_PCA --LDscoresMatchBp --qCovarCol=COV_{1:5} --statsFile=SampleData1/Fold_0/SampleData1.lmmForceNonInf_stat
    0001-4000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 7000001-10000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 8000001-11000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 9000001-12000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 82000001-85000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 83000001-86000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 84000001-87000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 85000001-88000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 86000001-89000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 87000001-90000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 88000001-91000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 89000001-92000001 doesn't exist
    [WARNING]  output file for chromosome 16 bp 90000001-93000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 1000001-4000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 7000001-10000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 8000001-11000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 9000001-12000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 79000001-82000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 80000001-83000001 doesn't exist
    [WARNING]  output file for chromosome 17 bp 81000001-84000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 1000001-4000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 7000001-10000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 8000001-11000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 9000001-12000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 63000001-66000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 64000001-67000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 65000001-68000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 66000001-69000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 67000001-70000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 68000001-71000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 69000001-72000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 70000001-73000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 71000001-74000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 72000001-75000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 73000001-76000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 74000001-77000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 75000001-78000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 76000001-79000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 77000001-80000001 doesn't exist
    [WARNING]  output file for chromosome 18 bp 78000001-81000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 1000001-4000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 7000001-10000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 8000001-11000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 9000001-12000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 19 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 1-3000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 1000001-4000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 2000001-5000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 3000001-6000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 4000001-7000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 5000001-8000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 6000001-9000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 7000001-10000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 8000001-11000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 9000001-12000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 10000001-13000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 11000001-14000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 51000001-54000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 52000001-55000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 53000001-56000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 54000001-57000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 55000001-58000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 56000001-59000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 57000001-60000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 58000001-61000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 59000001-62000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 60000001-63000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 61000001-64000001 doesn't exist
    [WARNING]  output file for chromosome 20 bp 62000001-65000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 12000001-15000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 13000001-16000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 21 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 14000001-17000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 15000001-18000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 16000001-19000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 17000001-20000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 18000001-21000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 19000001-22000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 20000001-23000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 21000001-24000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 22000001-25000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 23000001-26000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 24000001-27000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 25000001-28000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 26000001-29000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 27000001-30000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 28000001-31000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 29000001-32000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 30000001-33000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 31000001-34000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 32000001-35000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 33000001-36000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 34000001-37000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 35000001-38000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 36000001-39000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 37000001-40000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 38000001-41000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 39000001-42000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 40000001-43000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 41000001-44000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 42000001-45000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 43000001-46000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 44000001-47000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 45000001-48000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 46000001-49000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 47000001-50000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 48000001-51000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 49000001-52000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 50000001-53000001 doesn't exist
    [WARNING]  output file for chromosome 22 bp 51000001-54000001 doesn't exist
    [INFO]  Wrote aggregated results to SampleData1/Fold_0/SampleData1.finalpolypred
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
        --LDscoresMatchBp \
        --qCovarCol=COV_{1:5} \
        --statsFile=SampleData1/Fold_0/SampleData1.lmmForceNonInf_stat 
    
    Setting number of threads to 1
    fam: SampleData1/Fold_0/train_data.QC.clumped.pruned.fam
    bim(s): SampleData1/Fold_0/train_data.QC.clumped.pruned.bim
    bed(s): SampleData1/Fold_0/train_data.QC.clumped.pruned.bed
    
    === Reading genotype data ===
    
    Total indivs in PLINK data: Nbed = 380
    Total indivs stored in memory: N = 380
    Reading bim file #1: SampleData1/Fold_0/train_data.QC.clumped.pruned.bim
    

.. parsed-literal::

    2743it [00:01, 1929.71it/s]
    

.. parsed-literal::

        Read 172070 snps
    Total snps in PLINK data: Mbed = 172070
    
    Breakdown of SNP pre-filtering results:
      172070 SNPs to include in model (i.e., GRM)
      0 additional non-GRM SNPs loaded
      0 excluded SNPs
    Allocating 172070 x 380/4 bytes to store genotypes
    Reading genotypes and performing QC filtering on snps and indivs...
    Reading bed file #1: SampleData1/Fold_0/train_data.QC.clumped.pruned.bed
        Expecting 16346650 (+3) bytes for 380 indivs, 172070 snps
    

.. parsed-literal::

    WARNING: Genetic map appears to be in cM units; rescaling by 0.01
    

.. parsed-literal::

    Total indivs after QC: 380
    Total post-QC SNPs: M = 172070
      Variance component 1: 172070 post-QC SNPs (name: 'modelSnps')
    Time for SnpData setup = 0.738851 sec
    
    === Reading phenotype and covariate data ===
    
    Read data for 380 indivs (ignored 0 without genotypes) from:
      SampleData1/Fold_0/train_data.COV_PCA
    Read data for 380 indivs (ignored 0 without genotypes) from:
      SampleData1/Fold_0/train_data.PHENO
    Number of indivs with no missing phenotype(s) to use: 380
        Using quantitative covariate: COV_1
        Using quantitative covariate: COV_2
        Using quantitative covariate: COV_3
        Using quantitative covariate: COV_4
        Using quantitative covariate: COV_5
        Using quantitative covariate: CONST_ALL_ONES
    Number of individuals used in analysis: Nused = 380
    Singular values of covariate matrix:
        S[0] = 36.3836
        S[1] = 5.21992
        S[2] = 1
        S[3] = 1
        S[4] = 1
        S[5] = 0.992648
    Total covariate vectors: C = 6
    Total independent covariate vectors: Cindep = 6
    
    === Initializing Bolt object: projecting and normalizing SNPs ===
    
    

.. parsed-literal::

    NOTE: Using all-1s vector (constant term) in addition to specified covariates
    

.. parsed-literal::

    Number of chroms with >= 1 good SNP: 22
    Average norm of projected SNPs:           371.280532
    Dimension of all-1s proj space (Nused-1): 379
    Time for covariate data setup + Bolt initialization = 0.305559 sec
    
    Phenotype 1:   N = 380   mean = 170.14   std = 0.945865
    
    === Computing linear regression (LINREG) stats ===
    
    Time for computing LINREG stats = 0.0798419 sec
    
    === Estimating variance parameters ===
    
    Using CGtol of 0.005 for this step
    Using default number of random trials: 15 (for Nused = 380)
    
    Estimating MC scaling f_REML at log(delta) = 1.09131, h2 = 0.25...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.25  rNorms/orig: (0.02,0.02)  res2s: 868.646..147.756
      iter 2:  time=0.25  rNorms/orig: (0.0004,0.0005)  res2s: 869.426..147.897
      Converged at iter 2: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 31.7%, memory/overhead = 68.3%
      MCscaling: logDelta = 1.09, h2 = 0.250, f = -0.00114392
    
    Estimating MC scaling f_REML at log(delta) = 1.93861, h2 = 0.125...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.25  rNorms/orig: (0.009,0.01)  res2s: 2395.32..201.169
      iter 2:  time=0.25  rNorms/orig: (0.0001,0.0001)  res2s: 2395.87..201.217
      Converged at iter 2: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 31.0%, memory/overhead = 69.0%
      MCscaling: logDelta = 1.94, h2 = 0.125, f = -0.000187244
    
    Estimating MC scaling f_REML at log(delta) = 2.10445, h2 = 0.10796...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.25  rNorms/orig: (0.007,0.008)  res2s: 2888.69..209.088
      iter 2:  time=0.25  rNorms/orig: (7e-05,8e-05)  res2s: 2889.19..209.126
      Converged at iter 2: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 31.0%, memory/overhead = 69.0%
      MCscaling: logDelta = 2.10, h2 = 0.108, f = -3.94063e-05
    
    Estimating MC scaling f_REML at log(delta) = 2.14865, h2 = 0.103777...
      Batch-solving 16 systems of equations using conjugate gradient iteration
      iter 1:  time=0.25  rNorms/orig: (0.007,0.008)  res2s: 3035.09..211.056
      iter 2:  time=0.25  rNorms/orig: (7e-05,8e-05)  res2s: 3035.57..211.091
      Converged at iter 2: rNorms/orig all < CGtol=0.005
      Time breakdown: dgemm = 31.0%, memory/overhead = 69.0%
      MCscaling: logDelta = 2.15, h2 = 0.104, f = -2.15314e-06
    
    Secant iteration for h2 estimation converged in 2 steps
    Estimated (pseudo-)heritability: h2g = 0.104
    To more precisely estimate variance parameters and estimate s.e., use --reml
    Variance params: sigma^2_K = 0.073458, logDelta = 2.148652, f = -2.15314e-06
    
    Time for fitting variance components = 2.65798 sec
    
    === Computing mixed model assoc stats (inf. model) ===
    
    Selected 30 SNPs for computation of prospective stat
    Tried 31; threw out 1 with GRAMMAR chisq > 5
    Assigning SNPs to 22 chunks for leave-out analysis
    Each chunk is excluded when testing SNPs belonging to the chunk
      Batch-solving 52 systems of equations using conjugate gradient iteration
      iter 1:  time=0.44  rNorms/orig: (0.006,0.009)  res2s: 211.028..301.542
      iter 2:  time=0.44  rNorms/orig: (5e-05,9e-05)  res2s: 211.066..301.603
      Converged at iter 2: rNorms/orig all < CGtol=0.0005
      Time breakdown: dgemm = 60.1%, memory/overhead = 39.9%
    
    AvgPro: 0.955   AvgRetro: 0.955   Calibration: 1.000 (0.000)   (30 SNPs)
    Ratio of medians: 1.000   Median of ratios: 1.000
    
    Time for computing infinitesimal model assoc stats = 1.00979 sec
    
    === Estimating chip LD Scores using 400 indivs ===
    
    Reducing sample size to 376 for memory alignment
    

.. parsed-literal::

    WARNING: Only 380 indivs available; using all
    

.. parsed-literal::

    
    Time for estimating chip LD Scores = 0.471238 sec
    
    === Reading LD Scores for calibration of Bayesian assoc stats ===
    
    Looking up LD Scores...
      Looking for column header 'CHR': column number = 2
      Looking for column header 'BP': column number = 3
      Looking for column header 'LDSCORE': column number = 5
    Found LD Scores for 171152/172070 SNPs
    
    Estimating inflation of LINREG chisq stats using MLMe as reference...
    Filtering to SNPs with chisq stats, LD Scores, and MAF > 0.01
    # of SNPs passing filters before outlier removal: 171152/172070
    Masking windows around outlier snps (chisq > 20.0)
    # of SNPs remaining after outlier window removal: 171152/171152
    Intercept of LD Score regression for ref stats:   0.998 (0.005)
    Estimated attenuation: -4.747 (83.057)
    Intercept of LD Score regression for cur stats: 0.999 (0.005)
    Calibration factor (ref/cur) to multiply by:      0.999 (0.000)
    LINREG intercept inflation = 1.0011
    
    === Estimating mixture parameters by cross-validation ===
    
    Setting maximum number of iterations to 250 for this step
    Max CV folds to compute = 5 (to have > 10000 samples)
    
    ====> Starting CV fold 1 <====
    
        Using quantitative covariate: COV_1
        Using quantitative covariate: COV_2
        Using quantitative covariate: COV_3
        Using quantitative covariate: COV_4
        Using quantitative covariate: COV_5
        Using quantitative covariate: CONST_ALL_ONES
    Number of individuals used in analysis: Nused = 304
    Singular values of covariate matrix:
        S[0] = 32.5614
        S[1] = 4.66584
        S[2] = 0.964683
        S[3] = 0.911773
        S[4] = 0.894565
        S[5] = 0.887678
    Total covariate vectors: C = 6
    Total independent covariate vectors: Cindep = 6
    
    === Initializing Bolt object: projecting and normalizing SNPs ===
    
    

.. parsed-literal::

    NOTE: Using all-1s vector (constant term) in addition to specified covariates
    

.. parsed-literal::

    Number of chroms with >= 1 good SNP: 22
    Average norm of projected SNPs:           295.802788
    Dimension of all-1s proj space (Nused-1): 303
      Beginning variational Bayes
      iter 1:  time=0.65 for 18 active reps
      iter 2:  time=0.47 for 18 active reps  approxLL diffs: (0.07,0.08)
      iter 3:  time=0.46 for 18 active reps  approxLL diffs: (0.00,0.00)
      Converged at iter 3: approxLL diffs each have been < LLtol=0.01
      Time breakdown: dgemm = 30.3%, memory/overhead = 69.7%
    Computing predictions on left-out cross-validation fold
    Time for computing predictions = 0.309486 sec
    
    Average PVEs obtained by param pairs tested (high to low):
     f2=0.1, p=0.01: 0.000265
     f2=0.3, p=0.01: 0.000260
     f2=0.1, p=0.02: 0.000258
                ...
      f2=0.5, p=0.5: 0.000252
    
    Detailed CV fold results:
      Absolute prediction MSE baseline (covariates only): 0.57355
      Absolute prediction MSE using standard LMM:         0.573405
      Absolute prediction MSE, fold-best f2=0.1, p=0.01:  0.573398
        Absolute pred MSE using   f2=0.5, p=0.5: 0.573405
        Absolute pred MSE using   f2=0.5, p=0.2: 0.573405
        Absolute pred MSE using   f2=0.5, p=0.1: 0.573405
        Absolute pred MSE using  f2=0.5, p=0.05: 0.573405
        Absolute pred MSE using  f2=0.5, p=0.02: 0.573404
        Absolute pred MSE using  f2=0.5, p=0.01: 0.573403
        Absolute pred MSE using   f2=0.3, p=0.5: 0.573405
        Absolute pred MSE using   f2=0.3, p=0.2: 0.573405
        Absolute pred MSE using   f2=0.3, p=0.1: 0.573405
        Absolute pred MSE using  f2=0.3, p=0.05: 0.573404
        Absolute pred MSE using  f2=0.3, p=0.02: 0.573403
        Absolute pred MSE using  f2=0.3, p=0.01: 0.573401
        Absolute pred MSE using   f2=0.1, p=0.5: 0.573405
        Absolute pred MSE using   f2=0.1, p=0.2: 0.573405
        Absolute pred MSE using   f2=0.1, p=0.1: 0.573405
        Absolute pred MSE using  f2=0.1, p=0.05: 0.573404
        Absolute pred MSE using  f2=0.1, p=0.02: 0.573402
        Absolute pred MSE using  f2=0.1, p=0.01: 0.573398
    
    ====> End CV fold 1: 18 remaining param pair(s) <====
    
    Estimated proportion of variance explained using inf model: 0.000
    Relative improvement in prediction MSE using non-inf model: 0.000
    
    ====> Starting CV fold 2 <====
    
        Using quantitative covariate: COV_1
        Using quantitative covariate: COV_2
        Using quantitative covariate: COV_3
        Using quantitative covariate: COV_4
        Using quantitative covariate: COV_5
        Using quantitative covariate: CONST_ALL_ONES
    Number of individuals used in analysis: Nused = 304
    Singular values of covariate matrix:
        S[0] = 32.3248
        S[1] = 4.70352
        S[2] = 0.946981
        S[3] = 0.913355
        S[4] = 0.885557
        S[5] = 0.870892
    Total covariate vectors: C = 6
    Total independent covariate vectors: Cindep = 6
    
    === Initializing Bolt object: projecting and normalizing SNPs ===
    
    

.. parsed-literal::

    NOTE: Using all-1s vector (constant term) in addition to specified covariates
    

.. parsed-literal::

    Number of chroms with >= 1 good SNP: 22
    Average norm of projected SNPs:           295.830429
    Dimension of all-1s proj space (Nused-1): 303
      Beginning variational Bayes
      iter 1:  time=0.64 for 18 active reps
      iter 2:  time=0.46 for 18 active reps  approxLL diffs: (0.07,0.07)
      iter 3:  time=0.46 for 18 active reps  approxLL diffs: (0.00,0.00)
      Converged at iter 3: approxLL diffs each have been < LLtol=0.01
      Time breakdown: dgemm = 30.2%, memory/overhead = 69.8%
    Computing predictions on left-out cross-validation fold
    Time for computing predictions = 0.309555 sec
    
    Average PVEs obtained by param pairs tested (high to low):
      f2=0.5, p=0.5: -0.000038
      f2=0.3, p=0.5: -0.000038
      f2=0.5, p=0.2: -0.000038
                ...
     f2=0.1, p=0.01: -0.000043
    
    Detailed CV fold results:
      Absolute prediction MSE baseline (covariates only): 0.873917
      Absolute prediction MSE using standard LMM:         0.874204
      Absolute prediction MSE, fold-best  f2=0.5, p=0.5:  0.874204
        Absolute pred MSE using   f2=0.5, p=0.5: 0.874204
        Absolute pred MSE using   f2=0.5, p=0.2: 0.874204
        Absolute pred MSE using   f2=0.5, p=0.1: 0.874205
        Absolute pred MSE using  f2=0.5, p=0.05: 0.874205
        Absolute pred MSE using  f2=0.5, p=0.02: 0.874207
        Absolute pred MSE using  f2=0.5, p=0.01: 0.874210
        Absolute pred MSE using   f2=0.3, p=0.5: 0.874204
        Absolute pred MSE using   f2=0.3, p=0.2: 0.874205
        Absolute pred MSE using   f2=0.3, p=0.1: 0.874205
        Absolute pred MSE using  f2=0.3, p=0.05: 0.874206
        Absolute pred MSE using  f2=0.3, p=0.02: 0.874210
        Absolute pred MSE using  f2=0.3, p=0.01: 0.874216
        Absolute pred MSE using   f2=0.1, p=0.5: 0.874204
        Absolute pred MSE using   f2=0.1, p=0.2: 0.874205
        Absolute pred MSE using   f2=0.1, p=0.1: 0.874206
        Absolute pred MSE using  f2=0.1, p=0.05: 0.874208
        Absolute pred MSE using  f2=0.1, p=0.02: 0.874214
        Absolute pred MSE using  f2=0.1, p=0.01: 0.874224
    
    ====> End CV fold 2: 18 remaining param pair(s) <====
    
    ====> Starting CV fold 3 <====
    
        Using quantitative covariate: COV_1
        Using quantitative covariate: COV_2
        Using quantitative covariate: COV_3
        Using quantitative covariate: COV_4
        Using quantitative covariate: COV_5
        Using quantitative covariate: CONST_ALL_ONES
    Number of individuals used in analysis: Nused = 304
    Singular values of covariate matrix:
        S[0] = 32.9844
        S[1] = 4.58824
        S[2] = 0.907683
        S[3] = 0.883725
        S[4] = 0.872965
        S[5] = 0.627063
    Total covariate vectors: C = 6
    Total independent covariate vectors: Cindep = 6
    
    === Initializing Bolt object: projecting and normalizing SNPs ===
    
    

.. parsed-literal::

    NOTE: Using all-1s vector (constant term) in addition to specified covariates
    

.. parsed-literal::

    Number of chroms with >= 1 good SNP: 22
    Average norm of projected SNPs:           295.843895
    Dimension of all-1s proj space (Nused-1): 303
      Beginning variational Bayes
      iter 1:  time=0.64 for 18 active reps
      iter 2:  time=0.47 for 18 active reps  approxLL diffs: (0.07,0.07)
      iter 3:  time=0.47 for 18 active reps  approxLL diffs: (0.00,0.00)
      Converged at iter 3: approxLL diffs each have been < LLtol=0.01
      Time breakdown: dgemm = 30.2%, memory/overhead = 69.8%
    Computing predictions on left-out cross-validation fold
    Time for computing predictions = 0.343775 sec
    
    Average PVEs obtained by param pairs tested (high to low):
      f2=0.5, p=0.5: 0.000444
      f2=0.3, p=0.5: 0.000444
      f2=0.5, p=0.2: 0.000444
                ...
     f2=0.1, p=0.01: 0.000435
    
    Detailed CV fold results:
      Absolute prediction MSE baseline (covariates only): 0.830538
      Absolute prediction MSE using standard LMM:         0.829367
      Absolute prediction MSE, fold-best  f2=0.5, p=0.5:  0.829367
        Absolute pred MSE using   f2=0.5, p=0.5: 0.829367
        Absolute pred MSE using   f2=0.5, p=0.2: 0.829367
        Absolute pred MSE using   f2=0.5, p=0.1: 0.829368
        Absolute pred MSE using  f2=0.5, p=0.05: 0.829368
        Absolute pred MSE using  f2=0.5, p=0.02: 0.829369
        Absolute pred MSE using  f2=0.5, p=0.01: 0.829372
        Absolute pred MSE using   f2=0.3, p=0.5: 0.829367
        Absolute pred MSE using   f2=0.3, p=0.2: 0.829367
        Absolute pred MSE using   f2=0.3, p=0.1: 0.829368
        Absolute pred MSE using  f2=0.3, p=0.05: 0.829369
        Absolute pred MSE using  f2=0.3, p=0.02: 0.829371
        Absolute pred MSE using  f2=0.3, p=0.01: 0.829376
        Absolute pred MSE using   f2=0.1, p=0.5: 0.829367
        Absolute pred MSE using   f2=0.1, p=0.2: 0.829368
        Absolute pred MSE using   f2=0.1, p=0.1: 0.829368
        Absolute pred MSE using  f2=0.1, p=0.05: 0.829370
        Absolute pred MSE using  f2=0.1, p=0.02: 0.829374
        Absolute pred MSE using  f2=0.1, p=0.01: 0.829381
    
    ====> End CV fold 3: 18 remaining param pair(s) <====
    
    ====> Starting CV fold 4 <====
    
        Using quantitative covariate: COV_1
        Using quantitative covariate: COV_2
        Using quantitative covariate: COV_3
        Using quantitative covariate: COV_4
        Using quantitative covariate: COV_5
        Using quantitative covariate: CONST_ALL_ONES
    Number of individuals used in analysis: Nused = 304
    Singular values of covariate matrix:
        S[0] = 32.5614
        S[1] = 4.66525
        S[2] = 0.965453
        S[3] = 0.896721
        S[4] = 0.89062
        S[5] = 0.858094
    Total covariate vectors: C = 6
    Total independent covariate vectors: Cindep = 6
    
    === Initializing Bolt object: projecting and normalizing SNPs ===
    
    

.. parsed-literal::

    NOTE: Using all-1s vector (constant term) in addition to specified covariates
    

.. parsed-literal::

    Number of chroms with >= 1 good SNP: 22
    Average norm of projected SNPs:           295.834790
    Dimension of all-1s proj space (Nused-1): 303
      Beginning variational Bayes
      iter 1:  time=0.64 for 18 active reps
      iter 2:  time=0.47 for 18 active reps  approxLL diffs: (0.07,0.07)
      iter 3:  time=0.47 for 18 active reps  approxLL diffs: (0.00,0.00)
      Converged at iter 3: approxLL diffs each have been < LLtol=0.01
      Time breakdown: dgemm = 30.2%, memory/overhead = 69.8%
    Computing predictions on left-out cross-validation fold
    Time for computing predictions = 0.309342 sec
    
    Average PVEs obtained by param pairs tested (high to low):
      f2=0.5, p=0.5: 0.000576
      f2=0.3, p=0.5: 0.000576
      f2=0.5, p=0.2: 0.000576
                ...
     f2=0.1, p=0.01: 0.000561
    
    Detailed CV fold results:
      Absolute prediction MSE baseline (covariates only): 0.74444
      Absolute prediction MSE using standard LMM:         0.743717
      Absolute prediction MSE, fold-best  f2=0.5, p=0.5:  0.743717
        Absolute pred MSE using   f2=0.5, p=0.5: 0.743717
        Absolute pred MSE using   f2=0.5, p=0.2: 0.743717
        Absolute pred MSE using   f2=0.5, p=0.1: 0.743717
        Absolute pred MSE using  f2=0.5, p=0.05: 0.743718
        Absolute pred MSE using  f2=0.5, p=0.02: 0.743720
        Absolute pred MSE using  f2=0.5, p=0.01: 0.743724
        Absolute pred MSE using   f2=0.3, p=0.5: 0.743717
        Absolute pred MSE using   f2=0.3, p=0.2: 0.743717
        Absolute pred MSE using   f2=0.3, p=0.1: 0.743718
        Absolute pred MSE using  f2=0.3, p=0.05: 0.743720
        Absolute pred MSE using  f2=0.3, p=0.02: 0.743724
        Absolute pred MSE using  f2=0.3, p=0.01: 0.743732
        Absolute pred MSE using   f2=0.1, p=0.5: 0.743717
        Absolute pred MSE using   f2=0.1, p=0.2: 0.743718
        Absolute pred MSE using   f2=0.1, p=0.1: 0.743719
        Absolute pred MSE using  f2=0.1, p=0.05: 0.743722
        Absolute pred MSE using  f2=0.1, p=0.02: 0.743729
        Absolute pred MSE using  f2=0.1, p=0.01: 0.743742
    
    ====> End CV fold 4: 18 remaining param pair(s) <====
    
    ====> Starting CV fold 5 <====
    
        Using quantitative covariate: COV_1
        Using quantitative covariate: COV_2
        Using quantitative covariate: COV_3
        Using quantitative covariate: COV_4
        Using quantitative covariate: COV_5
        Using quantitative covariate: CONST_ALL_ONES
    Number of individuals used in analysis: Nused = 304
    Singular values of covariate matrix:
        S[0] = 32.2773
        S[1] = 4.71004
        S[2] = 0.929851
        S[3] = 0.919533
        S[4] = 0.893736
        S[5] = 0.87937
    Total covariate vectors: C = 6
    Total independent covariate vectors: Cindep = 6
    
    === Initializing Bolt object: projecting and normalizing SNPs ===
    
    

.. parsed-literal::

    NOTE: Using all-1s vector (constant term) in addition to specified covariates
    

.. parsed-literal::

    Number of chroms with >= 1 good SNP: 22
    Average norm of projected SNPs:           295.802180
    Dimension of all-1s proj space (Nused-1): 303
      Beginning variational Bayes
      iter 1:  time=0.65 for 18 active reps
      iter 2:  time=0.47 for 18 active reps  approxLL diffs: (0.07,0.08)
      iter 3:  time=0.48 for 18 active reps  approxLL diffs: (0.00,0.00)
      Converged at iter 3: approxLL diffs each have been < LLtol=0.01
      Time breakdown: dgemm = 30.0%, memory/overhead = 70.0%
    Computing predictions on left-out cross-validation fold
    Time for computing predictions = 0.309355 sec
    
    Average PVEs obtained by param pairs tested (high to low):
      f2=0.5, p=0.5: 0.000587
      f2=0.3, p=0.5: 0.000587
      f2=0.5, p=0.2: 0.000587
                ...
     f2=0.1, p=0.01: 0.000571
    
    Detailed CV fold results:
      Absolute prediction MSE baseline (covariates only): 0.548952
      Absolute prediction MSE using standard LMM:         0.548606
      Absolute prediction MSE, fold-best  f2=0.5, p=0.5:  0.548606
        Absolute pred MSE using   f2=0.5, p=0.5: 0.548606
        Absolute pred MSE using   f2=0.5, p=0.2: 0.548606
        Absolute pred MSE using   f2=0.5, p=0.1: 0.548606
        Absolute pred MSE using  f2=0.5, p=0.05: 0.548607
        Absolute pred MSE using  f2=0.5, p=0.02: 0.548607
        Absolute pred MSE using  f2=0.5, p=0.01: 0.548609
        Absolute pred MSE using   f2=0.3, p=0.5: 0.548606
        Absolute pred MSE using   f2=0.3, p=0.2: 0.548606
        Absolute pred MSE using   f2=0.3, p=0.1: 0.548607
        Absolute pred MSE using  f2=0.3, p=0.05: 0.548607
        Absolute pred MSE using  f2=0.3, p=0.02: 0.548609
        Absolute pred MSE using  f2=0.3, p=0.01: 0.548612
        Absolute pred MSE using   f2=0.1, p=0.5: 0.548606
        Absolute pred MSE using   f2=0.1, p=0.2: 0.548606
        Absolute pred MSE using   f2=0.1, p=0.1: 0.548607
        Absolute pred MSE using  f2=0.1, p=0.05: 0.548608
        Absolute pred MSE using  f2=0.1, p=0.02: 0.548611
        Absolute pred MSE using  f2=0.1, p=0.01: 0.548616
    
    ====> End CV fold 5: 18 remaining param pair(s) <====
    
    Optimal mixture parameters according to CV: f2 = 0.5, p = 0.5
    
    Time for estimating mixture parameters = 25.2059 sec
    
    === Computing Bayesian mixed model assoc stats with mixture prior ===
    
    Assigning SNPs to 22 chunks for leave-out analysis
    Each chunk is excluded when testing SNPs belonging to the chunk
      Beginning variational Bayes
      iter 1:  time=0.67 for 22 active reps
      iter 2:  time=0.50 for 22 active reps  approxLL diffs: (0.09,0.09)
      iter 3:  time=0.50 for 22 active reps  approxLL diffs: (0.00,0.00)
      Converged at iter 3: approxLL diffs each have been < LLtol=0.01
      Time breakdown: dgemm = 29.3%, memory/overhead = 70.7%
    Filtering to SNPs with chisq stats, LD Scores, and MAF > 0.01
    # of SNPs passing filters before outlier removal: 171152/172070
    Masking windows around outlier snps (chisq > 20.0)
    # of SNPs remaining after outlier window removal: 171152/171152
    Intercept of LD Score regression for ref stats:   0.998 (0.005)
    Estimated attenuation: -4.747 (83.057)
    Intercept of LD Score regression for cur stats: 0.998 (0.005)
    Calibration factor (ref/cur) to multiply by:      1.000 (0.000)
    
    Time for computing Bayesian mixed model assoc stats = 2.0642 sec
    
    Calibration stats: mean and lambdaGC (over SNPs used in GRM)
      (note that both should be >1 because of polygenicity)
    Mean BOLT_LMM_INF: 1.0003 (172070 good SNPs)   lambdaGC: 1.00764
    Mean BOLT_LMM: 1.00031 (172070 good SNPs)   lambdaGC: 1.00762
    
    === Streaming genotypes to compute and write assoc stats at all SNPs ===
    
    Time for streaming genotypes and writing output = 4.04339 sec
    
    Total elapsed time for analysis = 36.5768 sec
    *********************************************************************
    * PolyPred (POLYgenic risk PREDiction)
    * Version 1.0.0
    * (C) 2020-2024 Omer Weissbrod
    *********************************************************************
    
    [INFO]  Loading betas file SampleData1/Fold_0/SampleData1.lmmForceNonInf_stat...
    [INFO]  done in 0.26 seconds
    

.. parsed-literal::

    SampleData1/Fold_0/train_data.QC.clumped.pruned.bed:   0%|          | 0/1 [00:00<?, ?it/s]

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to /tmp/tmp8hzxhfr7/3kdpkmzd.log.
    Options in effect:
      --allow-no-sex
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --extract /tmp/tmp8hzxhfr7/hq25porf
      --keep SampleData1/Fold_0/train_data.PHENO
      --memory 2048
      --out /tmp/tmp8hzxhfr7/3kdpkmzd
      --score /tmp/tmp8hzxhfr7/hq25porf sum
      --threads 1
    
    63761 MB RAM detected; reserving 2048 MB for main workspace.
    172070 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 172070 variants remaining.
    --keep: 380 people remaining.
    Using 1 thread.
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is exactly 1.
    172070 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172070 valid predictors loaded.
    --score: Results written to /tmp/tmp8hzxhfr7/3kdpkmzd.profile .
    [INFO]  Loading betas file SampleData1/Fold_0/SampleData1.finalpolypred...
    [INFO]  done in 0.02 seconds
    

.. parsed-literal::

    SampleData1/Fold_0/train_data.QC.clumped.pruned.bed: 100%|██████████| 1/1 [00:00<00:00,  1.03it/s]
    SampleData1/Fold_0/train_data.QC.clumped.pruned.bed:   0%|          | 0/1 [00:00<?, ?it/s]

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to /tmp/tmpwgtdo4fq/r8kz_iaw.log.
    Options in effect:
      --allow-no-sex
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --extract /tmp/tmpwgtdo4fq/x8nfu8ab
      --keep SampleData1/Fold_0/train_data.PHENO
      --memory 2048
      --out /tmp/tmpwgtdo4fq/r8kz_iaw
      --score /tmp/tmpwgtdo4fq/x8nfu8ab sum
      --threads 1
    
    63761 MB RAM detected; reserving 2048 MB for main workspace.
    172070 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 3441 variants remaining.
    --keep: 380 people remaining.
    Using 1 thread.
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is exactly 1.
    3441 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 3441 valid predictors loaded.
    --score: Results written to /tmp/tmpwgtdo4fq/r8kz_iaw.profile .
    [INFO]  In-sample R2: 0.780
    [INFO]  Writing mixing weights to SampleData1/Fold_0/GWAS_polypred.mixweights
    
    

.. parsed-literal::

    SampleData1/Fold_0/train_data.QC.clumped.pruned.bed: 100%|██████████| 1/1 [00:00<00:00,  2.21it/s]
    

.. parsed-literal::

    *********************************************************************
    * PolyPred (POLYgenic risk PREDiction)
    * Version 1.0.0
    * (C) 2020-2024 Omer Weissbrod
    *********************************************************************
    
    [INFO]  Loading betas file SampleData1/Fold_0/SampleData1.lmmForceNonInf_stat...
    [INFO]  done in 0.23 seconds
    

.. parsed-literal::

    SampleData1/Fold_0/train_data.QC.clumped.pruned.bed:   0%|          | 0/1 [00:00<?, ?it/s]

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to /tmp/tmpcss56u2y/4edirvq6.log.
    Options in effect:
      --allow-no-sex
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --extract /tmp/tmpcss56u2y/i2yfjqob
      --memory 2048
      --out /tmp/tmpcss56u2y/4edirvq6
      --q-score-range /tmp/tmpcss56u2y/ranges.txt /tmp/tmpcss56u2y/cw3cappm
      --score /tmp/tmpcss56u2y/i2yfjqob sum
      --threads 1
    
    63761 MB RAM detected; reserving 2048 MB for main workspace.
    172070 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 172070 variants remaining.
    Using 1 thread.
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is exactly 1.
    172070 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172070 valid predictors loaded.
    --score: 200 ranges processed.
    Results written to /tmp/tmpcss56u2y/4edirvq6.*.profile.
    

.. parsed-literal::

    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    SampleData1/Fold_0/train_data.QC.clumped.pruned.bed: 100%|██████████| 1/1 [00:01<00:00,  1.78s/it]
    SampleData1/Fold_0/train_data.QC.clumped.pruned.bed:   0%|          | 0/1 [00:00<?, ?it/s]

.. parsed-literal::

    [INFO]  Loading betas file SampleData1/Fold_0/SampleData1.finalpolypred...
    [INFO]  done in 0.02 seconds
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to /tmp/tmpwh1tbrly/fri346er.log.
    Options in effect:
      --allow-no-sex
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --extract /tmp/tmpwh1tbrly/y1pykoe0
      --memory 2048
      --out /tmp/tmpwh1tbrly/fri346er
      --q-score-range /tmp/tmpwh1tbrly/ranges.txt /tmp/tmpwh1tbrly/xtmfed8h
      --score /tmp/tmpwh1tbrly/y1pykoe0 sum
      --threads 1
    
    63761 MB RAM detected; reserving 2048 MB for main workspace.
    172070 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 3441 variants remaining.
    Using 1 thread.
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is exactly 1.
    3441 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 3441 valid predictors loaded.
    --score: 200 ranges processed.
    Results written to /tmp/tmpwh1tbrly/fri346er.*.profile.
    

.. parsed-literal::

    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    SampleData1/Fold_0/train_data.QC.clumped.pruned.bed: 100%|██████████| 1/1 [00:01<00:00,  1.01s/it]
    

.. parsed-literal::

    [INFO]  Saving PRS to SampleData1/Fold_0/Train_PRS_polypred.prs
    
    *********************************************************************
    * PolyPred (POLYgenic risk PREDiction)
    * Version 1.0.0
    * (C) 2020-2024 Omer Weissbrod
    *********************************************************************
    
    [INFO]  Loading betas file SampleData1/Fold_0/SampleData1.lmmForceNonInf_stat...
    [INFO]  done in 0.23 seconds
    

.. parsed-literal::

    SampleData1/Fold_0/test_data.bed:   0%|          | 0/1 [00:00<?, ?it/s]

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to /tmp/tmp5gtaxcz0/siyzl_i2.log.
    Options in effect:
      --allow-no-sex
      --bfile SampleData1/Fold_0/test_data
      --extract /tmp/tmp5gtaxcz0/xe7gx0ls
      --memory 2048
      --out /tmp/tmp5gtaxcz0/siyzl_i2
      --q-score-range /tmp/tmp5gtaxcz0/ranges.txt /tmp/tmp5gtaxcz0/vby2lorq
      --score /tmp/tmp5gtaxcz0/xe7gx0ls sum
      --threads 1
    
    63761 MB RAM detected; reserving 2048 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    --extract: 172070 variants remaining.
    Using 1 thread.
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is exactly 1.
    172070 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 172070 valid predictors loaded.
    --score: 200 ranges processed.
    Results written to /tmp/tmp5gtaxcz0/siyzl_i2.*.profile.
    

.. parsed-literal::

    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    SampleData1/Fold_0/test_data.bed: 100%|██████████| 1/1 [00:02<00:00,  2.99s/it]
    SampleData1/Fold_0/test_data.bed:   0%|          | 0/1 [00:00<?, ?it/s]

.. parsed-literal::

    [INFO]  Loading betas file SampleData1/Fold_0/SampleData1.finalpolypred...
    [INFO]  done in 0.02 seconds
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to /tmp/tmpzoi6zf02/axg1dvbz.log.
    Options in effect:
      --allow-no-sex
      --bfile SampleData1/Fold_0/test_data
      --extract /tmp/tmpzoi6zf02/03gdqkf_
      --memory 2048
      --out /tmp/tmpzoi6zf02/axg1dvbz
      --q-score-range /tmp/tmpzoi6zf02/ranges.txt /tmp/tmpzoi6zf02/oymevenw
      --score /tmp/tmpzoi6zf02/03gdqkf_ sum
      --threads 1
    
    63761 MB RAM detected; reserving 2048 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    --extract: 9903 variants remaining.
    Using 1 thread.
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999897.
    9903 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --score: 9903 valid predictors loaded.
    --score: 200 ranges processed.
    Results written to /tmp/tmpzoi6zf02/axg1dvbz.*.profile.
    

.. parsed-literal::

    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    polyfun/polypred.py:153: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      df_prs[scoresum_colname] = df_jk[scoresum_colname]
    SampleData1/Fold_0/test_data.bed: 100%|██████████| 1/1 [00:01<00:00,  1.90s/it]
    

.. parsed-literal::

    [INFO]  Saving PRS to SampleData1/Fold_0/Test_PRS_polypred.prs
    
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
   python PolyPred.py 0
   python PolyPred.py 1
   python PolyPred.py 2
   python PolyPred.py 3
   python PolyPred.py 4

The following files should exist after the execution:

1. ``SampleData1/Fold_0/PolyPred/Results.csv``
2. ``SampleData1/Fold_1/PolyPred/Results.csv``
3. ``SampleData1/Fold_2/PolyPred/Results.csv``
4. ``SampleData1/Fold_3/PolyPred/Results.csv``
5. ``SampleData1/Fold_4/PolyPred/Results.csv``

Check the results file for each fold.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython2

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
    
    


.. parsed-literal::

    Fold_ 0 Yes, the file exists.
    Number of P-values processed:  1
    Fold_ 1 Yes, the file exists.
    Number of P-values processed:  1
    Fold_ 2 No, the file does not exist.
    Fold_ 3 Yes, the file exists.
    Number of P-values processed:  1
    Fold_ 4 Yes, the file exists.
    Number of P-values processed:  1
    

Sum the results for each fold.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython2

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
            'h2model',
            
            'model',
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
    Fold_ 1 Yes, the file exists.
    Fold_ 2 No, the file does not exist.
    Fold_ 3 Yes, the file exists.
    Fold_ 4 Yes, the file exists.
    Iteration 1:
    Unique rows in current common DataFrame: 1
    Unique rows in next DataFrame: 1
    Common rows after merge: 1
    
    Iteration 2:
    Unique rows in current common DataFrame: 1
    Unique rows in next DataFrame: 1
    Common rows after merge: 1
    
    Iteration 3:
    Unique rows in current common DataFrame: 1
    Unique rows in next DataFrame: 1
    Common rows after merge: 1
    
    DataFrame 1 with extracted common rows has 1 rows.
    DataFrame 2 with extracted common rows has 1 rows.
    DataFrame 3 with extracted common rows has 1 rows.
    DataFrame 4 with extracted common rows has 1 rows.
       clump_p1  clump_r2  clump_kb  p_window_size  p_slide_size  p_LD_threshold  \
    0       1.0       0.1     200.0          200.0          50.0            0.25   
    
       pvalue  numberofpca  tempalpha  l1weight  Train_pure_prs  Train_null_model  \
    0     1.0          6.0        0.1       0.1        0.775701          0.230336   
    
       Train_best_model  Test_pure_prs  Test_null_model  Test_best_model  
    0          0.997862      -0.010322         0.110899         0.125464  
    

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

.. code:: ipython2

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
    
    |                  |           0 |
    |:-----------------|------------:|
    | clump_p1         |   1         |
    | clump_r2         |   0.1       |
    | clump_kb         | 200         |
    | p_window_size    | 200         |
    | p_slide_size     |  50         |
    | p_LD_threshold   |   0.25      |
    | pvalue           |   1         |
    | numberofpca      |   6         |
    | tempalpha        |   0.1       |
    | l1weight         |   0.1       |
    | Train_pure_prs   |   0.775701  |
    | Train_null_model |   0.230336  |
    | Train_best_model |   0.997862  |
    | Test_pure_prs    |  -0.0103224 |
    | Test_null_model  |   0.110899  |
    | Test_best_model  |   0.125464  |
    


.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAABGUAAAKjCAYAAAC0i7LNAAAAAXNSR0IArs4c6QAAIABJREFUeF7snQm4TdX7x997jZlnujJXpoafWSQKkSQqU0SGQqJBkiFTiIiUDEmIpJQxMhMhREWUeZY5c9d07///Xezj3OPcc9Y595y999m+63l6qnvWXsPnXWfzfu/7visqPj4+XthIgARIgARIgARIgARIgARIgARIgARIgARMJRBFUcZU3pyMBEiABEiABEiABEiABEiABEiABEiABBQBijI8CCRAAiRAAiRAAiRAAiRAAiRAAiRAAiRgAQGKMhZA55QkQAIkQAIkQAIkQAIkQAIkQAIkQAIkQFGGZ4AESIAESIAESIAESIAESIAESIAESIAELCBAUcYC6JySBEiABEiABEiABEiABEiABEiABEiABCjK8AyQAAmQAAmQAAmQAAmQAAmQAAmQAAmQgAUEKMpYAJ1TkgAJkAAJkAAJkAAJkAAJkAAJkAAJkABFGZ4BEiABEiABEiABEiABEiABEiABEiABErCAAEUZC6BzShIgARIgARIgARIgARIgARIgARIgARKgKMMzQAIkQAIkQAIkQAIkQAIkQAIkQAIkQAIWEKAoYwF0TkkCJEACJEACJEACJEACJEACJEACJEACFGV4BkiABEiABEiABEiABEiABEiABEiABEjAAgIUZSyAzilJgARIgARIgARIgARIgARIgARIgARIgKIMzwAJkAAJkAAJkAAJkAAJkAAJkAAJkAAJWECAoowF0DklCZAACZAACZAACZAACZAACZAACZAACVCU4RkgARIgARIgARIgARIgARIgARIgARIgAQsIUJSxADqnJAESIAESIAESIAESIAESIAESIAESIAGKMjwDJEACJEACJEACJEACJEACJEACJEACJGABAYoyFkDnlCRAAiRAAiRAAiRAAiRAAiRAAiRAAiRAUYZngARIgARIgARIgARIgARIgARIgARIgAQsIEBRxgLonJIESIAESIAESIAESIAESIAESIAESIAEKMrwDJAACZAACZAACZAACZAACZAACZAACZCABQQoylgAnVOSAAmQAAmQAAmQAAmQAAmQAAmQAAmQAEUZngESIAESIAESIAESIAESIAESIAESIAESsIAARRkLoHNKEiABEiABEiABEiABEiABEiABEiABEqAowzNAAiRAAiRAAiRAAiRAAiRAAiRAAiRAAhYQoChjAXROSQIkQAIkQAIkQAIkQAIkQAIkQAIkQAIUZXgGSIAESIAESIAESIAESIAESIAESIAESMACAhRlLIDOKUmABEiABEiABEiABEiABEiABEiABEiAogzPAAmQAAmQAAmQAAmQAAmQAAmQAAmQAAlYQICijAXQOSUJkAAJkAAJkAAJkAAJkAAJkAAJkAAJUJThGSABEiABEiABEiABEiABEiABEiABEiABCwhQlLEAOqckARIgARIgARIgARIgARIgARIgARIgAYoyPAMkQAIkQAIkQAIkQAIkQAIkQAIkQAIkYAEBijIWQOeUJEACJEACJEACJEACJEACJEACJEACJEBRhmeABEiABEiABEiABEiABEiABEiABEiABCwgQFHGAuickgRIgARIgARIgARIgARIgARIgARIgAQoyvAMkAAJkAAJkAAJkAAJkAAJkAAJkAAJkIAFBCjKWACdU5IACZAACZAACZAACZAACZAACZAACZAARRmeARIgARIgARIgARIgARIgARIgARIgARKwgABFGQugc0oSIAESIAESIAESIAESIAESIAESIAESoCjDM0ACJEACJEACJEACJEACJEACJEACJEACFhCgKGMBdE5JAiRAAiRAAiRAAiRAAiRAAiRAAiRAAhRleAZIgARIgARIgARIgARIgARIgARIgARIwAICFGUsgM4pSYAESIAESIAESIAESIAESIAESIAESICiDM8ACZAACZAACZAACZAACZAACZAACZAACVhAgKKMBdA5JQmQAAmQAAmQAAmQAAmQAAmQAAmQAAlQlOEZIAESIAESIAESIAESIAESIAESIAESIAELCFCUsQA6pyQBEiABEiABEiABEiABEiABEiABEiABijI8AyRAAiRAAiRAAiRAAiRAAiRAAiRAAiRgAQGKMhZA55QkQAIkQAIkQAIkQAIkQAIkQAIkQAIkQFGGZ4AESIAESIAESIAESIAESIAESIAESIAELCBAUcYC6JySBEiABEiABEiABEiABEiABEiABEiABCjK8AyQAAmQAAmQAAmQAAmQAAmQAAmQAAmQgAUEKMpYAJ1TkgAJkAAJkAAJkAAJkAAJkAAJkAAJkABFGZ4BEiABEiABEiABEiABEiABEiABEiABErCAAEUZC6Bz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vpPZs2fLrFmzlDlw++23V5y8Zz3imjKDBw+WOnXqyIYNG+TDDz+UKVOmqIL+2Weflbp16yY2+1tuuUWtknjrrbeU+eNKS1p7+ZpXUVGRdOzYUY488ki59tprVRovvfSS/PTTT6UpzZkzRx577DE5//zzZa+99ir9d/B82GGHJZL6EUccIY0bN05r8sybN0+6d+8u9913n7Rv397qejrfZazPrCB2IpimjBM0+JMEb3p/uIqbKQuUuIixfxIIUHdJoMgxTBCg9kxQ8ydGp5ApbzYwOdq1PUB+XLVcXv9LkbRuqD/3ectFjhhXKLXqNpT35vx3s61MSRbG69atk5NPPlkVo/fff3/alTnA4dVXX5U///nP+hPIQ8+ff/5Ztt1221hXjmvKzJgxQxXRUfvXv/4lI0eOlMsuu0wV7zatuLhYfvnlF7UiCaskYPjA1KtatarNsKWxmzZtEvypVq2a8XhJas84iQQCX3/9dTnvvPPU6rB27dqlHRH4gwesZkliC1G6i5RnyqD/McccI82bN5c77rjDatY632Wsz6wgdiKYpowTNPiTBG96f7iKmykLlLiIsX8SCFB3SaDIMUwQoPZMUPMnRqeQKW821113nfztb3+TGeeL1gqZsmPN+Fjk2HtFME7023z0SbIwxrklf//73+Wmm26SLl26xCLnxRdfFMRjpQ0MhVatWsnFF1+8WZEbFbb4bf/cuXNl8uTJ8sMPPyjzB/PC9qjUBmPl7rvvlueff16++eYb2WGHHdS2kksvvVT9d9R69+4tS5YskcgUwaoGjIltJchnwoQJ8t5778m3336rDk7eb7/91DYj/J3abE2ZTz75RE466STp0aOH3HDDDWpofC9gq9O7774rmM/uu+8uffr0UaseorZ06VK1UqNfv37K5AGOmM/111+vjICybcCAATJw4ED1zzgoGvNbtGiRMnD+8Ic/yKBBg6RJkyalYVjhhGsiJ2xNi86mgcGAhs+GDx+u8gNmK1askH333beUE+A/btw4wT2w6667ytVXX61Wa0Tag3mElTxvvvmmfPXVV4LzWbCa5KyzzlJb3VLbVVddJVOnTlV9oTNsF4I59Mc//lFdL5VXxAEb4IctXKtXr1aGJLaJXXHFFZutGoKWHnnkEfn888+VcXXggQcqneyzT8XL0YYOHao0Bo0UFham1X0mU+b7778XbGfD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" width="1000">


.. parsed-literal::

    2. Reporting Generalized Performance:
    
    |                  |           0 |
    |:-----------------|------------:|
    | clump_p1         |   1         |
    | clump_r2         |   0.1       |
    | clump_kb         | 200         |
    | p_window_size    | 200         |
    | p_slide_size     |  50         |
    | p_LD_threshold   |   0.25      |
    | pvalue           |   1         |
    | numberofpca      |   6         |
    | tempalpha        |   0.1       |
    | l1weight         |   0.1       |
    | Train_pure_prs   |   0.775701  |
    | Train_null_model |   0.230336  |
    | Train_best_model |   0.997862  |
    | Test_pure_prs    |  -0.0103224 |
    | Test_null_model  |   0.110899  |
    | Test_best_model  |   0.125464  |
    | Difference       |   0.872398  |
    | Sum              |   1.12333   |
    3. Reporting the correlation of hyperparameters and the performance of 'Train_null_model', 'Train_pure_prs', 'Train_best_model', 'Test_pure_prs', 'Test_null_model', and 'Test_best_model':
    
    3. For string hyperparameters, we used one-hot encoding to find the correlation between string hyperparameters and 'Train_null_model', 'Train_pure_prs', 'Train_best_model', 'Test_pure_prs', 'Test_null_model', and 'Test_best_model'.
    3. We performed this analysis for those hyperparameters that have more than one unique value.
    


.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <img src="data:image/png;base64,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5CIYLe34yjQY1fnVX/oCRV9IKQzV/x+hL4x0L9K2F+Gi//9E90z1hVZF6L6pL2z0haf+lofHRKb3DI5DAAEEEEAAgbILEIyW3ZAzIIAAAggkUEAzH/WBXqGhZimmK8UF5lGTpD61PerrcX4E4i6g97NCKc1upiCAAAIIIIAAAgggIAGCUcYBAggggAACaQQ0o0mz2LSct6gHZhCMMnQQiIeAZhprRYBm8qU+vCoeLaCWCCCAAAIIIIAAAlEIEIxGoco5EUAAAQRiKaBl4Vo2r60YtGRZ4age2lVUIRiNZTdTaQQQQAABBBBAAAEEEEDACRCMMhAQQAABBBD4P4HgoR96iI9CT+1PqYf7EIwyRBBAAAEEEEAAAQSKEtDD1x588EEbM2aMe/iingFw2mmnuQc9aj/04ooekqe9yrUvtfZo1t7p2iO8qKL94u+55x734L6qVau6B97169cv7QonnVvPOfjyyy/dcwr0sE7t5a09kSkIIPBfAYJRRgICCCCAAAIIIIAAAggggAACCCBQSgGtMHrmmWesW7du7kGCb775pnsQat++fd2e9cUVrVDSaqU999zTBaMVK1YsMhh95513rGfPnu4hjyeccIJ7iKYeJKkHueohg1WqVCm4lB5yqsC0devW1qFDB3duPYhx3333teHDh5cY2JaSgpchEDsBgtHYdRkVRgABBBBAAAEEEEAAAQQQQACBfBDQzM3OnTvbWWedZf379y+okmZmvvrqq+6fZoEWVZYuXep+X6lSJbeN08KFC4sMRnWdlStX2vjx4wtWNWk//N69exd6YKj21tZMUgWmTz75pDu3ysiRI92KqEGDBtmRRx6ZD3zUAYGNLkAwutG7gAoggAACCCCAAAIIIIAAAggggEAcBe688063jH7y5MkuiAzK+++/b6eccopde+217r+ZlOKC0QULFtgxxxyTdhbqUUcd5ZbKa9aqimas9urVy2699VYX2gYleBjhIYccYnfffXcmVeIYBBIvQDCa+C6mgQgggAACCCCAAAIIIIAAAgggEIWAZorOmzfPhZHhohCyefPm1rVrVxswYEBGly4uGNUDQi+77DJ7+OGH7aCDDip0Pv385ZdftpkzZ7rZoYMHD7a77rrLJkyYYA0aNCh0bPfu3e2HH36wiRMnZlQnDkIg6QIEo0nvYdqHAAIIIIAAAggggAACCCCAAAJpBZo2bVqijPYALap07NjRNt10Uxs7duwGh7Rp08btHaowM5NSXDA6dOhQu+2220x7h2qP0XDRz/X7t956y7bffnu74YYb7IknnjDNWtVDRcPloosucrNb9bAnCgII8PAlxkAMBBoMvCMGtaSKCCCAAAIIIIAAAggggEA0AvMvvDSaE2+ks677ttFGuvKGl212xKYl1qW4YPSII45wYaT28kwtbdu2tTp16tiIESNKvIYOKC4Yve+++2zgwIHuoU7169cvdD49pf7+++93+5nWrl3brrrqKhszZox9/PHHLrQNl8svv9yee+45mzNnDg9gyqhXOCjpAswYTXoPJ6B9BKMJ6ESagAACCCCAAAIIIIAAAqUWIBgtNV2JL6y407wSjynuAGaMlomPFyOw0QUIRjd6F1CBkgQIRksS4vcIIIAAAggggAACCCCQZAGC0eh6t6zBaEl7jHbp0sVuuummjBpQlj1GNZN01qxZGe0x+v3339ukSZMyqhMHIZB0AYLRpPdwAtpHMJqATqQJCCCAAAIIIIAAAgggUGqBpAWjf37bsNQW5f3CTXb6okynvOOOO2zIkCFFPpX+mmuusVNPPTWjaxQXjM6fP9/atWtX5FPpt956axs9erS7ztSpU+3ss88u8qn0Bx98sGn5PQUBBNhjlDEQAwGC0Rh0ElVEAAEEEEAAAQQQQACByAQIRiOjtbIGo5999plpVqhmjvbv37+gohdffLGblal9P3fccUdbtWqVLVmyxLbddlurXr162gYVF4zqBZ06dbKff/7Zxo8fb5tvvrk7x5QpU6x3797u2qqDypo1a0z7m9aqVcueeuopq1ixovv5yJEj7frrr3d7lR599NHRoXJmBGIkwIzRGHWWr1UlGPW152k3AggggAACCCCAAAIISCBpwejqpQ3ypmM323l+meuihx3pqfTdunWzvfbayz0dfsKECXbBBRe4GZ4qM2bMsNNPP73Qz/Tzd9991/1T0YzPlStXWq9evdz/rlmzpnXu3LmgftOnT3fhZ5MmTezEE0+05cuX2/Dhw13wqoctVa1ateDY//znPy4sbdOmjbVv394WLlxojz76qLVo0cIee+wxHrxU5l7nBEkRIBhNSk8muB0EownuXJqGAAIIIIAAAggggAACJQoQjJZIVOoDyiMYXbt2rT344IMuHF22bJmbqanl85oB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qlqrhx+wv/4vAB3UrqObTaoQdea3SwuO1wzVd87u42aZdhg1osQ+4IDcCjBWcusd96sxXuLeg7mtP+Mlt95xvhpjJc69l9u6M1Zy400wGp0zwWh0tpwZgTgIxDIYXbRokSk0feutt2zlypVWq1YtO/744+3ss8924WFQXn75ZRs6dKjNnz/fFILusMMO1rp1a7vxxhvtm2++scMPP3yDPmrZsqWNGJFZUBQEt88++6yNGjXKBbirV6+2Vq1a2dVXX2116vxv35Tg2LlzC4dZQT1uvvlm69q1q6vP2LFj7corr7QhQ4bYe++9Z5pJu2zZMps0aZLVrl3bVqxY4dr/2muv2fLly1272rdvbxdeeKFtttlmGY+74Np9+vRxwbOut3jxYlfv8847zzp27FjoXNrq4Nhjj3U/HzhwoH3++efWu3dv69u3r82ePdvuuece++ijj+yXX36xatWqWfPmze2qq66ymjVrZlyndAcmJRgNlsCf/+LzNuGLeYWa2n63Rjao3bElzuAc0627tdippu370P22ImU/0RvaHm6nNfubtX/iUZuz/Ad3/hm9+ljVypWt+eB77b9zRf9Xhh3XxdrWa2DNBt9rv65ZU6Y+4sXlK8BYKV/PpJ+N8ZL0Hi7f9jFeytczyWdjrCS5d8u3bYyV8vUs6mwEo9E5E4xGZ8uZEYiDQOyC0YULF9pJJ51km2++uQtDtcT+nXfesfHjx7uf33DDDc797bfftp49e9p+++1nRx99tJtpqUB18uTJ9uKLL9rvv//ugkwt0993332tW7du7nXbb7+9HXDAARn1XRB2NmnSxLbcckt3ne+++84ef/xx22abbWzcuHEuIFQpTTC62267uaBTYaSCXQWnFSpUsBNPPNF+/fVX116Fjp9++qmNHj3a2rRpYw8//LA7JpMSBKOqv+p96qmnunYomJ03b57deeed1qFDh4JTKRht0KCBC2O7d+/urr3zzjtb06ZNrV27drb11lu7um277bYuyFVwfemll9o+++yTSXWKPCYpweh97Y+1dg0buWX0n6Y8aKnpDjVsXPce9sjMD+yGNyYXafHeOedZ5UoVrfngQRsc06vFPvbPg9raueOfs5e//MI2r1zZPjn3Qpvzw/fWfuRjGxx/7SGH2RnNW6StT5k6jBeXWYCxUmZCr07AePGqu8vcWMZLmQm9OQFjxZuu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YINzvfXWW/bcc8/ZFVdc4c4TFAWjlStXtvPOO6/Qa7bYYgs74ogjogPlzAjETIBgNGYd5mN1CUZ97HXajAACCCCAAAIIIIAAAoEAwWh0Y6GswahmeXbu3NnOOuss69+/f0FFL774Ynv11Vfdvxo1ahTZAL125cqVNn78eNt8883dcVOmTLHevXu7GaJnnnlmsY3X799991174403bLvttisUjO644442atSo6PA4MwIJECAYTUAnJr0JBKNJ72HahwACCCCAAAIIIIAAAsUJJC0YHftli7zp8K67flimutx5551uGf3kyZPdLM+gvP/++3bKKafYtdde6/6brixYsMCOOeaYtDNLjzrqKNtmm23cTNSiytKlS00zQw855BAbPHhwocP0cwWjjz32mK1evdq23HLLMrWTFyOQVAGC0aT2bILaRTCaoM6kKQgggAACCCCAAAIIIJC1AMFo1mQZv6Cswahmis6bN88tnw+XNWvWWPPmza1r1642YMCAtPV5/vnn7bLLLrOHH37YDjrooELH6Ocvv/yyzZw50ypVqpT29QpD77rrLhs4cKAdffTRGwSj2u90/fr1tnbtWtt2223t2GOPtUsuuaRgZmrGSByIQIIFCEYT3LlJaRrBaFJ6knYggAACCCCAAAIIIIBAaQQIRkujltlrrum0tsQDP/nkkyKP6dixo2266aamhyOlljZt2tiee+7pgs90ZejQoXbbbbfZuHHj3B6j4aKf6/faQ3T77bdP+3qFoT/99JNNnTrV1SFc/v73v1uLFi1s1113td9++81ef/11mzBhgvvZiBEj3P6jFAQQMCMYZRTkvQDBaN53ERVEAAEEEEAAAQQQQACBCAWSFow+8+U+EWpld+rrO60u8QXFBaN6kJGCyyeffHKD87Rt29bq1Knjgsh05b777nOzPfWgpvr16xc6RE+pv//++90epbVr197g5ZpJetJJJ7kn319zzTUltkEHaHapZpnefPPNbiYrBQEECEYZAzEQIBiNQSdRRQQQQAABBBBAAAEEEIhMgGA0Mlo7cdf3y3TyjTVj9LrrrnMPVtIepM2aNcuoDb/++qvts88+pjrfcccdGb2GgxBIugAzRpPewwloH8FoAjqRJiCAAAIIIIAAAggggECpBZIWjD71xX6ltijvF57U8N0ynbKkPUa7dOliN910U9prlLTHqGaS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fHkEo5qQ9eCDD9qYMWPs+++/t1q1atlpp51mp556akaTtN5++2275557bPbs2Va1alU79NBDrV+/fla9evVCDTnssMNs8eLFGzROE8s04Sy1vPjiizZkyBD78ssvbZtttrF27drZxRdfbFtssUXmQByJQMIF8jIYnTFjhp1++ul244032i+//GKjRo2yJUuW2PDhw61Vq1buv9o3VG9u/X7nnXe2Dh062HnnnVcoOE23lP6KK66wZ5991t566y279dZbbfLkyS5wPeSQQ+y6666zatWqZdXlujHtuOOOpvPefPPN7ka29dZb2/HHH28XXHBBocBQx7Zs2XKDLQKCcHju3LkF1+7Ro4ctXLjQHnnkEXfe999/3/bcc08bMWKEO+add96xwYMH26xZs2zt2rXWqFEj69Onjx1xxBFZ1T+4tkzkrIB69erVzvnqq6+2OnXqFJwv8NSN9b333jPN/l22mBvKSwAAIABJREFUbJlNmjTJ7eH65JNP2siRI23RokXuNTvttJMdc8wxdtFFF2VVp9SDkxKM6pv0ocd1tZEfz7KrJ08q1MwBhx1pWh7fcdSIYpfUX3HAwW5Z/IUTXrAXPv/feNHJJvXoadWqVLVWDz9gf/1fADqoXUc3m1Qh6sxvlxZcUzNU3zm7j5tl2mHUf8cUJX8EGCv50xdxqAnjJQ69lD91ZLzkT1/ke00YK/neQ/lTP8ZKbvoiacFoz3d75gYug6sM3294BkcVf4g+Oz/zzDPWrVs3a9asmb355pv20ksvWd++fV0uUFzRZ/uePXuath084YQTbMWKFTZs2DCrWbOmjR492qpUqVLwcmUKlStXdtlHuCjoTM0Cxo0b58LV1q1bu7zkq6++sscee8z23Xdfl6mwqrbM3c4JEiKQ18Hobrvt5kI/hYy6GRx44IHWoEED918Fmfq9ZpAqNNQ3Iccee6zddtttBV1TXDCqkFGhnwLABQsW2BNPPOFuFv/+97+z6lrdmHRDUUDbsWNH23XXXe2NN96w119/3U4++WS7/vrrC93EsglGFZTq26IDDjjAmjdv7mZmnnjiifbyyy/bJZdc4n521FFHWaVKlWz8+PE2c+ZMu/32251DpiUIRps0aWJbbrmlHX300fbdd9/Z448/7r5R0s00CIsDT7lvttlm7jr6Zqxr164uYL7qqqvs8MMPd/2jItc5c+YU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Ski3bt2kUaNGYZ8PPdGZVKNDhgyRpk2bhn1r//79pl9LliyRo0ePytVXXy29e/eWatWqRdQvTkLAJoE0C0Y//fRT6dixozz//PPy559/yrRp02TXrl3y1ltvSeXKlc3/6ntDv/vuO/P9AgUKmB/W+++/Pyw4TeoR7Mcff1zee+89WblypQwePNj8cGrgevPNN8uzzz4rcXFxUZlrgJk/f3555pln5IUXXpCNGzfKxRdfLC1btpQHH3xQMmbMmNCevkbggQcekF69eoVdQ/u0du1aWbx4ccLX0yoYDY13+fLlZrzLli0z12jQoIE89dRTcurUKRk2bJjMmzdPDh8+LNWrV5cBAwZI7ty5w/q4YcMGGTFihOj/6mdKly4t9913n3FLfOgKXQ271f7gwYPmpqd9GD78P483nPmOUV0tPG7cONm6daucPn1arr32WuNWsWLFhGZ5lP7skgw9At9z3gcS/+03YSc0vPIqGdWgSYorOGe1bisVLrtcKr75uuw/432i/WvVlfblykvDKRPl632/mfY/7dJDLsycWa4bM1L+s1b0f8f425tLraLFpdyYkXLoxImofoY42VsBasVb36C1Tr0EbUa9HQ/14q1vkFqnVoI0m96OhVrx1jfUOsGod87pFYxqNtKqVSu59NJLpX379nL8+HGZOHGiGej7779vvp7coX8X79q1q2gWo/uDFCtWTOLj42XVqlUmR2jWrNlZfz/v3r27lCxZMqzZ66+/PiyI1X5oNvLzzz9L586dTc4wc+ZM8/d/zQE0vOVAIEgCaR6MXnnllXLy5Em54447JFu2bFKjRg0pXry4+V8N5PT7uoL0888/N8FekyZNRH9DETqSC0bLlCljfmA1aP3+++9lypQpJlzV32REc2iAqeHnsWPHpGHDhuYGouGjhpwaMOrNKXSkZzCq49W+3XDDDfLZZ58ZL73h6aZOGnSq5zfffCPvvPOOGcfLL7+c0G/1veeee+SSSy6RO++808yF2urN96WXXjLnhw4NW999912pXbu21KxZ05zzwQcfmMD5sssuCwtGJ0yYYDaT0mvruTrXepP88ccfZfz48WbzKT0IRs+uyNcaNpEGJa8yj9FvPmOjpbJ588mcth1kwoYvpP+yJecs53Xd7pfMmTLKdWNGnXVOlwo3yJM1a8l9c/8lC777Vi7KnFm+uq+3fP3br9Jw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" width="1200">



.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <img src="data:image/png;base64,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nzx0o4GPXtDGB1U3fWrTxqYGVCKibY1ulNPVC4z96xEDY8wl/3GUNkIMcZbGoMjz8C6SNmNBAAIhEsBcakG1RYsWyZQpU+T555+XTTbZRLp27Spjx46VDh06ZNL6008/lauuukpmz54tK1askI4dO8rpp58uPXv23Kj98uXLk7HmzZsnn332mXTu3FlGjRol3bp12+jauXPnyvTp02Xx4sXSpk0b+da3viXnnXde0v+Gx+OPPy6/+c1v5P/+7//k888/T+bev39/Oemkk5I15T183FTmnUtI7diY6FELLey0CLEGoDnmEuaSXd7Xt7bJf/LQjQY+ekEbH1Td9KlNG5saUImItjW6UU5XLzD2r0cMjDGX/MdR2ggxxFkagyLPw7tI2owFAQiESABzqYJqr776amLCbL/99jJo0CBZs2aN3HrrrcnVM2fOTP7e0lFXVyfDhg2TZ555RgYPHix77rmnzJkzRxYsWCCTJ0+W3r17NzQ3fffr10+WLl0qQ4cOlW233VZmzJiRmEc33XRTYmrVH/fcc49ccMEFst9++0mvXr1k9erVcscddyTzM20aG1+zZs2S0aNHywEHHJBc27p1a3nsscfkiSeeSOZk+sl7+LipzDuXkNqxMdGjFlrYaRFiDUDzdZrDAQ522c+TS7b8tLanNmhVRl/d9rEHIP78xx+MYeyCAOaSC4p2fZDLdvyqbQ3vaolxPQQgEBsBzKUKio8YMUKeeuqpxBDaaaedkquWLFmSmEIDBw5MNWbM00Wmj4kTJybmlDm++OKLpO2bb74p8+fPl8022yz5uzGtLrnkErnhhhukR48eyd/M00vGEGrXrp0Yk8gc//rXv+TQQw+VHXfcUe677z7ZdNNNk7+/++67ctxxx8lhhx0mV155ZcOKjGH13nvvyaOPPtowljG9fvjDH8rf//735ImsvIePm8q8cwmpHRsTPWqhhZ0WIdYANF+nORzgYJf9mEu2/LS2pzZoVUZf3faxByD+/McfjGHsggDmkguKdn2Qy3b8qm0N72qJcT0EIBAbAcylZhRftWqVHHzwwXL88cfLpEmTmlxx2mmnJU8UmSeQWjrOPfdcMQbTwoULpW3btg2X1j9NZF5r17179+Tv5hV1y5YtS0ygxsf1118v06ZNk4ceeih55d3f/vY36dOnj/y///f/5Kc//WmTa8844wx5+umn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45I7zr+fNUAooN3UYy1XYwwVcykG1dbqLBNnvAdojXH93WV4l2+NkhCAAATSI4C5VKfZGmusISuttJIcd9xx7q+fffaZO2JO/znvvPOcubTpppvKXHPNJZdddpm32q2aS5988omsuOKKsskmm8ixxx7bqf2dd95ZnnvuOfcdpV69eslLL70kG264oeyzzz4ydOjQTveut956bmeWHqVX9orxUFm2LymVa/fC5IV/vCR7LzesIbKz//57WWCZH6SEtHRf261F6Y4bKZjiHIDmIswBkxKIeCg/mfjkP9zLc49dEm1iEy5fvzVtfOaARhSsjbG8WnZLjnr0Bdl3xeE8B0SUKIc4xlyKGEBNVl0mzngGaBJuN7eV4V2+NUpCAAIQSI8A5lKdZksssYSoSXPAAQe4v/7tb3+TXXbZxX27SA0ZvY466ii57bbbnInje03JXPriiy/ktddek+mnn17mmGOOjqb0SL7/+7//c99Sqv1dd1sNHDhQttlmG/nVr37Vca8aYR999JHccsst8u1vf9v9/YEHHpA99tjDfbtJx1b2ivFQWbYvKZVr58JEz1Y+addz5IFr/tYQ2ZpbrSwHXrCXTJvBN5jaqUVKMduorynOAblrzhzQOZpzjwefecgn/+HuQz5uWbSJy9endmva+MwBjThYG6OPXhbLsgaoRpUc4hhzqZpYmlIrrcYZ+e+nWau8/VqjNAQgAIH0CGAu1Wm28soru+PujjjiCPfXk046SS688EJnMs0000zubyeccIJcccUVot9nKntdfvnlzvTR7zlddNFFsvzyy7sdU3qtvfba7oiemvG02Wabdeyk0t9ffPFFGTRokMw+++yy/fbby/jx4+Xiiy92ZfX4vnojSk0oNZC0Pi3z7rvvyh//+EeZc8455frrr5fpppuu7BAkxkNl6c4kVLCdC5MP3/lIhm94tDz/+EsNiS24XH855tbhMuNsE7851pOvdmrRE7imOAfkrjlzQOfMyz0efOYhn/yHuw/5uGXRJi5fn9qtaeMzBzTiYG2MPnpZLMsaoBpVcohjzKVqYmlKrbQaZ+S/n2at8vZrjdIQgAAE0iOAuVSn2Q477CCvv/66jBw5UqaaairZeOON3RF4V199dcdd++23nzzzzDNy9913l1ZbDaQxY8Z0W16Pu9t8880bmktaSNs/8cQTncHVu3dvd1TewQcf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atSoFPUO1ucUXywHG7xHRSxMPOAFLooWfkBTnAPQHHMJc8kv7zGXwvCzWgtzpFVlMJfsKpNWz8jx+HrlwBhzKX4cFbWQQ5wVMajyd3hXSZu2IACBFAlgLtWptssuu8jDDz8s++67rwwcOFDmmGMO+da3vpWirtH7nOKL5ehQmmiAhUkTkCq6BS38QKc4B6A55hLmkl/eYy6F4We1FuZIq8pgLtlVJq2ekePx9cqBMeZS/DgqaiGHOCtiUOXv8K6SNm1BAAIpEsBcqlNtqaWWkrXWWktOOeWUFLWstM8pvliuFFCDxliYWFCBF+whVEhxDiD/iH3MpRDZL+KT/+RhGA1i1II2MaiGqdOaNj5zQCMi1sYYRjlbtcA4vh45MMZcih9HRS3kEGdFDKr8Hd5V0o7b1pdffil/+MMf5Prrr5e3335b+vXrJ9tvv71st912oidUTem69dZb5f7775d//vOf8sorr7jNCH/5y18mK/LBBx/IDTfcIPfdd5+MHj1a/ve//8l8880nG2+8sey4444yzTTTdCqzww47yKOPPjpZPbrRQT/X0ugaP368bLLJJq4ve+65pxxwwAGT3fr555/LBRdcIH/+85/lP//5j3z729+WhRZaSA4++GBZbLHFuq36vffekw033FB0HEcddZQMGjSo4z4dz4033ih//etf5bXXXnObMb7//e/LTjvtJBtssEGn+pTviSeeKE8//bS89dZb7pMzylvr1vu/853vdLpf+9XdtdJKK8nFF1/c6af333/f6aiM//vf/8oss8wiyy67rOy9997Sv3//Tvfq72eccYb83//9n7zzzjsy++yzyyqrrCJ77bWXzD333B33Pv7443LRRRe5U8reffddmXbaaWX++eeXbbfdVjbddNPJ4kNZnH766fLkk086VnPNNZf8+Mc/lt13311mmmmmjnr1hLQ//elPLm607nHjxskf//hHWXnllScb7j333CNXX321PP/886I6TD/99G48O++8s6y99tpTjM/TTjtNzj77bJlzzjk7xeXXX38tN910k2jdzzzzjGOg415hhRVkn332cTxCXJhLdRRXXHFFt2PpkEMOCcG2R9cR46GyRwP7ZnAsTOyojBZ+WqQ4B6D5RM3hAAe/7Mdc8uVntTxzg1Vl7M3bMdYAxF/8+IMxjEMQwFwKQdGvDnLZj1+rpeHdKjG79x9++OFy7bXXys9+9jNZYoklnEly++23u5fsQ4cOnWLH1QT617/+JYsuuqgzdHr37t2tuaSGx5AhQ5yBoe+Y1UR57LHHnMGzzDLLyGWXXdbphCyt97nnnpPDDjusU/tav5pHja6zzjrLGUdqXnVnLunf9XSuF154QbbccktZYIEF5JNPPnEGhxo8a665ZrdVDx8+XG677TZXb1dz6fe//70zStREWXLJJUWNCzXd9P8v/PznP5df/OIXHXVqu7/5zW9k6aWXdmZGzSwbMWKEM7i0nqmmmqrjfv2b3rvNNtt06peaeGow1S411TbffHNnlm299dbOfNH/+6qrrhI1D9VImXfeed3takIpwy+++MLVq+bWiy++6AycGWec0fW9ZnKp4XjXXXc5003bVGPuwQcflAceeEAGDx7cSR81jLbYYgtXh/ZBzS2Njeuuu04WXHBBZ17WTkHT8aq2asL17dvXmUyNzKVzzjnH6bPIIovIrLPOKp9++qnccccdju+hhx7qTLnurldffdWZl8pTDal601Pr0LjTuNXNNGqCaf+VwXTTTef6Wm+ylc1ezKU6ciq4OsMqfpFrXRa4j1OuTu+FF17oxB8zZoxzQ9ddd13nUM8wwwyduqTtaLDoxKmOsrquOpnstttusvrqq5ftfke5GA+V3p1KoAIWJnZEQgs/LVKcA9B8ouZwgINf9mMu+fKzWp65waoy9ubtGGsA4i9+/MEYxiEIYC6FoOhXB7nsx6/V0vBulVjn+8d9+pmMHzde+kzXR6brO61fZR6l9aW9biZQw2XYsGEdNe2///5uV4f+o6ZCo0t3wNQ+naKGkL7Q727n0uuvv+6q0N1K9VdtZ8mZZ57p3uPWrinV1agv2oaaCWpinXTSSd2aS8cee6wzO/Sd8A9+8IOmyOn8rru4lMnJJ588mbmkxojWVb/rSA0m3ZGlO3/UrFOjZUqXvtM+/vjjnTG22mqrddyq5pIaQbrbaUqXmj177LGHqFGo7GqX6qc7l3RX1q677ur+fMUVV8iRRx4patrU7/y59NJL5eijjxbVpOuOq65tq2mm49LdZWoO6aX913GokVW/40qZ6y4r5b744ou7e3VXk5o+yky9BjWJGplL3Y1bfQA103Q3lfaju6s2XjXeusal/k11W2655ToV1U8CqVmlDJWl74W5VEfw448/dtvNvve978lBBx3ktpOFvnyccu2TBq+6zLqFTU0jTRZ1RtWlnXrqqTu6W2tHJ5wf/ehHzvHUAFeHUrfurb/++l5Di/FQ6dWhRAqzMLEjFFr4aZHiHIDmEzWHAxz8sh9zyZef1fLMDVaVsTdvx1gDEH/x4w/GMA5BAHMpBEW/OshlP36tloZ3q8Qm3v/xe5/Im6+MlT8df6P8d/SbMnf/uWSrgzeVub43h0w/S+cj0cq10FopNUtqR6nNM888HYXVFNGjz4444gj372auMoaQ7k766U9/Kvvtt58zQWpXrS7d8aRHpqmBUbTZQQ0Pfcd73HHHyTrrrDOZuaQ7lFZddVW3W0eNNN18oLt3dKdKo6tmYgwYMMCZGbpbp+vOpUZldTeW3qvvp7uaGF3L6E4xZaCfo/nJT37S8XPNXDrmmGNcf/UIv+6uW265xe2Q6moMqYGiR/j9+te/dgaZXuedd54z3+rNHv277sxSA01/X2ONNaYouZpTOq6HHnpIZpttNnevtqE7r9SgqTfT1DA74YQTRPv4wx/+cLJ6y5hLWonqrW3pGLvj+ctf/tJ5BRrDjUzP7gapvsLCCy882bGDzeRA13swl+qIaFJqwunZkHrpbqCu50Dq3zXR77777pZ5+zjlusVOt93plrvf/va3HW3feeedbguLmGteAAAgAElEQVSnbjesbR/UiWT55Zd3k4yeLVm79MxG3bWk2zN1UvW5YjxU+vQnlbIsTOwohRZ+WqQ4B6D5RM3hAAe/7Mdc8uVntTxzg1Vl7M3bMdYAxF/8+IMxjEMQwFwKQdGvDnLZj1+rpeHdKrGJxtK1J98sVx0zYrLC2wzfXAb9YpPKDSbdsaTfsum6+0N3dugRb2qo6G6WZq4y5pIesaYnSen7XH2vW7u0Lp1XdXfLZ5995t5Br7feeqKGgR6N1vXSd9H77ruv+66TGlHdmUu13T1q+Oi3hvRoNX3XrUezqTGj9Xe9LrnkErcRQc0f3ZTQirmkBo4aNWradN0lpXz1HbX+Ww02NY/efPNN0XfZ9Rs61FxS80uPo9PdUPqbHl+oR/7VH5+nO8i0/3rEnZo8tWPx1GjTXUJ60lftZK+a4bTUUks5k612LJ4aRnoSmJpG9XUrEzXttA/6b2WnMaHfXtJvTdUuPU5PTxDTY+b0fbzqpN+W0v7o5g7l2N3VrLmkG19Urw8//NB5D2qMalv6TaX6S48u1M0nG220kdux1UpcqiZ6bKOy1Pp9L8ylOoJFH8iqh33vvfe2zN7HKdetdbrF7sorr3QfKqu/9FxKdRv1N73UHFOXWt1aDe7apQmqLrJuPVSX1+eK8VDp059UyrIwsaMUWvhpkeIcgOYTNYcDHPyyH3PJl5/V8swNVpWxN2/HWAMQf/HjD8YwDkEAcykERb86yGU/fq2WhnerxERe+MdLsvdyk46e61rD2X//vSywTHNHtbXeevcl9FSnPn36uKPJul76TR/9Jo3uPGnmauUlvtZXOzpODQg1C2o7YPQ3PSZNjRQ1VyZMmOB2qOhOGz1WT/9d/wkU3dmkRoK+u9bTqvRbQ92ZS7X3xzPPPLM7yk9P6NJL/64Gjx4Tp2ZF7Ro7dqwzKfS7U3rvI4880rS5pMe16Q4kNW66Y1szVGptqcGl76lXXnnlTqjVSNI+fPe733XfSlIDR3cL6albXc0a/TaSmnS1jSFakRpI+h2qerb6d31Pfuqppzqjpnbp2PX9fHe7ow455BBn3NUuNWDUYKp9x0n/rjrpRo6LLrrI7TarXToG3fhR+95S11hq1lzS+NJj+PTSuvQbVzpe1bP+0m9g6be81BBUo7GVuKwd09j1yMBm4r+7ezCXypIrUc7HKa9t59Mg14971V86EWpAP/HEEx3bJ3XC0e8yaWDrLiZ1XfVMSHWsdUJRZ97n0kWlJlTtzEmfunIqW5t4prQdNSce7RxrT9FCty2340pxDugpmvvqDYeJBHsKh3bMAT7531O4++ahxfJoY1GVKc9X7ch/7ZHPHNCIMvEXP/5g3PMYt2MO6Jr/xFX8uOraAsyrZW6Vdzvyvxny+o2lk3Y9Rx645m8Nb19zq5XlwAv2kmkr/AaTvqBX40G/Td/1WnPNNZ2Zo8e7NXO18hJf66ttNOj6naBGbel3kvReNXt0Z0zt0nr0N32vq6ZTI3NJd7ioeaC7c9TMmn766V0VuiNGOcw999wycuTIjnoPPPBA0ZO2dHeOfm6lWXNJd9foe259F61cF1tsscmGpMbViy++6HYv6X1at34zqehbR1qR7t66+eabRb+RpEe41S49yvD888+XJZZYwply+g0q/d+qr37PqP6oOt0Yon1TM0s1VnNN343rN5H0RK9pppmmU5+1r9rnd955R3QHmBpYhx12WKdvK2kB3SGlOqjRp+3q2HT3l45Lj8br7mjDZs0lPblMzTDthxpHvXv3ll/96ldSf5zjCy+84L4hpptQ9LhFvZqNSx2X7gjTWKg/7ayZ2G90D+aSD70Wy/o45Toh6Mfa1NXWj27VLg0orVcvTVKdPPQaPXq0S8Rnn322497ZZ5/dObm+xpJWGOOhskWcSd5udWGSJEzPTvcULdq1qExxDugpmnuGfo8xVeAwkUA75gCf/CcPfSM3Xnm0icfWt+ZG2rQj/2M9BxB/vlFSXB7GxYx876iacTvmAMwl3yjxL191nPn3OO0arPJuR/43o+SH73wkwzc8Wp5//KWGty+4XH855tbhMuNsMzRTZZB7fN7Hdu1Asy/xtdzll18uv/vd7yb7xEnRoNRM0V0+NTNMj6pTE0G/q6PfFtKrkbmkO2p0V4se9afmQ/1V25mjBo0ewadHv+24445uE4JuXNCrGXNJd2PpEXtqfqiZsskmmxQNyf2uu47UNFMTqOvupa4V1L5Tpd8c0rb00t1f+kkY3XGjp3LVLr13s802cyd4qRmkl75H1+87qZG2wAILdNxbOzaw6/v17gagRwvq94x0nDXTSs2u2hGC9Tul9Jg9PXJP37urcdP1atZc6lpOv9Gl7/71W066+04vjUHVQNusXc3E5ZNPPul2p+nxhTqOUBtGMJeaCv8wN/k45Xo+pW411O2B6mDrOY7qzmqg6we71DHWBJlrrrlcZ9966y3njutkoZOSusQ6qenko1s91eH1uWIch+HTn1TKsqXajlJo4adFinMAmk/UHA5w8Mt+jsXz5We1PHODVWXszdsx1gDEX/z4gzGMQxDgWLwQFP3qIJf9+LVaGt6tEbO6c6noJCk1JvR7QM1czbzE13rUTBg+fLh7l3viiSe6HSjNXtof/aaO7o7Ra6+99nLfQtL3ubVdMfrtIjVTtt9+e2cYqNEx7bTTuqPSdDfS7rvv7jYd1F/aD93lc99997mdMGpYzTjjjJ1MqKeeesqZOVqH9l2P7dMdTbVLT7HScen49Ig77UOzlxoiyyyzjDvOT7/VNKVL32PrZ2H0uDk16PRSQ0i/7VR/cletDj1WT3ci1XZlqU56bJ8esVd/af/18zKrrLKKM4KmdCmLWvv6b710t5J+60k51l+620nrVLNO+XS9yppLtW88qWmo9dcMOu17vcmsmr3xxhvOkNQj/+p3cGlfdHea9k3jRP2Brr83q2F392VtLmlQalKqAApX/3czl5ZpdtKpr8/XKX/55Zddcj/zzDOuWu2HTji6tVGD6+9//7vb7qhH4Glb+o/eX7tq53PqPfUfI2tmzF3vifFQWaYfqZVhYWJHMbTw0yLFOQDNJ2oOBzj4ZT/mki8/q+WZG6wqY2/ejrEGIP7ixx+MYRyCAOZSCIp+dZDLfvxaLQ3vVonZ/OaSGhn6uZGaqVIble7g2XbbbVsySZoxl3SXiRo7a6yxhpx55pky1VRTNQ1SDRjdJPDDH/5QrrrqKldu00037XiOb1TRueee676lpBsQ1ltvPbebSM2k+kv7pH3TcasBsdxyy7l3ylO69Gi6BRdcsOMW/f6PfstI3zfrEXetXLoxQo0d3SXV1ZzpWo+++9bdV3qE2wEHHOB+3nXXXd27bzWXupp1Ombd2aPmml76vSa9agZdrf6vvvqqow96NN6ULv32keqt3NSs00uP/1N99Hi9+ks3eay++urObFPTretV1lzSz+PojrNTTjnFmX21b2pNqd9dv1Wlx/2pCakbUHS3kxqGIa+szSV1+NSgURdQtxs2u61Uy6jj1+oVyil/5ZVX3LmPel6k7lRS91S/r6QfO9OrFnh6FmfXHUq6hVJdzJoR1eoYavfHeKgs25eUyrEwsaMWWvhpkeIcgOYTNYcDHPyyH3PJl5/V8swNVpWxN2/HWAMQf/HjD8YwDkEAcykERb86yGU/fq2WhnerxEQ+fu8Tufbkm+WqY0ZMVnjbwzaXLQ/YRKaf5TutV+xRQj8Zov9xvr6XHTZsWEdN+++/vzs+7Z577nEv3PU/ytfdHzPPPHPDnR1F5lLtODY9cUoNrdpRZl27rztzdEdQ1+/+6O4kPWpON0LokXB66fF1en/9pbty1MTQ7/yokbTUUku5jRN6bbHFFm6n05133in6iRS99D2yGjCLLLJIx3Fq999/v3z55Zed6n3++efdN5vUJNGj61ZccUVnSuh1/PHHO1Ol3vDpThbt26yzzjrZT/o+Wt9L139P6r333puMtRpAek/tm0lqSOmlp3fpt7FOPfVU0Z1KtUv/f5OahPodouOOO879WXd7qZn4pz/9qdPnYdRcU7b1fWjUX71P76//7pPu9nrttdeciTXvvPN29EG1VhNTjyJUU6zrVWQuddcH5aDxpuPTTSXqA2jbqlHXS5moxkcffbTMMcccHX6Amo2q5be+9S2ne32fPVKqU9GszaVQEJutJ6RTXmtTj8nTcyY1qXTy0UudVz0ST5O2loC1+/VIPTWdHn74Ya8tcDEeKpvlmPJ9LEzsqIcWflqkOAeg+UTN4QAHv+zHXPLlZ7U8c4NVZezN2zHWAMRf/PiDMYxDEMBcCkHRrw5y2Y9fq6Xh3SqxiferwfTmK2PlmhNulDdGvyXz9J9TfnbQpjLX9+ao3FiqjaB2lJv+B/qLL764+w/09Yi1epOh9r2h+r9p+ccee8z9o9d1110nH374odtFo5ceL6emhl7//Oc/3Yt8NY0OPvhgmW666ToB/O53v9vxnlbb0h05uhtF/64bGfRvaiLo5gfdtaS7ixpdjb65pPfrXD148GBnMug3ivTS+vToNjVnllxyyYb1NvrmkhosalzojqrudizpcXdqfuil9ykvfV+tZoaadrpbSsc2//zzyzXXXOOO49PrjDPOcLuL9Lg5ZfnRRx+5bxzpziX9vpQaSrVLjRU1bj777DPZaqut3I4q/ZvupKppo0fW1TNQc0+NJ+2bfptJ255pppnciV41A0yNRzXmlIuajPquXfuqeqp5p2Zb7dKYUVNS71e2+m/dSaX1KZvrr7++wzDUXWHKWy/dpKJmn/a/xkmNIz1dTC/dzaWGpJp/2q+xY8c6A0tNQo01jacpXd2ZnmpIqhmmG1LUENRvLdVf+s2l7r4PNcWGuvkRc6lVYh73h3TKa93Qo/w0gHVy0wDUq3b+oiaPOsK1Sye/jTbayLnm6v76XDEeKn36k0pZFiZ2lEILPy1SnAPQfKLmcICDX/ZjLvnys1qeucGqMvbm7RhrAOIvfvzBGMYhCGAuhaDoVwe57Mev1dLwbpVY5/s/+/Qz+XzceJlmuj4ybd9p/SrzLK1Hsul/jK87SPTFfb9+/ZwRpC/la98xamQuqQGix9t1dy2//PIdBkJtd0qjrqqJUdtZo+aQHlv3r3/9y5k+uktFjRjdXaTmjb74n9I1JXNJy+l8rbtZnn76aVeNbj5QU6TrCVdd22hkLunRbHpSVqOrfseOGndq+KhBpDty9Ag7NdDUQNptt906DBWtS+/94x//6IwfNXXUmFtggQWcsbTlllt2aFNr9/XXX5ezzz7bmVf63SnlpBrss88+nY7vq7370HuVge7qUVNJv1u03377OSOrdumOHjW0Ro8e7YxDNQW1D3ocoZqRXY/g09PA9BjCWp91d5geSbjvvvu6NmpXTaNGzHTHXG0nkcbXgw8+KHpamZpCults4YUXdu2rAVl0dWcuFbWvOeDrD2i/MJeK1An8u49TrudZzjDDDC7AddJRZ1eTSRNcP95Wu3TC1CRUV3Tdddd17qd+h0m3Ampg6TZGTRCfK8ZDpU9/UinLwsSOUmjhp0WKcwCaT9QcDnDwy37MJV9+VsszN1hVxt68HWMNQPzFjz8YwzgEAcylEBT96iCX/fi1WhrerRLjfghAIDcCmEvdKK7na+pHu9TJHj9+/GR3qKM9ZMiQUrHi45TrR7t0e526tOqaLrroos7x1Q/Edb3U5bzooovcljvd/qaX7mzS80XXWWedUn2vLxTjodK7UwlUwMLEjkho4adFinMAmmOq1Ec98VB+DvDJf7iX5x67JNrEJly+fmva+MwBjShYG2N5teyWhHF8bXJgjLkUP46KWsghzooYVPk7vKukTVsQgECKBDCX6lSbMGGCHHnkkW6Hz9dff+223unfalftf+u/dVdQzleMh8oceLIwsaMyWvhpkeIcgOaYS5hLfnlfK+2T/+RhGA1i1II2MaiGqdOaNj5zAOZSmJgoU4u1OCozButlcmCMudT+KMwhztpPeVIP4G1JDfoCAQhYJIC5VKeK7gzSsy/1PMOtt97afWRrxx13dN8p0g+PnXfeee7jWvoRrdqZiBZFraJPMR4qq+h3u9tgYdJuBVgkhlIgxTmA/MNcwlwKMwP45D95GEaDGLWgTQyqYeq0po3PHIC5FCYmytRiLY7KjMF6mRwYYy61PwpziLP2U+a9gSUN6AsEIGCbAOZSnT6bbLKJ+3CYfoBNrwEDBsjQoUPdP3q9/PLLssUWW7gPf6nplPMV46EyB54sBO2ojBZ+WqQ4B6A55hLmkl/e10r75D95GEaDGLWgTQyqYeq0po3PHIC5FCYmytRiLY7KjMF6mRwYYy61PwpziLP2U8ZcsqQBfYEABGwTwFyq02fJJZd0u5YOO+ww99eFF15Y9thjDznggAM67jrwwAPdkXi33nqrbWUj9y7GQ2XkLpuonoWgCRlcJ9DCT4sU5wA0n6g5HODgl/0iPvlP/PnSj1cebeKx9a3ZmjY+c0AjFtbG6KuZxfIwjq9KDowxl+LHUVELOcRZEYMqf4d3lbRpCwIQSJEA5lKdassvv7zbmTRs2DD3V/3f6667rhx99NEddx1//PFyxRVXyFNPPZWi3sH6HOOhMljnDFfEwsSOOGjhp0WKcwCaY6rURz3xUH4O8Ml/uJfnHrsk2sQmXL5+a9r4zAGYS+XjwLektTjyHY/F8jkwxlxqf+TlEGftpzypB/C2pAZ9gQAELBLAXKpTRb+xpN9SOv30091fd9hhB3njjTfktttukz59+ri/qfn00UcfyV133WVRz8r6FOOhsrLOt7EhFiZthN+labTw0yLFOQDNMZcwl/zyvlbaJ//JwzAaxKgFbWJQDVOnNW185oBGRKyNMYxytmqBcXw9cmCMuRQ/jopayCHOihhU+Tu8q6RNWxCAQIoEMJfqVDvxxBPlmmuukb/+9a/OTLr55pvloIMOkkUWWURWXHFFefLJJ+WJJ56QPffc0313KecrxkNlDjxZmNhRGS38tEhxDkBzzCXMJb+8x1wKw89qLcyRVpWxd5xpjDUA8Rc//mAM4xAEMJdCUPSrg1z249dqaXi3Soz7IQCB3AhgLtUpPmbMGHnwwQfdUXizzjqr++Xcc8+V888/Xz799FOZZpppZNCgQe7YvKmnnjq3WOk03hgPlTkAZWFiR2W08NMixTkAzTGXMJf88h5zKQw/q7UwR1pVBnPJrjJp9Ywcj69XDowxl+LHUVELOcRZEYMqf4d3lbRpCwIQSJEA5lITqn311Vfy/vvvyyyzzCK9e/duokTPvyXFF8sWVGFhYkEFXrCHUCHFOYD8I/Yxl0Jkv4hP/pOHYTSIUQvaxKAapk5r2vjMAY2IWBtjGOVs1QLj+HrkwBhzKX4cFbWQQ5wVMajyd3hXSZu2IACBFAlgLtWpduihh8qAAQNkxx13TFHLSvsc46Gy0gG0qTEWJm0C302zaOGnRYpzAJpjLmEu+eV9rbRP/pOHYTSIUQvaxKAapk5r2vjMAZhLYWKiTC3W4qjMGKyXyYEx5lL7ozCHOGs/5Uk9gLclNegLBCBgkQDmUp0qSy65pAwePFgOPPBAi1qZ6lOMh0pTA4zUGRYmkcCWqBYtSkCrK5LiHIDmmEuYS355j7kUhp/VWpgjrSrDsXh2lUmrZ+R4fL1yYIy5FD+OilrIIc6KGFT5O7yrpE1bEIBAigQwl+pU23zzzeX73/++nHTSSSlqWWmfU3yxXCmgBo2xMLGgAi/YQ6iQ4hxA/hH7mEshsp9j8cJQtFcLc6Q9TWo9sqZNjDWAtTHajYbyPYNxeXbNlsyBMeZSs9EQ774c4iwevdZrhnfrzCgBAQjkRQBzqU7vW2+9VfRovMsuu0yWWGKJvCKhxdHGeKhssQtJ3s7CxI5saOGnRYpzAJpjLmEu+eV9rbRP/pOHYTSIUQvaxKAapk5r2vjMAY2IWBtjGOVs1QLj+HrkwBhzKX4cFbWQQ5wVMajyd3hXSZu2IACBFAlgLtWpNnLkSLn55pvlkUcekfXXX18WXnhhmXXWWaVXr16TaTtw4MAU9Q7W5xgPlcE6Z7giFiZ2xEELPy1SnAPQHHMJc8kv7zGXwvCzWgtzpFVlOBbPrjJp9Ywcj69XDowxl+LHUVELOcRZEYMqf4d3lbRpCwIQSJFA9uaS7lT68Y9/LOuss44MGDDAGUkTJkzopGW9uaS/6f8eNWpUinoH63OKL5aDDd6jIhYmHvACF0ULP6ApzgFojrmEueSX95hLYfhZrYU50qoymEt2lUmrZ+R4fL1yYIy5FD+OilrIIc6KGFT5O7yrpE1bEIBAigSyN5fUUBo6dKj7Z8SIEd3uUupO2M022yxFvYP1OcUXy8EG71ERCxMPeIGLooUf0BTnADTHXMJc8st7zKUw/KzWwhxpVRnMJbvKpNUzcjy+XjkwxlyKH0dFLeQQZ0UMqvwd3lXSpi0IQCBFAphLdeZSigK2q88pvlhuFyteZlogP3kfWCT66ZLiHIDmmEvMx355j7kUhp/VWpgjrSqDuWRXmbR6Ro7H1ysHxphL8eOoqIUc4qyIQZW/w7tK2rQFAQikSABzCXOpVNym+GK51EADF2JhEhioR3Vo4QFPRFKcA9AccwlzyS/vMZfC8LNaC3OkVWUwl+wqk1bPyPH4euXAGHMpfhwVtZBDnBUxqPJ3eFdJm7YgAIEUCWAuYS6VitsUXyyXGmjgQixMAgP1qA4tPOBhLvnBa3NpYh+TzTcEfdYAxJ8v/Xjl0SYeW9+arWnjMwc0YmFtjL6aWSwP4/iq5MAYcyl+HBW1kEOcFTGo8nd4V0mbtiAAgRQJYC4NGCD9+vVz/zR79erVSy655JJmb++R98V4qOyRoLoMioWJHZXRwk+LFOcANMdUqY964qH8HOCT/3Avzz12SbSJTbh8/da08ZkDMJfKx4FvSWtx5Dsei+VzYIy51P7IyyHO2k95Ug/gbUkN+gIBCFgkgLk0YEDLuqi5NGrUqJbL9aQCMR4qexIfHpztq8ki0U+jFOcANMdcwlzyy/taaZ/8Jw/DaBCjFrSJQTVMnda08ZkDWCOHiYkytViLozJjsF4mB8aYS+2PwhzirP2UMZcsaUBfIAAB2wQwlwYMkB133FEGDx7cklKt7HRqqeJEbo7xUJnI0L26yULQC1/QwmjhhzPFOQDNMZcwl/zyHnMpDD+rtTBHWlWGby7ZVSatnpHj8fXKgTHmUvw4KmohhzgrYlDl7/CukjZtQQACKRLAXOKbS6XiNsUXy6UGGrgQC5PAQD2qQwsPeHxzyQ9em0sT+5hsviHoswYg/nzpxyuPNvHY+tZsTRufOaARC2tj9NXMYnkYx1clB8aYS/HjqKiFHOKsiEGVv8O7Stq0BQEIpEgAcwlzqVTcxnioLNWRxAqxMLEjGFr4aZHiHIDmmCr1UU88lJ8DfPIf7uW5xy6JNrEJl6/fmjY+cwDmUvk48C1pLY58x2OxfA6MMZfaH3k5xFn7KU/qAbwtqUFfIAABiwQwlzCXSsVljIfKUh1JrBALEzuCoYWfFinOAWiOuYS55Jf3tdI++U8ehtEgRi1oE4NqmDqtaeMzB2AuhYmJMrVYi6MyY7BeJgfGmEvtj8Ic4qz9lDGXLGlAXyAAAdsEMJcwl0pFaIyHylIdSawQC0E7gqGFnxYpzgFojrmEueSX95hLYfhZrYU50qoyfHPJrjJp9Ywcj69XDowxl+LHUVELOcRZEYMqf4d3lbRpCwIQSJFA9uZSiqJZ6HOKL5YtcGNhYkEFXrCHUCHFOYD8I/Yxl0Jkv4hP/pOHYTSIUQvaxKAapk5r2vjMAY2IWBtjGOVs1QLj+HrkwBhzKX4cFbWQQ5wVMajyd3hXSZu2IACBFAlgLqWomoE+x3ioNDCs6F1gYRIdcdMNoEXTqLq9McU5AM0xlzCX/PK+Vton/8nDMBrEqAVtYlANU6c1bXzmAMylMDFRphZrcVRmDNbL5MAYc6n9UZhDnLWf8qQewNuSGvQFAhCwSABzyaIqCfQpxkNlAsP27iILE2+EwSpACz+UKc4BaI65hLnkl/eYS2H4Wa2FOdKqMhyLZ1eZtHpGjsfXKwfGmEvx46iohRzirIhBlb/Du0ratAUBCKRIAHMpRdUM9DnFF8sGsAkLEwsq8II9hAopzgHkH7GPuRQi+zkWLwxFe7UwR9rTpNYja9rEWANYG6PdaCjfMxiXZ9dsyRwYYy41Gw3x7sshzuLRa71meLfOjBIQgEBeBDCX8tI72GhjPFQG65zhiliY2BEHLfy0SHEOQHPMJcwlv7yvlfbJf/IwjAYxakGbGFTD1GlNG585oBERa2MMo5ytWmAcX48cGGMuxY+johZyiLMiBlX+Du8qadMWBCCQIgHMpRRVM9DnGA+VBoYVvQssTKIjbroBtGgaVbc3pjgHoDnmEuaSX95jLoXhZ7UW5kirynAsnl1l0uoZOR5frxwYYy7Fj6OiFnKIsyIGVf4O7ypp0xYEIJAiAcylFFUz0OcUXywbwMaxeBZE+KYPLBL9xEhxDkBzzCXMJb+8x1wKw89qLcyRVpXBXLKrTFo9I8fj65UDY8yl+HFU1EIOcVbEoMrf4V0lbdqCAARSJIC5lKJqBvqc4otlA9gwlyyIgLkURIUU5wAeDDCXMJeCpL/45D95GEaDGLWgTQyqYeq0po3PHNCIiLUxhlHOVi0wjq9HDowxl+LHUVELOcRZEYMqf4d3lbRpCwIQSJEA5lKKqhnoc4yHSgPDit4FFibRETfdAFo0jarbG1OcA9AccwlzyS/va6V98p88DKNBjFrQJgbVMHVa08ZnDsBcChMTZWqxFkdlxmC9TA6MMZfaH4U5xFn7KU/qAbwtqUFfIAABiwQwlyyqkkCfYjxUJjBs7y6yMPFGGKwCtPBDmeIcgOaYS5hLfnmPuRSGn9VamCOtKsOxeGupc0UAACAASURBVHaVSatn5Hh8vXJgjLkUP46KWsghzooYVPk7vKukTVsQgECKBDCXUlTNQJ9TfLFsABvH4lkQ4Zs+sEj0EyPFOQDNMZcwl/zyHnMpDD+rtTBHWlUGc8muMmn1jByPr1cOjDGX4sdRUQs5xFkRgyp/h3eVtGkLAhBIkQDmUoqqGehzii+WDWDDXLIgAuZSEBVSnAN4MMBcwlwKkv58cykMRnO1MEeak6SjQ9a0ibEGsDZGu9FQvmcwLs+u2ZI5MMZcajYa4t2XQ5zFo9d6zfBunRklIACBvAhgLuWld7DRxnioDNY5wxWxMLEjDlr4aZHiHIDmmEuYS355Xyvtk//kYRgNYtSCNjGohqnTmjY+c0AjItbGGEY5W7XAOL4eOTDGXIofR0Ut5BBnRQyq/B3eVdKmLQhAIEUCmEspqmagzzEeKg0MK3oXWJhER9x0A2jRNKpub0xxDkBzzCXMJb+8x1wKw89qLcyRVpXhWDy7yqTVM3I8vl45MMZcih9HRS3kEGdFDKr8Hd5V0qYtCEAgRQKYSymqZqDPKb5YNoCNY/EsiPBNH1gk+omR4hyA5phLmEt+eY+5FIaf1VqYI60qg7lkV5m0ekaOx9crB8aYS/HjqKiFHOKsiEGVv8O7Stq0BQEIpEgAcylF1Qz0OcUXywawYS5ZEAFzKYgKKc4BPBhgLmEuBUl/vrkUBqO5WpgjzUnS0SFr2sRYA1gbo91oKN8zGJdn12zJHBhjLjUbDfHuyyHO4tFrvWZ4t86MEhCAQF4EMJcq1vvLL7+UP/zhD3L99dfL22+/Lf369ZPtt99etttuO+nVq9cUe/PVV1/JhRde6MqOGTNGZpppJll33XXlgAMOkBlmmGGysnr/VVddJdddd528/PLL0qdPH+nfv78MGTJEVlttNa+Rx3io9OpQIoVZmNgRCi38tEhxDkBzzCXMJb+8r5X2yX/yMIwGMWpBmxhUw9RpTRufOaAREWtjDKOcrVpgHF+PHBhjLsWPo6IWcoizIgZV/g7vKmnTFgQgkCIBzKWKVTv88MPl2muvlZ/97GeyxBJLyF//+le5/fbbZZ999pGhQ4dOsTcHHXSQ3HTTTbLhhhvKCiusIK+99ppcccUVsuCCCzoTaeqpp+4o//XXX8u+++4rDzzwgGy22Way+OKLy7hx4+TFF190//egQYO8Rh7jodKrQ4kUZmFiRyi08NMixTkAzTGXMJf88h5zKQw/q7UwR1pVhmPx7CqTVs/I8fh65cAYcyl+HBW1kEOcFTGo8nd4V0mbtiAAgRQJYC5VqNqoUaNk4MCBsssuu8iwYcM6Wt5///3lnnvucf/MMccc3fboX//6l2yxxRay9dZby29/+9uOe+68805nTP3mN7+RbbbZpuPvl156qfz+97+XSy65RJZbbrngo0zxxXJwCCUqZGFSAlqkImjhBzbFOQDNMZcwl/zyHnMpDD+rtTBHWlUGc8muMmn1jByPr1cOjDGX4sdRUQs5xFkRgyp/h3eVtGkLAhBIkQDmUoWqnXzyye5IvPvuu0/mmWeejpYff/xx2XbbbeWII45w/+7uuvjii+XYY4+VK6+8UpZddtlOtyy99NKy8MILu9/00l1L66yzjtuhdPrpp7v/rbuW+vbtG2y0Kb5YDjZ4j4pYmHjAC1wULfyApjgHoDnmEuaSX95jLoXhZ7UW5kirymAu2VUmrZ6R4/H1yoEx5lL8OCpqIYc4K2JQ5e/wrpI2bUEAAikSwFyqUDXdsfT888+7o/Dqr/Hjx8uSSy4pm2++uRx99NHd9ui8886Tk046SW644QZZZJFFOt2z0korOfPoiSeecN9t0qPvNtpoI/ctprFjx8qIESPc73PPPbfsueeebveT75Xii2XfMYcoz8IkBMUwdaCFH8cU5wA0x1zCXPLLe8ylMPys1sIcaVUZzCW7yqTVM3I8vl45MMZcih9HRS3kEGdFDKr8Hd5V0qYtCEAgRQKYSxWqtvHGG0ufPn2c2dP1UoNo0UUXlQsuuKDbHt19990yZMgQOfTQQ2WnnXbquOeFF14QrVevRx55RGaaaSap3TvzzDPLdNNN5wwl3bX0pz/9SR599FH59a9/Ldttt53XyHVROWHChKC7obw6lEhhNfn0Ul242kugp2gxYMCAtoBMcQ7oKZr7Cg6HiQR7Cod2zAE++d9TuPvmocXyaGNRlSnPV+3If+2RzxzQiDLxFz/+YNzzGLdjDuia/8RV/Ljq2gLMq2VulXc78r9a8rQGAQikQgBzqUKlfvzjH8tss80mV1999WStrrnmmjLffPPJZZdd1m2PdHfTT37yE3n//ffl8MMPlx/96Efy+uuvy1FHHSWvvvqqfPHFF/LAAw/IXHPNJTfeeKMcfPDBMvXUU8ttt93m6tXryy+/dEbUBx984HZPTTXVVKVHH+OhsnRnEipodWGSEMJgXe0pWrRrUZniHNBTNPdNAjhgLvnGkE/+E3++9OOVR5t4bH1rbqRNT1oDEH++UVJcHsbFjHzvqJpxO+YAzCXfKPEvX3Wc+fc47Rqs8m5H/qetJL2HAARiEcBcikW2m3p9di5pdS+//LIceOCB8swzz7ja9Qi8zTbbTD7++GO566675O9//7tMP/30cscdd8i+++4ryy+//GRm1RlnnCFnnnmm3HTTTbLQQguVHn2KR2KVHmzAgmypDgjTsyq08AOY4hyA5hM1hwMc/LJ/4q4FvZZZZpmWqyL+WkZWWQG0qQx1yw1Z08ZnDmg0eGtjbFmkBArAOL5IOTDmWLz4cVTUQg5xVsSgyt/hXSVt2oIABFIkgLlUoWpF31xSo+iYY44p7NErr7wib7/9ttuRpDuVfvazn8mYMWPkoYcecmX120v6XSXd6XTKKad0qu+qq66S3/zmN3L55Ze73U9lrxgPlWX7klI5FiZ21EILPy1SnAPQfKLmcICDX/ZjLvnys1qeucGqMvbm7RhrAOIvfvzBGMYhCGAuhaDoVwe57Mev1dLwbpUY90MAArkRwFyqUPGTTjpJzjvvPLnvvvtknnnm6Wj58ccfl2233bbUt5D0mLzVVltNNtxwQznhhBNcnZ9++qmsuOKKsthii4maSfWXmk3nnnuu3HrrrdK/f//So4/xUFm6MwkVZGFiRyy08NMixTkAzTFV6qOeeCg/B/jkP9zLc49dEm1iEy5fvzVtfOaARhSsjbG8WnZLwji+NjkwxlyKH0dFLeQQZ0UMqvwd3lXSpi0IQCBFAphLFar27LPPumPsdAfTsGHDOlref//95e6775Z77rlH5pxzTveR8TfeeENmnnlmmWWWWabYw0MPPdR9Y+m6666TRRZZpONePRZPj8q74YYbpHYWq9aru5n0OD1tS/9d9orxUFm2LymVY2FiRy208NMixTkAzTGXMJf88r5W2if/ycMwGsSoBW1iUA1TpzVtfOYAzKUwMVGmFmtxVGYM1svkwBhzqf1RmEOctZ/ypB7A25Ia9AUCELBIAHOpYlWGDx8uI0aMcEfZLb744u4ou9tuu02GDh0q++yzj+vNI488IoMHD+70N/27fm9phhlmkAUWWEC++uor922lxx57TA455BDZeeedO43k1VdflUGDBjkDSevq27eva/eFF16Q008/XdZdd12vkcd4qPTqUCKFWZjYEQot/LRIcQ5Ac8wlzCW/vMdcCsPPai3MkVaV4Vg8u8qk1TNyPL5eOTDGXIofR0Ut5BBnRQyq/B3eVdKmLQhAIEUCmEsVq/bFF1/IH/7wB2f0jB07Vvr16yfbbbed7LDDDh07iRqZSxdffLFcf/318vrrr0vv3r1l0UUXld12203WWGONbkcxevRoOfHEE50BNX78eLezaciQIe4YPd8rxRfLvmMOUZ6FSQiKYepACz+OKc4BaI65hLnkl/eYS2H4Wa2FOdKqMphLdpVJq2fkeHy9cmCMuRQ/jopayCHOihhU+Tu8q6RNWxCAQIoEMJdSVM1An1N8sWwAm7AwsaACL9hDqJDiHED+EfuYSyGyX8Qn/8nDMBrEqAVtYlANU6c1bXzmgEZErI0xjHK2aoFxfD1yYIy5FD+OilrIIc6KGFT5O7yrpE1bEIBAigQwl1JUzUCfYzxUGhhW9C6wMImOuOkG0KJpVN3emOIcgOaYS5hLfnlfK+2T/+RhGA1i1II2MaiGqdOaNj5zAOZSmJgoU4u1OCozButlcmCMudT+KMwhztpPeVIP4G1JDfoCAQhYJIC5ZFGVBPoU46EygWF7d5GFiTfCYBWghR/KFOcANMdcwlzyy3vMpTD8rNbCHGlVGY7Fs6tMWj0jx+PrlQNjzKX4cVTUQg5xVsSgyt/hXSVt2oIABFIkgLmUomoG+pzii2UD2DgWz4II3/SBRaKfGCnOAWiOuYS55Jf3mEth+FmthTnSqjKYS3aVSatn5Hh8vXJgjLkUP46KWsghzooYVPk7vKukTVsQgECKBDCXUlTNQJ9TfLFsABvmkgURMJeCqJDiHMCDAeYS5lKQ9OebS2EwmquFOdKcJB0dsqZNjDWAtTHajYbyPYNxeXbNlsyBMeZSs9EQ774c4iwevdZrhnfrzCgBAQjkRQBzKS+9g402xkNlsM4ZroiFiR1x0MJPixTnADTHXMJc8sv7Wmmf/CcPw2gQoxa0iUE1TJ3WtPGZAxoRsTbGMMrZqgXG8fXIgTHmUvw4KmohhzgrYlDl7/CukjZtQQACKRLAXEpRNQN9jvFQaWBY0bvAwiQ64qYbQIumUXV7Y4pzAJpjLmEu+eU95lIYflZrYY60qgzH4tlVJq2ekePx9cqBMeZS/DgqaiGHOCtiUOXv8K6SNm1BAAIpEsBcSlE1A31O8cWyAWwci2dBhG/6wCLRT4wU5wA0x1zCXPLLe8ylMPys1sIcaVUZzCW7yqTVM3I8vl45MMZcih9HRS3kEGdFDKr8Hd5V0qYtCEAgRQKYSymqZqDPKb5YNoANc8mCCJhLQVRIcQ7gwQBzCXMpSPrzzaUwGM3VwhxpTpKODlnTJsYawNoY7UZD+Z7BuDy7ZkvmwBhzqdloiHdfDnEWj17rNcO7dWaUgAAE8iKAuZSX3sFGG+OhMljnDFfEwsSOOGjhp0WKcwCaYy5hLvnlfa20T/6Th2E0iFEL2sSgGqZOa9r4zAGNiFgbYxjlbNUC4/h65MAYcyl+HBW1kEOcFTGo8nd4V0mbtiAAgRQJYC6lqJqBPsd4qDQwrOhdYGESHXHTDaBF06i6vTHFOQDNMZcwl/zyHnMpDD+rtTBHWlWGY/HsKpNWz8jx+HrlwBhzKX4cFbWQQ5wVMajyd3hXSZu2IACBFAlgLqWomoE+p/hi2QA2jsWzIMI3fWCR6CdGinMAmmMuYS755T3mUhh+VmthjrSqDOaSXWXS6hk5Hl+vHBhjLsWPo6IWcoizIgZV/g7vKmnTFgQgkCIBzKUUVTPQ5xRfLBvAhrlkQQTMpSAqpDgH8GCAuYS5FCT9+eZSGIzmamGONCdJR4esaRNjDWBtjHajoXzPYFyeXbMlc2CMudRsNMS7L4c4i0ev9Zrh3TozSkAAAnkRwFzKS+9go43xUBmsc4YrYmFiRxy08NMixTkAzTGXMJf88r5W2if/ycMwGsSoBW1iUA1TpzVtfOaARkSsjTGMcrZqgXF8PXJgjLkUP46KWsghzooYVPk7vKukTVsQgECKBDCXUlTNQJ9jPFQaGFb0LrAwiY646QbQomlU3d6Y4hyA5phLmEt+eY+5FIaf1VqYI60qw7F4dpVJq2fkeHy9cmCMuRQ/jopayCHOihhU+Tu8q6RNWxCAQIoEMJdSVM1An1N8sWwAG8fiWRDhmz6wSPQTI8U5AM0xlzCX/PIecykMP6u1MEdaVQZzya4yafWMHI+vVw6MMZfix1FRCznEWRGDKn+Hd5W0aQsCEEiRAOZSiqoZ6HOKL5YNYMNcsiAC5lIQFVKcA3gwwFzCXAqS/nxzKQxGc7UwR5qTpKND1rSJsQawNka70VC+ZzAuz67ZkjkwxlxqNhri3ZdDnMWj13rN8G6dGSUgAIG8CGAu5aV3sNHGeKgM1jnDFbEwsSMOWvhpkeIcgOaYS5hLfnlfK+2T/+RhGA1i1II2MaiGqdOaNj5zQCMi1sYYRjlbtcA4vh45MMZcih9HRS3kEGdFDKr8Hd5V0qYtCEAgRQKYSymqZqDPMR4qDQwrehdYmERH3HQDaNE0qm5vTHEOQHPMJcwlv7zHXArDz2otzJFWleFYPLvKpNUzcjy+XjkwxlyKH0dFLeQQZ0UMqvwd3lXSpi0IQCBFAphLKapmoM8pvlg2gI1j8SyI8E0fWCT6iZHiHIDmmEuYS355j7kUhp/VWpgjrSqDuWRXmbR6Ro7H1ysHxphL8eOoqIUc4qyIQZW/w7tK2rQFAQikSABzKUXVDPQ5xRfLBrBhLlkQAXMpiAopzgE8GGAuYS4FSX++uRQGo7lamCPNSdLRIWvaxFgDWBuj3Wgo3zMYl2fXbMkcGGMuNRsN8e7LIc7i0Wu9Zni3zowSEIBAXgQwl/LSO9hoYzxUBuuc4YpYmNgRBy38tEhxDkBzzCXMJb+8r5X2yX/yMIwGMWpBmxhUw9RpTRufOaAREWtjDKOcrVpgHF+PHBhjLsWPo6IWcoizIgZV/g7vKmnTFgQgkCIBzKUUVTPQ5xgPlQaGFb0LLEyiI266AbRoGlW3N6Y4B6A55hLmkl/eYy6F4We1FuZIq8pwLJ5dZdLqGTkeX68cGGMuxY+johZyiLMiBlX+Du8qadMWBCCQIgHMpRRVM9DnFF8sG8DGsXgWRPimDywS/cRIcQ5Ac8wlzCW/vMdcCsPPai3MkVaVwVyyq0xaPSPH4+uVA2PMpfhxVNRCDnFWxKDK3+FdJW3aggAEUiSAuZSiagb6nOKLZQPYMJcsiIC5FESFFOcAHgwwlzCXgqQ/31wKg9FcLcyR5iTp6JA1bWKsAayN0W40lO8ZjMuza7ZkDowxl5qNhnj35RBn8ei1XjO8W2dGCQhAIC8CmEt56R1stDEeKoN1znBFLEzsiIMWflqkOAegOeYS5pJf3tdK++Q/eRhGgxi1oE0MqmHqtKaNzxzQiIi1MYZRzlYtMI6vRw6MMZfix1FRCznEWRGDKn+Hd5W0aQsCEEiRAOZSiqoZ6HOMh0oDw4reBRYm0RE33QBaNI2q2xtTnAPQHHMJc8kv7zGXwvCzWgtzpFVlOBbPrjJp9Ywcj69XDowxl+LHUVELOcRZEYMqf4d3lbRpCwIQSJEA5lKKqhnoc4ovlg1g41g8CyJ80wcWiX5ipDgHoDnmEuaSX95jLoXhZ7UW5kirymAu2VUmrZ6R4/H1yoEx5lL8OCpqIYc4K2JQ5e/wrpI2bUEAAikSwFxKUTUDfU7xxbIBbJhLFkTAXAqiQopzAA8GmEuYS0HSn28uhcForhbmSHOSdHTImjYx1gDWxmg3Gsr3DMbl2TVbMgfGmEvNRkO8+3KIs3j0Wq8Z3q0zowQEIJAXAcylvPQONtoYD5XBOme4IhYmdsRBCz8tUpwD0BxzCXPJL+9rpX3ynzwMo0GMWtAmBtUwdVrTxmcOaETE2hjDKGerFhjH1yMHxphL8eOoqIUc4qyIQZW/w7tK2rQFAQikSABzKUXVDPQ5xkOlgWFF7wILk+iIm24ALZpG1e2NKc4BaI65hLnkl/eYS2H4Wa2FOdKqMhyLZ1eZtHpGjsfXKwfGmEvx46iohRzirIhBlb/Du0ratAUBCKRIAHMpRdUM9DnFF8sGsHEsngURvukDi0Q/MVKcA9AccwlzyS/vMZfC8LNaC3OkVWUwl+wqk1bPyPH4euXAGHMpfhwVtZBDnBUxqPJ3eFdJm7YgAIEUCWAupaiagT6n+GLZADbMJQsiYC4FUSHFOYAHA8wlzKUg6c83l8JgNFcLc6Q5STo6ZE2bGGsAa2O0Gw3lewbj8uyaLZkDY8ylZqMh3n05xFk8eq3XDO/WmVECAhDIiwDmUl56BxttjIfKYJ0zXBELEzvioIWfFinOAWiOuYS55Jf3tdI++U8ehtEgRi1oE4NqmDqtaeMzBzQiYm2MYZSzVQuM4+uRA2PMpfhxVNRCDnFWxKDK3+FdJW3aggAEUiSAuZSiagb6HOOh0sCwoneBhUl0xE03gBZNo+r2xhTnADTHXMJc8st7zKUw/KzWwhxpVRmOxbOrTFo9I8fj65UDY8yl+HFU1EIOcVbEoMrf4V0lbdqCAARSJIC5lKJqBvqc4otlA9g4Fs+CCN/0gUWinxgpzgFojrmEueSX95hLYfhZrYU50qoymEt2lUmrZ+R4fL1yYIy5FD+OilrIIc6KGFT5O7yrpE1bEIBAigQwl1JUzUCfU3yxbAAb5pIFETCXgqiQ4hzAgwHmEuZSkPTnm0thMJqrhTnSnCQdHbKmTYw1gLUx2o2G8j2DcXl2zZbMgTHmUrPREO++HOIsHr3Wa4Z368woAQEI5EUAcykvvYONNsZDZbDOGa6IhYkdcdDCT4sU5wA0x1zCXPLL+1ppn/wnD8NoEKMWtIlBNUyd1rTxmQMaEbE2xjDK2aoFxvH1yIEx5lL8OCpqIYc4K2JQ5e/wrpI2bUEAAikSwFxKUTUDfY7xUGlgWNG7wMIkOuKmG0CLplF1e2OKcwCaYy5hLvnlPeZSGH5Wa2GOtKoMx+LZVSatnpHj8fXKgTHmUvw4KmohhzgrYlDl7/CukjZtQQACKRLAXEpRNQN9fvzxx10vevXqZaA36XRhwoQJcDMiV0/RQnNw6aWXrpxqinNAT9HcV2w4TCTYUzi0Yw7wyf+ewt03Dy2WRxuLqkx5vmpH/muPfOaARpSJv/jxB+Oex7gdc0DX/Ceu4sdV1xZgXi1zq7zbkf/Vkqc1CEAgFQKYS6koZayfMR4qjQ2R7kAgCQLtWlQyByQRHnQyAwLtmAPI/wwCiyEmQaAd+R/LXEoCOJ2EgDEC7ZgDWAMYCwK6ky2BduR/trAZOAQgMEUCmEsECAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAQNMEMJeaRsWNEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACmEvEAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAQNMEMJeaRsWNEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACmEvEAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAQNMEMJeaRsWNEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACmEvEAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAQNMEMJeaRsWNEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACmEvEAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAQNMEMJeaRsWN9QSeeOIJ9z+XXnppwEAAAhkSYA7IUHSGDIFvCJD/hAIE8ibAHJC3/ow+bwLkf976M3oIQAACEIBAVwKYS8REKQL/+Mc/XLllllmmVPlcC/373/92Qx8wYECuCMyMGy38pEhxDkDziZrDAQ5+2S/ik//Eny/9eOXRJh5b35qtaeMzBzRiYW2MvppZLA/j+KrkwLhr/ucw5viR01oLMG+Nl+/d8PYlSHkIQKCnE8Bc6ukKRxpfjIfKSF01VS0LEztyoIWfFinOAWiOqVIf9cRD+TnAJ//hXp577JJoE5tw+fqtaeMzB2AulY8D35LW4sh3PBbL58AYc6n9kZdDnLWf8qQewNuSGvQFAhCwSABzyaIqCfQpxkNlAsP27iILE2+EwSpACz+UKc4BaI65hLnkl/e10j75Tx6G0SBGLWgTg2qYOq1p4zMHYC6FiYkytViLozJjsF4mB8aYS+2PwhzirP2UMZcsaUBfIAAB2wQwl2zrY7Z3MR4qzQ42YMdYCAaE6VkVWvgBTHEOQHPMJcwlv7zHXArDz2otzJFWlbF3nGmMNQDxFz/+YAzjEAQwl0JQ9KuDXPbj12ppeLdKjPshAIHcCGAu5aZ4oPHGeKgM1DXT1bAwsSMPWvhpkeIcgOaYS5hLfnmPuRSGn9VamCOtKoO5ZFeZtHpGjsfXKwfGmEvx46iohRzirIhBlb/Du0ratAUBCKRIAHMpRdUM9DnFF8sGsAkLEwsq8II9hAopzgHkH7GPuRQi+0V88p88DKNBjFrQJgbVMHVa08ZnDmhExNoYwyhnqxYYx9cjB8aYS/HjqKiFHOKsiEGVv8O7Stq0BQEIpEgAcylF1Qz0OcZDpYFhRe8CC5PoiJtuAC2aRtXtjSnOAWiOuYS55Jf3tdI++U8ehtEgRi1oE4NqmDqtaeMzB2AuhYmJMrVYi6MyY7BeJgfGmEvtj8Ic4qz9lCf1AN6W1KAvEICARQKYSxZVSaBPMR4qExi2dxdZmHgjDFYBWvihTHEOQHPMJcwlv7zHXArDz2otzJFWleFYPLvKpNUzcjy+XjkwxlyKH0dFLeQQZ0UMqvwd3lXSpi0IQCBFAphLKapmoM8pvlg2gI1j8SyI8E0fWCT6iZHiHIDmmEuYS355j7kUhp/VWpgjrSqDuWRXmbR6Ro7H1ysHxphL8eOoqIUc4qyIQZW/w7tK2rQFAQikSABzKUXVDPQ5xRfLBrBhLlkQAXMpiAopzgE8GGAuYS4FSX++uRQGo7lamCPNSdLRIWvaxFgDWBuj3Wgo3zMYl2fXbMkcGGMuNRsN8e7LIc7i0Wu9Zni3zowSEIBAXgQwl/LSO9hoYzxUBuuc4YpYmNgRBy38tEhxDkBzzCXMJb+8r5X2yX/yMIwGMWpBmxhUw9RpTRufOaAREWtjDKOcrVpgHF+PHBhjLsWPo6IWcoizIgZV/g7vKmnTFgQgkCIBzKUUVTPQ5xgPXfjO8wAAIABJREFUlQaGFb0LLEyiI266AbRoGlW3N6Y4B6A55hLmkl/eYy6F4We1FuZIq8pwLJ5dZdLqGTkeX68cGGMuxY+johZyiLMiBlX+Du8qadMWBCCQIgHMpRRVM9DnFF8sG8DGsXgWRPimDywS/cRIcQ5Ac8wlzCW/vMdcCsPPai3MkVaVwVyyq0xaPSPH4+uVA2PMpfhxVNRCDnFWxKDK3+FdJW3aggAEUiSAuZSiagb6nOKLZQPYMJcsiIC5FESFFOcAHgwwlzCXgqQ/31wKg9FcLcyR5iTp6JA1bWKsAayN0W40lO8ZjMuza7ZkDowxl5qNhnj3/T97bwJ1R1Xl7e83hEBIGhZ0OvkkCGLbkm4WUwQB6W6QWQQMCBEJJA3SQEMiMzFAi4jEoAQbDSiIQ4hMBhCbuQnyKSIi8lehV/KBC3GAaAPayBAwEPJfVXDfvEnemxr22VX73PPctViRVJ19dj2/vY+n7o+qm0Kd2dGrHhne1ZkxAgIQSIsA5lJaege7WoubymDJOQ7ExsSPOGih0yLGNQDNMZcwl3R93xmt6X/6MIwGFlHQxoJqmJjetNGsAd2IeLvGMMr5igJjez1SYIy5ZF9HRTOkUGdFDJo8Du8maTMXBCAQIwHMpRhVc5CzxU2lg8syT4GNiTni0hOgRWlUg54Y4xqA5phLmEu6vsdcCsPPaxTWSK/K8Fo8v8rElRk9bq9XCowxl+zrqGiGFOqsiEGTx+HdJG3mggAEYiSAuRSjag5yjvGLZQfYeC2eBxHeyoFNok6MGNcANMdcwlzS9T3mUhh+XqOwRnpVBnPJrzJxZUaP2+uVAmPMJfs6KpohhTorYtDkcXg3SZu5IACBGAlgLsWomoOcY/xi2QE2zCUPImAuBVEhxjWAGwPMJcylIO3Pby6FweguCmukO0n6E/KmjcUewNs1+q2G+pnBuD67siNTYIy5VLYa7M5Loc7s6FWPDO/qzBgBAQikRQBzKS29g12txU1lsOQcB2Jj4kcctNBpEeMagOaYS5hLur7vjNb0P30YRgOLKGhjQTVMTG/aaNaAbkS8XWMY5XxFgbG9Hikwxlyyr6OiGVKosyIGTR6Hd5O0mQsCEIiRAOZSjKo5yNniptLBZZmnwMbEHHHpCdCiNKpBT4xxDUBzzCXMJV3fYy6F4ec1CmukV2V4LZ5fZeLKjB631ysFxphL9nVUNEMKdVbEoMnj8G6SNnNBAAIxEsBcilE1BznH+MWyA2y8Fs+DCG/lwCZRJ0aMawCaYy5hLun6HnMpDD+vUVgjvSqDueRXmbgyo8ft9UqBMeaSfR0VzZBCnRUxaPI4vJukzVwQgECMBDCXYlTNQc4xfrHsABvmkgcRMJeCqBDjGsCNAeYS5lKQ9uc3l8JgdBeFNdKdJP0JedPGYg/g7Rr9VkP9zGBcn13ZkSkwxlwqWw1256VQZ3b0qkeGd3VmjIAABNIigLmUlt7BrtbipjJYco4DsTHxIw5a6LSIcQ1Ac8wlzCVd33dGa/qfPgyjgUUUtLGgGiamN200a0A3It6uMYxyvqLA2F6PFBhjLtnXUdEMKdRZEYMmj8O7SdrMBQEIxEgAcylG1RzkbHFT6eCyzFNgY2KOuPQEaFEa1aAnxrgGoDnmEuaSru8xl8Lw8xqFNdKrMrwWz68ycWVGj9vrlQJjzCX7OiqaIYU6K2LQ5HF4N0mbuSAAgRgJYC7FqJqDnGP8YtkBNl6L50GEt3Jgk6gTI8Y1AM0xlzCXdH2PuRSGn9corJFelcFc8qtMXJnR4/Z6pcAYc8m+jopmSKHOihg0eRzeTdJmLghAIEYCmEsxquYg5xi/WHaADXPJgwiYS0FUiHEN4MYAcwlzKUj785tLYTC6i8Ia6U6S/oS8aWOxB/B2jX6roX5mMK7PruzIFBhjLpWtBrvzUqgzO3rVI8O7OjNGQAACaRHAXEpL72BXa3FTGSw5x4HYmPgRBy10WsS4BqA55hLmkq7vO6M1/U8fhtHAIgraWFANE9ObNpo1oBsRb9cYRjlfUWBsr0cKjDGX7OuoaIYU6qyIQZPH4d0kbeaCAARiJIC5FKNqDnK2uKl0cFnmKbAxMUdcegK0KI1q0BNjXAPQHHMJc0nX95hLYfh5jcIa6VUZXovnV5m4MqPH7fVKgTHmkn0dFc2QQp0VMWjyOLybpM1cEIBAjAQwl2JUzUHOMX6x7AAbr8XzIMJbObBJ1IkR4xqA5phLmEu6vsdcCsPPaxTWSK/KYC75VSauzOhxe71SYIy5ZF9HRTOkUGdFDJo8Du8maTMXBCAQIwHMpRhVc5BzjF8sO8CGueRBBMylICrEuAZwY4C5hLkUpP35zaUwGN1FYY10J0l/Qt60sdgDeLtGv9VQPzMY12dXdmQKjDGXylaD3Xkp1JkdveqR4V2dGSMgAIG0CGAuBdD79ddfl8svv1xuvPFGefbZZ2Xs2LFyxBFHyKRJk6Svr6/rDM8//7x85zvfkXvvvVeeeOIJWbJkibz97W+X/fffX6ZMmSLrrLNO17FLly6VAw44QH7961/L8ccfL6eccspq595+++1yxRVX5LE32GAD+cAHPiAnn3yyjBgxQn3VFjeV6qQiCMDGxI9IaKHTIsY1AM0xlzCXdH3fGa3pf/owjAYWUdDGgmqYmN600awB3Yh4u8YwyvmKAmN7PVJgjLlkX0dFM6RQZ0UMmjwO7yZpMxcEIBAjAcylAKqdc845Mn/+fJk4caJsvfXW8sMf/lDuvPNOmTZtmkydOrXrDJmpdOKJJ8ouu+wiO+20k4wcOVIeeughufXWW2X8+PEyb948WWuttQYdf+mll8qVV16ZG1KDmUv/+Z//KWeccUYe94Mf/GBuQl111VWy/fbbyze+8Y01ml5lkFjcVJaZN/Zz2Jj4URAtdFrEuAagOeYS5pKu7zGXwvDzGoU10qsyvBbPrzJxZUaP2+uVAmPMJfs6KpohhTorYtDkcXg3SZu5IACBGAlgLilVW7RokUyYMEGOPvpomT59en+07Amhe+65J/9n9OjRg87yu9/9Lv/77GmlgZ9LLrlELrvsMpkzZ47stddeq43NxmVPN2XG1OzZs1czl7Knmt7//vfLxhtvLNddd12/QXXNNdfIeeed1zVuFRQxfrFc5fqszmVjYkW2ely0qM5s4IgY1wA0f1NBOMBB1/3Ca/G0AJ2OZ21wKozDddtiD0D92dcfjGEcggDmUgiKuhj0so5f1dHwrkqM8yEAgdQIYC4pFb/44ovzV+JlTyFlZk7n8/DDD8vhhx8u5557bv5nlc9jjz0mBx54oJx00klywgknrDb0uOOOk5dffllmzZole+yxx2rmUvbk1Mc+9jG58MILc+Or88lMpx133FF23XVX+Y//+I8qKa12rsVNpSqhSAazMfEjFFrotIhxDUBzTJWBVU891F8DNP0P9/rcrUeijTXh+vG9aaNZA7pR8HaN9dXyOxLG9tqkwBhzyb6OimZIoc6KGDR5HN5N0mYuCEAgRgKYS0rVsieWHn/88fxVeAM/mZGzzTbbyMEHHywXXHBBpVnuu+8+OeaYY/KnjA477LCVxi5YsEA+/vGP57/VlP120mDm0le+8hX5whe+IHfccYe8853vXGn8Rz/6UXnuuefk7rvvrpTTqidb3FSqEopkMBsTP0KhhU6LGNcANMdcwlzS9X1ntKb/6cMwGlhEQRsLqmFietNGswZ0I+LtGsMo5ysKjO31SIEx5pJ9HRXNkEKdFTFo8ji8m6TNXBCAQIwEMJeUqmWvpxs2bJjcdNNNq0XaeeedZcstt8x/G6ns54033pApU6bIo48+KpmRNGrUqP6hr7zySv77Sbvvvrtkv/P01FNPDWouffrTn5arr75asqenst9xGvjJnobKnrJ65JFHyqY06HnZpnL58uW5wcWnPIFMw+wzfPjw8oM404RAr2gxbtw4Ez5FQWNcA3pF8yJtio7D4U1CvcKhjTVA0/+9wr2oz2I8jjZ+VeumTRv9n1HSrAHdKFN/9vUH495j3MYasGr/U1f2dbXqDDBvlrlX3m30f7PkmQ0CEIiFAOaSUqk999wzN4Cy3zZa9bPbbrvlv6c0b9680rN0XrOXmUdHHnnkSuOyY/Pnz5e77rpL1l9//a7m0llnnSU33nhjblBlxtfAz5lnninf/e5389/c6OvrK53Xqida3FTWTiaigV43JhEhDJZqr2jR1qYyxjWgVzTXNgEcMJe0NaTpf+pPS99uPNrYsdVGxlzSEmR8RoAet6+Dphm3cR+AuWRfR0UzNF1nRfn0+nGvvNvo/17XmuuDAATqEcBcqsetf1TIJ5e+9a1vyfnnn5+/Ci97Jd7Az69+9av8d5iy33A69NBD80NtP7mU5TB+/HglwbSG80i1H73RQqeFxStxdBkVj0bzNxnBAQ7F3bLmMzT9T/1p6duNRxs7ttrI3rTRrAHdWHi7Rq1mHsfD2F6VFBjzWjz7OiqaIYU6K2LQ5HF4N0mbuSAAgRgJYC4pVSv6zaWDDjpIZs6cWThL9lq97Imj/fbbTy666CIZMmTISmP+7d/+TTKDKXvFXueJoz/84Q8yadIkOeKII+Soo47Kn6Bad911peg3l5599tn8lXuaj8VNpSafWMayMfGjFFrotIhxDUBzTJWBVU891F8DNP0P9/rcrUeijTXh+vG9aaNZAzCX6teBdqS3OtJej8fxKTDGXGq/8lKos/Ypr8gA3p7UIBcIQMAjAcwlpSqzZ8+WK664Iv8do4033rg/WvZ7R4cffrh88pOfzA2gNX1uu+02Of3002XXXXeVOXPmyNChQ1c7/UMf+lD/f23eLVZmKr3//e+X++67T4455hi58MILZcKECf2nL126VHbccUf553/+Z7nkkktUV25xU6lKKJLBbEz8CIUWOi1iXAPQHHMJc0nX953Rmv6nD8NoYBEFbSyohonpTRvNGtCNiLdrDKOcrygwttcjBcaYS/Z1VDRDCnVWxKDJ4/BukjZzQQACMRLAXFKqtnDhQsmeTsqeYJo+fXp/tJNPPjl/Ouiee+6RMWPG5O+4Xrx4sWy44Yay0UYb9Z+XnXPSSSfJDjvskJtUq/5GUufEH//4x/LSSy+tlO0f//jH3Lzad9995YADDpBtt902f3opM5Gy33saO3asXH/99f1PQV1zzTX56/a++MUvyj777KO6coubSlVCkQxmY+JHKLTQaRHjGoDmmEuYS7q+x1wKw89rFNZIr8r4e52pxR6A+rOvPxjDOAQBzKUQFHUx6GUdv6qj4V2VGOdDAAKpEcBcCqB49jq77LV2EydOlK222kruv/9+ueOOO2Tq1Kkybdq0fIYHH3xQJk+evNLfPfLII/lTTWuvvbaceeaZMnz48JWy2XTTTWW77bbrmmG331zKBtx888252bXzzjvnr9r7zW9+I3Pnzs3jXXXVVf2v1qt7+RY3lXVziWkcGxM/aqGFTosY1wA0x1zCXNL1PeZSGH5eo7BGelUGc8mvMnFlRo/b65UCY8wl+zoqmiGFOiti0ORxeDdJm7kgAIEYCWAuBVDttddek8svvzw3mJ555pn8iaHMNDryyCP7TZzBzKXs/BkzZnTNIHsiatasWbXMpWxQ9rq97Gmo7LeaNthgg/wJp+yJqpEjR6qvOsYvltUXHSAAG5MAEAOFQAsdyBjXADTHXMJc0vU95lIYfl6jsEZ6VQZzya8ycWVGj9vrlQJjzCX7OiqaIYU6K2LQ5HF4N0mbuSAAgRgJYC7FqJqDnGP8YtkBtv7fzRo3bpyHdJLOgU2iTv4Y1wA0x1zCXNL1PeZSGH5eo7BGelUGc8mvMnFlRo/b65UCY8wl+zoqmiGFOiti0ORxeDdJm7kgAIEYCWAuxaiag5xj/GLZATbMJQ8ivJUDm0SdGDGuAWiOuYS5pOt7zKUw/LxGYY30qgzmkl9l4sqMHrfXKwXGmEv2dVQ0Qwp1VsSgyePwbpI2c0EAAjESwFyKUTUHOcf4xbIDbJhLHkTAXAqiQoxrADcGmEuYS0HaXzT9Tx+G0cAiCtpYUA0T05s2mjWgGxFv1xhGOV9RYGyvRwqMMZfs66hohhTqrIhBk8fh3SRt5oIABGIkgLkUo2oOcra4qXRwWeYpsDExR1x6ArQojWrQE2NcA9AccwlzSdf3ndGa/qcPw2hgEQVtLKiGielNG80agLkUpibqRPFWR3WuwfuYFBhjLrVfhSnUWfuUV2QAb09qkAsEIOCRAOaSR1UiyMnipjKCy1anyMZEjTBYALTQoYxxDUBzzCXMJV3fYy6F4ec1CmukV2V4LZ5fZeLKjB631ysFxphL9nVUNEMKdVbEoMnj8G6SNnNBAAIxEsBcilE1BznH+MWyA2y8Fs+DCG/lwCZRJ0aMawCaYy5hLun6HnMpDD+vUVgjvSqDueRXmbgyo8ft9UqBMeaSfR0VzZBCnRUxaPI4vJukzVwQgECMBDCXYlTNQc4xfrHsABvmkgcRMJeCqBDjGsCNAeYS5lKQ9uc3l8JgdBeFNdKdJP0JedPGYg/g7Rr9VkP9zGBcn13ZkSkwxlwqWw1256VQZ3b0qkeGd3VmjIAABNIigLmUlt7BrtbipjJYco4DsTHxIw5a6LSIcQ1Ac8wlzCVd33dGa/qfPgyjgUUUtLGgGiamN200a0A3It6uMYxyvqLA2F6PFBhjLtnXUdEMKdRZEYMmj8O7SdrMBQEIxEgAcylG1RzkbHFT6eCyzFNgY2KOuPQEaFEa1aAnxrgGoDnmEuaSru8xl8Lw8xqFNdKrMrwWz68ycWVGj9vrlQJjzCX7OiqaIYU6K2LQ5HF4N0mbuSAAgRgJYC7FqJqDnGP8YtkBNl6L50GEt3Jgk6gTI8Y1AM0xlzCXdH2PuRSGn9corJFelcFc8qtMXJnR4/Z6pcAYc8m+jopmSKHOihg0eRzeTdJmLghAIEYCmEsxquYg5xi/WHaADXPJgwiYS0FUiHEN4MYAcwlzKUj785tLYTC6i8Ia6U6S/oS8aWOxB/B2jX6roX5mMK7PruzIFBhjLpWtBrvzUqgzO3rVI8O7OjNGQAACaRHAXEpL72BXa3FTGSw5x4HYmPgRBy10WsS4BqA55hLmkq7vO6M1/U8fhtHAIgraWFANE9ObNpo1oBsRb9cYRjlfUWBsr0cKjDGX7OuoaIYU6qyIQZPH4d0kbeaCAARiJIC5FKNqDnK2uKl0cFnmKbAxMUdcegK0KI1q0BNjXAPQHHMJc0nX95hLYfh5jcIa6VUZXovnV5m4MqPH7fVKgTHmkn0dFc2QQp0VMWjyOLybpM1cEIBAjAQwl2JUzUHOMX6x7AAbr8XzIMJbObBJ1IkR4xqA5phLmEu6vsdcCsPPaxTWSK/KYC75VSauzOhxe71SYIy5ZF9HRTOkUGdFDJo8Du8maTMXBCAQIwHMpRhVc5BzjF8sO8CGueRBBMylICrEuAZwY4C5hLkUpP35zaUwGN1FYY10J0l/Qt60sdgDeLtGv9VQPzMY12dXdmQKjDGXylaD3Xkp1JkdveqR4V2dGSMgAIG0CGAupaV3sKu1uKkMlpzjQGxM/IiDFjotYlwD0BxzCXNJ1/ed0Zr+pw/DaGARBW0sqIaJ6U0bzRrQjYi3awyjnK8oMLbXIwXGmEv2dVQ0Qwp1VsSgyePwbpI2c0EAAjESwFyKUTUHOVvcVDq4LPMU2JiYIy49AVqURjXoiTGuAWiOuYS5pOt7zKUw/LxGYY30qgyvxfOrTFyZ0eP2eqXAGHPJvo6KZkihzooYNHkc3k3SZi4IQCBGAphLMarmIOcYv1h2gI3X4nkQ4a0c2CTqxIhxDUBzzCXMJV3fYy6F4ec1CmukV2Uwl/wqE1dm9Li9Xikwxlyyr6OiGVKosyIGTR6Hd5O0mQsCEIiRAOZSjKo5yDnGL5YdYMNc8iAC5lIQFWJcA7gxwFzCXArS/vzmUhiM7qKwRrqTpD8hb9pY7AG8XaPfaqifGYzrsys7MgXGmEtlq8HuvBTqzI5e9cjwrs6MERCAQFoEMJfS0jvY1VrcVAZLznEgNiZ+xEELnRYxrgFojrmEuaTr+85oTf/Th2E0sIiCNhZUw8T0po1mDehGxNs1hlHOVxQY2+uRAmPMJfs6KpohhTorYtDkcXg3SZu5IACBGAlgLsWomoOcLW4qHVyWeQpsTMwRl54ALUqjGvTEGNcANMdcwlzS9T3mUhh+XqOwRnpVhtfi+VUmrszocXu9UmCMuWRfR0UzpFBnRQyaPA7vJmkzFwQgECMBzKUYVXOQc4xfLDvAxmvxPIjwVg5sEnVixLgGoDnmEuaSru8xl8Lw8xqFNdKrMphLfpWJKzN63F6vFBhjLtnXUdEMKdRZEYMmj8O7SdrMBQEIxEgAcylG1RzkHOMXyw6wYS55EAFzKYgKMa4B3BhgLmEuBWl/fnMpDEZ3UVgj3UnSn5A3bSz2AN6u0W811M8MxvXZlR2ZAmPMpbLVYHdeCnVmR696ZHhXZ8YICEAgLQKYS2npHexqLW4qgyXnOBAbEz/ioIVOixjXADTHXMJc0vV9Z7Sm/+nDMBpYREEbC6phYnrTRrMGdCPi7RrDKOcrCozt9UiBMeaSfR0VzZBCnRUxaPI4vJukzVwQgECMBDCXYlTNQc4WN5UOLss8BTYm5ohLT4AWpVENemKMawCaYy5hLun6HnMpDD+vUVgjvSrDa/H8KhNXZvS4vV4pMMZcsq+johlSqLMiBk0eh3eTtJkLAhCIkQDmUoyqOcg5xi+WHWDjtXgeRHgrBzaJOjFiXAPQHHMJc0nX95hLYfh5jcIa6VUZzCW/ysSVGT1ur1cKjDGX7OuoaIYU6qyIQZPH4d0kbeaCAARiJIC5FKNqDnKO8YtlB9gwlzyIgLkURIUY1wBuDDCXMJeCtD+/uRQGo7sorJHuJOlPyJs2FnsAb9fotxrqZwbj+uzKjkyBMeZS2WqwOy+FOrOjVz0yvKszYwQEIJAWAcyltPQOdrUWN5XBknMciI2JH3HQQqdFjGsAmmMuYS7p+r4zWtP/9GEYDSyioI0F1TAxvWmjWQO6EfF2jWGU8xUFxvZ6pMAYc8m+jopmSKHOihg0eRzeTdJmLghAIEYCmEsxquYgZ4ubSgeXZZ4CGxNzxKUnQIvSqAY9McY1AM0xlzCXdH2PuRSGn9corJFeleG1eH6ViSszetxerxQYYy7Z11HRDCnUWRGDJo/Du0nazAUBCMRIAHMpRtUc5BzjF8sOsPFaPA8ivJUDm0SdGDGuAWiOuYS5pOt7zKUw/LxGYY30qgzmkl9l4sqMHrfXKwXGmEv2dVQ0Qwp1VsSgyePwbpI2c0EAAjESwFyKUTUHOcf4xbIDbJhLHkTAXAqiQoxrADcGmEuYS0Han99cCoPRXRTWSHeS9CfkTRuLPYC3a/RbDfUzg3F9dmVHpsAYc6lsNdidl0Kd2dGrHhne1ZkxAgIQSIsA5lJaege7WoubymDJOQ7ExsSPOGih0yLGNQDNMZcwl3R93xmt6X/6MIwGFlHQxoJqmJjetNGsAd2IeLvGMMr5igJjez1SYIy5ZF9HRTOkUGdFDJo8Du8maTMXBCAQIwHMpRhVc5CzxU2lg8syT4GNiTni0hOgRWlUg54Y4xqA5phLmEu6vsdcCsPPaxTWSK/K8Fo8v8rElRk9bq9XCowxl+zrqGiGFOqsiEGTx+HdJG3mggAEYiTQc+bS5MmTa+nQ19cnc+fOrTU2xUExfrHsQSc2Jh5U4Av2ECrEuAbQf9Q+5lKI7hdeixcGo7sorJHuJOlPyJs2FnsAb9fotxrqZwbj+uzKjkyBMeZS2WqwOy+FOrOjVz0yvKszYwQEIJAWgZ4zl3bffffaCn7ve9+rPTa1gRY3lSkwZGPiR2W00GkR4xqA5phLmEu6vu+M1vQ/fRhGA4soaGNBNUxMb9po1oBuRLxdYxjlfEWBsb0eKTDGXLKvo6IZUqizIgZNHod3k7SZCwIQiJFAz5lLMYoQY84WN5UxcqiaMxuTqsTszkcLHdsY1wA0x1zCXNL1PeZSGH5eo7BGelWG1+L5VSauzOhxe71SYIy5ZF9HRTOkUGdFDJo8Du8maTMXBCAQIwHMpRhVc5BzjF8sO8AmbEw8qMAX7CFUiHENoP+ofcylEN3Pa/HCUPQXhTXSnyadjLxpY7EH8HaNfquhfmYwrs+u7MgUGGMula0Gu/NSqDM7etUjw7s6M0ZAAAJpEUjKXPrlL38pTz75pCxZskQmTJiQltKBr9bipjJwii7DsTHxIwta6LSIcQ1Ac8wlzCVd33dGa/qfPgyjgUUUtLGgGiamN200a0A3It6uMYxyvqLA2F6PFBhjLtnXUdEMKdRZEYMmj8O7SdrMBQEIxEggCXPppz/9qXzqU5+SJ554ol+jRYsW5f87O/axj31MZs+eLXvuuWeMGraSs8VNZSsX0vCkbEwaBr6G6dBCp0WMawCaYy5hLun6HnMpDD+vUVgjvSrDa/H8KhNXZvS4vV4pMMZcsq+johlSqLMiBk0eh3eTtJkLAhCIkUDPm0uPPPKIHHHEEbLuuuvKIYcckhtMP/jBD6RjLmWiZabSVlttJV/4whdi1LCVnGP8YrkVUKtMysbEgwp8wR5ChRjXAPqP2sdcCtH9vBYvDEV/UVgj/WnSycibNhZ7AG/X6Lca6mcG4/rsyo5MgTHmUtlqsDsvhTqzo1c9MryrM2MEBCCQFoGeN5f+9V//VTKD6eabb5a3ve1tMmfOHLn00ktXMpdOO+00efTRR+W//uu/0lJfcbUWN5WKdKIZysbEj1RoodMixjUAzTGXMJd0fd8Zrel/+jCMBhZR0MaCapiY3rTRrAHdiHi7xjDK+YoCY3s9UmCMuWRfR0UzpFBnRQyaPA7vJmkzFwQgECOBnjeXtt9+e9l3333lM5/5TK7PYObS5z//ebn66qvl5z+puFMJAAAgAElEQVT/eYwatpKzxU1lKxfS8KRsTBoGvobp0EKnRYxrAJpjLmEu6foecykMP69RWCO9KsNr8fwqE1dm9Li9Xikwxlyyr6OiGVKosyIGTR6Hd5O0mQsCEIiRQM+bS9tuu6185CMfkRkzZnQ1l/793/9d7rjjjvz3l/iUIxDjF8vlrsz2LDYmtnyrREeLKrRWPzfGNQDNMZcwl3R9j7kUhp/XKKyRXpXBXPKrTFyZ0eP2eqXAGHPJvo6KZkihzooYNHkc3k3SZi4IQCBGAj1vLh100EEydOhQmT9//qDm0rJly2S//faTjTbaSK699toYNWwl5xi/WG4F1CqTsjHxoAJfsIdQIcY1gP6j9jGXQnQ/v7kUhqK/KKyR/jTpZORNG4s9gLdr9FsN9TODcX12ZUemwBhzqWw12J2XQp3Z0aseGd7VmTECAhBIi0DPm0tXXXWVzJw5U0444QSZNm1a/ntLnd9cWrp0qcyaNSs3lc4//3w55JBD0lJfcbUWN5WKdKIZysbEj1RoodMixjUAzTGXMJd0fd8Zrel/+jCMBhZR0MaCapiY3rTRrAHdiHi7xjDK+YoCY3s9UmCMuWRfR0UzpFBnRQyaPA7vJmkzFwQgECOBnjeXsieTpk6dKvfee6+MGTNG1l13Xfntb38r73vf++Sxxx6T5557Tvbaay/50pe+FKN+reVscVPZ2sU0ODEbkwZhF0yFFjotYlwD0BxzCXNJ1/eYS2H4eY3CGulVGV6L51eZuDKjx+31SoEx5pJ9HRXNkEKdFTFo8ji8m6TNXBCAQIwEet5cykRZvny5XHPNNfkTSk888UT+79nnHe94h3z0ox+VyZMnS19fX4z6tZZzjF8stwZrwMRsTDyowBfsIVSIcQ2g/6h9zKUQ3c9r8cJQ9BeFNdKfJp2MvGljsQfwdo1+q6F+ZjCuz67syBQYYy6VrQa781KoMzt61SPDuzozRkAAAmkRSMJcGijpK6+8Ii+88IKMGDFCRo4cGUTt119/XS6//HK58cYb5dlnn5WxY8fKEUccIZMmTVqjafX888/Ld77znfypqsz0WrJkibz97W+X/fffX6ZMmSLrrLPOSvldeOGF8tBDD8nvfve7/NzsSaz3vve9+Sv/Ntlkk5XOPfLII+UnP/nJate31lprycKFC9XXbXFTqU4qggBsTPyIhBY6LWJcA9AccwlzSdf3ndGa/qcPw2hgEQVtLKiGielNG80a0I2It2sMo5yvKDC21yMFxphL9nVUNEMKdVbEoMnj8G6SNnNBAAIxEkjOXLIQ6ZxzzpH58+fLxIkTZeutt5Yf/vCHcuedd+a/8ZS9kq/bJzOVTjzxRNlll11kp512ys2uzDy69dZbZfz48TJv3jzJzKDOJzOrtthiC9l0001zcywzmW644QbJzK3M2MqMqc4nM5ey1/6dffbZK00/ZMgQOeCAA9QYLG4q1UlFEICNiR+R0EKnRYxrAJq/qTkc4KDrfp5c0vLzOp61wasy/tZtiz0A9WdffzCGcQgCmEshKOpi0Ms6flVHw7sqMc6HAARSI9Bz5tLixYtra7jxxhtXHrto0SKZMGGCHH300TJ9+vT+8SeffLLcc889+T+jR48eNG5mDmWfgaZQ9u+XXHKJXHbZZTJnzpz896DW9Hn00UflkEMOkWOPPVZOO+20/lMzc+k3v/mN/OAHP6h8TWUGWNxUlpk39nPYmPhREC10WsS4BqD5m5rDAQ667sdc0vLzOp61wasy/tZtiz0A9WdffzCGcQgCmEshKOpi0Ms6flVHw7sqMc6HAARSI9Bz5tK4ceNq/35SZhRV/Vx88cX5K/Gyp5AGmlMPP/ywHH744XLuuefmf1b5ZE8cHXjggXLSSSflr7xb0+dPf/qT7LzzznLYYYfJeeed139qx1zK8speBZg96RTyd6UsbiqrMIr1XDYmfpRDC50WMa4BaI6pMrDqqYf6a4Cm/+Fen7v1SLSxJlw/vjdtNGtANwrerrG+Wn5HwthemxQYYy7Z11HRDCnUWRGDJo/Du0nazAUBCMRIoOfMpS996UurmSg/+9nP5P7775d3vvOdst1228lf//Vfyx//+EfJNkZPPvlk/lq67O/X9Aq7buJmTyw9/vjj+avwBn6WLl0q22yzjRx88MFywQUXVKqN++67T4455pjcLMpMo4GfN954Q7Lfalq2bJk8/fTTcumll+ZPJ2V/7rnnniuZS9n1DR06VF599dX8lXt77723nH766fn1az8WN5XanGIYz8bEj0poodMixjUAzd/UHA5w0HU/Ty5p+Xkdz9rgVRl/67bFHoD6s68/GMM4BAHMpRAUdTHoZR2/qqPhXZUY50MAAqkR6DlzaVUBH3jggfyVcZ/5zGfkQx/60Gr63nzzzfLJT35Srrjiivx3j6p+9t9/fxk2bJjcdNNNqw3Nnijacsst5corrywdNjOPpkyZItnr7hYsWCCjRo1aaexTTz0le+yxR//fbbjhhnL88cfLv/zLv6x03owZM2TMmDH5bzQtX75cMg7Z7zNlr+DL/lx//fVL5zTYidmmMoubPRHFpzyB7Cmy7DN8+PDygzjThECvaJE9rdnGJ8Y1oFc01+oNhzcJ9gqHNtYATf/3CndtH3ocjzYeVVnzetVG/2cZadaAbpSpP/v6g3HvMW5jDVi1/6kr+7padQaYN8vcK+82+r9Z8swGAQjEQqDnzaXsyZ/sdXXZ6+u6fU455RT5/e9/L9ddd11l3bKnhTIDaLCxu+22W27mzJs3r3Tczmv2zjnnHMlebbfq5y9/+Ytkr9zLnox64okn5NZbb83Npuz1eUOGDFnjPPPnz5csbvaE1rRp00rnNNiJFjeVqoQiGex1YxIJvqBp9ooWbW0qY1wDekVzbSPAYc1f1mr5Nj2+jTVA0//UX9MVUn4+tCnPqukzu2nTRv9n165ZA7qxo/7sqwrGvce4jTUAc8m+jopmoJeLCIU97pV3G/0flizRIACBXiHQ8+bStttumz8JlBlI3T6ZoZMZQNnr86p+Qj659K1vfUvOP//81X4/aU05LV68WLIcMiNqTdfYibHjjjvK5ptvXstIG5iHxeswqrKP8XweqfajGlrotIhxDUDzNzWHAxx03c9r8bT8vI5nbfCqjL9122IPQP3Z1x+MYRyCAK/FC0FRF4Ne1vGrOhreVYlxPgQgkBqBnjeX/vEf/1E23XRTueaaa7pq+9GPflR++9vf5r/LVPVT9JtLBx10kMycObMwbPZavbPOOkv2228/ueiiiwqfQhoYMPt9pscee0yy32oq+mT5LFmyRO66666iU9d43OKmUpVQJIPZmPgRCi10WsS4BqA5psrAqqce6q8Bmv6He33u1iPRxppw/fjetNGsAd0oeLvG+mr5HQlje21SYIy5ZF9HRTOkUGdFDJo8Du8maTMXBCAQI4GeN5ey31q6+uqrZcKECfmr4LJX5HU+2VM/X/rSlyT73aVJkyblr4yr+pk9e3b+e0333nvvSrGzV9cdfvjh+e85ZbHX9Lntttvk9NNPl1133VXmzJkjQ4cOrZRG9tTSI488Ir/4xS/WOC77PafsyaV3vetdcu2111aaY9WTLW4qVQlFMpiNiR+h0EKnRYxrAJq/qTkc4KDrfp5c0vLzOp61wasy/tZtiz0A9WdffzCGcQgCmEshKOpi0Ms6flVHw7sqMc6HAARSI9Dz5tLLL78s//qv/5q/G3yttdaSv/mbv5GNNtpI/vSnP8mzzz4ry5Ytk/e85z3y1a9+VdZbb73K+i9cuFCyp4GyJ5imT5/eP/7kk0+WBQsWyD333CNjxozJfzg8M7M23HDDfP7OJzvnpJNOkh122CE3qYYNGzZoDi+88IIMHz5c1l577ZWOZ/9Hd+ihh0r2+r/Obzu99NJL+XnrrLPOSudeeeWV8vnPf15OPfVUOe644ypf68ABFjeVqoQiGczGxI9QaKHTIsY1AM0xVQZWPfVQfw3Q9D/c63O3Hok21oTrx/emjWYN6EbB2zXWV8vvSBjba5MCY8wl+zoqmiGFOiti0ORxeDdJm7kgAIEYCfS8uZSJkj2xc8MNN8gtt9wijz/+uGTmy8iRI+Xd7363HHjggfLhD3+40mvoVhU6e51d9lq7iRMnylZbbZW/Xu+OO+6QqVOn5k9LZZ8HH3xQJk+evNLfZU8bZU81ZUbQmWeemZtHAz/Z6/y22267/K8yE+pTn/qU7LvvvrLZZpvlRtkvf/nL/Kmr7DN37lzZeuut++fKfn8pe8VeFqOvry+f/+6775bsR/+yp5bqGGkDc7O4qYyxgarmzMakKjG789FCxzbGNQDN39QcDnDQdT9PLmn5eR3P2uBVGX/rtsUegPqzrz8YwzgEAcylEBR1MehlHb+qo+FdlRjnQwACqRFIwlyyFvW1116Tyy+/PDeYnnnmGRk7dmxuGmWvq8uMnewzmLmUnT9jxoyu6WVPRM2aNSs/nv0m1Je//OX8CaxsjmzO0aNH56+5O/bYY2XzzTfvj/PUU0/lv9v03//93/Lcc8/lT2dtsskmsvfee+fnjhgxQo3E4qZSnVQEAdiY+BEJLXRaxLgGoDmmysCqpx7qrwGa/od7fe7WI9HGmnD9+N600awB3Sh4u8b6avkdCWN7bVJgjLlkX0dFM6RQZ0UMmjwO7yZpMxcEIBAjgSTNpeXLl/ebPjGK5iFni5tKD9dlnQMbE2vC5eOjRXlWg50Z4xqA5phLmEu6vu+M1vQ/fRhGA4soaGNBNUxMb9po1gDMpTA1USeKtzqqcw3ex6TAGHOp/SpMoc7ap7wiA3h7UoNcIAABjwSSMJcyM+n666+X7373u/nrgF599VVZd91181fETZgwIX+dXecJI48ieczJ4qbS43WGzomNSWii9eOhRX122cgY1wA0x1zCXNL1PeZSGH5eo7BGelWG1+L5VSauzOhxe71SYIy5ZF9HRTOkUGdFDJo8Du8maTMXBCAQI4GeN5eWLl0qxx13nPz4xz/ODaS3ve1tMmrUqPx1cX/4wx/y32Paaaed8tfaDRs2LEYNW8k5xi+WWwG1yqRsTDyowBfsIVSIcQ2g/6h9zKUQ3a8zl+nDMBpYREEbC6phYnrTxmIP4O0awyjnKwqM7fVIgTHmkn0dFc2QQp0VMWjyOLybpM1cEIBAjAR63lyaM2eOZP988IMflFNOOSX/7aHOZ/HixTJ79my5/fbbZerUqXLiiSfGqGErOVvcVLZyIQ1PysakYeBrmA4tdFrEuAagOeYS5pKu7zujNf1PH4bRwCIK2lhQDRPTmzaaNaAbEW/XGEY5X1FgbK9HCowxl+zrqGiGFOqsiEGTx+HdJG3mggAEYiTQ8+bSPvvsIxtssIF8+9vf7qrPRz7yEXn++eflrrvuilHDVnK2uKls5UIanpSNScPAMZfMgMe4BtB/mEuYS2GWBE3/04dhNLCIgjYWVMPE9KaNZg3AXApTE3WieKujOtfgfUwKjDGX2q/CFOqsfcorMoC3JzXIBQIQ8Eig582lrbbaSo466ig59dRTu/K/+OKL5Rvf+IY8+uijHjVymZPFTaXLCw2cFBuTwEAV4dBCAY/fXNLBa3k0tY/Jpi1BzR6A+tPStxuPNnZstZG9aaNZAzCXtNVQf7y3Oqp/JX5HpsAYc6n9+kuhztqnjLnkSQNygQAEfBPoeXPpfe97n+yyyy7y+c9/vqsSZ5xxhtx///3yox/9yLdajrKzuKl0dHlmqbARNENbOTBaVEa20oAY1wA0x1QZWMTUQ/01QNP/cK/P3Xok2lgTrh/fmzaaNQBzqX4daEd6qyPt9XgcnwJjzKX2Ky+FOmufMuaSJw3IBQIQ8E2g582l7HeW7r77brn00ktl1113XU2N73//+/lvLe29996SPcHEpxwBi5vKcjPHfRYbQT/6oYVOixjXADTHXMJc0vV9Z7Sm/+nDMBpYREEbC6phYnrTRrMGYC6FqYk6UbzVUZ1r8D4mBcaYS+1XYQp11j5lzCVPGpALBCDgm0DPm0u//vWv5ZBDDpGXX35Ztt9+exk/fryMGjVKnnvuOck2Rj/96U/lr/7qr/LfZHrHO97hWy1H2VncVDq6PLNU2Aiaoa0cGC0qI1tpQIxrAJpjLmEu6foecykMP69RWCO9KiPiTRuLPYC3a/RbDfUzg3F9dmVHpsAYc6lsNdidl0Kd2dGrHhne1ZkxAgIQSItAz5tLmZyPPfaYfOpTn5Kf/exnq6n7nve8R84991x597vfnZbyyqu1uKlUphTFcDYmfmRCC50WMa4BaI65hLmk63vMpTD8vEZhjfSqDOaSX2Xiyowet9crBcaYS/Z1VDRDCnVWxKDJ4/BukjZzQQACMRJIwlzqCPP000/nRtNLL70kI0eOlC222ELGjh0bo26t5xzjF8utQxN/Xw54YNJWDmwSdeRjXAPQHHMJc0nX95hLYfh5jcIa6VUZf/tHiz0A9WdffzCGcQgCmEshKOpi0Ms6flVHw7sqMc6HAARSI5CUuZSauJbXa3FTaZmvl9hsTLwo4e+LGj9kymUS4xpA/2EuYS6V6++iszT9Tx8W0W3vONq0x75oZm/aaNaAbtfq7RqLNInxOIztVUuBMeaSfR0VzZBCnRUxaPI4vJukzVwQgECMBDCXYlTNQc4WN5UOLss8BTYm5ohLT4AWpVENemKMawCaYy5hLun6vjNa0//0YRgNLKKgjQXVMDG9aaNZAzCXwtREnSje6qjONXgfkwJjzKX2qzCFOmuf8ooM4O1JDXKBAAQ8EkjCXHrggQdk7ty58vjjj8szzzwjy5YtW02Lvr4+WbhwoUeNXOZkcVPp8kIDJ8XGJDBQRTi0UMATkRjXADTHXMJc0vU95lIYfl6jsEZ6Vcbf09YWewDqz77+YAzjEAQwl0JQ1MWgl3X8qo6Gd1VinA8BCKRGoOfNpeuuu07OO+88Wb58uWy22WYyatQoGTJkyKA6z5s3LzX9a1+vxU1l7WQiGsjGxI9YaKHTIsY1AM0xlzCXdH2PuRSGn9corJFelcFc8qtMXJnR4/Z6pcAYc8m+jopmSKHOihg0eRzeTdJmLghAIEYCPW8u7b777vKXv/xFvva1r8m4ceNi1MhlzjF+sewBJBsTDyrwBXsIFWJcA+g/ah9zKUT3655cpA/DaGARBW0sqIaJ6U0biz2At2sMo5yvKDC21yMFxphL9nVUNEMKdVbEoMnj8G6SNnNBAAIxEuh5c2mbbbaRiRMnytlnnx2jPm5ztripdHuxARNjYxIQpjIUWugAxrgGoDnmEuaSru87ozX9Tx+G0cAiCtpYUA0T05s2mjWgGxFv1xhGOV9RYGyvRwqMMZfs66hohhTqrIhBk8fh3SRt5oIABGIk0PPm0sEHHyzvfve7ZdasWTHq4zZni5tKtxcbMDE2JgFhKkOhhQ5gjGsAmmMuYS7p+h5zKQw/r1FYI70qw2vx/CoTV2b0uL1eKTDGXLKvo6IZUqizIgZNHod3k7SZCwIQiJFAz5tL//Vf/yWf+MQn5Prrr5e/+7u/i1EjlznH+MWyB5BsTDyowBfsIVSIcQ2g/6h9zKUQ3c9r8cJQ9BeFNdKfJp2MvGljsQfwdo1+q6F+ZjCuz67syBQYYy6VrQa781KoMzt61SPDuzozRkAAAmkR6HlzKZPz9ttvl5kzZ0r2+0tbbLGFjBgxYlCVJ0yYkJb6iqu1uKlUpBPNUDYmfqRCC50WMa4BaI65hLmk6/vOaE3/04dhNLCIgjYWVMPE9KaNZg3oRsTbNYZRzlcUGNvrkQJjzCX7OiqaIYU6K2LQ5HF4N0mbuSAAgRgJ9Ly5tGTJkvzJpbvvvluWL1+ea9TX17eSVtnfZ3+3aNGiGDVsJWeLm8pWLqThSdmYNAx8DdOhhU6LGNcANMdcwlzS9T3mUhh+XqOwRnpVhtfi+VUmrszocXu9UmCMuWRfR0UzpFBnRQyaPA7vJmkzFwQgECOBnjeXMmPp5ptvli233FL22msvGTVqlKy11lqDanXQQQfFqGErOcf4xXIroFaZlI2JBxX4gj2ECjGuAfQftY+5FKL7eS1eGIr+orBG+tOkk5E3bSz2AN6u0W811M8MxvXZlR2ZAmPMpbLVYHdeCnVmR696ZHhXZ8YICEAgLQI9by7tuOOOsvnmm8s111wjQ4YMSUtdw6u1uKk0TNdNaDYmbqQQtNBpEeMagOaYS5hLur7vjNb0P30YRgOLKGhjQTVMTG/aaNaAbkS8XWMY5XxFgbG9Hikwxlyyr6OiGVKosyIGTR6Hd5O0mQsCEIiRQM+bS+9973vlwx/+sEyfPj1GfdzmbHFT6fZiAybGxiQgTGUotNABjHENQHPMJcwlXd9jLoXh5zUKa6RXZXgtnl9l4sqMHrfXKwXGmEv2dVQ0Qwp1VsSgyePwbpI2c0EAAjES6Hlz6eMf/7j8+c9/lrlz58aoj9ucY/xi2QNMNiYeVOAL9hAqxLgG0H/UPuZSiO7ntXhhKPqLwhrpT5NORt60sdgDeLtGv9VQPzMY12dXdmQKjDGXylaD3Xkp1JkdveqR4V2dGSMgAIG0CPS8ufT73/9eJk2aJAceeKCccMIJMmzYsLQUNrpai5tKo1RdhWVj4kcOtNBpEeMagOaYS5hLur7vjNb0P30YRgOLKGhjQTVMTG/aaNaAbkS8XWMY5XxFgbG9Hikwxlyyr6OiGVKosyIGTR6Hd5O0mQsCEIiRQM+bS5MnT5YXXnhBHnvsMRkxYoRsttlm+Z+rfvr6+ni6qUIFW9xUVpg+2lPZmPiRDi10WsS4BqA55hLmkq7vMZfC8PMahTXSqzK8Fs+vMnFlRo/b65UCY8wl+zoqmiGFOiti0ORxeDdJm7kgAIEYCfS8uTRu3LhSumTm0qJFi0qdy0m6V+KkzI+NiR/10UKnBeaSjl+bo6l9TDZt/Wn6n/rT0rcbjzZ2bLWRvWmjWQO6sfB2jVrNPI6Hsb0qKTDGXLKvo6IZUqizIgZNHod3k7SZCwIQiJFAz5tLMYoSQ84WN5UxXLc2RzYmWoLhxqOFjmWMawCaY6oMrHrqof4aoOl/uNfnbj0SbawJ14/vTRvNGoC5VL8OtCO91ZH2ejyOT4Ex5lL7lZdCnbVPeUUG8PakBrlAAAIeCWAudVEl+z+Q7J8JEyZ41K31nCxuKlu/qAYSYGPSAOSSU6BFSVBdTotxDUBzzCXMJV3fd0Zr+p8+DKOBRRS0saAaJqY3bTRrAOZSmJqoE8VbHdW5Bu9jUmCMudR+FaZQZ+1TxlzypAG5QAACvglgLnXRZ86cOXLppZfyqrwe+mLZQyuyEfSgAl+wh1DB4oulEHmtKQb9R+1jLoXpMk3/04dhNLCIgjYWVMPE9KaNZg3AXApTE3WieKujOtfgfUwKjDGX2q/CFOqsfcqYS540IBcIQMA3AcwlzKVaFWpxU1krkcgGsRH0Ixha6LSIcQ1Ac8wlzCVd33dGa/qfPgyjgUUUtLGgGiamN200awDmUpiaqBPFWx3VuQbvY1JgjLnUfhWmUGftU8Zc8qQBuUAAAr4JYC5hLtWqUIubylqJRDaIjaAfwdBCp0WMawCaYy5hLun6HnMpDD+vUVgjvSoj+au6s8+4ceNcJGmxB/B2jS5AB04CxoGBDhIuBcaYS/Z1VDRDCnVWxKDJ4/BukjZzQQACMRLAXMJcqlW3FjeVtRKJbBAbEz+CoYVOixjXADTHXMJc0vU95lIYfl6jsEZ6VQZzya8ycWVGj9vrlQJjzCX7OiqaIYU6K2LQ5HF4N0mbuSAAgRgJYC5hLtWq2xi/WK51oYEHsTEJDFQRDi0U8EQkxjUAzTGXMJd0fY+5FIaf1yiskV6VwVzyq0xcmdHj9nqlwBhzyb6OimZIoc6KGDR5HN5N0mYuCEAgRgKYS5hLteo2xi+Wa11o4EFsTAIDVYRDCwU8zCUdvJZHU/uYbNoS1OwBqD8tfbvxaGPHVhvZmzaaNaAbC2/XqNXM43gY26uSAmPMJfs6KpohhTorYtDkcXg3SZu5IACBGAlgLmEu1apbi5vKWolENoiNiR/B0EKnRYxrAJpjqgyseuqh/hqg6X+41+duPRJtrAnXj+9NG80agLlUvw60I73VkfZ6PI5PgTHmUvuVl0KdtU95RQbw9qQGuUAAAh4JYC5hLtWqS4ubylqJRDaIjYkfwdBCp0WMawCaYy5hLun6vjNa0//0YRgNLKKgjQXVMDG9aaNZAzCXwtREnSje6qjONXgfkwJjzKX2qzCFOmufMuaSJw3IBQIQ8E0AcwlzqVaFWtxU1kokskFsBP0IhhY6LWJcA9AccwlzSdf3mEth+HmNwhrpVRl+c8mvMnFlRo/b65UCY8wl+zoqmiGFOiti0ORxeDdJm7kgAIEYCWAudVEt+z+QRYsWyUEHHRSjruY5x/jFsjmUEhOwMSkBqaFT0EIHOsY1AM0xlzCXdH2PuRSGn9corJFelcFc8qtMXJnR4/Z6pcAYc8m+jopmSKHOihg0eRzeTdJmLghAIEYCmEsxquYg5xi/WHaATdiYeFCBL9hDqBDjGkD/UfuYSyG6X0TT//RhGA0soqCNBdUwMb1po1kDuhHxdo1hlPMVBcb2eqTAGHPJvo6KZkihzooYNHkc3k3SZi4IQCBGAj1nLo0bN076+voqa5GNWbhwYeVxqQ6wuKlMgSUbEz8qo4VOixjXADTHXMJc0vV9Z7Sm/+nDMBpYREEbC6phYnrTRrMGYC6FqYk6UbzVUZ1r8D4mBcaYS+1XYQp11j7lFRnA25Ma5AIBCHgk0HPm0ic+8Yla5lImzmc/+1mPGrnMyeKm0uWFBk6KjUlgoIpwaKGAJ7onF3Qz1x+N5phLmEv1+2fgSM0egD4Mo4FFFLSxoBompjdtNGsA5lKYmqgTxVsd1bkG72NSYIy51H4VplBn7VPGXPKkAblAAI1zuKcAACAASURBVAK+CfScueQbd+9kZ3FT2Tt0ul8JG0E/KqOFTosY1wA0x1zCXNL1fWe0pv/pwzAaWERBGwuqYWJ600azBmAuhamJOlG81VGda/A+JgXGmEvtV2EKddY+ZcwlTxqQCwQg4JsA5pJvfdxmZ3FT6fZiAybmaSP4ysuvytJXlsqw4cNk+Ih1A15lHKE8aREHsZWzjHENQPMV5tLQoUNl7Ns2YQ0Qkex1unyqEdD0P31YjbXV2YPtAdDGirY+rjdtNGsA5pK+HupG8FZHda/D87gUGGMutV+B2jpL/XuAqgpqeVedj/MhAAEIxEYAcyk2xZzka3FT6eTSTNPwsDF58U8vyR9+/Yxc/7nvyu+f+IO87W//j3zkzA/J/3nHaPmrjUaaXr+n4B608MSjai4xrgFo/qbKTz35tLz47BK58eJbWQMwl6q2fn6+pv/pw1rIgw1a0x7g6WeeyufBcA2GO1ggb32jWQO6QfF2jcHEcxQIxvZipMAYc8m+jopmqFtnfA9QRHbw43V515uNURCAAATiI9Bz5tIee+xRS4W+vj5ZsGBBrbEpDrK4qUyBY9sbk2xDOf/iW+TamTethvujZx0sh556QDIGU9taxF7vMa4BaC6SrQHfvui7ct2sm1kD/t//44v0mguRpv/pw5rQAwwr2gPs/bF/lpf+8iLmUgDWoUN46xvNGoC5FLo6ysfzVkflM4/nzBQYYy61X4916qxoD5DS9wBVFazDu+ocnA8BCEAgZgI9Zy4deeSRtfWYN29erbGvv/66XH755XLjjTfKs88+K2PHjpUjjjhCJk2aJJlp1e3z/PPPy3e+8x2599575YknnpAlS5bI29/+dtl///1lypQpss4666w09MILL5SHHnpIfve73+XnjhkzRt773vfKCSecIJtssslq0zzwwANyySWXyKJFi2T48OHy/ve/X8444wzZaKONal3nwEEWN5XqpCII0PbG5Jf/36/khO2ndyV12U8vlL8b/84ISOpTbFsL/RW0GyHGNQDNRVgDVvQN9VB/DdH0P9zrc9eOLOr/Lz44U9besE/e9a53aadifGAC3vpGswZ0Q+PtGgNL6CIcjO1lSIEx5pJ9HRXNUKfOivYAKX0PUMR31eN1eFedg/MhAAEIxEyg58ylNsQ455xzZP78+TJx4kTZeuut5Yc//KHceeedMm3aNJk6dWrXlDJT6cQTT5RddtlFdtppJxk5cmRuHt16660yfvx4ycyutdZaq398ZlZtscUWsummm8qIESNyk+mGG26QzNzKjK3MmOp8fvKTn8hRRx2Vn3/IIYfIn/70J/n6178uG2+8cT5m3XV1v7FjcVPZhnZNz9nmxiR7t/Lsj31Zvv/tH3W97N0+8j457cp/k3UT+A2mNrVouu4s5otxDUhdc9aAlTsh9XrQrAua/oe7hnz9sWX7f8qsQ2WTzVb/D5bqz8zIEAS89Y1mDcBcClER9WJ4q6N6V+F7VAqMMZfar8GqdVZ2D5DK9wBVFazKu2p8zocABCAQOwHMJaWC2VNBEyZMkKOPPlqmT1/xRMjJJ58s99xzT/7P6NGjB50lM4eyz0BTKPv37Gmjyy67TObMmSN77bXXGjN89NFHc/Po2GOPldNOO63/3CynP//5z3LbbbfJeuutl//997///fy8GTNmyL/8y7+ortziplKVUCSD29yY/Pm5F+SsD1wgjz/8q6603r3938rM28+SDUatHwnR+mm2qUX9rP2MjHENSF1z1oCV+yf1etCsJpr+h7uGfP2xZft/+rUnyqZ/u+I/Vqo/IyNDEvDWN5o1oBsXb9cYUj8vsWBsr0QKjDGX7OuoaIaqdVZ2D5DK9wBFfFc9XpV31ficDwEIQCB2AphLSgUvvvji/JV42VNI2VNBnc/DDz8shx9+uJx77rn5n1U+jz32mBx44IFy0kkn5a+8W9MneyJp5513lsMOO0zOO++8/NQnn3xS9t1330GfnNp7771lgw02yJ+00nwsbio1+cQyts2NCf/F0spV0qYWsdTrmvKMcQ1IXXPWANaAUGuPpv9T78NQGlSNU7b/eXKpKtlmzvfWN5o1oBsxb9fYjLLNzgJje94pMMZcsq+johmq1lnZPQBPLg1OvirvIv04DgEIQKDXCPS8uTR58uRSmmW/jTR37txS5w48KXti6fHHH89fhTfws3TpUtlmm23k4IMPlgsuuKBS3Pvuu0+OOeaY3CzKTKOBnzfeeEOy32patmyZPP3003LppZfKD37wg/zPPffcMz/1lltukdNPP12uvPJK+ad/+qeVxmd/f9ddd8nPf/7zlV65VylBEbG4qayaQ4znt70x4V3LK6qmbS1irN+BOce4BqA5v7k0sIaph/qrkKb/4V6fu3Zk0R6A31zSErYb761vNGtAN0rertFOzfYiw9iefQqMMZfs66hohjp1VrQH4DeXulOvw7tIQ45DAAIQ6CUCPW8u7b777oPq9fLLL+evjctMpVGjRsnaa68t3/ve9ypru//++8uwYcPkpptuWm1s9kTRlltumZs8ZT+ZeTRlyhTJXne3YMGCPLeBn6eeekr22GOP/r/acMMN5fjjj1/pNXdf+9rX5HOf+5z853/+Z/6bSwM/2d9nx++///7VYpfNMTsv21QuX748/+0nPuUJvPLKK/nJw4cPLz8o4Jkj1/krueur/1eum3XzalEPm3GQ7HPMrvLSX14MOKPfUG1rEYrMuHHjQoWqFCfGNaBXNK8k1ConswasANIr9dDGGqDp/17hrunDtsYW9f+eR/2j/OnFP7a2R2mLSwzzduubNvrf6j6AtcG+EmHce4zbWANW3QNQV/Z1teoMdZgX7QFS+h6gqmJ1eFedo875bfR/nTwZAwEI9D6BnjeX1iRh9uTPhRdeKP/zP/8jX//612sZJdnTQpkBdN1116021W677Zb/ntK8efNKV1LnNXvnnHOOHHnkkauN+8tf/iLZK/eyJ6OeeOIJufXWW3OzKXt93pAhQ/Lzs6eYvvjFL8qdd94pm2+++UoxOr/nlP0W1Cab1P/BZs0XS6Vh9OCJHjYm2cbyxeeWyE0X3yqLn/gf2fhvx8jBp+4vfzVqvWSMpay0PGgRosTb2lTGuAb0iubauhk+dD1Z8qdX5aYv3MYa0KLZr9WxM76NNUDT//RhKOXrxVnTHuDZ55/Jg7b1H8DUu6I0RmEupaGz9VWy/loTbv7+wsMegLqyr6tVZ6jLnO8B6mlVl3e92cqPaqP/y2fHmRCAQEoEkjaXMqFff/11Oeigg2T77bfPfx+p6ifkk0vf+ta35Pzzz1/p95OK8lm8eLFkOWRG1CmnnJKf3tSTS9lc48ePL0qR4wMIeHqk+tWXX5W/vLJU1hk+TNYdsW5yOnnSIkb4Fq/EseaA5m8SzjgMHTpUNnnbJqwBIsKNWfXO0/Q/fVidt8WIwfYAaGNBOkxMb9po1oBuRLxdYxjlfEWBsb0eKTDmtXj2dVQ0g7bOUv8eoIjvqse1vKvOx/kQgAAEYiOQvLmUCZb9JtJtt90mP/rRjyrrV/SbS5lxNXPmzMK42Wv1zjrrLNlvv/3koosu6n8KqXCgSP77TI899phkv9WUfYp+cyl7oukXv/gFv7lUBm7gc9iYBAaqCIcWCniR/u4amr+pORzgoOt+3e8uUn9a+nbj0caOrTayN20wl7SKtjPeWx21Q8F21hQYYy7Z1lCZ6CnUWRkOTZ0D76ZIMw8EIBArAcwlEfnEJz4hd9xxR264VP3Mnj1brrjiCrn33ntl44037h+evbru8MMPl09+8pMyadKkNYbNjK3TTz9ddt11V5kzZ07+X5RX+WRPLT3yyCP9+f/qV7+SD3zgAzJt2jSZOnXqSqH23ntvWX/99eWGG26oMsVq51rcVKoSimQwGxM/QqGFTosY1wA0x1QZWPXUQ/01QNP/cK/P3Xok2lgTrh/fmzaaNaAbBW/XWF8tvyNhbK9NCowxl+zrqGiGFOqsiEGTx+HdJG3mggAEYiSQtLm0fPny/DeLZsyYIdtss41cffXVlTVcuHBh/lq97Amm6dOn948/+eSTZcGCBZL9ttGYMWPy33fJXmG34YYbykYbbdR/XnbOSSedJDvssENuUg0bNmzQHF544YX8Hfhrr732Ssez/6M79NBDZdttt13pt50+9KEPSTYmM67WW2+9fMz3v/99OfbYY/M8s3w1H4ubSk0+sYxlY+JHKbTQaRHjGoDmb2oOBzjoup8nl7T8vI5nbfCqjL9122IPQP3Z1x+MYRyCAOZSCIq6GPSyjl/V0fCuSozzIQCB1Aj0vLm0xx57DKrpsmXL5I9//GP+m0sjRoyQb3zjG7LVVlvV0j97nV32WruJEyfmMe6///78SajsqaHs6aHs8+CDD8rkyZNX+rvsaaPsqabMMDrzzDNX+wHlTTfdVLbbbrt8fGZCfepTn5J9991XNttss/yVdr/85S/l5ptvzo/PnTtXtt566/78f/zjH+cGUvZbEpn5lF1rdo2Z0XXjjTeqf6zZ4qayFvzIBrEx8SMYWui0iHENQHNMlYFVTz3UXwM0/Q/3+tytR6KNNeH68b1po1kDulHwdo311fI7Esb22qTAGHPJvo6KZkihzooYNHkc3k3SZi4IQCBGAj1vLmWvjBvsM2TIkPz1cFtuuaUcfPDBMnr06Nr6vfbaa3L55ZfnBtMzzzwjY8eOzU2jbO6+vr487mDmUnZ+9tRUt0/2RNSsWbPyw7/97W/ly1/+smSbuWyObM4s5x133DF/GmnzzTdfLUz2G1KXXHKJLFq0KDeTdtttNznjjDNk1KhRta+1M9DiplKdVAQB2Jj4EQktdFrEuAag+ZuawwEOuu7nySUtP6/jWRu8KuNv3bbYA1B/9vUHYxiHIIC5FIKiLga9rONXdTS8qxLjfAhAIDUCPW8upSZoU9drcVPZVO5tzsPGpE36K8+NFjotYlwD0BxTZWDVUw/11wBN/8O9PnfrkWhjTbh+fG/aaNaAbhS8XWN9tfyOhLG9Nikwxlyyr6OiGVKosyIGTR6Hd5O0mQsCEIiRAOZSjKo5yNniptLBZZmnwMbEHHHpCdCiNKpBT4xxDUBzzCXMJV3fd0Zr+p8+DKOBRRS0saAaJqY3bTRrAOZSmJqoE8VbHdW5Bu9jUmCMudR+FaZQZ+1TXpEBvD2pQS4QgIBHAphLHlWJICeLm8oILludIhsTNcJgAdBChzLGNQDNMZcwl3R9j7kUhp/XKKyRXpXhtXh+lYkrM3rcXq8UGGMu2ddR0Qwp1FkRgyaPw7tJ2swFAQjESCAJc+nee++VefPmycKFC+XFF1+UN954YzWtst9Gyo7zKUcgxi+Wy12Z7VlsTGz5VomOFlVorX5ujGsAmmMuYS7p+h5zKQw/r1FYI70qg7nkV5m4MqPH7fVKgTHmkn0dFc2QQp0VMWjyOLybpM1cEIBAjAR63ly68cYb5ZxzzpG11lpL3vOe98jo0aNl6NChg2r12c9+NkYNW8k5xi+WWwG1yqRsTDyowBfsIVSIcQ2g/6h9zKUQ3S+i6X/6MIwGFlHQxoJqmJjetNGsAd2IeLvGMMr5igJjez1SYIy5ZF9HRTOkUGdFDJo8Du8maTMXBCAQI4GeN5f22WcfWbJkiVx77bWyySabxKiRy5wtbipdXmjgpNiYBAaqCIcWCnii+3JZN3P90WiOuYS5VL9/Bo7U7AHowzAaWERBGwuqYWJ600azBmAuhamJOlG81VGda/A+JgXGmEvtV2EKddY+5RUZwNuTGuQCAQh4JNDz5tJWW20lhx12mJx99tke+Uebk8VNZbQwKiTOxqQCLONT0UIHOMY1AM0xlzCXdH3fGa3pf/owjAYWUdDGgmqYmN600awBmEthaqJOFG91VOcavI9JgTHmUvtVmEKdtU8Zc8mTBuQCAQj4JtDz5tJ+++0n2267rcycOdO3EpFlZ3FTGRmCWumyEayFzWQQWuiwxrgGoDnmEuaSru8xl8Lw8xqFNdKrMvzmkl9l4sqMHrfXKwXGmEv2dVQ0Qwp1VsSgyePwbpI2c0EAAjES6Hlz6dvf/rbMnj1bbr75Znnb294Wo0Yuc47xi2UPINmYeFCBL9hDqBDjGkD/UfuYSyG6X/daTPowjAYWUdDGgmqYmN60sdgDeLvGMMr5igJjez1SYIy5ZF9HRTOkUGdFDJo8Du8maTMXBCAQI4GeN5ceeugh+cY3viE///nPZfLkybLFFlvIyJEjB9Vqhx12iFHDVnK2uKls5UIanpSNScPA1zAdWui0iHENQHPMJcwlXd93Rmv6nz4Mo4FFFLSxoBompjdtNGtANyLerjGMcr6iwNhejxQYYy7Z11HRDCnUWRGDJo/Du0nazAUBCMRIoOfNpXHjxklfX58sX748/3NNn0WLFsWoYSs5W9xUtnIhDU/KxqRh4JhLZsBjXAPoP8wlzKUwS4Km/+nDMBpYREEbC6phYnrTRrMGYC6FqYk6UbzVUZ1r8D4mBcaYS+1XYQp11j7lFRnA25Ma5AIBCHgk0PPm0pe+9KVCU6kjzNSpUz1q5DIni5tKlxcaOCk2JoGBKsKhhQKe6F6LpZu5/mg0x1zCXKrfPwNHavYA9GEYDSyioI0F1TAxvWmjWQMwl8LURJ0o3uqozjV4H5MCY8yl9qswhTprnzLmkicNyAUCEPBNoOfNJd/4483O4qYyXhrlM2cjWJ6V9ZlooSMc4xqA5phLmEu6vu+M1vQ/fRhGA4soaGNBNUxMb9po1gDMpTA1USeKtzqqcw3ex6TAGHOp/SpMoc7ap4y55EkDcoEABHwTwFzqos/cuXPlqquuknvuuce3gi1lZ3FT2dKlNDotG8FGca9xMrTQaRHjGoDmmEuYS7q+x1wKw89rFNZIr8qIeNPGYg/g7Rr9VkP9zGBcn13ZkSkwxlwqWw1256VQZ3b0qkeGd3VmjIAABNIigLnURe85c+bIpZdeKvwO0+CALG4qU2g9NiZ+VEYLnRYxrgFojrmEuaTre8ylMPy8RmGN9KoM5pJfZeLKjB631ysFxphL9nVUNEMKdVbEoMnj8G6SNnNBAAIxEsBcwlyqVbcxfrFc60IDD2JjEhioIhxaKODxm0s6eC2PpvYx2bQlqNkDUH9a+nbj0caOrTayN200a0A3Ft6uUauZx/EwtlclBcaYS/Z1VDRDCnVWxKDJ4/BukjZzQQACMRLAXMJcqlW3FjeVtRKJbBAbEz+CoYVOixjXADTHVBlY9dRD/TVA0/9wr8/deiTaWBOuH9+bNpo1AHOpfh1oR3qrI+31eByfAmPMpfYrL4U6a5/yigzg7UkNcoEABDwSwFzCXKpVlxY3lbUSiWwQGxM/gqGFTosY1wA0x1zCXNL1fWe0pv/pwzAaWERBGwuqYWJ600azBmAuhamJOlG81VGda/A+JgXGmEvtV2EKddY+ZcwlTxqQCwQg4JsA5hLmUq0KtbiprJVIZIPYCPoRDC10WsS4BqA55hLmkq7vMZfC8PMahTXSqzL85pJfZeLKjB631ysFxphL9nVUNEMKdVbEoMnj8G6SNnNBAAIxEsBcwlyqVbcxfrFc60IDD2JjEhioIhxaKODxm0s6eC2PpvYx2bQlqNkDUH9a+nbj0caOrTayN200a0A3Ft6uUauZx/EwtlclBcaYS/Z1VDRDCnVWxKDJ4/BukjZzQQACMRLAXMJcqlW3FjeVtRKJbBAbEz+CoYVOixjXADTHVBlY9dRD/TVA0/9wr8/deiTaWBOuH9+bNpo1AHOpfh1oR3qrI+31eByfAmPMpfYrL4U6a5/yigzg7UkNcoEABDwSwFzCXKpVlxY3lbUSiWwQGxM/gqGFTosY1wA0x1zCXNL1fWe0pv/pwzAaWERBGwuqYWJ600azBmAuhamJOlG81VGda/A+JgXGmEvtV2EKddY+ZcwlTxqQCwQg4JsA5lIXfb75zW/KVVddJd/73vd8K9hSdhY3lS1dSqPTshFsFPcaJ0MLnRYxrgFojrmEuaTre8ylMPy8RmGN9KoMv7nkV5m4MqPH7fVKgTHmkn0dFc2QQp0VMWjyOLybpM1cEIBAjAQwl2JUzUHOMX6x7ACbsDHxoAJfsIdQIcY1gP6j9jGXQnS/iKb/6cMwGlhEQRsLqmFietNGswZ0I+LtGsMo5ysKjO31SIEx5pJ9HRXNkEKdFTFo8ji8m6TNXBCAQIwEkjCXXn31Vbn77rtl4cKF8uKLL8qyZctW06qvr09mzpwZo4at5GxxU9nKhTQ8KRuThoGvYTq00GkR4xqA5phLmEu6vu+M1vQ/fRhGA4soaGNBNUxMb9po1gDMpTA1USeKtzqqcw3ex6TAGHOp/SpMoc7ap7wiA3h7UoNcIAABjwR63lz61a9+JR/72MfkD3/4gyxfvryrBpm5tGjRIo8auczJ4qbS5YUGToqNSWCginBooYAnuicXdDPXH43mmEuYS/X7Z+BIzR6APgyjgUUUtLGgGiamN200awDmUpiaqBPFWx3VuQbvY1JgjLnUfhWmUGftU8Zc8qQBuUAAAr4J9Ly5dPTRR8sDDzwgH//4x2XChAkyevRoWWuttXyrEkF2FjeVEVy2OkU2gmqEwQKghQ5ljGsAmmMuYS7p+r4zWtP/9GEYDSyioI0F1TAxvWmjWQMwl8LURJ0o3uqozjV4H5MCY8yl9qswhTprnzLmkicNyAUCEPBNoOfNpW233Vbe//73yxe+8AXfSkSWncVNZWQIaqXLRrAWNpNBaKHDGuMagOaYS5hLur7HXArDz2sU1kivyoi73+y02ANQf/b1B2MYhyCAuRSCoi4GvazjV3U0vKsS43wIQCA1Aj1vLu200075E0uf+MQnUtPW9HotbipNE3YSnI2JEyHE3xc1fsiUyyTGNYD+w1zCXCrX30VnafqfPiyi295xtGmPfdHM3rTRrAHdrtXbNRZpEuNxGNurlgJjzCX7OiqaIYU6K2LQ5HF4N0mbuSAAgRgJ9Ly5dPbZZ8vChQvlpptukux3lfiEIWBxUxkmM99R2Jj40QctdFrEuAagOeYS5pKu7zujNf1PH4bRwCIK2lhQDRPTmzaaNQBzKUxN1InirY7qXIP3MSkwxlxqvwpTqLP2Ka/IAN6e1CAXCEDAI4GeN5defPFFOeqoo+Qd73iHnHHGGTJmzBiPOkSXk8VNZXQQaiTMxqQGNKMhaKEDG+MagOaYS5hLur7HXArDz2sU1kivyvh72tpiD0D92dcfjGEcggDmUgiKuhj0so5f1dHwrkqM8yEAgdQI9Ly5tMcee8hrr70mzz77bK7t+uuvLyNHjlxN5+yppgULFqSmf+3rtbiprJ1MRAPZmPgRCy10WsS4BqA55hLmkq7vMZfC8PMahTXSqzKYS36ViSszetxerxQYYy7Z11HRDCnUWRGDJo/Du0nazAUBCMRIoOfNpd133720Lt/73vdKn5v6iTF+sexBMzYmHlTgC/YQKsS4BtB/1D7mUojuF9H0P30YRgOLKGhjQTVMTG/aaNaAbkS8XWMY5XxFgbG9Hikwxlyyr6OiGVKosyIGTR6Hd5O0mQsCEIiRQM+bSzGKEkPOFjeVMVy3Nkc2JlqC4cajhY5ljGsAmmMuYS7p+r4zWtP/9GEYDSyioI0F1TAxvWmjWQMwl8LURJ0o3uqozjV4H5MCY8yl9qswhTprn/KKDODtSQ1ygQAEPBLAXPKoSgQ5WdxURnDZ6hTZmKgRBguAFjqUMa4BaI65hLmk63vMpTD8vEZhjfSqDK/F86tMXJnR4/Z6pcAYc8m+jopmSKHOihg0eRzeTdJmLghAIEYCmEsxquYg5xi/WHaATdiYeFCBL9hDqBDjGkD/UfuYSyG6n9fihaHoLwprpD9NOhl508ZiD+DtGv1WQ/3MYFyfXdmRKTDGXCpbDXbnpVBndvSqR4Z3dWaMgAAE0iLQc+bSjBkzpK+vT0499VQZNWqUZP9e5pONmTlzZplTOUd0XyylDJCNiR/10UKnhcUXS7qMikejOeYS5lJxn5Q5Q9P/9GEZwu2cgzbtcC8zqzdtNGtAt+v1do1ldIntHBjbK5YCY8wl+zoqmiGFOiti0ORxeDdJm7kgAIEYCfScuTRu3LjcXLr99ttl8803l+zfy3yyMYsWLSpzKudgLtWuATYmtdEFH4gWOqQWXyzpMioejeaYS5hLxX1S5gxN/9OHZQi3cw7atMO9zKzetNGsAZhLZRS3OcdbHdlcZbtRU2CMudRujWWzp1Bn7VNekQG8PalBLhCAgEcCPWcueYTcizlZ3FT2IqdVr4mNiR+V0UKnRYxrAJpjLmEu6fq+M1rT//RhGA0soqCNBdUwMb1po1kDMJfC1ESdKN7qqM41eB+TAmPMpfarMIU6a58y5pInDcgFAhDwTQBzybc+brOzuKl0e7EBE2MjGBCmMhRa6ADGuAagOeYS5pKu7zGXwvDzGoU10qsy/v4rdYs9APVnX38whnEIAphLISjqYtDLOn5VR8O7KjHOhwAEUiOAuZSa4oGu1+KmMlBqrsOwMfEjD1rotIhxDUBzzCXMJV3fYy6F4ec1CmukV2Uwl/wqE1dm9Li9Xikwxlyyr6OiGVKosyIGTR6Hd5O0mQsCEIiRQDLm0uLFi+UnP/mJPPPMM7J06dLVtMp+c+nEE0+MUcNWco7xi+VWQK0yKRsTDyrwBXsIFWJcA+g/ah9zKUT3i2j6nz4Mo4FFFLSxoBompjdtNGtANyLerjGMcr6iwNhejxQYYy7Z11HRDCnUWRGDJo/Du0nazAUBCMRIoOfNpeXLl8unP/1puf766+WNN96QzETK/q7z6fx79ueiRYti1LCVnC1uKlu5kIYnZWPSMPA1TIcWOi1iXAPQHHMJc0nX953Rmv6nD8NoYBEFbSyohonpTRvNGoC5FKYm6kTxVkd1rsH7mBQYYy61X4Upa4TtkQAAIABJREFU1Fn7lFdkAG9PapALBCDgkUDPm0vf/OY3ZdasWTJx4kQ57LDD5OCDD5YpU6bIBz/4QXn44YfliiuukB122EHOPPNM2WSTTTxq5DIni5tKlxcaOCk2JoGBKsKhhQKe6J5c0M1cfzSaYy5hLtXvn4EjNXsA+jCMBhZR0MaCapiY3rTRrAGYS2Fqok4Ub3VU5xq8j0mBMeZS+1WYQp21TxlzyZMG5AIBCPgm0PPm0gEHHCBrr7223HTTTbkS48aNk6lTp+b/ZJ8nn3xSPvzhD8tJJ52Um058yhGwuKksN3PcZ7ER9KMfWui0iHENQHPMJcwlXd93Rmv6nz4Mo4FFFLSxoBompjdtNGsA5lKYmqgTxVsd1bkG72NSYIy51H4VplBn7VPGXPKkAblAAAK+CfS8ubTNNtvkTy2dffbZuRJ///d/L8cee6yccsop/cqcdtpp+Svxbr/9dt9qOcrO4qbS0eWZpcJG0Axt5cBoURnZSgNiXAPQHHMJc0nX95hLYfh5jcIa6VUZEW/aWOwBvF2j32qonxmM67MrOzIFxphLZavB7rwU6syOXvXI8K7OjBEQgEBaBHreXHrve9+bP5k0ffr0XNns3/faay+54IIL+pX+3Oc+J1dffbX84he/SEt9xdVa3FQq0olmKBsTP1KhhU6LGNcANMdcwlzS9T3mUhh+XqOwRnpVBnPJrzJxZUaP2+uVAmPMJfs6KpohhTorYtDkcXg3SZu5IACBGAn0vLmU/cZS9ltKX/ziF3N9jjzySFm8eLHccccdMmzYsPzvMvPphRdekLvvvruWhq+//rpcfvnlcuONN8qzzz4rY8eOlSOOOEImTZokfX19XWM+//zz8p3vfEfuvfdeeeKJJ2TJkiXy9re/Xfbff//8FX3rrLNO/9gq53au8yc/+clqc6+11lqycOHCWtc5cFCMXyyrLzpAADYmASAGCoEWOpAxrgFo/qbmcICDrvt1v7lG/Wnp241HGzu22sjetLHYA3i7Rq1mHsfD2F6VFBhjLtnXUdEMKdRZEYMmj8O7SdrMBQEIxEig582liy66SL797W/LD3/4w9xMuuWWW+SMM86Qf/iHf5CddtpJfv7zn8vPfvYzOf744/PfXarzOeecc2T+/Pn56/e23nrrfK4777xTpk2b1v/bToPFzUylE088UXbZZZc8l5EjR8pDDz0kt956q4wfP17mzZsnmRmUfaqcm52fmWiPPfZY/+sAO/MPGTJEst+h0n4sbiq1OcUwno2JH5XQQqdFjGsAmmOqDKx66qH+GqDpf7jX5249Em2sCdeP700bzRrQjYK3a6yvlt+RMLbXJgXGmEv2dVQ0Qwp1VsSgyePwbpI2c0EAAjES6Hlz6emnn5b77rsvfxXeX//1X+cafeUrX5GvfvWr8vLLL+dPBx166KH5a/PWXnvtyhpmv9U0YcIEOfroo/tfvZcFOfnkk+Wee+7J/xk9evSgcX/3u9/lf589rTTwc8kll8hll10mc+bMyfPOPlXOzc7PzKXf/OY38oMf/KDyNZUZYHFTWWbe2M9hY+JHQbTQaRHjGoDmb2oOBzjoup8nl7T8vI5nbfCqjL9122IPQP3Z1x+MYRyCAOZSCIq6GPSyjl/V0fCuSozzIQCB1Aj0vLnUTdBly5bJ//7v/8pGG20k2dM8dT8XX3xx/kq87MmijTfeuD/Mww8/LIcffrice+65+Z9VPtkTRwceeGD+JNUJJ5ywxqHdzu2YS1ler7zyiowYMWKNr+irkl92rsVNZdUcYjyfjYkf1dBCp0WMawCaY6oMrHrqof4aoOl/uNfnbj0SbawJ14/vTRvNGtCNgrdrrK+W35EwttcmBcaYS/Z1VDRDCnVWxKDJ4/BukjZzQQACMRLoeXNpxowZMm7cuPw3jCw+2RNLjz/+eP4qvIGfpUuXyjbbbCPZbz5dcMEFlabOnrQ65phj5LzzzpPDDjtsjWO7nZuZS9nGb+jQofLqq6/mr9zbe++95fTTT+9/gqtSUqucnMVevnx5blrxKU8gM/qyz/Dhw8sP4kwTAr2iRba+tfGJcQ3oFc21esPhTYK9wqGNNUDT/73CXduHHsejjUdV1rxetdH/WUaaNaAbZerPvv5g3HuM21gDVu1/6sq+rladAebNMvfKu43+b5Y8s0EAArEQ6HlzKTN4Jk+eLKeddpqJJvvvv3/+W0433XTTavF33nln2XLLLeXKK68sPfcbb7yRG2GPPvqoLFiwQEaNGtV17JrOzUy1MWPGyBZbbJGbQA888IDccMMN+Sv4sj/XX3/90jkNdqLFTaUqoUgGe92YRIIvaJq9okVbm8oY14Be0VzbCHDAXNLWkKb/qT8tfbvxaGPHVhu5mza9tAeg/rRVUjwexsWMtGc0zbiNNQBzSVsl+vFN15k+47gjeOXdRv/HrSTZQwACVgR63lzKnhzafPPNZfbs2SYM99xzz9wAuu6661aLv9tuu+Vmzrx580rP3XnN3jnnnJP/btKaPlXOzeLMnz9fsrhTp06VadOmlc6pm7mU/f348eNVcVIbzCPVfhRHC50WFq/E0WVUPBrN32QEBzgUd8uaz9D0P/WnpW83Hm3s2Goje9NGswZ0Y+HtGrWaeRwPY3tVUmDMa/Hs66hohhTqrIhBk8fh3SRt5oIABGIk0PPm0u233y7ZUzyZwbP11lsH1yjkk0vf+ta35Pzzz89fhZe9Em9NnyrnDoyz44475mbbYGZYFTgWN5VV5o/1XDYmfpRDC50WMa4BaI6pMrDqqYf6a4Cm/+Fen7v1SLSxJlw/vjdtNGsA5lL9OtCO9FZH2uvxOD4FxphL7VdeCnXWPuUVGcDbkxrkAgEIeCTQ8+bSzTffLLfccos8+OCDss8++8jf//3f57851NfXt5oeEyZMqKxR0W8uHXTQQTJz5szCuNlr9c466yzZb7/95KKLLpIhQ4Z0HVPl3FWDZPksWbJE7rrrrsKc1nSCxU2lKqFIBrMx8SMUWui0iHENQHPMJcwlXd93Rmv6nz4Mo4FFFLSxoBompjdtNGsA5lKYmqgTxVsd1bkG72NSYIy51H4VplBn7VPGXPKkAblAAAK+CfSkuZQ9qZS9rm6PPfaQ7D2kmZGU/e7QwM9Acyk7lv37okWLKquVvW7viiuukHvvvVc23njj/vEPP/ywHH744fLJT35SJk2atMa4t912m5x++umy6667ypw5c2To0KFdz69y7qpBst9oyp5cete73iXXXntt5WsdOMDiplKVUCSD2Qj6EQotdFrEuAag+ZuawwEOuu4X0fQ/9aelbzcebezYaiN700azBnRj4e0atZp5HA9je1VSYIy5ZF9HRTOkUGdFDJo8Du8maTMXBCAQI4GeNJcyQyn7XaHsn+wpn8GeUhpMrOypnqqfhQsXSjYue4Jp+vTp/cNPPvlkWbBggdxzzz0yZswYyX4EcPHixbLhhhvKRhtt1H9eds5JJ50kO+ywQ25SDRs2rGsKZc996aWXZO2115Z11llnpVhXXnmlfP7zn5dTTz1VjjvuuKqXutL5FjeVqoQiGczGxI9QaKHTIsY1AM0xVQZWPfVQfw3Q9D/c63O3Hok21oTrx/emjWYN6EbB2zXWV8vvSBjba5MCY8wl+zoqmiGFOiti0ORxeDdJm7kgAIEYCfS8udSEKNnr7DITa+LEibLVVlvJ/fffL3fccUdubk2bNi1PIXst3+TJk1f6u0ceeSR/qikzgs4880wZPnz4Suluuummst122+V/V+XcbK5TTjklf8VeFiMz17K/u/vuu/MnubKnltZbbz0VGoubSlVCkQxmY+JHKLTQaRHjGoDmmEuYS7q+74zW9D99GEYDiyhoY0E1TExv2mjWAMylMDVRJ4q3OqpzDd7HpMAYc6n9KkyhztqnvCIDeHtSg1wgAAGPBDCXAqjy2muvyeWXX54bTM8884yMHTs2N42OPPLI/qemBjOXsvOzV/h1+2RPRM2aNSs/XOXcp556Kv/dpv/+7/+W5557TpYtWyabbLKJ7L333nLsscfKiBEj1FdtcVOpTiqCAGxM/IiEFjotYlwD0BxzCXNJ1/eYS2H4eY3CGulVGX+vM7XYA1B/9vUHYxiHIIC5FIKiLga9rONXdTS8qxLjfAhAIDUCmEupKR7oei1uKgOl5joMGxM/8qCFTosY1wA0x1zCXNL1PeZSGH5eo7BGelUGc8mvMnFlRo/b65UCY8wl+zoqmiGFOiti0ORxeDdJm7kgAIEYCfSsuZQ9PZT9U/aTvTpu7ty5ZU9P/rwYv1j2IBobEw8q8AV7CBViXAPoP2ofcylE94to+p8+DKOBRRS0saAaJqY3bTRrQDci3q4xjHK+osDYXo8UGGMu2ddR0Qwp1FkRgyaPw7tJ2swFAQjESKBnzaWqYmTm0qJFi6oOS/Z8i5vKFGCyMfGjMlrotIhxDUBzzCXMJV3fd0Zr+p8+DKOBRRS0saAaJqY3bTRrAOZSmJqoE8VbHdW5Bu9jUmCMudR+FaZQZ+1TXpEBvD2pQS4QgIBHAj1rLk2ZMkUmT55ciXmVJ50qBe7Bky1uKnsQ02qXxMbEj8poodMixjUAzTGXMJd0fY+5FIaf1yiskV6V4bV4fpWJKzN63F6vFBhjLtnXUdEMKdRZEYMmj8O7SdrMBQEIxEigZ82lqVOnSvYPHxsCMX6xbEOiWlQ2JtV4WZ6NFjq6Ma4BaI65hLmk63vMpTD8vEZhjfSqDOaSX2Xiyowet9crBcaYS/Z1VDRDCnVWxKDJ4/BukjZzQQACMRLAXIpRNQc5x/jFsgNswsbEgwp8wR5ChRjXAPqP2sdcCtH9/OZSGIr+orBG+tOkk5E3bSz2AN6u0W811M8MxvXZlR2ZAmPMpbLVYHdeCnVmR696ZHhXZ8YICEAgLQKYS2npHexqLW4qgyXnOBAbEz/ioIVOixjXADTHXMJc0vV9Z7Sm/+nDMBpYREEbC6phYnrTRrMGdCPi7RrDKOcrCozt9UiBMeaSfR0VzZBCnRUxaPI4vJukzVwQgECMBDCXYlTNQc4WN5UOLss8BTYm5ohLT4AWpVENemKMawCaYy5hLun6HnMpDD+vUVgjvSrDa/H8KhNXZvS4vV4pMMZcsq+johlSqLMiBk0eh3eTtJkLAhCIkUBPmksxChFbzjF+seyBMRsTDyrwBXsIFWJcA+g/ah9zKUT381q8MBT9RWGN9KdJJyNv2ljsAbxdo99qqJ8ZjOuzKzsyBcaYS2Wrwe68FOrMjl71yPCuzowREIBAWgQwl9LSO9jVWtxUBkvOcSA2Jn7EQQudFjGuAWiOuYS5pOv7zmhN/9OHYTSwiII2FlTDxPSmjWYN6EbE2zWGUc5XFBjb65ECY8wl+zoqmiGFOiti0ORxeDdJm7kgAIEYCWAuxaiag5wtbiodXJZ5CmxMzBGXngAtSqMa9MQY1wA0x1zCXNL1PeZSGH5eo7BGelWG1+L5VSauzOhxe71SYIy5ZF9HRTOkUGdFDJo8Du8maTMXBCAQIwHMpRhVc5BzjF8sO8AmbEw8qMAX7CFUiHENoP+ofcylEN3Pa/HCUPQXhTXSnyadjLxpY7EH8HaNfquhfmYwrs+u7MgUGGMula0Gu/NSqDM7etUjw7s6M0ZAAAJpEcBcSkvvYFdrcVMZLDnHgdiY+BEHLXRaxLgGoDnmEuaSru87ozX9Tx+G0cAiCtpYUA0T05s2mjWgGxFv1xhGOV9RYGyvRwqMMZfs66hohhTqrIhBk8fh3SRt5oIABGIkgLkUo2oOcra4qXRwWeYpsDExR1x6ArQojWrQE2NcA9AccwlzSdf3mEth+HmNwhrpVRlei+dXmbgyo8ft9UqBMeaSfR0VzZBCnRUxaPI4vJukzVwQgECMBDCXYlTNQc4xfrHsABuvxfMgwls5sEnUiRHjGoDmmEuYS7q+x1wKw89rFNZIr8pgLvlVJq7M6HF7vVJgjLlkX0dFM6RQZ0UMmjwO7yZpMxcEIBAjAcylGFVzkHOMXyw7wIa55EEEzKUgKsS4BnBjgLmEuRSk/UXT//RhGA0soqCNBdUwMb1po1kDuhHxdo1hlPMVBcb2eqTAGHPJvo6KZkihzooYNHkc3k3SZi4IQCBGAphLMarmIGeLm0oHl2WeAhsTc8SlJ0CL0qgGPTHGNQDNMZcwl3R93xmt6X/6MIwGFlHQxoJqmJjetNGsAZhLYWqiThRvdVTnGryPSYEx5lL7VZhCnbVPeUUG8PakBrlAAAIeCWAueVQlgpwsbiojuGx1imxM1AiDBUALHcoY1wA0x1zCXNL1PeZSGH5eo7BGelWG1+L5VSauzOhxe71SYIy5ZF9HRTOkUGdFDJo8Du8maTMXBCAQIwHMpRhVc5BzjF8sO8DGa/E8iPBWDmwSdWLEuAagOeYS5pKu7zGXwvDzGoU10qsymEt+lYkrM3rcXq8UGGMu2ddR0Qwp1FkRgyaPw7tJ2swFAQjESABzKUbVHOQc4xfLDrBhLnkQAXMpiAoxrgHcGGAuYS4FaX9+cykMRndRWCPdSdKfkDdtLPYA3q7RbzXUzwzG9dmVHZkCY8ylstVgd14KdWZHr3pkeFdnxggIQCAtAphLaekd7GotbiqDJec4EBsTP+KghU6LGNcANMdcwlzS9X1ntKb/6cMwGlhEQRsLqmFietNGswZ0I+LtGsMo5ysKjO31SIEx5pJ9HRXNkEKdFTFo8ji8m6TNXBCAQIwEMJdiVM1BzhY3lQ4uyzwFNibmiEtPgBalUQ16YoxrAJpjLmEu6foecykMP69RWCO9KsNr8fwqE1dm9Li9Xikwxlyyr6OiGVKosyIGTR6Hd5O0mQsCEIiRAOZSjKo5yDnGL5YdYOO1eB5EeCsHNok6MWJcA9AccwlzSdf3mEth+HmNwhrpVRnMJb/KxJUZPW6vVwqMMZfs66hohhTqrIhBk8fh3SRt5oIABGIkgLkUo2oOco7xi2UH2DCXPIiAuRREhRjXAG4MMJcwl4K0P7+5FAajuyiske4k6U/ImzYWewBv1+i3GupnBuP67MqOTIEx5lLZarA7L4U6s6NXPTK8qzNjBAQgkBYBzKW09A52tRY3lcGScxyIjYkfcdBCp0WMawCaYy5hLun6vjNa0//0YRgNLKKgjQXVMDG9aaNZA7oR8XaNYZTzFQXG9nqkwBhzyb6OimZIoc6KGDR5HN5N0mYuCEAgRgKYSzGq5iBni5tKB5dlngIbE3PEpSdAi9KoBj0xxjUAzTGXMJd0fY+5FIaf1yiskV6V4bV4fpWJKzN63F6vFBhjLtnXUdEMKdRZEYMmj8O7SdrMBQEIxEgAcylG1RzkHOMXyw6w8Vo8DyK8lQObRJ0YMa4BaI65hLmk63vMpTD8vEZhjfSqDOaSX2Xiyowet9crBcaYS/Z1VDRDCnVWxKDJ4/BukjZzQQACMRLAXIpRNQc5x/jFsgNsmEseRMBcCqJCjGsANwaYS5hLQdqf31wKg9FdFNZId5L0J+RNG4s9gLdr9FsN9TODcX12ZUemwBhzqWw12J2XQp3Z0aseGd7VmTECAhBIiwDmUlp6B7tai5vKYMk5DsTGxI84aKHTIsY1AM0xlzCXdH3fGa3pf/owjAYWUdDGgmqYmN600awB3Yh4u8YwyvmKAmN7PVJgjLlkX0dFM6RQZ0UMmjwO7yZpMxcEIBAjAcylGFVzkLPFTaWDyzJPgY2JOeLSE6BFaVSDnhjjGoDmmEuYS7q+x1wKw89rFNZIr8rwWjy/ysSVGT1ur1cKjDGX7OuoaIYU6qyIQZPH4d0kbeaCAARiJIC5FKNqDnKO8YtlB9h4LZ4HEd7KgU2iTowY1wA0x1zCXNL1PeZSGH5eo7BGelUGc8mvMnFlRo/b65UCY8wl+zoqmiGFOiti0ORxeDdJm7kgAIEYCWAuxaiag5xj/GLZATbMJQ8iYC4FUSHGNYAbA8wlzKUg7c9vLoXB6C4Ka6Q7SfoT8qaNxR7A2zX6rYb6mcG4PruyI1NgjLlUthrszkuhzuzoVY8M7+rMGAEBCKRFAHMpLb2DXa3FTWWw5BwHYmPiRxy00GkR4xqA5phLmEu6vu+M1vQ/fRhGA4soaGNBNUxMb9po1oBuRLxdYxjlfEWBsb0eKTDGXLKvo6IZUqizIgZNHod3k7SZCwIQiJEA5lKMqjnI2eKm0sFlmafAxsQccekJ0KI0qkFPjHENQHPMJcwlXd9jLoXh5zUKa6RXZXgtnl9l4sqMHrfXKwXGmEv2dVQ0Qwp1VsSgyePwbpI2c0EAAjESwFyKUTUHOcf4xbIDbLwWz4MIb+XAJlEnRoxrAJpjLmEu6foecykMP69RWCO9KoO55FeZuDKjx+31SoEx5pJ9HRXNkEKdFTFo8ji8m6TNXBCAQIwEMJdiVM1BzjF+sewAG+aSBxEwl4KoEOMawI0B5hLmUpD25zeXwmB0F4U10p0k/Ql508ZiD+DtGv1WQ/3MYFyfXdmRKTDGXCpbDXbnpVBndvSqR4Z3dWaMgAAE0iKAuZSW3sGu1uKmMlhyjgOxMfEjDlrotIhxDUBzzCXMJV3fd0Zr+p8+DKOBRRS0saAaJqY3bTRrQDci3q4xjHK+osDYXo8UGGMu2ddR0Qwp1FkRgyaPw7tJ2swFAQjESABzKUbVHORscVPp4LLMU2BjYo649ARoURrVoCfGuAagOeYS5pKu7zGXwvDzGoU10qsyvBbPrzJxZUaP2+uVAmPMJfs6KpohhTorYtDkcXg3SZu5IACBGAlgLsWomoOcY/xi2QE2XovnQYS3cmCTqBMjxjUAzTGXMJd0fY+5FIaf1yiskV6VwVzyq0xcmdHj9nqlwBhzyb6OimZIoc6KGDR5HN5N0mYuCEAgRgKYSzGq5iDnGL9YdoANc8mDCJhLQVSIcQ3gxgBzCXMpSPvzm0thMLqLwhrpTpL+hLxpY7EH8HaNfquhfmYwrs+u7MgUGGMula0Gu/NSqDM7etUjw7s6M0ZAAAJpEcBcSkvvYFdrcVMZLDnHgdiY+BEHLXRaxLgGoDnmEuaSru87ozX9Tx+G0cAiCtpYUA0T05s2mjWgGxFv1xhGOV9RYGyvRwqMMZfs66hohhTqrIhBk8fh3SRt5oIABGIkgLkUo2oOcra4qXRwWeYpsDExR1x6ArQojWrQE2NcA9AccwlzSdf3mEth+HmNwhrpVRlei+dXmbgyo8ft9UqBMeaSfR0VzZBCnRUxaPI4vJukzVwQgECMBDCXYlTNQc4xfrHsANv/z965QF9VVft//nyQiNmVS5qgPa7dZOTFB6lpXm/mKzNTMCXjZWpZKSjeVNS8aqM0K6GhYje8pleptAR8lKKFoiY+Mq6JDki7vnpgYfYX83HDB/+xD/5+/IDfYT/mmnvPddbnjMEw2XvNtdbnO+ds7fP1nMPX4nkQ4c01cEjUiRFjD0BzzCXMJV3dYy6F4ec1Cj3SqzKYS36ViWtl1Li9Xikwxlyyz6O8GVLIszwGdV6Hd520mQsCEIiRAOZSjKo5WHOMbyw7wIa55EEEzKUgKsTYA3gwwFzCXApS/vzmUhiM7qLQI91J0rMgb9pYnAG87dFvNlRfGYyrsys6MgXGmEtFs8HuvhTyzI5e+cjwLs+MERCAQFoEMJfS0jvYbi0eKoMtznEgDiZ+xEELnRYx9gA0x1zCXNLVffdoTf1Th2E0sIiCNhZUw8T0po2mB7Qj4m2PYZTzFQXG9nqkwBhzyT6P8mZIIc/yGNR5Hd510mYuCEAgRgKYSzGq5mDNFg+VDrZlvgQOJuaIC0+AFoVR9XljjD0AzTGXMJd0dY+5FIaf1yj0SK/K8LV4fpWJa2XUuL1eKTDGXLLPo7wZUsizPAZ1Xod3nbSZCwIQiJEA5lKMqjlYc4xvLDvAxtfieRDhzTVwSNSJEWMPQHPMJcwlXd1jLoXh5zUKPdKrMphLfpWJa2XUuL1eKTDGXLLPo7wZUsizPAZ1Xod3nbSZCwIQiJEA5lIA1V577TWZPn26zJo1S5599lkZMmSIjB07VsaMGSNdXV1tZ3j++efluuuuk3nz5snjjz8uL7/8smy99dZy0EEHyZFHHilvectbesaWubd70L333isXXnihLF68WPr37y8f+chH5JRTTpGBAweqdx3jG8vqTQcIwMEkAMRAIdBCBzLGHoDmmEuYS7q6x1wKw89rFHqkV2Uwl/wqE9fKqHF7vVJgjLlkn0d5M6SQZ3kM6rwO7zppMxcEIBAjAcylAKqdeeaZcu2118qoUaNk++23l7vvvltuueUWmThxokyYMKHtDJmpdPzxx8see+whu+22m2yyySbywAMPyE9/+lMZPny4zJgxQ9Zff/3W+DL3Zvf/8pe/lKOOOkq23XZbOeyww+Svf/2rXH755TJ48GCZOXOmbLTRRqqdx/jGsmrDgQZzMAkEMkAYtNBBjLEHoDnmEuaSru4xl8Lw8xqFHulVGcwlv8rEtTJq3F6vFBhjLtnnUd4MKeRZHoM6r8O7TtrMBQEIxEgAc0mpWvapoBEjRsjRRx8tkydP7ok2adIkue2221p/Nt988z5n+f3vf9/6++zTSr1f2aeNvvOd78i0adNkv/32a10qc292f7amZcuWyU033SQbb7xxK8add94pxx57rJx++unymc98RrXzGN9YVm040GAOJoFABgiDFjqIMfYANF+pORzgoKt+EU39k39a+nbj0caOrTayN200PaAdC2971GrmcTyM7VVJgTHmkn0e5c2QQp7lMajzOrzrpM1cEIBAjAQwl5SqTZ06tfWVeNkni7JPBXW/FixYIKNHj5azzz679c8yr0cffVQOPvhgOfHEE+W4445b59C+7n3yySflgAMO6POTU/vvv7+87W1va33SSvOyeKjUrCeWsRxM/CiFFjotYuwBaI6p0jvryYfqPUBT/3Cvzt16JNrFEPKyAAAgAElEQVRYE64e35s2mh6AuVQ9D7QjveWRdj8ex6fAGHOp+cxLIc+ap7xqBfD2pAZrgQAEPBLAXFKqkn1i6bHHHmt9FV7v1/Lly2WHHXaQQw89VM4999xSs/ziF7+Qz372s/KVr3xFjjjiiHWO7even/zkJ3LyySfLZZddJnvuuedq47O/v/XWW+XXv/51z1fulVrcmzdnh8oVK1bIgAEDqgxPdswrr7zS2nv2G1i8miXQKVoMHTq0EZAx9oBO0VwrOBxWEuwUDk30AE39dwp3bR16HI82HlVZd79qov6zFWl6QDvK5J99/sG48xg30QPWrH/yyj6v1pwB5vUy98q7ifqvlzyzQQACsRDAXFIqddBBB0m/fv1k9uzZa0XafffdZbvttmuZPEVfb7zxhhx55JHy8MMPy9y5c2XQoEFth7a793vf+55885vflBtvvLH1m0u9X9nfZ9fnz5+/zth567V4qMybsxOuez2YdALbsnvoFC2aOlTG2AM6RfOyuc4Dad/EOiUfmugBmvrvFO7aOvQ4Hm08qoK55FeV+FZGjdtrVjdjD2eAuvdsr6L/GWBer0ZeeTdR//WSZzYIQCAWAphLSqX23XfflklzzTXXrBVpr732av2e0owZMwrP0v01e2eeeaaMGzdunePa3XvJJZfIRRddJLfccou85z3vWS1G9+85Zb8FtdVWWxVe15o3WnwdRuXFRDSQj1T7EQstdFrE2APQfKXmcICDrvr5zSUtP6/j6Q1elfHXty3OAOSfff7BGMYhCPC1eCEo6mJQyzp+ZUfDuywx7ocABFIjgLmkVDzkJ5e+//3vy1e/+tXWV+FlX4m3rte67q3rk0vZ+oYPH64kmNZwDiZ+9EYLnRYWbyzpVpQ/Gs0xVXpnCfmQXzPt7tDUP9yrc7ceiTbWhKvH96aNpge0o+Btj9XV8jsSxvbapMAYc8k+j/JmSCHP8hjUeR3eddJmLghAIEYCmEtK1fJ+c2nkyJFy3nnn5c6Sfa3eGWecIQceeKBccMEFst5667Udk3dv3m8uZZ9oeuihh9S/uYS5lCvrWjdwMCnPzGoEWujIWryxpFtR/mg0x1zCXMqvkyJ3aOqfOixCuJl70KYZ7kVm9aaNpgdgLhVR3OYeb3lks8tmo6bAGHOp2RzLZk8hz5qnvGoF8PakBmuBAAQ8EsBcUqoyZcoUufTSS2XevHkyePDgnmgLFiyQ0aNHy1lnnSVjxoxZ5yw33XSTnHzyyfLhD39Ypk2bJhtssEHb+4vc+8QTT8jHPvYxmThxokyYMGG1WPvvv79suummMnPmTNXOLR4qVQuKZDAHEz9CoYVOixh7AJpjLmEu6eq+e7Sm/qnDMBpYREEbC6phYnrTRtMDMJfC5ESVKN7yqMoevI9JgTHmUvNZmEKeNU8Zc8mTBqwFAhDwTQBzSanPokWLJPt0UvYJpsmTJ/dEmzRpksydO1ey3zbaYostJPsRwCVLlshmm20mAwcO7Lkvu+fEE0+UXXbZpWVS9evXr+2Kytx7yCGHyAsvvCCZGbXxxhu3Yt55551y7LHHttaZrVfzsnio1KwnlrEcBP0ohRY6LWLsAWi+UnM4wEFX/fzmkpaf1/H0Bq/K+OvbFmcA8s8+/2AM4xAEMJdCUNTFoJZ1/MqOhndZYtwPAQikRgBzKYDi2dfZZV9VN2rUKBk2bJjMnz9f5syZ0/rUUPbpoex1//33y/jx41f7u4ULF7Y+1bThhhvKqaeeKv37919tNe985ztlp512av1dmXuz+++7776WgTR06FA5/PDD5bnnnpMrrriiZXTNmjVrrbnKYrB4qCy7hhjv52DiRzW00GkRYw9Ac0yV3llPPlTvAZr6h3t17tYj0caacPX43rTR9IB2FLztsbpafkfC2F6bFBhjLtnnUd4MKeRZHoM6r8O7TtrMBQEIxEgAcymAaq+++qpMnz69ZTAtXbpUhgwZ0jKNxo0bJ11dXW3Npez+008/ve0Ksk9EnX/++a3rZe7tDnjPPffIhRdeKIsXL26ZSXvttZeccsopMmjQIPWuLR4q1YuKIAAHEz8ioYVOixh7AJpjLmEu6eq+e7Sm/qnDMBpYREEbC6phYnrTRtMDMJfC5ESVKN7yqMoevI9JgTHmUvNZmEKeNU951Qrg7UkN1gIBCHgkgLnkUZUI1mTxUBnBttVL5GCiRhgsAFroUMbYA9AccwlzSVf3mEth+HmNQo/0qgxfi+dXmbhWRo3b65UCY8wl+zzKmyGFPMtjUOd1eNdJm7kgAIEYCWAuxaiagzXH+MayA2z81okHEd5cA4dEnRgx9gA0x1zCXNLVPeZSGH5eo9AjvSqDueRXmbhWRo3b65UCY8wl+zzKmyGFPMtjUOd1eNdJm7kgAIEYCWAuxaiagzXH+MayA2yYSx5EwFwKokKMPYAHA8wlzKUg5S+a+qcOw2hgEQVtLKiGielNG00PaEfE2x7DKOcrCozt9UiBMeaSfR7lzZBCnuUxqPM6vOukzVwQgECMBDCXYlTNwZotHiodbMt8CRxMzBEXngAtCqPq88YYewCaYy5hLunqvnu0pv6pwzAaWERBGwuqYWJ600bTAzCXwuRElSje8qjKHryPSYEx5lLzWZhCnjVPedUK4O1JDdYCAQh4JIC55FGVCNZk8VAZwbbVS+RgokYYLABa6FDG2APQHHMJc0lX95hLYfh5jUKP9KoMX4vnV5m4VkaN2+uVAmPMJfs8ypshhTzLY1DndXjXSZu5IACBGAlgLsWomoM1x/jGsgNsfC2eBxHeXAOHRJ0YMfYANMdcwlzS1T3mUhh+XqPQI70qg7nkV5m4VkaN2+uVAmPMJfs8ypshhTzLY1DndXjXSZu5IACBGAlgLsWomoM1x/jGsgNsmEseRMBcCqJCjD2ABwPMJcylIOXPby6FweguCj3SnSQ9C/KmjcUZwNse/WZD9ZXBuDq7oiNTYIy5VDQb7O5LIc/s6JWPDO/yzBgBAQikRQBzKS29g+3W4qEy2OIcB+Jg4kcctNBpEWMPQHPMJcwlXd13j9bUP3UYRgOLKGhjQTVMTG/aaHpAOyLe9hhGOV9RYGyvRwqMMZfs8yhvhhTyLI9BndfhXSdt5oIABGIkgLkUo2oO1mzxUOlgW+ZL4GBijrjwBGhRGFWfN8bYA9AccwlzSVf3mEth+HmNQo/0qgxfi+dXmbhWRo3b65UCY8wl+zzKmyGFPMtjUOd1eNdJm7kgAIEYCWAuxaiagzXH+MayA2x8LZ4HEd5cA4dEnRgx9gA0x1zCXNLVPeZSGH5eo9AjvSqDueRXmbhWRo3b65UCY8wl+zzKmyGFPMtjUOd1eNdJm7kgAIEYCWAuxaiagzXH+MayA2yYSx5EwFwKokKMPYAHA8wlzKUg5c9vLoXB6C4KPdKdJD0L8qaNxRnA2x79ZkP1lcG4OruiI1NgjLlUNBvs7kshz+zolY8M7/LMGAEBCKRFAHMpLb2D7dbioTLY4hwH4mDiRxy00GkRYw9Ac8wlzCVd3XeP1tQ/dRhGA4soaGNBNUxMb9poekA7It72GEY5X1FgbK9HCowxl+zzKG+GFPIsj0Gd1+FdJ23mggAEYiSAuRSjag7WbPFQ6WBb5kvgYGKOuPAEaFEYVZ83xtgD0BxzCXNJV/eYS2H4eY1Cj/SqDF+L51eZuFZGjdvrlQJjzCX7PMqbIYU8y2NQ53V410mbuSAAgRgJYC7FqJqDNcf4xrIDbHwtngcR3lwDh0SdGDH2ADTHXMJc0tU95lIYfl6j0CO9KoO55FeZuFZGjdvrlQJjzCX7PMqbIYU8y2NQ53V410mbuSAAgRgJYC7FqJqDNcf4xrIDbJhLHkTAXAqiQow9gAcDzCXMpSDlz28uhcHoLgo90p0kPQvypo3FGcDbHv1mQ/WVwbg6u6IjU2CMuVQ0G+zuSyHP7OiVjwzv8swYAQEIpEUAcyktvYPt1uKhMtjiHAfiYOJHHLTQaRFjD0BzzCXMJV3dd4/W1D91GEYDiyhoY0E1TExv2mh6QDsi3vYYRjlfUWBsr0cKjDGX7PMob4YU8iyPQZ3X4V0nbeaCAARiJIC5FKNqDtZs8VDpYFvmS+BgYo648ARoURhVnzfG2APQHHMJc0lX95hLYfh5jUKP9KoMX4vnV5m4VkaN2+uVAmPMJfs8ypshhTzLY1DndXjXSZu5IACBGAlgLsWomoM1x/jGsgNsfC2eBxHeXAOHRJ0YMfYANMdcwlzS1T3mUhh+XqPQI70qg7nkV5m4VkaN2+uVAmPMJfs8ypshhTzLY1DndXjXSZu5IACBGAlgLsWomoM1x/jGsgNsmEseRMBcCqJCjD2ABwPMJcylIOXPby6FweguCj3SnSQ9C/KmjcUZwNse/WZD9ZXBuDq7oiNTYIy5VDQb7O5LIc/s6JWPDO/yzBgBAQikRQBzKS29g+12wYIFrVhdXV3BYqYQaMWKFXBzInSnaJHV4E477VQ71Rh7QKdorhUbDisJdgqHJnqApv47hbu2Dj2ORxuPqqy7XzVR/9mKND2gHWXyzz7/YNx5jJvoAWvWP3lln1drzgDzepl75d1E/ddLntkgAIFYCGAuxaKUs3VaPFQ62yLLgUAUBJo6VNIDokgPFpkAgSZ6APWfQGKxxSgINFH/VuZSFMBZJAScEWiiB3AGcJYELCdZAk3Uf7Kw2TgEILBOAphLJAgEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgEBhAphLhVFxIwQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAOYSOQABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIFCYAOZSYVTcCAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQggLlEDkAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCBQmgLlUGBU3QgACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIYC6RAxCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAoUJYC4VRsWNEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACmEvkAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAQGECmEuFUXFjbwIPPvhg61932mknwEAAAgkSoAckKDpbhsCbBKh/UgECaROgB6StP7tPmwD1n7b+7B4CEIAABCCwJgHMJXKiEoH/+Z//aY0bPnx4pfGpDvrNb37T2vrQoUNTReBm32ihkyLGHoDmKzWHAxx01S+iqX/yT0vfbjza2LHVRvamjaYHtGPhbY9azTyOh7G9KikwXrP+U9izfeaUmwHm5Xhp74a3liDjIQCBTieAudTpChvtz+Kh0miprsJyMPEjB1rotIixB6A5pkrvrCcfqvcATf3DvTp365FoY024enxv2mh6AOZS9TzQjvSWR9r9eByfAmPMpeYzL4U8a57yqhXA25MarAUCEPBIAHPJoyoRrMnioTKCbauXyMFEjTBYALTQoYyxB6A55hLmkq7uu0dr6p86DKOBRRS0saAaJqY3bTQ9AHMpTE5UieItj6rswfuYFBhjLjWfhSnkWfOUMZc8acBaIAAB3wQwl3zr43Z1Fg+VbjcbcGEcBAPCVIZCCx3AGHsAmmMuYS7p6h5zKQw/r1HokV6V8fd1phZnAPLPPv9gDOMQBDCXQlDUxaCWdfzKjoZ3WWLcDwEIpEYAcyk1xQPt1+KhMtDSXIfhYOJHHrTQaRFjD0BzzCXMJV3dYy6F4ec1Cj3SqzKYS36ViWtl1Li9Xikwxlyyz6O8GVLIszwGdV6Hd520mQsCEIiRAOZSjKo5WHOMbyw7wCYcTDyowBvsIVSIsQdQf+Q+5lKI6hfR1D91GEYDiyhoY0E1TExv2mh6QDsi3vYYRjlfUWBsr0cKjDGX7PMob4YU8iyPQZ3X4V0nbeaCAARiJIC5FKNqDtZs8VDpYFvmS+BgYo648ARoURhVnzfG2APQHHMJc0lX992jNfVPHYbRwCIK2lhQDRPTmzaaHoC5FCYnqkTxlkdV9uB9TAqMMZeaz8IU8qx5yqtWAG9ParAWCEDAIwHMJY+qRLAmi4fKCLatXiIHEzXCYAHQQocyxh6A5phLmEu6usdcCsPPaxR6pFdl+Fo8v8rEtTJq3F6vFBhjLtnnUd4MKeRZHoM6r8O7TtrMBQEIxEgAcylG1RysOcY3lh1g42vxPIjw5ho4JOrEiLEHoDnmEuaSru4xl8Lw8xqFHulVGcwlv8rEtTJq3F6vFBhjLtnnUd4MKeRZHoM6r8O7TtrMBQEIxEgAcylG1RysOcY3lh1gw1zyIALmUhAVYuwBPBhgLmEuBSl/fnMpDEZ3UeiR7iTpWZA3bSzOAN726Dcbqq8MxtXZFR2ZAmPMpaLZYHdfCnlmR698ZHiXZ8YICEAgLQKYS2npHWy3Fg+VwRbnOBAHEz/ioIVOixh7AJpjLmEu6eq+e7Sm/qnDMBpYREEbC6phYnrTRtMD2hHxtscwyvmKAmN7PVJgjLlkn0d5M6SQZ3kM6rwO7zppMxcEIBAjAcylGFVzsGaLh0oH2zJfAgcTc8SFJ0CLwqj6vDHGHoDmmEuYS7q6x1wKw89rFHqkV2X4Wjy/ysS1MmrcXq8UGGMu2edR3gwp5Fkegzqvw7tO2swFAQjESABzKUbVHKw5xjeWHWDja/E8iPDmGjgk6sSIsQegOeYS5pKu7jGXwvDzGoUe6VUZzCW/ysS1MmrcXq8UGGMu2edR3gwp5Fkegzqvw7tO2swFAQjESABzKUbVHKw5xjeWHWDDXPIgAuZSEBVi7AE8GGAuYS4FKX9+cykMRndR6JHuJOlZkDdtLM4A3vboNxuqrwzG1dkVHZkCY8ylotlgd18KeWZHr3xkeJdnxggIQCAtAphLaekdbLcWD5XBFuc4EAcTP+KghU6LGHsAmmMuYS7p6r57tKb+qcMwGlhEQRsLqmFietNG0wPaEfG2xzDK+YoCY3s9UmCMuWSfR3kzpJBneQzqvA7vOmkzFwQgECMBzKUYVXOwZouHSgfbMl8CBxNzxIUnQIvCqPq8McYegOaYS5hLurrHXArDz2sUeqRXZfhaPL/KxLUyatxerxQYYy7Z51HeDCnkWR6DOq/Du07azAUBCMRIAHMpRtUcrDnGN5YdYONr8TyI8OYaOCTqxIixB6A55hLmkq7uMZfC8PMahR7pVRnMJb/KxLUyatxerxQYYy7Z51HeDCnkWR6DOq/Du07azAUBCMRIAHMpRtUcrDnGN5YdYMNc8iAC5lIQFWLsATwYYC5hLgUpf35zKQxGd1Hoke4k6VmQN20szgDe9ug3G6qvDMbV2RUdmQJjzKWi2WB3Xwp5ZkevfGR4l2fGCAhAIC0CmEtp6R1stxYPlcEW5zgQBxM/4qCFTosYewCaYy5hLunqvnu0pv6pwzAaWERBGwuqYWJ600bTA9oR8bbHMMr5igJjez1SYIy5ZJ9HeTOkkGd5DOq8Du86aTMXBCAQIwHMpRhVc7Bmi4dKB9syXwIHE3PEhSdAi8Ko+rwxxh6A5phLmEu6usdcCsPPaxR6pFdl+Fo8v8rEtTJq3F6vFBhjLtnnUd4MKeRZHoM6r8O7TtrMBQEIxEgAcymAaq+99ppMnz5dZs2aJc8++6wMGTJExo4dK2PGjJGurq51znDzzTfLHXfcIQsXLpSnnnpKNt98c7nrrrv6HFPm3izA66+/LldffbXMnDlTnnzySenXr59ss802cvzxx8uee+6p2nmMbyyrNhxoMAeTQCADhEELHcQYewCaYy5hLunqHnMpDD+vUeiRXpXBXPKrTFwro8bt9UqBMeaSfR7lzZBCnuUxqPM6vOukzVwQgECMBDCXAqh25plnyrXXXiujRo2S7bffXu6++2655ZZbZOLEiTJhwoR1zjBu3Dh55JFHZLvttmuZS+utt15bc6nMvW+88YaccMIJcuedd8rIkSNl2LBh8sorr8j//u//tv734Ycfrtp5jG8sqzYcaDAHk0AgA4RBCx3EGHsAmmMuYS7p6h5zKQw/r1HokV6VwVzyq0xcK6PG7fVKgTHmkn0e5c2QQp7lMajzOrzrpM1cEIBAjAQwl5SqLV68WEaMGCFHH320TJ48uSfapEmT5Lbbbmv9yT6N1O71zDPPtK6vv/76kplHTz/9dFtzqcy9V111lXzjG9+QK6+8UnbeeWflLtceHuMby8EhVAjIwaQCNKMhaKEDG2MPQHPMJcwlXd1jLoXh5zUKPdKrMphLfpWJa2XUuL1eKTDGXLLPo7wZUsizPAZ1Xod3nbSZCwIQiJEA5pJStalTp7a+Em/evHkyePDgnmgLFiyQ0aNHy9lnn936Z5FXnrnUO8a67s0+tbTPPvu0PqF00UUXSfbv2aeWBgwYUGQZhe6J8Y3lQhszvomDiTHgEuHRogSsPm6NsQegOeYS5pKu7jGXwvDzGoUe6VUZzCW/ysS1MmrcXq8UGGMu2edR3gwp5Fkegzqvw7tO2swFAQjESABzSala9omlxx57rPVVeL1fy5cvlx122EEOPfRQOffccwvNEspcyr767uMf/7icdNJJsnTpUpk9e3bLXNpyyy3lC1/4ghxxxBGF1rOum7JD5YoVK4IaVupFRRAg0yF79e/fP4LVdvYSO0WLoUOHNiJUjD2gUzTXCg6HlQQ7hUMTPUBT/53CXVuHHsejjUdV1t2vmqj/bEWaHtCOMvlnn38w7jzGTfSANeufvLLPqzVngHm9zL3ybqL+6yXPbBCAQCwEMJeUSh100EHSr1+/loGz5mv33Xdv/ZbSZZddVmiWUObS3Llz5fjjj5fNNtusZWJkhlL2qaUf/ehH8stf/lLOOussGTNmTKE1tbvJ4qFStaBIBns9mESCL+gyO0WLpg6VMfaATtFcWwhwwFzS5pCm/sk/LX278Whjx1YbuZ02nXQGIP+0WZI/Hsb5jLR31M24iR6AuaTNEv34uvNMv+K4I3jl3UT9x60kq4cABKwIYC4pye67774yaNAgueaaa9aKtNdee8nWW28tM2bMKDRLKHPphhtukFNPPVU23HBDmTNnTmsN2eu1116TzAx7/vnnW5+02mCDDQqtq6+bYvxKrMqbDTiQj1QHhKkMhRY6gDH2ADRfqTkc4KCr/pWfWshew4cPLx2K/CuNrLYBaFMb6tITedNG0wPabd7bHkuLFMEAGNuLlAJjvhbPPo/yZkghz/IY1Hkd3nXSZi4IQCBGAphLStU8fnLp1ltvlRNOOEF23XXXtYytiy++WKZNmyY33nijbLvttpV3b/FQWXkxEQ3kYOJHLLTQaRFjD0BzTJXeWU8+VO8BmvqHe3Xu1iPRxppw9fjetNH0AMyl6nmgHektj7T78Tg+BcaYS81nXgp51jzlVSuAtyc1WAsEIOCRAOaSUpW831waOXKknHfeeYVmCfXJpQcffLD1u0oHHnigfPvb315t7quvvlrOOecc+f73vy+77LJLoXX1dZPFQ2XlxUQ0kIOJH7HQQqdFjD0AzTGXMJd0dd89WlP/1GEYDSyioI0F1TAxvWmj6QGYS2FyokoUb3lUZQ/ex6TAGHOp+SxMIc+ap4y55EkD1gIBCPgmgLmk1GfKlCly6aWXyrx582Tw4ME90RYsWCCjR48u9ftGocyll156SXbbbTf5l3/5F8nMpN6vzGz67ne/KzfffLNss802lXdv8VBZeTERDeQg6EcstNBpEWMPQHPMJcwlXd1jLoXh5zUKPdKrMv6+ztTiDED+2ecfjGEcggDmUgiKuhjUso5f2dHwLkuM+yEAgdQIYC4pFV+0aJFkn07KPsE0efLknmiTJk2SuXPnym233SZbbLGFZD8CuGTJEtlss81k4MCBfc4aylzKgmdfi/fzn/9crrvuOun+ob9sDdmnmbq6ulrryv5Z9WXxUFl1LTGN42DiRy200GkRYw9Ac8wlzCVd3WMuheHnNQo90qsymEt+lYlrZdS4vV4pMMZcss+jvBlSyLM8BnVeh3edtJkLAhCIkQDmUgDVzjjjDJk9e7aMGjVKhg0bJvPnz5c5c+bIhAkTZOLEia0Z7r//fhk/fvxqf5f9/QMPPND6k71mzpwpy5Ytk2OOOab179knoUaMGNGzwjL3Pv3003L44Ye3DKRs3gEDBrTW+Nvf/lYuuugi2W+//VQ7j/GNZdWGAw3mYBIIZIAwaKGDGGMPQPOVmsMBDrrqF9HUP/mnpW83Hm3s2Goje9NG0wPasfC2R61mHsfD2F6VFBhjLtnnUd4MKeRZHoM6r8O7TtrMBQEIxEgAcymAaq+++qpMnz69Zd4sXbpUhgwZImPGjJHsk0jdnw5qZy5dfPHFMm3atD5Xseuuu8qMGTN6rpW5Nxv0+OOPywUXXNAyr5YvXy7vf//75fjjj5c999xTvWuLh0r1oiIIwMHEj0hoodMixh6A5pgqvbOefKjeAzT1D/fq3K1Hoo014erxvWmj6QGYS9XzQDvSWx5p9+NxfAqMMZeaz7wU8qx5yqtWAG9ParAWCEDAIwHMJY+qRLAmi4fKCLatXiIHEzXCYAHQQocyxh6A5phLmEu6uu8eral/6jCMBhZR0MaCapiY3rTR9ADMpTA5USWKtzyqsgfvY1JgjLnUfBamkGfNU8Zc8qQBa4EABHwTwFzyrY/b1Vk8VLrdbMCFcRAMCFMZCi10AGPsAWiOuYS5pKt7zKUw/LxGoUd6Vcbf15lanAHIP/v8gzGMQxDAXApBUReDWtbxKzsa3mWJcT8EIJAaAcyl1BQPtF+Lh8pAS3MdhoOJH3nQQqdFjD0AzTGXMJd0dY+5FIaf1yj0SK/KYC75VSaulVHj9nqlwBhzyT6P8mZIIc/yGNR5Hd510mYuCEAgRiY+MoQAACAASURBVAKYSzGq5mDNMb6x7ACbcDDxoAJvsIdQIcYeQP2R+5hLIapfRFP/1GEYDSyioI0F1TAxvWmj6QHtiHjbYxjlfEWBsb0eKTDGXLLPo7wZUsizPAZ1Xod3nbSZCwIQiJEA5lKMqjlYs8VDpYNtmS+Bg4k54sIToEVhVH3eGGMPQHPMJcwlXd13j9bUP3UYRgOLKGhjQTVMTG/aaHoA5lKYnKgSxVseVdmD9zEpMMZcaj4LU8iz5imvWgG8PanBWiAAAY8EMJc8qhLBmiweKiPYtnqJHEzUCIMFQAsdyhh7AJpjLmEu6eoecykMP69R6JFeleFr8fwqE9fKqHF7vVJgjLlkn0d5M6SQZ3kM6rwO7zppMxcEIBAjAcylGFVzsOYY31h2gI2vxfMgwptr4JCoEyPGHoDmmEuYS7q6x1wKw89rFHqkV2Uwl/wqE9fKqHF7vVJgjLlkn0d5M6SQZ3kM6rwO7zppMxcEIBAjAcylGFVzsOYY31h2gA1zyYMImEtBVIixB/BggLmEuRSk/PnNpTAY3UWhR7qTpGdB3rSxOAN426PfbKi+MhhXZ1d0ZAqMMZeKZoPdfSnkmR298pHhXZ4ZIyAAgbQIYC6lpXew3Vo8VAZbnONAHEz8iIMWOi1i7AFojrmEuaSr++7RmvqnDsNoYBEFbSyohonpTRtND2hHxNsewyjnKwqM7fVIgTHmkn0e5c2QQp7lMajzOrzrpM1cEIBAjAQwl2JUzcGaLR4qHWzLfAkcTMwRF54ALQqj6vPGGHsAmmMuYS7p6h5zKQw/r1HokV6V4Wvx/CoT18qocXu9UmCMuWSfR3kzpJBneQzqvA7vOmkzFwQgECMBzKUYVXOw5hjfWHaAja/F8yDCm2vgkKgTI8YegOaYS5hLurrHXArDz2sUeqRXZTCX/CoT18qocXu9UmCMuWSfR3kzpJBneQzqvA7vOmkzFwQgECMBzKUYVXOw5hjfWHaADXPJgwiYS0FUiLEH8GCAuYS5FKT8+c2lMBjdRaFHupOkZ0HetLE4A3jbo99sqL4yGFdnV3RkCowxl4pmg919KeSZHb3ykeFdnhkjIACBtAhgLqWld7DdWjxUBluc40AcTPyIgxY6LWLsAWiOuYS5pKv77tGa+qcOw2hgEQVtLKiGielNG00PaEfE2x7DKOcrCozt9UiBMeaSfR7lzZBCnuUxqPM6vOukzVwQgECMBDCXYlTNwZotHiodbMt8CRxMzBEXngAtCqPq88YYewCaYy5hLunqHnMpDD+vUeiRXpXha/H8KhPXyqhxe71SYIy5ZJ9HeTOkkGd5DOq8Du86aTMXBCAQIwHMpRhVc7DmGN9YdoCNr8XzIMKba+CQqBMjxh6A5phLmEu6usdcCsPPaxR6pFdlMJf8KhPXyqhxe71SYIy5ZJ9HeTOkkGd5DOq8Du86aTMXBCAQIwHMpRhVc7DmGN9YdoANc8mDCJhLQVSIsQfwYIC5hLkUpPz5zaUwGN1FoUe6k6RnQd60sTgDeNuj32yovjIYV2dXdGQKjDGXimaD3X0p5JkdvfKR4V2eGSMgAIG0CGAupaV3sN1aPFQGW5zjQBxM/IiDFjotYuwBaI65hLmkq/vu0Zr6pw7DaGARBW0sqIaJ6U0bTQ9oR8TbHsMo5ysKjO31SIEx5pJ9HuXNkEKe5TGo8zq866TNXBCAQIwEMJdiVM3Bmi0eKh1sy3wJHEzMEReeAC0Ko+rzxhh7AJpjLmEu6eoecykMP69R6JFeleFr8fwqE9fKqHF7vVJgjLlkn0d5M6SQZ3kM6rwO7zppMxcEIBAjAcylGFVzsOYY31h2gI2vxfMgwptr4JCoEyPGHoDmmEuYS7q6x1wKw89rFHqkV2Uwl/wqE9fKqHF7vVJgjLlkn0d5M6SQZ3kM6rwO7zppMxcEIBAjAcylGFVzsOYY31h2gA1zyYMImEtBVIixB/BggLmEuRSk/PnNpTAY3UWhR7qTpGdB3rSxOAN426PfbKi+MhhXZ1d0ZAqMMZeKZoPdfSnkmR298pHhXZ4ZIyAAgbQIYC6lpXew3Vo8VAZbnONAHEz8iIMWOi1i7AFojrmEuaSr++7RmvqnDsNoYBEFbSyohonpTRtND2hHxNsewyjnKwqM7fVIgTHmkn0e5c2QQp7lMajzOrzrpM1cEIBAjAQwl2JUzcGaLR4qHWzLfAkcTMwRF54ALQqj6vPGGHsAmmMuYS7p6h5zKQw/r1HokV6V4Wvx/CoT18qocXu9UmCMuWSfR3kzpJBneQzqvA7vOmkzFwQgECMBzKUYVXOw5hjfWHaAja/F8yDCm2vgkKgTI8YegOaYS5hLurrHXArDz2sUeqRXZTCX/CoT18qocXu9UmCMuWSfR3kzpJBneQzqvA7vOmkzFwQgECMBzKUYVXOw5hjfWHaADXPJgwiYS0FUiLEH8GCAuYS5FKT8+c2lMBjdRaFHupOkZ0HetLE4A3jbo99sqL4yGFdnV3RkCowxl4pmg919KeSZHb3ykeFdnhkjIACBtAhgLqWld7DdWjxUBluc40AcTPyIgxY6LWLsAWiOuYS5pKv77tGa+qcOw2hgEQVtLKiGielNG00PaEfE2x7DKOcrCozt9UiBMeaSfR7lzZBCnuUxqPM6vOukzVwQgECMBDCXYlTNwZotHiodbMt8CRxMzBEXngAtCqPq88YYewCaYy5hLunqHnMpDD+vUeiRXpXha/H8KhPXyqhxe71SYIy5ZJ9HeTOkkGd5DOq8Du86aTMXBCAQIwHMpRhVc7DmGN9YdoCNr8XzIMKba+CQqBMjxh6A5phLmEu6usdcCsPPaxR6pFdlMJf8KhPXyqhxe71SYIy5ZJ9HeTOkkGd5DOq8Du86aTMXBCAQIwHMpRhVc7DmGN9YdoANc8mDCJhLQVSIsQfwYIC5hLkUpPz5zaUwGN1FoUe6k6RnQd60sTgDeNuj32yovjIYV2dXdGQKjDGXimaD3X0p5JkdvfKR4V2eGSMgAIG0CGAupaV3sN1aPFQGW5zjQBxM/IiDFjotYuwBaI65hLmkq/vu0Zr6pw7DaGARBW0sqIaJ6U0bTQ9oR8TbHsMo5ysKjO31SIEx5pJ9HuXNkEKe5TGo8zq866TNXBCAQIwEMJdiVM3Bmi0eKh1sy3wJHEzMEReeAC0Ko+rzxhh7AJpjLmEu6eoecykMP69R6JFeleFr8fwqE9fKqHF7vVJgjLlkn0d5M6SQZ3kM6rwO7zppMxcEIBAjAcylGFVzsOYY31h2gI2vxfMgwptr4JCoEyPGHoDmmEuYS7q6x1wKw89rFHqkV2Uwl/wqE9fKqHF7vVJgjLlkn0d5M6SQZ3kM6rwO7zppMxcEIBAjAcylGFVzsOYY31h2gA1zyYMImEtBVIixB/BggLmEuRSk/PnNpTAY3UWhR7qTpGdB3rSxOAN426PfbKi+MhhXZ1d0ZAqMMZeKZoPdfSnkmR298pHhXZ4ZIyAAgbQIYC6lpXew3Vo8VAZbnONAHEz8iIMWOi1i7AFojrmEuaSr++7RmvqnDsNoYBEFbSyohonpTRtND2hHxNsewyjnKwqM7fVIgTHmkn0e5c2QQp7lMajzOrzrpM1cEIBAjAQwl2JUzcGaLR4qHWzLfAkcTMwRF54ALQqj6vPGGHsAmmMuYS7p6h5zKQw/r1HokV6V4Wvx/CoT18qocXu9UmCMuWSfR3kzpJBneQzqvA7vOmkzFwQgECMBzKUYVXOw5hjfWHaAja/F8yDCm2vgkKgTI8YegOaYS5hLurrHXArDz2sUeqRXZTCX/CoT18qocXu9UmCMuWSfR3kzpJBneQzqvA7vOmkzFwQgECMBzKUYVXOw5hjfWHaADXPJgwiYS0FUiLEH8GCAuYS5FKT8+c2lMBjdRaFHupOkZ0HetLE4A3jbo99sqL4yGFdnV3RkCowxl4pmg919KeSZHb3ykeFdnhkjIACBtAhgLgXQ+7XXXpPp06fLrFmz5Nlnn5UhQ4bI2LFjZcyYMdLV1bXOGW6++Wa54447ZOHChfLUU0/J5ptvLnfddVefY8rc2zvA8uXL5ROf+EQr/he+8AU56aST1Lu2eKhULyqCABxM/IiEFjotYuwBaI65hLmkq/vu0Zr6pw7DaGARBW0sqIaJ6U0bTQ9oR8TbHsMo5ysKjO31SIEx5pJ9HuXNkEKe5TGo8zq866TNXBCAQIwEMJcCqHbmmWfKtddeK6NGjZLtt99e7r77brnllltk4sSJMmHChHXOMG7cOHnkkUdku+22a5k/6623Xltzqcy9vSe95JJL5LLLLpOXX34ZcymA3poQHEw09MKORQsdT4s3lnQryh+N5phLmEv5dVLkDk39U4dFCDdzD9o0w73IrN600fQAzKUiitvc4y2PbHbZbNQUGGMuNZtj2ewp5FnzlFetAN6e1GAtEICARwKYS0pVFi9eLCNGjJCjjz5aJk+e3BNt0qRJctttt7X+ZJ9Gavd65plnWtfXX399ycyjp59+uq25VObe7vl+//vfy0EHHSTHH3+8TJkyBXNJqbd2OAcTLcFw49FCx9LijSXdivJHo/lKRnCAQ361rPsOTf2Tf1r6duPRxo6tNrI3bTQ9oB0Lb3vUauZxPIztVUmBMeaSfR7lzZBCnuUxqPM6vOukzVwQgECMBDCXlKpNnTq19ZV48+bNk8GDB/dEW7BggYwePVrOPvvs1j+LvPLMpd4xit77+c9/Xl566SU5//zzZZ999sFcKiKE4T0cTAzhlgyNFiWBrXG7xRtLuhXlj0bzlYzgAIf8aln3HZr6J/+09O3Go40dW21kb9poekA7Ft72qNXM43gY26uSAmPMJfs8ypshhTzLY1DndXjXSZu5IACBGAlgLilVyz6x9Nhjj7W+Cq/3K/udox122EEOPfRQOffccwvNUtQwyoIVuXfu3LlywgknyHXXXScDBgzAXCqkgu1NHExs+ZaJjhZlaK19r8UbS7oV5Y9G85WM4ACH/GpZ9x2a+if/tPTtxqONHVttZG/aaHpAOxbe9qjVzON4GNurkgJjzCX7PMqbIYU8y2NQ53V410mbuSAAgRgJYC4pVcu+cq5fv34ye/bstSLtvvvurd9Syn7vqMiriGHUHSfv3ldeeUU+/vGPy9577y3Zb0L94Q9/CG4urVixomVa8SpOINMle/Xv37/4IO40IdApWgwdOtSET17Q7MEyth7QKZrnaZN3HQ4rCXUKhyZ6gKb+O4V7Xp3FeB1t/KrWTpsm6j+jpOkB7SiTf/b5B+POY9xED1iz/skr+7xacwaY18vcK+8m6r9e8swGAQjEQgBzSanUvvvuK4MGDZJrrrlmrUh77bWXbL311jJjxoxCs+QZRr2D5N2bfV3ftddeK7feeqtsuummmEuFFLC/yevBxH7n/mboFC2aOlRavLFknSWdormWExwwl7Q5pKl/8k9L32482tix1UbGXNISZHxGgBq3z4O6GTfxHIC5ZJ9HeTPUnWd56+n06155N1H/na41+4MABKoRwFyqxq1nlMdPLj3xxBNy8MEHt37v6fDDD2+t1eKTS1nc4cOHKwmmNZyPVPvRGy10Wlh8JY5uRfmj0XwlIzjAIb9a1n2Hpv7JPy19u/FoY8dWG9mbNpoe0I6Ftz1qNfM4Hsb2qqTAmK/Fs8+jvBlSyLM8BnVeh3edtJkLAhCIkQDmklK1vN9cGjlypJx33nmFZsn7NFLvIOu694tf/KJkBlP2dXxdXV2tYX/6059kzJgxMnbsWDnqqKNan7baaKONCq2rr5ssHiorLyaigRxM/IiFFjotYuwBaI6p0jvryYfqPUBT/3Cvzt16JNpYE64e35s2mh6AuVQ9D7QjveWRdj8ex6fAGHOp+cxLIc+ap7xqBfD2pAZrgQAEPBLAXFKqMmXKFLn00ktl3rx5Mnjw4J5oCxYskNGjR8tZZ53VMnWKvEKZS4ccckjPf5nebt7vfve78pGPfKTIsvq8x+KhsvJiIhrIwcSPWGih0yLGHoDmmEuYS7q67x6tqX/qMIwGFlHQxoJqmJjetNH0AMylMDlRJYq3PKqyB+9jUmCMudR8FqaQZ81TxlzypAFrgQAEfBPAXFLqs2jRIsk+nZR9gmny5Mk90SZNmiRz586V2267TbbYYovWd1wvWbJENttsMxk4cGCfs4Yyl+677z558cUXV5vjueeeaxldBxxwgHziE5+QHXfcsfXppaovi4fKqmuJaRwHQT9qoYVOixh7AJpjLmEu6eoecykMP69R6JFelfH3daYWZwDyzz7/YAzjEAQwl0JQ1MWglnX8yo6Gd1li3A8BCKRGAHMpgOJnnHGGzJ49W0aNGiXDhg2T+fPny5w5c2TChAkyceLE1gz333+/jB8/frW/y/7+gQceaP3JXjNnzpRly5bJMccc0/r37JNQI0aM6FlhmXvX3Ba/uRRA6AAhOJgEgBgoBFroQFq8saRbUf5oNF/JCA5wyK+Wdd+hqX/yT0vfbjza2LHVRvamjaYHtGPhbY9azTyOh7G9Kikwxlyyz6O8GVLIszwGdV6Hd520mQsCEIiRQMeZS5mBU+WV/TbRlVdeWWWovPrqqzJ9+vSWwbR06VIZMmRI66vwsk8idf/mUTtz6eKLL5Zp06b1Oe+uu+4qM2bM6LlW5t41A2IuVZI2+CAOJsGRVg6IFpXRtQZavLGkW1H+aDRfyQgOcMivlnXfoal/8k9L32482tix1Ub2po2mB7Rj4W2PWs08joexvSopMMZcss+jvBlSyLM8BnVeh3edtJkLAhCIkUDHmUt77713ZR1uv/32ymNTG2jxUJkCQw4mflRGC50WMfYANMdU6Z315EP1HqCpf7hX5249Em2sCVeP700bTQ/AXKqeB9qR3vJIux+P41NgjLnUfOalkGfNU161Anh7UoO1QAACHgl0nLnkEXInrsniobITOa25Jw4mflRGC50WMfYANMdcwlzS1X33aE39U4dhNLCIgjYWVMPE9KaNpgdgLoXJiSpRvOVRlT14H5MCY8yl5rMwhTxrnjLmkicNWAsEIOCbAOaSb33crs7iodLtZgMujINgQJjKUGihAxhjD0BzzCXMJV3dYy6F4ec1Cj3SqzL+vs7U4gxA/tnnH4xhHIIA5lIIiroY1LKOX9nR8C5LjPshAIHUCCRlLv32t7+VJ598Ul5++WUZMWJEaloH3a/FQ2XQBToNxsHEjzBoodMixh6A5phLmEu6usdcCsPPaxR6pFdlMJf8KhPXyqhxe71SYIy5ZJ9HeTOkkGd5DOq8Du86aTMXBCAQI4EkzKVf/epXcs4558jjjz/eo9HixYtb/zu7dswxx8iUKVNk3333jVHDRtYc4xvLjYBaY1IOJh5U4A32ECrE2AOoP3IfcylE9Yto6p86DKOBRRS0saAaJqY3bTQ9oB0Rb3sMo5yvKDC21yMFxphL9nmUN0MKeZbHoM7r8K6TNnNBAAIxEuh4c2nhwoUyduxY2WijjeSwww5rGUx33XWXdJtLmWiZqTRs2DD59re/HaOGjazZ4qGykY3UPCkHk5qBr2M6tNBpEWMPQHPMJcwlXd13j9bUP3UYRgOLKGhjQTVMTG/aaHoA5lKYnKgSxVseVdmD9zEpMMZcaj4LU8iz5imvWgG8PanBWiAAAY8EOt5c+tznPieZwXT99dfLlltuKdOmTZNLLrlkNXPpS1/6kjz88MPys5/9zKNGLtdk8VDpcqOBF8XBJDBQRTi0UMAT3ScXdDNXH43mmEuYS9Xrp/dIzRmAOgyjgUUUtLGgGiamN200PQBzKUxOVIniLY+q7MH7mBQYYy41n4Up5FnzlDGXPGnAWiAAAd8EOt5c2nnnneWAAw6Qr33tay0l+jKXvvWtb8kPfvAD+fWvf+1bLUers3iodLQ9s6VwEDRDWzowWpRGttqAGHsAmmMuYS7p6r57tKb+qcMwGlhEQRsLqmFietNG0wMwl8LkRJUo3vKoyh68j0mBMeZS81mYQp41TxlzyZMGrAUCEPBNoOPNpR133FE+9alPyemnn97WXPqP//gPmTNnTuv3l3gVI2DxUFls5rjv4iDoRz+00GkRYw9Ac8wlzCVd3WMuheHnNQo90qsyIt60sTgDeNuj32yovjIYV2dXdGQKjDGXimaD3X0p5JkdvfKR4V2eGSMgAIG0CHS8uTRy5EjZYIMN5Nprr+3TXHr99dflwAMPlIEDB8rVV1+dlvqK3Vo8VCqWE81QDiZ+pEILnRYx9gA0x1zCXNLVPeZSGH5eo9AjvSqDueRXmbhWRo3b65UCY8wl+zzKmyGFPMtjUOd1eNdJm7kgAIEYCXS8uXTVVVfJeeedJ8cdd5xMnDix9XtL3b+5tHz5cjn//PNbptJXv/pVOeyww2LUsJE1x/jGciOg1piUg4kHFXiDPYQKMfYA6o/cx1wKUf2631yjDsNoYBEFbSyohonpTRuLM4C3PYZRzlcUGNvrkQJjzCX7PMqbIYU8y2NQ53V410mbuSAAgRgJdLy5lH0yacKECTJv3jzZYostZKONNpLf/e538qEPfUgeffRR+ctf/iL77befXHzxxTHq19iaLR4qG9tMjRNzMKkRds5UaKHTIsYegOaYS5hLurrvHq2pf+owjAYWUdDGgmqYmN600fSAdkS87TGMcr6iwNhejxQYYy7Z51HeDCnkWR6DOq/Du07azAUBCMRIoOPNpUyUFStWyA9/+MPWJ5Qef/zx1r9nr3e/+93y6U9/WsaPHy9dXV0x6tfYmi0eKhvbTI0TczCpETbmkinsGHsA9Ye5hLkUpi1o6p86DKOBRRS0saAaJqY3bTQ9AHMpTE5UieItj6rswfuYFBhjLjWfhSnkWfOUV60A3p7UYC0QgIBHAkmYS73Bv/LKK/LCCy/IgAEDZJNNNvGoSRRrsniojGLjykVyMFECDDgcLXQwY+wBaI65hLmkq/vu0Zr6pw7DaGARBW0sqIaJ6U0bTQ/AXAqTE1WieMujKnvwPiYFxphLzWdhCnnWPGXMJU8asBYIQMA3geTMJd9yxLM6i4fKeHZffaUcBKuzCz0SLXREY+wBaI65hLmkq3vMpTD8vEahR3pVRsSbNhZnAG979JsN1VcG4+rsio5MgTHmUtFssLsvhTyzo1c+MrzLM2MEBCCQFoGOM5eWLFlSWcHBgwdXHpvaQIuHyhQYcjDxozJa6LSIsQegOeYS5pKu7jGXwvDzGoUe6VUZzCW/ysS1MmrcXq8UGGMu2edR3gwp5Fkegzqvw7tO2swFAQjESKDjzKWhQ4dW/v2kxYsXx6hhI2uO8Y3lRkCtMSkHEw8q8AZ7CBVi7AHUH7mPuRSi+kU09U8dhtHAIgraWFANE9ObNpoe0I6Itz2GUc5XFBjb65ECY8wl+zzKmyGFPMtjUOd1eNdJm7kgAIEYCXScuXTxxRevZS49+OCDMn/+fPmnf/on2WmnneQf//Ef5bnnnmu9OfLkk0/KHnvs0fr7CRMmxKhhI2u2eKhsZCM1T8rBpGbg65gOLXRaxNgD0BxzCXNJV/fdozX1Tx2G0cAiCtpYUA0T05s2mh6AuRQmJ6pE8ZZHVfbgfUwKjDGXms/CFPKsecqrVgBvT2qwFghAwCOBjjOX1oR87733yrHHHitf+9rX5JBDDllLg+uvv17OOussufTSS2W33XbzqJHLNVk8VLrcaOBFcTAJDFQRDi0U8ET3yQXdzNVHoznmEuZS9frpPVJzBqAOw2hgEQVtLKiGielNG00PwFwKkxNVonjLoyp78D4mBcaYS81nYQp51jxlzCVPGrAWCEDAN4GON5eOOOIIyX5LaerUqW2VOOmkk+SZZ56Ra665xrdajlZn8VDpaHtmS+EgaIa2dGC0KI1stQEx9gA0x1zCXNLVffdoTf1Th2E0sIiCNhZUw8T0po2mB2AuhcmJKlG85VGVPXgfkwJjzKXmszCFPGueMuaSJw1YCwQg4JtAx5tLO+64oxx55JGSGUjtXpnxNGPGDMm+Po9XMQIWD5XFZo77Lg6CfvRDC50WMfYANMdcwlzS1T3mUhh+XqPQI70qI+JNG4szgLc9+s2G6iuDcXV2RUemwBhzqWg22N2XQp7Z0SsfGd7lmTECAhBIi0DHm0v/+q//Ku985zvlhz/8YVtlP/3pT8vvfve71u8y8SpGwOKhstjMcd/FwcSPfmih0yLGHoDmmEuYS7q6x1wKw89rFHqkV2Uwl/wqE9fKqHF7vVJgjLlkn0d5M6SQZ3kM6rwO7zppMxcEIBAjgY43l7LfWvrBD34gI0aMkIkTJ7a+Iq/7tWTJErn44osl+92lMWPGyJlnnhmjho2sOcY3lhsBtcakHEw8qMAb7CFUiLEHUH/kPuZSiOrX/eYadRhGA4soaGNBNUxMb9pYnAG87TGMcr6iwNhejxQYYy7Z51HeDCnkWR6DOq/Du07azAUBCMRIoOPNpZdeekk+97nPSXYIWn/99eXtb3+7DBw4UP7617/Ks88+K6+//rp84AMfkP/6r/+SjTfeOEYNG1mzxUNlIxupeVIOJjUDX8d0aKHTIsYegOaYS5hLurrvHq2pf+owjAYWUdDGgmqYmN600fSAdkS87TGMcr6iwNhejxQYYy7Z51HeDCnkWR6DOq/Du07azAUBCMRIoOPNpUyUN954Q2bOnCk/+clP5LHHHpMXX3xRNtlkE3nf+94nBx98sHzyk5+U9dZbL0b9GluzxUNlY5upcWIOJjXCzpkKLXRaxNgD0BxzCXNJV/eYS2H4eY1Cj/SqDF+L51eZuFZGjdvrlQJjzCX7PMqbIYU8y2NQ53V410mbuSAAgRgJJGEuxSiM9zXH+MayB6YcTDyowBvsIVSIsQdQf+Q+5lKI6udr8cJQ9BeFHulPk+4VedPG4gzgbY9+s6H6ymBcnV3RkSkwxlwqmg1296WQZ3b0ykeGd3lmjIAAWSzfcwAAIABJREFUBNIikKS5tGLFCunq6kpL6cC7tXioDLxEl+E4mPiRBS10WsTYA9AccwlzSVf33aM19U8dhtHAIgraWFANE9ObNpoe0I6Itz2GUc5XFBjb65ECY8wl+zzKmyGFPMtjUOd1eNdJm7kgAIEYCSRhLmVm0o9+9CO54YYbJPs/hv/7v/+TjTbaSIYOHSojRoyQUaNGYTaVzF6Lh8qSS4jydg4mfmRDC50WMfYANMdcwlzS1T3mUhh+XqPQI70qw9fi+VUmrpVR4/Z6pcAYc8k+j/JmSCHP8hjUeR3eddJmLghAIEYCHW8uLV++XD7/+c/Lfffd1zKQttxySxk0aJD85S9/kT/96U+t32PabbfdZPr06dKvX78YNWxkzTG+sdwIqDUm5WDiQQXeYA+hQow9gPoj9zGXQlQ/X4sXhqK/KPRIf5p0r8ibNhZnAG979JsN1VcG4+rsio5MgTHmUtFssLsvhTyzo1c+MrzLM2MEBCCQFoGON5emTZsm2Z+Pf/zjctJJJ8lWW23Vo/CSJUtkypQpcvPNN8uECRPk+OOPT0t9xW4tHioVy4lmKAcTP1KhhU6LGHsAmmMuYS7p6r57tKb+qcMwGlhEQRsLqmFietNG0wPaEfG2xzDK+YoCY3s9UmCMuWSfR3kzpJBneQzqvA7vOmkzFwQgECOBjjeXPvrRj8rb3vY2+fGPf9xWn0996lPy/PPPy6233hqjho2s2eKhspGN1DwpB5Oaga9jOrTQaRFjD0BzzCXMJV3dYy6F4ec1Cj3SqzJ8LZ5fZeJaGTVur1cKjDGX7PMob4YU8iyPQZ3X4V0nbeaCAARiJNDx5tKwYcPkqKOOkn//939vq8/UqVPliiuukIcffjhGDRtZc4xvLDcCao1JOZh4UIE32EOoEGMPoP7IfcylENXP1+KFoegvCj3SnybdK/KmjcUZwNse/WZD9ZXBuDq7oiNTYIy5VDQb7O5LIc/s6JWPDO/yzBgBAQikRaDjzaUPfehDsscee8i3vvWttsqecsopMn/+fLnnnnvSUl+xW4uHSsVyohnKwcSPVGih0yLGHoDmmEuYS7q67x6tqX/qMIwGFlHQxoJqmJjetNH0gHZEvO0xjHK+osDYXo8UGGMu2edR3gwp5Fkegzqvw7tO2swFAQjESKDjzaXsd5Z+/vOfyyWXXCIf/vCH19LozjvvbP3W0v777y/ZJ5h4FSNg8VBZbOa47+Jg4kc/tNBpEWMPQHPMJcwlXd1jLoXh5zUKPdKrMnwtnl9l4loZNW6vVwqMMZfs8yhvhhTyLI9BndfhXSdt5oIABGIk0PHm0lNPPSWHHXaYvPTSS7LzzjvL8OHDZdCgQfKXv/xFsoPRr371K3nrW9/a+k2md7/73TFq2MiaY3xjuRFQa0zKwcSDCrzBHkKFGHsA9UfuYy6FqH6+Fi8MRX9R6JH+NOlekTdtLM4A3vboNxuqrwzG1dkVHZkCY8ylotlgd18KeWZHr3xkeJdnxggIQCAtAh1vLmVyPvroo3LOOefIgw8+uJa6H/jAB+Tss8+W973vfWkpr9ytxUOlcklRDOdg4kcmtNBpEWMPQHPMJcwlXd13j9bUP3UYRgOLKGhjQTVMTG/aaHpAOyLe9hhGOV9RYGyvRwqMMZfs8yhvhhTyLI9BndfhXSdt5oIABGIkkIS51C3MH//4x5bR9OKLL8omm2wi2267rQwZMiRG3Rpfs8VDZeObqmEBHExqgFxwCrQoCKrNbTH2ADTHXMJc0tU95lIYfl6j0CO9KsPX4vlVJq6VUeP2eqXAGHPJPo/yZkghz/IY1Hkd3nXSZi4IQCBGAkmZS1YCvfbaazJ9+nSZNWuWPPvssy3DauzYsTJmzBjp6upa57Q333yz3HHHHbJw4ULJvsJv8803l7vuuqvPMUXvff755+W6666TefPmyeOPPy4vv/yybL311nLQQQfJkUceKW95y1vUKGJ8Y1m96QABOJgEgBgoBFroQMbYA9AccwlzSVf3mEth+HmNQo/0qgzmkl9l4loZNW6vVwqMMZfs8yhvhhTyLI9BndfhXSdt5oIABGIkgLkUQLUzzzxTrr32Whk1apRsv/32cvfdd8stt9wiEydOlAkTJqxzhnHjxskjjzwi2223XctcWm+99dqaS0XvzUyl448/XvbYYw/ZbbfdWp/SeuCBB+SnP/1p6zenZsyYIeuvv75q5zG+sazacKDBHEwCgQwQBi10EGPsAWiOuYS5pKt7zKUw/LxGoUd6VQZzya8yca2MGrfXKwXGmEv2eZQ3Qwp5lsegzuvwrpM2c0EAAjESSMJcuvfee+XKK6+Uxx57TJYuXSqvv/76WlplnzBatGhRaQ0XL14sI0aMkKOPPlomT57cM37SpEly2223tf5kn0Zq93rmmWda1zOzJzOPnn766bbmUtF7f//737emyz6t1Pt14YUXyne+8x2ZNm2a7LfffqX32ntAjG8sqzYcaDAHk0AgA4RBCx3EGHsAmq/UHA5w0FW/iKb+yT8tfbvxaGPHVhvZmzaaHtCOhbc9ajXzOB7G9qqkwBhzyT6P8mZIIc/yGNR5Hd510mYuCEAgRgIdby5dc8018pWvfEVWrFgh73rXu2TQoEGtTwf19co+0VP2NXXq1NZX4mWfFho8eHDP8AULFsjo0aPl7LPPbv2zyCvPXOodo8y93eOy35s6+OCD5cQTT5TjjjuuyJLa3mPxUKlaUCSDOZj4EQotdFrE2APQfKXmcICDrvoxl7T8vI6nN3hVxl/ftjgDkH/2+QdjGIcggLkUgqIuBrWs41d2NLzLEuN+CEAgNQIdby7tvffe8ve//12+973vydChQ4Prm31iKftEVPZVeL1fy5cvlx122EEOPfRQOffccwvNW8YwKnNv9+S/+MUv5LOf/WzLbDviiCMKrandTRYPlaoFRTKYg4kfodBCp0WMPQDNMVV6Zz35UL0HaOof7tW5W49EG2vC1eN700bTA9pR8LbH6mr5HQlje21SYIy5ZJ9HeTOkkGd5DOq8Du86aTMXBCAQI4GON5cygyf7LaQvf/nLJvocdNBB0q9fP5k9e/Za8XfffffWbylddtllheYuYxiVuTeb/I033pAjjzxSHn74YZk7d27rE1yaV3aozD4NNmDAAE2Y5Ma+8sorrT33798/ub1723CnaGFhmhfRKsYe0CmaF9FnXffAYSWdTuHQRA/Q1H+ncNfWocfxaONRlXX3qybqP1uRpge0o0z+2ecfjDuPcRM9YM36J6/s82rNGWBeL3OvvJuo/3rJMxsEIBALgY43l7JPDr3vfe+T888/30STfffdt2XUZF+/t+Zrr732av3uUdGv2ytjGJW5N1tX99f3nXnmma3fdtK+LB4qtWuKYbzXg0kM7EKvsVO0aOpQGWMP6BTNtbUAB8wlbQ5p6p/809K3G482dmy1kdtp00lnAPJPmyX542Gcz0h7R92Mm+gBmEvaLNGPrzvP9CuOO4JX3k3Uf9xKsnoIQMCKQMebSz/72c/ktNNOkx/96Efyz//8z8E5xvDJpe9///vy1a9+tfVVeNlX4oV4WXwdRoh1eY/BR6r9KIQWOi1i7AFovlJzOMBBV/385pKWn9fx9Aavyvjr2xZnAPLPPv9gDOMQBPhavBAUdTGoZR2/sqPhXZYY90MAAqkR6HhzKRP05ptvlvPOO0+y31/adttt236V24gRI0rrn/ebSyNHjmzNXeRV5tNIRe/Nvq7vjDPOkAMPPFAuuOACWW+99YosJfcei4fK3Ek74AYOJn5ERAudFjH2ADTHVOmd9eRD9R6gqX+4V+duPRJtrAlXj+9NG00PaEfB2x6rq+V3JIzttUmBMeaSfR7lzZBCnuUxqPM6vOukzVwQgECMBDreXHr55Zdbn1z6+c9/3vqNoOzV1dW1mlbZ32d/t3jx4tIaTpkyRS699FKZN2+eDB48uGf8ggULZPTo0XLWWWfJmDFjCsUtahhlwYrce9NNN8nJJ58sH/7wh2XatGmywQYbFFpHkZssHiqLzBv7PRxM/CiIFjotYuwBaL5SczjAQVf9fHJJy8/reHqDV2X89W2LMwD5Z59/MIZxCAKYSyEo6mJQyzp+ZUfDuywx7ocABFIj0PHmUmYsXX/99bLddtvJfvvt1/p9pPXXX79PnbNPGZV9LVq0SLJx2SeYJk+e3DN80qRJMnfuXLnttttkiy22aP1w+JIlS2SzzTaTgQMH9jlNEcOoe2DevdncJ554ouyyyy4t86tfv35lt7bO+y0eKoMu0GkwDiZ+hEELnRYx9gA0x1TpnfXkQ/UeoKl/uFfnbj0SbawJV4/vTRtND2hHwdseq6vldySM7bVJgTHmkn0e5c2QQp7lMajzOrzrpM1cEIBAjAQ63lz64Ac/KO95z3vkhz/8YbCvhFtT6Oxr57Kvnxs1apQMGzZM5s+fL3PmzJEJEybIxIkTW7fff//9Mn78+NX+Lvv7Bx54oPUne82cOVOWLVsmxxxzTOvfs09C9f6qvqL3Lly4sPVpqQ033FBOPfVU6d+//2pLfuc73yk77bSTKl8tHipVC4pkMAcTP0KhhU6LGHsAmq/UHA5w0FU/n1zS8vM6nt7gVRl/fdviDED+2ecfjGEcggDmUgiKuhjUso5f2dHwLkuM+yEAgdQIdLy5tOuuu8onP/nJ1T5VFFrkV199VaZPn94ymJYuXSpDhgxpmTvZp4u6v4Kvnbl08cUXt76yrq9XtvYZM2b0XCp6b7aO008/ve02s09anX/++SoMFg+VqgVFMpiDiR+h0EKnRYw9AM0xVXpnPflQvQdo6h/u1blbj0Qba8LV43vTRtMD2lHwtsfqavkdCWN7bVJgjLlkn0d5M6SQZ3kM6rwO7zppMxcEIBAjgY43l0444YTWp4GuvPLKGPVxu2aLh0q3mw24MA4mAWEqQ6GFDmCMPQDNMZcwl3R13z1aU//UYRgNLKKgjQXVMDG9aaPpAZhLYXKiShRveVRlD97HpMAYc6n5LEwhz5qnvGoF8PakBmuBAAQ8Euh4c+mZZ55pfYro4IMPluOOOy74bw95FLWONVk8VNax7qbn4GDStAIcEkMpEGMPoP4wlzCXwnQATf1Th2E0sIiCNhZUw8T0po2mB2AuhcmJKlG85VGVPXgfkwJjzKXmszCFPGueMu8beNKAtUAAAr4JdLy5lP3O0QsvvCCPPvqoDBgwQN71rne1/rnmK/v6Oj7dVDxZLR4qi88e750cBP1ohxY6LWLsAWiOuYS5pKv77tGa+qcOw2hgEQVtLKiGielNG00PwFwKkxNVonjLoyp78D4mBcaYS81nYQp51jxlzCVPGrAWCEDAN4GON5eGDh1aSIHMXFq8eHGhe7lJ92PeKfPjIOhHfbTQaWHxxpJuRfmj0RxzCXMpv06K3KGpf+qwCOFm7kGbZrgXmdWbNpoegLlURHGbe7zlkc0um42aAmPMpWZzLJs9hTxrnjLmkicNWAsEIOCbQMebS77xx7s6i4fKeGkUXzkHweKsrO9ECx3hGHsAmmMuYS7p6r57tKb+qcMwGlhEQRsLqmFietNG0wMwl8LkRJUo3vKoyh68j0mBMeZS81mYQp41TxlzyZMGrAUCEPBNAHOpjT7Z/2Fnf0aMGOFbwYZWZ/FQ2dBWap2Wg2CtuNc5GVrotIixB6A55hLmkq7uMZfC8PMahR7pVRl//5W6xRmA/LPPPxjDOAQBzKUQFHUxqGUdv7Kj4V2WGPdDAAKpEcBcaqP4tGnT5JJLLuGr8trwsXioTKH4OJj4URktdFrE2APQHHMJc0lX95hLYfh5jUKP9KoM5pJfZeJaGTVur1cKjDGX7PMob4YU8iyPQZ3X4V0nbeaCAARiJIC5hLlUKW9jfGO50kYDD+JgEhioIhxaKOBJnL+7huaYS5hLurrHXArDz2sUeqRXZTCX/CoT18qocXu9UmCMuWSfR3kzpJBneQzqvA7vOmkzFwQgECMBzCXMpUp5i7lUCRs/vlkNm8koDok6rDH2ADTHXMJc0tU95lIYfl6j0CO9KoO55FeZuFZGjdvrlQJjzCX7PMqbIYU8y2NQ53V410mbuSAAgRgJYC5hLlXK2xjfWK600cCDOJgEBqoIhxYKeHxySQev4dHkPiabNgU1ZwDyT0vfbjza2LHVRvamjaYHtGPhbY9azTyOh7G9Kikwxlyyz6O8GVLIszwGdV6Hd520mQsCEIiRAOYS5lKlvLV4qKy0kMgGcTDxIxha6LSIsQegOaZK76wnH6r3AE39w706d+uRaGNNuHp8b9poegDmUvU80I70lkfa/XgcnwJjzKXmMy+FPGue8qoVwNuTGqwFAhDwSABzCXOpUl5aPFRWWkhkgziY+BEMLXRaxNgD0BxzCXNJV/fdozX1Tx2G0cAiCtpYUA0T05s2mh6AuRQmJ6pE8ZZHVfbgfUwKjDGXms/CFPKsecqYS540YC0QgIBvAphLmEuVMtTiobLSQiIbxEHQj2BoodMixh6A5phLmEu6usdcCsPPaxR6pFdl+M0lv8rEtTJq3F6vFBhjLtnnUd4MKeRZHoM6r8O7TtrMBQEIxEgAcwlzqVLexvjGcqWNBh7EwSQwUEU4tFDA4zeXdPAaHk3uY7JpU1BzBiD/tPTtxqONHVttZG/aaHpAOxbe9qjVzON4GNurkgJjzCX7PMqbIYU8y2NQ53V410mbuSAAgRgJYC5hLlXKW4uHykoLiWwQBxM/gqGFTosYewCaY6r0znryoXoP0NQ/3Ktztx6JNtaEq8f3po2mB2AuVc8D7UhveaTdj8fxKTDGXGo+81LIs+Ypr1oBvD2pwVogAAGPBDCX2qiS/R/I4sWLZeTIkR51a3xNFg+VjW+qhgV4Opi88tL/yfJXlku//v2k/4CNati9ryk8aeGLTLHVxNgD0HyVubTBBhvIkC23ogeIyNChQ4slPXf1ENDUP3XoI5H6OgOgjQ9t+lqFN200PQBzqbk8651HqT8HWKngrVYt9om5ZEG1XMwU8qwcEdu74W3Ll+gQgED8BDCX4tewkR1YPFQ2spGaJ/VwMPnbX1+UPz21VH70zRvkmcf/JFtu8w751KmHyDvevbm8deAmNRNpbjoPWjS3e/3MMfYANF+p+x+e/KP87dmXZdbUn9IDMJcqNQNN/VOHlZAHG7SuM8Afl/6hNQ+GazDcwQJ5qxtND8BcCpYWpQNlebTJW94qr/y/vyf/HFAaXsEB3mq14LJL3Ya5VAqXyc0p5JkJuIpB4V0RHMMgAIFkCHScuZQ9EHd1dZUWMBuzaNGi0uNSHWDxUJkCy6YPJtmbStdO/Ylcfd7stXB/+oxD5fB//0QyBlPTWsSe7zH2ADQXyXrAjy+4Qa45/3p6wG9+wxvpFRuRpv6pw4rQAwzLOwPsf8y/yYt//xvmUgDWoUN4qxtND8BcCp0dxeNl/3HJrf91B2eA4shK3+mtVktvoMAAzKUCkIxvSSHPjBGWCg/vUri4GQIQSJBAx5lLp512WiVzKdP+61//eoIpUG3LFg+V1VYS16imDya//Z8n5LidJ7eF9p1ffUP+efg/xQW14mqb1qList0Mi7EHoLkIPWBVCZEP1duJpv7hXp27dmRe/V90/3my4WZd8t73vlc7FeMDE/BWN5oe0A6Ntz0GltBFuMW//K2csNsZPAcYqpFCHmMuGSZQwdAp5FlBFLXcBu9aMDMJBCAQMYGOM5ci1iKqpVs8VEYFoOJimzyYZN+tPuWY/5Q7f3xP29Xv9akPyZcu+6JslMBvMDWpRcX0cTUsxh6Quub0gNVLKPV80DQUTf3DXUO++tii9X/k+YfLVu/aqvpEjDQh4K1uND0Ac8kkRXKDFu0BqTwH5AKreIO3Wq24jXUOw1yyoFouZgp5Vo6I7d3wtuVLdAhAIH4CmEvxa9jIDiweKhvZSM2TNnkwWfaXF+SMj50rjy14ou2u37fzNnLezWfI2wZtWjOZ+qdrUov6dxt+xhh7QOqa0wMwl0J1Ak39p16HoTQoG6do/U+++nh55zZblw3P/cYEvNWNpgdgLhknS5vwRXtAKs8BVip4q1WLfWIuWVAtFzOFPCtHxPZueNvyJToEIBA/Acyl+DVsZAcWD5WNbKTmSZs8mPBfLPLGcsh0j7EHNFl/IdlXjUUPoAdUzZ01x2nqP/U6DKVB2ThF659PLpUlW8/93upG0wMwl+rJmTVnKdoD+OSSTh9vtarbTd+jMZcsqJaLmUKelSNieze8bfkSHQIQiJ9Ax5lL++yzTyVVurq6ZO7cuZXGpjjI4qEyBY5NH0zyfm+B31xKIQvD7DHGHtB0/YUhr4tCD1jFj3yonkua+od7de7akXn1z28uaQnbjfdWN5oegLlklyd5kfnNpTxC+uvealW/o7UjYC5ZUC0XM4U8K0fE9m542/IlOgQgED+BjjOXxo0bV1mVGTNmVB6b2kCLh8oUGDZ9MPnbX1+Ua6f+RK4+b/ZauEd/+VA57KRPyFsHbpKCFNK0FrFDjrEHoLlI1gN+fMENcs3519MDfvObFoOhQ4fGXo61r19T/9Rh7XL1TJh3Btjv6H+TF//+N2qiOYnazuytbjQ9AHOpuQT7w5N/lFv/6w7OAIYSeKtVi61iLllQLRczhTwrR8T2bnjb8iU6BCAQP4GOM5filySOHVg8VMaxc90qPRxMsjeX/vTUUvnxt26QJY//WQZvs4WMOuUQece7N0/GWMpU9KCFLpuaHR1jD0DzlTmTvbn0t7+8LLOn/pQegLlUqZFo6p86rIQ82KB1nQH+uPQPrXkwXIPhDhbIW91oegDmUrC0KB0oy6NN3vJWeeX//T3554DS8AoO8FarBZdd6jbMpVK4TG5OIc9MwFUMCu+K4BgGAQgkQwBzKRmpw27U4qEy7Ap9RvN0MPm/l/5P/v7KcnlL/36y0YCNfAIzXJUnLQy3aRY6xh6A5ivTIeOwwQYbyFZbbkUP4I30Sj1CU//UYSXkwQf1dQZAm+CYgwX0po2mB2AuBUuL0oF651HqzwGl4RUc4K1WCy671G2YS6VwmdycQp6ZgKsYFN4VwTEMAhBIhgDmUjJSh92oxUNl2BX6jMbBxI8uaKHTIsYegOarzKXsf6X+6QTyoXoP0NQ/3Ktztx6JNtaEq8f3po2mB2AuVc8D7UhveaTdj8fxKTDGXGo+81LIs+Ypr1oBvD2pwVogAAGPBDreXBo/fnwh7l1dXXLllVcWupebRCweKlPgysHEj8poodMixh6A5phLvbOefKjeAzT1D/fq3K1Hoo014erxvWmj6QGYS9XzQDvSWx5p9+NxfAqMMZeaz7wU8qx5yphLnjRgLRCAgG8CHW8u7b333n0q8NJLL8myZcskM5UGDRokG264odx+++2+1XK0OouHSkfbM1sKB0EztKUDo0VpZKsNiLEHoDnmEuaSru67R2vqnzoMo4FFFLSxoBompjdtND0AcylMTlSJ4i2PquzB+5gUGGMuNZ+FKeRZ85QxlzxpwFogAAHfBDreXFoX/j/+8Y/yjW98Q/785z/L5ZdfLgMGDPCtlqPVWTxUOtqe2VI4CJqhLR0YLUojw1zSIXMzmtzHZNMmo+YMQP5p6duNRxs7ttrI3rTR9ADMJW02VB/vLY+q78TvyBQYYy41n38p5FnzlDGXPGnAWiAAAd8EkjaXMmlee+01GTlypOy8885y9tln+1bL0eosHiodbc9sKRwEzdCWDowWpZFhLumQuRlN7mMuaZNRcwYg/7T07cajjR1bbWRv2mh6AOaSNhuqj/eWR9V34ndkCowxl5rPvxTyrHnKmEueNGAtEICAbwLJm0uZPOeee67cdNNNcs899/hWy9HqLB4qHW3PbCkcBM3Qlg6MFqWRYS7pkLkZTe5jLmmTUXMGIP+09O3Go40dW21kb9poegDmkjYbqo/3lkfVd+J3ZAqMMZeaz78U8qx5yphLnjRgLRCAgG8CmEsictppp8mcOXPkoYce8q2Wo9VZPFQ62p7ZUjgImqEtHRgtSiPDXNIhczOa3Mdc0iaj5gxA/mnp241HGzu22sjetNH0AMwlbTZUH+8tj6rvxO/IFBhjLjWffynkWfOUMZc8acBaIAAB3wSSNpdWrFghP/3pT+X000+XHXbYQX7wgx/4VsvR6iweKh1tz2wpHATN0JYOjBalkWEu6ZC5GU3uYy5pk1FzBiD/tPTtxqONHVttZG/aaHoA5pI2G6qP95ZH1Xfid2QKjDGXms+/FPKsecqYS540YC0QgIBvAh1vLu2zzz59KvD666/Lc8891/rNpQEDBsgVV1whw4YN862Wo9VZPFQ62p7ZUjgImqEtHRgtSiPDXNIhczOa3Mdc0iaj5gxA/mnp241HGzu22sjetNH0AMwlbTZUH+8tj6rvxO/IFBhjLjWffynkWfOUMZc8acBaIAAB3wQ63lwaN25cnwqst956summm8p2220nhx56qGy++ea+lXK2OouHSmdbNFkOB0ETrJWCokUlbD2DYuwBaI6p0jvryYfqPUBT/3Cvzt16JNpYE64e35s2mh6AuVQ9D7QjveWRdj8ex6fAGHOp+cxLIc+ap4y55EkD1gIBCPgm0PHmkm/88a7O4qEyXhrFV85BsDgr6zvRQkc4xh6A5phLmEu6uu8eral/6jCMBhZR0MaCapiY3rTR9ADMpTA5USWKtzyqsgfvY1JgjLnUfBamkGfNU8Zc8qQBa4EABHwTwFwKoE/21XrTp0+XWbNmybPPPitDhgyRsWPHypgxY6Srq2udM9x8881yxx13yMKFC+Wpp55qfYLqrrvu6nNMmXuzAPfee69ceOGFsnjxYunfv7985CMfkVNOOUUGDhyo3rXFQ6V6UREE4CDoRyS00GkRYw9Ac8wlzCVd3WMuheHnNQo90qsyIt60sTgDeNuj32yovjIYV2dXdGQKjDGXimZithx8AAAgAElEQVSD3X0p5JkdvfKR4V2eGSMgAIG0CGAuBdD7zDPPlGuvvVZGjRol22+/vdx9991yyy23yMSJE2XChAnrnCH72r5HHnmk9fV8mbmUfV1fO3OpzL2//OUv5aijjpJtt91WDjvsMPnrX/8ql19+uQwePFhmzpwpG220kWrnFg+VqgVFMpiDiR+h0EKnRYw9AM0xlzCXdHWPuRSGn9co9EivymAu+VUmrpVR4/Z6pcAYc8k+j/JmSCHP8hjUeR3eddJmLghAIEYCSZhL8+bNkxkzZsiiRYvkb3/7m7zxxhtraZV9wii7XvaVfSpoxIgRcvTRR8vkyZN7hk+aNEluu+221p91/Z7TM88807q+/vrrS2YePf30023NpTL3ZmtatmyZ3HTTTbLxxhu31nXnnXfKscceK6effrp85jOfKbvV1e6P8Y1l1YYDDeZgEghkgDBooYMYYw9Ac8wlzCVd3WMuheHnNQo90qsymEt+lYlrZdS4vV4pMMZcss+jvBlSyLM8BnVeh3edtJkLAhCIkUDHm0vZV9VlnyzKzJsPfOADLSNngw026FOrr3/966U1nDp1ausr8TIDK/tUUPdrwYIFMnr0aDn77LNb/yzyyjOXesdY171PPvmkHHDAAX1+cmr//feXt73tba1PWmleMb6xrNlvqLEcTEKR1MdBCx3DGHsAmq/UHA5w0FW/iKb+yT8tfbvxaGPHVhvZmzaaHtCOhbc9ajXzOB7G9qqkwBhzyT6P8mZIIc/yGNR5Hd510mYuCEAgRgIdby599KMflZdfflmuvvpq2WqrrYJrlH1i6bHHHmt9FV7v1/Lly2WHHXaQQw89VM4999xC84Yyl37yk5/IySefLJdddpnsueeeq82d/f2tt94qv/71r1uGW9VXdqhcsWKFDBgwoGqIJMe98sorrX1nv4HFq1kCnaLF0KFDGwEZYw/oFM21gsNhJcFO4dBED9DUf6dw19ahx/Fo41GVdferJuo/W5GmB7SjTP7Z5x+MO49xEz1gzfonr+zzas0ZYF4vc6+8m6j/eskzGwQgEAuBjjeXhg0bJkcccYR8+ctfNtHkoIMOkn79+sns2bPXir/77ru3fkspM3mKvEKZS9/73vfkm9/8ptx4442t31zq/cr+Prs+f/58GTRoUJFl9XmPxUNl5cVENNDrwSQihMGW2ilaNHWojLEHdIrm2iKAA+aSNoc09U/+aenbjUcbO7bayO206aQzAPmnzZL88TDOZ6S9o27GTfQAzCVtlujH151n+hXHHcEr7ybqP24lWT0EIGBFoOPNpQMPPFB23HFHOe+880wY7rvvvi2T5pprrlkr/l577SVbb7116/eeirxCmUuXXHKJXHTRRXLLLbfIe97zntWmvvDCC+U73/lO67egNJ/ksvg6jCKMYr+Hj1T7URAtdFrE2APQfKXmcICDrvr5WjwtP6/j6Q1elfHXty3OAOSfff7BGMYhCPC1eCEo6mJQyzp+ZUfDuywx7ocABFIj0PHm0o9//GOZMmWKXH/99bLlllsG1zflTy5lMIcPHx6caScH5GDiR1200Glh8caSbkX5o9EcU6V3lpAP+TXT7g5N/cO9OnfrkWhjTbh6fG/aaHpAOwre9lhdLb8jYWyvTQqMMZfs8yhvhhTyLI9BndfhXSdt5oIABGIk0PHm0gMPPCBXXHFF6zeGxo8f3/qauE022aRPrXbZZZfSGub95tLIkSMLf2oq1CeX8n5zKftE00MPPaT+zSXMpdLpwicGyiMzG8EhUYfW4o0l3YryR6M55hLmUn6dFLlDU//UYRHCzdyDNs1wLzKrN200PQBzqYjiNvd4yyObXTYbNQXGmEvN5lg2ewp51jzlVSuAtyc1WAsEIOCRQMebS9n3kHZ1dcmKFSta/1zXa/HixaU1yj4Vdemll8q8efNk8ODBPeMXLFggo0ePlrPOOkvGjBlTKG4oc+mJJ56Qj33sYzJx4kSZMGHCanPvv//+summm8rMmTMLrandTRYPlaoFRTKYg4kfodBCp0WMPQDNV2oOBzjoqp+vxdPy8zqe3uBVGX992+IMQP7Z5x+MYRyCAOZSCIq6GNSyjl/Z0fAuS4z7IQCB1Ah0vLl08cUX55pK3aKvacQUSYZFixZJ9umk7BNMkydP7hkyadIkmTt3buu3jbbYYgvJfgRwyZIlstlmm8nAgQP7DB3KXMqCH3LIIfLCCy/ITTfdJBtvvHFrvjvvvFOOPfbY1jqz9WpeFg+VmvXEMpaDiR+l0EKnRYw9AM0xVXpnPflQvQdo6h/u1blbj0Qba8LV43vTRtMD2lHwtsfqavkdCWN7bVJgjLlkn0d5M6SQZ3kM6rwO7zppMxcEIBAjgY43l+oQ5YwzzpDZs2fLqFGjZNiwYTJ//nyZM2dO61ND2aeHstf999/f+lq+3n+X/X32tX3Zn+yVfZpo2bJlcswxx7T+Pfsk1IgRI3q2UObe++67r2UgZZ/cOvzww+W5555rfT1gZnTNmjVL+vfvr0Jj8VCpWlAkgzmY+BEKLXRaxNgD0Hyl5nCAg676+eSSlp/X8fQGr8r469sWZwDyzz7/YAzjEAQwl0JQ1MWglnX8yo6Gd1li3A8BCKRGAHOpjeJXXnmlXHXVVa1PHuW9Xn31VZk+fXrLYFq6dKkMGTKk9VV42SeRur+Kr525lH2yatq0aX1Oseuuu8qMGTN6rpW5Nxt0zz33yIUXXijZ1/1lZtJee+0lp5xyigwaNChvS7nXLR4qcyftgBs4mPgRES10WsTYA9AcU6V31pMP1XuApv7hXp279Ui0sSZcPb43bTQ9oB0Fb3usrpbfkTC21yYFxphL9nmUN0MKeZbHoM7r8K6TNnNBAAIxEsBcaqNaZvhccsklLWOG19oELB4qU+DMwcSPymih0yLGHoDmmEuYS7q67x6tqX/qMIwGFlHQxoJqmJjetNH0AMylMDlRJYq3PKqyB+9jUmCMudR8FqaQZ81TXrUCeHtSg7VAAAIeCWAuYS5VykuLh8pKC4lsEAcTP4KhhU6LGHsAmmMuYS7p6h5zKQw/r1HokV6V4Wvx/CoT18qocXu9UmCMuWSfR3kzpJBneQzqvA7vOmkzFwQgECMBzCXMpUp5G+Mby5U2GngQB5PAQBXh0EIBT3S/uaKbufpoNMdcwlyqXj+9R2rOANRhGA0soqCNBdUwMb1po+kB7Yh422MY5XxFgbG9Hikwxlyyz6O8GVLIszwGdV6Hd520mQsCEIiRAOYS5lKlvLV4qKy0kMgGcTDxIxha6LSIsQegOeYS5pKu7rtHa+qfOgyjgUUUtLGgGiamN200PQBzKUxOVIniLY+q7MH7mBQYYy41n4Up5FnzlFetAN6e1GAtEICARwKYS5hLlfLS4qGy0kIiG8TBxI9gaKHTIsYegOaYS5hLurrHXArDz2sUeqRXZfhaPL/KxLUyatxerxQYYy7Z51HeDCnkWR6DOq/Du07azAUBCMRIAHMJc6lS3sb4xnKljQYexMEkMFBFOLRQwONr8XTwGh5N7mOyaVNQcwYg/7T07cajjR1bbWRv2mh6QDsW3vao1czjeBjbq5ICY8wl+zzKmyGFPMtjUOd1eNdJm7kgAIEYCWAuYS5VyluLh8pKC4lsEAcTP4KhhU6LGHsAmmOq9M568qF6D9DUP9yrc7ceiTbWhKvH96aNpgdgLlXPA+1Ib3mk3Y/H8SkwxlxqPvNSyLPmKa9aAbw9qcFaIAABjwQwlzCXKuWlxUNlpYVENoiDiR/B0EKnRYw9AM0xlzCXdHXfPVpT/9RhGA0soqCNBdUwMb1po+kBmEthcqJKFG95VGUP3sekwBhzqfksTCHPmqeMueRJA9YCAQj4JoC51Eaf//7v/5arrrpKbr/9dt8KNrQ6i4fKhrZS67QcBGvFvc7J0EKnRYw9AM0xlzCXdHWPuRSGn9co9EivyvCbS36ViWtl1Li9Xikwxlyyz6O8GVLIszwGdV6Hd520mQsCEIiRQMebS+PHj5dDDz1URowY0VafG264QWbNmtUyk3gVIxDjG8vFdmZ7FwcTW75loqNFGVpr3xtjD0BzzCXMJV3dYy6F4ec1Cj3SqzKYS36ViWtl1Li9Xikwxlyyz6O8GVLIszwGdV6Hd520mQsCEIiRQMebS0OHDpUJEya0/rR7/ed//qdcdNFFsnjx4hg1bGTNMb6x3AioNSblYOJBBd5gD6FCjD2A+iP3MZdCVL+Ipv6pwzAaWERBGwuqYWJ600bTA9oR8bbHMMr5igJjez1SYIy5ZJ9HeTOkkGd5DOq8Du86aTMXBCAQIwHMJRE5//zz5eqrr5aHHnooRg0bWbPFQ2UjG6l5Ug4mNQNfx3RoodMixh6A5phLmEu6uu8eral/6jCMBhZR0MaCapiY3rTR9ADMpTA5USWKtzyqsgfvY1JgjLnUfBamkGfNU161Anh7UoO1QAACHgl0pLl0/fXX97A+7bTTZN999239WfP1xhtvyJ///GfJfl9p8ODBct1113nUyOWaLB4qXW408KI4mAQGqgiHFgp4ovvkgm7m6qPRfCU7OMChehWtHKk5A5B/Wvp249HGjq02sjdtND2gHQtve9Rq5nE8jO1VSYEx5pJ9HuXNkEKe5TGo8zq866TNXBCAQIwEOtJcyr4Kr6urq5AeK1askH79+snUqVP7NKAKBUnwJouHyhQwcjDxozJa6LSIsQegOaZK76wnH6r3AE39w706d+uRaGNNuHp8b9poegDmUvU80I70lkfa/XgcnwJjzKXmMy+FPGue8qoVwNuTGqwFAhDwSKAjzaXuTyBlxtEZZ5zRMo322WeftfhnBtQ//MM/yA477CCbbbaZR33crsniodLtZgMujINJQJjKUGihAxhjD0BzzCXMJV3dd4/W1D91GEYDiyhoY0E1TExv2mh6AOZSmJyoEsVbHlXZg/cxKTDGXGo+C1PIs+YpYy550oC1QAACvgl0pLnUG/npp5/e1lzyLY3v1Vk8VPrecZjVcRAMwzFEFLTQUYyxB6A55hLmkq7uMZfC8PMahR7pVRl/X2dqcQYg/+zzD8YwDkEAcykERV0MalnHr+xoeJclxv0QgEBqBDreXEpN0Lr2a/FQWdfam5yHg0mT9FefGy10WsTYA9AccwlzSVf3mEth+HmNQo/0qgzmkl9l4loZNW6vVwqMMZfs8yhvhhTyLI9BndfhXSdt5oIABGIk0PHm0hNPPNH68fJ/+7d/k0022aSl0fLly+Wiiy6SefPmtX5v6TOf+YwccsghMerX2JpjfGO5MVi9JuZg4kEF3mAPoUKMPYD6I/cxl0JUv4im/qnDMBpYREEbC6phYnrTRtMD2hHxtscwyvmKAmN7PVJgjLlkn0d5M6SQZ3kM6rwO7zppMxcEIBAjgY43l0466SR54IEH5K677pL11luvpdE3vvENueKKK2TjjTeWV199VV577TW5/PLLZffdd49Rw0bWbPFQ2chGap6Ug0nNwNcxHVrotIixB6A55hLmkq7uu0dr6p86DKOBRRS0saAaJqY3bTQ9AHMpTE5UieItj6rswfuYFBhjLjWfhSnkWfOUV60A3p7UYC0QgIBHAh1vLu29994yfPhwueCCC1r8MzMpM5He+973ypVXXikvvPCCjBw5Ut7//vfLpZde6lEjl2uyeKh0udHAi+JgEhioIhxaKOCJ7pMLupmrj0ZzzCXMper103uk5gxAHYbRwCIK2lhQDRPTmzaaHoC5FCYnqkTxlkdV9uB9TAqMMZeaz8IU8qx5yphLnjRgLRCAgG8CHW8u7bjjjjJu3Dj50pe+1FLiV7/6lYwdO1a++c1vysEHH9z6u3POOUduv/321qebeBUjYPFQWWzmuO/iIOhHP7TQaRFjD0DzlZrDAQ666teZy+Sflr7deLSxY6uN7E0bizOAtz1qNfM4Hsb2qqTAGHPJPo/yZkghz/IY1Hkd3nXSZi4IQCBGAh1vLn3wgx9smUhf/vKXW/pMmzZNLrnkkpaR9Pa3v731d1OmTGl9imnhwoUxatjImi0eKhvZSM2TcjCpGfg6pkMLnRYx9gA0x1TpnfXkQ/UeoKl/uFfnbj0SbawJV4/vTRtND2hHwdseq6vldySM7bVJgTHmkn0e5c2QQp7lMajzOrzrpM1cEIBAjAQ63lz61Kc+JcuWLZMbb7xRurq6WkbThhtu2Pr37lf2qabskDRv3rwYNWxkzRYPlY1spOZJOZjUDBxzyQx4jD2A+sNcwlwK0xI09U8dhtHAIgraWFANE9ObNpoegLkUJieqRPGWR1X24H1MCowxl5rPwhTyrHnKq1YAb09qsBYIQMAjgY43lzIT6dRTT5V3vOMdLVPpD3/4g5x11lny6U9/ukeP/fbbT7bZZhv57ne/61Ejl2uyeKh0udHAi+JgEhioIhxaKODxm0s6eA2PJvcx2bQpqDkDkH9a+nbj0caOrTayN200PQBzSZsN1cd7y6PqO/E7MgXGmEvN518KedY8ZcwlTxqwFghAwDeBjjeXMvxXXXWVXH/99S0lDjjgADn22GN7VMl+g+mLX/xi6zeZjjjiCN9qOVqdxUOlo+2ZLYWDoBna0oHRojSy1QbE2APQHFOldxKTD9V7gKb+4V6d+/9v70zApSjONfwjhrC5AeJV3FARUAExXuMSruCCGhcWEUUEXKImiglRQTGISxRwjQhoJAoiEVFRUBYBF4SgAm4ohhNXxD0gRhAVceE+f+mcHA5nzvRMdfVUTb/1PDwmp7uq/nq/qpru/rqqXedEG9eECy/fN21s5gDMpcL7gW1O3/qRbXt8zJ8GxphLxe95aehnxaeMueSTBsQCAQj4TSAV5pLfEoQZnYubyjBJ5Bc1F4L58XJ5NlrY0Q1xDkBzzCXMJbtxn8ltM/4Zh/Fo4KIUtHFBNZ4yfdPGZg7AXIqnTxRSim/9qJA2+J4nDYwxl4rfC9PQz4pPGXPJJw2IBQIQ8JsA5pLf+ngbnYubSm8bG2NgXAjGCNOyKLSwAxjiHIDmmEuYS3bjHnMpHn6+lsIc6asyIr5p4+IawLc2+tsbCo8MxoWzi5ozDYwxl6L2BnfnpaGfuaOXf8nwzp8ZOSAAgXQRSIW5tH79ehk3bpzMnDlTli1bJuvWrZOlS5capfWHYuLEidK7d2/Zbbfd0qW+RWtd3FRahBNMVi5M/JEKLey0CHEOQHPMJcwlu3GPuRQPP19LYY70VRnMJX+VCSsyxrh7vdLAGHPJfT/KVUMa+lkuBkkeh3eStKkLAhAIkUDJm0tr1641xpGaSQ0bNpSaNWvKypUrpayszOilx9u1ayc9e/aUiy++OEQNixJziA+WiwKqUqVcmPigAg/Y41AhxDmA8Uffx1yKY/SL2Ix/xmE8GrgoBW1cUI2nTN+0sZkDshHxrY3xKOdXKTB2r0caGGMuue9HuWpIQz/LxSDJ4/BOkjZ1QQACIRIoeXNp6NChZtXSoEGDjIE0cuRIue2228rNJRXt3HPPlRUrVsjkyZND1LAoMbu4qSxKQxKulAuThIFXUx1a2GkR4hyA5phLmEt24z6T22b8Mw7j0cBFKWjjgmo8Zfqmjc0cgLkUT58opBTf+lEhbfA9TxoYYy4VvxemoZ8Vn/J/I4C3T2oQCwQg4COBkjeXOnToIM2aNZPRo0cb/moujRo1aiNz6ZprrpFp06bJggULfNTIy5hc3FR62dCYg+LCJGagFsWhhQU8sVu5YFdz4bnRHHMJc6nw8VMxp801AOMwHg1clII2LqjGU6Zv2tjMAZhL8fSJQkrxrR8V0gbf86SBMeZS8XthGvpZ8SljLvmkAbFAAAJ+Eyh5c6lVq1ZmW7z+/ftnNZeGDRsmEyZMkFdffdVvtTyKzsVNpUfNcxYKF4LO0OZdMFrkjWyjDCHOAWiOuYS5ZDfuM7ltxj/jMB4NXJSCNi6oxlOmb9rYzAGYS/H0iUJK8a0fFdIG3/OkgTHmUvF7YRr6WfEpYy75pAGxQAACfhMoeXOpffv2su+++8ott9yS1Vw666yz5IMPPpBZs2b5rZZH0bm4qfSoec5C4ULQGdq8C0aLvJFhLtkh8yY3fR+TzbYz2lwD0P9s6bvLjzbu2NqW7Js2NnMA5pJtbyg8v2/9qPCW+JszDYwxl4rf/9LQz4pPGXPJJw2IBQIQ8JtAyZtL+q2lRx55RCZNmiTNmzffZFu8RYsWSZ8+fczqpoEDB/qtlkfRubip9Kh5zkLhQtAZ2rwLRou8kWEu2SHzJjd9H3PJtjPaXAPQ/2zpu8uPNu7Y2pbsmzY2cwDmkm1vKDy/b/2o8Jb4mzMNjDGXit//0tDPik8Zc8knDYgFAhDwm0BJmktqEh1xxBFy+OGHyyeffCJdunSRdevWSa9evWT58uUye/ZsufHGG2Xx4sUyceJE2WqrrYwB1bBhQ7/V8ig6FzeVHjXPWShcCDpDm3fBaJE3MswlO2Te5KbvYy7ZdkabawD6ny19d/nRxh1b25J908ZmDsBcsu0Nhef3rR8V3hJ/c6aBMeZS8ftfGvpZ8SljLvmkAbFAAAJ+EyhJc6lFixbSt29f80/T22+/LQMGDJB//vOf5WrUqFFDNmzYIC1btjRG0+677+63Up5F5+Km0rMmOgmHC0EnWAsqFC0KwlaeKcQ5AM0xVSr2evpD4XOAzfiHe+HcXedEG9eECy/fN21s5gDMpcL7gW1O3/qRbXt8zJ8GxphLxe95aehnxaeMueSTBsQCAQj4TSAV5lJGgtdee01effVVWbNmjdSrV09atWplvsdEyp+Ai5vK/KMILwcXgv5ohhZ2WoQ4B6A55hLmkt24z+S2Gf+Mw3g0cFEK2rigGk+ZvmljMwdgLsXTJwopxbd+VEgbfM+TBsaYS8XvhWnoZ8WnjLnkkwbEAgEI+E0gVeaSKym+++47ueOOO+Shhx6SlStXSpMmTeS0006Tnj17iq6Qqi7NmDFDnn76aWN6vfvuu9K4cWOZN29e1izPPfecDB8+XMrKyqROnTrSoUMH6d+/vzRo0GCjPF9//bWMGTNGpk+fLh9++KFsscUWss8++8jvfvc7adOmjTUKFzeV1kEFUAAXgv6IhBZ2WoQ4B6A55hLmkt24x1yKh5+vpTBH+qqMiG/auLgG8K2N/vaGwiODceHsouZMA2PMpai9wd15aehn7ujlXzK882dGDghAIF0EMJdi0HvQoEHy4IMPSvfu3aV169Yyf/58mTlzplxwwQXlW/Nlq0a/A6Urqvbee29jLm222WZZzaVFixbJGWecIc2bN5du3brJZ599ZgykHXbYQSZNmiS1a9cur+bMM8+UBQsWmPPUVFq1apX5vpT+995777U2mFzcVMYghfdFcGHij0RoYadFiHMAmmMuYS7ZjXvMpXj4+VoKc6SvymAu+atMWJExxt3rlQbGmEvu+1GuGtLQz3IxSPI4vJOkTV0QgECIBErWXNLVQ/ovatIVRuPGjYt6evl5uoKoc+fOombOJZdcUv73fv36yZNPPmn+6WqkbOnjjz82x2vWrClqNC1fvjyruaT1rF692qxGqlu3rily7ty5cs4558jAgQPl9NNPN39766235NhjjzVG1KWXXlpetf4odurUyayoGjx4cN5trZghxAfLVg2OKTMXJjGBjKEYtLCDGOIcgOY/ag4HONiNfhGb8U//s6XvLj/auGNrW7Jv2tjMAdlY+NZGW818zA9j96qkgTHmkvt+lKuGNPSzXAySPA7vJGlTFwQgECKBkjWX8hVDzSU1ivJNN998s9kSb86cOWYFUSa9+OKLcuqpp8oVV1xh/hslVWcuLVu2TI4++ugqV0N17NhRttpqK7N6StPixYvl5JNPlgEDBshZZ51VXrWudDrooIPM3/SYTXJxU2kTTyh5uTDxRym0sNMixDkAzX/UHA5wsBv9mEu2/HzNz9zgqzL+zdsurgHof+77H4xhHAcBzKU4KNqVwVi245dvbnjnS4zzIQCBtBEoWXOpT58+0rt377z0zGelU6ZgXbH0xhtvmK3wKqb169ebree6du0q1157baQ4qjOXpk6dKhdffLHceeed0q5du43K07/PmjXLmEq6Auqr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" width="1499.5555555555557">


