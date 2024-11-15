HAIL
====

In this notebook, we used `Hail <https://github.com/hail-is/hail>`__ to
calculate Polygenic Risk Scores (PRS). Hail does not calculate new betas
but instead uses the existing weights from the GWAS and applies a custom
formula for the calculation.

Possible issues
---------------

1. Computation Time: Kindly note that HAIL is computationally expensive
   because it requires data in a specific format (.mt / matrix format),
   and then it loads the data for the computation.

| We followed this tutorial to calculate the PRS:
| https://nbviewer.org/github/ddbj/imputation-server-wf/blob/main/Notebooks/hail-prs-tutorial.ipynb

Basic Process
-------------

1. Read the genotype data.
2. Convert it to VCF.
3. Use Beagle to convert the data to Beagle format.
4. Convert the Beagle format to Hail format.
5. Pass the data to Hail, GWAS, and genotype.
6. Calculate PRS using Hail.

Genotype Data Processing
~~~~~~~~~~~~~~~~~~~~~~~~

1. | **Convert genotype data to VCF format.**
   | Hail requires data in ``beagle.vcf.gz`` format. The first step is
     to convert the ``bed``, ``bim``, and ``fam`` files to VCF format. A
     simple approach is to extract genotype data for each chromosome and
     then convert the data to VCF.

   .. code:: bash

      plink --bfile traindirec/newtrainfilename.clumped.pruned --chr 22 --make-bed --out traindirec/newtrainfilename.clumped.pruned.22
      plink --bfile traindirec/newtrainfilename.clumped.pruned.22 --chr 22 --recode vcf --out traindirec/hail.train.22

2. | **Download the phased reference panel.**
   | To run the following commands, download the phased reference panel
     from 1000 Genomes and place it in the current working directory:
   | `1000 Genomes Reference
     Panel <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/>`__

   .. code:: bash

      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

3. **Alternatively, download the reference panel from Hail.**

   -  | Reference panels:
      | `Hail Reference
        Panel <https://sc.ddbj.nig.ac.jp/en/advanced_guides/imputation_server_tutorial>`__

   -  Files (2.6 GB):

      -  `test-data.GRCh37.vcf.gz <https://zenodo.org/records/6650681/files/test-data.GRCh37.vcf.gz?download=1>`__
         (1.3 GB, ``md5:aff8bca4689cc70f6dbc1c3296590458``)
      -  `test-data.GRCh38.vcf.gz <https://zenodo.org/records/6650681/files/test-data.GRCh38.vcf.gz?download=1>`__
         (1.3 GB, ``md5:d28a741e820444ca926f7b0d5ac2e196``)

4. | **Download genetic distances from the Beagle website.**
   | `Beagle Genetic
     Maps <https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/>`__

5. **Download Beagle.**

   -  Beagle documentation:
      `Beagle
      Documentation <https://faculty.washington.edu/browning/beagle/beagle_5.4_18Mar22.pdf>`__
   -  Beagle download page:
      `Download
      Beagle <https://faculty.washington.edu/browning/beagle/beagle.html>`__

6. | **Run the Beagle command.**
   | Run the following command to perform phasing and imputation:

   .. code:: bash

      java -Xmx50g -jar beagle gt=traindirec/hail.train.22.vcf ref=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr22.GRCh37.map out=traindirec/beagle.hail.train.22

7. | **Follow the remaining process.**
   | The rest of the process is straightforward and can be followed from
     this GitHub repository:
   | `PRS on Hail GitHub
     Repository <https://github.com/hacchy1983/prs-on-hail-public>`__

GWAS file processing for HAIL for Binary Phenotypes.
----------------------------------------------------

When the effect size relates to disease risk and is thus given as an
odds ratio (OR) rather than BETA (for continuous traits), the PRS is
computed as a product of ORs. To simplify this calculation, take the
natural logarithm of the OR so that the PRS can be computed using
summation instead.

.. code:: ipython3

    import os
    import pandas as pd
    import numpy as np
    import sys
    
    #filedirec = sys.argv[1]
    
    filedirec = "SampleData1"
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
    


.. parsed-literal::

    Length of DataFrame! 499617
    Length of DataFrame! 499617
    |    |   chr_name |   chr_position | effect_allele   | other_allele   |   effect_weight |
    |---:|-----------:|---------------:|:----------------|:---------------|----------------:|
    |  0 |          1 |         756604 | A               | G              |     -0.00211532 |
    |  1 |          1 |         768448 | A               | G              |      0.00068708 |
    |  2 |          1 |         779322 | G               | A              |     -0.00239932 |
    |  3 |          1 |         801536 | G               | T              |      0.00203363 |
    |  4 |          1 |         808631 | G               | A              |      0.00130747 |
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
    

Execute HAIL
------------

.. code:: ipython3

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
    


.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC
      --geno 0.001
      --hwe 1e-05
      --indep-pairwise 200 50 0.25
      --maf 0.2
      --out SampleData1/Fold_0/train_data
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    491952 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999894.
    4954 variants removed due to missing genotype data (--geno).
    --hwe: 0 variants removed due to Hardy-Weinberg exact test.
    346375 variants removed due to minor allele threshold(s)
    (--maf/--max-maf/--mac/--max-mac).
    140623 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    Pruned 6339 variants from chromosome 1, leaving 4777.
    Pruned 6277 variants from chromosome 2, leaving 4595.
    Pruned 5191 variants from chromosome 3, leaving 3992.
    Pruned 4996 variants from chromosome 4, leaving 3727.
    Pruned 4537 variants from chromosome 5, leaving 3581.
    Pruned 6868 variants from chromosome 6, leaving 3561.
    Pruned 4379 variants from chromosome 7, leaving 3276.
    Pruned 4071 variants from chromosome 8, leaving 2898.
    Pruned 3287 variants from chromosome 9, leaving 2823.
    Pruned 3822 variants from chromosome 10, leaving 3054.
    Pruned 4158 variants from chromosome 11, leaving 2849.
    Pruned 3544 variants from chromosome 12, leaving 2868.
    Pruned 2560 variants from chromosome 13, leaving 2246.
    Pruned 2424 variants from chromosome 14, leaving 1989.
    Pruned 2615 variants from chromosome 15, leaving 2060.
    Pruned 2867 variants from chromosome 16, leaving 2231.
    Pruned 2782 variants from chromosome 17, leaving 2164.
    Pruned 2223 variants from chromosome 18, leaving 2047.
    Pruned 2593 variants from chromosome 19, leaving 1818.
    Pruned 1980 variants from chromosome 20, leaving 1779.
    Pruned 1254 variants from chromosome 21, leaving 992.
    Pruned 1403 variants from chromosome 22, leaving 1126.
    Pruning complete.  80170 of 140623 variants removed.
    Marker lists written to SampleData1/Fold_0/train_data.prune.in and
    SampleData1/Fold_0/train_data.prune.out .
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC
      --clump SampleData1/SampleData1.txt
      --clump-field P
      --clump-kb 200
      --clump-p1 1
      --clump-r2 0.1
      --clump-snp-field SNP
      --extract SampleData1/Fold_0/train_data.prune.in
      --out SampleData1/Fold_0/train_data
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    491952 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 60453 variants remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is exactly 1.
    60453 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --clump: 38646 clumps formed from 60453 top variants.
    Results written to SampleData1/Fold_0/train_data.clumped .
    

.. parsed-literal::

    Warning: 'rs3134762' is missing from the main dataset, and is a top variant.
    Warning: 'rs3132505' is missing from the main dataset, and is a top variant.
    Warning: 'rs3130424' is missing from the main dataset, and is a top variant.
    439161 more top variant IDs missing; see log file.
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC
      --extract SampleData1/Fold_0/train_data.valid.snp
      --indep-pairwise 200 50 0.25
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    491952 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 77998 variants remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999939.
    77998 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
    Pruned 2 variants from chromosome 1, leaving 6245.
    Pruned 2 variants from chromosome 2, leaving 5993.
    Pruned 1 variant from chromosome 3, leaving 5253.
    Pruned 0 variants from chromosome 4, leaving 4966.
    Pruned 4 variants from chromosome 5, leaving 4812.
    99%

.. parsed-literal::

    Warning: At least 5268 duplicate IDs in --extract file.
    

.. parsed-literal::

    Pruned 3 variants from chromosome 6, leaving 4414.
    Pruned 4 variants from chromosome 7, leaving 4212.
    Pruned 2 variants from chromosome 8, leaving 3866.
    Pruned 1 variant from chromosome 9, leaving 3555.
    Pruned 3 variants from chromosome 10, leaving 3949.
    Pruned 1 variant from chromosome 11, leaving 3758.
    Pruned 1 variant from chromosome 12, leaving 3715.
    Pruned 2 variants from chromosome 13, leaving 2866.
    Pruned 1 variant from chromosome 14, leaving 2590.
    Pruned 2 variants from chromosome 15, leaving 2604.
    Pruned 2 variants from chromosome 16, leaving 2870.
    Pruned 0 variants from chromosome 17, leaving 2668.
    Pruned 0 variants from chromosome 18, leaving 2603.
    Pruned 0 variants from chromosome 19, leaving 2136.
    Pruned 1 variant from chromosome 20, leaving 2239.
    Pruned 0 variants from chromosome 21, leaving 1331.
    Pruned 0 variants from chromosome 22, leaving 1321.
    Pruning complete.  32 of 77998 variants removed.
    Marker lists written to
    SampleData1/Fold_0/train_data.QC.clumped.pruned.prune.in and
    SampleData1/Fold_0/train_data.QC.clumped.pruned.prune.out .
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data
      --extract SampleData1/Fold_0/train_data.valid.snp
      --indep-pairwise 200 50 0.25
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    551892 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    --extract: 77998 variants remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999939.
    77998 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    Pruned 402 variants from chromosome 1, leaving 5845.
    Pruned 358 variants from chromosome 2, leaving 5637.
    Pruned 333 variants from chromosome 3, leaving 4921.
    Pruned 317 variants from chromosome 4, leaving 4649.
    Pruned 329 variants from chromosome 5, leaving 4487.
    Pruned 270 variants from chromosome 6, leaving 4147.
    Pruned 284 variants from chromosome 7, leaving 3932.
    Pruned 255 variants from chromosome 8, leaving 3613.
    Pruned 208 variants from chromosome 9, leaving 3348.
    Pruned 229 variants from chromosome 10, leaving 3723.
    Pruned 236 variants from chromosome 11, leaving 3523.
    Pruned 247 variants from chromosome 12, leaving 3469.
    81%

.. parsed-literal::

    Warning: At least 5268 duplicate IDs in --extract file.
    

.. parsed-literal::

    Pruned 170 variants from chromosome 13, leaving 2698.
    Pruned 124 variants from chromosome 14, leaving 2467.
    Pruned 137 variants from chromosome 15, leaving 2469.
    Pruned 170 variants from chromosome 16, leaving 2702.
    Pruned 124 variants from chromosome 17, leaving 2544.
    Pruned 134 variants from chromosome 18, leaving 2469.
    Pruned 98 variants from chromosome 19, leaving 2038.
    Pruned 147 variants from chromosome 20, leaving 2093.
    Pruned 68 variants from chromosome 21, leaving 1263.
    Pruned 57 variants from chromosome 22, leaving 1264.
    Pruning complete.  4697 of 77998 variants removed.
    Marker lists written to SampleData1/Fold_0/test_data.clumped.pruned.prune.in
    and SampleData1/Fold_0/test_data.clumped.pruned.prune.out .
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --extract SampleData1/Fold_0/train_data.valid.snp
      --out SampleData1/Fold_0/test_data
      --pca 6
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    77998 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    --extract: 77998 variants remaining.
    Using up to 8 threads (change this with --threads).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999939.
    77998 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    53700 markers complete.

.. parsed-literal::

    Warning: At least 5268 duplicate IDs in --extract file.
    

.. parsed-literal::

    Relationship matrix calculation complete.
    --pca: Results saved to SampleData1/Fold_0/test_data.eigenval and
    SampleData1/Fold_0/test_data.eigenvec .
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --extract SampleData1/Fold_0/train_data.valid.snp
      --out SampleData1/Fold_0/train_data
      --pca 6
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    77998 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    --extract: 77998 variants remaining.
    Using up to 8 threads (change this with --threads).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999939.
    77998 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    5700 markers complete.

.. parsed-literal::

    Warning: At least 5268 duplicate IDs in --extract file.
    

.. parsed-literal::

    Relationship matrix calculation complete.
    --pca: Results saved to SampleData1/Fold_0/train_data.eigenval and
    SampleData1/Fold_0/train_data.eigenvec .
    Deleted directory: SampleData1/Fold_0/hail.train.1.mt
    Deleted directory: SampleData1/Fold_0/hail.train.2.mt
    Deleted directory: SampleData1/Fold_0/hail.train.3.mt
    Deleted directory: SampleData1/Fold_0/hail.train.4.mt
    Deleted directory: SampleData1/Fold_0/hail.train.5.mt
    Deleted directory: SampleData1/Fold_0/hail.train.6.mt
    Deleted directory: SampleData1/Fold_0/hail.train.7.mt
    Deleted directory: SampleData1/Fold_0/hail.train.8.mt
    Deleted directory: SampleData1/Fold_0/hail.train.9.mt
    Deleted directory: SampleData1/Fold_0/hail.train.10.mt
    Deleted directory: SampleData1/Fold_0/hail.train.11.mt
    Deleted directory: SampleData1/Fold_0/hail.train.12.mt
    Deleted directory: SampleData1/Fold_0/hail.train.13.mt
    Deleted directory: SampleData1/Fold_0/hail.train.14.mt
    Deleted directory: SampleData1/Fold_0/hail.train.15.mt
    Deleted directory: SampleData1/Fold_0/hail.train.16.mt
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.1.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 1
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.1
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    6247 out of 77998 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999944.
    6247 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.1.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.1.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.1.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.1.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.1
      --chr 1
      --out SampleData1/Fold_0/hail.train.1
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    6247 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999944.
    6247 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.1.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 01:49 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.1.vcf
      ref=ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr1.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.1
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [1:10177-18871145]
    Reference markers:              579,450
    Study     markers:                  554
    
    Burnin  iteration 1:           4 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  7735
    Estimated err:                 6.6e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               33 seconds
    
    Window 2 [1:18270328-53836747]
    Reference markers:            1,004,286
    Study     markers:                  721
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  6790
    Estimated err:                 8.6e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               22 seconds
    
    Window 3 [1:49070327-87760656]
    Reference markers:            1,094,455
    Study     markers:                  758
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  7514
    Estimated err:                 7.9e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               28 seconds
    
    Window 4 [1:85578740-144962397]
    Reference markers:            1,037,157
    Study     markers:                  705
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  6710
    Estimated err:                 2.0e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               3 minutes 29 seconds
    
    Window 5 [1:144224289-175819109]
    Reference markers:              830,857
    Study     markers:                  578
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  7663
    Estimated err:                 8.2e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               1 minute 15 seconds
    
    Window 6 [1:173093312-213583529]
    Reference markers:            1,151,258
    Study     markers:                  766
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  7897
    Estimated err:                 7.5e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               20 seconds
    
    Window 7 [1:211080854-240239399]
    Reference markers:              854,129
    Study     markers:                  715
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  5728
    Estimated err:                 7.5e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               15 seconds
    
    Window 8 [1:239251841-249240543]
    Reference markers:              306,271
    Study     markers:                  296
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  32903
    Estimated err:                 1.2e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    

.. parsed-literal::

    Initializing Hail with default parameters...
    

.. parsed-literal::

    
    Imputation time:               9 seconds
    
    Cumulative Statistics:
    
    Reference markers:            6,468,094
    Study     markers:                4,870
    
    Haplotype phasing time:        1 minute 34 seconds
    Imputation time:               6 minutes 50 seconds
    Total time:                    9 minutes 13 seconds
    
    End time: 01:58 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.1.vcf ref=ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr1.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.1
    

.. parsed-literal::

    SLF4J: Failed to load class "org.slf4j.impl.StaticLoggerBinder".
    SLF4J: Defaulting to no-operation (NOP) logger implementation
    SLF4J: See http://www.slf4j.org/codes.html#StaticLoggerBinder for further details.
    Running on Apache Spark version 3.5.1
    SparkUI available at http://login02:4040
    Welcome to
         __  __     <>__
        / /_/ /__  __/ /
       / __  / _ `/ / /
      /_/ /_/\_,_/_/_/   version 0.2.132-678e1f52b999
    LOGGING: writing to /data/ascher01/uqmmune1/BenchmarkingPGSTools/hail-20240927-1358-0.2.132-678e1f52b999.log
    2024-09-27 13:58:45.333 Hail: INFO: scanning VCF for sortedness...
    SLF4J: Failed to load class "org.slf4j.impl.StaticMDCBinder".
    SLF4J: Defaulting to no-operation MDCAdapter implementation.
    SLF4J: See http://www.slf4j.org/codes.html#no_static_mdc_binder for further details.
    2024-09-27 13:59:05.274 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 14:04:57.046 Hail: INFO: wrote matrix table with 6468094 rows and 380 columns in 2 partitions to SampleData1/Fold_0/hail.train.1.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.2.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 2
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.2
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8998 out of 114118 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999929.
    8998 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.2.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.2.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.2.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.2.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.2
      --chr 2
      --out SampleData1/Fold_0/hail.train.2
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8998 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999929.
    8998 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.2.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 02:04 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.2.vcf
      ref=ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr2.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.2
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [2:10179-18145025]
    Reference markers:              573,772
    Study     markers:                  880
    
    Burnin  iteration 1:           5 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  7544
    Estimated err:                 1.1e-03
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               25 seconds
    
    Window 2 [2:16417586-49114274]
    Reference markers:            1,053,726
    Study     markers:                1,173
    
    Burnin  iteration 1:           4 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           3 seconds
    
    Estimated ne:                  10035
    Estimated err:                 8.6e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               32 seconds
    
    Window 3 [2:47323299-88186309]
    Reference markers:            1,255,428
    Study     markers:                1,184
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           3 seconds
    
    Estimated ne:                  8361
    Estimated err:                 8.7e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               33 seconds
    
    Window 4 [2:86643697-134624574]
    Reference markers:            1,221,295
    Study     markers:                1,092
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  8660
    Estimated err:                 2.2e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               45 seconds
    
    Window 5 [2:133390491-177333395]
    Reference markers:            1,235,397
    Study     markers:                1,202
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  8841
    Estimated err:                 8.2e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           2 seconds
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               34 seconds
    
    Window 6 [2:175084768-220694589]
    Reference markers:            1,278,696
    Study     markers:                1,258
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           3 seconds
    Burnin  iteration 3:           3 seconds
    
    Estimated ne:                  9252
    Estimated err:                 2.1e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           2 seconds
    Phasing iteration 6:           2 seconds
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          2 seconds
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               35 seconds
    
    Window 7 [2:218829611-242575162]
    Reference markers:              734,807
    Study     markers:                1,103
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  9018
    Estimated err:                 8.8e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          0 seconds
    
    Imputation time:               16 seconds
    
    Window 8 [2:241170836-243188367]
    Reference markers:               65,148
    Study     markers:                   65
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  245134
    Estimated err:                 1.2e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               4 seconds
    
    Cumulative Statistics:
    
    Reference markers:            7,081,600
    Study     markers:                7,568
    
    Haplotype phasing time:        2 minutes 22 seconds
    Imputation time:               3 minutes 44 seconds
    Total time:                    7 minutes 49 seconds
    
    End time: 02:12 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.2.vcf ref=ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr2.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.2
    

.. parsed-literal::

    2024-09-27 14:12:47.026 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 14:13:05.331 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 14:18:18.220 Hail: INFO: wrote matrix table with 7081600 rows and 380 columns in 2 partitions to SampleData1/Fold_0/hail.train.2.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.3.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 3
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.3
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    7649 out of 114118 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999916.
    7649 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.3.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.3.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.3.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.3.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.3
      --chr 3
      --out SampleData1/Fold_0/hail.train.3
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    7649 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999916.
    7649 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.3.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 02:18 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.3.vcf
      ref=ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr3.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.3
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [3:60069-21643990]
    Reference markers:              718,343
    Study     markers:                1,014
    
    Burnin  iteration 1:           5 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  7964
    Estimated err:                 9.1e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               29 seconds
    
    Window 2 [3:19789059-60278775]
    Reference markers:            1,194,873
    Study     markers:                1,173
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  8022
    Estimated err:                 9.3e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               37 seconds
    
    Window 3 [3:59462914-105372608]
    Reference markers:            1,325,897
    Study     markers:                1,253
    
    Burnin  iteration 1:           4 seconds
    Burnin  iteration 2:           3 seconds
    Burnin  iteration 3:           3 seconds
    
    Estimated ne:                  8492
    Estimated err:                 8.9e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               36 seconds
    
    Window 4 [3:103398092-143445761]
    Reference markers:            1,153,131
    Study     markers:                1,258
    
    Burnin  iteration 1:           4 seconds
    Burnin  iteration 2:           3 seconds
    Burnin  iteration 3:           3 seconds
    
    Estimated ne:                  8060
    Estimated err:                 8.4e-04
    
    Phasing iteration 1:           2 seconds
    Phasing iteration 2:           2 seconds
    Phasing iteration 3:           2 seconds
    Phasing iteration 4:           1 second
    Phasing iteration 5:           2 seconds
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               32 seconds
    
    Window 5 [3:141741447-182238285]
    Reference markers:            1,155,993
    Study     markers:                1,208
    
    Burnin  iteration 1:           4 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           3 seconds
    
    Estimated ne:                  7660
    Estimated err:                 8.4e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               27 seconds
    
    Window 6 [3:179790776-197962381]
    Reference markers:              545,816
    Study     markers:                  848
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  8667
    Estimated err:                 8.6e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           1 second
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               13 seconds
    
    Cumulative Statistics:
    
    Reference markers:            5,832,276
    Study     markers:                6,419
    
    Haplotype phasing time:        2 minutes 1 second
    Imputation time:               2 minutes 53 seconds
    Total time:                    6 minutes 14 seconds
    
    End time: 02:24 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.3.vcf ref=ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr3.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.3
    

.. parsed-literal::

    2024-09-27 14:24:32.872 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 14:24:57.642 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 14:33:22.829 Hail: INFO: wrote matrix table with 5832276 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.3.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.4.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 4
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.4
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    7312 out of 114118 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999935.
    7312 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.4.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.4.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.4.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.4.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.4
      --chr 4
      --out SampleData1/Fold_0/hail.train.4
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    7312 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999935.
    7312 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.4.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 02:33 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.4.vcf
      ref=ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr4.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.4
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [4:10005-22520039]
    Reference markers:              786,334
    Study     markers:                  973
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  8066
    Estimated err:                 9.2e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          1 second
    Phasing iteration 12:          0 seconds
    
    Imputation time:               18 seconds
    
    Window 2 [4:20544690-60086601]
    Reference markers:            1,127,621
    Study     markers:                1,186
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  9301
    Estimated err:                 8.4e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               19 seconds
    
    Window 3 [4:57715343-105764064]
    Reference markers:            1,443,530
    Study     markers:                1,220
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  8989
    Estimated err:                 2.0e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               23 seconds
    
    Window 4 [4:102709973-150650017]
    Reference markers:            1,376,146
    Study     markers:                1,281
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  9225
    Estimated err:                 8.8e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               22 seconds
    
    Window 5 [4:148086441-182155049]
    Reference markers:            1,012,842
    Study     markers:                1,191
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  8541
    Estimated err:                 7.8e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          1 second
    Phasing iteration 12:          0 seconds
    
    Imputation time:               13 seconds
    
    Window 6 [4:181312927-191043881]
    Reference markers:              316,589
    Study     markers:                  612
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  10400
    Estimated err:                 9.5e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               6 seconds
    
    Cumulative Statistics:
    
    Reference markers:            5,732,585
    Study     markers:                6,121
    
    Haplotype phasing time:        1 minute 13 seconds
    Imputation time:               1 minute 41 seconds
    Total time:                    4 minutes 6 seconds
    
    End time: 02:37 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.4.vcf ref=ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr4.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.4
    

.. parsed-literal::

    2024-09-27 14:37:29.240 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 14:37:52.712 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 14:45:55.265 Hail: INFO: wrote matrix table with 5732585 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.4.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.5.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 5
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.5
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    10632 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999882.
    10632 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.5.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.5.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.5.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.5.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.5
      --chr 5
      --out SampleData1/Fold_0/hail.train.5
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    10632 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999882.
    10632 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.5.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 02:45 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.5.vcf
      ref=ALL.chr5.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr5.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.5
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [5:10043-24839896]
    Reference markers:              785,833
    Study     markers:                1,656
    
    Burnin  iteration 1:           4 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  11105
    Estimated err:                 2.0e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               25 seconds
    
    Window 2 [5:22775078-69306446]
    Reference markers:            1,297,643
    Study     markers:                1,974
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  19142
    Estimated err:                 2.2e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               29 seconds
    
    Window 3 [5:67534999-111982027]
    Reference markers:            1,223,992
    Study     markers:                2,012
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           3 seconds
    
    Estimated ne:                  19277
    Estimated err:                 2.0e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               30 seconds
    
    Window 4 [5:109533719-149269607]
    Reference markers:            1,201,670
    Study     markers:                1,968
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           3 seconds
    
    Estimated ne:                  18121
    Estimated err:                 1.9e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               30 seconds
    
    Window 5 [5:148226432-174766607]
    Reference markers:              787,384
    Study     markers:                1,787
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  12549
    Estimated err:                 2.1e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               15 seconds
    
    Window 6 [5:174167120-180904689]
    Reference markers:              201,966
    Study     markers:                  480
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  7643
    Estimated err:                 9.3e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               4 seconds
    
    Cumulative Statistics:
    
    Reference markers:            5,265,763
    Study     markers:                9,450
    
    Haplotype phasing time:        1 minute 33 seconds
    Imputation time:               2 minutes 13 seconds
    Total time:                    4 minutes 54 seconds
    
    End time: 02:50 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.5.vcf ref=ALL.chr5.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr5.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.5
    

.. parsed-literal::

    2024-09-27 14:50:49.927 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 14:51:13.675 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 14:59:00.308 Hail: INFO: wrote matrix table with 5265763 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.5.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.6.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 6
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.6
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    10068 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999902.
    10068 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.6.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.6.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.6.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.6.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.6
      --chr 6
      --out SampleData1/Fold_0/hail.train.6
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    10068 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999902.
    10068 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.6.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 02:59 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.6.vcf
      ref=ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr6.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.6
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [6:63854-21294216]
    Reference markers:              658,963
    Study     markers:                1,608
    
    Burnin  iteration 1:           4 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  11789
    Estimated err:                 2.2e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          0 seconds
    
    Imputation time:               24 seconds
    
    Window 2 [6:20005714-56895330]
    Reference markers:            1,115,840
    Study     markers:                1,961
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  15522
    Estimated err:                 6.4e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               25 seconds
    
    Window 3 [6:53556519-110284141]
    Reference markers:            1,605,180
    Study     markers:                2,128
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  18780
    Estimated err:                 1.9e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               31 seconds
    
    Window 4 [6:108297729-149276979]
    Reference markers:            1,157,048
    Study     markers:                2,014
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  19637
    Estimated err:                 2.0e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               22 seconds
    
    Window 5 [6:148481212-170715177]
    Reference markers:              694,787
    Study     markers:                1,580
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  11587
    Estimated err:                 2.1e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               12 seconds
    
    Window 6 [6:169508794-171053406]
    Reference markers:               46,342
    Study     markers:                   86
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  209543
    Estimated err:                 5.0e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               1 second
    
    Cumulative Statistics:
    
    Reference markers:            5,024,119
    Study     markers:                8,904
    
    Haplotype phasing time:        1 minute 17 seconds
    Imputation time:               1 minute 55 seconds
    Total time:                    4 minutes 2 seconds
    
    End time: 03:03 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.6.vcf ref=ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr6.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.6
    

.. parsed-literal::

    2024-09-27 15:03:03.199 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 15:03:23.648 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 15:10:43.593 Hail: INFO: wrote matrix table with 5024119 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.6.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.7.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 7
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.7
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    9496 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999897.
    9496 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.7.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.7.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.7.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.7.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.7
      --chr 7
      --out SampleData1/Fold_0/hail.train.7
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    9496 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999897.
    9496 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.7.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 03:10 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.7.vcf
      ref=ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr7.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.7
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [7:14808-22507950]
    Reference markers:              834,932
    Study     markers:                1,527
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  10288
    Estimated err:                 2.0e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               23 seconds
    
    Window 2 [7:21712814-53284120]
    Reference markers:              941,559
    Study     markers:                1,889
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  20441
    Estimated err:                 2.2e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               22 seconds
    
    Window 3 [7:51371519-102207161]
    Reference markers:            1,373,374
    Study     markers:                1,945
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  17409
    Estimated err:                 6.5e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               28 seconds
    
    Window 4 [7:98972040-141848611]
    Reference markers:            1,207,792
    Study     markers:                2,088
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           3 seconds
    
    Estimated ne:                  19500
    Estimated err:                 2.0e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               21 seconds
    
    Window 5 [7:140215356-159128567]
    Reference markers:              589,076
    Study     markers:                1,342
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  8786
    Estimated err:                 9.1e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               11 seconds
    
    Cumulative Statistics:
    
    Reference markers:            4,716,715
    Study     markers:                8,421
    
    Haplotype phasing time:        1 minute 17 seconds
    Imputation time:               1 minute 45 seconds
    Total time:                    3 minutes 42 seconds
    
    End time: 03:14 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.7.vcf ref=ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr7.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.7
    

.. parsed-literal::

    2024-09-27 15:14:26.167 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 15:14:46.616 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 15:21:32.811 Hail: INFO: wrote matrix table with 4716715 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.7.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.8.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 8
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.8
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8867 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    8867 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.8.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.8.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.8.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.8.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.8
      --chr 8
      --out SampleData1/Fold_0/hail.train.8
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8867 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    8867 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.8.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 03:21 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.8.vcf
      ref=ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr8.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.8
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [8:11740-18799855]
    Reference markers:              877,773
    Study     markers:                1,473
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  9532
    Estimated err:                 2.1e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               24 seconds
    
    Window 2 [8:17678413-59729351]
    Reference markers:            1,232,933
    Study     markers:                1,985
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  20310
    Estimated err:                 2.0e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               27 seconds
    

.. parsed-literal::

    java.lang.IllegalArgumentException: Duplicate marker: 8	89329148	.	G	A
    	at vcf.Markers.markerSet(Markers.java:129)
    	at vcf.Markers.<init>(Markers.java:83)
    	at vcf.Markers.create(Markers.java:62)
    	at vcf.RefGT.<init>(RefGT.java:90)
    	at vcf.RefTargSlidingWindow$Reader.window(RefTargSlidingWindow.java:326)
    	at vcf.RefTargSlidingWindow$Reader.readWindow(RefTargSlidingWindow.java:303)
    	at vcf.RefTargSlidingWindow$Reader.run(RefTargSlidingWindow.java:246)
    	at java.base/java.lang.Thread.run(Thread.java:829)
    
    Terminating program.
    

.. parsed-literal::

    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.8.vcf ref=ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr8.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.8
    

.. parsed-literal::

    2024-09-27 15:24:02.637 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 15:24:11.134 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 15:27:05.495 Hail: INFO: wrote matrix table with 2039669 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.8.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.9.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 9
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.9
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    7768 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999892.
    7768 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.9.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.9.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.9.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.9.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.9
      --chr 9
      --out SampleData1/Fold_0/hail.train.9
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    7768 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999892.
    7768 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.9.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 03:27 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.9.vcf
      ref=ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr9.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.9
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [9:10163-19603557]
    Reference markers:              766,446
    Study     markers:                1,563
    
    Burnin  iteration 1:           4 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  11323
    Estimated err:                 2.2e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               24 seconds
    
    Window 2 [9:18636347-80220415]
    Reference markers:            1,009,966
    Study     markers:                1,724
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  16522
    Estimated err:                 1.9e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               27 seconds
    
    Window 3 [9:79151500-112369214]
    Reference markers:              978,438
    Study     markers:                1,899
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  18394
    Estimated err:                 2.1e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               24 seconds
    
    Window 4 [9:110879571-137041122]
    Reference markers:              775,425
    Study     markers:                1,621
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  10878
    Estimated err:                 2.2e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               13 seconds
    
    Window 5 [9:136539237-141148854]
    Reference markers:              166,298
    Study     markers:                  363
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  13390
    Estimated err:                 8.8e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               21 seconds
    
    Cumulative Statistics:
    
    Reference markers:            3,560,687
    Study     markers:                6,877
    
    Haplotype phasing time:        1 minute 5 seconds
    Imputation time:               1 minute 49 seconds
    Total time:                    3 minutes 31 seconds
    
    End time: 03:30 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.9.vcf ref=ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr9.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.9
    

.. parsed-literal::

    2024-09-27 15:30:37.280 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 15:30:52.044 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 15:35:56.524 Hail: INFO: wrote matrix table with 3560687 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.9.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.10.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 10
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.10
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8824 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999901.
    8824 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.10.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.10.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.10.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.10.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.10
      --chr 10
      --out SampleData1/Fold_0/hail.train.10
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8824 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999901.
    8824 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.10.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 03:35 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.10.vcf
      ref=ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr10.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.10
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [10:60494-18538847]
    Reference markers:              622,290
    Study     markers:                1,543
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  10772
    Estimated err:                 2.2e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               21 seconds
    
    Window 2 [10:17068251-54798179]
    Reference markers:            1,000,400
    Study     markers:                1,560
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  16442
    Estimated err:                 8.9e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               15 seconds
    
    Window 3 [10:53792357-90734239]
    Reference markers:            1,092,804
    Study     markers:                1,942
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  17386
    Estimated err:                 7.0e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          0 seconds
    
    Imputation time:               14 seconds
    
    Window 4 [10:88784362-123900556]
    Reference markers:            1,013,192
    Study     markers:                1,904
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  21656
    Estimated err:                 2.4e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          0 seconds
    
    Imputation time:               10 seconds
    
    Window 5 [10:123093921-135524373]
    Reference markers:              406,013
    Study     markers:                1,118
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  10684
    Estimated err:                 9.8e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               5 seconds
    
    Cumulative Statistics:
    
    Reference markers:            3,992,219
    Study     markers:                7,763
    
    Haplotype phasing time:        49 seconds
    Imputation time:               1 minute 5 seconds
    Total time:                    2 minutes 29 seconds
    
    End time: 03:38 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.10.vcf ref=ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr10.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.10
    

.. parsed-literal::

    2024-09-27 15:38:26.409 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 15:38:42.411 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 15:44:28.695 Hail: INFO: wrote matrix table with 3992219 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.10.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.11.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 11
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.11
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8420 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.9999.
    8420 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.11.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.11.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.11.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.11.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.11
      --chr 11
      --out SampleData1/Fold_0/hail.train.11
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8420 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.9999.
    8420 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.11.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 03:44 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.11.vcf
      ref=ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr11.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.11
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [11:61395-22535739]
    Reference markers:              711,691
    Study     markers:                1,629
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  11525
    Estimated err:                 2.0e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               22 seconds
    
    Window 2 [11:21349458-70803324]
    Reference markers:            1,428,012
    Study     markers:                2,035
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  18794
    Estimated err:                 1.9e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          0 seconds
    
    Imputation time:               33 seconds
    
    Window 3 [11:69604836-112290497]
    Reference markers:            1,283,964
    Study     markers:                2,125
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  21225
    Estimated err:                 1.9e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               25 seconds
    
    Window 4 [11:110168424-133043804]
    Reference markers:              692,681
    Study     markers:                1,740
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  18833
    Estimated err:                 2.1e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           1 second
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          1 second
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               12 seconds
    
    Window 5 [11:132350788-134946453]
    Reference markers:               88,218
    Study     markers:                  273
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  30941
    Estimated err:                 7.5e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               1 second
    
    Cumulative Statistics:
    
    Reference markers:            4,045,628
    Study     markers:                7,460
    
    Haplotype phasing time:        54 seconds
    Imputation time:               1 minute 34 seconds
    Total time:                    3 minutes 23 seconds
    
    End time: 03:47 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.11.vcf ref=ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr11.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.11
    

.. parsed-literal::

    2024-09-27 15:47:52.751 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 15:48:09.331 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 15:53:57.158 Hail: INFO: wrote matrix table with 4045628 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.11.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.12.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 12
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.12
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8198 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    8198 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.12.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.12.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.12.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.12.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.12
      --chr 12
      --out SampleData1/Fold_0/hail.train.12
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8198 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    8198 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.12.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 03:53 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.12.vcf
      ref=ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr12.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.12
      nthreads=8
    

.. parsed-literal::

    java.lang.IllegalArgumentException: Duplicate marker: 12	8400000	.	T	G
    	at vcf.Markers.markerSet(Markers.java:129)
    	at vcf.Markers.<init>(Markers.java:83)
    	at vcf.Markers.create(Markers.java:62)
    	at vcf.RefGT.<init>(RefGT.java:90)
    	at vcf.RefTargSlidingWindow$Reader.window(RefTargSlidingWindow.java:326)
    	at vcf.RefTargSlidingWindow$Reader.readWindow(RefTargSlidingWindow.java:303)
    	at vcf.RefTargSlidingWindow$Reader.run(RefTargSlidingWindow.java:246)
    	at java.base/java.lang.Thread.run(Thread.java:829)
    
    Terminating program.
    

.. parsed-literal::

    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.12.vcf ref=ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr12.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.12
    

.. parsed-literal::

    2024-09-27 15:54:16.909 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 15:54:17.167 Hail: INFO: Coerced sorted VCF - no additional import work to do
    2024-09-27 15:54:17.658 Hail: INFO: wrote matrix table with 0 rows and 380 columns in 0 partitions to SampleData1/Fold_0/hail.train.12.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.13.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 13
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.13
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    6350 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999901.
    6350 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.13.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.13.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.13.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.13.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.13
      --chr 13
      --out SampleData1/Fold_0/hail.train.13
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    6350 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999901.
    6350 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.13.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 03:54 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.13.vcf
      ref=ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr13.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.13
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [13:19020047-40892639]
    Reference markers:              671,822
    Study     markers:                1,483
    
    Burnin  iteration 1:           4 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  9869
    Estimated err:                 2.2e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               20 seconds
    
    Window 2 [13:39849087-81726083]
    Reference markers:            1,198,543
    Study     markers:                1,992
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  19344
    Estimated err:                 2.0e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               25 seconds
    
    Window 3 [13:79491658-110702359]
    Reference markers:              937,844
    Study     markers:                1,855
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  12534
    Estimated err:                 2.0e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           1 second
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           1 second
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          1 second
    Phasing iteration 12:          0 seconds
    
    Imputation time:               14 seconds
    
    Window 4 [13:109580155-115109852]
    Reference markers:              180,096
    Study     markers:                  461
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  8263
    Estimated err:                 1.3e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               3 seconds
    
    Cumulative Statistics:
    
    Reference markers:            2,857,916
    Study     markers:                5,502
    
    Haplotype phasing time:        49 seconds
    Imputation time:               1 minute 2 seconds
    Total time:                    2 minutes 22 seconds
    
    End time: 03:56 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.13.vcf ref=ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr13.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.13
    

.. parsed-literal::

    2024-09-27 15:56:40.644 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 15:56:53.103 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 16:01:03.478 Hail: INFO: wrote matrix table with 2857916 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.13.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.14.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 14
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.14
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5742 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999852.
    5742 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.14.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.14.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.14.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.14.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.14
      --chr 14
      --out SampleData1/Fold_0/hail.train.14
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5742 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999852.
    5742 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.14.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 04:01 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.14.vcf
      ref=ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr14.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.14
      nthreads=8
    

.. parsed-literal::

    java.lang.IllegalArgumentException: Duplicate marker: 14	21649957	.	C	T
    	at vcf.Markers.markerSet(Markers.java:129)
    	at vcf.Markers.<init>(Markers.java:83)
    	at vcf.Markers.create(Markers.java:62)
    	at vcf.RefGT.<init>(RefGT.java:90)
    	at vcf.RefTargSlidingWindow$Reader.window(RefTargSlidingWindow.java:326)
    	at vcf.RefTargSlidingWindow$Reader.readWindow(RefTargSlidingWindow.java:303)
    	at vcf.RefTargSlidingWindow$Reader.run(RefTargSlidingWindow.java:246)
    	at java.base/java.lang.Thread.run(Thread.java:829)
    
    Terminating program.
    

.. parsed-literal::

    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.14.vcf ref=ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr14.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.14
    

.. parsed-literal::

    2024-09-27 16:01:21.924 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 16:01:22.171 Hail: INFO: Coerced sorted VCF - no additional import work to do
    2024-09-27 16:01:22.527 Hail: INFO: wrote matrix table with 0 rows and 380 columns in 0 partitions to SampleData1/Fold_0/hail.train.14.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.15.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 15
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.15
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5569 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999887.
    5569 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.15.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.15.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.15.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.15.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.15
      --chr 15
      --out SampleData1/Fold_0/hail.train.15
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5569 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999887.
    5569 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.15.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 04:01 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.15.vcf
      ref=ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr15.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.15
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [15:20000041-32890838]
    Reference markers:              325,223
    Study     markers:                  628
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  6297
    Estimated err:                 3.1e-03
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               2 minutes 42 seconds
    
    Window 2 [15:32564155-62934668]
    Reference markers:              905,216
    Study     markers:                1,635
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  18800
    Estimated err:                 2.3e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               13 seconds
    
    Window 3 [15:61592413-93656456]
    Reference markers:              943,933
    Study     markers:                1,797
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  19673
    Estimated err:                 2.1e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           1 second
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          1 second
    Phasing iteration 12:          0 seconds
    
    Imputation time:               14 seconds
    
    Window 4 [15:92986651-102521131]
    Reference markers:              318,635
    Study     markers:                  976
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  10402
    Estimated err:                 1.0e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Cumulative Statistics:
    
    Reference markers:            2,424,689
    Study     markers:                4,877
    
    Haplotype phasing time:        38 seconds
    Imputation time:               3 minutes 17 seconds
    Total time:                    4 minutes 7 seconds
    
    End time: 04:05 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.15.vcf ref=ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr15.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.15
    

.. parsed-literal::

    2024-09-27 16:05:29.868 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 16:05:40.131 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 16:09:07.754 Hail: INFO: wrote matrix table with 2424689 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.15.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.16.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 16
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.16
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    6069 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.99988.
    6069 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.16.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.16.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.16.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.16.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.16
      --chr 16
      --out SampleData1/Fold_0/hail.train.16
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    6069 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.99988.
    6069 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.16.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 04:09 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.16.vcf
      ref=ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr16.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.16
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [16:60086-19243846]
    Reference markers:              737,912
    Study     markers:                1,487
    
    Burnin  iteration 1:           5 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  9948
    Estimated err:                 9.0e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               22 seconds
    
    Window 2 [16:18062715-58867225]
    Reference markers:              828,530
    Study     markers:                1,630
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  12060
    Estimated err:                 2.1e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    

.. parsed-literal::

    java.util.zip.ZipException: invalid distance too far back
    	at java.base/java.util.zip.InflaterInputStream.read(InflaterInputStream.java:165)
    	at java.base/java.util.zip.GZIPInputStream.read(GZIPInputStream.java:118)
    	at java.base/java.io.FilterInputStream.read(FilterInputStream.java:107)
    	at blbutil.BGZipIt.inflateBlock(BGZipIt.java:314)
    	at blbutil.BGZipIt.lambda$readAndInflateBlocks$3(BGZipIt.java:263)
    	at java.base/java.util.stream.ReferencePipeline$3$1.accept(ReferencePipeline.java:195)
    	at java.base/java.util.ArrayList$ArrayListSpliterator.forEachRemaining(ArrayList.java:1655)
    	at java.base/java.util.stream.AbstractPipeline.copyInto(AbstractPipeline.java:484)
    	at java.base/java.util.stream.AbstractPipeline.wrapAndCopyInto(AbstractPipeline.java:474)
    	at java.base/java.util.stream.Nodes$SizedCollectorTask.compute(Nodes.java:1886)
    	at java.base/java.util.concurrent.CountedCompleter.exec(CountedCompleter.java:746)
    	at java.base/java.util.concurrent.ForkJoinTask.doExec(ForkJoinTask.java:290)
    	at java.base/java.util.concurrent.ForkJoinPool$WorkQueue.topLevelExec(ForkJoinPool.java:1020)
    	at java.base/java.util.concurrent.ForkJoinPool.scan(ForkJoinPool.java:1656)
    	at java.base/java.util.concurrent.ForkJoinPool.runWorker(ForkJoinPool.java:1594)
    	at java.base/java.util.concurrent.ForkJoinWorkerThread.run(ForkJoinWorkerThread.java:183)
    
    Terminating program.
    

.. parsed-literal::

    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.16.vcf ref=ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr16.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.16
    

.. parsed-literal::

    2024-09-27 16:10:30.430 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 16:10:33.858 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 16:11:36.417 Hail: INFO: wrote matrix table with 725178 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.16.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.17.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 17
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.17
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5723 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999884.
    5723 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.17.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.17.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.17.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.17.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.17
      --chr 17
      --out SampleData1/Fold_0/hail.train.17
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5723 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999884.
    5723 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.17.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 04:11 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.17.vcf
      ref=ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr17.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.17
      nthreads=8
    

.. parsed-literal::

    java.lang.IllegalArgumentException: Duplicate marker: 17	1144632	.	C	CT
    	at vcf.Markers.markerSet(Markers.java:129)
    	at vcf.Markers.<init>(Markers.java:83)
    	at vcf.Markers.create(Markers.java:62)
    	at vcf.RefGT.<init>(RefGT.java:90)
    	at vcf.RefTargSlidingWindow$Reader.window(RefTargSlidingWindow.java:326)
    	at vcf.RefTargSlidingWindow$Reader.readWindow(RefTargSlidingWindow.java:303)
    	at vcf.RefTargSlidingWindow$Reader.run(RefTargSlidingWindow.java:246)
    	at java.base/java.lang.Thread.run(Thread.java:829)
    
    Terminating program.
    2024-09-27 16:11:51.161 Hail: INFO: scanning VCF for sortedness...
    

.. parsed-literal::

    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.17.vcf ref=ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr17.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.17
    

.. parsed-literal::

    2024-09-27 16:11:51.378 Hail: INFO: Coerced sorted VCF - no additional import work to do
    2024-09-27 16:11:52.022 Hail: INFO: wrote matrix table with 0 rows and 380 columns in 0 partitions to SampleData1/Fold_0/hail.train.17.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.18.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 18
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.18
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5578 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999906.
    5578 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.18.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.18.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.18.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.18.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.18
      --chr 18
      --out SampleData1/Fold_0/hail.train.18
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5578 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999906.
    5578 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.18.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 04:11 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.18.vcf
      ref=ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr18.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.18
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [18:10083-15185474]
    Reference markers:              478,762
    Study     markers:                1,287
    
    Burnin  iteration 1:           4 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  9156
    Estimated err:                 8.7e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               19 seconds
    
    Window 2 [18:13460439-56841603]
    Reference markers:            1,173,645
    Study     markers:                2,085
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  17072
    Estimated err:                 6.9e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               22 seconds
    
    Window 3 [18:55786342-76632429]
    Reference markers:              654,683
    Study     markers:                1,587
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  8853
    Estimated err:                 2.3e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           1 second
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           1 second
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               11 seconds
    
    Window 4 [18:75886666-78017156]
    Reference markers:               77,128
    Study     markers:                  181
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  92867
    Estimated err:                 3.4e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               1 second
    
    Cumulative Statistics:
    
    Reference markers:            2,267,185
    Study     markers:                4,912
    
    Haplotype phasing time:        43 seconds
    Imputation time:               53 seconds
    Total time:                    2 minutes 5 seconds
    
    End time: 04:13 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.18.vcf ref=ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr18.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.18
    

.. parsed-literal::

    2024-09-27 16:13:57.798 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 16:14:07.440 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 16:17:28.490 Hail: INFO: wrote matrix table with 2267185 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.18.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.19.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 19
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.19
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    4364 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.99987.
    4364 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.19.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.19.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.19.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.19.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.19
      --chr 19
      --out SampleData1/Fold_0/hail.train.19
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    4364 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.99987.
    4364 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.19.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 04:17 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.19.vcf
      ref=ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr19.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.19
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [19:60842-15654475]
    Reference markers:              522,636
    Study     markers:                1,252
    
    Burnin  iteration 1:           4 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  8340
    Estimated err:                 9.4e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               18 seconds
    
    Window 2 [19:14657353-48467440]
    Reference markers:              980,357
    Study     markers:                1,650
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  18607
    Estimated err:                 2.0e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               15 seconds
    
    Window 3 [19:47339481-59118924]
    Reference markers:              395,275
    Study     markers:                1,025
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  9736
    Estimated err:                 1.0e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           1 second
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Cumulative Statistics:
    
    Reference markers:            1,832,506
    Study     markers:                3,788
    
    Haplotype phasing time:        38 seconds
    Imputation time:               41 seconds
    Total time:                    1 minute 45 seconds
    
    End time: 04:19 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.19.vcf ref=ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr19.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.19
    

.. parsed-literal::

    2024-09-27 16:19:13.692 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 16:19:22.149 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 16:21:59.026 Hail: INFO: wrote matrix table with 1832506 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.19.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.20.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 20
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.20
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    4916 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999887.
    4916 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.20.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.20.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.20.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.20.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.20
      --chr 20
      --out SampleData1/Fold_0/hail.train.20
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    4916 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999887.
    4916 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.20.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 04:21 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.20.vcf
      ref=ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr20.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.20
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [20:60343-17038340]
    Reference markers:              507,580
    Study     markers:                1,478
    
    Burnin  iteration 1:           5 seconds
    Burnin  iteration 2:           2 seconds
    Burnin  iteration 3:           2 seconds
    
    Estimated ne:                  9314
    Estimated err:                 8.3e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               19 seconds
    
    Window 2 [20:16012786-51450735]
    Reference markers:              956,124
    Study     markers:                1,813
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  18916
    Estimated err:                 2.1e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               14 seconds
    
    Window 3 [20:50160268-62965354]
    Reference markers:              419,983
    Study     markers:                1,207
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  10473
    Estimated err:                 2.3e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Cumulative Statistics:
    
    Reference markers:            1,812,841
    Study     markers:                4,346
    
    Haplotype phasing time:        40 seconds
    Imputation time:               40 seconds
    Total time:                    1 minute 42 seconds
    
    End time: 04:23 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.20.vcf ref=ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr20.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.20
    

.. parsed-literal::

    2024-09-27 16:23:41.589 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 16:23:49.312 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 16:26:24.612 Hail: INFO: wrote matrix table with 1812841 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.20.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.21.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 21
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.21
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    2811 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999906.
    2811 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.21.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.21.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.21.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.21.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.21
      --chr 21
      --out SampleData1/Fold_0/hail.train.21
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    2811 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999906.
    2811 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.21.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 04:26 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.21.vcf
      ref=ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr21.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.21
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                  380
    
    Window 1 [21:9411239-36416332]
    Reference markers:              719,753
    Study     markers:                1,469
    
    Burnin  iteration 1:           3 seconds
    Burnin  iteration 2:           2 seconds
    
    Estimated ne:                  17989
    Estimated err:                 2.3e-04
    
    Phasing iteration 1:           1 second
    Phasing iteration 2:           1 second
    Phasing iteration 3:           1 second
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           1 second
    Phasing iteration 7:           1 second
    Phasing iteration 8:           1 second
    Phasing iteration 9:           1 second
    Phasing iteration 10:          1 second
    Phasing iteration 11:          1 second
    Phasing iteration 12:          1 second
    
    Imputation time:               21 seconds
    
    Window 2 [21:35704550-48119740]
    Reference markers:              407,769
    Study     markers:                1,061
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  9530
    Estimated err:                 1.1e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Cumulative Statistics:
    
    Reference markers:            1,105,538
    Study     markers:                2,458
    
    Haplotype phasing time:        21 seconds
    Imputation time:               29 seconds
    Total time:                    1 minute 16 seconds
    
    End time: 04:27 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.21.vcf ref=ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr21.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.21
    

.. parsed-literal::

    2024-09-27 16:27:40.991 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 16:27:45.921 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 16:29:24.159 Hail: INFO: wrote matrix table with 1105538 rows and 380 columns in 1 partition to SampleData1/Fold_0/hail.train.21.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/train_data.QC.clumped.pruned.22.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned
      --chr 22
      --make-bed
      --out SampleData1/Fold_0/train_data.QC.clumped.pruned.22
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    2831 out of 172878 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    2831 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/train_data.QC.clumped.pruned.22.bed +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.22.bim +
    SampleData1/Fold_0/train_data.QC.clumped.pruned.22.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.train.22.log.
    Options in effect:
      --bfile SampleData1/Fold_0/train_data.QC.clumped.pruned.22
      --chr 22
      --out SampleData1/Fold_0/hail.train.22
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    2831 variants loaded from .bim file.
    380 people (183 males, 197 females) loaded from .fam.
    380 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 380 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999896.
    2831 variants and 380 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.train.22.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 04:29 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.train.22.vcf
      ref=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr22.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.train.22
      nthreads=8
    

.. parsed-literal::

    java.util.zip.ZipException: invalid distance too far back
    	at java.base/java.util.zip.InflaterInputStream.read(InflaterInputStream.java:165)
    	at java.base/java.util.zip.GZIPInputStream.read(GZIPInputStream.java:118)
    	at java.base/java.io.FilterInputStream.read(FilterInputStream.java:107)
    	at blbutil.BGZipIt.inflateBlock(BGZipIt.java:314)
    	at blbutil.BGZipIt.lambda$readAndInflateBlocks$3(BGZipIt.java:263)
    	at java.base/java.util.stream.ReferencePipeline$3$1.accept(ReferencePipeline.java:195)
    	at java.base/java.util.ArrayList$ArrayListSpliterator.forEachRemaining(ArrayList.java:1655)
    	at java.base/java.util.stream.AbstractPipeline.copyInto(AbstractPipeline.java:484)
    	at java.base/java.util.stream.AbstractPipeline.wrapAndCopyInto(AbstractPipeline.java:474)
    	at java.base/java.util.stream.Nodes$SizedCollectorTask.compute(Nodes.java:1886)
    	at java.base/java.util.concurrent.CountedCompleter.exec(CountedCompleter.java:746)
    	at java.base/java.util.concurrent.ForkJoinTask.doExec(ForkJoinTask.java:290)
    	at java.base/java.util.concurrent.ForkJoinTask.doInvoke(ForkJoinTask.java:408)
    	at java.base/java.util.concurrent.ForkJoinTask.invoke(ForkJoinTask.java:736)
    	at java.base/java.util.stream.Nodes.collect(Nodes.java:333)
    	at java.base/java.util.stream.ReferencePipeline.evaluateToNode(ReferencePipeline.java:109)
    	at java.base/java.util.stream.AbstractPipeline.evaluate(AbstractPipeline.java:545)
    	at java.base/java.util.stream.AbstractPipeline.evaluateToArrayNode(AbstractPipeline.java:260)
    	at java.base/java.util.stream.ReferencePipeline.toArray(ReferencePipeline.java:517)
    	at blbutil.BGZipIt.readAndInflateBlocks(BGZipIt.java:264)
    	at blbutil.BGZipIt.fillBuffer(BGZipIt.java:148)
    	at blbutil.BGZipIt.next(BGZipIt.java:137)
    	at blbutil.BGZipIt.next(BGZipIt.java:53)
    	at blbutil.BlockLineReader.lambda$startFileReadingThread$0(BlockLineReader.java:107)
    	at java.base/java.util.concurrent.Executors$RunnableAdapter.call(Executors.java:515)
    	at java.base/java.util.concurrent.FutureTask.run(FutureTask.java:264)
    	at java.base/java.util.concurrent.ThreadPoolExecutor.runWorker(ThreadPoolExecutor.java:1128)
    	at java.base/java.util.concurrent.ThreadPoolExecutor$Worker.run(ThreadPoolExecutor.java:628)
    	at java.base/java.lang.Thread.run(Thread.java:829)
    
    Terminating program.
    2024-09-27 16:29:36.872 Hail: INFO: scanning VCF for sortedness...
    

.. parsed-literal::

    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.train.22.vcf ref=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr22.GRCh37.map out=SampleData1/Fold_0/beagle.hail.train.22
    

.. parsed-literal::

    2024-09-27 16:29:37.076 Hail: INFO: Coerced sorted VCF - no additional import work to do
    2024-09-27 16:29:37.469 Hail: INFO: wrote matrix table with 0 rows and 380 columns in 0 partitions to SampleData1/Fold_0/hail.train.22.mt
    


.. raw:: html

    <table><thead><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="4"><div style="text-align: left;"></div></td></tr><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="4"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">info</div></td></tr><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">locus</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">alleles</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">rsid</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">qual</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">filters</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">AF</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">DR2</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">IMP</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">END</div></td></tr><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">locus&lt;GRCh37&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">array&lt;str&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">str</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">float64</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">set&lt;str&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">array&lt;float64&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">array&lt;float64&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">bool</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">int32</td></tr>
    </thead><tbody><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10177</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;A&quot;,&quot;AC&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[3.95e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[9.40e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10235</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;T&quot;,&quot;TA&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10352</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;T&quot;,&quot;TA&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[4.20e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[9.40e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10505</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;A&quot;,&quot;T&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10506</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;C&quot;,&quot;G&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    </tbody></table><p style="background: #fdd; padding: 0.4em;">showing top 5 rows</p>
    


.. parsed-literal::

    None
    1 SampleData1/Fold_0/hail.train.1.mt (6468094, 380) (6468094, 380)
    2 SampleData1/Fold_0/hail.train.2.mt (7081600, 380) (13549694, 380)
    3 SampleData1/Fold_0/hail.train.3.mt (5832276, 380) (19381970, 380)
    4 SampleData1/Fold_0/hail.train.4.mt (5732585, 380) (25114555, 380)
    5 SampleData1/Fold_0/hail.train.5.mt (5265763, 380) (30380318, 380)
    6 SampleData1/Fold_0/hail.train.6.mt (5024119, 380) (35404437, 380)
    7 SampleData1/Fold_0/hail.train.7.mt (4716715, 380) (40121152, 380)
    8 SampleData1/Fold_0/hail.train.8.mt (2039669, 380) (42160821, 380)
    9 SampleData1/Fold_0/hail.train.9.mt (3560687, 380) (45721508, 380)
    10 SampleData1/Fold_0/hail.train.10.mt (3992219, 380) (49713727, 380)
    11 SampleData1/Fold_0/hail.train.11.mt (4045628, 380) (53759355, 380)
    12 SampleData1/Fold_0/hail.train.12.mt (0, 380) (53759355, 380)
    13 SampleData1/Fold_0/hail.train.13.mt (2857916, 380) (56617271, 380)
    14 SampleData1/Fold_0/hail.train.14.mt (0, 380) (56617271, 380)
    15 SampleData1/Fold_0/hail.train.15.mt (2424689, 380) (59041960, 380)
    16 SampleData1/Fold_0/hail.train.16.mt (725178, 380) (59767138, 380)
    17 SampleData1/Fold_0/hail.train.17.mt (0, 380) (59767138, 380)
    18 SampleData1/Fold_0/hail.train.18.mt (2267185, 380) (62034323, 380)
    19 SampleData1/Fold_0/hail.train.19.mt (1832506, 380) (63866829, 380)
    20 SampleData1/Fold_0/hail.train.20.mt (1812841, 380) (65679670, 380)
    21 SampleData1/Fold_0/hail.train.21.mt (1105538, 380) (66785208, 380)
    22 SampleData1/Fold_0/hail.train.22.mt (0, 380) (66785208, 380)
    (66785208, 380)
    

.. parsed-literal::

    WARNING: An illegal reflective access operation has occurred
    WARNING: Illegal reflective access by org.apache.spark.util.SizeEstimator$ (file:/data/ascher01/uqmmune1/miniconda3/envs/genetics/lib/python3.10/site-packages/pyspark/jars/spark-core_2.12-3.5.1.jar) to field java.util.concurrent.locks.ReentrantReadWriteLock.readerLock
    WARNING: Please consider reporting this to the maintainers of org.apache.spark.util.SizeEstimator$
    WARNING: Use --illegal-access=warn to enable warnings of further illegal reflective access operations
    WARNING: All illegal access operations will be denied in a future release
    2024-09-27 16:34:45.048 Hail: INFO: wrote matrix table with 66785208 rows and 380 columns in 20 partitions to SampleData1/Fold_0/HAILTRAIN.mt
    


.. raw:: html

    <table><thead><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="4"><div style="text-align: left;"></div></td></tr><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="4"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">info</div></td></tr><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">locus</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">alleles</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">rsid</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">qual</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">filters</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">AF</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">DR2</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">IMP</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">END</div></td></tr><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">locus&lt;GRCh37&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">array&lt;str&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">str</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">float64</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">set&lt;str&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">array&lt;float64&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">array&lt;float64&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">bool</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">int32</td></tr>
    </thead><tbody><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10177</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;A&quot;,&quot;AC&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[3.95e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[9.40e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10235</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;T&quot;,&quot;TA&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10352</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;T&quot;,&quot;TA&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[4.20e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[9.40e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10505</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;A&quot;,&quot;T&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10506</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;C&quot;,&quot;G&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    </tbody></table><p style="background: #fdd; padding: 0.4em;">showing top 5 rows</p>
    


.. parsed-literal::

    None
    

.. parsed-literal::

    2024-09-27 16:34:46.492 Hail: INFO: Reading table to impute column types
    2024-09-27 16:34:48.264 Hail: INFO: Finished type imputation        (0 + 1) / 1]
      Loading field 'chr_name' as type int32 (imputed)
      Loading field 'chr_position' as type int32 (imputed)
      Loading field 'effect_allele' as type str (imputed)
      Loading field 'other_allele' as type str (imputed)
      Loading field 'effect_weight' as type float64 (imputed)
    2024-09-27 16:34:48.331 Hail: WARN: cols(): Resulting column table is sorted by 'col_key'.
        To preserve matrix table column order, first unkey columns with 'key_cols_by()'
    2024-09-27 16:36:27.301 Hail: INFO: Ordering unsorted dataset with network shuffle
    2024-09-27 16:36:29.066 Hail: INFO: Ordering unsorted dataset with network shuffle
    2024-09-27 16:39:14.336 Hail: INFO: Ordering unsorted dataset with network shuffle
    2024-09-27 16:48:09.426 Hail: INFO: Coerced sorted dataset=====>  (19 + 1) / 20]
    2024-09-27 16:48:09.727 Hail: INFO: merging 9 files totalling 10.1K...
    2024-09-27 16:48:09.771 Hail: INFO: while writing:
        SampleData1/Fold_0/HAIL/train_PRS.txt
      merge time: 44.263ms
    2024-09-27 16:48:10.122 Hail: INFO: Reading table to impute column types
    2024-09-27 16:48:10.372 Hail: INFO: Finished type imputation
      Loading field 's' as type str (imputed)
      Loading field 'prs' as type float64 (imputed)
    2024-09-27 16:48:10.574 Hail: INFO: Coerced sorted dataset
    

.. parsed-literal::

             subjectID  PGS002724
    0  HG00097_HG00097     18.058
    1  HG00099_HG00099    -35.684
    2  HG00101_HG00101     23.943
    3  HG00102_HG00102    -20.056
    4  HG00103_HG00103     17.401
             subjectID  PGS002724
    0  HG00097_HG00097     18.058
    1  HG00099_HG00099    -35.684
    2  HG00101_HG00101     23.943
    3  HG00102_HG00102    -20.056
    4  HG00103_HG00103     17.401
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.1.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 1
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.1
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    14013 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999881.
    14013 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.1.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.1.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.1.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.1.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.1
      --chr 1
      --out SampleData1/Fold_0/hail.test.1
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    14013 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999881.
    14013 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.1.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 04:48 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.1.vcf
      ref=ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr1.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.1
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [1:10177-18871145]
    Reference markers:              579,450
    Study     markers:                1,344
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  12600
    Estimated err:                 2.5e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 2 [1:18270328-53836747]
    Reference markers:            1,004,286
    Study     markers:                1,861
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  12483
    Estimated err:                 8.0e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Window 3 [1:49070327-87760656]
    Reference markers:            1,094,455
    Study     markers:                1,963
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  15371
    Estimated err:                 1.2e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           1 second
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Window 4 [1:85578740-144962397]
    Reference markers:            1,037,157
    Study     markers:                1,824
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  16088
    Estimated err:                 8.7e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               18 seconds
    
    Window 5 [1:144224289-175819109]
    Reference markers:              830,857
    Study     markers:                1,537
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  14501
    Estimated err:                 2.0e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               10 seconds
    
    Window 6 [1:173093312-213583529]
    Reference markers:            1,151,258
    Study     markers:                1,960
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  15330
    Estimated err:                 6.6e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Window 7 [1:211080854-240239399]
    Reference markers:              854,129
    Study     markers:                1,793
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  8274
    Estimated err:                 7.2e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               6 seconds
    
    Window 8 [1:239251841-249240543]
    Reference markers:              306,271
    Study     markers:                  716
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  16635
    Estimated err:                 1.2e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               1 second
    
    Cumulative Statistics:
    
    Reference markers:            6,468,094
    Study     markers:               12,391
    
    Haplotype phasing time:        42 seconds
    Imputation time:               1 minute 5 seconds
    Total time:                    3 minutes 49 seconds
    
    End time: 04:51 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.1.vcf ref=ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr1.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.1
    

.. parsed-literal::

    2024-09-27 16:52:00.610 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 16:52:14.012 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 16:54:52.578 Hail: INFO: wrote matrix table with 6468094 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.1.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.2.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 2
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.2
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    13813 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999898.
    13813 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.2.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.2.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.2.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.2.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.2
      --chr 2
      --out SampleData1/Fold_0/hail.test.2
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    13813 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999898.
    13813 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.2.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 04:54 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.2.vcf
      ref=ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr2.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.2
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [2:10179-18145025]
    Reference markers:              573,772
    Study     markers:                1,369
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  7930
    Estimated err:                 2.6e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 2 [2:16417586-49114274]
    Reference markers:            1,053,726
    Study     markers:                1,900
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  16394
    Estimated err:                 8.2e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Window 3 [2:47323299-88186309]
    Reference markers:            1,255,428
    Study     markers:                1,915
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  16727
    Estimated err:                 6.6e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Window 4 [2:86643697-134624574]
    Reference markers:            1,221,295
    Study     markers:                1,730
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  12237
    Estimated err:                 1.4e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               11 seconds
    
    Window 5 [2:133390491-177333395]
    Reference markers:            1,235,397
    Study     markers:                2,049
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  11427
    Estimated err:                 1.2e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           1 second
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 6 [2:175084768-220694589]
    Reference markers:            1,278,696
    Study     markers:                2,058
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  11824
    Estimated err:                 1.3e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 7 [2:218829611-242575162]
    Reference markers:              734,807
    Study     markers:                1,751
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  17532
    Estimated err:                 2.1e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               4 seconds
    
    Window 8 [2:241170836-243188367]
    Reference markers:               65,148
    Study     markers:                  108
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  164880
    Estimated err:                 4.3e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               0 seconds
    
    Cumulative Statistics:
    
    Reference markers:            7,081,600
    Study     markers:               12,251
    
    Haplotype phasing time:        40 seconds
    Imputation time:               53 seconds
    Total time:                    3 minutes 43 seconds
    
    End time: 04:58 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.2.vcf ref=ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr2.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.2
    

.. parsed-literal::

    2024-09-27 16:58:36.213 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 16:58:51.021 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:01:45.909 Hail: INFO: wrote matrix table with 7081600 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.2.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.3.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 3
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.3
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    11785 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999887.
    11785 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.3.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.3.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.3.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.3.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.3
      --chr 3
      --out SampleData1/Fold_0/hail.test.3
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    11785 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999887.
    11785 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.3.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:01 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.3.vcf
      ref=ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr3.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.3
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [3:60069-21643990]
    Reference markers:              718,343
    Study     markers:                1,639
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  9255
    Estimated err:                 2.3e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 2 [3:19789059-60278775]
    Reference markers:            1,194,873
    Study     markers:                1,958
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  14946
    Estimated err:                 7.7e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 3 [3:59462914-105372608]
    Reference markers:            1,325,897
    Study     markers:                1,972
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  13192
    Estimated err:                 1.4e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           1 second
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 4 [3:103398092-143445761]
    Reference markers:            1,153,131
    Study     markers:                2,033
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  15223
    Estimated err:                 1.3e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Window 5 [3:141741447-182238285]
    Reference markers:            1,155,993
    Study     markers:                2,010
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  16136
    Estimated err:                 8.0e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Window 6 [3:179790776-197962381]
    Reference markers:              545,816
    Study     markers:                1,364
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  9843
    Estimated err:                 2.4e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               3 seconds
    
    Cumulative Statistics:
    
    Reference markers:            5,832,276
    Study     markers:               10,457
    
    Haplotype phasing time:        35 seconds
    Imputation time:               40 seconds
    Total time:                    3 minutes 33 seconds
    
    End time: 05:05 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.3.vcf ref=ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr3.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.3
    

.. parsed-literal::

    2024-09-27 17:05:19.815 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 17:05:31.288 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:07:53.762 Hail: INFO: wrote matrix table with 5832276 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.3.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.4.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 4
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.4
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    11041 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999902.
    11041 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.4.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.4.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.4.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.4.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.4
      --chr 4
      --out SampleData1/Fold_0/hail.test.4
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    11041 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999902.
    11041 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.4.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:07 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.4.vcf
      ref=ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr4.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.4
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [4:10005-22520039]
    Reference markers:              786,334
    Study     markers:                1,498
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  12681
    Estimated err:                 2.2e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               5 seconds
    
    Window 2 [4:20544690-60086601]
    Reference markers:            1,127,621
    Study     markers:                1,878
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  14831
    Estimated err:                 7.3e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               4 seconds
    
    Window 3 [4:57715343-105764064]
    Reference markers:            1,443,530
    Study     markers:                2,053
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  14756
    Estimated err:                 1.2e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               5 seconds
    
    Window 4 [4:102709973-150650017]
    Reference markers:            1,376,146
    Study     markers:                2,059
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  13702
    Estimated err:                 1.3e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               5 seconds
    
    Window 5 [4:148086441-182155049]
    Reference markers:            1,012,842
    Study     markers:                1,896
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  13437
    Estimated err:                 7.9e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               4 seconds
    
    Window 6 [4:181312927-191043881]
    Reference markers:              316,589
    Study     markers:                  921
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  9458
    Estimated err:                 8.9e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               2 seconds
    
    Cumulative Statistics:
    
    Reference markers:            5,732,585
    Study     markers:                9,762
    
    Haplotype phasing time:        20 seconds
    Imputation time:               24 seconds
    Total time:                    2 minutes 8 seconds
    
    End time: 05:10 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.4.vcf ref=ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr4.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.4
    

.. parsed-literal::

    2024-09-27 17:10:02.238 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 17:10:13.123 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:12:31.502 Hail: INFO: wrote matrix table with 5732585 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.4.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.5.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 5
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.5
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    10632 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999882.
    10632 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.5.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.5.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.5.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.5.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.5
      --chr 5
      --out SampleData1/Fold_0/hail.test.5
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    10632 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999882.
    10632 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.5.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:12 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.5.vcf
      ref=ALL.chr5.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr5.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.5
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [5:10043-24839896]
    Reference markers:              785,833
    Study     markers:                1,647
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  16474
    Estimated err:                 2.4e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               5 seconds
    
    Window 2 [5:22775078-69306446]
    Reference markers:            1,297,643
    Study     markers:                1,975
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  15101
    Estimated err:                 1.5e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               5 seconds
    
    Window 3 [5:67534999-111982027]
    Reference markers:            1,223,992
    Study     markers:                2,011
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  13613
    Estimated err:                 6.9e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               4 seconds
    
    Window 4 [5:109533719-149269607]
    Reference markers:            1,201,670
    Study     markers:                1,971
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  12773
    Estimated err:                 7.3e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               4 seconds
    
    Window 5 [5:148226432-174766607]
    Reference markers:              787,384
    Study     markers:                1,791
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  14571
    Estimated err:                 2.1e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               3 seconds
    
    Window 6 [5:174167120-180904689]
    Reference markers:              201,966
    Study     markers:                  483
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  10030
    Estimated err:                 2.9e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               1 second
    
    Cumulative Statistics:
    
    Reference markers:            5,265,763
    Study     markers:                9,449
    
    Haplotype phasing time:        18 seconds
    Imputation time:               22 seconds
    Total time:                    1 minute 44 seconds
    
    End time: 05:14 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.5.vcf ref=ALL.chr5.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr5.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.5
    

.. parsed-literal::

    2024-09-27 17:14:15.812 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 17:14:25.845 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:16:34.038 Hail: INFO: wrote matrix table with 5265763 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.5.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.6.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 6
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.6
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    10068 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999902.
    10068 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.6.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.6.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.6.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.6.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.6
      --chr 6
      --out SampleData1/Fold_0/hail.test.6
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    10068 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999902.
    10068 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.6.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:16 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.6.vcf
      ref=ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr6.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.6
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [6:63854-21294216]
    Reference markers:              658,963
    Study     markers:                1,608
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  12514
    Estimated err:                 2.4e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 2 [6:20005714-56895330]
    Reference markers:            1,115,840
    Study     markers:                1,951
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  12421
    Estimated err:                 1.3e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           1 second
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 3 [6:53556519-110284141]
    Reference markers:            1,605,180
    Study     markers:                2,135
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  17115
    Estimated err:                 6.6e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           1 second
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               9 seconds
    
    Window 4 [6:108297729-149276979]
    Reference markers:            1,157,048
    Study     markers:                2,011
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  14130
    Estimated err:                 1.2e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           1 second
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 5 [6:148481212-170715177]
    Reference markers:              694,787
    Study     markers:                1,585
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  16045
    Estimated err:                 8.2e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               4 seconds
    
    Window 6 [6:169508794-171053406]
    Reference markers:               46,342
    Study     markers:                   84
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  168185
    Estimated err:                 7.2e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               0 seconds
    
    Cumulative Statistics:
    
    Reference markers:            5,024,119
    Study     markers:                8,898
    
    Haplotype phasing time:        31 seconds
    Imputation time:               36 seconds
    Total time:                    2 minutes 55 seconds
    
    End time: 05:19 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.6.vcf ref=ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr6.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.6
    

.. parsed-literal::

    2024-09-27 17:19:29.657 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 17:19:39.682 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:21:41.817 Hail: INFO: wrote matrix table with 5024119 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.6.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.7.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 7
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.7
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    9496 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999897.
    9496 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.7.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.7.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.7.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.7.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.7
      --chr 7
      --out SampleData1/Fold_0/hail.test.7
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    9496 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999897.
    9496 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.7.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:21 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.7.vcf
      ref=ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr7.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.7
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [7:14808-22507950]
    Reference markers:              834,932
    Study     markers:                1,538
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  13893
    Estimated err:                 1.9e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 2 [7:21712814-53284120]
    Reference markers:              941,559
    Study     markers:                1,883
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  16574
    Estimated err:                 1.6e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Window 3 [7:51371519-102207161]
    Reference markers:            1,373,374
    Study     markers:                1,942
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  12995
    Estimated err:                 1.1e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               9 seconds
    
    Window 4 [7:98972040-141848611]
    Reference markers:            1,207,792
    Study     markers:                2,091
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  14555
    Estimated err:                 6.7e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 5 [7:140215356-159128567]
    Reference markers:              589,076
    Study     markers:                1,337
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  7130
    Estimated err:                 2.2e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               3 seconds
    
    Cumulative Statistics:
    
    Reference markers:            4,716,715
    Study     markers:                8,420
    
    Haplotype phasing time:        28 seconds
    Imputation time:               35 seconds
    Total time:                    2 minutes 43 seconds
    
    End time: 05:24 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.7.vcf ref=ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr7.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.7
    

.. parsed-literal::

    2024-09-27 17:24:25.990 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 17:24:35.374 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:26:29.775 Hail: INFO: wrote matrix table with 4716715 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.7.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.8.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 8
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.8
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8867 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    8867 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.8.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.8.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.8.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.8.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.8
      --chr 8
      --out SampleData1/Fold_0/hail.test.8
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8867 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    8867 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.8.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:26 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.8.vcf
      ref=ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr8.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.8
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [8:11740-18799855]
    Reference markers:              877,773
    Study     markers:                1,470
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  8676
    Estimated err:                 9.2e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               9 seconds
    
    Window 2 [8:17678413-59729351]
    Reference markers:            1,232,933
    Study     markers:                1,987
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  15278
    Estimated err:                 7.0e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    

.. parsed-literal::

    java.lang.IllegalArgumentException: Duplicate marker: 8	89329148	.	G	A
    	at vcf.Markers.markerSet(Markers.java:129)
    	at vcf.Markers.<init>(Markers.java:83)
    	at vcf.Markers.create(Markers.java:62)
    	at vcf.RefGT.<init>(RefGT.java:90)
    	at vcf.RefTargSlidingWindow$Reader.window(RefTargSlidingWindow.java:326)
    	at vcf.RefTargSlidingWindow$Reader.readWindow(RefTargSlidingWindow.java:303)
    	at vcf.RefTargSlidingWindow$Reader.run(RefTargSlidingWindow.java:246)
    	at java.base/java.lang.Thread.run(Thread.java:829)
    
    Terminating program.
    2024-09-27 17:28:22.413 Hail: INFO: scanning VCF for sortedness...
    

.. parsed-literal::

    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.8.vcf ref=ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr8.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.8
    

.. parsed-literal::

    2024-09-27 17:28:26.744 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:29:17.542 Hail: INFO: wrote matrix table with 2039669 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.8.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.9.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 9
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.9
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    7768 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999892.
    7768 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.9.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.9.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.9.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.9.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.9
      --chr 9
      --out SampleData1/Fold_0/hail.test.9
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    7768 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999892.
    7768 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.9.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:29 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.9.vcf
      ref=ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr9.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.9
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [9:10163-19603557]
    Reference markers:              766,446
    Study     markers:                1,563
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  8930
    Estimated err:                 2.0e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 2 [9:18636347-80220415]
    Reference markers:            1,009,966
    Study     markers:                1,725
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  12947
    Estimated err:                 8.2e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Window 3 [9:79151500-112369214]
    Reference markers:              978,438
    Study     markers:                1,894
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  12438
    Estimated err:                 7.9e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Window 4 [9:110879571-137041122]
    Reference markers:              775,425
    Study     markers:                1,625
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  15593
    Estimated err:                 8.7e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               4 seconds
    
    Window 5 [9:136539237-141148854]
    Reference markers:              166,298
    Study     markers:                  361
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  12562
    Estimated err:                 7.8e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               5 seconds
    
    Cumulative Statistics:
    
    Reference markers:            3,560,687
    Study     markers:                6,879
    
    Haplotype phasing time:        24 seconds
    Imputation time:               31 seconds
    Total time:                    2 minutes 2 seconds
    
    End time: 05:31 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.9.vcf ref=ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr9.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.9
    

.. parsed-literal::

    2024-09-27 17:31:20.202 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 17:31:27.445 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:32:53.929 Hail: INFO: wrote matrix table with 3560687 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.9.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.10.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 10
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.10
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8824 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999901.
    8824 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.10.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.10.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.10.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.10.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.10
      --chr 10
      --out SampleData1/Fold_0/hail.test.10
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8824 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999901.
    8824 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.10.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:32 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.10.vcf
      ref=ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr10.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.10
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [10:60494-18538847]
    Reference markers:              622,290
    Study     markers:                1,540
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  9145
    Estimated err:                 1.9e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Window 2 [10:17068251-54798179]
    Reference markers:            1,000,400
    Study     markers:                1,557
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  11245
    Estimated err:                 1.6e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               6 seconds
    
    Window 3 [10:53792357-90734239]
    Reference markers:            1,092,804
    Study     markers:                1,949
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  13861
    Estimated err:                 1.4e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 4 [10:88784362-123900556]
    Reference markers:            1,013,192
    Study     markers:                1,914
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  16709
    Estimated err:                 7.7e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           1 second
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               6 seconds
    
    Window 5 [10:123093921-135524373]
    Reference markers:              406,013
    Study     markers:                1,114
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  11500
    Estimated err:                 2.8e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               2 seconds
    
    Cumulative Statistics:
    
    Reference markers:            3,992,219
    Study     markers:                7,769
    
    Haplotype phasing time:        26 seconds
    Imputation time:               30 seconds
    Total time:                    2 minutes 14 seconds
    
    End time: 05:35 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.10.vcf ref=ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr10.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.10
    

.. parsed-literal::

    2024-09-27 17:35:08.632 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 17:35:16.600 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:36:57.724 Hail: INFO: wrote matrix table with 3992219 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.10.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.11.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 11
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.11
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8420 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.9999.
    8420 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.11.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.11.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.11.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.11.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.11
      --chr 11
      --out SampleData1/Fold_0/hail.test.11
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8420 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.9999.
    8420 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.11.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:36 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.11.vcf
      ref=ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr11.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.11
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [11:61395-22535739]
    Reference markers:              711,691
    Study     markers:                1,628
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  17448
    Estimated err:                 2.3e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 2 [11:21349458-70803324]
    Reference markers:            1,428,012
    Study     markers:                2,039
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  15187
    Estimated err:                 1.4e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               10 seconds
    
    Window 3 [11:69604836-112290497]
    Reference markers:            1,283,964
    Study     markers:                2,121
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  16129
    Estimated err:                 1.2e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 4 [11:110168424-133043804]
    Reference markers:              692,681
    Study     markers:                1,732
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  19162
    Estimated err:                 2.1e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               3 seconds
    
    Window 5 [11:132350788-134946453]
    Reference markers:               88,218
    Study     markers:                  272
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  20805
    Estimated err:                 1.2e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               0 seconds
    
    Cumulative Statistics:
    
    Reference markers:            4,045,628
    Study     markers:                7,450
    
    Haplotype phasing time:        25 seconds
    Imputation time:               29 seconds
    Total time:                    2 minutes 31 seconds
    
    End time: 05:39 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.11.vcf ref=ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr11.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.11
    

.. parsed-literal::

    2024-09-27 17:39:29.012 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 17:39:37.312 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:41:17.495 Hail: INFO: wrote matrix table with 4045628 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.11.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.12.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 12
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.12
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8198 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    8198 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.12.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.12.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.12.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.12.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.12
      --chr 12
      --out SampleData1/Fold_0/hail.test.12
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    8198 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999891.
    8198 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.12.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:41 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.12.vcf
      ref=ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr12.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.12
      nthreads=8
    

.. parsed-literal::

    java.lang.IllegalArgumentException: Duplicate marker: 12	8400000	.	T	G
    	at vcf.Markers.markerSet(Markers.java:129)
    	at vcf.Markers.<init>(Markers.java:83)
    	at vcf.Markers.create(Markers.java:62)
    	at vcf.RefGT.<init>(RefGT.java:90)
    	at vcf.RefTargSlidingWindow$Reader.window(RefTargSlidingWindow.java:326)
    	at vcf.RefTargSlidingWindow$Reader.readWindow(RefTargSlidingWindow.java:303)
    	at vcf.RefTargSlidingWindow$Reader.run(RefTargSlidingWindow.java:246)
    	at java.base/java.lang.Thread.run(Thread.java:829)
    
    Terminating program.
    2024-09-27 17:41:36.003 Hail: INFO: scanning VCF for sortedness...
    

.. parsed-literal::

    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.12.vcf ref=ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr12.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.12
    

.. parsed-literal::

    2024-09-27 17:41:36.194 Hail: INFO: Coerced sorted VCF - no additional import work to do
    2024-09-27 17:41:36.688 Hail: INFO: wrote matrix table with 0 rows and 95 columns in 0 partitions to SampleData1/Fold_0/hail.test.12.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.13.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 13
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.13
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    6350 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999901.
    6350 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.13.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.13.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.13.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.13.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.13
      --chr 13
      --out SampleData1/Fold_0/hail.test.13
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    6350 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999901.
    6350 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.13.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:41 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.13.vcf
      ref=ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr13.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.13
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [13:19020047-40892639]
    Reference markers:              671,822
    Study     markers:                1,479
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  13287
    Estimated err:                 2.3e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               5 seconds
    
    Window 2 [13:39849087-81726083]
    Reference markers:            1,198,543
    Study     markers:                1,989
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  13447
    Estimated err:                 1.3e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               4 seconds
    
    Window 3 [13:79491658-110702359]
    Reference markers:              937,844
    Study     markers:                1,859
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  18191
    Estimated err:                 7.0e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               3 seconds
    
    Window 4 [13:109580155-115109852]
    Reference markers:              180,096
    Study     markers:                  464
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  11961
    Estimated err:                 1.4e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               1 second
    
    Cumulative Statistics:
    
    Reference markers:            2,857,916
    Study     markers:                5,503
    
    Haplotype phasing time:        12 seconds
    Imputation time:               13 seconds
    Total time:                    1 minute 10 seconds
    
    End time: 05:42 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.13.vcf ref=ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr13.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.13
    

.. parsed-literal::

    2024-09-27 17:42:47.366 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 17:42:52.923 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:44:02.178 Hail: INFO: wrote matrix table with 2857916 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.13.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.14.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 14
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.14
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5742 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999852.
    5742 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.14.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.14.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.14.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.14.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.14
      --chr 14
      --out SampleData1/Fold_0/hail.test.14
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5742 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999852.
    5742 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.14.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:44 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.14.vcf
      ref=ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr14.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.14
      nthreads=8
    

.. parsed-literal::

    java.lang.IllegalArgumentException: Duplicate marker: 14	21649957	.	C	T
    	at vcf.Markers.markerSet(Markers.java:129)
    	at vcf.Markers.<init>(Markers.java:83)
    	at vcf.Markers.create(Markers.java:62)
    	at vcf.RefGT.<init>(RefGT.java:90)
    	at vcf.RefTargSlidingWindow$Reader.window(RefTargSlidingWindow.java:326)
    	at vcf.RefTargSlidingWindow$Reader.readWindow(RefTargSlidingWindow.java:303)
    	at vcf.RefTargSlidingWindow$Reader.run(RefTargSlidingWindow.java:246)
    	at java.base/java.lang.Thread.run(Thread.java:829)
    
    Terminating program.
    2024-09-27 17:44:22.736 Hail: INFO: scanning VCF for sortedness...
    

.. parsed-literal::

    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.14.vcf ref=ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr14.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.14
    

.. parsed-literal::

    2024-09-27 17:44:22.955 Hail: INFO: Coerced sorted VCF - no additional import work to do
    2024-09-27 17:44:23.429 Hail: INFO: wrote matrix table with 0 rows and 95 columns in 0 partitions to SampleData1/Fold_0/hail.test.14.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.15.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 15
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.15
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5569 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999887.
    5569 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.15.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.15.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.15.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.15.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.15
      --chr 15
      --out SampleData1/Fold_0/hail.test.15
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5569 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999887.
    5569 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.15.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:44 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.15.vcf
      ref=ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr15.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.15
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [15:20000041-32890838]
    Reference markers:              325,223
    Study     markers:                  632
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  8658
    Estimated err:                 5.1e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               50 seconds
    
    Window 2 [15:32564155-62934668]
    Reference markers:              905,216
    Study     markers:                1,640
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  7920
    Estimated err:                 2.1e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               6 seconds
    
    Window 3 [15:61592413-93656456]
    Reference markers:              943,933
    Study     markers:                1,788
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  13708
    Estimated err:                 7.3e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               4 seconds
    
    Window 4 [15:92986651-102521131]
    Reference markers:              318,635
    Study     markers:                  973
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  7332
    Estimated err:                 1.1e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               2 seconds
    
    Cumulative Statistics:
    
    Reference markers:            2,424,689
    Study     markers:                4,873
    
    Haplotype phasing time:        17 seconds
    Imputation time:               1 minute 3 seconds
    Total time:                    1 minute 32 seconds
    
    End time: 05:45 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.15.vcf ref=ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr15.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.15
    

.. parsed-literal::

    2024-09-27 17:45:55.571 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 17:46:00.603 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:47:01.284 Hail: INFO: wrote matrix table with 2424689 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.15.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.16.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 16
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.16
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    6069 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.99988.
    6069 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.16.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.16.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.16.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.16.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.16
      --chr 16
      --out SampleData1/Fold_0/hail.test.16
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    6069 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.99988.
    6069 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.16.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:47 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.16.vcf
      ref=ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr16.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.16
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [16:60086-19243846]
    Reference markers:              737,912
    Study     markers:                1,489
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  7912
    Estimated err:                 2.3e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               10 seconds
    
    Window 2 [16:18062715-58867225]
    Reference markers:              828,530
    Study     markers:                1,634
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  9672
    Estimated err:                 2.0e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           1 second
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               6 seconds
    

.. parsed-literal::

    java.util.zip.ZipException: invalid distance too far back
    	at java.base/java.util.zip.InflaterInputStream.read(InflaterInputStream.java:165)
    	at java.base/java.util.zip.GZIPInputStream.read(GZIPInputStream.java:118)
    	at java.base/java.io.FilterInputStream.read(FilterInputStream.java:107)
    	at blbutil.BGZipIt.inflateBlock(BGZipIt.java:314)
    	at blbutil.BGZipIt.lambda$readAndInflateBlocks$3(BGZipIt.java:263)
    	at java.base/java.util.stream.ReferencePipeline$3$1.accept(ReferencePipeline.java:195)
    	at java.base/java.util.ArrayList$ArrayListSpliterator.forEachRemaining(ArrayList.java:1655)
    	at java.base/java.util.stream.AbstractPipeline.copyInto(AbstractPipeline.java:484)
    	at java.base/java.util.stream.AbstractPipeline.wrapAndCopyInto(AbstractPipeline.java:474)
    	at java.base/java.util.stream.Nodes$SizedCollectorTask.compute(Nodes.java:1886)
    	at java.base/java.util.concurrent.CountedCompleter.exec(CountedCompleter.java:746)
    	at java.base/java.util.concurrent.ForkJoinTask.doExec(ForkJoinTask.java:290)
    	at java.base/java.util.concurrent.ForkJoinPool$WorkQueue.helpCC(ForkJoinPool.java:1115)
    	at java.base/java.util.concurrent.ForkJoinPool.externalHelpComplete(ForkJoinPool.java:1957)
    	at java.base/java.util.concurrent.ForkJoinTask.tryExternalHelp(ForkJoinTask.java:378)
    	at java.base/java.util.concurrent.ForkJoinTask.externalAwaitDone(ForkJoinTask.java:323)
    	at java.base/java.util.concurrent.ForkJoinTask.doInvoke(ForkJoinTask.java:412)
    	at java.base/java.util.concurrent.ForkJoinTask.invoke(ForkJoinTask.java:736)
    	at java.base/java.util.stream.Nodes.collect(Nodes.java:333)
    	at java.base/java.util.stream.ReferencePipeline.evaluateToNode(ReferencePipeline.java:109)
    	at java.base/java.util.stream.AbstractPipeline.evaluate(AbstractPipeline.java:545)
    	at java.base/java.util.stream.AbstractPipeline.evaluateToArrayNode(AbstractPipeline.java:260)
    	at java.base/java.util.stream.ReferencePipeline.toArray(ReferencePipeline.java:517)
    	at blbutil.BGZipIt.readAndInflateBlocks(BGZipIt.java:264)
    	at blbutil.BGZipIt.fillBuffer(BGZipIt.java:148)
    	at blbutil.BGZipIt.next(BGZipIt.java:137)
    	at blbutil.BGZipIt.next(BGZipIt.java:53)
    	at blbutil.BlockLineReader.lambda$startFileReadingThread$0(BlockLineReader.java:107)
    	at java.base/java.util.concurrent.Executors$RunnableAdapter.call(Executors.java:515)
    	at java.base/java.util.concurrent.FutureTask.run(FutureTask.java:264)
    	at java.base/java.util.concurrent.ThreadPoolExecutor.runWorker(ThreadPoolExecutor.java:1128)
    	at java.base/java.util.concurrent.ThreadPoolExecutor$Worker.run(ThreadPoolExecutor.java:628)
    	at java.base/java.lang.Thread.run(Thread.java:829)
    
    Terminating program.
    2024-09-27 17:48:18.737 Hail: INFO: scanning VCF for sortedness...
    

.. parsed-literal::

    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.16.vcf ref=ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr16.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.16
    

.. parsed-literal::

    2024-09-27 17:48:22.556 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:49:01.126 Hail: INFO: wrote matrix table with 1522700 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.16.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.17.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 17
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.17
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5723 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999884.
    5723 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.17.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.17.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.17.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.17.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.17
      --chr 17
      --out SampleData1/Fold_0/hail.test.17
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5723 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999884.
    5723 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.17.vcf ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:49 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.17.vcf
      ref=ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr17.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.17
      nthreads=8
    

.. parsed-literal::

    java.lang.IllegalArgumentException: Duplicate marker: 17	1144632	.	C	CT
    	at vcf.Markers.markerSet(Markers.java:129)
    	at vcf.Markers.<init>(Markers.java:83)
    	at vcf.Markers.create(Markers.java:62)
    	at vcf.RefGT.<init>(RefGT.java:90)
    	at vcf.RefTargSlidingWindow$Reader.window(RefTargSlidingWindow.java:326)
    	at vcf.RefTargSlidingWindow$Reader.readWindow(RefTargSlidingWindow.java:303)
    	at vcf.RefTargSlidingWindow$Reader.run(RefTargSlidingWindow.java:246)
    	at java.base/java.lang.Thread.run(Thread.java:829)
    
    Terminating program.
    2024-09-27 17:49:17.082 Hail: INFO: scanning VCF for sortedness...
    

.. parsed-literal::

    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.17.vcf ref=ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr17.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.17
    

.. parsed-literal::

    2024-09-27 17:49:17.277 Hail: INFO: Coerced sorted VCF - no additional import work to do
    2024-09-27 17:49:17.685 Hail: INFO: wrote matrix table with 0 rows and 95 columns in 0 partitions to SampleData1/Fold_0/hail.test.17.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.18.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 18
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.18
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5578 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999906.
    5578 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.18.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.18.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.18.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.18.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.18
      --chr 18
      --out SampleData1/Fold_0/hail.test.18
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    5578 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999906.
    5578 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.18.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:49 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.18.vcf
      ref=ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr18.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.18
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [18:10083-15185474]
    Reference markers:              478,762
    Study     markers:                1,290
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  9011
    Estimated err:                 2.2e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Window 2 [18:13460439-56841603]
    Reference markers:            1,173,645
    Study     markers:                2,084
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  16300
    Estimated err:                 7.5e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           1 second
    Phasing iteration 5:           1 second
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               8 seconds
    
    Window 3 [18:55786342-76632429]
    Reference markers:              654,683
    Study     markers:                1,588
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  16706
    Estimated err:                 8.8e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               3 seconds
    
    Window 4 [18:75886666-78017156]
    Reference markers:               77,128
    Study     markers:                  182
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  84265
    Estimated err:                 5.4e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               0 seconds
    
    Cumulative Statistics:
    
    Reference markers:            2,267,185
    Study     markers:                4,917
    
    Haplotype phasing time:        17 seconds
    Imputation time:               18 seconds
    Total time:                    1 minute 24 seconds
    
    End time: 05:50 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.18.vcf ref=ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr18.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.18
    

.. parsed-literal::

    2024-09-27 17:50:41.984 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 17:50:46.764 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:51:43.016 Hail: INFO: wrote matrix table with 2267185 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.18.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.19.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 19
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.19
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    4364 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.99987.
    4364 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.19.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.19.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.19.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.19.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.19
      --chr 19
      --out SampleData1/Fold_0/hail.test.19
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    4364 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.99987.
    4364 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.19.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:51 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.19.vcf
      ref=ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr19.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.19
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [19:60842-15654475]
    Reference markers:              522,636
    Study     markers:                1,251
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  12625
    Estimated err:                 2.3e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               6 seconds
    
    Window 2 [19:14657353-48467440]
    Reference markers:              980,357
    Study     markers:                1,641
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  16421
    Estimated err:                 2.6e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Window 3 [19:47339481-59118924]
    Reference markers:              395,275
    Study     markers:                1,022
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  8805
    Estimated err:                 2.5e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               3 seconds
    
    Cumulative Statistics:
    
    Reference markers:            1,832,506
    Study     markers:                3,773
    
    Haplotype phasing time:        13 seconds
    Imputation time:               16 seconds
    Total time:                    1 minute 11 seconds
    
    End time: 05:52 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.19.vcf ref=ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr19.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.19
    

.. parsed-literal::

    2024-09-27 17:52:54.163 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 17:52:58.073 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:53:45.037 Hail: INFO: wrote matrix table with 1832506 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.19.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.20.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 20
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.20
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    4916 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999887.
    4916 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.20.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.20.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.20.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.20.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.20
      --chr 20
      --out SampleData1/Fold_0/hail.test.20
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    4916 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999887.
    4916 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.20.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:53 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.20.vcf
      ref=ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr20.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.20
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [20:60343-17038340]
    Reference markers:              507,580
    Study     markers:                1,470
    
    Burnin  iteration 1:           2 seconds
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           1 second
    
    Estimated ne:                  9206
    Estimated err:                 1.9e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Window 2 [20:16012786-51450735]
    Reference markers:              956,124
    Study     markers:                1,807
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    
    Estimated ne:                  15582
    Estimated err:                 7.2e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               7 seconds
    
    Window 3 [20:50160268-62965354]
    Reference markers:              419,983
    Study     markers:                1,208
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    
    Estimated ne:                  15124
    Estimated err:                 2.4e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               2 seconds
    
    Cumulative Statistics:
    
    Reference markers:            1,812,841
    Study     markers:                4,333
    
    Haplotype phasing time:        15 seconds
    Imputation time:               16 seconds
    Total time:                    1 minute 8 seconds
    
    End time: 05:54 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.20.vcf ref=ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr20.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.20
    

.. parsed-literal::

    2024-09-27 17:54:53.622 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 17:54:57.517 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:55:43.685 Hail: INFO: wrote matrix table with 1812841 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.20.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.21.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 21
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.21
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    2811 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999906.
    2811 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.21.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.21.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.21.fam ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.21.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.21
      --chr 21
      --out SampleData1/Fold_0/hail.test.21
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    2811 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 10111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989 done.
    Total genotyping rate is 0.999906.
    2811 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.21.vcf ... 101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:55 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.21.vcf
      ref=ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr21.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.21
      nthreads=8
    
    Reference samples:                2,504
    Study     samples:                   95
    
    Window 1 [21:9411239-36416332]
    Reference markers:              719,753
    Study     markers:                1,463
    
    Burnin  iteration 1:           1 second
    Burnin  iteration 2:           1 second
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  11785
    Estimated err:                 2.4e-04
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               9 seconds
    
    Window 2 [21:35704550-48119740]
    Reference markers:              407,769
    Study     markers:                1,059
    
    Burnin  iteration 1:           0 seconds
    Burnin  iteration 2:           0 seconds
    Burnin  iteration 3:           0 seconds
    
    Estimated ne:                  8340
    Estimated err:                 1.1e-03
    
    Phasing iteration 1:           0 seconds
    Phasing iteration 2:           0 seconds
    Phasing iteration 3:           0 seconds
    Phasing iteration 4:           0 seconds
    Phasing iteration 5:           0 seconds
    Phasing iteration 6:           0 seconds
    Phasing iteration 7:           0 seconds
    Phasing iteration 8:           0 seconds
    Phasing iteration 9:           0 seconds
    Phasing iteration 10:          0 seconds
    Phasing iteration 11:          0 seconds
    Phasing iteration 12:          0 seconds
    
    Imputation time:               3 seconds
    
    Cumulative Statistics:
    
    Reference markers:            1,105,538
    Study     markers:                2,450
    
    Haplotype phasing time:        8 seconds
    Imputation time:               12 seconds
    Total time:                    43 seconds
    
    End time: 05:56 PM AEST on 27 Sep 2024
    beagle.06Aug24.a91.jar finished
    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.21.vcf ref=ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr21.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.21
    

.. parsed-literal::

    2024-09-27 17:56:26.726 Hail: INFO: scanning VCF for sortedness...
    2024-09-27 17:56:29.161 Hail: INFO: Coerced prefix-sorted VCF, requiring additional sorting within data partitions on each query.
    2024-09-27 17:56:56.949 Hail: INFO: wrote matrix table with 1105538 rows and 95 columns in 1 partition to SampleData1/Fold_0/hail.test.21.mt
    

.. parsed-literal::

    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/test_data.clumped.pruned.22.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned
      --chr 22
      --make-bed
      --out SampleData1/Fold_0/test_data.clumped.pruned.22
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    2831 out of 172878 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999896.
    2831 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --make-bed to SampleData1/Fold_0/test_data.clumped.pruned.22.bed +
    SampleData1/Fold_0/test_data.clumped.pruned.22.bim +
    SampleData1/Fold_0/test_data.clumped.pruned.22.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
    PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
    (C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to SampleData1/Fold_0/hail.test.22.log.
    Options in effect:
      --bfile SampleData1/Fold_0/test_data.clumped.pruned.22
      --chr 22
      --out SampleData1/Fold_0/hail.test.22
      --recode vcf
    
    63761 MB RAM detected; reserving 31880 MB for main workspace.
    2831 variants loaded from .bim file.
    95 people (44 males, 51 females) loaded from .fam.
    95 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 95 founders and 0 nonfounders present.
    Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    Total genotyping rate is 0.999896.
    2831 variants and 95 people pass filters and QC.
    Phenotype data is quantitative.
    --recode vcf to SampleData1/Fold_0/hail.test.22.vcf ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
    beagle.06Aug24.a91.jar (version 5.4)
    Copyright (C) 2014-2022 Brian L. Browning
    Enter "java -jar beagle.06Aug24.a91.jar" to list command line argument
    Start time: 05:56 PM AEST on 27 Sep 2024
    
    Command line: java -Xmx51200m -jar beagle.06Aug24.a91.jar
      gt=SampleData1/Fold_0/hail.test.22.vcf
      ref=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      map=plink.chr22.GRCh37.map
      out=SampleData1/Fold_0/beagle.hail.test.22
      nthreads=8
    

.. parsed-literal::

    java.util.zip.ZipException: invalid distance too far back
    	at java.base/java.util.zip.InflaterInputStream.read(InflaterInputStream.java:165)
    	at java.base/java.util.zip.GZIPInputStream.read(GZIPInputStream.java:118)
    	at java.base/java.io.FilterInputStream.read(FilterInputStream.java:107)
    	at blbutil.BGZipIt.inflateBlock(BGZipIt.java:314)
    	at blbutil.BGZipIt.lambda$readAndInflateBlocks$3(BGZipIt.java:263)
    	at java.base/java.util.stream.ReferencePipeline$3$1.accept(ReferencePipeline.java:195)
    	at java.base/java.util.ArrayList$ArrayListSpliterator.forEachRemaining(ArrayList.java:1655)
    	at java.base/java.util.stream.AbstractPipeline.copyInto(AbstractPipeline.java:484)
    	at java.base/java.util.stream.AbstractPipeline.wrapAndCopyInto(AbstractPipeline.java:474)
    	at java.base/java.util.stream.Nodes$SizedCollectorTask.compute(Nodes.java:1886)
    	at java.base/java.util.concurrent.CountedCompleter.exec(CountedCompleter.java:746)
    	at java.base/java.util.concurrent.ForkJoinTask.doExec(ForkJoinTask.java:290)
    	at java.base/java.util.concurrent.ForkJoinTask.doInvoke(ForkJoinTask.java:408)
    	at java.base/java.util.concurrent.ForkJoinTask.invoke(ForkJoinTask.java:736)
    	at java.base/java.util.stream.Nodes.collect(Nodes.java:333)
    	at java.base/java.util.stream.ReferencePipeline.evaluateToNode(ReferencePipeline.java:109)
    	at java.base/java.util.stream.AbstractPipeline.evaluate(AbstractPipeline.java:545)
    	at java.base/java.util.stream.AbstractPipeline.evaluateToArrayNode(AbstractPipeline.java:260)
    	at java.base/java.util.stream.ReferencePipeline.toArray(ReferencePipeline.java:517)
    	at blbutil.BGZipIt.readAndInflateBlocks(BGZipIt.java:264)
    	at blbutil.BGZipIt.fillBuffer(BGZipIt.java:148)
    	at blbutil.BGZipIt.next(BGZipIt.java:137)
    	at blbutil.BGZipIt.next(BGZipIt.java:53)
    	at blbutil.BlockLineReader.lambda$startFileReadingThread$0(BlockLineReader.java:107)
    	at java.base/java.util.concurrent.Executors$RunnableAdapter.call(Executors.java:515)
    	at java.base/java.util.concurrent.FutureTask.run(FutureTask.java:264)
    	at java.base/java.util.concurrent.ThreadPoolExecutor.runWorker(ThreadPoolExecutor.java:1128)
    	at java.base/java.util.concurrent.ThreadPoolExecutor$Worker.run(ThreadPoolExecutor.java:628)
    	at java.base/java.lang.Thread.run(Thread.java:829)
    
    Terminating program.
    2024-09-27 17:57:11.246 Hail: INFO: scanning VCF for sortedness...
    

.. parsed-literal::

    java -Xmx50g -jar beagle gt=SampleData1/Fold_0/hail.test.22.vcf ref=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz map=plink.chr22.GRCh37.map out=SampleData1/Fold_0/beagle.hail.test.22
    

.. parsed-literal::

    2024-09-27 17:57:11.446 Hail: INFO: Coerced sorted VCF - no additional import work to do
    2024-09-27 17:57:12.056 Hail: INFO: wrote matrix table with 0 rows and 95 columns in 0 partitions to SampleData1/Fold_0/hail.test.22.mt
    


.. raw:: html

    <table><thead><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="4"><div style="text-align: left;"></div></td></tr><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="4"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">info</div></td></tr><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">locus</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">alleles</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">rsid</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">qual</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">filters</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">AF</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">DR2</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">IMP</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">END</div></td></tr><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">locus&lt;GRCh37&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">array&lt;str&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">str</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">float64</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">set&lt;str&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">array&lt;float64&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">array&lt;float64&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">bool</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">int32</td></tr>
    </thead><tbody><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10177</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;A&quot;,&quot;AC&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[4.58e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[9.90e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10235</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;T&quot;,&quot;TA&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10352</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;T&quot;,&quot;TA&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[4.69e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[9.90e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10505</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;A&quot;,&quot;T&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10506</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;C&quot;,&quot;G&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    </tbody></table><p style="background: #fdd; padding: 0.4em;">showing top 5 rows</p>
    


.. parsed-literal::

    None
    1 SampleData1/Fold_0/hail.test.1.mt (6468094, 95) (6468094, 95)
    2 SampleData1/Fold_0/hail.test.2.mt (7081600, 95) (13549694, 95)
    3 SampleData1/Fold_0/hail.test.3.mt (5832276, 95) (19381970, 95)
    4 SampleData1/Fold_0/hail.test.4.mt (5732585, 95) (25114555, 95)
    5 SampleData1/Fold_0/hail.test.5.mt (5265763, 95) (30380318, 95)
    6 SampleData1/Fold_0/hail.test.6.mt (5024119, 95) (35404437, 95)
    7 SampleData1/Fold_0/hail.test.7.mt (4716715, 95) (40121152, 95)
    8 SampleData1/Fold_0/hail.test.8.mt (2039669, 95) (42160821, 95)
    9 SampleData1/Fold_0/hail.test.9.mt (3560687, 95) (45721508, 95)
    10 SampleData1/Fold_0/hail.test.10.mt (3992219, 95) (49713727, 95)
    11 SampleData1/Fold_0/hail.test.11.mt (4045628, 95) (53759355, 95)
    12 SampleData1/Fold_0/hail.test.12.mt (0, 95) (53759355, 95)
    13 SampleData1/Fold_0/hail.test.13.mt (2857916, 95) (56617271, 95)
    14 SampleData1/Fold_0/hail.test.14.mt (0, 95) (56617271, 95)
    15 SampleData1/Fold_0/hail.test.15.mt (2424689, 95) (59041960, 95)
    16 SampleData1/Fold_0/hail.test.16.mt (1522700, 95) (60564660, 95)
    17 SampleData1/Fold_0/hail.test.17.mt (0, 95) (60564660, 95)
    18 SampleData1/Fold_0/hail.test.18.mt (2267185, 95) (62831845, 95)
    19 SampleData1/Fold_0/hail.test.19.mt (1832506, 95) (64664351, 95)
    20 SampleData1/Fold_0/hail.test.20.mt (1812841, 95) (66477192, 95)
    21 SampleData1/Fold_0/hail.test.21.mt (1105538, 95) (67582730, 95)
    22 SampleData1/Fold_0/hail.test.22.mt (0, 95) (67582730, 95)
    (67582730, 95)
    

.. parsed-literal::

    2024-09-27 17:59:42.287 Hail: INFO: wrote matrix table with 67582730 rows and 95 columns in 18 partitions to SampleData1/Fold_0/HAILTEST.mt
    


.. raw:: html

    <table><thead><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="4"><div style="text-align: left;"></div></td></tr><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;"></div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="4"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">info</div></td></tr><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">locus</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">alleles</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">rsid</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">qual</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">filters</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">AF</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">DR2</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">IMP</div></td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; " colspan="1"><div style="text-align: left;border-bottom: solid 2px #000; padding-bottom: 5px">END</div></td></tr><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">locus&lt;GRCh37&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">array&lt;str&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">str</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">float64</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">set&lt;str&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">array&lt;float64&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">array&lt;float64&gt;</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">bool</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; text-align: left;">int32</td></tr>
    </thead><tbody><tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10177</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;A&quot;,&quot;AC&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[4.58e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[9.90e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10235</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;T&quot;,&quot;TA&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10352</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;T&quot;,&quot;TA&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[4.69e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[9.90e-01]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10505</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;A&quot;,&quot;T&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    <tr><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">1:10506</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[&quot;C&quot;,&quot;G&quot;]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">-1.00e+01</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">{}</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">[0.00e+00]</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">True</td><td style="white-space: nowrap; max-width: 500px; overflow: hidden; text-overflow: ellipsis; ">NA</td></tr>
    </tbody></table><p style="background: #fdd; padding: 0.4em;">showing top 5 rows</p>
    


.. parsed-literal::

    None
    

.. parsed-literal::

    2024-09-27 17:59:43.085 Hail: INFO: Reading table to impute column types
    2024-09-27 17:59:44.498 Hail: INFO: Finished type imputation        (0 + 1) / 1]
      Loading field 'chr_name' as type int32 (imputed)
      Loading field 'chr_position' as type int32 (imputed)
      Loading field 'effect_allele' as type str (imputed)
      Loading field 'other_allele' as type str (imputed)
      Loading field 'effect_weight' as type float64 (imputed)
    2024-09-27 18:00:29.825 Hail: INFO: Ordering unsorted dataset with network shuffle
    2024-09-27 18:00:31.442 Hail: INFO: Ordering unsorted dataset with network shuffle
    2024-09-27 18:04:23.664 Hail: INFO: Ordering unsorted dataset with network shuffle
    2024-09-27 18:14:10.295 Hail: INFO: Coerced sorted dataset====>   (17 + 1) / 18]
    2024-09-27 18:14:10.414 Hail: INFO: merging 9 files totalling 2.5K...
    2024-09-27 18:14:10.424 Hail: INFO: while writing:
        SampleData1/Fold_0/HAIL/test_PRS.txt
      merge time: 9.174ms
    2024-09-27 18:14:10.864 Hail: INFO: Reading table to impute column types
    2024-09-27 18:14:11.204 Hail: INFO: Finished type imputation
      Loading field 's' as type str (imputed)
      Loading field 'prs' as type float64 (imputed)
    2024-09-27 18:14:11.486 Hail: INFO: Coerced sorted dataset
    

.. parsed-literal::

             subjectID  PGS002724
    0  HG00096_HG00096    -40.344
    1  HG00109_HG00109    -23.098
    2  HG00112_HG00112     2.3007
    3  HG00119_HG00119     -4.378
    4  HG00125_HG00125     9.9591
             subjectID  PGS002724
    0  HG00096_HG00096   -40.3440
    1  HG00109_HG00109   -23.0980
    2  HG00112_HG00112     2.3007
    3  HG00119_HG00119    -4.3780
    4  HG00125_HG00125     9.9591
    Continous Phenotype!
               subjectID  PGS002724
    0    HG00097_HG00097   18.05800
    1    HG00099_HG00099  -35.68400
    2    HG00101_HG00101   23.94300
    3    HG00102_HG00102  -20.05600
    4    HG00103_HG00103   17.40100
    ..               ...        ...
    375  NA20818_NA20818    6.77010
    376  NA20826_NA20826    5.06360
    377  NA20827_NA20827    9.35300
    378  NA20828_NA20828   22.78500
    379  NA20832_NA20832    0.96845
    
    [380 rows x 2 columns]
          SCORE      FID      IID
    0  -40.3440  HG00096  HG00096
    1  -23.0980  HG00109  HG00109
    2    2.3007  HG00112  HG00112
    3   -4.3780  HG00119  HG00119
    4    9.9591  HG00125  HG00125
    ..      ...      ...      ...
    90  17.7130  NA20814  NA20814
    91   1.5582  NA20815  NA20815
    92  19.4650  NA20819  NA20819
    93   4.6435  NA20821  NA20821
    94  17.9070  NA20822  NA20822
    
    [95 rows x 3 columns]
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
    
              PC6     SCORE  
    0    0.032581  18.05800  
    1    0.009107 -35.68400  
    2    0.019296  23.94300  
    3   -0.003620 -20.05600  
    4   -0.030565  17.40100  
    ..        ...       ...  
    375 -0.047403   6.77010  
    376  0.010779   5.06360  
    377  0.003873   9.35300  
    378  0.015832  22.78500  
    379 -0.023086   0.96845  
    
    [380 rows x 10 columns]
      clump_p1 clump_r2 clump_kb p_window_size p_slide_size p_LD_threshold  \
    0        1      0.1      200           200           50           0.25   
    
      numberofpca tempalpha l1weight  Train_pure_prs  Train_null_model  \
    0           6       0.1      0.1     -399.690143          0.227477   
    
       Train_best_model  Test_pure_prs  Test_null_model  Test_best_model  
    0          0.243729    -447.487354          0.14297         0.196619  
    

Repeat the process for each fold.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Change the ``foldnumber`` variable.

.. code:: python

   #foldnumber = sys.argv[1]
   foldnumber = "0"  # Setting 'foldnumber' to "0"

Or uncomment the following line:

.. code:: python

   # foldnumber = sys.argv[1]
   python HAIL.py 0
   python HAIL.py 1
   python HAIL.py 2
   python HAIL.py 3
   python HAIL.py 4

The following files should exist after the execution:

1. ``SampleData1/Fold_0/HAIL/Results.csv``
2. ``SampleData1/Fold_1/HAIL/Results.csv``
3. ``SampleData1/Fold_2/HAIL/Results.csv``
4. ``SampleData1/Fold_3/HAIL/Results.csv``
5. ``SampleData1/Fold_4/HAIL/Results.csv``

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
    Fold_ 2 Yes, the file exists.
    Number of P-values processed:  1
    Fold_ 3 Yes, the file exists.
    Number of P-values processed:  1
    Fold_ 4 Yes, the file exists.
    Number of P-values processed:  1
    

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
    
     


.. parsed-literal::

    We have to ensure when we sum the entries across all Folds, the same rows are merged!
    Fold_ 0 Yes, the file exists.
    Fold_ 1 Yes, the file exists.
    Fold_ 2 Yes, the file exists.
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
    
    Iteration 4:
    Unique rows in current common DataFrame: 1
    Unique rows in next DataFrame: 1
    Common rows after merge: 1
    
    DataFrame 1 with extracted common rows has 1 rows.
    DataFrame 2 with extracted common rows has 1 rows.
    DataFrame 3 with extracted common rows has 1 rows.
    DataFrame 4 with extracted common rows has 1 rows.
    DataFrame 5 with extracted common rows has 1 rows.
       clump_p1  clump_r2  clump_kb  p_window_size  p_slide_size  p_LD_threshold  \
    0       1.0       0.1     200.0          200.0          50.0            0.25   
    
       numberofpca  tempalpha  l1weight  Train_pure_prs  Train_null_model  \
    0          6.0        0.1       0.1     -408.190965          0.233856   
    
       Train_best_model  Test_pure_prs  Test_null_model  Test_best_model  
    0          0.252479    -425.091606         0.137656         0.156992  
    

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
    
    divided_result['pvalue'] = 1
    
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
    | clump_p1         |    1        |
    | clump_r2         |    0.1      |
    | clump_kb         |  200        |
    | p_window_size    |  200        |
    | p_slide_size     |   50        |
    | p_LD_threshold   |    0.25     |
    | numberofpca      |    6        |
    | tempalpha        |    0.1      |
    | l1weight         |    0.1      |
    | Train_pure_prs   | -408.191    |
    | Train_null_model |    0.233856 |
    | Train_best_model |    0.252479 |
    | Test_pure_prs    | -425.092    |
    | Test_null_model  |    0.137656 |
    | Test_best_model  |    0.156992 |
    | pvalue           |    1        |
    


.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAABGUAAAKjCAYAAAC0i7LNAAAAAXNSR0IArs4c6QAAIABJREFUeF7snQmcTeX/xz8z1qzZ1+wSWiUkIvuWyi8lJFtZSiJCSGRNicguRBGttiRb2SKEUtnLlj27aWjm//88OtOd69459965c5eZz/f18ipzn/Ms7+ecY57P/S4RsbGxsZCJgAiIgAiIgAiIgAiIgAiIgAiIgAiIgAgElECERJmA8tZgIiACIiACIiACIiACIiACIiACIiACImAISJTRjSACIiACIiACIiACIiACIiACIiACIiACQSAgUSYI0DWkCIiACIiACIiACIiACIiACIiACIiACEiU0T0gAiIgAiIgAiIgAiIgAiIgAiIgAiIgAkEgIFEmCNA1pAiIgAiIgAiIgAiIgAiIgAiIgAiIgAhIlNE9IAIiIAIiIAIiIAIiIAIiIAIiIAIiIAJBICBRJgjQNaQIiIAIiIAIiIAIiIAIiIAIiIAIiIAISJTRPSACIiACIiACIiACIiACIiACIiACIiACQSAgUSYI0DWkCIiACIiACIiACIiACIiACIiACIiACEiU0T0gAiIgAiIgAiIgAiIgAiIgAiIgAiIgAkEgIFEmCNA1pAiIgAiIgAiIgAiIgAiIgAiIgAiIgAhIlNE9IAIiIAIiIAIiIAIiIAIiIAIiIAIiIAJBICBRJgjQNaQIiIAIiIAIiIAIiIAIiIAIiIAIiIAISJTRPSACIiACIiACIiACIiACIiACIiACIiACQSAgUSYI0DWkCIiACIiACIiACIiACIiACIiACIiACEiU0T0gAiIgAiIgAiIgAiIgAiIgAiIgAiIgAkEgIFEmCNA1pAiIgAiIgAiIgAiIgAiIgAiIgAiIgAhIlNE9IAIiIAIiIAIiIAIiIAIiIAIiIAIiIAJBICBRJgjQNaQIiIAIiIAIiIAIiIAIiIAIiIAIiIAISJTRPSACIiACIiACIiACIiACIiACIiACIiACQSAgUSYI0DWkCIiACIiACIiACIiACIiACIiACIiACEiU0T0gAiIgAiIgAiIgAiIgAiIgAiIgAiIgAkEgIFEmCNA1pAiIgAiIgAiIgAiIgAiIgAiIgAiIgAhIlNE9IAIiIAIiIAIiIAIiIAIiIAIiIAIiIAJBICBRJgjQNaQIiIAIiIAIiIAIiIAIiIAIiIAIiIAISJTRPSACIiACIiACIiACIiACIiACIiACIiACQSAgUSYI0DWkCIiACIiACIiACIiACIiACIiACIiACEiU0T0gAiIgAiIgAiIgAiIgAiIgAiIgAiIgAkEgIFEmCNA1pAiIgAiIgAiIgAiIgAiIgAiIgAiIgAhIlNE9IAIiIAIiIAIiIAIiIAIiIAIiIAIiIAJBICBRJgjQNaQIiIAIiIAIiIAIiIAIiIAIiIAIiIAISJTRPSACIiACIiACIiACIiACIiACIiACIiACQSAgUSYI0DWkCIiACIiACIiACIiACIiACIiACIiACEiU0T0gAiIgAiIgAiIgAiIgAiIgAiIgAiIgAkEgIFEmCNA1pAiIgAiIgAiIgAiIgAiIgAiIgAiIgAhIlNE9IAIiIAIiIAIiIAIiIAIiIAIiIAIiIAJBICBRJgjQNaQIiIAIiIAIiIAIiIAIiIAIiIAIiIAISJTRPSACIiACIiACIiACIiACIiACIiACIiACQSAgUSYI0DWkCIiACIiACIiACIiACIiACIiACIiACEiU0T0gAiIgAiIgAiIgAiIgAiIgAiIgAiIgAkEgIFEmCNA1pAiIgAiIgAiIgAiIgAiIgAiIgAiIgAhIlNE9IAIiIAIiIAIiIAIiIAIiIAIiIAIiIAJBICBRJgjQNaQIiIAIiIAIiIAIiIAIiIAIiIAIiIAISJTRPSACIiACIiACIiACIiACIiACIiACIiACQSAgUSYI0DWkCIiACIiACIiACIiACIiACIiACIiACEiU0T0gAiIgAiIgAiIgAiIgAiIgAiIgAiIgAkEgIFEmCNA1pAiIgAiIgAiIgAiIgAiIgAiIgAiIgAhIlNE9IAIiIAIiIAIiIAIiIAIiIAIiIAIiIAJBICBRJgjQNaQIiIAIiIAIiIAIiIAIiIAIiIAIiIAISJTRPSACIiACIiACIiACIiACIiACIiACIiACQSAgUSYI0DWkCIiACIiACIiACIiACIiACIiACIiACEiU0T0gAiIgAiIgAiIgAiIgAiIgAiIgAiIgAkEgIFEmCNA1pAiIgAiIgAiIgAiIgAiIgAiIgAiIgAhIlNE9IAIiIAIiIAIiIAIiIAIiIAIiIAIiIAJBICBRJgjQNaQIiIAIiIAIiIAIiIAIiIAIiIAIiIAISJTRPSACIiACIiACIiACIiACIiACIiACIiAC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LVHG8TNWnsqVK5cJ32FVLYpqzkZRhVW3KGCQEz2LyJKCG+dFLxgaw8TYhiISRRyKHnfeeacJX2NbR6NQN3DgQBNqRD4UxijyJGSJFWXYN/PycFzueVRUFEqWLGlCn1588UWz/47GHDU9evQwa7I8lhiexyTQFHCceXu6zxzDWZTx973ozX2rtiIgAiIgAiKQVAQkyiQVWfUrAiIgAiIgAiIgAiIgAiIgAiIgAiIgAgkQkCij20MEREAEREAEREAEREAEREAEREAEREAEgkBAokwQoGtIERABERABERABERABERABERABERABEZAoo3tABERABERABEQg2RHY/dln2DFpEo5v3Yorp07hiVWrUMipMlJSLPqzrZ9h0reTsPXgVpy6eAqreqxC9VLxKzIlxbjqUwREQAREQAREIDwJSJQJz33TrEVABERABERABBIgsHPWLJzbvx9ZihTB0tatAybKzNowC/tP7keRnEXQenpriTK6S0VABERABERABBIkIFFGN4gIiIAIiIAIiEDQCcytXh3ZSpZE2ixZsHPGDMRcu4YyLVvioXfeQaq0aX2e36VjxzAhXz63okz1kdVRMk9JZEmfBTPWz8C1mGtoWbEl3nnyHaRN7fu4x84dQ74e+STK+LxzulAEREAEREAEUgYBiTIpY5+1ShEQAREQAREIaQIUZU5s3YrSzZujXNeu+GvPHnzdrh3uaN8eDw4bhu+HDsXGoUMTXEPtSZNQpkWLeG08EWUYatS8QnN0rdUVe47vQbuZ7dC+ansMazIMQxcPxdCvEh53UstJaFEp/rgSZUL6dtPkREAEREAERCBkCEiUCZmt0EREQAREQAREIOUSoChz4dAhtN+zBxGRkQbEj+PH49uePdHl7FlEX7iAqDNnEgSUMU8epM2c2WtR5tBfh7Bn8B5E/jvu+FXj0fOTnjg75iwu/H0BZy4lPG6eLHmQOX38cSXKpNx7WSsXAREQAREQAW8ISJTxhpbaIioqCps3b0bevHmROnVqEREBERABERABvxBY3awZbsqTBxXHjInr7+wvv2B5w4aou2IFMhcr5tM4USdPYlGFCnhwzhzkrlTphj6azWqGPJnzYMyj/437y/Ff0HBqQ6zouALFcvg27smLJ1FhTAXMaTkHlQrfOK5Pi9FFIiACIiACIuBE4Nq1azh27BjKly+P9OnTi08Y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3O6MuWRSpUoV9HlqAuFPwO69ovNZ+O+xRJnw38OArkAPfUBxazAREAEREAEREIEwImB3eAqjpQRtqtbvmglNgIJIsI15b1iO2p0xCTJz7chEILEE7N4rOp8llnDwr5coE/w9CKsZ6KEPq+3SZEVABERABERABAJIwO7wFMCphO1Q0dHRJkltQla+fPmgr49Jcw8cOOB2HpkzZ0apUqWCPk9NIPwJ2L1XdD4L/z2WKBP+exjQFeihDyhuDSYCIiACIiACIhBGBOwOT2G0FE1VBEQgRAjYvVd0PguRjUrENCTKJAJeSrxUD31K3HWtWQREQAREQAREwBMCdocnT/pQGxEQARFwJGD3XtH5LPzvF4ky4b+HAV2BHvqA4tZgIiACIiACIiACYUTA7vAURkvRVEVABEKEgN17ReezENmoRExDokwi4KXES/XQp8Rd15pFQAREQAREQAQ8IWB3ePKkD7URAREQAUcCdu8Vnc/C/36RKBP+exjQFeihDyhuDSYCIiACIiACIhBGBOwOT2G0FE1VBEQgRAjYvVd0PguRjUrENCTKJAJeSrxUD31K3HWtWQREQAREQAREwBMCdocnT/pQGxEQARFwJGD3XtH5LPzvF4ky4b+HAV2BHvqA4tZgIiACIiACIiACYUTA7vDk61KioqJw4cIFsMxy+vTpfe1G14mACIQhAbv3is5nYbipTlOWKBP+exjQFeihDyhuDSYCIiACIiACIhBGBOwOT94s5fTp05g+fTomTpiAffv3x11avFgxdOzUCW3atEGOHDm86VJtRUAEwpCA3XtF57Mw3FSJMuG/acFcgR76YNLX2CIgAiIgAiIgAqFMwO7w5MncY2JiMHDgQIwYPhx/R0cjX0QEbo2NBf1jogDsiojAsdhYpEubFr379MGAAQMQERHhSddet3n99dfNXCxLnTo1brnlFjRu3NiMmy1bNq/7tLuAY9aoUQMPPvigXVNUr14d3377rWlHBvnz58d9992HwYMHo2zZsrbXe9IgNjYWr7zyCj788EMcO3bMrP2LL77w5FK1+ZfAG2+8geXLl5u9ctwzd4AKFy4M61nyBWKRIkVQq1YtTJ061avLX331VWzduhVLly716rqkbmz3XtH5LKl3IOn7l6dM0jNOViPooU9W26nFiIAIiIAIiIAI+JGA3eHJbigKMu3atcOMGTNwK4Ca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NWrUyOpCYUKZ//znP+5+gy049s9//rPdddddWV2vpJ2C9j799NN2wgknuDBBA/bN3cI4bO614nD8k08+aWeffbYtXbrUtthii02apGdYgU1pM2MKCgpMfzT7Jdeb+tULL7xgH374YamXyvR7hfFZriuV+/MTyuTeOFFXoNMnqpzcDAIIIIAAAgiUo0CmwVO6S2kmimakWE/LboZM8RMtMrOp5ma0DB06tNzuqHv37vbJJ5/Y119/nfUg9d1337W///3v9sYbb9hvv/1mBx54oA0fPnyjmSTbb7+9devWzf3s+uuvt2+//dbatGlj99577yYzTv75z3+6v1cYUr9+fTvxxBPthhtuKApdpk+fbl26dLHx48fbmDFjbNKkSdayZUv32o9CpTvuuMPef/99W7lype2444524YUXulkUxYOVbF5fKh7K6Bx6jWn//fd319WmYOWWW24xzXzacsst3YLMt956q9WrV8/9PPhvab2+opkx2l+zJjQzJvX1KO2rWqqmy5cvt8suu8yee+4598+aDXLppZfaSSedVHQf2k+hg/x1j2+//badeeaZdtRRRzmfiRMnulk+WjxagdXFF1/szik3hQNq1x//+Ed76KGHNgqEVEvdm/xr1apl++23n7u/1NflNNNDNlqUWrOmdF+qsdpz9NFHb/Q8vvTSS67maqde29F5NONIz4M2LXD9j3/8wx577DE3G0kLW5911ll2+eWXuxkppW16XlUPtaOkraRQJngW99lnH7v99tudg+x0Xd37v//9b/cqWdOmTU3n18wwPYfBVvz1Jc3GOfXUU905tK+eQe2v1/yuuuqqjV5V+vjjj90rcO+8845zTbdl+r3C+KzcfuVV2IkIZSqM3s8L0+n9rButRgABBBBAAIHcC2QaPJXUAg1CW2zXwpYWLLXCvmVf2yJvQp41qdLEvvryq3L5KpMCldq1a9txxx1nCiyy2TS47NSpkwtbzjvvPKtSpYrdfPPN9tlnn5kGoBroatNAWLMVWrRoYZdccombmaBXdjTo/u9//1sUAOnvtO6G9tHgV68KKUDo2bOnPfXUU+5cQSijc2uNG71mtWHDBrfPPffcY5otoldL9KWq119/3QUFGnyrfdo2Z6bMzz//bA0aNHADboUsWqD5ggsusHPOOceFMd99950NGTLE3e/MmTPdgDy43jbbbOPCEq1jozbqNaUnnnjChTDPPvus6ed6VUz/K1OFPAozFCz961//cqGFrqngRZuOU7jRqlUrd2/77ruvC1F++eUXdx21QYGWHBXuKOiSr9ZCkak2BSraTwFLsCnA6tq1qzVv3tx0v3r9R46qk9qmTaGMwh39u86pfUeMGGHTpk1z+6nN2hQKHX/88W5mldZ70etZemYUeARBmZ43vRqnNimoeOutt+yaa66xiy66yIUc6ba1a9faVltt5a4b1Lb4vulCGYViavPgwYNdmxTQ6BW9Bx54wDp37uxq/OWXX7owUD9Xm4ItXSij4Ozkk0+2Aw44wKZMmeKCyYcffnijV5XUBxTWKWDTc5Juy/R7hfFZNr+d4r0PoUy86xO71tHpY1cSGoQAAggggAACMRHINHgqqZnPPPOMCz7sIDPbaTNu5HMzm2Gm8x1zzDGbcaLfD12yZIkbLGs2hWYuBJsGkgo9gk1BQ7BQqQbvWghYA+3g9Y8VK1a4mRcadN90003uMA3816xZY3r9ScGPNgURavebb75pbdu2tYULF7rBvEKdv/71r0XX0ywIBRkKKXbfffeiUEazRhRUpNuC11I0i0Shg2ZzaAsbyihk0KZZHJptojBCIYIW+9XAXgPx1FeZFGAopJo8ebJ7FSy4ngb7M2bM2Ki5CmV0H7p3GWnTDJe+ffu6WUDuOfnfduihh7qgS7M4FGYFs61GjRrlQqJgC0Kr1FesVD+FWApZVIMgLNOsJIVKmrmk2hffdJxCRP1Msz50/9qCRZCD2unv9Eqe9tOzo9BBz812223nQiPdd0mzXl577TU76KCDNnmGFUYpcFK7ghlHxdumRajbt2/vwi95l7SlC2V++uknF7psvfXWaZ8fhZR6rnUNzfLRzC5t6UIZPeu672DTjBiFbHoOUjfdr8KkCRMmpL12pt8rjM/S0nnzA0IZb0oVj4bS6eNRB1qBAAIIIIAAAvETyDR4KqnFem1l4osTreCE/332uqy3pc9lP5VvfXv2dTMhNndLF8oEg/zg/IcccoibbaEZGZpFoMGzZkukbnqFRjMP9EqTNgUOmgUxbty4ot30tZxdd93VzabQjBfNUtBrK4sWLbImTZoU7acBtF5R0YwN/TxoT3Bc6nUXL17swgOFJhrQB2GSXt/RzAptYUKZ4q8XaQaFBvqaGaNZIZoBohBGoVLqptdXFIpopkVwPb0CpBlAqVtJoUwwW0i+qV/pCV6TkZtmZQShjMKQhg0bFp028NGgv0+fPkV/rxBJM0QUNASbXi1S2KNXb/Qqkzb9ne5Rr7HJPti0dst9993n/lWhjMIj1Th108wZBW0KexRmacaS2n3KKaeU+HhqdoxeN9N5NMsq2PT6mUxVx3TrF+mZ16tSmpWl56ikLV0oo69caW2X1E0hkmY+6TlTcKUQMdg0U0lr/2hLF8oo9EtdpFmLZCtEk2PqJh8FQql1KN72TL9XGJ+VWG6v/pJQxqtyVXxj6fQVXwNagAACCCCAAALxFMg0eCqp1Rpsvv3l22Z9y+Gexpu1bdXWzTbZ3C3d60tam0VBgDYNzDVzQaGMZmxoJkC6TcFBcFywjseDDz5YtHvxcESzI7SmR7pNAYfWGQlCB82y0Gs+waaZMVrrRYu+6tUQDdQVGumVHwU+wWKvYUIZfX1Jg3/N8lAwpJkxQVASzOBJ114FSBrgB9fTrJ7UNWF0XEmhjF7zURihdX1St2ARZs0Q0SsywetLqbOYtH/gU3wGSUmfci6+bzAzREGIZjpp5kv16tVduKPwRgGLtmBNmdRFkPX3qXUOZgwpvArWjylupVexUp+J4j/XQr4KN0raFJTo9SwFKJqNU9KWLpRRsKgZRqmbwiG9MqU/un8FcHodTa+lpb6GlC6UUb1S+0M6owEDBrj1jzTzK92W6fcK47PN/W1X8ccTylR8DbxqAZ3eq3LRWAQQQAABBBCIUCDT4Kmkpuy62642d8Vcs8PLoaGTzXapu4v997PfX7HZ3E2DZw0W0y30mzogXb16tdWtW9e90qJ1Q4pvNWvWNM1I0JZNKKM1T7Q2iIICLZhbfFMgopkY6UKHzz//3M0gSZ3VoHMoXNBitmUJZUpa6Ddol9YNOfzww931gjVUUtvcqFEjtwBxaSFQupkyeh1KMzVKmimjdXZ22mknF8oodFCYlrptTiijmSta00czZILX0RT6aK0ahSNhQplgJlRpM2W0potm3yj8KmnTa3Cpi+ym7qPgSusI6fPkWk+npK20hX6Lh0EdO3Z0IZ7Cr2BTsKXXzsozlNEzo75T/FW21PZn+r3C+Gxzf9NV/PGEMhVfA69aQKf3qlw0FgEEEEAAAQQiFMg0eCqpKXGdKaO2BkFDuk9iF58loAGrBuypA9mS7jmbUEazLnbZZRc3eyTd7AidO13ooM8Ma8FWvSIVfAFI69toFoXWvSnvUEbnVlCkWTmawZNuCxvKBGvKFF8rSOvT6B5T15Qp71DmL3/5iwtefvjhh6LXiYIZQXoFKUwoI2/VXX/SBRCvvvqqadaKXnPT2i1htsA1td7Fjw8TyijY0ULUqWu9BDN5yjOU0StdmuGlGVylPTP6WbDOUPH9GJ+FeVLiuS+hTDzrEttW0eljWxoahgACCCCAAAIVLFCWUCaua8oElFqTRevE6KtGej1Er2ToFSa92qKvGGmh0mANG/2d/l1f+9GgXa/4aG0avbqi2RxaVyUYXGoWTmmvL2k/zZzQuh56hUSLt2qdERlr/Q/N4NDMk3ShjNZL0YyVLbbYwi0wrBkeWnRWsz70ikt5hzJqr772pLZqjRm98qKASgsCa10W3bvWcQkbyqjdGrTPmTOn6OtL+vKUAhG9hqXXm7TlYqZMEMrp1Rstrqx1Um677TY3a0eLD4cJZdRGBUv9+vVzM1o0Y0kzqzSzRbOIgsWJNctK4YzW29G6Q1pYeP78+W5BZYV9qWvNFO/uejb1SljqwtSp+4QJZbTAtdb90RefFO4pHNOnwfXslFcoo36kRX71Kpt8CWUq+Bd4BV6eUKYC8X28NKGMj1WjzQgggAACCCAQhUBZQpm4fn0p1UuhghZr1folCjU0mN57773d14A0mNbCucGmxUwVECgs0WsZWoekXbt27rWmYAHcbGbKBOfTq0YKO/QaldYz0bEa1F955ZXuy03pQhkdrwG/whCt2aGBv/5Zn5/WsbkIZXRNzazQYF6L02rTbAsFUJpBo9etwoYyOsfy5ctdQKXwS19M0mtZWgBYX3oKtlyEMjq3PjGt8E1r82j2iP5Z4ZxCsrChjM6noEfhiHz0mfI99tjDhX6aIaNNIZSuqeBDs6X06prCtV69erlFm0v6alNgoGdM6xsFX9Yq3ufDhDIKnhSw6atgCvi0ho5qoHWKyiuUefrpp+3UU091i1ArnCGUieK3dDyvQSgTz7rEtlWEMrEtDQ1DAAEEEEAAgQoWKEsoo5kALbZrYUsLllph38Iy30HehDxrUqWJffXlV26wy4ZAZRPQujV6HUivdaV++SiuDlo0WAtlB+EWoUxcK5X7dhHK5N44UVcglElUObkZBBBAAAEEEChHgbKEMrp88Dlj62lm6T9glL6li8zs/61HqvMMHTq0HO+IUyHgl4BeX9JsGr0SFOdNn+7WzCPNAPvDH/5QalMz/V5hfBbnSmfXNkKZ7JzY638CdHoeBQQQQAABBBBAoGSBTIOndG7ff/+9tdmvjS1eutgKDi8waxBCeJlZ/uR8a9a4mb337u/rc7AhUFkF9KWwRx991LQmTPDFqDhaaN0cvbakT2Jn2jL9XmF8lkkw/j8nlIl/jWLVQjp9rMpBYxBAAAEEEEAgRgKZBk+lNVX/j3mHjh1s1dpVVtClILsZM4vM8v+db7Vr1rY3Xn/Drc/BhgACyRLI9HuF8Zn/9SaU8b+Gkd4BnT5Sbi6GAAIIIIAAAh4JZBo8ZboVBTM9D+tpi75eZHlN8qxw10KzVmZWNeXI38xsoVneZ3lWuLTQtm2xrU2dMpVAJhMuP0fAU4FMv1cYn3la2JRmE8r4X8NI74BOHyk3F0MAAQQQQAABjwQyDZ6yuRW9yqSvDd1737225Lslll8z3wrqFJhp7d5fzfJX5lvB2gJruk1TO+fsc+y8887jlaVsYNkHAU8FMv1eYXzmaWEJZfwvXEXdAZ2+ouS5LgIIIIAAAgjEXSDT4ClM+/VVpvHjx9vo0aNt8eLF9vOKn22rultZs2bNbODAgda3b1++shQGlH0R8FQg0+8VxmeeFpZQxv/CVdQd0OkrSp7rIoAAAggggEDcBTINnuLeftqHAALxE8j0e4XxWfxqFrZFvL4UVqyS70+nr+QPALePAAIIIIAAAmkFMg2ewtBppszzzz9fNFNm1coVVrtO3aKZMkcccQQzZcKAsi8Cngpk+r3C+MzTwqY0m1DG/xpGegd0+ki5uRgCCCCAAAIIeCSQafCUza1oTZm7777b7h95n3373RKrVyvfdm5YYHVqmK1cZzbvh3z7aU2BbdO0iZ119jl2/vnns6ZMNrDsg4CnApl+rzA+87SwhDL+F66i7oBOX1HyXBcBBBBAAAEE4i6QafCUqf36+tJhPQ+1rxd9Y+23z7Pz2hfasXuZ1dQiv//b1v5qNvYjs3veyLPZXxbadi22tcl8fSkTLT9HwFuBTL9XGJ95W9qihjNTxv8aRnoHdPpIubkYAggggAACCHgkkGnwVNqtKJA5sGMHK1y/ysaeVGDdd85849PmmR33eL7lVa9t/3n9jZx9FvuVV16xO++802bPnm0//vij1a1b1/bee2879thj7bTTTrOaNWtmbmxM9wj+2/bxxx93Cyin2/70pz/Zo48+WvTjxo0b21577WXDhg2zDh06lNvd3Xbbbc76m2++sT333NM++OCDcjt3ZTiR6nj11Vfb3Llz7YwzztioZunuv7CwcLNoVCO9anjxxRe7vhFs//3vf61Nmzb22WefWcuWLct8jUy/VxiflZk2NgcSysSmFH40hE7vR51oJQIIIIAAAghEL5Bp8JSuRXpl6Y/7tbGff1hsr51bYHs1y77tHy0263xvvm3VsJm98+575f4qk0IHDXJ79+5t/fv3t2233daWLVtmL774ohvw3njjjTZo0KDsGxyzPcOEMgqnxo4d6+7gq6++smuvvdbmzZtn7733XrkEYh999JELu/72t7/ZkUceaXXq1LHWrVvHTCy+zVm/fr3tvPPOduWVV9rpp59u8+fPN/WtYHvkkUds5MiRNmvWrI1uol27dpt1Uzrvqaeeal9//bXrH6lbv379rFatWqZ9yrpl+r3C+KyssvE5jlAmPrXwoiV0ei/KRCMRQAABBBBAoAIEMg2e0jVJoYfCj5fOsqxmyBQ/z0tzzQ59wFx4MnTo0HK786lTp9phhx1mV1xxhV0wH6DOAAAgAElEQVR33XWbnHfBggX2xRdfWI8ePcrtmuVxonXr1lmNGjWyOlWYUOY///mPu99gC47985//bHfddVdW1ytpp6C9Tz/9tJ1wwgkuTNhhhx3KfL7gwDAOm32xGJzgySeftLPPPtuWLl1qW2yxxSYt0jOswGZzZ8YUP3Fpocy0adOsV69etmjRItPsqrJsmX6vMD4ri2q8jiGUiVc9Yt8aOn3sS0QDEUAAAQQQQKCCBDINnkpqlr6y1HK7FtayxlKbdUHZX6Nod1eefb2+if2fL78qt68yde/e3T755BM3A6Bq1apZqb777rv297//3d544w377bff7MADD7Thw4dvNJNk++23t27durmfXX/99fbtt9+61zzuvffeTWac/POf/3R/rzCkfv36duKJJ9oNN9xQFLpMnz7dunTpYuPHj7cxY8bYpEmT3KsieqVEodIdd9xh77//vq1cudJ23HFHu/DCC90siuLBSjavLxUPZXQODbT3339/d11tClZuueUW0+toW265pekrWbfeeqvVq1fP/Tz4b+n777/fNDNG+2uGh2bGpL4epX0VsCloW758uV12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MgU6Wpx4++23z7pt7IgAAggggAACCCCAAAIIIFC+AoQy5etZEWcjlKkI9ZRrllco8+abb1rXrl3tyCOPtNGjR6e9q2z3C06gsEgLEBffCGUq+MHh8ggggAACCCCAAAIIIFDpBQhl/H8ECGUquIbl8fqSPonduXNn22effWzy5Mnus9olbdnul3osM2Uq+AHh8ggggAACCCCAAAIIIIBAGgFCGf8fDUKZMtbwxRdfdOu9LFu2zL1+pM9Tz54927RAr4KWbLeDDz7YhSgvvfTSRofoPN27d7f777+/1FOpE+pLTc2aNbN///vfVrt27RL3z3a/TO2m02cS4ucIIIAAAggggAACCCCAQDQCjM+icc7lVQhlQupq5kifPn3s9ddftzp16tiqVavs7bfftjZt2rh1XBo2bGh33HFH1mfVl5L0JSS9DqRgRZvCnfbt29uECRPctdJtS5YssQMPPNB9vUkBka5d0pbtftk0mk6fjRL7IIAAAggggAACCCCAAAK5F2B8lnvjXF+BUCak8BlnnGFTp061sWPH2v777+9mubzzzjsulHnkkUfslltuCbU4r76w1Lp1a2vcuLFb7HfNmjU2ePBga9q0qQt+8vLyXAuvueYa92f+/PnWsmVL++WXX6xDhw42d+5ce/TRR91MndRt3333tRo1amS9X7YMdPpspdgPAQQQQAABBBBAAAEEEMitAOOz3PpGcXZCmZDKjRo1suHDh7uvHG3YsMGqVatWFMro9aEjjjjCfco6zKagZdCgQTZ9+nR3Ps2OGTFihDVo0KDoNMGCu8ECu0HnS3edsPtl2146fbZS7IcAAggggAACCCCAAAII5FaA8VlufaM4O6FMSOVatWrZ+PHj3XovxUOZSZMmWf/+/UOHMiGbUKG70+krlJ+LI4AAAggggAACCCCAAAJFAozP/H8YCGVC1lCvLLVr187uuuuuTUIZzXb54IMPbMaMGSHP6s/udHp/akVLEUAAAQQQQAABBBBAINkCjM/8ry+hTMgaPv3003biiSfaBRdcYCeccIJb10V/t2DBAhs6dKg988wz1rt375Bn9Wd3Or0/taKlCCCAAAIIIIAAAgggkGwBxmf+15dQpgw1vPvuu90Xk1auXGmFhYXuDPoU9U033WTnnntuGc7ozyF0en9qRUsRQAABBBBAAAEEEEAg2QKMz/yvL6FMGWu4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+/HHH03rzNSqVWszzuTPoXR6f2pFSxFAAAEEEEAAAQQQQCDZAozP/K8voUyWNRw1apTdeeed7nPXzZs3d68uXXnllW5Nmcq00ekrU7W5VwQQQAABBBBAAAEEEIizAOOzOFcnu7YRymTh9Pjjj9spp5xiO+64o1vQd+HChW5x3wsuuMBuv/32LM5Q+i4636BBg+zVV1+16tWrW9++fd2rUPXr1097oF6b0j6TJk2yOXPm2Lp169w6N3qt6uijj97kuFdeecUuu+wy++STT6xRo0Z25plnuvVwqlSpEqr9dPpQXOyMAAIIIIAAAggggAACCORMgPFZzmgjOzGhTBbU+++/v7Vs2dKefvrpohDjmmuusRtuuMFWrVplVatWzeIsJe+ycuVKa926tTVs2NCGDRtmq1evtsGDB7vZOFqzJi8vr8QDdd0WLVq4sKhbt25Wo0YNGzNmjD344IN2zz33uE9zB5vWwNGaN8ccc4ydccYZLpjRNS6++GK7/vrrQ7WdTh+Ki50RQAABBBBAAAEEEEAAgZwJMD7LGW1kJyaUyYJaa8YokDn00EOL9v7hhx+scePG9vnnn9sf/vCHLM5S8i633Xabm7ESvBalvd544w0XokycONF69+5d4oGaKbNixQqrV6/eRj/v2bOnzZs3z50v2PSp7i+//NI++OADtyixNgVKCpa++eYba9CgQdbtp9NnTcWOCCCAAAIIIIAAAggggEBOBRif5ZQ3kpMTymTBrCBDn7pu27Zt0d4KRbSejGahtGnTJouzlLxLly5d3EybadOmbbRDq1atrEePHjZy5MhQ577iiivs1ltvda8zaVu/fr3VrVvXrX+jnwWbQprtt9/ennzySevfv3/W16DTZ03FjggggAACCCCAAAIIIIBATgUYn+WUN5KTE8pkwaxQ5oUXXrB99tmnaO/ffvvNhRpTpkxxrx+lbs2aNcvirL/v0qRJExeKFF+bplevXm4mjF5hCrN16tTJ9EqUZsVo++yzz2z33Xe35557zo488siNTrXlllu6V5iuvfbarC9Bp8+aih0RQAABBBBAAAEEEEAAgZwKMD7LKW8kJyeUyYJZoUxJa7sUFhaW+PeaRZPtpoV99frS1VdfvdEhAwcOtPfff98+/fTTbE9lTzzxhJ100knufwcMGOCOC16F+ve//20HH3zwRufadttt7YgjjrC777477TWWL19u+hNsixYtMgU/WpxYoRQbAggggAACCCCAAAIIIIBAxQgQylSMe3lelVAmC81HH300i73+/y5afDfbrbxCmTfffNO6du3qZsOMHj266PKbG8ooLNICxMU3QplsK8x+CCCAAAIIIIAAAggggEBuBAhlcuMa5VkJZaLULuFa5fH6kj6J3blzZ/d61eTJk91ntYNtc19fYqZMBT8gXB4BBBBAAAEEEEAAAQQQSCNAKOP/o0EoU8E11CtFClFeeumljVqihX67d+9u999/f6ktVCfUl5q0jo1eUapdu/ZG+2uh3zp16tjQoUPda1LBFiz0q1k1J554YtYKdPqsqdgRAQQQQAABBBBAAAEEEMipAOOznPJGcnJCmUiY019EX0rSV5H0OlCwQLC+9NS+fXubMGGC6XPW6bYlS5bYgQce6L7epAWBGzZsWOKu+qz2119/7daoCT6JfeONN7qgRp/ETndcSSej01fwA8PlEUAAAQQQQAABBBBAAIH/CTA+8/9RIJSp4BrqC0v6elPjxo3dYr9r1qyxwYMHW9OmTe31118vWkj4mmuuMf2ZP3++tWzZ0n755Rfr0KGDzZ0717TmTYsWLTa6k3333ddq1Kjh/u6tt95y4U2/fv3s9NNPt08++cRd48ILLzSFM2E2On0YLfZFAAEEEEAAAQQQQAABBHInwPgsd7ZRnZlQJirpUq6joGXQoEE2ffp0q1atmpsdM2LECGvQoEHRUcGCu8ECu0HnS3fa4gvxTps2zS6//HIXyGhmzJlnnml///vfrUqVKqEE6PShuNgZAQQQQAABBBBAAAEEEMiZAOOznNFGdmJCmciok3EhOn0y6shdIIAAAggggAACCCCAgP8CjM/8ryGhjP81jPQO6PSRcnMxBBBAAAEEEEAAAQQQQCCtAOMz/x8OQpksanjWWWdlsdfvu+Tl5dnIkSOz3t+3Hen0vlWM9iKAAAIIIIAAAggggEBSBRif+V9ZQpksarj99tsXLbibaXeFMgsWLMi0m7c/p9N7WzoajgACCCCAAAIIIIAAAgkTYHzmf0EJZfyvYaR3QKePlJuLIYAAAggggAACCCCAAAJpBRif+f9wEMr4X8NI74BOHyk3F0MAAQQQQAABBBBAAAEECGUS/AwQypShuL/88os99thjNnPmTFu2bJndddddtuOOO9pzzz1ne+65p+20005lOKsfhxDK+FEnWokAAggggAACCCCAAALJF2B85n+NCWVC1nDRokXWtWtX+/LLL23XXXe1Tz75xN5++21r06aNaUHgwsJCe+CBB0Ke1Z/d6fT+1IqWIoAAAggggAACCCCAQLIFGJ/5X19CmZA1PO6442zOnDk2efJka968uVWvXt3eeecdF8r861//sqFDh9q8efNCntWf3en0/tSKliKAAAIIIIAAAggggECyBRif+V9fQpmQNdx6663twQcftGOPPdY2bNhg1apVKwplZsyYYYcffritXr065Fn92Z1O70+taCkCCCCAAAIIIIAAAggkW4Dxmf/1JZQJWcPatWvb2LFj7bDDDtsklBk/frz96U9/sp9++inkWf3ZnU7vT61oKQIIIIAAAggggAACCCRbgPGZ//UllAlZw4MOOsi23XZbGz169CahzCmnnGJLly61KVOmhDyrP7vT6f2pFS1FAAEEEEAAAQQQQACBZAswPvO/voQyIWv48ssvW8+ePa1Xr17Wv39/GzBggA0fPtz++9//2kMPPWSvvvqqdezYMeRZ/dmdTu9PrWgpAggggAACCCCAAAIIJFuA8Zn/9SWUKUMNJ0yYYBdffLEtWLCg6OiWLVu6T2P37t27DGf05xA6vT+1oqUIIIAAAggggAACCCCQbAHGZ/7Xl1BmM2r4xRdfuNeV6tev7z6PXRk2On1lqDL3iAACCCCAAAIIIIAAAj4IMD7zoUqlt5FQxv8aRnoHdPpIubkYAggggAACCCCAAAIIIJBWgPGZ/w8HoUwWNbzmmmuy2Ov3XfLy8uzKK6/Men/fdqTT+1Yx2osAAggggAACCCCAAAJJFWB85n9lCWWyqGG1atU22WvDhg2b/F2VKlVcKLN+/foszurnLnR6P+tGqxFAAAEEEEAAAQQQQCB5AozP/K8poUzIGn744Yd2zDHH2KBBg+zYY4+1xo0bu3Vlxo4da3feeac999xz1rp165Bn9Wd3Or0/taKlCCCAAAIIIIAAAgggkGwBxmf+15dQJmQNDzzwQDvuuONcKFN8u/32223cuHE2c+bMkGf1Z3c6vT+1oqUIIIAAAggggAACCCCQbAHGZ/7Xl1AmZA1r1aplzz//vPXo0WOTI1988UU76qijbM2aNSHP6s/udHp/akVLEUAAAQQQQAABBBBAINkCjM/8ry+hTMga/uEPf7AOHTrY448/vsmRAwYMsNmzZ9v8+fNDntWf3en0/tSKliKAAAIIIIAAAggggECyBRif+V9fQpmQNXzooYfsjDPOsAMOOMCOPPLIojVlNHvmzTfftFGjRtmpp54a8qz+7E6n96dWtBQBBBBAAAEEEEAAAQSSLcD4zP/6EsqUoYaTJk2y6667zt5991377bffrGrVqrbffvu5T2EffvjhZTijP4fQ6f2pFS1FAAEEEEAAAQQQQACBZAswPvO/voQym1HDgoIC+/77761Ro0aWn5+/GWfy51A6vT+1oqUIIIAAAggggAACCCCQbAHGZ/7Xl1BmM2q4fv16+/nnn22rrbay6tWrb8aZ/DmUTu9PrWgpAggggAACCCCAAAIIJFuA8Zn/9SWUKUMNp02bZkOHDrW3337bNFtGs2Tatm1rw4YNs27dupXhjP4cQqf3p1a0FAEEEEAAAQQQQAABBJItwPjM//oSyoSs4dSpU6137962yy67WL9+/axp06b27bff2tixY23evHn2wgsv2KGHHhryrP7sTqf3p1a0FAEEEEAAAQQQQAABBJItwPjM//oSyoSsob661KRJExs/frzl5eUVHV1YWGh9+/Z1a8zos9hJ3ej0Sa0s94UAAggggAACCCCAAAK+CTA+861im7aXUCZkDWvVqmXjxo2zww47bJMjJ0+ebMcee6ytWbMm5Fn92Z1O70+taCkCCCCAAAIIIIAAAggkW4Dxmf/1JZQJWcMGDRrY7bffbieddNImRz722GP2l7/8xZYtWxbyrP7sTqf3p1a0FAEEEEAAAQQQQAABBJItwPjM//oSyoSsYf/+/W3mzJk2YcIEa9OmTdHR77//vnt9qXPnzjZ69OiQZ/Vndzq9P7WipQgggAACCCCAAAIIIJBsAcZn/teXUCZkDb/55hs76KCDbOHChbbjjju6hX6/++47++KLL2yHHXawGTNmWLNmzUKe1Z/d6fT+1IqWIoAAAggggAACCCCAQLIFGJ/5X19CmTLUcNWqVfbwww+7GTM//fST1atXz82QOfXUU23LLbcswxn9OYRO70+taCkCCCCAAAIIIIAAAggkW4Dxmf/1JZTxv4aR3gGdPlJuLoYAAggggAACCCCAAAIIpBVgfOb/w0Eo438NI70DOn2k3FwMAQQQQAABBBBAAAEEECCUSfAzQCiTRXG1Vky2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RiiTC1XOiQACCCCAAAIIIIAAAgiEF2B8Ft4sbkcQysStIjFvD50+5gWieQgggAACCCCAAAIIIFBpBBif+V9qQhn/axjpHdDpI+XmYggggAACCCCAAAIIIIBAWgHGZ/4/HIQy/tcw0jug00fKzcUQQAABBBBAAAEEEEAAAUKZBD8DhDIJLm4ubo1QJheqnBMBBBBAAAEEEEAAAQQQCC/A+Cy8WdyOIJSJW0Vi3h46fcwLRPMQQAABBBBAAAEEEECg0ggwPvO/1IQyMajhwoULbdCgQfbqq69a9erVrW/fvjZ8+HCrX79+qa37z3/+Yw888IC99dZbNnfuXOvcubNNnz69xGNGjRpld911l33xxRdWp04d69Spk91www224447hhKg04fiYmcEEEAAAQQQQAABBBBAIGcCjM9yRhvZiQllIqMu+UIrV6601q1bW8OGDW3YsGG2evVqGzx4sDVv3txmzpxpeXl5aVt4zTXX2OOPP25t27a1WbNm2XbbbVdiKHP//ffb2Wef7YKfPn362NKlS+3qq6+2devW2ccff+xCmmw3On22UuyHAAIIIIAAAggggAACCORWgPFZbn2jODuhTBTKpVzjtttusyFDhtiCBQtcEKPtjTfesI4dO9rEiROtd+/eaY8uKCiw/Px89/Nu3brZb7/9VmIoo1kx2m/GjBlF53rttdfsoIMOsilTpljPnj2zVqDTZ03FjggggAACCCCAAAIIIIBATgUYn+WUN5KTE8pEwpz+Il26dLGqVavatGnTNtqpVatW1qNHDxs5cmRWLSwtlGnXrp01atTIhTzB9uGHH9o+++xjkyZNssMPPzyra2gnOn3WVOyIAAIIIIAAAggggAACCORUgPFZTnkjOTmhTCTM6S/SpEkT69+/v91+++0b7dSrVy9bsWKFe4Upm620UOahhx5yry89/PDDbr0avb503nnn2eLFi+29995z69hku9Hps5ViPwQQQAABBBBAAAEEEEAgtwKMz3LrG8XZCWWiUC7lGgpE9PqS1nhJ3QYOHGjvv/++ffrpp1m1sLRQRifQgsB//vOfbf369e58e++9t3t1aZtttin1/MuXLzf9CbZFixa5RYK1OPH222+fVdvYCQEEEEAAAQQQQAABBBBAoPwFCGXK3zTqMxLKRC1e7HpRhDLPP/+8DRgwwC30q1eilixZYlokuEaNGm4mzpZbbplWQWGRFiAuvhHKVPCDw+URQAABBBBAAAEEEECg0gsQyvj/CBDKVHANc/36UmFhoTVt2tSOOuoou++++4ru9ssvvzStW3PnnXe6GTTpNmbKVPADwuURQAABBBBAAAEEEEAAgTQChDL+PxqEMhVcw4MPPtit6fLSSy9t1BIFJt27dzd9zjqbLd3rS1o/RsGPFgw+66yzNjqVFv898cQT7Y477sjmEm4fOn3WVOyIAAIIIIAAAggggAACCORUgPFZTnkjOTmhTCTM6S9y66232hVXXOHWaGnWrJnbcfbs2da+fXubMGGC9enTJ6sWpgtlNFOmTp067vWl1C85qfPusMMONnz4cLvooouyugahTNZM7IgAAggggAACCCCAAAII5FyAUCbnxDm/AKFMzolLv4C+sNS6dWtr3LixW+x3zZo1NnjwYPfK0euvv255eXnuBFoDRn/mz59vLVu2dH/3/fff24wZM4p+vmHDhqL1X/bff/+i/S699FIXvvztb38zhTdaU+a6665zx2shYV07241On60U+yGAAAIIIIAAAggggAACuRVgfJZb3yjOTigThXKGayho0SK806dPt2rVqrnZMSNGjLAGDRoUHRksuJu6wK7279KlS4ln1+ev//SnP7mf/frrr+4VJf2dOu3WW29tBxxwgF1//fW22267hRKg04fiYmcEEEAAAQQQQAABBBBAIGcCjM9yRhvZiQllIqNOxoXo9MmoI3eBAAIIIIAAAggggAAC/gswPvO/hoQy/tcw0jug00fKzcUQQAABBBBAAAEEEEAAgbQCjM/8fzgIZfyvYaR3QKePlJuLIYAAAggggAACCCCAAAKEMgl+BghlElzcXNwaoUwuVDknAggggAACCCCAAAIIIBBegPFZeLO4HUEoE7eKxLw9dPqYF4jmIYAAAggggAACCCCAQKURYHzmf6kJZfyvYaR3QKePlJuLIYAAAggggAACCCCAAAJpBRif+f9wEMr4X8NI74BOHyk3F0MAAQQQQAABBBBAAAEECGUS/AwQyiS4uLm4NUKZXKhyTgQQQAABBBBAAAEEEEAgvADjs/BmcTuCUCZuFYl5e+j0MS8QzUMAAQQQQAABBBBAAIFKI8D4zP9SE8r4X8NI74BOHyk3F0MAAQQQQAABBBBAAAEE0gowPvP/4SCU8b+Gkd4BnT5Sbi6G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" width="1000">


.. parsed-literal::

    2. Reporting Generalized Performance:
    
    |                  |            0 |
    |:-----------------|-------------:|
    | clump_p1         |    1         |
    | clump_r2         |    0.1       |
    | clump_kb         |  200         |
    | p_window_size    |  200         |
    | p_slide_size     |   50         |
    | p_LD_threshold   |    0.25      |
    | numberofpca      |    6         |
    | tempalpha        |    0.1       |
    | l1weight         |    0.1       |
    | Train_pure_prs   | -408.191     |
    | Train_null_model |    0.233856  |
    | Train_best_model |    0.252479  |
    | Test_pure_prs    | -425.092     |
    | Test_null_model  |    0.137656  |
    | Test_best_model  |    0.156992  |
    | pvalue           |    1         |
    | Difference       |    0.0954873 |
    | Sum              |    0.409471  |
    3. Reporting the correlation of hyperparameters and the performance of 'Train_null_model', 'Train_pure_prs', 'Train_best_model', 'Test_pure_prs', 'Test_null_model', and 'Test_best_model':
    
    3. For string hyperparameters, we used one-hot encoding to find the correlation between string hyperparameters and 'Train_null_model', 'Train_pure_prs', 'Train_best_model', 'Test_pure_prs', 'Test_null_model', and 'Test_best_model'.
    3. We performed this analysis for those hyperparameters that have more than one unique value.
    


.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAABUYAAAOECAYAAABtsevkAAAAAXNSR0IArs4c6QAAIABJREFUeF7s3QecFEX+//8PIGaUE1AkCBgwYhZJBlQwgVlUDBjAjKinYs6YFfOdgUNUQDGLngFFPQQBURAwE0xgjqgY8f9/1/dX6+wwuzu9bM92d73q8eDhsfRUVz+rpvfmPdVVdf7666+/jIIAAggggAACCCCAAAIIIIAAAggggAACCAQkUIdgNKDe5lIRQAABBBBAAAEEEEAAAQQQQAABBBBAwAkQjDIQEEAAAQQQQAABBBBAAAEEEEAAAQQQQCA4AYLR4LqcC0YAAQQQQAABBBBAAAEEEEAAAQQQQAABglHGAAIIIIAAAggggAACCCCAAAIIIIAAAggEJ0AwGlyXc8EIIIAAAggggAACCCCAAAIIIIAAAgggQDDKGEAAAQQQQAABBBBAAAEEEEAAAQQQQACB4AQIRoPrci4YAQQQQAABBBBAAAEEEEAAAQQQQAABBAhGGQMIIIAAAggggAACCCCAAAIIIIAAAgggEJwAwWhwXc4FI4AAAggggAACCCCAAAIIIIAAAggggADBKGMAAQQQQAABBBBAAAEEEEAAAQQQQAABBIITIBgNrsu5YAQQQAABBBBAAAEEEEAAAQQQQAABBBAgGGUMIIAAAggggAACCCCAAAIIIIAAAggggEBwAgSjwXU5F4wAAggggAACCCCAAAIIIIAAAggggAACBKOMAQQQQAABBBBAAAEEEEAAAQQQQAABBBAIToBgNLgu54IRQAABBBBAAAEEEEAAAQQQQAABBBBAgGCUMYAAAggggAACCCCAAAIIIIAAAggggAACwQkQjAbX5VwwAggggAACCCCAAAIIIIAAAggggAACCBCMMgYQQAABBBBAAAEEEEAAAQQQQAABBBBAIDgBgtHgupwLRgABBBBAAAEEEEAAAQQQQAABBBBAAAGCUcYAAggggAACCCCAAAIIIIAAAggggAACCAQnQDAaXJdzwQgggAACCCCAAAIIIIAAAggggAACCCBAMMoYQAABBBBAAAEEEEAAAQQQQAABBBBAAIHgBAhGg+tyLhgBBBBAAAEEEEAAAQQQQAABBBBAAAEECEYZAwgggAACCCCAAAIIIIAAAggggAACCCAQnADBaHBdzgUjgAACCCCAAAIIIIAAAggggAACCCCAAMEoYwABBBBAAAEEEEAAAQQQQAABBBBAAAEEghMgGA2uy7lgBBBAAAEEEEAAAQQQQAABBBBAAAEEECAYZQwggAACCCCAAAIIIIAAAggggAACCCCAQHACBKPBdTkXjAACCCCAAAIIIIAAAggggAACCCCAAAIEo4wBBBBAAAEEEEAAAQQQQAABBBBAAAEEEAhOgGA0uC7nghFAAAEEEEAAAQQQQAABBBBAAAEEEECAYJQxgAACCCCAAAIIIIAAAggggAACCCCAAALBCRCMBtflXDACCCCAAAIIIIAAAggggAACCCCAAAIIEIwyBhBAAAEEEEAAAQQQQAABBBBAAAEEEEAgOAGC0eC6nAtGAAEEEEAAAQQQQAABBBBAAAEEEEAAAYJRxgACCCCAAAIIIIAAAggggAACCCCAAAIIBCdAMBpcl3PBCCCAAAIIIIAAAggggAACCCCAAAIIIEAwyhhAAAEEEEAAAQQQQAABBBBAAAEEEEAAgeAECEaD63IuGAEEEEAAAQQQQAABBBBAAAEEEEAAAQQIRhkDCCCAAAIIIIAAAggggAACCCCAAAIIIBCcAMFocF3OBSOAAAIIIIAAAggggAACCCCAAAIIIIAAwShjAAEEEEAAAQQQQAABBBBAAAEEEEAAAQSCEyAYDa7LuWAEEEAAAQQQQAABBBBAAAEEEEAAAQQQIBhlDCCAAAIIIIAAAggggAACCCCAAAIIIIBAcAIEo8F1OReMAAIIIIAAAggggAACCCCAAAIIIIAAAgSjjAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQCE6AYDS4LueCEUAAAQQQQAABBBBAAAEEEEAAAQQQQIBglDGAAAIIIIAAAggggAACCCCAAAIIIIAAAsEJEIwG1+VcMAIIIIAAAggggAACCCCAAAIIIIAAAggQjDIGEEAAAQQQQAABBBBAAAEEEEAAAQQQQCA4AYLR4LqcC0YAAQQQQAABBBBAAAEEEEAAAQQQQAABglHGAAIIIIAAAggggAACCCCAAAIIIIAAAggEJ0AwGlyXc8EIIIAAAggggAACCCCAAAIIIIAAAgggQDDKGEAAAQQQQAABBBBAAAEEEEAAAQQQQACB4AQIRoPrci4YAQQQQAABBBBAAAEEEEAAAQQQQAABBAhGGQMIIIAAAggggAACCCCAAAIIIIAAAgggEJwAwWhwXc4FI4AAAggggAACCCCAAAIIIIAAAggggADBKGMAAQQQQAABBBBAAAEEEEAAAQQQQAABBIITIBgNrsu5YAQQQAABBBBAAAEEEEAAAQQQQAABBBAgGGUMIIAAAggggAACCCCAAAIIIIAAAggggEBwAgSjwXU5F4wAAggggAACCCCAAAIIIIAAAggggAACBKOMAQQQQAABBBBAAAEEEEAAAQQQQAABBBAIToBgNLgu54IRQAABBBBAAAEEEEAAAQQQQAABBBBAgGCUMYAAAggggAACCCCAAAIIIIAAAggggAACwQkQjAbX5VwwAggggAACCCCAAAIIIIAAAggggAACCBCMMgYQQAABBBBAAAEEEEAAAQQQQAABBBBAIDgBgtHgupwLRgABBBBAAAEEEEAAAQQQQAABBBBAAAGCUcYAAggggAACCCCAAAIIIIAAAggggAACCAQnQDAaXJdzwQgggAACCCCAAAIIIIAAAggggAACCCBAMMoYQAABBBBAAAEEEEAAAQQQQAABBBBAAIHgBAhGg+tyLhgBBBBAAAEEEEAAAQQQQAABBBBAAAEECEYZAwgggAACCCCAAAIIIIAAAggggAACCCAQnADBaHBdzgUjgAACCCCAAAIIIIAAAggggAACCCCAAMEoYwABBBBAAAEEEEAAAQQQQAABBBBAAAEEghMgGA2uy7lgBBBAAAEEEEAAAQQQQAABBBBAAAEEECAYZQwggAACCCCAAAIIIIAAAggggAACCCCAQHACBKPBdTkXjAACCCCAAAIIIIAAAggggAACCCCAAAIEo4wBBBBAAAEEEEAAAQQQQAABBBBAAAEEEAhOgGA0uC7nghFAAAEEEEAAAQQQQAABBBBAAAEEEECAYJQxgAACCCCAAAIIIIAAAggggAACCCCAAALBCRCMBtflXDACCCCAAAIIIIAAAggggAACCCCAAAIIEIwyBhBAAAEEEEAAAQQQQAABBBBAAAEEEEAgOAGC0eC6nAtGAAEEEEAAAQQQQAABBBBAAAEEEEAAAYJRxgACCCCQIIG77rrLjjjiCHvhhRds++23r/GWxV1/jTc4p8L//e9/dtZZZ9nMmTPthx9+sAsuuMAuvPDCOE9J3QikQqBOnTrWp08f0/u7mKL3zUUXXRTbfaaYNlTnmDTfv6pzvUl9DffipPYM7UIAAQQQQACB6ggQjFZHjdcggEAiBf766y97+OGH7e6777bJkyfb119/bcsvv7xtsMEGtscee9jRRx9tq6yySiLb7htVEx/8Fai+9NJLJo/8UhP11wbgt99+a2uttZY7tQKglVde2QXHlYXHCou23nprmzhxYsEm//vf/7bjjjvOhg4daocffnhtXBbnjFlAY2C77bazF198MeYz1W71hYLR1q1bu0Z98MEHizVuSYNRfx8p9qpr6ouemr5/nX/++XbFFVfYhhtuaNOmTSvqclq1alXQtKgXFzhI/dOmTZtIwbaq0T1r2LBhZTVqDOi+uOmmm7r7Wq9evarbpEpfV517cSwNoVIEEEAAAQQQQKCGBAhGawiSahBAoHYFFixY4D4IPv3009awYUPbZZddTMHAzz//bJMmTXJBqX7+zTff1G5Dqzh7TXzwrywY/f777+3TTz+1NdZYw4XGaSnPPvus7bzzzi7EGDhwYFHNJhgtiinTB4USjL7zzjsuFFt99dXL+jPOYFQh4qOPPlpu7Ch81hcye+65pwvncotCPN+eJRlwNX3/2myzzewf//iHHXbYYYuFnddff73pfJqZnlv0e+Tkk09eksso99olDUZPOOEEa9y4sf3xxx/2/vvvuy8H9b8vu+wyN8O+pkt17sU13QbqQwABBBBAAAEEalKAYLQmNakLAQRqTUAfxh9//HHbd999bciQIS4kyC1Tp061Y445xgWkSS5xB6NJvvbK2qZZwJopGmV2J8FoWnu75todSjBaSCzOYLTQ+fws1Cjv0Zrr6eg1zZs3z1q0aGHXXnutnXrqqYtVIL8PP/yw4Mz76Ger+BVLGoy+/fbbtt5665WdQLNzd9ppJ1tmmWXsyy+/tBVWWKEmm+ueyIh6L67RBlAZAggggAACCCBQwwIEozUMSnUIIFB6Ac0S3XXXXa1du3Y2ZcoUW3rppQs24rfffiv3b7/++qtdddVVNnz4cDdbqEGDBta1a1e7+OKLy33QVGU+ZFD9Z5xxhj355JPuQ+ecOXPcY7paF1SBgGZhXnnllfbWW2/ZAQccULbmn/5+ySWXuDX99Chiy5Yt7ZBDDnEzevQB1pdCwahmvd566602evRoe/fdd93rNTNMywNonUDNePJFQVCh4sOKioLX6lho1thpp53mZo6pjZ07d7abbrppMbvKRsSIESPsxhtvtBkzZli9evVMM7jOPPNM159VXdPcuXMrnYUWNRhduHChNW3a1NZZZx03jvKL1jbVGNPsN3nmjotXX33VTjnlFHvqqadMll26dLHrrrvOLeOQX1555RW7/PLLbfz48fbTTz/Z2muvbf369bOTTjrJcvsv93FnjZ+bb77ZZs2aZWeffbZbW9U/SqtZYgor1LdffPGFO6dmuenLgtyi86nd48aNs08++aTMW+N59913L3esxrTeC6pnm222cefTlwtbbrmlG+8a+//617/sv//9r82ePdut+apHjDXmzznnHFt22WXL1effP6+//robMxrLctphhx3c2G7WrJmbbXjuueeajllppZXs2GOPNT3qnD+mNdauueYau++++9z7T8GP6hk0aJC1bdvWnde3v9DYyx03mkGu12kMy0Tvpd12280uvfRS1yZfcsOrAQMGuDGqJRq0NIfqU9F9ROP/vffes19++cWaNGliHTp0cHbrr79+hW8DjYdOnTotNhu6f//+rs8PPPBAGzlyZNnr9cVP3759nWGPHj3cz3MfpfdtLXRC/0h77tjS+NFYlaX6UNd25JFHVva2LfhvFQWjfgb7jz/+aOedd56NGjXKzVp//vnn3VIYGkP33nuvm9mvsHK55ZZzbjpWLrml0P0rd6zuuOOO7p7qx9DBBx/sXOvXr79Ym2+77TY3xnRP9eMm96CKgtFixp+vR32r8aX2aGkXjZeNNtrINIbUd5UtSVBoKZTc9vn3f34wqmN0D9DP9UXgVltt5V5WE/cd/b6p7D21aNEiu+WWW+zOO+907wPdBzp27OjuI1rWJLdUNi50nL//6L+65+n3jd6fmiGrPta51Lc61/z58917TL9LdL/KLdUZX7pPamaw/r+A7r2rrrqqW56g0AxcjQe9f+6//353L9Q1r7vuuu73u9qaW+655x5339TvO/Wvft/pnqenISgIIIAAAgggUHsCBKO1Z8+ZEUCghgQOOuggF5L85z//cQFlMUUfqvS4/ZgxY2yLLbZwM2wUjOhDuz6YK6TJfRxUH5IV5Cjs0AcafQDX4/sKO/Vooc7bvXt397q99trLBQwK2RSW6WcKWxRe6N/0c4UqCqr08yeeeKIs/Cn0wV+PyW688cYuRFBop/BVH6x03k022cTV5cNVhROqQzOdch8B1Xl1PYXqr46FQubmzZvb77//7j7AKtB48MEH3c/U3mJm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GoIAAqUXIBgtvTlnjChAMBoRjMMRQAABBBBAAAEEEEAgUwIEo/F1J8FofLbUjEAaBAhG09BLgbeRYDTwAcDlI4AAAggggAACCCAQuADBaHwDgGA0PltqRiANAgSjaeilwNtIMBr4AODyEUAAAQQQQAABBBAIXCBrwegbH7VMTI9ussbHiWkLDUEAgdILZCoYrVOnTpWCffr0sbvuuqvK4yo7QK+vW7euHXbYYUtUDy8uTiCrwWiLlVay87fdwTq2aGm/L/rTnpsz2waNe8m+//WXKmE6tVjDzui8jbVt1Mi+WbjQ7ps5w26dMskW/fVXudf2WGddO6F9B2u9ckP79McFNmTqazZ8xhtV1s8ByRNgvCSvT5LaIsZKUnsmme1ivCSzX5LYKsZKEnsluW1ivNR83xCM1rypr5FgND5bakYgDQKZCkYnTpxYznz//fe3jTfe2M4777yynzdp0sTWWmutJeqb7bff3pZaail77rnnlqgeXlycQBaD0RXq17enDu5j3y5caNdPmmDL1a9vAztta5//9KP1evC+SmHarbqajdr/QHtm1vt2/5szrG2jxjaw8zYu9Lz2lfFlr92xzVp2R8+9bOjU12zM3NnWvlkL69++g537wnPudZT0CDBe0tNXtd1Sxkpt90C6zs94SVd/1WZrGSu1qZ++czNe4umzrAWjr3+0RjxQ1ah18zU+qsareAkCCGRFIFPBaH6nrL322talS5clniGaX29Sg9FFixaZ/ii0rc3y66+/2jLLLFNjTchiMHrUZlvYaZ262PZ3DXFhqMrmTZvZg70Osr6PP2JjP5hToZ/CzuYNVrLdR9xtfn7o8Vu2t/7tO1rH/9xm3/3yfzNOn+p9mM1fsMCOGv1IWV2DduhmO7VZyx2XP7u0xjqMimpcgPFS46SZrZCxktmujeXCGC+xsGayUsZKJrs1totivMRDSzAaj6v7HEYwGh8uNSOQAoHggtHXXnvNzj33XJswYYL98ccfLji97rrrbMMNNyzrrieffNIuvvhie+utt0yP57dq1cpOOukk69evnykUfemll8p1bbGP5/tA9cgjj7QLL7zQPvroIzej9eabb7b27duX1dm6dWs7/PDD3TG+fPDBB9amTRu755577JBDDnE/1nE77bSTbbrppnb99debjnn11Vdts802c7NZdQ26XgWl3bp1s8GDB1vLlsWt5eLPd/vtt9uUKVNs1KhRLnTdb7/97IYbbrAVV1zRteHFF1+0rl272mOPPeaOkZ28pk2bZq+88oqdffbZNnXqVPv999+tRYsWduihhzr/KCWLwejwffa3Pxf9ZYc9+mA5ipf69LVxH33gZnUWKvXr1rXpx/a3myZPdI/O+9KsQQN7+YijbcDTT9ro994x//dTnnnSHnv3nbLj2jdvYffte4DtO2qETf3s0yjdwLG1KMB4qUX8lJ2asZKyDqvl5jJearkDUnR6xkqKOisBTWW8xNMJBKPxuKpWgtH4bKkZgTQIBBWMKuDbZpttXBh6/PHHW7169eyqq66yt99+22bMmGHNmjWz2bNn2/rrr28HHHCAC/G0lqgC0t9++81OO+00978VTCpsvPHGG10fF/t4voLRd99911ZYYQW75JJLbNlll7VLL73UZs2a5f6oHh94FhuMql3Nmze3gQMHWoMGDVxIOmnSJNt7771NSwnoGn766ScXsurYN954w52/quKDUdUtL7XnnXfecUHnnnvuaSNHjnRV+GBUdgpNe/ToYX/++ad17tzZ1lhjDevYsaOdcMIJttxyy7lrnDt3rl1++eVVnb7cv2cxGJ3c91gb/e47dsm4F8td65Cee9uKSy9tBzx0f0Gjtf6xio059Ag75onHbMycWeWOmXncSTZk6hQbPHGCbdeqtQ3dc183q/Ttr74sO26V5ZazKf2OtzOfe8ZGvTUzUj9wcO0JMF5qzz5tZ2aspK3Hare9jJfa9U/T2Rkraeqt2m8r4yWePshaMPrqR63jgapGrVut8UE1XsVLEEAgKwJBBaM77LCDffPNN24GpH/c/IcffrA111zTjjrqKLvyyivtwQcfdIHi999/byuttFLBfq7uo/R+tqnCSc0UVfniiy/cDMtTTjnFLrvsMvezKDNGv/32W/vwww+tYcOGZW3VGqoKd5944omynymQXHfddd1sz+OOO67K8euD0S222MJ5+aIw+OSTT3YB8XrrrVcWjCqAvfvuu8uO02u22morF8T6a63ypBUckMVg9J0TTrZ/TZlkN0x6pdxVX9d9V9ugyaq2y/BhBTX84/YHPXS/TZr3Sbljxh95tNvA6YIXn7c92q5n1++yu20z9A6bt+CHsuPq1alj7/c/1S5/+SW74/W/+7W6fcPrSiPAeCmNcxbOwljJQi+W7hoY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upNgND5bakYgDQIEo2nopcDbSDAa+ADg8hFAAAEEEEAAAQQQCFwga8HoY3M2TUyP7rnmtMhtmTt3rg0YMMDGjh1rSy+9tO2xxx523XXX2SqrrFJlXc8//7ydeeaZNnPmTGvSpIn169fPzj77bKtXr5577QcffGBt2rSpsJ5XXnnFOnTo4P798MMPt2HDhi127E033WQnnnhilW3hAAQQMCMYZRQkXoBgNPFdRAMRQAABBBBAAAEEEEAgRgGC0fhwowajCxYssHbt2lnjxo3toosusp9++skGDhxozZs3t3HjxlmdOnUqbOyUKVOsc+fOtu+++1rfvn1dOKrXnnrqqTZo0CD3ul9//dWmTp26WB3HH3+8zZ8/3+bNm1cWoioYVdD6wAMPlDtewepqq60WHxo1I5AhAYLRDHVmVi+FYDSrPct1IYAAAggggAACCCCAQDECBKPFKFXvmKjB6LXXXutmeM6ZM8eFoSoTJkxwgefo0aOtR48eFTakZ8+e9uGHH9q0adOsbt267rjLLrvMLr74Yhd4NmrUqOBrP//8c3eu/v372+DBg8uOUTD68ssv26xZs6p38bwKAQSYMcoYSL4AwWjy+4gWIoAAAggggAACCCCAQHwCWQtGH569WXxYEWveZ63FZ2dWVkXXrl1tqaWWsjFjxpQ7TLM0u3fvbrfddlvBl//222+20kor2XnnnWfnnHNO2TEKSlu3bm0jRoywgw46qOBrFYZqVunrr79um232tx3BaMTO5nAECggwY5RhkXgBgtHEdxENRAABBBBAAAEEEEAAgRgFCEbjw40ajOoRdQWY119/fblG7b777vbDDz+4x+kLlbfffts22GADe+SRR2yvvfYqd8gKK6zggs9LLrmk4Gs333xzU7CqR+9zi4LRkSNHml6vR/zXXXddt/ap1i2lIIBAcQIEo8U5cVQtChCM1iI+p0YAAQQQQAABBBBAAIFaFyAYja8Ldmj0gn333XflTtCwYUPTn0JFmy3pUfoLL7yw3D8fcsghbm3QN998s+Dr/OP2L7zwgm2//fbljmnRooXtueeedssttyz2WtW30UYb2ZVXXmlnnHFGuX+/4YYb3OzVDTfc0IWyw4cPt1GjRrnH888666z40KgZgQwJEIxmqDOzeikEo1ntWa4LAQQQQAABBBBAAAEEihHIWjD6wOwtirnskhzz5j093SZKueWCCy5YLPj0/17qYFQ72F999dX20Ucfla1pWhmMNnZ65pln7KuvvrJll122JIacBIE0CxCMprn3Amk7wWggHc1lIoAAAggggAACCCCAQEEBgtH4Bka3Rs9HmjFaykfpFy1aZK1atXKPyD/33HNFIWjG6AEHHLDYeqRFvZiDEAhQgGA0wE5P2yUTjKatx2gvAggggAACCCCAAAII1KRA1oLR+2dtVZM8S1TXAWu/Gun1egxes0afffbZcq/T5kvdunWz22+/vWB9WiO0QYMGptmoehTfF7/5kh6D7927d7nXjh071nbccUcbNmyYHXbYYUW10wejeqx/0003Leo1HIRAyAIEoyH3fkqunWA0JR1FMxFAAAEEEEAAAQQQQCAWAYLRWFhdpVGD0WuuucbtKj937lxr1qyZq2PixInWsWNHe/zxx61nz54VNrZHjx728ccfu7VI69at64674oorXFg6b948a9y4cbnXHnHEEfbAAw/YZ599ZiuuuGJRCHvvvbebXfrll1/yKH1RYhwUugDBaOgjIAXXTzCagk6iiQgggAACCCCAAAIIIBCbAMFobLSRg1FtctSuXTtbddVV3TqkP//8sw0cONCaNm1q48ePtzp16rjGXnzxxe7P7Nmz3ePwKpMnT7YuXbpYr1697KijjnK7zOu1J510kgtIc8vChQtNj+1rU6Z77rlnMQDNNO3Tp48deOCBtvbaa7td6e+99157+OGH3Zqkp512Wnxo1IxAhgQIRjPUmVm9FILRrPYs14UAAggggAACCCCAAALFCGQtGB0xa+tiLrskx/Ree1Lk8yjsHDBggL344otWv359N0t08ODB1qhRo7K6FJpqUyfNLG3dunXZz8eMGeN2jFcoqhmi/fr1s3PPPdfq1atXrh0jR450j9ZrI6Xu3bsv1sZvvvnGjjzySLeW6BdffOFer8C2f//+dvDBB0e+Jl6AQKgCBKOh9nyKrptgNEWdRVMRQAABBBBAAAEEEECgxgUIRmuctKzC6gSj8bWGmhFAoNQCBKOlFud8kQUIRiOT8QIEEEAAAQQQQAABBBDIkADBaHydSTAany01I5AGAYLRNPRS4G0kGA18AHD5CCCAAAIIIIAAAggELpC1YPSe9zskpkcPXWdiYtpCQxBAoPQCBKOlN+eMEQUIRiOCcTgCCCCAAAIIIIAAAghkSoBg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WsZKbeqn79yMl3j6jGA0HlfVWp1gdO7cuTZgwAAbO3asLb300rbHHnvYddddZ6usskqVDX3++eftzDPPtJkzZ5omi/Xr18/OPvtsq1evXtlrDz/8cBs2bPFNoW666SY78cQTy53jvvvus0svvdRmzZplLVu2tFNPPdWOO+64KtvBAQgg8H8CtRaM+kDwscces1GjRtmTTz5prVq1smnTptnTTz9tN9xwg5sRumDBAjcr9KSTTrKjjjqqrN8KPUrvH9U+5ZRT7KyzznI3hvXXX9+uv/5622abbYrucwWjF110kf3yyy82ZMgQ+/HHH22XXXaxW2+91VZbbbWyen7//Xe7/PLL7e6777aPPvrImjVrZkcffbQ7t19K4Pvvv7eBAwfaE088YV9++aW7UW6xxRbuJjdjxgw3WzS/6CbbunXrStvr/UaPHm0PPPCAPfrooy5U1k1QN1m5nnPOOSanLbfc0v7zn//YmmuuWVbnd99954575JFHTP+7bdu2dsYZZ9ihhx5a7rxqt34+Z84cW3fddd3N/phjjlnsUfrXXnvNLT0wYcIE++OPP9y/69gNN9ywrD4epf+b9qjNtrDTOnWx7e8a4sJQlc2bNrMHex1kfR9/xMZ+MKfC/lfY2bzBSrb7iLvNzw89fsv21r99R+v4n9vsu1/+b8bpU70Ps/kLFthRox8pq2vQDt1spzZruePyZ5cW/QbhwJILMF5KTp7aEzJWUtt1tdJwxkutsKfypIyVVHZbrTWa8RIPfdaC0SNePSIeqGrUOnSroZFepYyiXbt27slV5QY//fST+8zfvHlz9yRoZcsKTpkyxTp37mz77ruv9e3b14Wjeq0+xw8aNKisHQpGFaDqs35u0XKCuZnE448/bnvuuacLafXf//3vf3bxxRfbv//9bxe4UhBAoGqBWg9GFSbut99+1qNHD/vzzz/LAkjNHlWoWb9+fRs/frz7BkQB5/HHH++uqqJg9L333nPfuiicXGmlleyCCy6w999/3x3fsGHDqkWUFtep425qOr9uMJ9//rm7WSk8VPDni5YHeOaZZ1wAqbBz8uTJ7iZ08skn2xVXXOEOU5ir0FcBqpYN+OKLL9wamwoR1b57773XTjjhBLvtttvcUgMqm222WZUzZ30wqgC1d+/eptBRIee//vUvO/3009051C4VXYOOe/nll93f5ayg+M0333Q3XwXPI0eOdAHv7bffXnYDnT59uruuHXbYwfr372+fffaZW3v1hx9+cDddvz6sbu6qT2Go+kffdF111VX29ttvu/BXfaxCMPr38Bu+z/7256K/7LBHHyw3Jl/q09fGffSBm9VZqNSvW9emH9vfbpo80T0670uzBg3s5SOOtgFPP2mj33vH/N9PeeZJe+zdd8qOa9+8hd237wG276gRNvWzT4t6P3BQ7QswXmq/D9LSAsZKWnoqGe1kvCSjH9LQCsZKGnopOW1kvMTTFwSj8biq1qjB6LXXXutmeGrykHIDFeUECjw1cUnZRkWlZ8+e9uGHH7oJYXXr1nWHXXbZZS5HmDdvnjVq1Mj9TMGoPr9rsldlRQGt9ipR5uCLJmspMFV9ubNQ4xOkZgTSLVDrwahmKCqQq6gsWrTI9EczRnVjUFinUlEwqrVNFY76jYxef/11F+7df//91qtXr6J6S8HoXHm9SwAAIABJREFU6quvbpq56R/t99/EKAjt3r27+yZmu+22swcffNB92+OLgkaFh59++qn94x//sI022sgdr9mThUp1H6X3r9M0ek2nV1HgqRBSs1R1k/aB5M033+yCTbWpadOm7matqf6aqatw15edd97ZBZm6gcrgwAMPdOaqS48HqHiH3I2zFJx+8803poBUmyupKDzVDFUFw1deeaX7GcHo3yNgct9jbfS779gl414sNyyG9NzbVlx6aTvgofsLjpe1/rGKjTn0CDvmicdszJzyvyRnHneSDZk6xQZPnGDbtWptQ/fc180qffurL8vqWmW55WxKv+PtzOeesVFvzSzq/cBBtS/AeKn9PkhLCxgraempZLST8ZKMfkhDKxgraeil5LSR8RJPXxCMxuOqWqMGo3rqU597x4wZU65Rms2pz/6a9FSo/Pbbb25ylPZd8ZOYdJyCUk1kGjFihB100EHupcUEo/51mmx18MEHl53ypZdecp+9FdZ27NgxPjhqRiAjArUejGpquGaM5pb58+fb+eef72ZjKsxT4KeikFKPt6tUFIxqWrse6/ZFM0+XXXZZu/rqq+20004rqtsUCupbltwb2l9//WXLLbec+2ZIbdONTI/7f/311+W+hdHj/+3bt3dt101RNzSFiZpxqr9vuumm5abWL2kwqrr1rZMvnTp1Mt1wFVL68uyzz5pCz1dffdU9Vq8ZpTfeeKMt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g3zarpM8vnK5jNu4Pk2zz0mJF2BeEt8DXypgVnzpVDDqZF6C0QcfqmBWfOhScGpkXtz0gmDUjatZlWDUnS0rI+CDQCCDUXNruLli9OTjiy++kN69e8vChQttwGdCQHOY28eP/veZiikFo+bWe3Ord+QwV56efvrp8vTTT0vnzp3T1CcTFJorRIcMGRI939RirmZ9//337e3qibiV3gSVNWvWjNZUtmxZ+e2332xIGTnefvttMaFnpE5zRemzzz4rR44ciYa85tzILf87duywt9YXLVpUqlWrZh8XEDlOnDgh2bJls6GpOd+skT17dnslqVn35OOuu+6S7777zt5ebw6C0f/p8P8wpunHjpP+K8C8MAppFWBW0irFeUaAeWEO0irArKRVivP4u8XdDIQtGC065XF3WDGuvKfxYzF+gtMRQCBMAoEMRt99990kz/80V4def/318s0330j37t3lkksusWGcue173LhxYsI6c/zdM0bNreUnH39+xmdqTU3u/MhzPlesWGFvI09EMBr57kj9yYWPf66zZcuWNmA2zzw9+TDPcb311lvFXGF7ww032PDYhJ0m9Dz5KFSokNx44402GDVXl0aesZqcoQlYTdBqjrQEo+Z5p+afk49K01+RTGdkS61FXv05txh51a6EF8u8JLwF3hTArHjTqkAUyrwEog1eFMGseNGmwBTJvLhpBcGoG1ezKsGoO1tWRsAHgUAGo38O+8zzQE3AZp6Daa5UjBzmNm3zUqGgBKPmGZznn3++DR3NlZaRw1xlaq5MdfHypfQEoybsNLff//LLL8leMWoeO3DxxRen6YrRn3/+2T631VxNW79+/b/MvAlXr7jiijQHo3379pV+/folWSdn9ZvlnFtv8eHnKc01Tq1dT/7v+DFp8vrsJJ8xL19a+fk+6fHOomTXsi9feqidjFz7noxevy56TuRlSw+/tUDmffKxfWv9imatJPLryImRly/dPXOabPzqizTXy4mJFWBeEuvv07czKz51K/G1Mi+J74EvFTArvnQqGHUyL276QDDqxpVg1J0rKyPgi4AXwegHH3xgX3o0e/ZsqV27trU1z7A0t3qbZ1sGJRiN3GY+ePBgGxRGDvMMTvOm+qAEo5FnjM6aNUvq1KkTrdNcLWqsY33GaMWKFeWMM84Qc8Xp3x1cMfo/nZbXXCePlC0nlV4aL9/8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" width="1200">



.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAABpcAAAaYCAYAAABirITVAAAAAXNSR0IArs4c6QAAIABJREFUeF7s3Qu0nVV16PF5gsGEhIchgheBXsRK0EYhA1FrKW8kCpIgWFLyQKRoIalQCBAEFUolIEEDCVDKQ6CIXIIEAggtD4O8JLVAEFNpARVsFUpqIs8UOXd8ew9ScrLPzlpzrTW/9a39P2N03FtYc61v/+ac+6zD7D6nr7+/v1/4QgABBBBAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQMBBoI/hkoMSSxBAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABBFoCDJcoBAQQQAABBBBAAAEEEEAAAQQQQAABBBBAAAEEEEAAAWcBhkvOVCxEAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABBBBguEQNIIAAAggggAACCCCAAAIIIIAAAggggAACCCCAAAIIOAswXHKmYiECCCCAAAIIIIAAAggggAACCCCAAAIIIIAAAggggADDJWoAAQQQQAABBBBAAAEEEEAAAQQQQAABBBBAAAEEEEDAWYDhkjMVCxFAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABBBguUQMIIIAAAggggAACCCCAAAIIIIAAAggggAACCCCAAALOAgyXnKlYiAACCCCAAAIIIIAAAggggAACCCCAAAIIIIAAAgggwHCJGkAAAQQQQAABBBBAAAEEEEAAAQQQQAABBBBAAAEEEHAWYLjkTMVCBBBAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABhkvUAAIIIIAAAggggAACCCCAAAIIIIAAAggggAACCCCAgLMAwyVnKhYigAACCCCAAAIIIIAAAggggAACCCCAAAIIIIAAAggwXKIGEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABBBBAAAEEnAUYLjlTsRABBBBAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQIDhEjWAAAIIIIAAAggggAACCCCAAAIIIIAAAggggAACCCDgLMBwyZmKhQgggAACCCCAAAIIIIAAAggggAACCCCAAAIIIIAAAgyXqAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABBBBAAAFnAYZLzlQsRAABBBBAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQYLhEDSCAAAIIIIAAAggggAACCCCAAAIIIIAAAggggAACCDgLMFxypmIhAggggAACCCCAAAIIIIAAAggggAACCCCAAAIIIIAAwyVqAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABBBBAwFmA4ZIzFQsRQAABBBBAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQYLlEDCCCAAAIIIIAAAggggAACCCCAAAIIIIAAAggggAACzgIMl5ypWIgAAggggAACCCCAAAIIIIAAAggggAACCCCAAAIIIMBwiRpAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABBBBwFmC45EzVrIUPP/xw64F33HHHZj04T4sAAsEC9H8wIRsg0FgB+r+xqePBEYgiwHtAFEY2QaCRAvR/I9PGQyOAAAIIINBoAYZLjU7f4A//L//yL61/OW7cuOxf4b/+67+2nnHMmDHZP2svPSB5aW62ffufXKfJNa64phHovqtv/9fxjE06kz62zxbmYealvwdQH2H1YRlNriy122eV3v+pRKlVvSx22OkFiEQAgVIEGC6VkskBr6NJF0suJHkWIXnJMy8uT+Xb/+TaRdV/Da7+Zi4RuDJccqmTWGuot1iS7vtg7m7VaaXvHSDsNPto6sPeXHsiudLK6eNK73+9TPdIalUvix12egEiEUCgFAGGS6VkkuFSoZms72VxUazPPvRk3x8syXWoeOd4XHFNI8BwydKVPrbUbp+FeZi57x0g7DT7aOrD3lx7IrnSyunjSu9/vQzDJexSCej35T1Sb0ckAgjkJcBwKa98RHuaJl0s+aYaLe1RNyIvUTlNN/Ptf3KdJj244ppGgOGSpSt9bKnNcCmGtu8dIMaZlnvQk5baYWeRqzA/TXTp/a8xcYmhVl2UOq/BDju9AJEIIFCKAMOlUjI54HU06WLJhSTPIiQveebF5al8+59cu6j6r8HV38wlAleGSy51EmsN9RZL0n0fzN2tOq30vQOEnWYfTX3Ym2tPJFdaOX1c6f2vl+keSa3qZbHDTi9AJAIIlCLAcKmUTDJcKjST9b0sLor12Yee7PuDJbkOFe8cjyuuaQQYLlm60seW2u2zMA8z970DhJ1mH0192JtrTyRXWjl9XOn9r5dhuIRdKgH9vrxH6u2IRACBvAQYLuWVj2hP06SLJd9Uo6U96kbkJSqn6Wa+/U+u06QHV1zTCDBcsnSljy21GS7F0Pa9A8Q403IPetJSO+wschXmp4kuvf81Ji4x1KqLUuc12GGnFyASAQRKEWC4VEomB7yOJl0suZDkWYTkJc+8uDyVb/+TaxdV/zW4+pu5RODKcMmlTmKtod5iSbrvg7m7VaeVvneAsNPso6kPe3PtieRKK6ePK73/9TLdI6lVvSx22OkFiEQAgVIEGC6VkkmGS4Vmsr6XxUWxPvvQk31/sCTXoeKd43HFNY0AwyVLV/rYUrt9FuZh5r53gLDT7KOpD3tz7YnkSiunjyu9//UyDJewSyWg35f3SL0dkQggkJcAw6W88hHtaZp0seSbarS0R92IvETlNN3Mt//JdZr04IprGgGGS5au9LGlNsOlGNq+d4AYZ1ruQU9aaoedRa7C/DTRpfe/xsQlhlp1Ueq8Bjvs9AJEIoBAKQIMl0rJ5IDX0aSLJReSPIuQvOSZF5en8u1/cu2i6r8GV38zlwhcGS651EmsNdRbLEn3fTB3t+q00vcOEHaafTT1YW+uPZFcaeX0caX3v16meyS1qpfFDju9AJEIIFCKAMOlUjLJcKnQTNb3srgo1mcferLvD5bkOlS8czyuuKYRYLhk6UofW2q3z8I8zNz3DhB2mn009WFvrj2RXGnl9HGl979ehuESdqkE9PvyHqm3IxIBBPISYLiUVz6iPU2TLpZ8U42W9qgbkZeonKab+fY/uU6THlxxTSPAcMnSlT621Ga4FEPb9w4Q40zLPehJS+2ws8hVmJ8muvT+15i4xFCrLkqd12CHnV6ASAQQKEWA4VIpmRzwOpp0seRCkmcRkpc88+LyVL79T65dVP3X4Opv5hKBK8MllzqJtYZ6iyXpvg/m7ladVvreAcJOs4+mPuzNtSeSK62cPq70/tfLdI+kVvWy2GGnFyASAQRKEWC4VEomGS4Vmsn6XhYXxfrsQ0/2/cGSXIeKd47HFdc0AgyXLF3pY0vt9lmYh5n73gHCTrOPpj7szbUnkiutnD6u9P7XyzBcwi6VgH5f3iP1dkQigEBeAgyX8spHtKdp0sWSb6rR0h51I/ISldN0M9/+J9dp0oMrrmkEGC5ZutLHltoMl2Jo+94BYpxpuQc9aakddha5CvPTRJfe/xoTlxhq1UWp8xrssNMLEIkAAqUIMFwqJZMDXkeTLpZcSPIsQvKSZ15cnsq3/8m1i6r/Glz9zVwicGW45FInsdZQb7Ek3ffB3N2q00rfO0DYafbR1Ie9ufZEcqWV08eV3v96me6R1KpeFjvs9AJEIoBAKQIMl0rJJMOlQjNZ38violiffejJvj9YkutQ8c7xuOKaRoDhkqUrfWyp3T4L8zBz3ztA2Gn20dSHvbn2RHKlldPHld7/ehmGS9ilEtDvy3uk3o5IBBDIS4DhUl75iPY0TbpY8k01WtqjbkReonKabubb/+Q6TXpwxTWNAMMlS1f62FKb4VIMbd87QIwzLfegJy21w84iV2F+mujS+19j4hJDrboodV6DHXZ6ASIRQKAUAYZLpWRywOto0sWSC0meRUhe8syLy1P59j+5dlH1X4Orv5lLBK4Ml1zqJNYa6i2WpPs+mLtbdVrpewcIO80+mvqwN9eeSK60cvq40vtfL9M9klrVy2KHnV6ASAQQKEWA4VIpmWS4VGgm63tZXBTrsw892fcHS3IdKt45Hldc0wgwXLJ0pY8ttdtnYR5m7nsHCDvNPpr6sDfXnkiutHL6uNL7Xy/DcAm7VAL6fXmP1NsRiQACeQkwXMorH9GepkkXS76pRkt71I3IS1RO0818+59cp0kPrrimEWC4ZOlKH1tqM1yKoe17B4hxpuUe9KSldthZ5CrMTxNdev9rTFxiqFUXpc5rsMNOL0AkAgiUIsBwqZRMDngdTbpYciHJswjJS555cXkq3/4n1y6q/mtw9TdzicCV4ZJLncRaQ73FknTfB3N3q04rfe8AYafZR1Mf9ubaE8mVVk4fV3r/62W6R1KrelnssNMLEIkAAqUIMFwqJZMMlwrNZH0vi4tiffahJ/v+YEmuQ8U7x+OKaxoBhkuWrvSxpXb7LMzDzH3vAGGn2UdTH/bm2hPJlVZOH1d6/+tlGC5hl0pAvy/vkXo7IhFAIC8Bhkt55SPa0zTpYsk31Whpj7oReYnKabqZb/+T6zTpwRXXNAIMlyxd6WNLbYZLMbR97wAxzrTcg5601A47i1yF+WmiS+9/jYlLDLXqotR5DXbY6QWIRACBUgQYLpWSyQGvo0kXSy4keRYheckzLy5P5dv/5NpF1X8Nrv5mLhG4MlxyqZNYa6i3WJLu+2DubtVppe8dIOw0+2jqw95ceyK50srp40rvf71M90hqVS+LHXZ6ASIRQKAUAYZLpWSS4VKhmazvZXFRrM8+9GTfHyzJdah453hccU0jwHDJ0pU+ttRun4V5mLnvHSDsNPto6sPeXHsiudLK6eNK73+9DMMl7FIJ6PflPVJvRyQCCOQlwHApr3xEe5omXSz5phot7VE3Ii9ROU038+1/cp0mPbjimkaA4ZKlK31sqc1wKYa27x0gxpmWe9CTltphZ5GrMD9NdOn9rzFxiaFWXZQ6r8EOO70AkQggUIoAw6VSMjngdTTpYsmFJM8iJC955sXlqXz7n1y7qPqvwdXfzCUCV4ZLLnUSaw31FkvSfR/M3a06rfS9A4SdZh9Nfdiba08kV1o5fVzp/a+X6R5JreplscNOL0AkAgiUIsBwqZRMMlwqNJP1vSwuivXZh57s+4MluQ4V7xyPK65pBBguWbrSx5ba7bMwDzP3vQOEnWYfTX3Ym2tPJFdaOX1c6f2vl2G4hF0qAf2+vEfq7YhEAIG8BBgu5ZWPaE/TpIsl31SjpT3qRuQlKqfpZr79T67TpAdXXNMIMFyydKWPLbUZLsXQ9r0DxDjTcg960lI77CxyFeaniS69/zUmLjHUqotS5zXYYacXIBIBBEoRYLhUSiYHvI4mXSy5kORZhOQlz7y4PJVv/5NrF1X/Nbj6m7lE4MpwyaVOYq2h3mJJuu+DubtVp5W+d4Cw0+yjqQ97c+2J5Eorp48rvf/1Mt0jqVW9LHbY6QWIRACBUgQYLpWSSYZLhWayvpfFRbE++9CTfX+wJNeh4p3jccU1jQDDJUtX+thSu30W5mHmvneAsNPso6kPe3PtieRKK6ePK73/9TIMl7BLJaDfl/dIvR2RCCCQlwDDpbzyEe1pmnSx5JtqtLRH3Yi8ROU03cy3/8l1mvTgimsaAYZLlq70saU2w6UY2r53gBhnWu5BT1pqh51FrsL8NNGl97/GxCWGWnVR6rwGO+z0AkQigEApAgyXSsnkgNfRpIslF5I8i5C85JkXl6fy7X9y7aLqvwZXfzOXCFwZLrnUSaw11FssSfd9MHe36rTS9w4Qdpp9NPVhb649kVxp5fRxpfe/XqZ7JLWql8UOO70AkQggUIoAw6VSMslwqdBM1veyuCjWZx96su8PluQ6VLxzPK64phFguGTpSh9barfPwjzM3PcOEHaafTT1YW+uPZFcaeX0caX3v16G4RJ2qQT0+/IeqbcjEgEE8hJguJRXPqI9TZMulnxTjZb2qBuRl6icppv59j+5TpMeXHFNI8BwydKVPrbUZrgUQ9v3DhDjTMs96ElL7bCzyFWYnya69P7XmLjEUKsuSp3XYIedXoBIBBAoRYDhUimZHPA6mnSx5EKSZxGSlzzz4vJUvv1Prl1U/dfg6m/mEoErwyWXOom1hnqLJem+D+buVp1W+t4Bwk6zj6Y+7M21J5IrrZw+rvT+18t0j6RW9bLYYacXIBIBBEoRYLhUSiYZLhWayfpeFhfF+uxDT/b9wZJch4p3jscV1zQCDJcsXeljS+32WZiHmfveAcJOs4+mPuzNtSeSK62cPq70/tfLMFzCLpWAfl/eI/V2RCKAQF4CDJfyyke0p2nSxZJvqtHSHnUj8hKV03Qz3/4n12nSgyuuaQQYLlm60seW2gyXYmj73gFinGm5Bz1pqR12FrkK89NEl97/GhOXGGrVRanzGuyw0wsQiQACpQgwXColkwNeR5MullxI8ixC8pJnXlyeyrf/ybWLqv8aXP3NXCJwZbjkUiex1lBvsSTd98Hc3arTSt87QNhp9tHUh7259kRypZXTx5Xe/3qZ7pHUql4WO+z0AkQigEApAgyXSskkw6VCM1nfy+KiWJ996Mm+P1iS61DxzvG44ppGgOGSpSt9bKndPgvzMHPfO0DYafbR1Ie9ufZEcqWV08eV3v96GYZL2KUS0O/Le6TejkgEEMhLgOFSXvmI9jRNuljyTTVa2qNuRF6icppu5tv/5DpNenDFNY0AwyVLV/rYUpvhUgxt3ztAjDMt96AnLbXDziJXYX6a6NL7X2PiEkOtuih1XoMddnoBIhFAoBQBhkulZHLA62jSxZILSZ5FSF7yzIvLU/n2P7l2UfVfg6u/mUsErgyXXOok1hrqLZak+z6Yu1t1Wul7Bwg7zT6a+rA3155IrrRy+rjS+18v0z2SWtXLYoedXoBIBBAoRYDhUimZZLhUaCbre1lcFOuzDz3Z9wdLch0q3jkeV1zTCDBcsnSljy2122dhHmbuewcIO80+mvqwN9eeSK60cvq40vtfL8NwCbtUAvp9eY/U2xGJAAJ5CTBcyisf0Z6mSRdLvqlGS3vUjchLVE7TzXz7n1ynSQ+uuKYRYLhk6UofW2ozXIqh7XsHiHGm5R70pKV22FnkKsxPE116/2tMXGKoVRelzmuww04vQCQCCJQiwHDJM5OTJ0+WJUuWyP777y/nnHPOGtG33nqr/OAHP5ClS5fKz3/+c9lss83knnvuGfSEBx54QObOnSvLli2T4cOHy+677y4zZ86UUaNGeT7V2subdLHkQhKc7iQbkJckrCab+vY/uU6TFlxxTSPQfVff/q/jGZt0Jn1sny3Mw8xLfw+gPsLqwzKaXFlqt88qvf9TiVKrelnssNMLEIkAAqUIMFzyyOTChQvltNNOk5dffrnjcGnKlCnyk5/8RD7wgQ+0hktDhgwZdLj00EMPyec+9znZbrvt5KCDDpLly5fLZZddJltssYUsWLBAhg0b5vFkDJeCsAjuKMBFsbmF4fuDJblOk2tccU0jwHDJ0pU+ttRun4V5mLnvHSDsNPto6sPeXHsiudLK6eNK73+9TPdIalUvix12egEiEUCgFAGGS46ZXLlypYwfP16mTZsmc+bM6Thc+s///M/Wp5XWW289qQZNv/jFLwYdLk2YMEFWrFght9xyi2ywwQatp1i8eLEceeSRMmvWLDnssMMcn6zzsiZdLLmQBKU6WTB5SUabfGPf/ifXaVKCK65pBLrv6tv/dTxjk86kj+2zhXmYeenvAdRHWH1YRpMrS+32WaX3fypRalUvix12egEiEUCgFAGGS46ZPP300+W+++6TRYsWydixYzsOl966Vbfh0tNPPy377ruvzJgxQ6ZPn77GE+yzzz6y8cYby3XXXef4ZAyXgqAIHlSAi2Jzi8P3B0tynSbXuOKaRoDhkqUrfWyp3T4L8zBz3ztA2Gn20dSHvbn2RHKlldPHld7/epnukdSqXhY77PQCRCKAQCkCDJccMvn444+3fnXdRRddJLvuumvrV9l1+ptLb92q23CpGlAdf/zxcskll8guu+yyxhNU//z222+XRx55pPUJKO1Xky6WXEi0WU4bR17S+qbc3bf/yXWabOCKaxqB7rv69n8dz9ikM+lj+2xhHmZe+nsA9RFWH5bR5MpSu31W6f2fSpRa1ctih51egEgEEChFgOHSOjL5xhtvyCGHHCKjRo1qDZeqr9Dh0qWXXipnn3223HTTTa293vpV/fPq31efkho9erS6zqqLZX9/v4wYMUK9h1XgK6+80jpq+PDhVkdyjoMAeXFASrRkzJgxQTv79j+5DuIeNBhXXDUC1v2vecZeiqGP7bPdy+ah/f/mf1xuys8Amurq5frQeNUZQ6789UPfA3x/BvB/wjIjqFV9XrGLZxfa//onIRIBBBAIE+iZ4VL1Q9aqVauctIYMGSJDhw5trb322mvljDPOaP1tpK233rr1z0KHS/Pnz5fzzjtPbrvtNtlmm23WeKa5c+fKBRdcIHfeeadsueWWTs/baVGTLpZcSNRpThpIXpLydt089GLp2//kOk2uccVVI2Dd/5pn7KUY+tg+271sHtr/VbZ87wD2GQ47sZfrI0zOPppc+ZuHvgeU3v/+om4R1KqbU6dV2MWzC+1//ZMQiQACCIQJ9Mxw6Yknnmj9KjuXr4kTJ8rs2bNl+fLlMn78eJk0aZIcc8wxq0NDh0tWn1yqHnjcuHEuL7nWNXyUulb+QQ8nL3nmxeWpfH8lBrl2UfVfg6u/mUsErt2VfPvfxbyX11Bv9tnHPMy89PcA6iOsPiyjyZWldvus0vs/lSi1qpfFDju9AJEIIFCKQM8Ml1auXCl33HGHU96qTyjttNNOrU8sVX8f6eqrr5Zhw4atjt1zzz1lr732klmzZskmm2wiI0eOXGvfkL+5VH2i6dFHH+VvLjlli0WpBLgoppJNv6/vD5bkOk1OcMU1jQDDJUtX+thSu30W5mHmvneAsNPso6kPe3PtieRKK6ePK73/9TLdI6lVvSx22OkFiEQAgVIEema4pEnYUUcd1fr1dN2+Tj31VJk8efJaS7oNl5566qnWJ6JmzJgh06dPXyN2n332kY022kgWLFigeeTVMU26WHIhCUp1smDykow2+ca+/U+u06QEV1zTCHTf1bf/63jGJp1JH9tnC/Mw89LfA6iPsPqwjCZXltrts0rv/1Si1KpeFjvs9AJEIoBAKQIMl7pksvr00PPPP7/WiqOPPlp23nlnmTZtWuvvL2211VZew6Vq8QEHHCDVp6mqv+W0wQYbtOIXL14sRx55pJx44oly+OGHB9VYky6WXEiCUp0smLwko02+sW//k+s0KcEV1zQC3Xf17f86nrFJZ9LH9tnCPMy89PcA6iOsPiyjyZWldvus0vs/lSi1qpfFDju9AJEIIFCKAMMlRSYH+5tLS5Yskep/qq/qk0crVqyQz3/+863/fYsttpAJEyasPu3BBx9sDZCqP9p38MEHywsvvCCXX365bL755nL99dfL8OHDFU/2vyFNulhyIQlKdbJg8pKMNvnGvv1PrtOkBFdc0wh039W3/+t4xiadSR/bZwvzMPPS3wOoj7D6sIwmV5ba7bNK7/9UotSqXhY77PQCRCKAQCkCDJcUmRxsuHT++efLvHnzOu5YfdLpqquuWuPf3X///TJ37lxZtmxZa5i02267ycyZM2X06NGKp1ozpEkXSy4kwelOsgF5ScJqsqlv/5PrNGnBFdc0AgyXLF3pY0vt9lmYh5n73gHCTrOPpj7szbUnkiutnD6u9P7Xy3SPpFb1sthhpxcgEgEEShFguFRKJge8jiZdLLmQ5FmE5CXPvLg8lW//k2sXVf81uPqbuUTgynDJpU5iraHeYkm674O5u1Wnlb53gLDT7KOpD3tz7YnkSiunjyu9//UyDJewSyWg35f3SL0dkQggkJcAw6W88hHtaZp0seSbarS0R92IvETlNN3Mt//JdZr04IprGgGGS5au9LGldvsszMPMfe8AYafZR1Mf9ubaE8mVVk4fV3r/62UYLmGXSkC/L++RejsiEUAgLwGGS3nlI9rTNOliyTfVaGmPuhF5icppuplv/5PrNOnBFdc0AgyXLF3pY0tthksxtH3vADHOtNyDnrTUDjuLXIX5aaJL73+NiUsMteqi1HkNdtjpBYhEAIFSBBgulZLJAa+jSRdLLiR5FiF5yTMvLk/l2//k2kXVfw2u/mYuEbgyXHKpk1hrqLdYku77YO5u1Wml7x0g7DT7aOrD3lx7IrnSyunjSu9/vUz3SGpVL4sddnoBIhFAoBQBhkulZJLhUqGZrO9lcVGszz70ZN8fLMl1qHjneFxxTSPAcMnSlT621G6fhXmYue8dIOw0+2jqw95ceyK50srp40rvf70MwyXsUgno9+U9Um9HJAII5CXAcCmvfER7miZdLPmmGi3tUTciL1E5TTfz7X9ynSY9uOKaRoDhkqUrfWypzXAphrbvHSDGmZZ70JOW2mFnkaswP0106f2vMXGJoVZdlDqvwQ47vQCRCCBQigDDpVIyOeB1NOliyYUkzyIkL3nmxeWpfPufXLuo+q/B1d/MJQJXhksudRJrDfUWS9J9H8zdrTqt9L0DhJ1mH0192JtrTyRXWjl9XOn9r5fpHkmt6mWxw04vQCQCCJQiwHCplEwyXCo0k/W9LC6K9dmHnuz7gyW5DhXvHI8rrmkEGC5ZutLHltrtszAPM/e9A4SdZh9Nfdiba08kV1o5fVzp/a+XYbiEXSoB/b68R+rtiEQAgbwEGC7llY9oT9OkiyXfVKOlPepG5CUqp+lmvv1PrtOkB1dc0wgwXLJ0pY8ttRkuxdD2vQPEONNyD3rSUjvsLHIV5qeJLr3/NSYuMdSqi1LnNdhhpxcgEgEEShFguFRKJge8jiZdLLmQ5FmE5CXPvLg8lW//k2sXVf81uPqbuUTgynDJpU5iraHeYkm674O5u1Wnlb53gLDT7KOu9qskAAAgAElEQVSpD3tz7YnkSiunjyu9//Uy3SOpVb0sdtjpBYhEAIFSBBgulZJJhkuFZrK+l8VFsT770JN9f7Ak16HineNxxTWNAMMlS1f62FK7fRbmYea+d4Cw0+yjqQ97c+2J5Eorp48rvf/1MgyXsEsloN+X90i9HZEIIJCXAMOlvPIR7WmadLHkm2q0tEfdiLxE5TTdzLf/yXWa9OCKaxoBhkuWrvSxpTbDpRjavneAGGda7kFPWmqHnUWuwvw00aX3v8bEJYZadVHqvAY77PQCRCKAQCkCDJdKyeSA19GkiyUXkjyLkLzkmReXp/Ltf3Ltouq/Bld/M5cIXBkuudRJrDXUWyxJ930wd7fqtNL3DhB2mn009WFvrj2RXGnl9HGl979epnsktaqXxQ47vQCRCCBQigDDpVIyyXCp0EzW97K4KNZnH3qy7w+W5DpUvHM8rrimEWC4ZOlKH1tqt8/CPMzc9w4Qdpp9NPVhb649kVxp5fRxpfe/XobhEnapBPT78h6ptyMSAQTyEmC4lFc+oj1Nky6WfFONlvaoG5GXqJymm/n2P7lOkx5ccU0jwHDJ0pU+ttRmuBRD2/cOEONMyz3oSUvtsLPIVZifJrr0/teYuMRQqy5Knddgh51egEgEEChFgOFSKZkc8DqadLHkQpJnEZKXPPPi8lS+/U+uXVT91+Dqb+YSgSvDJZc6ibWGeosl6b4P5u5WnVb63gHCTrOPpj7szbUnkiutnD6u9P7Xy3SPpFb1sthhpxcgEgEEShFguFRKJhkuFZrJ+l4WF8X67ENP9v3BklyHineOxxXXNAIMlyxd6WNL7fZZmIeZ+94Bwk6zj6Y+7M21J5IrrZw+rvT+18swXMIulYB+X94j9XZEIoBAXgIMl/LKR7SnadLFkm+q0dIedSPyEpXTdDPf/ifXadKDK65pBBguWbrSx5baDJdiaPveAWKcabkHPWmpHXYWuQrz00SX3v8aE5cYatVFqfMa7LDTCxCJAAKlCDBcKiWTA15Hky6WXEjyLELykmdeXJ7Kt//JtYuq/xpc/c1cInBluORSJ7HWUG+xJN33wdzdqtNK3ztA2Gn20dSHvbn2RHKlldPHld7/epnukdSqXhY77PQCRCKAQCkCDJdKySTDpUIzWd/L4qJYn33oyb4/WJLrUPHO8bjimkaA4ZKlK31sqd0+C/Mwc987QNhp9tHUh7259kRypZXTx5Xe/3oZhkvYpRLQ78t7pN6OSAQQyEuA4VJe+Yj2NE26WPJNNVrao25EXqJymm7m2//kOk16cMU1jQDDJUtX+thSm+FSDG3fO0CMMy33oCcttcPOIldhfpro0vtfY+ISQ626KHVegx12egEiEUCgFAGGS6VkcsDraNLFkgtJnkVIXvLMi8tT+fY/uXZR9V+Dq7+ZSwSuDJdc6iTWGuotlqT7Ppi7W3Va6XsHCDvNPpr6sDfXnkiutHL6uNL7Xy/TPZJa1ctih51egEgEEChFgOFSKZlkuFRoJut7WVwU67MPPdn3B0tyHSreOR5XXNMIMFyydKWPLbXbZ2EeZu57Bwg7zT6a+rA3155IrrRy+rjS+18vw3AJu1QC+n15j9TbEYkAAnkJMFzKKx/RnqZJF0u+qUZLe9SNyEtUTtPNfPufXKdJD664phFguGTpSh9bajNciqHteweIcablHvSkpXbYWeQqzE8TXXr/a0xcYqhVF6XOa7DDTi9AJAIIlCLAcKmUTA54HU26WHIhybMIyUueeXF5Kt/+J9cuqv5rcPU3c4nAleGSS53EWkO9xZJ03wdzd6tOK33vAGGn2UdTH/bm2hPJlVZOH1d6/+tlukdSq3pZ7LDTCxCJAAKlCGQxXHrttdfkxRdflE022UTWW2+91bb33HOP3H333bL++uvLZz/7Wdl2221LcU/+Opp0seRCkrwcVAeQFxVbFkG+/U+u06QNV1zTCDBcsnSljy2122dhHmbuewcIO80+mvqwN9eeSK60cvq40vtfL8NwCbtUAvp9eY/U2xGJAAJ5CWQxXPra174mN954o9x7770yYsSIltB1110nX/nKV6S/v7/1v1f/fMGCBbLNNtvkJZjp0zTpYsk31TyLiLzkmReXp/Ltf3Ltouq/Bld/M5cIXLsr+fa/i3kvr6He7LOPeZh56e8B1EdYfVhGkytL7fZZpfd/KlFqVS+LHXZ6ASIRQKAUgSyGS+PHj28NjS644ILVrrvvvrv09fXJN77xDXnhhRfkhBNOkGrdmWeeWYp90tfRpIslF5KkpaDenLyo6WoP9O1/cp0mZbjimkaA4ZKlK31sqd0+C/Mwc987QNhp9tHUh7259kRypZXTx5Xe/3qZ7pHUql4WO+z0AkQigEApAlkMl3baaSc5+OCD5cQTT2y5PvHEE/LpT39aZs2aJdOmTWv9s+OPP14eeeQRueOOO0qxT/o6mnSx5EKStBTUm5MXNV3tgb79T67TpAxXXNMIdN/Vt//reMYmnUkf22cL8zDz0t8DqI+w+rCMJleW2u2zSu//VKLUql4WO+z0AkQigEApAlkMl3bccUeZNGlS69NJ1ddVV10lX//61+Xmm29e/XeWzj33XLniiivk0UcfLcU+6eto0sWSC0nSUlBvTl7UdLUH+vY/uU6TMlxxTSPAcMnSlT621G6fhXmYue8dIOw0+2jqw95ceyK50srp40rvf71M90hqVS+LHXZ6ASIRQKAUgSyGS/vtt5+84x3vaA2Vqq+pU6fKL37xC1m8ePFq55NPPlnuueee1t9l4mvdAk26WHIhWXc+61hBXupQj3Omb/+T6zjuA3fBFdc0At139e3/Op6xSWfSx/bZwjzMvPT3AOojrD4so8mVpXb7rNL7P5UotaqXxQ47vQCRCCBQikAWw6ULL7xQ5s6dK/vss48MGzZMFi1aJIcffrjMnDlztfNBBx0kb3/72+Xqq68uxT7p62jSxZILSdJSUG9OXtR0tQf69j+5TpMyXHFNI8BwydKVPrbUbp+FeZi57x0g7DT7aOrD3lx7IrnSyunjSu9/vUz3SGpVL4sddnoBIhFAoBSBLIZLr776autvKt15553S398vH//4x+X888+XDTbYoOX85JNPyqc+9SmZPn1663/4WrdAky6WXEjWnc86VpCXOtTjnOnb/+Q6jvvAXXDFNY1A9119+7+OZ2zSmfSxfbYwDzMv/T2A+girD8tocmWp3T6r9P5PJUqt6mWxw04vQCQCCJQikMVw6U3M3/3ud9LX1ycjR45cw3f58uXy3HPPybvf/W7ZcMMNS7FP+jqadLHkQpK0FNSbkxc1Xe2Bvv1PrtOkDFdc0wgwXLJ0pY8ttdtnYR5m7nsHCDvNPpr6sDfXnkiutHL6uNL7Xy/TPZJa1ctih51egEgEEChFIIvh0qxZs2TMmDEybdq0Ulxrfx1NulhyIam9XDo+AHnJMy8uT+Xb/+TaRdV/Da7+Zi4RuDJccqmTWGuot1iS7vtg7m7VaaXvHSDsNPto6sPeXHsiudLK6eNK73+9DMMl7FIJ6PflPVJvRyQCCOQlkMVw6UMf+pBMnTpVjjvuuLx0Gvw0TbpY8k01z0IjL3nmxeWpfPufXLuo+q/B1d/MJQJXhksudRJrDfUWS9J9H8zdrRguhVkRnVaAXk7r24v9n0qUWtXLYoedXoBIBBAoRSCL4dKBBx4o22yzjcyZM6cU19pfh+9/XK7zgbmQ1Kk/+NnkJc+8uDyVb/+TaxdV/zW4+pu5RODKcMmlTmKtod5iSbrvg7m7VS/+x2XqI6w+LKPJlaV2+yzfnwHsnzDPE6lVfV6ww04vQCQCCJQikMVw6dZbb5XqV+NdddVV8sEPfrAU21pfR5MullxIai2VQQ8nL3nmxeWpfPufXLuo+q/B1d/MJQJXhksudRJrDfUWS9J9H8zdrRguhVkRnVaAXk7r24v9n0qUWtXLYoedXoBIBBAoRSCL4dLChQtl0aJF8qMf/Ug+8YlPyPbbby+bbrqp9PX1reU8YcKEUuyTvg7f/7ic9GHWsTkXkjr1Bz+bvOSZF5en8u1/cu2i6r8GV38zlwhcuyv59r+LeS+vod7ss495mHnp7wHUR1h9WEaTK0vt9lml938qUWpVL4sddnoBIhFAoBSBLIZLY8aMaQ2S+vv713B963Cp+nfV/75s2bJS7JO+jiZdLLmQJC0F9ebkRU1Xe6Bv/5PrNCnDFdc0AgyXLF3pY0vt9lmYh5n73gHCTrOPpj7szbUnkiutnD6u9P7Xy3SPpFb1sthhpxcgEgEEShHIYrh0ww03OHtOnDjReW0vL2zSxZILSZ6VSl7yzIvLU/n2P7l2UfVfg6u/mUsErgyXXOok1hrqLZak+z6Yu1t1Wul7Bwg7zT6a+rA3155IrrRy+rjS+18vw3AJu1QC+n15j9TbEYkAAnkJZDFcyoukjKdp0sWSb6p51hx5yTMvLk/l2//k2kXVfw2u/mYuEbgyXHKpk1hrqLdYku77YO5uxXApzIrotAL0clrfXuz/VKLUql4WO+z0AkQigEApAgyXSsnkgNfh+x+X62TgQlKn/uBnk5c88+LyVL79T65dVP3X4Opv5hKBK8MllzqJtYZ6iyXpvg/m7la9+B+XqY+w+rCMJleW2u2zfH8GsH/CPE+kVvV5wQ47vQCRCCBQikBWw6W7775bbr75ZnniiSfkxRdflJEjR8r73vc+2X///WW33XYrxdzkdTTpYsmFxKQkvA8hL95k2QT49j+5TpM6XHFNI8BwydKVPrbUbp+FeZi57x0g7DT7aOrD3lx7IrnSyunjSu9/vUz3SGpVL4sddnoBIhFAoBSBLIZLq1atkmOPPVbuuusu6e/vl7e97W2yySabyG9/+1t5/fXXpa+vT/bYYw/55je/Keuvv34p9klfR5MullxIkpaCenPyoqarPdC3/8l1mpThimsaAYZLlq70saU2w6UY2r53gBhnWu5BT1pqh51FrsL8NNGl97/GxCWGWnVR6rwGO+z0AkQigEApAlkMl77xjW/IpZdeKrvssovMmDFDxo4d2xooVYOmxx57TM4//3y599575YgjjpDjjjuuFPukr6NJF0suJElLQb05eVHT1R7o2//kOk3KcMU1jQDDJUtX+thSm+FSDG3fO0CMMy33oCcttcPOIldhfpro0vtfY+ISQ626KDFc0ithF9uO/RBAIC+BLIZL1VDpne98p3zve9/rqFMNmT7zmc/I888/Lz/84Q/zEsz0aZp0seQyl2cRkZc88+LyVL79T65dVP3X4Opv5hKBK8MllzqJtYZ6iyXpvg/m7ladVvreAcJOs4+mPuzNtSeSK62cPq70/tfLdI+kVvWy2GGnFyASAQRKEchiuLTDDjvI1KlT5a//+q8HdT333HPlyiuvlEceeaQU+6Svo0kXSy4kSUtBvTl5UdPVHujb/+Q6TcpwxTWNAMMlS1f62FK7fRbmYea+d4Cw0+yjqQ97c+2J5Eorp48rvf/1MgyXsEsloN+X90i9HZEIIJCXQBbDpc9+9rOy9dZbyznnnDOoTvXr8J599lm59tpr8xLM9GmadLHkm2qeRURe8syLy1P59j+5dlH1X4Orv5lLBK7dlXz738W8l9dQb/bZxzzMvPT3AOojrD4so8mVpXb7rNL7P5UotaqXxQ47vQCRCCBQikAWw6X7779fvvjFL8rs2bPlk5/85Fq2t9xyi5x88sly0UUXycc+9rFS7JO+jiZdLLmQJC0F9ebkRU1Xe6Bv/5PrNCnDFdc0AgyXLF3pY0vt9lmYh5n73gHCTrOPpj7szbUnkiutnD6u9P7Xy3SPpFb1sthhpxcgEgEEShHIYrg0b948efjhh6UaMr33ve+V6tfkjRo1SpYvX976NXj//u//Ln/8x38sO+644xrufX19cvTRR5vmYvLkybJkyRLZf//91/ik1W9/+1u54YYb5O6775Ynn3xSXn75Zdlqq61kv/32k2nTpsnb3/72tZ7zgQcekLlz58qyZctk+PDhsvvuu8vMmTNbrz30q0kXSy4kodlOE09e0rha7Orb/+Q6TVZwxTWNQPddffu/jmds0pn0sX22MA8zL/09gPoIqw/LaHJlqd0+q/T+TyVKreplscNOL0AkAgiUIpDFcGnMmDEqz2q4VA1mrL4WLlwop512WmtwNHC4VA2VqkHXxz/+cfnoRz8qI0eObA2hbr75Zhk3bpxcddVVst56661+1Iceekg+97nPyXbbbScHHXRQa5B22WWXyRZbbCELFiyQYcOGBb2sJl0suZAEpTpZMHlJRpt8Y9/+J9dpUoIrrmkEGC5ZutLHltrtszAPM/e9A4SdZh9Nfdiba08kV1o5fVzp/a+X6R5JreplscNOL0AkAgiUIpDFcKkatGi/dt55Z22oV9zKlStl/PjxrU8hzZkzZ63h0jPPPNPar/q00lu/qk8mXXDBBVJ9Omvvvfde/a8mTJggK1askOpX/m2wwQatf7548WI58sgjZdasWXLYYYd5Pd/AxU26WHIhCUp1smDykow2+ca+/U+u06QEV1zTCHTf1bf/63jGJp1JH9tnC/Mw89LfA6iPsPqwjCZXltrts0rv/1Si1KpeFjvs9AJEIoBAKQJZDJe0mC+++KJUQ5/q0z6pv04//XS57777ZNGiRTJ27Ni1hkuDnf+zn/1MPv3pT8uXvvQlOeqoo1rLnn76adl3331lxowZMn369DVC99lnH9l4443luuuuC3pJTbpYciEJSnWyYPKSjDb5xr79T67TpARXXNMIdN/Vt//reMYmnUkf22cL8zDz0t8DqI+w+rCMJleW2u2zSu//VKLUql4WO+z0AkQigEApAo0eLlWfBpo/f37yX433+OOPt3513UUXXSS77rpr61fZDfy1eIMVxA9/+EM54ogjWr9O75BDDmktqwZUxx9/vFxyySWyyy67rBFa/fPbb7+99bem3vpr9HwLrrpY9vf3y4gRI3xDzde/8sorrTOrvzvFVz4C5KW+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Ti9AJAIIlCLAcKmUTDJcKjST9b0sLor12Yee7PuDJbkOFe8cjyuuaQQYLlm60seW2u2zMA8z970DhJ1mH0192JtrTyRXWjl9XOn9r5dhuIRdKgH9vrxH6u2IRACBvAQYLuWVj2hP06SLJd9Uo6U96kbkJSqn6Wa+/U+u06QHV1zTCDBcsnSljy21GS7F0Pa9A8Q403IPetJSO+wschXmp4kuvf81Ji4x1KqLUuc12GGnFyASAQRKEWC45JnJyZMny5IlS2T//feXc845Z43os846q/XvnnnmGXn55Zdl8803l5133lmOOuoo2XLLLdc66YEHHpC5c+fKsmXLZPjw4bL77rvLzJkzZdSoUZ5PtfbyJl0suZAEpzvJBuQlCavJpr79T67TpAVXXNMIdN/Vt//reMYmnUkf22cL8zDz0t8DqI+w+rCMJleW2u2zSu//VKLUql4WO+z0AkQigEApAgyXPDK5cOFCOe2001qDo07DpUMPPVS222472XrrrWXEiBGtIdOCBQvk9ddfl+uvv1622mqr1ac99NBD8rnPfa61/qCDDpLly5fLZZddJltssUUrZtiwYR5PxnApCIvgjgJcFJtbGL4/WJLrNLnGFdc0AgyXLF3pY0vt9lmYh5n73gHCTrOPpj7szbUnkiutnD6u9P7Xy3SPpFb1sthhpxcgEgEEShFo9HDp29/+tlx55ZVy1113Jc/HypUrZfz48TJt2jSZM2dOx+FSp4d47LHHWsOjI488Uo477rjVSyZMmCArVqyQW265RTbYYIPWP1+8eHFr3axZs+Swww4Lek1NulhyIQlKdbJg8pKMNvnGvv1PrtOkBFdc0wh039W3/+t4xiadSR/bZwvzMPPS3wOoj7D6sIwmV5ba7bNK7/9UotSqXhY77PQCRCKAQCkCWQyXpk6dKgceeKBUA5fBvm688cbWp3+qYVIdX6effrrcd999smjRIhk7dqzzcKn6RNLHPvYxOeSQQ1qfeqq+nn76adl3331lxowZMn369DVezj777CMbb7yxXHfddUEvs0kXSy4kQalOFkxektEm39i3/8l1mpTgimsage67+vZ/Hc/YpDPpY/tsYR5mXvp7APURVh+W0eTKUrt9Vun9n0qUWtXLYoedXoBIBBAoRSCL4dKYMWNaQ5aBg5a3Il944YVy3nnntf4+kfXX448/3vr00UUX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bdo0mTNnTtfh0jPPPCP77befHH300a21A4dL1aeadt99d9liiy3ku9/9rqy33nqtp/jOd77TGl7NmzdP9t57b8cn67ysSRdLLiRBqU4WTF6S0Sbf2Lf/yXWalOCKaxqB7rv69n8dz9ikM+lj+2xhHmZe+nsA9RFWH5bR5MpSu31W6f2fSpRa1ctih51egEgEEChFgOGSYyZPP/10ue+++2TRokUyduzYrsOlL3zhC/LSSy/J7NmzZc8991xruHTvvffK5z//eTnrrLNkwoQJq5+gGjp95CMfkV133VW+9a1vOT4Zw6UgKIIHFeCi2Nzi8P3BklynyTWuuKYRYLhk6UofW2q3z8I8zNz3DhB2mn009WFvrj2RXGnl9HGl979epnsktaqXxQ47vQCRCCBQigDDJYdMPv7443LQQQfJRRdd1Br8bLfddoMOl+644w75q7/6K7nhhhtafzup03Cp2ueb3/ymfP/735f3vOc9azzBpEmT5L/+67/kn/7pnxyebPAl1cWyv78/yt9vCnoQh+BXXnmltWr48OEOq1liJUBerKTXPmfMmDFBh/v2P7kO4h40GFdcNQLW/a95xl6KoY/ts93L5qH9X2XL9w5gn+GwE3u5PsLk7KPJlb956HtA6f3vL+oWQa26OXVahV08u9D+1z8JkQgggECYAMOldfi98cYbcsghh8ioUaNaw6Xqa7DhUvWNtfr7SXvssYeccsop8uyzz3YcLlWfgrr66qvlxz/+sYwcOXKNJ/jSl74kd999tyxdujQos026WHIhCUp1smDykox2nRuHXix9+59crzMlqgW4qtjWGVS6q3X/rxO8xxeUXm85preXzUP7v8qn7x0gxxro9ky9XB/kqmkC/s8b+h5Qev/7i7pF8L7i5tRpFXbx7EL7X/8kRCKAAAJhAj0zXKo+xVP92jmXryFDhsjQoUNbS6+99lo544wz5JZbbpGtt9669c8GGy6de+65ct1118ntt98uG2200aDDpZNPPlmuv/56eeyxx2T99ddf45FOOOEEufHGG1u/EqSvr+//s3cvYFdV1eL/B6go4s8bHC3RPGal5RXM2zFTS9AML3iCPCKQHh8tAcH+KmGkede8FIieMLIUC/uJoKmpCXG8Zcrh4C04WuSPn+W/QA09pBy88H/W2n/euOx3M+eYa4611tzf/Tw9dWCNNdf6jDEW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EFeNq2//a9Zopxiej/bZxjzMPPVnAPURVh+W0eTKUpvhUog2tarXww47vQCRCCCQikBbDJfee+89mTx5stx1112ydOlS6d27t5xyyikyZMgQ6dKlS6e5XLZsmcycOVPmzJkjixYtkrffflt22mknGTBggAwfPlw23XTT9WLff/99mTZtmkyfPl1efvll6datm+y6664yYsQIOfTQQ9c6/sknn5QJEybIwoULpXv37nLEEUfIeeedJ9tuu21wfdXpB0s2JMHpjnIC8hKF1eSkvv1PruOkBVdc4wi0Pqtv/5dxjXVakz62zxbmYeapPwOoj7D6sIwmV5baDJdCtKlVvR522OkFiEQAgVQE2mK4NH78eLnzzjtl8ODBsvfee8vjjz8uDz74oIwaNUpGjhzZaS6zoVI2FDrkkEPkoIMOki222ELmzp0r9913n/Tt21emTp0qG220UUf8Bx98IGeffbY88sgjMnDgQNlrr73knXfekd///vf5/x40aFDHsU8//bSceuqpsttuu8mXvvQleeONN+SWW26RHXbYIR9MbbbZZkE1VqcfLNmQBKU6WjB5iUYb/cS+/U+u46QEV1zjCDBcsnSljy21G2thHmbuuwcIW80+mvqwN9euSK60cvq41PtfL9M6klrVy2KHnV6ASAQQSEUg+eFS9q6gE044QU477TQZO3ZsR97GjBkjs2fPzv+z3XbbNc3nK6+8kv959m6lNV/Zu41uuukmmTRpkvTr16/jr2677Ta5+uqr5dZbb5VPf/rTLWsku6Y333xT7r//ftl8883zY7Oh1BlnnCHjxo2Tr3zlK0E1VqeNJRuSoFRHCyYv0Wijn9i3/8l1nJTgimscgdZn9e3/Mq6xTmvSx/bZwjzMPPVnAPURVh+W0eTKUruxVur9H0uUWtXLYoedXoBIBBBIRSD54dL111+ffyRe9i6k7F1Bq1/z5s2Tk08+WS666KL8v31eL774ohx33HEyevRoOeuss/LQ7F1Ln//85/N3KE2cODH/v7N3LfXo0WO9U2cfl3f00Uc3fedU//79ZauttsrfaRXyqtPGkg1JSKbjxZKXeLaxz+zb/+Q6TkZwxTWOQOuz+vZ/GddYpzXpY/tsYR5mnvozgPoIqw/LaHJlqd1YK/X+jyVKreplscNOL0AkAgikIpD8cCl7x9JLL72UfxTemq+VK1fKPvvs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3j4eV1zTCDBcsnSljy21W3thHmbuewYI280+mvqwN9fuSK60cvq43PtfL8NwCbtUAvp1eY3U2xGJAAL1EmC4VK98RLubJh0s+aEaLe1RFyIvUTlNF/Ptf3KdJj244ppGgOGSpSt9bKnNcCmGtu8ZIMaelmvQk5baYXuRqzA/TXTu/a8xcYmhVl2U2l+DHXZ6ASIRQCAXgR4xXHrnnXdkwYIFcuutt8qqVatkyJAhMnHiRJkwYYL06tWr21y+/vrrcvvtt8v9998vzz77rLz11luy4447yuGHHy5TpkyRzTffvNvY1atXy5gxY+T555+Xk08+WWbOnLnRtXfffbcsXLiwXHubbbaRz33uczJjxgzZaqutguurSQdLDiTB6U6yAHlJwmqyqG//k+s0acEV1zQCDJcsXeljS22GSzG0fc8AMfa0XIOetNQO24tchflponPvf42JSwy16qLEcEmvhF1sO9ZDAIF6CfSI4dI555wjixcvlnHjxsmuu+4qDz30kNxzzz0ybdo0mTp1arcZKYZKp512muy7774ycuRIGTBggDz22GNy5513yogRI2TRokXSp0+ftvHz58+Xq6++uhxItRsu/eQnP5EzzzyzXPewww4rh1DXXXed7LnnnnLNNdd0HHq5lFCTDpYc5lwyan8NebE3j7Wjb/+T61jyG66DK65pBDqv6tv/Vdxjk/akj+2zhXmYee6vAdRHWH1YRpMrS+3WXrn3fypRalUvix12egEiEUAgF4Hsh0srVqyQsWPHygknnCCzZ89en7fiHULLli0r/9l+++3b5vOFF14o/7x4t9L7H1dccYVceeWVMm/ePDnooIM2ii3iinc3FYOpuXPnbjRcKt7VtP/++8vgwYPl5ptvXj+guvHGG+X888/vdl2fomvSwZIDiU9m7a4lL3bWsXfy7X9yHTsDrfVwxTWNQOdVffu/ints0p70sX22MA8zz/01gPoIqw/LaHJlqd3aK/f+TyVKreplscNOL0AkAgjkIpD9cOmyyy4rPxKveBdSMcxZ91i+fLmMHz9ezjvvvPJ/fR5PP/20HHHEETJ9+nQ59dRTNwo96aST5M0335Q5c+bIAQccsNFwqXjn1Je+9CW55JJLysHXukcxdNp7771l1KhR8p3vfMfnlja6tkkHSw4kQalOFkxektEmX9i3/8l1mpTgimsagc6r+vZ/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NVnUt//JdZq04IprGoHOq/r2fxX32KQ96WP7bGEeZp77awD1EVYfltHkylK7tVfu/Z9KlFrVy2KHnV6ASAQQyEUg++FS8fF0/fr1k9tuu22jnO2zzz6y8847l9+N5Pp47733ZMqUKfLkk09KMUgaOHDg+tC3335bDjvsMBk9erQU3/P04osvdjtcKt7dNHz4cPnoRz9avsvpgQcekCVLlpR/VnyXU/F9TCGP4mBZDK6KAVfdH4Vb8ejfv3/db7VH3R95qS7dw4YNC9rct//JdRB3t8G44qoRsO5/zT32pBj62D7bPdk8tP/X/cflpvwOoKmunlwfGq8qY8iVv37oa4Dv7wD+d5hnBLWqzyt28exC+19/J0QigAACYQLZD5cOPPDAcgB08803byS13377ld+nVAxzXB/rPmavGB5NmjRpg7Di7xYvXiz33nuvbL311h2HS+32u/zyy+Wqq66Sb33rW+U7qkIeTTpYciAJyXS6WPKSznZTK4ceLH37n1xvKiO6v8dV57apqNxdrft/U949/e9zr7c65rcnm4f2f5FP3zNAHWug0z315PogV00T8L/f0NeA3PvfX9QtgtcVN6d2V2EXzy60//V3QiQCCCAQJpD9cCnmO5euv/56ufDCC8uPwis+Eu/9j+eee678HqbiO5yOOeaY8q86vXOpXdreeOMN2WOPPaS457lz5wZltklvieet1EGpThZMXpLRJl/Yt//JdZqU4IprGoHOq/r2fxX32KQ96WP7bGEeZp77awD1EVYfltHkylK7tVfu/Z9KlFrVy2KHnV6ASAQQyEUg++HSpr5z6aijjpKLL754k/ksPlbv7LPPlkMPPVQuvfRS6d279wYxp5xyihQDpuIj9nr16lX+3csvvywTJkyQiRMnyvHHH1++g2qLLbbouNcnP/lJ+cQnPiE/+MEPNnlPnS5o0sGSA0lQqpMFk5dktMkX9u1/cp0mJbjimkag86q+/V/FPTZpT/rYPluYh5nn/hpAfYTVh2U0ubLUbu2Ve/+nEqVW9bLYYacXIBIBBHIRyH64VLwDaOHChXL//ffL4MGD1+dt+fLlMn78eDn33HPLAVCnx1133SVnnHGGjBo1SubNmyebbbbZRpcfeeSRsu4Ha3drFR95t//++3e71euvvy577723jBkzphxghTyadLDkQBKS6XSx5CWdbeqVffuf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Jl/Yt//JdZqU4IprGoHOq/r2fxX32KQ96WP7bGEeZp77awD1EVYfltHkylK7tVfu/Z9KlFrVy2KHnV6ASAQQyEUg++FS8Y6lZ555pvwovPc/ikHObrvtJkcffbRcdNFFXvl88MEH5cQTTyzfZXTsscduELt06VI5/fTTy+9qKr47qd1w6aqrrpLLL79clixZIh/5yEc2iD/uuOPk1Vdflfvuu8/rnrpe3KSDJQeSoFQnCyYvyWiTL+zb/+Q6TUpwxTWNQOdVffu/ints0p70sX22MA8zz/01gPoIqw/LaHJlqd3aK/f+TyVKreplscNOL0AkAgjkIlCb4dL//M//yO9+97vyO4re/26eUOji4+n69esnt91220ZL7bPPPrLzzjuX343k+njvvfdkypQp8uSTT0oxSBo4cOD60Lfffrv8/qTi4/+K73l68cUX2w6XLrjgArnhhhukePdU8T1O738U74Yq3mX1xBNPuN5S2+uKg+XatWvLAVfdH4Vb8ejfv3/db7VH3R95qS7dw4YNC9rct//JdRB3t8G44qoRsO5/zT32pBj62D7bPdk8tP+LbPmeAewzHLZjT66PMDn7aHLlbx76GpB7//uLukVQq25O7a7CLp5daP/r74RIBBBAIEyg8uHSr371K/nGN74hzz777PpnUnxPUvEo/q74+Ljie4sOPPBA1TMt4ooBUPHdRl0f++23X/l9SosWLXJee93H7BXDo0mTJm0QV/zd4sWL5d5775Wtt9662+HS2WefLbfeems5oCoGX+9/zJo1S3784x+X39XRq1cv5/vqemGTDpYcSNRpThpIXpLydlw89GDp2//kOk2uccVVI2Dd/5p77Ekx9LF9tnuyeWj/F9nyPQPYZzhsx55cH2Fy9tHkyt889DUg9/73F3WLoFbdnNpdhV08u9D+198JkQgggECYQKXDpeLdORMnTpQttthCvvCFL5QDpn//93+XdcOl4qkVw6Fddtml/Bg5zSPmO5euv/56ufDCC8uPwis+Eu/9j+eee678HqbiO5yOOeaY8q+qfudScQ8jRozQsJnG8FZqU27nzciLM1XtLvT9SAxynSaFuOKaRqDzqr79X8U9NmlP+tg+W5iHmef+GkB9hNWHZTS5stRu7ZV7/6cSpVb1sthhpxcgEgEEchGodLj05S9/ufz4tzvuuEM+9KEPybx582T+/PkbDJe+9rWvle/w+bd/+zeV+aa+c+moo46Siy++eJNrFx+rV7zj6NBDD5VLL71UevfuvUHMKaecIsWAqfiIvXXvOHr55ZdlwoQJ5QDt+OOPL99BVQzSNvWdS6tWrSo/ci/k0aSDJQeSkEyniyUv6WxTr+zb/+Q6TUZwxTWNAMMlS1f62FK7tRfmYea+Z4Cw3eyjqQ97c+2O5Eorp4/Lvf/1Mp0jqVW9LHbY6QWIRACBXAQqHS7tueeecsghh8g3v/nN0rPdcOnb3/52+f1Ev/71r1XmxUfqLVy4sPweo8GDB69fo/i+o/Hjx8u5555bDoA6Pe666y4544wzZNSoUeU9brbZZhtdfuSRR67/Zbi7tYqh0v777y8PPvignHjiiXLJJZds8P1Sq1evlr333ls+85nPyBVXXKF6vuuCmnSw5EASlOpkweQlGW3yhX37n1ynSQmuuKYR6Lyqb/9XcY9N2pM+ts8W5mHmub8GUB9h9WEZTa4stVt75d7/qUSpVb0sdtjpBYhEAIFcBCodLu2+++7yxS9+Uc4666zSs91w6Z/+6Z9kyZIl5fcvaR5PPfWUFO9OKt7BNHv27PVLzJgxo3x30LJly2TQoEFSfFbsypUrZdttt5Xttttu/XXFNdOnT5e99tqrHFJ1/Y6kdRc+8sgj8sYbb2xwi6+99lo5vCoGaGPGjJHi+RbvXiqGSMX3PQ0ZMkRuueWW9e+CuvHGG8uP2/vud78rBx98sObpro9p0sGSA0lQqpMFk5dktMkX9u1/cp0mJbjimkag86q+/V/FPTZpT/rYPluYh5nn/hpAfYTVh2U0ubLUbu2Ve/+nEqVW9bLYYacXIBIBBHIRqHS4VAx9incBLV68uPTsOlx69913y4+hK4Y9N910k9q8+Di74mPtxo0bV35/08MPP1wOrKZOnSrTpk0r13300Udl8uTJG/xZ8ZF9xbua+vbtK7NmzZL+/ftvcA9Dhw6V4cOHd3tf3X3nUhFQfBRgMezaZ599yuf4+9//Xq699tpyveuuu279R+tpn3STDpYcSLRZThtHXtL6plzdt//JdZps4IprGoHOq/r2fxX32KQ96WP7bGEeZp77awD1EVYfltHkylK7tVfu/Z9KlFrVy2KHnV6ASAQQyEWg0uFSMUQpvu/o1FNPLYc8xfctrfvOpeLdPXPmzCmHShdeeKF84QtfUJuvWbNGFixYUA6YXnnllfIdQ8XQaNKkSeuHOO2GS8X1695V1W7zYjhW3GN3j07DpSKm+Li94t1QxXc1bbPNNuU7nIp3VA0YMED9XNcFNulgyYEkON1JFiAvSVhNFvXtf3KdJi244ppGoPOqvv1fxT02aU/62D5bmIeZ5/4aQH2E1YdlNLmy1G7tlXv/pxKlVvWy2GGnFyASAQRyEah0uFS8M6l491DxfUjFR9NtscUW8oc//EE+9alPydNPPy2vvvqqHHTQQfK9730vF2+z59GkgyUHErOy8NqIvHhx1epi3/4n12nShyuuaQQYLlm60seW2q29MA8z9z0DhO1mH0192JtrdyRXWjl9XO79r5fpHEmt6mWxw04vQCQCCOQiUOlwqUBcu3atFN81VLxD6dlnny3/vXjstNNOctxxx5UfVderV69cvM2eR5MOlhxIzMrCayPy4sVVq4t9+59cp0kfrrimEei8qm//V3GPTdqTPrbPFuZh5rm/BlAfYfVhGU2uLLVbe+Xe/6lEqVW9LHbY6QWIRACBXAQqHy69H/Ltt9+WP/3pT7LVVltF+Wi4XJKkeR5NOlhyINFkOH0MeUlvnGoH3/4n12kygSuuaQQYLlm60seW2q29MA8z9z0DhO1mH0192JtrdyRXWjl9XO79r5fpHEmt6mWxw04vQCQCCOQiUKvhUi6odXgeTTpYciCpQ8VsfA/kpZ55cbkr3/4n1y6q/tfg6m/mEoErwyWXOol1DfUWS9J9Hczdrdpd6XsGCNvNPpr6sDfX7kiutHL6uNz7Xy/DcAm7VAL6dXmN1NsRiQAC9RIwHS6tXLlS/ewHDx6sju2JgU06WPJDtZ4VSl7qmReXu/Ltf3Ltoup/Da7+Zi4RuHZW8u1/F/OefA31Zp99zMPMc38NoD7C6sMymlxZarf2yr3/U4lSq3pZ7LDTCxCJAAK5CJgOl4YNG6b+/qQVK1bkYm7yPJp0sORAYlIS3puQF2+y2gT49j+5TpM6XHFNI8BwydKVPrbUbu2FeZi57xkgbDf7aOrD3ly7I7nSyunjcu9/vUznSGpVL4sddnoBIhFAIBcB0+HS9773vY2GS48//rg8/PDD8pGPfESGDx8uH/zgB+W1114r/183v/vd72Tfffct/3zq1Km5mJs8jyYdLDmQmJSE9ybkxZusNgG+/U+u06QOV1zTCDBcsnSljy21GS7F0PY9A8TY03INetJSO2wvchXmp4nOvf81Ji4x1KqLUvtrsMNOL0AkAgjkImA6XOqK9stf/lK+8pWvyDe/+U058sgjNzK944475Nxzz5WFCxfKyJEjczE3eR5NOlhyIDEpCe9NyIs3WW0CfPufXKdJHa64phFguGTpSh9bajNciqHtewaIsaflGvSkpXbYXuQqzE8TnXv/a0xcYqhVFyWGS3ol7GLbsR4CCNRLoNLh0rHHHivFdylddtll3arMnDlT/vd//1duvvnmesnV/G6adLDkMFfPYiIv9cyLy1359j+5dlH1vwZXfzOXCFwZLrnUSaxrqLdYku7rYO5u1e5K3zNA2G720dSHvbl2R3KlldPH5d7/epnOkdSqXhY77PQCRCKAQC4ClQ6Xdt99d5kyZYoUA6TuHsXgadGiRVJ8fB4Pd4EmHSw5kLjn1fJK8mKpHXcv3/4n13H9162GK65pBDqv6tv/Vdxjk/akj+2zhXmYee6vAdRHWH1YRpMrS+3WXrn3fypRalUvix12egEiEUAgF4FKh0uf/vSnZejQoXLjjTd263ncccfJH/7wh/J7mXi4CzTpYMmBxD2vlleSF0vtuHv59j+5juvPcCmNJ65urr7977Zqz72K10f73GMeZp77awD1EVYfltHkylKb4VKINrWq18MOO70AkQggkItApcOl4ruWbrjhBhk7dqxMmzat/Ii8dY+VK1fK9773PSm+d2nChAlyzjnn5GJu8jya9IslBxKTkvDehLx4k9UmwLf/yXWa1OGKaxqBzqv69n8V99ikPelj+2xhHmae+2sA9RFWH5bR5MpSm+FSiDa1qtfDDju9AJEIIJCLQKXDpTfffFO+/OUvl2/f7tOnj/zN3/yNbLfddvLHP/5RVq1aJe+++67sscce8v3vf1+23HLLXMxNnkeTfrHkQGJSEt6bkBdvstoE+PY/uU6TOlxxTSPAcMnSlT621G7thXmYue8ZIGw3+2jqw95cuyO50srp43Lvf71M50hqVS+LHXZ6ASIRQCAXgUqHSwXie++9J//6r/8qP/3pT+WZZ56RN954QwYMGCAf+9jH5IgjjpDPf/7z0rt371y8zZ5Hkw6WHEjMysJrI/LixVWri337n1ynSR+uuKYRYLhk6UofW2ozXIqh7XsGiLGn5Rr0pKV22F7kKsxPE517/2tMXGKoVRel9tdgh51egEgEEMhFoPLhUi6QdXseTTpYciCpW/XwH3fqmRH3u/Ltf3rQ3dbnSlx9tNyvxZXhknu1hF9JvYUb+q6Aua/Yhtf7ngHCdrOPpj7szbU7kiutnD4u9/7Xy3SOpFb1sthhpxcgEgEEchGo3XBp7dq10qtXr1x8K3seTTpYciCprEw6bkxe6pkXl7vy7X9y7aLqfw2u/mYuEbgyXHKpk1jXUG+xJN3Xwdzdqt2VvmeAsN3so6kPe3PtjuRKK6ePy73/9TIMl7BLJaBfl9dIvR2RCCBQL4HKh0vFMOmWW26RH//4x+VnrP/lL3+RLbbYQoYNGyZjx46VcePGMWxS1EyTDpb8UFUk2CCEvBggJ9rCt//JdZpE4IprGgGGS5au9LGldmsvzMPMfc8AYbvZR1Mf9ubaHcmVVk4fl3v/62UYLmGXSkC/Lq+RejsiEUCgXgKVDpdWr14tJ510kjzyyCPlAOlDH/qQDBw4UF599VV5+eWXy+9jGjlypCxYsED69etXL7ma302TDpb8UK1nMZGXeubF5a58+59cu6j6X4Orv5lLBK4Ml1zqJNY11FssSfd1MHe3anel7xkgbDf7aOrD3ly7I7nSyunjcu9/vQzDJexSCejX5TVSb0ckAgjUS6DS4dK8efOk+Oewww6TmTNnyg477LBeZ+XKlTJ37ly5++67ZerUqXLaaafVS67md9OkgyU/VOtZTOSlnnlxuSvf/ifXLqr+1+Dqb+YSgWtnJd/+dzHvyddQb/bZxzzMPPfXAOojrD4so8mVpXZrr9z7P5UotaqXxQ47vQCRCCCQi0Clw6WDDz5YttlmG/nRj37UrecXv/hFef311+Xee+/NxdzkeTTpYMmBxKQkvDchL95ktQnw7X9ynSZ1uOKaRoDhkqUrfWyp3doL8zBz3zNA2G720dSHvbl2R3KlldPH5d7/epnOkdSqXhY77PQCRCKAQC4ClQ6XdtllFzn++OPlq1/9areel112mVxzzTXy5JNP5mJu8jyadLDkQGJSEt6bkBdvstoE+PY/uU6TOlxxTSPAcMnSlT621Ga4FEPb9wwQY0/LNehJS+2wvchVmJ8mOvf+15i4xFCrLkrtr8EOO70AkQggkItApcOlT33qU7LvvvvKt7/97W49zzzzTHn44YflF7/4RS7mJs+jSQdLDiQmJeG9CXnxJqtNgG//k+s0qcMV1zQCDJcsXeljS22GSzG0fc8AMfa0XIOetNQO24tchflponPvf42JSwy16qLEcEmvhF1sO9ZDAIF6CVQ6XCq+Z+m+++6T+fPny6hRozaS+fnPf15+19JnP/tZKd7BxMNdoEkHSw5z7nm1vJK8WGrH3cu3/8l1XP91q+GKaxqBzqv69n8V99ikPelj+2xhHmae+2sA9RFWH5bR5MpSu7VX7v2fSpRa1ctih51egEgEEMhFoNLh0vPPPy9f+MIX5M0335Q999xTRowYIQMHDpRXX321PBj96le/kg984APldzLttNNOuZibPI8mHSw5kJiUhPcm5MWbrDYBvv1PrtOkDldc0wgwXLJ0pY8ttVt7YR5m7nsGCNvNPpr6sDfX7kiutHL6uNz7Xy/TOZJa1ctih51egEgEEMhFoNLhUoH49NNPyze+8Q15/PHHNzLdY4895LzzzpOPfexjuXibPY8mHSw5kJiVhddG5MWLq1YX+/Y/uU6TPlxxTSPQeVXf/q/iHpu0J31sny3Mw8xzfw2gPsLqwzKaXFlqt/bKvf9TiVKrelnssNMLEIkAArkIVD5cWgf50ksvlYOmN954QwYMGCB///d/L0OGDMnF2fx5NOlgyYHEvDycNiQvTky1vMi3/8l1mjTiimsaAYZLlq70saV2ay/Mw8x9zwBhu9lHUx/25todyZVWTh+Xe//rZTpHUqt6Weyw0wsQiQACuQjUZriUC2hdnkeTDpYcSOpSNRveB3mpZ15c7sq3/8m1i6r/Nbj6m7lE4MpwyaVOYl1DvcWSdF8Hc3erdlf6ngHCdrOPpj7szbU7kiutnD4u9/7XyzBcwi6VgH5dXiP1dkQigEC9BBgu1Ssf0e6mSQdLfqhGS3vUhchLVE7TxXz7n1ynSQ+uuKYRYLhk6UofW2q39sI8zNz3DBC2m3009WFvrt2RXGnl9HG5979ehuESdqkE9OvyGqm3IxIBBOolUPlw6Ze//KVce+218swzz8grr7wi77777kZCvXr1kqeeeqpecjW/myYdLPmhWs9iIi/1zIvLXfn2P7l2UfW/Bld/M5cIXDsr+fa/i3lPvoZ6s88+5mHmub8GUB9h9WEZTa4stVt75d7/qUSpVb0sdtjpBYhEAIFcBCodLt18881y/vnny9q1a+Vv//ZvZeDAgdK7d++2tosWLcrF3OR5NOlgyYHEpCS8NyEv3mS1CfDtf3KdJnW44ppGgOGSpSt9bKnd2gvzMHPfM0DYbvbR1Ie9uXZHcqWV08fl3v96mc6R1KpeFjvs9AJEIoBALgKVDpdGjx4tf/3rX+Vf/uVfZNiwYbmY1uJ5NOlgyYGkFiWz0U2Ql/+PvTOBtqOo9vdOgAQMgwyCEmQQkQACgggIKiCTAyAEAU1IoqDgg4QZQoAl8nwEkMGnMgjmqRgmZZZRwiSDPkCeMpj8AYMo8yAyhkFC/mv3fSfv3ptzz63a1VXdXefrtbIiN7Wrqr/f3n3qnp9VXU9dXGblW/9o7ULVvw1c/Zm5RMAVc8klT8pqQ76VRdK9H5i7s2rX0ncNEDZa+mjyIz1z64hoZSVnj8u9/u1kMJdgF4uAvV+ekXZ2REIAAvUiUKm5tP7668vuu+8uRx99dL2oZDCbJi0s+VCtZ8KhSz11cZmVb/2jtQtV/zZw9WfmEgFXzCWXPCmrDflWFkn3fmDuzgpzKYwV0XEJUMtx+XZj/cciSq7aycIOdnYCREIAArkQqNRcGj16tHzkIx+RE088MReetbkP3y+Xq5w4C5Iq6Q88NrrUUxeXWfnWP1q7UPVvA1d/Zi4RcMVccsmTstqQb2WRdO8H5u6suvHLZfIjLD9SRqNVSto9Y/n+DpB+hvUckVy16wI72NkJEAkBCORCoFJz6YYbbpAjjzxSfvnLX8oaa6yRC9Na3EeTFpYsSGqRMgtMAl3qqYvLrHzrH61dqPq3gas/M5cIuGIuueRJWW3It7JIuvcDc3dWmEthrIiOS4Bajsu3G+s/FlFy1U4WdrCzEyASAhDIhUCl5pJCvPbaa2Xq1Kmi719ac801ZcSIEW3Z7rzzzrkwT3Ifvl8uJ5nUAIOwIKmS/sBjo0s9dXGZlW/9o7ULVf82cPVn5hIB186UfOvfhXk3tyHf0qsP8zDmuT8DyI+w/EgZjVYpafeMlXv9xyJKrtrJwg52dgJEQgACuRCo1FyaM2dOsXNpxowZMm/evILpkCFD+rDVn+vPZs2alQvzJPfRpIUlC5IkKeE9CLp4I6tNgG/9o3Uc6eAK1zgEMJdScqWOU9LuGQvmYcx91wBho6WPJj/SM7eOiFZWcva43OvfTqZzJLlqJws72NkJEAkBCORCoFJzSY2lK664QtZZZx3ZdtttZbnllpOFFlqoLdtddtklF+ZJ7qNJC0sWJElSwnsQdPFGVpsA3/pH6zjSwRWucQhgLqXkSh2npI25VAZt3zVAGWOm7IOaTEk7bCy0CuNnic69/i1MXGLIVRdK7dvADnZ2AkRCAAK5EKjUXNpkk01ktdVWkwsuuECGDh2aC9Na3EeTFpYsSGqRMgtMAl3qqYvLrHzrH61dqPq3gas/M5cIuGIuueRJWW3It7JIuvcDc3dW7Vr6rgHCRksfTX6kZ24dEa2s5Oxxude/nUznSHLVThZ2sLMTIBICEMiFQKXm0sYbbyy77rqrTJ48OReetbmPJi0sWZDUJm36TARd6qmLy6x86x+tXaj6t4GrPzOXCLhiLrnkSVltyLeySLr3A3N3VphLYayIjkuAWo7LtxvrPxZRctVOFnawsxMgEgIQyIVApebSAQccIC+//LKce+65ufCszX34frlc5cRZkFRJf+Cx0aWeurjMyrf+0dqFqn8buPozc4mAK+aSS56U1YZ8K4ukez8wd2fVjV8ukx9h+ZEyGq1S0u4Zy/d3gPQzrOeI5KpdF9jBzk6ASAhAIBcClZpLTz/9tIwdO1Z22mkn2W+//WTYsGG5cK38Ppq0sGRBUnm6tJ0AutRTF5dZ+dY/WrtQ9W8DV39mLhFwxVxyyZOy2pBvZZF07wfm7qwwl8JYER2XALUcl2831n8souSqnSzsYGcnQCQEIJALgUrNpfHjx8srr7wiDz30kIwYMUJWWWWV4u/+15AhQ9jd5Jlxvl8ue3ZfanMWJKXiLK0zdCkNZfKOfOsfreNIBFe4xiGAuZSSK3WcknbPWDAPY+67BggbLX00+ZGeuXVEtLKSs8flXv92Mp0jyVU7WdjBzk6ASAhAIBcClZpLo0aNcuKo5tKsWbOc2tKoh0CTFpYsSOqZtehST11cZuVb/2jtQtW/DVz9mblEwLUzJd/6d2HezW3It/TqwzyMee7PAPIjLD9SRqNVStrN+w4gPZ2BRyRX7WrADnZ2AkRCAAK5EKjUXMoFYh3vo0m/WLIgqWMG8f8crqcqbrPyrX9q0I2rbyu4+hJzaw9XzCW3TCmnFflWDkefXmDuQ2vBtr5rgLDR0keTH+mZW0dEKys5e1zu9W8n0zmSXLWThR3s7ASIhAAEciHQSHNJP8D0z84775yLDqXfR5MWlixISpe/lA7RpRSMlXTiW/9oHUcmuMI1DgHMpZRcqeOUtHvGgnkYc981QNho6aPJj/TMrSOilZWcPS73+reTwVyCXSwC9n55RtrZEQkBCNSLQCPNpdNPP13OOOMMp6Py3nnnHTn77LPl0ksvleeff15Gjhwpe+65p4wdO1b0uL2Brpdeekkuv/xyueWWW2T27NkyZ84c+eAHPyg77LCDTJgwQYYPH94n9KSTTpJ77rlHHn/88aLtCiusIBtvvLHst99+stJKK/VpO27cOLn77rsXGHqhhRaSmTNnlpIhTVpY8qFaiuSld4IupSNN1qFv/aN1HGngCtc4BDr36lv/VcyxSWNSx+nVgnkY89yfAeRHWH6kjEarlLR7xsq9/mMRJVftZGEHOzsBIiEAgVwIZG8uHXPMMXLxxSFOhvIAACAASURBVBfL7rvvLuutt57ccccdcv3118ukSZNk4sSJA+qoptL+++8vm2++uWy66aay+OKLF+bR1VdfLRtuuKFMnz5d1AxqXWpWrbnmmrLyyivLiBEjCpPpkksuETW31NhSY6p1qbn00EMPydFHH91n/KFDh8qOO+5YSm41aWHJgqQUyUvvBF1KR5qsQ9/6R+s40sAVrnEIdO7Vt/6rmGOTxqSO06sF8zDmuT8DyI+w/EgZjVYpafeMlXv9xyJKrtrJwg52dgJEQgACuRDI2lyaNWtWcXTeXnvtJZMnT56v2UEHHSQ33XRT8Wf55Zdvq6WaQ3r1NoX0v3/wgx/ImWeeKbp7atttt+2YBw888IB8+ctfln322UcOPfTQ+W3VXPrb3/4mt912W7Q8atLCkgVJtDQI6hhdgvBVGuxb/2gdRy64wjUOgc69+tZ/FXNs0pjUcXq1YB7GPPdnAPkRlh8po9EqJe2esXKv/1hEyVU7WdjBzk6ASAhAIBcCWZtLp512WnEknu5CWnHFFedrdu+998qYMWPk2GOPLf72uXTH0U477SQHHnhgceRdp+vFF1+UT37yk/KVr3xFjjvuuPlNW+aSzuuNN94odjp1OqLPZ36ttk1aWLIgsSgcPwZd4jOONYJv/aN1HCXgCtc4BDr36lv/VcyxSWNSx+nVgnkY89yfAeRHWH6kjEarlLR7xsq9/mMRJVftZGEHOzsBIiEAgVwIZG0u6Y6lhx9+uDgKr/f19ttvy/rrry+jR4+W448/3kvL22+/Xb7xjW8UZpGaRr2vd999V/RdTXPnzpUnn3yyeC+U7k7Sv7fZZpv5TdVc0oXfwgsvLG+++WZx5N52220nhx12mCy77LJe8xmosfY/b968wriq+6UGm16LLbZY3afaVfNDl+rkHjVqVNDgvvWP1kG4BwyGK1wtBFLXv2WO3RRDHadXu5uZh9Z/68vlpvwOYMmubs4PC68qY9DKn37oM8D3dwD/GeYZQa7adYVdeexC698+EyIhAAEIhBHI2lzaYYcdZNiwYXLZZZctQEl3FK2zzjoybdo0Z4JqHk2YMEH0uLsbb7xRlltuuT6xTzzxhGy99dbzf7b00kvLt771Lfna177Wp92UKVNkhRVWKN7RpL/8/f73vy/ez6RH8OnfSy65pPOcBmrYpIUlC5JguaN0gC5RsDp1Grqw9K1/tHaSxbsRXL2ROQXkzjV1/TtB7+JGuedbHaXtZuah9a96+q4B6pgDnebUzfmBVk0j4D/f0GdA7vXvT9QtgueKG6d2rWBXHrvQ+rfPhEgIQAACYQSyNpd0t5AaQBdddNEClLbccsvCzJk+fbozwdYxe8ccc4zo7qP+11tvvSV65J7ujJo9e7ZcffXVhdmkx+cNHTq04zgXX3yxaL8TJ06USZMmOc9poIZN2hLPVupguaN0gC5RsCbp1Lf+0TqOLHCFaxwCnXv1rf8q5tikManj9GrBPIx57s8A8iMsP1JGo1VK2j1j5V7/sYiSq3aysIOdnQCREIBALgSyNpfK3Ll03nnnyXe/+90F3p/UKRGeeuop0TmoEXXwwQcPmjObbLKJrLbaam3NsEGD+zVo0sKSBYmvumnao0sazjFG8a1/tI6hgghc4RqHAOZSSq7UcUraPWPBPIy57xogbLT00eRHeubWEdHKSs4el3v928l0jiRX7WRhBzs7ASIhAIFcCGRtLg32zqVddtlFpk6dOqiWeqzeUUcdJV/4whfklFNOGXQXUu8O9f1MDz30kOi7mga7dD5z5syR3/zmN4M1HfTfm7SwZEEyqJyVNECXSrCXMqhv/aN1KdgX6ASucI1DoHOvvvVfxRybNCZ1nF4tmIcxz/0ZQH6E5UfKaLRKSbtnrNzrPxZRctVOFnawsxMgEgIQyIVAI80l/QCbNWuWqBnT6Tr11FPlnHPOkVtuuUVWXHHF+U316LoxY8bIt7/9bRk7dmzHPq655ho57LDDZIsttpDTTz9dFl54YS/tddfS/fffL/fdd1/HOH2fk+5c+vCHPywXXnih1xjtGjdpYcmCJFjuKB2gSxSsSTr1rX+0jiMLXOEah0DnXn3rv4o5NmlM6ji9WjAPY577M4D8CMuPlNFolZJ2z1i5138souSqnSzsYGcnQCQEIJALgUaaS67wZ86cWRhQuoNp8uTJ88MOOuggufHGG+Wmm26SFVZYQfQlhHqE3dJLLy3LLLPM/Hba5sADD5RPfOIThUk1bNiwtkO/8sorsthii8kiiyzS59/1g3a33XaTj33sY/Pf7fTaa68V7YYPH96n7bRp0+Tkk0+WQw45RPbdd1/XWxywXZMWlixIguWO0gG6RMGapFPf+kfrOLLAFa5xCHTu1bf+q5hjk8akjtOrBfMw5rk/A8iPsPxIGY1WKWn3jJV7/cciSq7aycIOdnYCREIAArkQSGoujRo1SoYMGeLNTmPUKLJcepydHmu3++67y7rrrit33nmnXHfddTJx4kSZNGlS0eVdd90l48eP7/Mz3W2ku5rUCDriiCMK86j3tfLKK8sGG2xQ/EhNqO985zvyuc99TlZZZRVZaKGF5JFHHpErrrii+Pdzzz1X1ltvvflj6fuX9Ig97UPvTcefMWOGKB/dtfSe97zHcqt9Ypq0sGRBEix3lA7QJQrWJJ361j9ax5EFrnCNQ6Bzr771X8UcmzQmdZxeLZiHMc/9GUB+hOVHymi0Skm7Z6zc6z8WUXLVThZ2sLMTIBICEMiFQFJz6cgjjzSZSwr7hBNOMDH/17/+JWeffXZhMD333HMycuTIwjTS4+paRlc7c0nbT5kyZcAxdUfUiSeeWPz73//+dznrrLOKxZyOoWMuv/zyxTF3++yzj6y22mrz+3niiSeK9zY9+OCD8sILL8jcuXNlpZVWku22265oO2LECNN99g9q0sKSBUkpkpfeCbqUjjRZh771j9ZxpIErXOMQ6Nyrb/1XMccmjUkdp1cL5mHMc38GkB9h+ZEyGq1S0u4ZK/f6j0WUXLWThR3s7ASIhAAEciGQ1FzKBVoT7qNJC0sWJH4Z9cbrb8rbb7wtwxYbJouNWNQv2KM1unjAqllT3/pH6zgCxuKa6hkQh0p4r7G4hs+sHj341n89Zl3fWZBv6bWBeRjz3J8B3ZwfTfv872atwqrYHp17/dvJdI5sSq7W8RnQFHaxciekX9iF0CMWAhCoEwHMpTqpUeJcmrSw5EPVTfhXX3xNnnnsOfnl966Up2c/Ix9Y/f2yxxFfkvevurwssczibp14tEIXD1g1a+pb/2gdR8CyuaZ+BsShEt5r2VzDZ1SvHnzrv16zr99syLf0msA8jHnuz4BuzI+mfv53o1Zh1RsenXv9hxNq30Pdc7XOz4C6s4uVM2X0C7syKNIHBCBQBwKYS3VQIcIcmrSw5EN18ATQBeXFp10lF069bIHGXz1qtOx2yI6lG0zoMrgudW3hW/9oHUfJMrlW8QyIQyW81zK5hs+mfj341n/97qBeMyLf0usB8zDmuT8Dui0/mvz5321ahVVuOdG51385lBbspc65WvdnQJ3ZxcqXsvqFXVkk6QcCEKiaQFJzaeuttzbdr74b6cYbbzTFdmtQkxaWfKgOnqWP/M+jst9GkwdseOYfTpI1NvzQ4B15tEAXD1g1a+pb/2gdR8AyuVbxDIhDJbzXMrmGz6Z+PfjWf/3uoF4zIt/S6wHzMOa5PwO6LT+a/PnfbVqFVW450bnXfzmUmmUu1f0ZQJ3bsxJ2dnZEQgAC9SKQ1FwaN26c+e6nT59uju3GwCYtLPlQ7ZyherbyqXufJb/91e8GbLjlHpvJodP+TRYt8R1M6NLcJ4dv/aN1HK3L4lrVMyAOlfBey+IaPpN69uBb//W8i/rMinxLrwXMw5jn/gzopvxo+ud/N2kVVrXlRede/+WR6ttTXXO1Cc+AurKLlStl9gu7MmnSFwQgUCWBpOZSlTfabWM3aWHJh2rn7Hz5hVfkqM8fLw/f++iADT+y0eoy9dqjZKnlliwt1dGlNJTJO/Ktf7SOI1FZXKt6BsShEt5rWVzDZ1LPHnzrv553UZ9ZkW/ptYB5GPPcnwHdlB9N//zvJq3Cqra86NzrvzxSzTCXmvAMoM7tWQk7OzsiIQCBehHAXKqXHqXNpkkLSz5UO8te1f9jCV1KK8fkHfnWP1rHkagsrlU9A+JQCe+1LK7hM6lnD771X8+7qM+syLf0WsA8jHnuz4Buyo+mf/53k1ZhVVtedO71Xx6pZphLTXgGUOf2rISdnR2REIBAvQhgLtVLj9Jm06SFJR+qg8texVnL6DK4LnVt4Vv/aB1HyTK5VvEMiEMlvNcyuYbPpn49+NZ//e6gXjMi39LrAfMw5rk/A7otP5r8+d9tWoVVbjnRudd/OZQW7KXOuVr3Z0Cd2cXKl7L6hV1ZJOkHAhComkCl5tL48eOd7n/IkCFy7rnnOrWlUQ+BJi0s+VAdPGtfffE1ufi0q+TCqZct0HjM0aPlywfvKEsss/jgHXm0QBcPWDVr6lv/aB1HwDK5VvEMiEMlvNcyuYbPpn49+NZ//e6gXjMi39LrAfMw5rk/A7otP5r8+d9tWoVVbjnRudd/OZSaZS7V/RlAnduzEnZ2dkRCAAL1IlCpufTZz362LY3XX39dXn75ZVFTabnllpNFFllEbr755nqRq/lsmrSw5EPVLZl0YfnMY8/Jr06+Up6a/aysuPoKsvvhX5L3r7p86caSzghd3HSpYyvf+kfrOCqWzTX1MyAOlfBey+YaPqN69eBb//Waff1mQ76l1wTmYcxzfwZ0Y3409fO/G7UKq97w6NzrP5xQ+x7qnqt1fgbUnV2snCmjX9iVQZE+IACBOhCo1FzqBODJJ5+Uk046SZ599ln56U9/KiNGjKgDr8bMoUkLSz5U/dLqzdfflLfeeFuGLzZMFh2xqF+wR2t08YBVs6a+9Y/WcQSMxTXVMyAOlfBeY3ENn1k9evCt/3rMur6zIN/SawPzMOa5PwO6OT+a9vnfzVqFVbE9Ovf6t5PpHNmUXK3jM6Ap7GLlTki/sAuhRywEIFAnArU1lxTSO++8I7vssotstNFGcuyxx9aJW+3n0qSFJR+q9UwndKmnLi6z8q1/tHah6t8Grv7MXCLg2pmSb/27MO/mNuRbevVhHsY892cA+RGWHymj0Sol7Z6xcq//WETJVTtZ2MHOToBICEAgFwK1NpcU8vHHHy/XXHON/O53v8uFeZL7aNLCkgVJkpTwHgRdvJHVJsC3/tE6jnRwhWscAphLKblSxylp94wF8zDmvmuAsNHSR5Mf6ZlbR0QrKzl7XO71byfTOZJctZOFHezsBIiEAARyIVB7c+nII4+U6667Tu67775cmCe5jyYtLFmQJEkJ70HQxRtZbQJ86x+t40gHV7jGIYC5lJIrdZySNuZSGbR91wBljJmyD2oyJe2wsdAqjJ8lOvf6tzBxiSFXXSi1bwM72NkJEAkBCORCoLbm0rx58+Tqq6+WKVOmyPrrry/nn39+LsyT3EeTFpYsSJKkhPcg6OKNrDYBvvWP1nGkgytc4xDAXErJlTpOSRtzqQzavmuAMsZM2Qc1mZJ22FhoFcbPEp17/VuYuMSQqy6UMJfslGBXNjv6gwAE6kWgUnNp6623bktj7ty58o9//KN459KIESPkZz/7may77rr1Ilfz2TRpYclirp7JhC711MVlVr71j9YuVP3bwNWfmUsEXDGXXPKkrDbkW1kk3fuBuTurdi191wBho6WPJj/SM7eOiFZWcva43OvfTqZzJLlqJws72NkJEAkBCORCoFJzady4cW05Dh06VJZccklZZ511ZPTo0bL88svnwjvZfTRpYcmCJFlaeA2ELl64atXYt/7ROo58cIVrHAKYSym5UscpafeMBfMw5r5rgLDR0keTH+mZW0dEKys5e1zu9W8ng7kEu1gE7P3yjLSzIxICEKgXgUrNpXqhyGs2TVpY8qFaz9xDl3rq4jIr3/pHaxeq/m3g6s/MJQKumEsueVJWG/KtLJLu/cDcnVW7lr5rgLDR0keTH+mZW0dEKys5e1zu9W8ng7kEu1gE7P3yjLSz6/bIH/3oR3L66afLL37xC9lkk026HUep9//Zz35WRo4cKdOnTzf1e9ddd8n48ePlhBNOKDbLdMuFuZSp0k1aWPKhWs8kRJd66uIyK9/6R2sXqv5t4OrPzCUCrp0p+da/C/NubkO+pVcf5mHMc38GkB9h+ZEyGq1S0u4ZK/f6j0WUXLWThR3s7ASqi1xzzTWdB99ll13kxBNPdG4/WMPLLrtMpkyZEmQ+1MFcapkoEydOlEmTJg122435d8wlm1SYSzZutY9q0sKSBUk90wld6qmLy6x86x+tXaj6t4GrPzOXCLhiLrnkSVltyLeySLr3A3N3Vu1a+q4BwkZLH01+pGduHRGtrOTscbnXv51M50hy1U4WdrCzE6guUs2Z3tcrr7xS7ALSHStqJvW+1lprLdlmm21KmyzmUmkoo3SEuWTDWrm5dMsttxTbzWbOnCmvvvqqvPvuuwvcyZAhQ4p/53In0KSFJQsSd11TtkSXlLTLHcu3/tG6XP6t3uAK1zgEOvfqW/9VzLFJY1LH6dWCeRjz3J8B5EdYfqSMRquUtHvGyr3+YxElV+1kYQc7O4H6RD7xxBOy9dZby8Ybb2w+Ds31bjCXXElV0w5zyca9UnPp0ksvlWOOOUYWWmgh+fjHPy7LL7+8LLzwwm3vRM8r5HIn0KSFJQsSd11TtkSXlLTLHcu3/tG6XP6YS3F4wtWNq2/9u/Xava14PqbXHuZhzHN/BpAfYfmRMhqtUtLGXAqhTa7a6cEOdnYC9YnsZC49/fTTcsYZZ8jtt98u//jHP2S55ZaT7bbbrjgKboklluhzE7/+9a/lvPPOk8cee0zeeustWWaZZWT99dcv2q6++upy5JFHyuWXX77AjfuaWr2PxXvooYeKMXWeK620kuy1116y2267LTDGSy+9JD/+8Y/lxhtvlGeeeUaWWmop+cxnPiMHHXSQrLDCCn3a//a3v5Vp06bJX/7yF3nttddk6aWXlrXXXlv22Wcf2XDDDaU1fv9BdOfXzTff7Cxsb6NtxIgRctZZZ8mjjz5azGffffeVL3/5ywXH//zP/5RrrrlG/vnPf8p6660n3/nOd2SNNdZYYJyrrrqqMAcffvhhGTp0aDHnb37zm7LFFlss0PbBBx+U733ve3L//ffL8OHDZcstt5TJkycXY7Z755JrHvDOJWf5y2u4/fbby5w5c+TCCy8sioCrPAJN+sWSBUl5upfZE7qUSTNtX771j9Zx9IErXOMQ6Nyrb/1XMccmjUkdp1cL5mHMc38GkB9h+ZEyGq1S0u4ZK/f6j0WUXLWThR3s7ATqEzmQuTR79mzZc889i1O2WjtaHnnkkcJoUuPioosuKowJvfRYveOPP15WXXVV+fSnP138XE2c3//+93L00UfLF7/4xcLYUUPlpptuKnZK6ZF7eqmZMXr0aGcgLXNHzaH77ruv6HuRRRaR6667Tp577rnCxPr6178+v7/nn39exo4dK3//+98LQ0mNrieffFJmzJgh73//+0U3fqgRppfOcf/99y8Mnq222kqWXHJJefbZZ+UPf/iDfOUrXykMJjVQ1CTTP2qM6R+91Gz72te+5nwfLXNJx9E+1bRbfPHF5xtJZ555ZsFY9dl8882Le/vNb35T8Lrhhhv6bE5RA/CHP/xhcT/qNbzzzjty7bXXippq//Ef/1GYRq1LT0ZTHv/6178Kdu973/vktttuEz01TQ3E1VZbrc8ONp88wFxylr+8huuuu26RnFpoXOUSaNLCkgVJudqX1Ru6lEUyfT++9Y/WcTSCK1zjEOjcq2/9VzHHJo1JHadXC+ZhzHN/BpAfYfmRMhqtUtLuGSv3+o9FlFy1k4Ud7OwE6hM5kLm06667ihoLuiGiZQTprFtG0qGHHlqYLXrpu5rUmFDTY9FFF51/c2pyvPHGG/N3OZV5LN5iiy0mV155payyyirFeC+++KJ86Utfkpdffln0FTTLLrts8fMDDjigMJLOPvvswlxqXWpy7bfffrLHHnvIv//7vxc/njhxotx6662F2dIynPTn8+bNK/p973vfW7RrmSjaXndmWa4Wi2HDhsnFF18so0aNKrqZNWuW7LzzzgUzNfF+8pOfzDfx1MBT/rqb6fOf/3zRXnc7qUmkm1YuueSSYleWXrrbSPvR3U/KQ3dg6aU+xB//+Ef5+c9/Lp/85CeLn82dO7fY5XTnnXcucDyiTx5gLlkyITDmC1/4gnzsYx+TqVOnBvZEeH8CTVpYsiCpZ/6iSz11cZmVb/2jtQtV/zZw9WfmEgHXzpR869+FeTe3Id/Sqw/zMOa5PwPIj7D8SBmNVilp94yVe/3HIkqu2snCDnZ2AvWJbGcuPfDAA8VuFzUcDjvssD6Tfffdd+VTn/pUsUtGDRK91Fx65ZVXit1DapYMdJVpLn31q18tjojrfakRc8opp8hxxx1XmChqeOlcP/e5z8n3v//9Baalxonev5oieqlZpAaLHo2nu5YGuso0l/QYP91d1PvSXUx/+9vf5IILLiheo9O67r33XhkzZkyxu0pNM71aO7m0j/5HAp5++unFv6t5piZaS+tNNtmkMKl6X3/605+KNr2PKfTNA8ylCur6V7/6lZx66qlyxRVXyAc+8IEKZpDvkE1aWLIgqWceoks9dXGZlW/9o7ULVf82cPVn5hIB186UfOvfhXk3tyHf0qsP8zDmuT8DyI+w/EgZjVYpafeMlXv9xyJKrtrJwg52dgL1iWxnLp1//vmFIaFHrLV7v49+n63vI9IdMHrprqDTTjutOFJthx12kE033bR4P1B/o6lMc0nfGaQ7lXpfenydHvmmBsyxxx5bmES6u2qzzTYr3pfU/7r++uuLdyvp8X26U+nqq68W3ZGlxtmOO+5Y7OzRON0l1fsq01w65phjZNy4cX361/mrkaR/9Ki81qVH+2277baFidQypNQQ051Zei/Kv/d19913F323eLR2a33rW9+Sgw8+uE9bNQ31hDW9X313k16+eYC5VEFd33PPPfKzn/1M1B0cP368rLnmmn2SpveUPvGJT1Qww+YO2aSFJQuSeuYZutRTF5dZ+dY/WrtQ9W8DV39mLhFw7UzJt/5dmHdzG/ItvfowD2Oe+zOA/AjLj5TRaJWSds9Yudd/LKLkqp0s7GBnJ1CfyHbm0llnnVUcvTbY9dBDDxVN9Ng43WWj7wh6+OGHi5+pKaImyCGHHDLfZCrTXJo2bVrxfqfelx4Rp8fF7bTTTnLyySfLr3/9azn88MMHuw25+eabi3cZ6aXvKtIj4+6///7ivvT9UXrs3JQpU+bvZirTXDrhhBMWeOeUGkJqDLX4tm6gpZXuFDvxxBOLH+t7ntQcU2NNj9LrfemxhnpqWouHHiN4xBFHyFFHHSUTJkxYgIvu8ur9ziXfPMBcGjTVym+g5ynqC7M0WfXvTpeeucjlTqBJC0sWJO66pmyJLilplzuWb/2jdbn8W73BFa5xCHTu1bf+q5hjk8akjtOrBfMw5rk/A8iPsPxIGY1WKWn3jJV7/cciSq7aycIOdnYC9YlsZy613qukBo0aEz7Xs88+W5gd+q4m3Uyx9957F4aGXmWaSy47l1o7dQ466CD5t3/7N5/bKN7hpAaPvsfo9ttvL47W+8EPflD0USdzqcydS7rbbIMNNpi/c8k3DzCXvFKsnMZ67uFgplJrJE0WLncCTVpYsiBx1zVlS3RJSbvcsXzrH63L5d/qDa5wjUOgc6++9V/FHJs0JnWcXi2YhzHP/RlAfoTlR8potEpJu2es3Os/FlFy1U4WdrCzE6hPZDtzSY+703cW6Q4aPbbNcr355pvFsXJ6xJy+i0kvfS3M5MmT5fjjjy/e6WS5Wu8Y6vTOJX0Xk/67Gl1bbLGFbLnllvLjH//YMpzocXG6G0r7Ui76Pb6eRLbnnnvKfvvtJwceeKCp305Gm8/OpdZ7ldoxPfPMMwtDzPrOJd88wFwypQJBdSXQpIUlC5J6ZhG61FMXl1n51j9au1D1bwNXf2YuEXDtTMm3/l2Yd3Mb8i29+jAPY577M4D8CMuPlNFolZJ2z1i5138souSqnSzsYGcnUJ/IduaSnrClR6/p+4j+67/+SzbZZJM+E9b3LT3++OOy1lprFT9Xs6X/61yef/552WqrrWT11VcXPY5Nr1tuuUX0fT8HHHCA7L///iYILXNJ34Ok/a6yyipFP7rTaOedd5aXXnqpGGfZZZctfq4GkO5garcL66233iqOntMdO3rp0XL6zqGhQ4fOn9ucOXNku+22E22r96mXHv2n72TaddddZerUqab7KMtc+utf/1ocfbfyyisXu6xaR+OpGabvpFKTT3ksvfTSxTzVNFTTSI/+U/NPr7lz58o3v/lNufPOO2XjjTeev3PJNw8wl0ypUE3QueeeK7o1TYuDqz2BJi0sWZDUM4vRpZ66uMzKt/7R2oWqfxu4+jNziYBrZ0q+9e/CvJvbkG/p1Yd5GPPcnwHkR1h+pIxGq5S0e8bKvf5jESVX7WRhBzs7gfpEtjOXdHb6vp7x48fLP/7xD9lss81kjTXWkHfeeacwlfS4OD0uT3fE6LXRRhvJe9/73sKkWXHFFeXVV18tvrNWg0nfDaRGlV5qAOlOIjWG9GdLLbVU0V5NIderZS595jOfkfvuu694H9IiiyxS7I567rnn5Mgjj5Svf/3r87t74YUXZOzYsfLYY48VR7599KMfLcyjJ598srgPnbMaaHqpGaP3qwbTSiutJG+//bbceuutxT1PmjRJWqeKKQe9DzXZRo8eLe973/uK9zHpbibXqyxzScdr7V76wAc+INtvv32h0zXXXFMYbd/97neLd1+1rpkzZ8qYMWOKphbLeQAAIABJREFUNspO567H/uml9977nUu+eYC55Kp+Ddpp0pxxxhnCe5gGFqNJC0sWJDUoqjZTQJd66uIyK9/6R2sXqv5t4OrPzCUCrp0p+da/C/NubkO+pVcf5mHMc38GkB9h+ZEyGq1S0u4ZK/f6j0WUXLWThR3s7ATqEzmQuaQzVHPonHPOKQyWp59+WkaMGFEcc6dmkx5rp7uS9LrggguKNroLSA0KNZrWXHNN2WuvvWTzzTfvc7MzZswozJBHH320MG9675RxodIyl3TThY533nnnyVNPPVWYQfp+p95GSqs/Nbt++tOfyg033FAYRWpG6X3o2GpytXYuqSFz/fXXy5///GdRU0rv90Mf+lBhTunuoN6X7nI65ZRTRJ8Db7zxhowcOVJuvvlml1so2pRpLml/v/71r4sdR4888khxdN/aa68t++yzT2GC9b8eeOCBYieXmnOLLrpo0UZNOdVU70P76X255gHmkrP81TfEXBpcgyYtLFmQDK5nFS3QpQrq5YzpW/9oXQ73/r3AFa5xCHTu1bf+q5hjk8akjtOrBfMw5rk/A8iPsPxIGY1WKWn3jJV7/cciSq7aycIOdnYCREIAArkQGDJPDxBs2IW5NLhgTVpYsiAZXM8qWqBLFdTLGdO3/tG6HO6YS3E4wtWPq2/9+/Xefa15PqbXHOZhzHN/BpAfYfmRMhqtUtLGXAqhTa7a6cEOdnYCREIAArkQwFzKRcl+99GkXyxZkNQzCdGlnrq4zMq3/tHahap/G7j6M3OJgGtnSr7178K8m9uQb+nVh3kY89yfAeRHWH6kjEarlLQxl0Jok6t2erCDnZ0AkRCAQC4EMJdyURJzKVMlq7stForVsQ8d2feLJbQOJd4+Hq5wjUMAcyklV+o4Je2esWAextx3DRA2Wvpo8iM9c+uIaGUlZ4/Lvf7tZDpHkqt2srCDnZ0Akb0J/PznPxd9L1Kna4kllpCvfe1rtQan71F68sknB53jpEmTBm1Dg+YQwFxqjlZeM23SwpIFiZe0yRqjSzLUpQ/kW/9oXboEfEEaBylcHbj61r9Dl13dhOdjevlhHsY892cA+RGWHymj0Sol7Z6xcq//WETJVTtZ2MHOToDI3gQ++9nPDmrKjBw5Um6++eZagxs3bpzcfffdg87xoYceGrQNDZpDAHOpOVp5zbRJC0sWJF7SJmuMLslQlz6Qb/2jdekSYILEQQpXB66+9e/QZVc34fmYXn6YhzHP/RlAfoTlR8potEpJG3MphDa5aqcHO9jZCRAJAQjkQgBzKRcl+91Hk36xZEFSzyREl3rq4jIr3/pHaxeq/m3g6s/MJQKunSn51r8L825uQ76lVx/mYcxzfwaQH2H5kTIarVLSxlwKoU2u2unBDnZ2AkRCAAK5EGikuaRnUf7iF7+o/XbAKpOkSb9YsiCpMlMGHhtd6qmLy6x86x+tXaj6t4GrPzOXCLhiLrnkSVltyLeySLr3A3N3Vu1a+q4BwkZLH01+pGduHRGtrOTscbnXv51M50hy1U4WdrCzEyASAhDIhUAjzaVc4Me8jyYtLFmQxMwEe9/oYmdXdaRv/aN1HMXgCtc4BDCXUnKljlPS7hkL5mHMfdcAYaOljyY/0jO3johWVnL2uNzr304Gcwl2sQjY++UZaWdHJAQgUC8ClZtLb775psyYMUNmzpwpr776qsydO3cBQkOGDJGpU6fWi1zNZ9OkhSUfqvVMJnSppy4us/Ktf7R2oerfBq7+zFwi4NqZkm/9uzDv5jbkW3r1YR7GPPdnAPkRlh8po9EqJe2esXKv/1hEyVU7WdjBzk6ASAhAIBcClZpLjz76qOy9997yzDPPyLx58wZkqubSrFmzcmGe5D6atLBkQZIkJbwHQRdvZLUJ8K1/tI4jHVzhGocA5lJKrtRxSto9Y8E8jLnvGiBstPTR5Ed65tYR0cpKzh6Xe/3byXSOJFftZGEHOzsBIiEAgVwIVGou7bXXXvL73/9eDjjgANl5551l+eWXl4UWWigXtpXeR5MWlixIKk2VAQdHl3rq4jIr3/pHaxeq/m3g6s/MJQKumEsueVJWG/KtLJLu/cDcnVW7lr5rgLDR0keTH+mZW0dEKys5e1zu9W8ng7kEu1gE7P3yjLSzIxICEKgXgUrNpY997GOy1VZbyfe///2oVN555x05++yz5dJLL5Xnn39eRo4cKXvuuaeMHTtWdFfUQNdLL70kl19+udxyyy0ye/ZsmTNnjnzwgx+UHXbYQSZMmCDDhw/vE3rSSSfJPffcI48//njRdoUVVpCNN95Y9ttvP1lppZUWGEaNtR/84AfFrqzFFlusYHH44YfLMsssE8yjSQtLPlSD5Y7SAbpEwZqkU9/6R+s4ssAVrnEIdO7Vt/6rmGOTxqSO06sF8zDmuT8DyI+w/EgZjVYpafeMlXv9xyJKrtrJwg52dgJEQgACuRCo1FzadNNNix1LRx55ZFSexxxzjFx88cWy++67y3rrrSd33HGHXH/99TJp0iSZOHHigGOrqbT//vvL5ptvLjrXxRdfvDCPrr76atlwww1l+vTpfXZaqVm15pprysorrywjRowoTKZLLrlE1NxSY0uNqdZ19913y9e//vWi/Ze//GV58cUX5ac//amsuOKKRcyiiy4axKRJC0sWJEFSRwtGl2hoo3fsW/9oHUcSuMI1DgHMpZRcqeOUtHvGgnkYc981QNho6aPJj/TMrSOilZWcPS73+reT6RxJrtrJwg52dgJEQgACuRCo1Fw6+uijZebMmXLZZZd13EEUAlt3BamBpUfwTZ48eX5XBx10kNx0003FHz2Or92l5pBevU0h/W/dbXTmmWfK6aefLttuu23H6T3wwAOFebTPPvvIoYceOr+tzunll1+Wa665Rt7znvcUP//tb39btJsyZYp87WtfC7ntRv2/lliQBEkdLRhdoqGN3rHvL5ZoHUcSuMI1DoHOvfrWfxVzbNKY1HF6tWAexjz3ZwD5EZYfKaPRKiXtnrFyr/9YRMlVO1nYwc5OgEgIQCAXApWaS6+++mqxe2fVVVctjoPTY+TKvk477bTiSDzdhaS7glrXvffeK2PGjJFjjz22+Nvneuihh2SnnXaSAw88sDjyrtOlO5I++clPyle+8hU57rjjiqZ//etf5XOf+1zbnVPbbbedLLXUUsVOq5CrSQtLFiQhSseLRZd4bGP37Fv/aB1HEbjCNQ6Bzr361n8Vc2zSmNRxerVgHsY892cA+RGWHymj0Sol7Z6xcq//WETJVTtZ2MHOTqC+kdZXq+gdXXvttXLrrbfK/fffL4899lixmeG2225re7Pjxo0TPdWq/7XQQgsVGzH6X3PnzpULL7ywOO1Kv1ceNmyYrL766sWJW5/+9KfnN9d2//Vf/1WcoPXkk0/Ke9/73mJjxMEHHyxLLrlkn25d2/q+NkYHcX0VjN7TXXfdVTDT+epJYfqz/pduTNHNGANdm222mfzsZz+b/8+f/exni/76X3ra2IwZM/r82HUOraC33npLpk2bVpxs9sQTTxSbRvRksiOOOEI++tGPmubQe0LqJXz+858X5f4f//Efsttuu/WZ7z//+c/5XsfTTz9dvF7n4x//eOFRaE70v5599tlig4xuatG+tf36668vJ5xwQnFKW+hVqbm09dZby7/+9a/iPUh6aZK3uyl9L9KNN95oulfdsfTwww8XR+H1vt5+++0C5OjRo+X444/36vv222+Xb3zjG4VZpKZR7+vdd98txNcC1SQ+44wzigeJ/r3NNtsUTa+66io57LDDikTs/QDQf9Of/+Y3v5E//elPfY7c85rg/y4s582bVxzPV/frjTfeKKao753iqg8BdKlOi1GjRgUNrr9Y+tQ/WgfhHjAYrnC1EEhd/5Y5dlMMdZxe7W5mHlr/qpbvGiC9wmEjdnN+hJFLH41W/sxDnwG5178/UbcIctWNU7tWsCuPXWj922dCZH8C1leraD9qGD344IOyzjrrFObS0KFDO5pLunlBT/XqfWnMjjvu2Odn+l3zAQccUJgDu+yyi6y77rqi9feXv/yl+N+9zQfdvPHrX/+6MCc22WQT+fvf/y7nn3++fOQjHylMm0UWWWR+365tfV8b4/MqGDWBdPOJ3scf//hH0VpoZy7p6WKt/xNFbzh33nmnXHnllcUrd3QDS+vSfvVe+28K0e/JW9/P927rMgdtP2fOnOJ0tEceeaQ4qWyNNdaQ1157TfTkNGW+5ZZbmubQ+56OOuooue6664qx+ptL6meol6GmlnoSaibp/1Zmaoyq9iuttNL87mbPnl3kpd63nqKmG3vUYFKW3/3ud+V973tf8EOgUnNJhXa9br75ZtemfdrtsMMOhZurDmf/S3cUacGryeN6aUFPmDBB9Lg7NbyWW265PqEqqJpmrWvppZeWb33rW32OuVMH+Xvf+14huDqbvS/9uf67Fkf/vl3nqO2atLBkQeKjbLq26JKOdf+RQheWvvWP1nG0hitcLQRS179ljt0UQx2nV7ubmYfWf9N+B7BkVzfnh4VXlTFo5U8/9Bng+zuA/wzzjCBX7brCrjx2ofVvn0l9It94/U15+423Zdhiw2SxEWHvobfeVcirVXRM3UWiu5V095F+of+3v/2to7nU6d9738MvfvELOemkk+Tcc8+VjTbaaMDbU2Nr11137XN6lja+4YYbitOzvvOd78hXv/rVIt6nre9rY3xeBaMbM/SkMd1Yoj6Bmh/tzKWBblpfK3PPPfcUnJdddtn5zXz68pmD7vbR3WN64tiHPvShjqnmM4dWR/pZPnbsWNHX+ehpbP3NpdYrddQE1RxrXfraHzXSdPfU3nvvXfxY/0/nmg96TZ8+PdoGlErNJWux+8SpG6kmzUUXXbRAmLqJ+j4lBex6tY7Z6y9iK163xumRe+okqjuoW+TUbFKB1X3WS3cx/fCHP5Trr79eVltttT5Dt97npEnR22l0nV/vZNT/rdsJ636xlbqeCqFLPXVxmZXvkRho7ULVvw1c/Zm5RMC1MyXf+ndh3s1tyLf06sM8jHnuzwDyIyw/UkajVUraPWPlXv+xiJKrdrKwg52dwP9Fvvria/LMY8/JL793pTw9+xn5wOrvlz2O+JK8f9XlZYllwo/s8pljma9WcTWXdFeQGrW6s0QNlv6XbnLQ75V1Z49+l6z/3Wrfv+3Pf/7z4qizCy64oDgmrfe1wQYbyFprrVX8m14+bQdi2O61MSGvgvE1Y9TM05gttthCfvzjH/eZZqsvNeb0u3rX4986zUF3KH3qU58qDLrJkycXO4X0RLaBTuLynYOegqa7ktRs1r/Hjx+/gLl0zTXXyCGHHCLqH+grd1qXHiuoO9i+/e1vF+aUXno0oZpvymarrbaSN998szA+e+9e86mPgdpmby6VuXPpvPPOK7aM9X5/0mAiPPXUU6Jz0IeKnm+pV6qdSzoW5tJgCvHvAxFgodjc3PD9xRKt42gNV7jGIdC5V9/6r2KOTRqTOk6vFszDmOf+DCA/wvIjZTRapaTdM1bu9R+LKLlqJws72NkJ9ESqsXTxaVfJhVMXPG3qq0eNlt0O2TGpwVTmq1VczCV9bi+88MLFl/5qfmy33XbF61J678DRo++++MUvFt8pP/fcc8XJXGoufeADHyhOyur9upZzzjlHTj31VLn88stl7bXX7iOPnt6lcXr0nJpYPm0H0rnda2NCXgXjay6pafL973+/MN223377PtPUvvQ1PLp7Rw0gPVlMjxtUjvqOpIGuTnNo7RrS3UT//d//XbzWRvvWjSNq+Kh+vS/fOejOtNZmlEcffbStuaSGmo4zcuTIwkhqHYt34oknFq/p0Xdttd6t1TodTTfVqHGq2uvGF939pptm+p+mZq3n7M2lwR4Melbl1KlTB+WnxatnHn7hC1+QU045Zf4upEEDRYr3M6mbq0Wn12CFpjua7rvvvuB3LulYmEsuCtGmHQEWis3NC99fLNE6jtZwhWscAp179a3/KubYpDGp4/RqwTyMee7PAPIjLD9SRqNVSto9Y+Ve/7GIkqt2srCDnZ1AT+Qj//Oo7LfR5AG7OfMPJ8kaG3Y+eix0Dr3jy9ygMJi5NGXKlOIIOP2CXw0Q3WWix63pCVv6d8sg0Fey7L///oU5ojtk1FDSXU6//OUvRd9t1HunSqut9q07VlqXvh9I702vu+66S9773vcWr3rRfl3atmM80GtjQjZU+JpLaiipoaLft+srcXpf++67r+huLTVfXn/9dbn11luL9xjpz9RsGWj3Tqc5tHZ7qRZ6/GHrHU/6c/3e/6yzzip2CLUunzmocajvbJo4cWLRr+rUbueS9j1jxgw57rjjCvOsdX3sYx8rTkrr/YodPUVNT0bT+aqhpCaljnPmmWcWOaev63n/+98fXEJJzSVNWHVH1c3Tm9X/drk0xsUAateXOrbqxuo2Qz3DsXXp0XVjxozpU4QDzUW3nKlzrNvsTj/99MJV9rn0gaLb09Qw0kvdR00YPe9Sk6b3pe6jPkD0QRJyNWlhyYIkROl4segSj23snn3rH63jKAJXuMYh0LlX3/qvYo5NGpM6Tq8WzMOY5/4MID/C8iNlNFqlpN0zVu71H4souWonCzvY2QmI6DuWTt37LPntr343YDdb7rGZHDrt32TRRO9gKvPVKoOZS+1uWt/joztK9Lti/c5YryuvvLJ4j46aIWqOqPmklx7JpoaRmit33HFH8V21vqJFN0X885//LPr5xCc+Ifq+JN1po+930l02uvtGDQWftu3mOtBrY0JeBeNjLv3pT3+SPfbYozgCTg02l0t3OeluJz06UI+da3d1moOaMnocXcucW2KJJYouXn31VdHc0d1kV1xxRcepDDSHQw89VPSdX6q3at3JXFJP4yc/+Ymst956hTmpGut/q9fys5/9TJZZZpliDmowqmmpu9bUAGtdf/jDHwpuEyZMKDbShF5JzSU9M1CNomuvvbbYMub6wjqNUcCWa+bMmaK7k3QHk56H2Lr0xVjq0qqDp06xbg3UI+zUzWuJoG21zYEHHlgUpJpU/Z3QVn+vvPJK4SD3dz71w1bPPFQHsfe7nb70pS+Jxqhx1dqO19pep/PU+YZcTVpYsiAJUTpeLLrEYxu7Z9/6R+s4isAVrnEIdO7Vt/6rmGOTxqSO06sF8zDmuT8DyI+w/EgZjVYpafeMlXv9xyJKrtrJwg52dgIiL7/wihz1+ePl4XsfHbCbj2y0uky99ihZarklQ4Zyjk25c2mgSW2yySbFd+YXXXRR0USPXjvggANk44037vO9sv7bj370o2IThO5AaR1xpu88UqPiz3/+cxGv36nr9+JqgOiOFzUWWqaIT9ve8+302phUO5e+853vyIUXXihqyKnJ4nLpO5P0XVSqs25GaXd1Mpd++tOfykknnVQYU2pQ9b6OPPLI4jhCNX46vd+p3Rz0iD01etQAUiNIr4HMpQceeKB455Pukvr0pz89fwq6c0p1VtPo6KOPLn6uu9x0s40emaf/1vvS+1STTE9qC72Smkuhk7XGqwunsHbffffiBWh33nln4fb2doJbovX+me42UlHUMFKXuP8LulZeeeViO51eakJpYuvLtFZZZZXiSDvddthyLPXcxN7JromjBpIabGo+/eMf/yjcRTW69HzEgV4G5sqgSQtLFiSuqqZthy5peZc5mm/9o3WZ9P+vL7jCNQ6Bzr361n8Vc2zSmNRxerVgHsY892cA+RGWHymj0Sol7Z6xcq//WETJVTtZ2MHOTqCeO5fKerWKcrHsXNI4NQHmzJlTmEp66Xty9L1KuiNJd730vtRc0e+i1ezRTRG9r8cee6w4Nk13OulOJf1O/Mknnyy+E+9/+bQd7LUxIa+Ccd25pLuuPvWpT8n73ve+YtOGz6Um3Uc/+lFRo6jd1WkOV199dWHcffOb3yxOOOt96St0dPdQ/5PT2o3Rfw477bSTLLXUUn0MKz39TE9+0/FUe/UL1J/QU+DU02i9O6t3/3pK2vDhw+d7EbqjS49P1Hl95jOf6TMVzQf1InTTTejVFeaSbvs7++yzC4NJzxbUl16paaSFrg6uXu3MJW3f6eg+LXh1//T6+9//XriGuqDTMXRMPX9RHed99tmncJ37X7/73e+K7XS6K0vNpC233FIOP/zwPucjWgVu0sKSBYlV5bhx6BKXb8zefesfreOoAVe4xiHQuVff+q9ijk0akzpOrxbMw5jn/gwgP8LyI2U0WqWk3TNW7vUfiyi5aicLO9jZCfRE1u2dS2W8WqXFxGIu6XuM9HvkD3/4w8WuHL30fUGbbrppYYi0ftYao3XEmp4Qpu8WGujSY/J0l4uaDyeffHJH2Tq1dXltTMirYFzNpeuvv744ZUwNHjV6XC89QlD57rjjjqJmULur0xz0aEF9nU27eJ2L8tGdS60Tytr1324O+j4k3VnW6VLT7iMf+Yjsvffexe4zNZeGDh3aJ0TnpieuqQmmV+uYRT0WUTe29L7UbFLDStuEXl1hLoVCamJ8kxaWLEjqmWHoUk9dXGblW/9o7ULVvw1c/Zm5RMC1MyXf+ndh3s1tyLf06sM8jHnuzwDyIyw/UkajVUraPWPlXv+xiJKrdrKwg52dQE/kqy++JhefdpVcOHXBo7nGHD1avnzwjrLEMouHDuMcH/pqld4DdTKX9Gg03YWiu0x6X9OmTSvMH92xsu+++87/Jz0WT4+002PXWq+Y0de76I4W3TShu09amyfa3axunNB3+VxyySWy9tprd+QxUFvX18Zo59ZXwbiaS3rc22233Sa33nprsbGj/6UGjh79p6eK9b7+/d//Xc4///ziSDw9Gq/dNdgcdt11V1ED7YYbbih2TumlO8TU2FG22r9ePnPQ+9B3aPW+Hn744WJDim6O2WyzzQqDUY/bU6NIX7vzn//5n4VZ2Lp0DTBmzBjZeeed52+EefHFF2WrrbYqjkxUY7LFQ3dXKcP99tuvMOlCr1qYS/quo7vvvrvY8aNb2/pfWiD7779/6L12VXyTFpYsSOqZmuhST11cZuVb/2jtQtW/DVz9mblEwLUzJd/6d2HezW3It/TqwzyMee7PAPIjLD9SRqNVSto9Y+Ve/7GIkqt2srCDnZ3A/0WqwfTMY8/Jr06+Up6a/aysuPoKsvvhX5L3r7p8UmOpNSPrq1U0/p577in+6KVGzssvv1zsNNFrxRVXLL7410tPzzr44IMLc0hfuaLfe+vP1EBS80iNgN67X3THjO480Xbjx4+XESNGFKdz6etYfvjDH8q22247H6geo7bkkkvKGmusIXPnzi2O19M56TuBvv71r/eRzLWtz2tjdACfV8HcfPPN0nqW6Oti1EBp7bJRFmr29L70KDfddaPvJlIzrt2lbPRdVNtvv72stNJKhdegBo7OS4/TO+ecc/oYTz5z0M9a1UBNLX33kV6q1wsvvFCYPuuvv37xM9859L+Pgd65pCen6Tuf3nzzTdljjz2K3Uz6swsuuGB+3vXexdZ6T5TujlIz6tlnny3mueyyyxZz1OP4Qq9KzaV58+aJuoZ6/p9u/dMi0Z+1rtZ/6996dByXO4EmLSxZkLjrmrIluqSkXe5YvvWP1uXyb/UGV7jGIdC5V9/6r2KOTRqTOk6vFszDmOf+DCA/wvIjZTRapaTdM1bu9R+LKLlqJws72NkJLBj55utvyltvvC3DFxsmi45YtMyuvfqyvlpFB/nRj35UmBrtLn3Pjn6pr9cTTzxRHMv24IMPFqaEmkBqgujuF321ippH/a/Zs2cXMWoUqVmiu2R0I4Yed9f7+vnPfy6XXnqpPP7448Wxaeuss4584xvfkC222GKBPl3b+rw2pjWI66tg1PTSHVntrt6vo2n9u875hBNOKN4/peZcu+vPf/6znHHGGaI70XT3jvoKq666anGc3YQJE4pdY70v3zno563uHHrggQeKbjbYYAM56KCDZL311pvfre8c+t/HQOaStlNtzzzzzCIXnnnmmSJfNL8mTZpUmE39L9VPuf31r38t2upreXR3XLtdX17F8r+NKzWX9Mb0nUX6Eil9OZk6byryF7/4xeKMQnUS9YVkRxxxRFFkXO4EmrSwZEHirmvKluiSkna5Y/nWP1qXy7/VG1zhGodA515967+KOTZpTOo4vVowD2Oe+zOA/AjLj5TRaJWSds9Yudd/LKLkqp0s7GBnJ0AkBCCQC4FKzSV1DNUtVAdNL93uNnHixOKPXuqo6VmGev6fmk5c7gSatLBkQeKua8qW6JKSdrlj+dY/WpfLv9UbXOEah0DnXn3rv4o5NmlM6ji9WjAPY577M4D8CMuPlNFolZJ2z1i5138souSqnSzsYGcnQCQEIJALgUrNJT2HUHctHX300QXPtdZaq9j+p+dOti49/1GPxLv22mtzYZ7kPpq0sGRBkiQlvAdBF29ktQnwrX+0jiMdXOEahwDmUkqu1HFK2j1jwTyMue8aIGy09NHkR3rm1hHRykrOHpd7/dvJdI4kV+1kYQc7OwEiIQCBXAhUai7peYC6M2ny5MkFT/1vfQnZ8ccfP5/v9773PTn//PPlvvvuy4V5kvto0sKSBUmSlPAeBF28kdUmwLf+0TqOdHCFaxwCnXv1rf8q5tikManj9GrBPIx57s8A8iMsP1JGo1VK2j1j5V7/sYiSq3aysIOdnQCREIBALgQqNZf0HUv6LqUf/vCHBc9x48bJU089Jdddd50MGzas+JmaT6+88orMmDEjF+ZJ7qNJC0sWJElSwnsQdPFGVpsA3/pH6zjSwRWucQhgLqXkSh2npN0zFszDmPuuAcJGSx9NfqRnbh0Rrazk7HG517+dTOdIctVOFnawsxMgEgIQyIVApebSKaecIr/61a/kjjvuKMykq666Sg4//HBZe+21ZdNNN5U//elP8sc//lG+9a1vFe9d4nIn0KSWx22lAAAgAElEQVSFJQsSd11TtkSXlLTLHcu3/tG6XP6t3uAK1zgEOvfqW/9VzLFJY1LH6dWCeRjz3J8B5EdYfqSMRquUtHvGyr3+YxElV+1kYQc7OwEiIQCBXAhUai49+eSTcvvttxdH4S277LIF0x//+Mfyk5/8RF5//XUZPny47LbbbsWxeYssskguzJPcR5MWlixIkqSE9yDo4o2sNgG+9Y/WcaSDK1zjEMBcSsmVOk5Ju2csmIcx910DhI2WPpr8SM/cOiJaWcnZ43KvfzuZzpHkqp0s7GBnJ0AkBCCQC4FKzaWBIM6dO1f++c9/yjLLLCNDhw7NhXXS+2jSwpIFSdLUcB4MXZxR1a6hb/2jdRwJ4QrXOAQwl1JypY5T0sZcKoO27xqgjDFT9kFNpqQdNhZahfGzROde/xYmLjHkqgul9m1gBzs7ASIhAIFcCFRqLk2ZMkVGjRolEyZMyIVnbe6jSQtLFiS1SZs+E0GXeuriMivf+kdrF6r+beDqz8wlAq6YSy55UlYb8q0sku79wNydVbuWvmuAsNHSR5Mf6ZlbR0QrKzl7XO71byfTOZJctZOFHezsBIiEAARyIVCpubT++uvL+PHj5dBDD82FZ23uo0kLSxYktUkbzKV6SuE9K9/6pwa9ETsFwNUJk3cjuGIueSdNQAD5FgDPGApzI7j/DfNdA4SNlj6a/EjP3DoiWlnJ2eNyr387Gcwl2MUiYO+XZ6SdHZEQgEC9CFRqLo0ePVpWW201OfXUU+tFJYPZNGlhyYdqPRMOXeqpi8usfOsfrV2o+reBqz8zlwi4Yi655ElZbci3ski69wNzd1btWvquAcJGSx9NfqRnbh0Rrazk7HG517+dDOYS7GIRsPfLM9LOjkgIQKBeBCo1l6699lrRo/GmT58u6623Xr3INHw2TVpY8qFaz2RDl3rq4jIr3/pHaxeq/m3g6s/MJQKumEsueVJWG/KtLJLu/cDcnRXmUhgrouMSoJbj8u3G+o9FlFy1k4Ud7OwEiIQABHIhUKm5dMUVV8hVV10ld911l2y//fay1lprybLLLitDhgxZgO/OO++cC/Mk9+H75XKSSQ0wCAuSKukPPDa61FMXl1n51j9au1D1bwNXf2YuEXDtTMm3/l2Yd3Mb8i29+jAPY577M4D8CMuPlNFolZJ2z1i5138souSqnSzsYGcnQCQEIJALgeTmku5U2mabbWTrrbeWUaNGFUbSvHnz+vDsbS7pv+l/z5o1KxfmSe6jSQtLFiRJUsJ7EHTxRlabAN/6R+s40sEVrnEIYC6l5Eodp6TdMxbMw5j7rgHCRksfTX6kZ24dEa2s5Oxxude/nUznSHLVThZ2sLMTIBICEMiFQHJzSQ2liRMnFn8uu+yytruU2sHdZZddcmGe5D6atLBkQZIkJbwHQRdvZLUJ8K1/tI4jHVzhGocA5lJKrtRxStqYS2XQ9l0DlDFmyj6oyZS0w8ZCqzB+lujc69/CxCWGXHWh1L4N7GBnJ0AkBCCQC4FKzaVcINbxPpq0sGRBUscM4v85XE9V3GblW//UoBtX31Zw9SXm1h6umEtumVJOK/KtHI4+vcDch9aCbX3XAGGjpY8mP9Izt46IVlZy9rjc699OpnMkuWonCzvY2QkQCQEI5EIAcykXJfvdR5MWlixI6pmE6FJPXVxm5Vv/aO1C1b8NXP2ZuUTAFXPJJU/KakO+lUXSvR+Yu7Nq19J3DRA2Wvpo8iM9c+uIaGUlZ4/Lvf7tZDCXYBeLgL1fnpF2dkRCAAL1IoC5VC89SptNkxaWfKiWJnupHaFLqTiTduZb/2gdRx64wjUOAcyllFyp45S0e8aCeRhz3zVA2Gjpo8mP9MytI6KVlZw9Lvf6t5PBXIJdLAL2fnlG2tkRCQEI1ItAJebSyJEjRf+4XkOGDJFzzz3XtTntRKRJC0s+VOuZsuhST11cZuVb/2jtQtW/DVz9mblEwLUzJd/6d2HezW3It/TqwzyMee7PAPIjLD9SRqNVSto9Y+Ve/7GIkqt2srCDnZ0AkRCAQC4EKjGXfOGpuTRr1izfsK5u36SFJQuSeqYqutRTF5dZ+dY/WrtQ9W8DV39mLhFwxVxyyZOy2pBvZZF07wfm7qzatfRdA4SNlj6a/EjP3DoiWlnJ2eNyr387mc6R5KqdLOxgZydAJAQgkAuBSsylCRMmyPjx470Y+ux08uo408ZNWliyIKlnEqJLPXVxmZVv/aO1C1X/NnD1Z+YSAVfMJZc8KasN+VYWSfd+YO7OCnMpjBXRcQlQy3H5dmP9xyJKrtrJwg52dgJEQgACuRCoxFyaOHGi6B+ueAR8v1yON5PBe2ZBMjijKlqgSxXUyxnTt/7Ruhzu/XuBK1zjEOjcq2/9VzHHJo1JHadXC+ZhzHN/BpAfYfmRMhqtUtLuGSv3+o9FlFy1k4Ud7OwEiIQABHIhgLmUi5L97qNJC0sWJPVMQnSppy4us/Ktf7R2oerfBq7+zFwi4Iq55JInZbUh38oi6d4PzN1ZtWvpuwYIGy19NPmRnrl1RLSykrPH5V7/djKdI8lVO1nYwc5OgEgIQCAXAphLuSiJuZSpktXdFgvF6tiHjuz7iyVahxJvHw9XuMYhgLmUkit1nJJ2z1gwD2PuuwYIGy19NPmRnrl1RLSykrPH5V7/djKYS7CLRcDeL89IOzsiIQCBehHAXKqXHqXNpkkLSz5US5O91I7QpVScSTvzrX+0jiMPXOEahwDmUkqu1HFK2phLZdD2XQOUMWbKPqjJlLTDxkKrMH6W6Nzr38LEJYZcdaHUvg3sYGcnQCQEIJALgeTmUi7g6n4fTVpYsiCpZzahSz11cZmVb/2jtQtV/zZw9WfmEgFXzCWXPCmrDflWFkn3fmDuzqpdS981QNho6aPJj/TMrSOilZWcPS73+reT6RxJrtrJwg52dgJEQgACuRDAXMpFyX730aSFJQuSeiYhutRTF5dZ+dY/WrtQ9W8DV39mLhFwxVxyyZOy2pBvZZF07wfm7qwwl8JYER2XALUcl2831n8souSqnSzsYGcnQCQEIJALAcylXJTEXMpUyepui4VidexDR8ZcCiVYTjw1VA7H/r3AFXMpTma175V8S0m7ZyyYhzH3XQOEjZY+mvxIz9w6IlpZydnjcq9/O5nOkeSqnSzsYGcnQCQEIJALAcylXJTsdx9NWliyIKlnEqJLPXVxmZVv/aO1C1X/NnD1Z+YSAVfMJZc8KasN+VYWSfd+YO7Oql1L3zVA2Gjpo8mP9MytI6KVlZw9Lvf6t5PBXIJdLAL2fnlG2tkRCQEI1IsA5lK99ChtNk1aWPKhWprspXaELqXiTNqZb/2jdRx54ArXOAQwl1JypY5T0u4ZC+ZhzH3XAGGjpY8mP9Izt46IVlZy9rjc699OBnMJdrEI2PvlGWlnRyQEIFAvAphL9dKjtNk0aWHJh2ppspfaEbqUijNpZ771j9Zx5IErXOMQwFxKyZU6Tkkbc6kM2r5rgDLGTNkHNZmSdthYaBXGzxKde/1bmLjEkKsulNq3gR3s7ASIhAAEciGAuZSLkv3uo0kLSxYk9UxCdKmnLi6z8q1/tHah6t8Grv7MXCLgirnkkidltSHfyiLp3g/M3Vm1a+m7BggbLX00+ZGeuXVEtLKSs8flXv92Mp0jyVU7WdjBzk6ASAhAIBcCmEu5KIm5lKmS1d0WC8Xq2IeO7PuLJVqHEm8fD1e4xiGAuZSSK3WcknbPWDAPY+67BggbLX00+ZGeuXVEtLKSs8flXv92MphLsItFwN4vz0g7OyIhAIF6EcBcqpcepc2mSQtLPlRLk73UjtClVJxJO/Otf7SOIw9c4RqHAOZSSq7UcUramEtl0PZdA5QxZso+qMmUtMPGQqswfpbo3OvfwsQlhlx1odS+DexgZydAJAQgkAsBzKVclOx3H01aWLIgqWcSoks9dXGZlW/9o7ULVf82cPVn5hIBV8wllzwpqw35VhZJ935g7s6qXUvfNUDYaOmjyY/0zK0jopWVnD0u9/q3k+kcSa7aycIOdnYCREIAArkQwFzKRUnMpUyVrO62WChWxz50ZN9fLNE6lHj7eLjCNQ4BzKWUXKnjlLR7xoJ5GHPfNUDYaOmjyY/0zK0jopWVnD0u9/q3k8Fcgl0sAvZ+eUba2REJAQjUiwDmUr30KG02TVpY8qFamuyldoQupeJM2plv/aN1HHngCtc4BDCXUnKljlPSxlwqg7bvGqCMMVP2QU2mpB02FlqF8bNE517/FiYuMeSqC6X2bWAHOzsBIiEAgVwIYC7lomS/+2jSwpIFST2TEF3qqYvLrHzrH61dqPq3gas/M5cIuGIuueRJWW3It7JIuvcDc3dW7Vr6rgHCRksfTX6kZ24dEa2s5Oxxude/nUznSHLVThZ2sLMTIBICEMiFAOZSLkpiLmWqZHW3xUKxOvahI/v+YonWocTbx8MVrnEIYC6l5Eodp6TdMxbMw5j7rgHCRksfTX6kZ24dEa2s5Oxxude/nQzmEuxiEbD3yzPSzo5ICECgXgQwl+qlR2mzadLCkg/V0mQvtSN0KRVn0s586x+t48gDV7jGIYC5lJIrdZySNuZSGbR91wBljJmyD2oyJe2wsdAqjJ8lOvf6tzBxiSFXXSi1bwM72NkJEAkBCORCAHMpFyX73UeTFpYsSOqZhOhST11cZuVb/2jtQtW/DVz9mblEwBVzySVPympDvpVF0r0fmLuzatfSdw0QNlr6aPIjPXPriGhlJWePy73+7WQ6R5KrdrKwg52dAJEQgEAuBLI3l9555x05++yz5dJLL5Xnn39eRo4cKXvuuaeMHTtWhgwZMqCOL730klx++eVyyy23yOzZs2XOnDnywQ9+UHbYYQeZMGGCDB8+fH6sT1sNGjdunNx9990LjL3QQgvJzJkzS8mtJi0sWZCUInnpnaBL6UiTdehb/2gdRxq4wjUOgc69+tZ/FXNs0pjUcXq1YB7GPPdnAPkRlh8po9EqJe2esXKv/1hEyVU7WdjBzk6ASAhAIBcC2ZtLxxxzjFx88cWy++67y3rrrSd33HGHXH/99TJp0iSZOHHigDqqqbT//vvL5ptvLptuuqksvvjics8998jVV18tG264oUyfPl3UDNLLp23LXHrooYfk6KOP7jP+0KFDZccddywlt5q0sGRBUorkpXeCLqUjTdahb/2jdRxp4ArXOAQwl1JypY5T0u4ZC+ZhzH3XAGGjpY8mP9Izt46IVlZy9rjc699OpnMkuWonCzvY2QkQCQEI5EIga3Np1qxZsvPOO8tee+0lkydPnq/ZQQcdJDfddFPxZ/nll2+r5eOPP178XHcr9b5+8IMfyJlnnimnn366bLvttsU/+bTV9rpz6W9/+5vcdttt0fKoSQtLFiTR0iCoY3QJwldpsG/9o3UcueAK1zgEOvfqW/9VzLFJY1LH6dWCeRjz3J8B5EdYfqSMRquUtHvGyr3+YxElV+1kYQc7OwEiIQCBXAhkbS6ddtppxZF4urNoxRVXnK/ZvffeK2PGjJFjjz22+Nvn0h1HO+20kxx44IGy3377dQwdqG3LXNJ5vfHGGzJixIiOR/T5zK/VtkkLSxYkFoXjx6BLfMaxRvCtf7SOowRc4RqHAOZSSq7UcUraPWPBPIy57xogbLT00eRHeubWEdHKSs4el3v928l0jiRX7WRhBzs7ASIhAIFcCGRtLumOpYcffrg4Cq/39fbbb8v6668vo0ePluOPP95Ly9tvv12+8Y1vyHHHHSdf+cpXOsYO1FbNJV34LbzwwvLmm28WR+5tt912cthhh8myyy7rNZ+BGmv/8+bNK4yrul9qsOm12GKL1X2qXTU/dKlO7lGjRgUN7lv/aB2Ee8BguMLVQiB1/Vvm2E0x1HF6tbuZeWj9q1q+a4D0CoeN2M35EUYufTRa+TMPfQbkXv/+RN0iyFU3Tu1awa48dqH1b58JkRCAAATCCGRtLu2www4ybNgwueyyyxag9MlPflLWWWcdmTZtmjPBd999VyZMmCAPPPCA3HjjjbLccssNGNup7ZQpU2SFFVaQNddcszCAfv/738sll1xSHMGnfy+55JLOcxqoYZMWlixIguWO0gG6RMHq1GnowtK3/tHaSRbvRnD1RuYUkDvX1PXvBL2LG+Web3WUtpuZh9a/6um7BqhjDnSaUzfnB1o1jYD/fEOfAbnXvz9RtwieK26c2rWCXXnsQuvfPhMiIQABCIQRyNpc2mabbQoD6KKLLlqA0pZbblmYOdOnT3cm2Dpm75hjjinem9Tp8mmr/Vx88cWi/U6cOFEmTZrkPKeBGjZpSzxbqYPljtIBukTBmqRT3/pH6ziywBWucQh07tW3/quYY5PGpI7TqwXzMOa5PwPIj7D8SBmNVilp94yVe/3HIkqu2snCDnZ2AkRCAAK5EMjaXCpz59J5550n3/3ud4uj8PRIvE6XT9ve/WyyySay2mqrtTXDfBOuSQtLFiS+6qZpjy5pOMcYxbf+0TqGCrw3JA5VuA7G1bf+B+uv2/+d52P6DIB5GPPcnwHkR1h+pIxGq5S0MZdCaJOrdnqwg52dAJEQgEAuBLI2lwZ759Iuu+wiU6dOHVRLPVbvqKOOki984QtyyimnyNChQweM8WnbvxOdz5w5c+Q3v/nNoHMarEGTfrFkQTKYmtX8O7pUw72MUX3rH63LoL5gH3CFaxwCnXv1rf8q5tikManj9GrBPIx57s8A8iMsP1JGo1VK2phLIbTJVTs92MHOToBICEAgFwJZm0unnnqqnHPOOXLLLbfIiiuuOF+ze++9V8aMGSPf/va3ZezYsR21vOaaa+Swww6TLbbYQk4//XRZeOGFB2zv07Z/J/qOJt259OEPf1guvPDC4Pxq0i+WLEiC5Y7SAbpEwZqkU9/6R+s4ssAVrnEIYC6l5Eodp6TdMxbMw5j7rgHCRksfTX6kZ24dEa2s5Oxxude/nUznSHLVThZ2sLMTIBICEMiFQNbm0syZM0V3A+kOpsmTJ8/X7KCDDpIbb7xRbrrpJllhhRVEX0L41FNPydJLLy3LLLPM/Hba5sADD5RPfOIThUk1bNiwAXV3bfvaa6/JIossIsOHD+/T17Rp0+Tkk0+WQw45RPbdd9/g/GrSwpIFSbDcUTpAlyhYk3TqW/9oHUcWuMI1DoHOvfrWfxVzbNKY1HF6tWAexjz3ZwD5EZYfKaPRKiXtnrFyr/9YRMlVO1nYwc5OgEgIQCAXAlmbSyqSHmenR9Xtvvvusu6668qdd94p1113nUycOFEmTZpU6HjXXXfJ+PHj+/zs/vvvL3Y1qRF0xBFHyGKLLdZH85VXXlk22GCD4mc+bXWsgw8+uDhiT/sYMmRIMf6MGTNk1KhRxa6l97znPcH51aSFJQuSYLmjdIAuUbAm6dS3/tE6jixwhWscAp179a3/KubYpDGp4/RqwTyMee7PAPIjLD9SRqNVSto9Y+Ve/7GIkqt2srCDnZ0AkRCAQC4EsjeX/vWvf8nZZ59dGEzPPfecjBw5sjCNxo0bVxg7A5lL2n7KlCkD6qw7ok488cTi333aPvHEE8V7mx588EF54YUXZO7cubLSSivJdtttJ/vss4+MGDGilNxq0sKSBUkpkpfeCbqUjjRZh771j9ZxpIErXOMQwFxKyZU6Tkm7ZyyYhzH3XQOEjZY+mvxIz9w6IlpZydnjcq9/O5nOkeSqnSzsYGcnQCQEIJALgezNpVyE8r2PJi0sWZD4qpumPbqk4RxjFN/6R+sYKvAFaRyqcB2Mq2/9D9Zft/87z8f0GQDzMOa5PwPIj7D8SBmNVilp94yVe/3HIkqu2snCDnZ2AkRCAAK5EMBcykXJfvfRpIUlC5J6JiG61FMXl1n51j9au1D1bwNXf2YuEXDtTMm3/l2Yd3Mb8i29+jAPY577M4D8CMuPlNFolZI25lIIbXLVTg92sLMTIBICEMiFAOZSLkpiLmWqZHW3xUKxOvahI/t+sYTWocTbx8MVrnEIYC6l5Eodp6TdMxbMw5j7rgHCRksfTX6kZ24dEa2s5Oxxude/nUznSHLVThZ2sLMTIBICEMiFAOZSLkpiLmWqZHW3xUKxOvahI/v+YonWocQxl+IQhKuFq2/9W8bophiej+nVhnkY89yfAeRHWH6kjEarlLR7xsq9/mMRJVftZGEHOzsBIiEAgVwIYC7loiTmUqZKVndbLBSrYx86su8vlmgdShwTJA5BuFq4+ta/ZYxuiuH5mF5tmIcxz/0ZQH6E5UfKaLRKSRtzKYQ2uWqnBzvY2QkQCQEI5EIAcykXJTGXMlWyuttioVgd+9CRfb9YQutQ4pggcQjC1cLVt/4tY3RTDM/H9GrDPIx57s8A8iMsP1JGo1VK2phLIbTJVTs92MHOToBICEAgFwKYS7koibmUqZLV3RYLxerYh47s+8USWocSxwSJQxCuFq6+9W8Zo5tieD6mVxvmYcxzfwaQH2H5kTIarVLSxlwKoU2u2unBDnZ2AkRCAAK5EMBcykVJzKVMlazutlgoVsc+dGTfL5bQOpQ4JkgcgnC1cPWtf8sY3RTD8zG92jAPY577M4D8CMuPlNFolZI25lIIbXLVTg92sLMTIBICEMiFAOZSLkpiLmWqZHW3xUKxOvahI/t+sYTWocQxQeIQhKuFq2/9W8bophiej+nVhnkY89yfAeRHWH6kjEarlLQxl0Jok6t2erCDnZ0AkRCAQC4EMJdyURJzKVMlq7stForVsQ8d2feLJbQOJY4JEocgXC1cfevfMkY3xfB8TK82zMOY5/4MID/C8iNlNFqlpI25FEKbXLXTgx3s7ASIhAAEciGAuZSLkphLmSpZ3W2xUKyOfejIvl8soXUocUyQOAThauHqW/+WMbophudjerVhHsY892cA+RGWHymj0SolbcylENrkqp0e7GBnJ0AkBCCQCwHMpVyUxFzKVMnqbouFYnXsQ0f2/WIJrUOJY4LEIQhXC1ff+reM0U0xPB/Tqw3zMOa5PwPIj7D8SBmNVilpYy6F0CZX7fRgBzs7ASIhAIFcCGAu5aIk5lKmSlZ3WywUq2MfOrLvF0toHUocEyQOQbhauPrWv2WMborh+ZhebZiHMc/9GUB+hOVHymi0SkkbcymENrlqpwc72NkJEAkBCORCAHMpFyUxlzJVsrrbYqFYHfvQkX2/WELrUOKYIHEIwtXC1bf+LWN0UwzPx/RqwzyMee7PAPIjLD9SRqNVStqYSyG0yVU7PdjBzk6ASAhAIBcCmEu5KIm5lKmS1d0WC8Xq2IeO7PvFElqHEscEiUMQrhauvvVvGaObYng+plcb5mHMc38GkB9h+ZEyGq1S0sZcCqFNrtrpwQ52dgJEQgACuRDAXMpFScylTJWs7rZYKFbHPnRk3y+W0DqUOCZIHIJwtXD1rX/LGN0Uw/MxvdowD2Oe+zOA/AjLj5TRaJWSNuZSCG1y1U4PdrCzEyASAhDIhQDmUi5KYi5lqmR1t8VCsTr2oSP7frGE1qHEMUHiEISrhatv/VvG6KYYno/p1YZ5GPPcnwHkR1h+pIxGq5S0MZdCaJOrdnqwg52dAJEQgEAuBDCXclEScylTJau7LRaK1bEPHdn3iyW0DiWOCRKHIFwtXH3r3zJGN8XwfEyvNszDmOf+DCA/wvIjZTRapaSNuRRCm1y104Md7OwEiIQABHIhgLmUi5KYS5kqWd1tsVCsjn3oyL5fLKF1KHFMkDgE4Wrh6lv/ljG6KYbnY3q1YR7GPPdnAPkRlh8po9EqJW3MpRDa5KqdHuxgZydAJAQgkAsBzKVclMRcylTJ6m6LhWJ17ENH9v1iCa1DiWOCxCEIVwtX3/q3jNFNMTwf06sN8zDmuT8DyI+w/EgZjVYpaWMuhdAmV+30YAc7OwEiIQCBXAhgLuWiJOZSpkpWd1ssFKtjHzqy7xdLaB1KHBMkDkG4Wrj61r9ljG6K4fmYXm2YhzHP/RlAfoTlR8potEpJG3MphDa5aqcHO9jZCRAJAQjkQgBzKRclMZcyVbK622KhWB370JF9v1hC61DimCBxCMLVwtW3/i1jdFMMz8f0asM8jHnuzwDyIyw/UkajVUramEshtMlVOz3Ywc5OgEgIQCAXAphLuSiJuZSpktXdFgvF6tiHjuz7xRJahxLHBIlDEK4Wrr71bxmjm2J4PqZXG+ZhzHN/BpAfYfmRMhqtUtLGXAqhTa7a6cEOdnYCREIAArkQwFzKRUnMpUyVrO62WChWxz50ZN8vltA6lDgmSByCcLVw9a1/yxjdFMPzMb3aMA9jnvszgPwIy4+U0WiVkjbmUghtctVOD3awsxMgEgIQyIUA5lIuSmIuZapkdbfFQrE69qEj+36xhNahxDFB4hCEq4Wrb/1bxuimGJ6P6dWGeRjz3J8B5EdYfqSMRquUtDGXQmiTq3Z6sIOdnQCREIBALgQwl3JREnMpUyWruy0WitWxDx3Z94sltA4ljgkShyBcLVx9698yRjfF8HxMrzbMw5jn/gwgP8LyI2U0WqWkjbkUQptctdODHezsBIiEAARyIYC5lIuSmEuZKlndbbFQrI596Mi+XyyhdShxTJA4BOFq4epb/5YxuimG52N6tWEexjz3ZwD5EZYfKaPRKiVtzKUQ2uSqnR7sYGcnQCQEIJALAcylXJTEXMpUyepui4VidexDR/b9YgmtQ4ljgsQhCFcLV9/6t4zRTTE8H9OrDfMw5rk/A8iPsPxIGY1WKWljLoXQJlft9GAHOzsBIiEAgVwIYC7loiTmUqZKVndbLBSrYx86su8XS2gdShwTJA5BuFq4+ta/ZYxuiuH5mF5tmIcxz/0ZQH6E5UfKaLRKSRtzKYQ2uWqnBzvY2QkQCQEI5EIAcykXJTGXMlWyuttioVgd+9CRfb9YQutQ4pggcQjC1cLVt/4tY3RTDM/H9GrDPIx57s8A8i0Aoo0AACAASURBVCMsP1JGo1VK2phLIbTJVTs92MHOToBICEAgFwKYS7koibmUqZLV3RYLxerYh47s+8USWocSxwSJQxCuFq6+9W8Zo5tieD6mVxvmYcxzfwaQH2H5kTIarVLSxlwKoU2u2unBDnZ2AkRCAAK5EMBcykVJzKVMlazutlgoVsc+dGTfL5bQOpQ4JkgcgnC1cPWtf8sY3RTD8zG92jAPY577M4D8CMuPlNFolZI25lIIbXLVTg92sLMTIBICEMiFAOZSLkpiLmWqZHW3xUKxOvahI/t+sYTWocQxQeIQhKuFq2/9W8bophiej+nVhnkY89yfAeRHWH6kjEarlLQxl0Jok6t2erCDnZ0AkRCAQC4EMJdyURJzKVMlq7stForVsQ8d2feLJbQOJY4JEocgXC1cfevfMkY3xfB8TK82zMOY5/4MID/C8iNlNFqlpI25FEKbXLXTgx3s7ASIhAAEciGAuZSLkphLmSpZ3W2xUKyOfejIvl8soXUocUyQOAThauHqW/+WMbophudjerVhHsY892cA+RGWHymj0SolbcylENrkqp0e7GBnJ0AkBCCQCwHMpVyUxFzKVMnqbouFYnXsQ0f2/WIJrUOJY4LEIQhXC1ff+reM0U0xPB/Tqw3zMOa5PwPIj7D8SBmNVilpYy6F0CZX7fRgBzs7ASIhAIFcCGAu5aIk5lKmSlZ3WywUq2MfOrLvF0toHUocEyQOQbhauPrWv2WMborh+ZhebZiHMc/9GUB+hOVHymi0SkkbcymENrlqpwc72NkJEAkBCORCAHMpFyX73ce9995b/GTIkCG1v8N58+Y1Zq61h1niBNGlRJieXWndbrDBBp5R/9fct/7R2oy6YyBc4WohkLr+LXPsphjqOL3a3cw8tP5VLd81QHqFw0bs5vwII5c+Gq38mYc+A3Kvf3+ibhHkqhundq1gVx670Pq3z4RICEAAAmEEMJfC+NU2moVlbaVhYhAYlEDowpL6HxQxDSBQWwLUf22lYWIQiE4gtP67wVyKLgIDQKBCAqHPAH4HqFA8hoZAIIHQ+g8cnnAIQAACZgKYS2Z0BEIAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIACB7iOAudR9mnPHEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQMBMAHPJjI5ACEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEINB9BDCXuk9z7hgCEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEImAlgLpnREQgBCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEuo8A5lL3ac4dQwACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEzAcwlMzoCIQABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgED3EcBc6j7NuWMIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgYCaAuWRGRyAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQ6D4CmEvdpzl3DAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAATMBDCXzOgIhAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAALdRwBzqfs0544hAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAgJkA5pIZHYEQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAoPsIYC51n+bcMQQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhAwE8BcMqMjEAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQh0HwHMpe7TnDuGAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAmYCmEtmdARCAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAge4jgLnUfZpzxxCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEDATABzyYyOQAhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCDQfQQwl7pPc+4YAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCJgJYC6Z0REIAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABLqPAOZSppr/8Y9/LO5sgw02yPQOuS0IQGAgAtQ/uQGB7iVA/Xev9tw5BJQAzwDyAALdS4D6717tuXMIQAACEIBAVQQwl6oiH3nc//mf/ylG2HDDDSOPFN79//t//6/oZNSoUeGd0UNpBNClNJTJO/Ktf7SOIxFc4RqHQOdefeu/ijk2aUzqOL1aMA9jnvszgPwIy4+U0WiVknbPWLnXfyyi5KqdLOxgZydAJAQgkAsBzKVclOx3H01aWLIgqWcSoks9dXGZlW/9o7ULVf82cPVn5hIBV8wllzwpqw35VhZJ935g7s6qXUvfNUDYaOmjyY/0zK0jopWVnD0u9/q3k+kcSa7aycIOdnYCREIAArkQwFzKRUnMpUyVrO62WChWxz50ZN9fLNE6lHj7eLjCNQ4BzKWUXKnjlLR7xoJ5GHPfNUDYaOmjyY/0zK0jopWVnD0u9/q3k8Fcgl0sAvZ+eUba2REJAQjUiwDmUr30KG02TVpY8qFamuyldoQupeJM2plv/aN1HHngCtc4BDCXUnKljlPSxlwqg7bvGqCMMVP2QU2mpB02FlqF8bNE517/FiYuMeSqC6X2bWAHOzsBIiEAgVwIYC7lomS/+2jSwpIFST2TEF3qqYvLrHzrH61dqPq3gas/M5cIuGIuueRJWW3It7JIuvcDc3dW7Vr6rgHCRksfTX6kZ24dEa2s5Oxxude/nUznSHLVThZ2sLMTIBICEMiFAOZSLkpiLmWqZHW3xUKxOvahI/v+YonWocTbx8MVrnEIYC6l5Eodp6TdMxbMw5j7rgHCRksfTX6kZ24dEa2s5Oxxude/nQzmEuxiEbD3yzPSzo5ICECgXgQwl+qlR2mzadLCkg/V0mQvtSN0KRVn0s586x+t48gDV7jGIYC5lJIrdZySNuZSGbR91wBljJmyD2oyJe2wsdAqjJ8lOvf6tzBxiSFXXSi1bwM72NkJEAkBCORCAHMpFyX73UeTFpYsSOqZhOhST11cZuVb/2jtQtW/DVz9mblEwBVzySVPympDvpVF0r0fmLuzatfSdw0QNlr6aPIjPXPriGhlJWePy73+7WQ6R5KrdrKwg52dAJEQgEAuBDCXclEScylTJau7LRaK1bEPHdn3F0u0DiXePh6ucI1DAHMpJVfqOCXtnrFgHsbcdw0QNlr6aPIjPXPriGhlJWePy73+7WQwl2AXi4C9X56RdnZEQgAC9SKAuVQvPUqbTZMWlnyoliZ7qR2hS6k4k3bmW/9oHUceuMI1DgHMpZRcqeOUtDGXyqDtuwYoY8yUfVCTKWmHjYVWYfws0bnXv4WJSwy56kKpfRvYwc5OgEgIQCAXAphLuSjZ7z6atLBkQVLPJESXeuriMivf+kdrF6r+beDqz8wlAq6YSy55UlYb8q0sku79wNydVbuWvmuAsNHSR5Mf6ZlbR0QrKzl7XO71byfTOZJctZOFHezsBIiEAARyIYC5lIuSmEuZKlndbbFQrI596Mi+v1iidSjx9vFwhWscAphLKblSxylp94wF8zDmvmuAsNHSR5Mf6ZlbR0QrKzl7XO71byeDuQS7WATs/fKMtLMjEgIQqBcBzKV66VHabJq0sORDtTTZS+0IXUrFmbQz3/pH6zjywBWucQhgLqXkSh2npI25VAZt3zVAGWOm7IOaTEk7bCy0CuNnic69/i1MXGLIVRdK7dvADnZ2AkRCAAK5EMBcykXJfvfRpIUlC5J6JiG61FMXl1n51j9au1D1bwNXf2YuEXDFXHLJk7LakG9lkXTvB+burNq19F0DhI2WPpr8SM/cOiJaWcnZ43KvfzuZzpHkqp0s7GBnJ0AkBCCQCwHMpVyUxFzKVMnqbouFYnXsQ0f2/cUSrUOJt4+HK1zjEMBcSsmVOk5Ju2csmIcx910DhI2WPpr8SM/cOiJaWcnZ43KvfzsZzCXYxSJg75dnpJ0dkRCAQL0IYC7VS4/SZtOkhSUfqqXJXmpH6FIqzqSd+dY/WseRB65wjUMAcyklV+o4JW3MpTJo+64ByhgzZR/UZEraYWOhVRg/S3Tu9W9h4hJDrrpQat8GdrCzEyASAhDIhQDmUi5K9ruPJi0sWZDUMwnRpZ66uMzKt/7R2oWqfxu4+jNziYAr5pJLnpTVhnwri6R7PzB3Z9Wupe8aIGy09NHkR3rm1hHRykrOHpd7/dvJdI4kV+1kYQc7OwEiIQCBXAhgLuWiJOZSpkpWd1ssFKtjHzqy7y+WaB1KvH08XOEahwDmUkqu1HFK2j1jwTyMue8aIGy09NHkR3rm1hHRykrOHpd7/dvJYC7BLhYBe788I+3siIQABOpFAHOpXnqUNpsmLSz5UC1N9lI7QpdScSbtzLf+0TqOPHCFaxwCmEspuVLHKWljLpVB23cNUMaYKfugJlPSDhsLrcL4WaJzr38LE5cYctWFUvs2/5+9N4G2q6rStmeCgTR0hggWAiVFfRK0aAsRSyma0CuYIOEDTUNjMfxJEBAChKYsKKQzgEBoCulToEjA0OMggJFWKH+6wjvgV8CysMNECT0i+cfeZxDJzcnJWnOtNffa6z5nDAZy73rXXPt559xnnjs958AOdnoCKCEAgVIIMFwqxcl+19GmxpKGJM8kxJc8fXE5lW/947ULVf81cPVn5qKAK8MllzyJtYZ8i0XSfR+Yu7PqttK3BwiLZq8mP+yZayPilZacXld6/evJ9FaSq3qysIOdngBKCECgFAIMl0pxkuFSoU42d1k0is2xD43s+8ISr0OJd9fDFa5pCDBcsuRKHVvS7sSCeRhz3x4gLJq9mvywZ66NiFdacnpd6fWvJ8NwCXapCOj35R6pZ4cSAhDIiwDDpbz8iHaaNjWWPKlGsz3qRvgSFafpZr71j9dp7IErXNMQYLhkyZU6tqTNcCkGbd8eIEZMyz2oSUvaYbHwKoyfRl16/WuYuGjIVRdK3dfADnZ6AighAIFSCDBcWo6T3/nOd+See+6RF154QV555RX50Ic+JJtssokccsghsuGGGy6lvv7662XWrFny/PPPy/Dhw2WbbbaRadOmyVprrbXU2r6+PvnWt74ljz32mKywwgqy9dZbyzHHHCPrrrtucH61qbGkIQm2O8kG+JIEq8mmvvWP12lsgStc0xBguGTJlTq2pM1wKQZt3x4gRkzLPahJS9phsfAqjJ9GXXr9a5i4aMhVF0oMl/SUYBebHftBAAJ5EWC4tBw/vva1r8kqq6wiG2ywgay66qry29/+Vm688UZ56aWX6iHSZptttniHb3/723LRRRfJZz/7Wdlhhx1k/vz59ZpKd8MNN8jqq6++eO0vfvELGT9+vIwaNUomTJggb731llx11VX17+fMmVP/POTRpsaSZi7E6XRafEnHNvXOvvWP12kcgStc0xDovatv/TdxxjbFpI7t3YJ5GPPS7wHkR1h+WKrxypJ2J1bp9Z+KKLmqJws72OkJoIQABEohwHBJ4WQ1WNpuu+1kl112kbPPPrveofrZ9ttvL5/+9KelerfTe4+nnnpK9tlnHznooIPkqKOOWvzzKVOmyEMPPSR33HHH4nc1PfvsszJ27Fj50pe+JCeccILiZH+VtKmxpCEJsjqZGF+SoU2+sW/943UaS+AK1zQEeu/qW/9NnLFNMalje7dgHsa89HsA+RGWH5ZqvLKk3YlVev2nIkqu6snCDnZ6AighAIFSCDBcUjj57rvvypZbbilbbLGFXHrppfUOd911l0ydOrX+mLs999xziV0/97nPyauvvirz5s2rf/7aa6/Jpz71Kdljjz3ktNNOW2LtAQccIM8884w8+OCDipP9VdKmxpKGJMjqZGJ8SYY2+ca+9Y/XaSyBK1zTEOi9q2/9N3HGNsWkju3dgnkY89LvAeRHWH5YqvHKknYnVun1n4oouaonCzvY6QmghAAESiHAcMnRyQULFsiiRYvk97//vVx55ZX1R9dNnz5d9t9//3qH2267Tb7+9a/LzJkzZaeddlpi17333luqdzA98MAD9cfdVU3ffvvtJyeddJLsu+++S6w955xz5OKLL64HUR/+8IcdT7f0sjY1ljQkapuTCvElKd6km/vWP16nsQOucE1DoPeuvvXfxBnbFJM6tncL5mHMS78HkB9h+WGpxitL2p1Ypdd/KqLkqp4s7GCnJ4ASAhAohQDDJUcnN9xww8Urhw8fLhMnTpTDDjtMVlhhhfrnfX199UfaTZo0SY4//vjFa6uhVPX9S2+88Ub9XU2f+MQn5M4776y11fczVb97/+Oaa66Rk08+Wb7//e/Lpptu6ni6pZdVjWU1DBsxYoR6DythxaZ6DBs2zCokcRwI4IsDpERLRo8eHbSzb/3jdRDuZYrhClcNAev615xxIGmoY3u3BzLz0Pqv3PLtAewdDos4kPMjjJy9Gq/8mYfeA0qvf3+ibgpy1Y1Tt1Wwi8cutP71J0EJAQhAIIzAgBkuVYOWt99+24nW4MGDZciQIUusrT6m7p133pFf/epX9buWqsHPtGnTZKWVVlq8rvqupCeffLL+bqUxY8bI/Pnz5cwzz6x/9uc//1mqwVH1cXqV/phjjpHLLrtMPvvZzy4RZ/bs2fVw6uqrr64/Ok/7aFNjSUOidTmtDl/S8u21e2hj6Vv/eJ3Ga7jCVUPAuv41ZxxIGurY3u2BzDy0/iu3fHsAe4fDIg7k/AgjZ6/GK3/mofeA0uvfn6ibglx149RtFezisQutf/1JUEIAAhAIIzBghkvPPvts/R1HLo9x48bJ6aefvsylr7zySr1XNSiaMWPG4nV/+MMf6oHT+78vadttt5W/+Zu/ke9973ty0003SfWEYfXOpepg1fdC5f7grdR5OoQvefricirfj8TAaxeq/mvg6s/MRQHX3pR869+F+UBeQ77Zuw/zMOal3wPIj7D8sFTjlSXtTqzS6z8VUXJVTxZ2sNMTQAkBCJRCYMAMlxYuXChz58518m299darB0e9Ht/4xjfqj657/PHHl3j3UqV58cUX5de//rWstdZaUu11xBFH1AOl//qv/6o/ps7lO5d+9KMf1UMp7aNNjSUNidbltDp8Scs35e6+9Y/XadyAK1zTEGC4ZMmVOrak3YkF8zDmvj1AWDR7Nflhz1wbEa+05PS60utfT6a3klzVk4Ud7PQEUEIAAqUQGDDDpdiGTZ8+vf4OpepdSmusscYyt68+Dq9699JHP/pRufbaa+t1r776qmy99db1u59OO+20JbQHHHCAPPPMM/LAAw/IoEGD1MduU2NJQ6K2OakQX5LiTbq5b/3jdRo74ArXNAR67+pb/02csU0xqWN7t2Aexrz0ewD5EZYflmq8sqTdiVV6/aciSq7qycIOdnoCKCEAgVIIMFzq4eTrr79e/3b48OFLrPrd734ne+65Z/0upHvuuadnLsycOVPOP/98ueiii2SHHXZYvPaQQw6Rhx9+uH5H05prrln/vProvrFjx8p+++0nJ554YlCOtamxpCEJsjqZGF+SoU2+sW/943UaS+AK1zQEeu/qW/9NnLFNMalje7dgHsa89HsA+RGWH5ZqvLKk3YlVev2nIkqu6snCDnZ6AighAIFSCDBc6uFkX1+fTJ48WXbddVdZf/316yHTCy+8UL9jqfrepWpoNGbMmMU7nHHGGfKnP/1JPv7xj8sHPvABmTdvntx7770yceJEOeGEE5aI9POf/1zGjx8vH/rQh2TChAny9ttvy5VXXlmvqfZ/b+CkTbQ2NZY0JFqX0+rwJS3flLv71j9ep3EDrnBNQ4DhkiVX6tiSdicWzMOY+/YAYdHs1eSHPXNtRLzSktPrSq9/PZneSnJVTxZ2sNMTQAkBCJRCgOFSDycXLFgg5513Xv1dSb/5zW/kzTffrD8C7x//8R+l+vi6TTbZZAn1LbfcIldccUU9gHr33XflYx/7mHzpS1+q343U7fH000/LjBkz6u9tGjx4cP1ReUcffbT87d/+bXB+tamxpCEJtjvJBviSBKvJpr71j9dpbIErXNMQ6L2rb/03ccY2xaSO7d2CeRjz0u8B5EdYfliq8cqSdidW6fWfiii5qicLO9jpCaCEAARKIcBwqRQn+11HmxpLGpI8kxBf8vTF5VS+9Y/XLlT918DVn5mLAq4Ml1zyJNYa8i0WSfd9YO7OqttK3x4gLJq9mvywZ66NiFdacnpd6fWvJ9NbSa7qycIOdnoCKCEAgVIIMFwqxUmGS4U62dxl0Sg2xz40su8LS7wOJd5dD1e4piHAcMmSK3VsSbsTC+ZhzH17gLBo9mryw565NiJeacnpdaXXv54MwyXYpSKg35d7pJ4dSghAIC8CDJfy8iPaadrUWPKkGs32qBvhS1Scppv51j9ep7EHrnBNQ4DhkiVX6tiSNsOlGLR9e4AYMS33oCYtaYfFwqswfhp16fWvYeKiIVddKHVfAzvY6QmghAAESiHAcKkUJ/tdR5saSxqSPJMQX/L0xeVUvvWP1y5U/dfA1Z+ZiwKuDJdc8iTWGvItFkn3fWDuzqrbSt8eICyavZr8sGeujYhXWnJ6Xen1ryfTW0mu6snCDnZ6AighAIFSCDBcKsVJhkuFOtncZdEoNsc+NLLvC0u8DiXeXQ9XuKYhwHDJkit1bEm7EwvmYcx9e4CwaPZq8sOeuTYiXmnJ6XWl17+eDMMl2KUioN+Xe6SeHUoIQCAvAgyX8vIj2mna1FjypBrN9qgb4UtUnKab+dY/XqexB65wTUOA4ZIlV+rYkjbDpRi0fXuAGDEt96AmLWmHxcKrMH4aden1r2HioiFXXSh1XwM72OkJoIQABEohwHCpFCf7XUebGksakjyTEF/y9MXlVL71j9cuVP3XwNWfmYsCrgyXXPIk1hryLRZJ931g7s6q20rfHiAsmr2a/LBnro2IV1pyel3p9a8n01tJrurJwg52egIoIQCBUggwXCrFSYZLhTrZ3GXRKDbHPjSy7wtLvA4l3l0PV7imIcBwyZIrdWxJuxML5mHMfXuAsGj2avLDnrk2Il5pyel1pde/ngzDJdilIqDfl3uknh1KCEAgLwIMl/LyI9pp2tRY8qQazfaoG+FLVJymm/nWP16nsQeucE1DgOGSJVfq2JI2w6UYtH17gBgxLfegJi1ph8XCqzB+GnXp9a9h4qIhV10odV8DO9jpCaCEAARKIcBwqRQn+11HmxpLGpI8kxBf8vTF5VS+9Y/XLlT918DVn5mLAq4Ml1zyJNYa8i0WSfd9YO7OqttK3x4gLJq9mvywZ66NiFdacnpd6fWvJ9NbSa7qycIOdnoCKCEAgVIIMFwqxUmGS4U62dxl0Sg2xz40su8LS7wOJd5dD1e4piHAcMmSK3VsSbsTC+ZhzH17gLBo9mryw565NiJeacnpdaXXv54MwyXYpSKg35d7pJ4dSghAIC8CDJfy8iPaadrUWPKkGs32qBvhS1Scppv51j9ep7EHrnBNQ4DhkiVX6tiSNsOlGLR9e4AYMS33oCYtaYfFwqswfhp16fWvYeKiIVddKHVfAzvY6QmghAAESiHAcKkUJ/tdR5saSxqSPJMQX/L0xeVUvvWP1y5U/dfA1Z+ZiwKuDJdc8iTWGvItFkn3fWDuzqrbSt8eICyavZr8sGeujYhXWnJ6Xen1ryfTW0mu6snCDnZ6AighAIFSCDBcKsVJhkuFOtncZdEoNsc+NLLvC0u8DiXeXQ9XuKYhwHDJkit1bEm7EwvmYcx9e4CwaPZq8sOeuTYiXmnJ6XWl17+eDMMl2KUioN+Xe6SeHUoIQCAvAgyX8vIj2mna1FjypBrN9qgb4UtUnKab+dY/XqexB65wTUOA4ZIlV+rYkjbDpRi0fXuAGDEt96AmLWmHxcKrMH4aden1r2HioiFXXSh1XwM72OkJoIQABEohwHCpFCf7XUebGksakjyTEF/y9MXlVL71j9cuVP3XwNWfmYsCrgyXXPIk1hryLRZJ931g7s6q20rfHiAsmr2a/LBnro2IV1pyel3p9a8n01tJrurJwg52egIoIQCBUggwXCrFSYZLhTrZ3GXRKDbHPjSy7wtLvA4l3l0PV7imIcBwyZIrdWxJuxML5mHMfXuAsGj2avLDnrk2Il5pyel1pde/ngzDJdilIqDfl3uknh1KCEAgLwIMl/LyI9pp2tRY8qQazfaoG+FLVJymm/nWP16nsQeucE1DgOGSJVfq2JI2w6UYtH17gBgxLfegJi1ph8XCqzB+GnXp9a9h4qIhV10odV8DO9jpCaCEAARKIcBwqRQn+11HmxpLGpI8kxBf8vTF5VS+9Y/XLlT918DVn5mLAq4Ml1zyJNYa8i0WSfd9YO7OqttK3x4gLJq9mvywZ66NiFdacnpd6fWvJ9NbSa7qycIOdnoCKCEAgVIIMFwqxUmGS4U62dxl0Sg2xz40su8LS7wOJd5dD1e4piHAcMmSK3VsSbsTC+ZhzH17gLBo9mryw565NiJeacnpdaXXv54MwyXYpSKg35d7pJ4dSghAIC8CDJfy8iPaadrUWPKkGs32qBvhS1Scppv51j9ep7EHrnBNQ4DhkiVX6tiSNsOlGLR9e4AYMS33oCYtaYfFwqswfhp16fWvYeKiIVddKHVfAzvY6QmghAAESiHAcKkUJ/tdR5saSxqSPJMQX/L0xeVUvvWP1y5U/dfA1Z+ZiwKuDJdc8iTWGvItFkn3fWDuzqrbSt8eICyavZr8sGeujYhXWnJ6Xen1ryfTW0mu6snCDnZ6AighAIFSCDBcKsVJhkuFOtncZdEoNsc+NLLvok/emAAAIABJREFUC0u8DiXeXQ9XuKYhwHDJkit1bEm7EwvmYcx9e4CwaPZq8sOeuTYiXmnJ6XWl17+eDMMl2KUioN+Xe6SeHUoIQCAvAgyX8vIj2mna1FjypBrN9qgb4UtUnKab+dY/XqexB65wTUOA4ZIlV+rYkjbDpRi0fXuAGDEt96AmLWmHxcKrMH4aden1r2HioiFXXSh1XwM72OkJoIQABEohwHCpFCf7XUebGksakjyTEF/y9MXlVL71j9cuVP3XwNWfmYsCrgyXXPIk1hryLRZJ931g7s6q20rfHiAsmr2a/LBnro2IV1pyel3p9a8n01tJrurJwg52egIoIQCBUggwXCrFSYZLhTrZ3GXRKDbHPjSy7wtLvA4l3l0PV7imIcBwyZIrdWxJuxML5mHMfXuAsGj2avLDnrk2Il5pyel1pde/ngzDJdilIqDfl3uknh1KCEAgLwIMl/LyI9pp2tRY8qQazfaoG+FLVJymm/nWP16nsQeucE1DgOGSJVfq2JI2w6UYtH17gBgxLfegJi1ph8XCqzB+GnXp9a9h4qIhV10odV8DO9jpCaCEAARKIcBwqRQn+11HmxpLGpI8kxBf8vTF5VS+9Y/XLlT918DVn5mLAq4Ml1zyJNYa8i0WSfd9YO7OqttK3x4gLJq9mvywZ66NiFdacnpd6fWvJ9NbSa7qycIOdnoCKCEAgVIIMFwqxUmGS4U62dxl0Sg2xz40su8LS7wOJd5dD1e4piHAcMmSK3VsSbsTC+ZhzH17gLBo9mryw565NiJeacnpdaXXv54MwyXYpSKg35d7pJ4dSghAIC8CDJfy8iPaadrUWPKkGs32qBvhS1Scppv51j9ep7EHrnBNQ4DhkiVX6tiSNsOlGLR9e4AYMS33oCYtaYfFwqswfhp16fWvYeKiIVddKHVfAzvY6QmghAAESiHAcKkUJ/tdR5saSxqSPJMQX/L0xeVUvvWP1y5U/dfA1Z+ZiwKuDJdc8iTWGvItFkn3fWDuzqrbSt8eICyavZr8sGeujYhXWnJ6Xen1ryfTW0mu6snCDnZ6AighAIFSCDBcKsVJhkuFOtncZdEoNsc+NLLvC0u8DiXeXQ9XuKYhwHDJkit1bEm7EwvmYcx9e4CwaPZq8sOeuTYiXmnJ6XWl17+eDMMl2KUioN+Xe6SeHUoIQCAvAgyX8vIj2mna1FjypBrN9qgb4UtUnKab+dY/XqexB65wTUOA4ZIlV+rYkjbDpRi0fXuAGDEt96AmLWmHxcKrMH4aden1r2HioiFXXSh1XwM72OkJoIQABEohwHCpFCf7XUebGksakjyTEF/y9MXlVL71j9cuVP3XwNWfmYsCrgyXXPIk1hryLRZJ931g7s6q20rfHiAsmr2a/LBnro2IV1pyel3p9a8n01tJrurJwg52egIoIQCBUggwXCrFSYZLhTrZ3GXRKDbHPjSy7wtLvA4l3l0PV7imIcBwyZIrdWxJuxML5mHMfXuAsGj2avLDnrk2Il5pyel1pde/ngzDJdilIqDfl3uknh1KCEAgLwIMl/LyI9pp2tRY8qQazfaoG+FLVJymm/nWP16nsQeucE1DgOGSJVfq2JI2w6UYtH17gBgxLfegJi1ph8XCqzB+GnXp9a9h4qIhV10odV8DO9jpCaCEAARKIcBwqRQn+11HmxpLGpI8kxBf8vTF5VS+9Y/XLlT918DVn5mLAq4Ml1zyJNYa8i0WSfd9YO7OqttK3x4gLJq9mvywZ66NiFdacnpd6fWvJ9NbSa7qycIOdnoCKCEAgVIIMFwqxUmGS4U62dxl0Sg2xz40su8LS7wOJd5dD1e4piHAcMmSK3VsSbsTC+ZhzH17gLBo9mryw565NiJeacnpdaXXv54MwyXYpSKg35d7pJ4dSghAIC8CDJfy8iPaadrUWPKkGs32qBvhS1Scppv51j9ep7EHrnBNQ4DhkiVX6tiSNsOlGLR9e4AYMS33oCYtaYfFwqswfhp16fWvYeKiIVddKHVfAzvY6QmghAAESiHAcKkUJ/tdR5saSxqSPJMQX/L0xeVUvvWP1y5U/dfA1Z+ZiwKuDJdc8iTWGvItFkn3fWDuzqrbSt8eICyavZr8sGeujYhXWnJ6Xen1ryfTW0mu6snCDnZ6AighAIFSCDBcKsVJhkuFOtncZdEoNsc+NLLvC0u8DiXeXQ9XuKYhwHDJkit1bEm7EwvmYcx9e4CwaPZq8sOeuTYiXmnJ6XWl17+eDMMl2KUioN+Xe6SeHUoIQCAvAgyX8vIj2mna1FjypBrN9qgb4UtUnKab+dY/XqexB65wTUOA4ZIlV+rYkjbDpRi0fXuAGDEt96AmLWmHxcKrMH4aden1r2HioiFXXSh1XwM72OkJoIQABEohwHCpFCf7XUebGksakjyTEF/y9MXlVL71j9cuVP3XwNWfmYsCrgyXXPIk1hryLRZJ931g7s6q20rfHiAsmr2a/LBnro2IV1pyel3p9a8n01tJrurJwg52egIoIQCBUggwXCrFSYZLhTrZ3GXRKDbHPjSy7wtLvA4l3l0PV7imIcBwyZIrdWxJuxML5mHMfXuAsGj2avLDnrk2Il5pyel1pde/ngzDJdilIqDfl3uknh1KCEAgLwIMl/LyI9pp2tRY8qQazfaoG+FLVJymm/nWP16nsQeucE1DgOGSJVfq2JI2w6UYtH17gBgxLfegJi1ph8XCqzB+GnXp9a9h4qIhV10odV8DO9jpCaCEAARKIcBwqRQn+11HmxpLGpI8kxBf8vTF5VS+9Y/XLlT918DVn5mLAq4Ml1zyJNYa8i0WSfd9YO7OqttK3x4gLJq9mvywZ66NiFdacnpd6fWvJ9NbSa7qycIOdnoCKCEAgVIIZDFceuutt+TVV1+V1VdfXVZYYYXFbH/84x/LvffeKyuuuKLss88+ssEGG5TCPfl1tKmxpCFJng6qAPiiwpaFyLf+8TqNbXCFaxoCDJcsuVLHlrQ7sWAexty3BwiLZq8mP+yZayPilZacXld6/evJMFyCXSoC+n25R+rZoYQABPIikMVw6d/+7d/kpptukvvvv19GjBhRE7r++uvlX//1X2XRokX1f1c/nz17tqy//vp5Ecz0NG1qLHlSzTOJ8CVPX1xO5Vv/eO1C1X8NXP2ZuSjg2puSb/27MB/Ia8g3e/dhHsa89HsA+RGWH5ZqvLKk3YlVev2nIkqu6snCDnZ6AighAIFSCGQxXNptt93qodGFF164mOv2228vgwYNkm9961syf/58Ofroo6Vad9ppp5XCPul1tKmxpCFJmgrqzfFFja5xoW/943Uay+AK1zQEGC5ZcqWOLWl3YsE8jLlvDxAWzV5Nftgz10bEKy05va70+teT6a0kV/VkYQc7PQGUEIBAKQSyGC5tueWWMn78eDnmmGNqrs8++6zsueeeMn36dJk8eXL9s6OOOkoef/xxmTt3binsk15HmxpLGpKkqaDeHF/U6BoX+tY/XqexDK5wTUOg966+9d/EGdsUkzq2dwvmYcxLvweQH2H5YanGK0vanVil138qouSqnizsYKcngBICECiFQBbDpc0331z222+/+t1J1WPWrFly6qmnyq233rr4e5bOPvtsueqqq+SJJ54ohX3S62hTY0lDkjQV1Jvjixpd40Lf+sfrNJbBFa5pCDBcsuRKHVvS7sSCeRhz3x4gLJq9mvywZ66NiFdacnpd6fWvJ9NbSa7qycIOdnoCKCEAgVIIZDFc+vznPy8f/OAH66FS9Zg0aZL88pe/lHnz5i3mfNxxx8mPf/zj+nuZeCyfQJsaSxqS5fvZxAp8aYJ6nJi+9Y/Xcbj33wWucE1DoPeuvvXfxBnbFJM6tncL5mHMS78HkB9h+WGpxitL2p1Ypdd/KqLkqp4s7GCnJ4ASAhAohUAWw6WLLrpIzj33XNl5551l6NChcsstt8iBBx4o06ZNW8x57733lpVWWkmuueaaUtgnvY42NZY0JElTQb05vqjRNS70rX+8TmMZXOGahgDDJUuu1LEl7U4smIcx9+0BwqLZq8kPe+baiHilJafXlV7/ejK9leSqnizsYKcngBICECiFQBbDpTfffLP+TqW7775bFi1aJJ/5zGfk/PPPl+HDh9ecf/GLX8jnPvc5mTp1av1Pk48JEybIo48+KnvssYfMmDFjqaM89NBD9aCsr69Phg0bJttvv309JBs5cuRSa2+//Xa55JJL6utbbbXVZLfddpPDDz9cRowYEXyJbWosaUiC7U6yAb4kwWqyqW/943UaW+AK1zQEeu/qW/9NnLFNMalje7dgHsa89HsA+RGWH5ZqvLKk3YlVev2nIkqu6snCDnZ6AighAIFSCGQxXHoP5iuvvCKDBg2SlVdeeQm+CxYskN///vfykY98RFZZZZXG2M+ZM0dOOukkef3117sOlx555BE54IADZMMNN5TqnVbVuS+//HJZe+21Zfbs2fW7st573HzzzfXQaeutt64HZy+88IJcffXVsuWWW8oVV1xRcwh5tKmxpCEJcTqdFl/SsU29s2/943UaR+AK1zQEeu/qW/9NnLFNMalje7dgHsa89HsA+RGWH5ZqvLKk3YlVev2nIkqu6snCDnZ6AighAIFSCGQxXJo+fbqMHj1aJk+enC3XhQsX1u8sqs541llndR0ujR07Vl5++WW57bbbFr/rqvreqIMPPliqa9x///3r63v77bfrdzRVQ6fvfe97ssIKK9Q/v/baa+vh1cyZM2WnnXYKYtGmxpKGJMjqZGJ8SYY2+ca+9Y/XaSyBK1zTEGC4ZMmVOrak3YkF8zDmvj1AWDR7Nflhz1wbEa+05PS60utfT6a3klzVk4Ud7PQEUEIAAqUQyGK4tOmmm8qkSZPkyCOPzJbrySefLA888ED9fVAbb7zxUsOl559/XnbddVc59NBDl/rovuq7pKqPvbv++uvr67v//vvloIMOkjPOOEOqgdR7j2ro9KlPfUq23XZb+fa3vx3Eok2NJQ1JkNXJxPiSDG3yjX3rH6/TWAJXuKYh0HtX3/pv4oxtikkd27sF8zDmpd8DyI+w/LBU45Ul7U6s0us/FVFyVU8WdrDTE0AJAQiUQiCL4dJee+0l66+/fv2OoBwfTz/9dP0xdxdffHE9+Kk+9q7/dy5VQ6fqe6MuvfRS2WabbZa4jOrnP/zhD+Xxxx+v36VU7XPOOefIHXfcIX/3d3+3xNr99ttP/vCHP8hdd90VhKJqLKvvr4rx/U1BB3EQv/HGG/Wq6juqeORDAF+a86J6J2fIw7f+8TqE9rK1cIWrhoB1/WvOOJA01LG92wOZeWj9v/fH5ba8BtBk10DODw2vJjV45U8/9B7g+xrA/4RlKshVva+wi8cutP71J0EJAQhAIIxAFsOl22+/vf7YuFmzZskmm2wSdkWR1e+++67su+++MnLkyHooVD26DZcuu+wyOfPMM6X6LqXq9+9/VD+vfl+982nUqFFSvQvqmmuukZ/+9KdLfb/UYYcdJvfee688+eSTQVfSpsaShiTI6mRifEmGdrkbhzaWvvWP18u1RLUAripsyxWVztW6/pcLfIAvKD3fcrR3IDMPrf/KT98eIMcc6HWmgZwfeNU2Av7nDb0HlF7//kTdFNxX3Dh1WwW7eOxC619/EpQQgAAEwghkMVyaM2dO/XFzP/nJT2SXXXaRjTbaSNZYYw0ZNGjQUlf3/o+R87n06v/BV33snMtj8ODBMmTIkHrpddddJ6ecckr9PUrrrbde/bNuw6ULLrhAzjvvPLnzzjvrd2G9/3HuuefKhRdeKHfffbess846ctxxx8kNN9wgTz31lKy44opLrD366KPlpptuqj9vvtv1u5z/vReW1b+32GILV0lj63grdWPoewbGlzx9cTmV70di4LULVf81cPVn5qKAa29KvvXvwnwgryHf7N2HeRjz0u8B5EdYfliq8cqSdidW6fWfiii5qicLO9jpCaCEAARKIZDFcKma0FeDlGoA9P7H+4cr1e+q/+7r61Oxf/bZZ+uPsnN5jBs3Tk4//XRZsGCB7LbbblJ9VN3hhx++WNqWdy5VB2a45OI4a7oRoFFsb174vrDE6zRewxWuaQj03tW3/ps4Y5tiUsf2bsE8jHnp9wDyIyw/LNV4ZUm7E6v0+k9FlFzVk4Ud7PQEUEIAAqUQyGK49IMf/MCZZzX40TwWLlwoc+fOdZJW71Dacsst63csVe+oqj7CbujQoYu1Y8aMkR133LH+KL/VV1+9/mi75X3nUvWOpieeeMLpO5deeukl57Mu64La1FjSkDilpfkifDFHHi2gb/3jdTT0S2wEV7imIdB7V9/6b+KMbYpJHdu7BfMw5qXfA8iPsPywVOOVJe1OrNLrPxVRclVPFnaw0xNACQEIlEIgi+FSrjAPOeSQ+qPsej1OPPFEmTBhgjz33HP1u5wOPfRQmTp16hKSnXfeWVZddVWZPXt2/fP77rtPvvKVr8gZZ5wh7/+Yv+pj+z71qU/JP//zP0v1UXohjzY1ljQkIU6n0+JLOrapd/atf7xO4whc4ZqGQO9dfeu/iTO2KSZ1bO8WzMOYl34PID/C8sNSjVeWtDuxSq//VETJVT1Z2MFOTwAlBCBQCgGGSz2crN5pVL2LqP9jypQpstVWW8nkyZPr719ad9116yVf+MIXpHqHVPX9TMOHD69/Nm/ePDn44IPlmGOOkQMPPLD+WTVE2m677eQjH/lI/Z1O1Xc8VY9rr71WTjrppPq7m6rvngp5tKmxpCEJcTqdFl/SsU29s2/943UaR+AK1zQEeu/qW/9NnLFNMalje7dgHsa89HsA+RGWH5ZqvLKk3YlVev2nIkqu6snCDnZ6AighAIFSCGQ1XLr33nvl1ltvler7kV599dX64+Y+9rGP1d+VVA1jcnl0+86l6mwPP/xwPUCqvkNq/PjxMn/+fLniiitkrbXWkhtuuEGGDRu2+BLmzJlTD5w+/elPy+677y6//OUv5aqrrpLNN99crr766vr7pUIebWosaUhCnE6nxZd0bFPv7Fv/eJ3GEbjCNQ0BhkuWXKljS9qdWDAPY+7bA4RFs1eTH/bMtRHxSktOryu9/vVkeivJVT1Z2MFOTwAlBCBQCoEshkvVO3mOOOIIueeee2TRokXygQ98oP4uoz/96U/yzjvv1IOWHXbYQc455xxZccUVG2e/rOFSdbAHH3yw/ki7vr6+ephUDcWmTZsmo0aNWurc1TucLrnkkvoj9VZbbTXZdddd5fDDD6+HaqGPNjWWNCShbqfR40sarha7+tY/XqdxBa5wTUOg966+9d/EGdsUkzq2dwvmYcxLvweQH2H5YanGK0vanVil138qouSqnizsYKcngBICECiFQBbDpW9961ty2WWXyTbbbFN/Z9HGG29cD5SqQdNTTz0l559/vtx///319xQdeeSRpbBPeh1taixpSJKmgnpzfFGja1zoW/94ncYyuMI1DQGGS5ZcqWNL2p1YMA9j7tsDhEWzV5Mf9sy1EfFKS06vK73+9WR6K8lVPVnYwU5PACUEIFAKgSyGS9VQ6UMf+pDceOONXblWQ6YvfvGL9fcf3XfffaWwT3odbWosaUiSpoJ6c3xRo2tc6Fv/eJ3GMrjCNQ0BhkuWXKljS9oMl2LQ9u0BYsS03IOatKQdFguvwvhp1KXXv4aJi4ZcdaHUfQ3sYKcngBICECiFQBbDpc0220wmTZokX//615fJ9eyzz66/i+jxxx8vhX3S62hTY0lDkjQV1Jvjixpd40Lf+sfrNJbBFa5pCDBcsuRKHVvSZrgUg7ZvDxAjpuUe1KQl7bBYeBXGT6Muvf41TFw05KoLJYZLekqwi82O/SAAgbwIZDFc2meffWS99daTGTNmLJNO9XF4//u//yvXXXddXgQzPU2bGkuauTyTCF/y9MXlVL71j9cuVP3XwNWfmYsCrgyXXPIk1hryLRZJ931g7s6q20rfHiAsmr2a/LBnro2IV1pyel3p9a8n01tJrurJwg52egIoIQCBUghkMVx68MEH5atf/aqcfvrpsvvuuy/F9rbbbpPjjjtOLr74Yvn0pz9dCvuk19GmxpKGJGkqqDfHFzW6xoW+9Y/XaSyDK1zTEGC4ZMmVOrak3YkF8zDmvj1AWDR7Nflhz1wbEa+05PS60utfT4bhEuxSEdDvyz1Szw4lBCCQF4EshkszZ86Uxx57TKoh09///d9L9TF5I0eOlAULFtQfg/fzn/9c/umf/kk233zzJegNGjRIpkyZkhfRTE7TpsaSJ9VMkqbfMfAlT19cTuVb/3jtQtV/DVz9mbko4Nqbkm/9uzAfyGvIN3v3YR7GvPR7APkRlh+WaryypN2JVXr9pyJKrurJwg52egIoIQCBUghkMVwaPXq0imc1XOrr61NpSxe1qbGkIckzG/ElT19cTuVb/3jtQtV/DVz9mbko4MpwySVPYq0h32KRdN8H5u6suq307QHCotmryQ975tqIeKUlp9eVXv96Mr2V5KqeLOxgpyeAEgIQKIVAFsOlRx55RM1zq622UmtLFrapsaQhyTMT8SVPX1xO5Vv/eO1C1X8NXP2ZuSjgynDJJU9irSHfYpF03wfm7qwYLoWxQp2WALWclu9ArP9URMlVPVnYwU5PACUEIFAKgSyGS1qYr776qixcuFDWXntt7RbF6nz/uNwkCBqSJukvOza+5OmLy6l86x+vXaj6r4GrPzMXBVwZLrnkSaw15Fssku77wNyd1UD84zL5EZYflmq8sqTdieX7GsD+hHlGJFf1vsAOdnoCKCEAgVIItHq4VH1X0wUXXMBH43XJxjY1ljQked5O8CVPX1xO5Vv/eO1C1X8NXP2ZuSjgynDJJU9irSHfYpF03wfm7qwYLoWxQp2WALWclu9ArP9URMlVPVnYwU5PACUEIFAKAYZLpTjZ7zp8/7jcJAYakibpLzs2vuTpi8upfOsfr12o+q+Bqz8zFwVcGS655EmsNeRbLJLu+8DcndVA/OMy+RGWH5ZqvLKk3Ynl+xrA/oR5RiRX9b7ADnZ6AighAIFSCDBcKsXJftfRpsaShiTPJMSXPH1xOZVv/eO1C1X/NXD1Z+aigCvDJZc8ibWGfItF0n0fmLuzYrgUxgp1WgLUclq+A7H+UxElV/VkYQc7PQGUEIBAKQQYLpXiJMOlQp1s7rJoFJtjHxqZ4VIowTh6aigOx/67wJXhUprM6r4r+WZJuxML5mHMfXuAsGj2avLDnrk2Il5pyel1pde/nkxvJbmqJws72OkJoIQABEohwHCpFCf7XUebGksakjyTEF/y9MXlVL71j9cuVP3XwNWfmYsCrgyXXPIk1hryLRZJ931g7s6q20rfHiAsmr2a/LBnro2IV1pyel3p9a8nw3AJdqkI6PflHqlnhxICEMiLAMOlvPyIdpo2NZY8qUazPepG+BIVp+lmvvWP12nsgStc0xBguGTJlTq2pN2JBfMw5r49QFg0ezX5Yc9cGxGvtOT0utLrX0+G4RLsUhHQ78s9Us8OJQQgkBcBhkt5+RHtNG1qLHlSjWZ71I3wJSpO08186x+v09gDV7imIcBwyZIrdWxJm+FSDNq+PUCMmJZ7UJOWtMNi4VUYP4269PrXMHHRkKsulLqvgR3s9ARQQgACpRBguFSKk/2uo02NJQ1JnkmIL3n64nIq3/rHaxeq/mvg6s/MRQFXhksueRJrDfkWi6T7PjB3Z9VtpW8PEBbNXk1+2DPXRsQrLTm9rvT615PprSRX9WRhBzs9AZQQgEApBBguleIkw6VCnWzusmgUm2MfGtn3hSVehxLvrocrXNMQYLhkyZU6tqTdiQXzMOa+PUBYNHs1+WHPXBsRr7Tk9LrS619PhuES7FIR0O/LPVLPDiUEIJAXgVYPl+bOnSt33323nHbaaXlRzeA0bWoseVLNIGG6HAFf8vTF5VS+9Y/XLlT918DVn5mLAq69KfnWvwvzgbyGfLN3H+ZhzEu/B5AfYflhqcYrS9qdWKXXfyqi5KqeLOxgpyeAEgIQKIVAVsOlt956Sx566CF5/vnn5fXXX5cpU6bUnKufv/rqq/LBD35QBg8eXAr7pNfRpsaShiRpKqg3xxc1usaFvvWP12ksgytc0xBguGTJlTq2pN2JBfMw5r49QFg0ezX5Yc9cGxGvtOT0utLrX0+mt5Jc1ZOFHez0BFBCAAKlEMhmuHT77bfLySefLC+//LIsWrRIBg0aJH19fTXnp556SvbZZ5/6HUpjx44thX3S62hTY0lDkjQV1Jvjixpd40Lf+sfrNJbBFa5pCPTe1bf+mzhjm2JSx/ZuwTyMeen3APIjLD8s1XhlSbsTq/T6T0WUXNWThR3s9ARQQgACpRDIYrh0//33y7/8y7/IuuuuK5MnT5bHHntMbrvttsXDpQr25z//+fr3F110USnsk15HmxpLGpKkqaDeHF/U6BoX+tY/XqexDK5wTUOA4ZIlV+rYknYnFszDmPv2AGHR7NXkhz1zbUS80pLT60qvfz2Z3kpyVU8WdrDTE0AJAQiUQiCL4dKXv/xl+dWvflUPlFZZZRWZOXOmXHDBBUsMl44++mj56U9/Wn/HEo/lE2hTY0lDsnw/m1iBL01QjxPTt/7xOg73/rvAFa5pCPTe1bf+mzhjm2JSx/ZuwTyMeen3APIjLD8s1XhlSbsTq/T6T0WUXNWThR3s9ARQQgACpRDIYri0+eab1x93941vfKPm2m24dNZZZ8nVV18tTzzxRCnsk15HmxpLGpKkqaDeHF/U6BoX+tY/XqexDK5wTUOA4ZIlV+rYknYnFszDmPv2AGHR7NXkhz1zbUS80pLT60qvfz2Z3kpyVU8WdrDTE0AJAQiUQiCL4dIWW2whe+21l5xwwgk1127DpenTp8s999wjP/nJT0phn/Q62tRY0pAkTQX15viiRte40Lf+8TqNZXCFaxrjL5fKAAAgAElEQVQCvXf1rf8mztimmNSxvVswD2Ne+j2A/AjLD0s1XlnS7sQqvf5TESVX9WRhBzs9AZQQgEApBLIYLu27777yyiuvyC233CKDBw9earj05ptvyi677CIbbLCBXH755aWwT3odbWosaUiSpoJ6c3xRo2tc6Fv/eJ3GMrjCNQ0BhkuWXKljS9qdWDAPY+7bA4RFs1eTH/bMtRHxSktOryu9/vVkeivJVT1Z2MFOTwAlBCBQCoEshktz5syRY489VsaNG1e/e+mKK65Y/J1Lf/zjH+Vf//VfZe7cuXLeeefJTjvtVAr7pNfRpsaShiRpKqg3xxc1usaFvvWP12ksgytc0xBguGTJlTq2pM1wKQZt3x4gRkzLPahJS9phsfAqjJ9GXXr9a5i4aMhVF0rd18AOdnoCKCEAgVIIZDFcqmCefPLJcu2118qQIUNklVVWkWqotN5668mLL74o77zzjkycOFGOP/74Urgnv442NZY0JMnTQRUAX1TYshD51j9ep7ENrnBNQ4DhkiVX6tiSNsOlGLR9e4AYMS33oCYtaYfFwqswfhp16fWvYeKiIVddKDFc0lOCXWx27AcBCORFIJvhUoXlRz/6kXz3u9+Vp556ShYuXCgjRoyQf/iHf5D99ttPdtxxx7zIZX6aNjWWNHN5JhO+5OmLy6l86x+vXaj6r4GrPzMXBVwZLrnkSaw15Fssku77wNydVbeVvj1AWDR7Nflhz1wbEa+05PS60utfT6a3klzVk4Ud7PQEUEIAAqUQyGK49Oijj8rKK68sG220USlcG7+ONjWWNCSNp0vXA+BLnr64nMq3/vHahar/Grj6M3NRwJXhkkuexFpDvsUi6b4PzN1ZMVwKY4U6LQFqOS3fgVj/qYiSq3qysIOdngBKCECgFAJZDJc+/vGP1+9OOvHEE0vh2vh1+P5xuckD05A0SX/ZsfElT19cTuVb/3jtQtV/DVz9mbko4MpwySVPYq0h32KRdN8H5u6sBuIfl8mPsPywVOOVJe1OLN/XAPYnzDMiuar3BXaw0xNACQEIlEIgi+HStttuW3/sHcOleGnVpsaShiSe7zF3wpeYNG338q1/vE7jD1zhmoYAwyVLrtSxJe1OLJiHMfftAcKi2avJD3vm2oh4pSWn15Ve/3oyvZXkqp4s7GCnJ4ASAhAohUAWw6XTTjtN7r33Xrn55ptl6NChpbBt9Dra1FjSkDSaKssMji95+uJyKt/6x2sXqv5r4OrPzEUBV4ZLLnkSaw35Fouk+z4wd2fVbaVvDxAWzV5Nftgz10bEKy05va70+teTYbgEu1QE9Ptyj9SzQwkBCORFIIvh0htvvCFTp06V1157rf539TF5I0eOzItUy07TpsaSJ9U8kwtf8vTF5VS+9Y/XLlT918DVn5mLAq4Ml1zyJNYa8i0WSfd9YO7OiuFSGCvUaQlQy2n5DsT6T0WUXNWThR3s9ARQQgACpRDIYri00UYb1TwXLVokgwYNWibb6nc/+9nPSmGf9Dp8/7ic9DDL2ZyGpEn6y46NL3n64nIq3/rHaxeq/mvg6s/MRQHX3pR869+F+UBeQ77Zuw/zMOal3wPIj7D8sFTjlSXtTqzS6z8VUXJVTxZ2sNMTQAkBCJRCIIvh0sSJE515zpo1y3ntQF7YpsaShiTPTMWXPH1xOZVv/eO1C1X/NXD1Z+aigCvDJZc8ibWGfItF0n0fmLuz6rbStwcIi2avJj/smWsj4pWWnF5Xev3ryfRWkqt6srCDnZ4ASghAoBQCWQyXSoGZ03W0qbGkIckpc/56FnzJ0xeXU/nWP167UPVfA1d/Zi4KuDJccsmTWGvIt1gk3feBuTsrhkthrFCnJUAtp+U7EOs/FVFyVU8WdrDTE0AJAQiUQoDhUilO9rsO3z8uN4mBhqRJ+suOjS95+uJyKt/6x2sXqv5r4OrPzEUBV4ZLLnkSaw35Fouk+z4wd2c1EP+4TH6E5YelGq8saXdi+b4GsD9hnhHJVb0vsIOdngBKCECgFAIMl0pxst91tKmxpCHJMwnxJU9fXE7lW/947ULVfw1c/Zm5KODKcMklT2KtId9ikXTfB+burBguhbFCnZYAtZyW70Cs/1REyVU9WdjBTk8AJQQgUAqBLIZLkyZNcuI5aNAgueqqq5zWDvRFvn9cbpIXDUmT9JcdG1/y9MXlVL71j9cuVP3XwNWfmYsCrgyXXPIk1hryLRZJ931g7s5qIP5xmfwIyw9LNV5Z0u7E8n0NYH/CPCOSq3pfYAc7PQGUEIBAKQSyGC7tsMMOXXm+9tpr8vLLL0s1VBo1apQMGTJE7rnnnlLYJ72ONjWWNCRJU0G9Ob6o0TUu9K1/vE5jGVzhmoYAwyVLrtSxJe1OLJiHMfftAcKi2avJD3vm2oh4pSWn15Ve/3oyvZXkqp4s7GCnJ4ASAhAohUAWw6VeMF988UU544wz5He/+51cfvnlMmLEiFLYJ72ONjWWNCR+qfDGa2/K22+8LSsOW1GGjRjqJ/ZYjS8esDJb6lv/eJ3GwFRcre4BaaiE75qKa/jJ8tjBt/7zOHW+p8gt3wZC/efGPN/s7H6y0u8B5Ed7MhKv7L0qvf5TEW1LrubYA7SFXarcCdkXdiH00EIAAjkRyH64VMF65513ZNy4cbLlllvKN77xjZz4ZXuWNjWWPKm6pdErC16V377we7nuzJvkN7/4rfzNBh+W/3v0F+TDH11TVhm5stsmHqvwxQNWZkt96x+v0xgYm6v1PSANlfBdY3MNP1FeO/jWf16nz+80ueTbQKr/XJjnl41uJyr9HkB+uOVBDqvwyt6F0us/FdHcczXnHiB3dqlyJsa+sItBkT0gAIEcCLRiuFSB+uY3vym33XabPPjggzlwy/4MbWoseVJdfjpVDeX1Z98i3z31xqUW73fcXjL+63tEHzDhy/J9yXWFb/3jdRonY3Jt4h6Qhkr4rjG5hp8mvx186z+/K8jrRDnk20Cr/xyY55WFfqcp/R5AfvjlQ5Or8cqefun1n4pozrmaew+QM7tU+RJrX9jFIsk+EIBA0wRaM1w69thj5Y477pAnnniiaWatiN+mxpIn1eWn1P/3/z4nh2x5zDIXXvhfZ8j/2eLvlr+Rxwp88YCV2VLf+sfrNAbG5NrEPSANlfBdY3INP01+O/jWf35XkNeJcsi3gVb/OTDPKwv9TlP6PYD88MuHJlfjlT390us/FdGcczX3HiBndqnyJda+sItFkn0gAIGmCWQ/XFq0aJHceuutMn36dNl0003lmmuuaZpZK+K3qbHkSbV3SlWfrXzWQRfJvO8v+1172/3ff5IjL/1/ZGjE72DCl1aUetdD+tY/XqfxOhbXpu4BaaiE7xqLa/hJ8tzBt/7zvIp8TtV0vg3E+m+aeT7ZpztJ6fcA8kOXF02o8Mqeeun1n4porrnahh4gV3apciXmvrCLSZO9IACBJglkMVwaM2ZMVwZ/+ctfZP78+fV3Lo0YMUKuuOIK2XjjjZvk1ZrYbWoseVLtnVYv/2GhHLfbN+XZnz63zIUf23IDOfX242S1UatGy1F8iYbSfCPf+sfrNBbF4trUPSANlfBdY3ENP0meO/jWf55Xkc+pms63gVj/TTPPJ/t0Jyn9HkB+6PKiCRVe2VMvvf5TEc01V9vQA+TKLlWuxNwXdjFpshcEINAkgSyGSxMnTuzKYPDgwbLqqqvKJz7xCdlrr71kzTXXbJJVq2K3qbHkSbV3ajX1/1jCl1aV/BKH9a1/vE7jdSyuTd0D0lAJ3zUW1/CT5LmDb/3neRX5nKrpfBuI9d8083yyT3eS0u8B5IcuL5pQ4ZU99dLrPxXRXHO1DT1AruxS5UrMfWEXkyZ7QQACTRLIYrjUJIBSY7epseRJdflZ2MRnLePL8n3JdYVv/eN1Gidjcm3iHpCGSviuMbmGnya/HXzrP78ryOtEOeTbQKv/HJjnlYV+pyn9HkB++OVDk6vxyp5+6fWfimjOuZp7D5Azu1T5Emtf2MUiyT4QgEDTBBguNe1Aovhtaix5Ul1+Eryy4FW5/uxb5Lun3rjU4i8dv5fsfcQessrIlZe/kccKfPGAldlS3/rH6zQGxuTaxD0gDZXwXWNyDT9Nfjv41n9+V5DXiXLIt4FW/zkwzysL/U5T+j2A/PDLhyZX45U9/dLrPxXRnHM19x4gZ3ap8iXWvrCLRZJ9IACBpglkNVxasGCBzJ07V5555hl59dVXZeWVV5aPfexjstNOO8nIkSObZtWq+G1qLHlSdUutqrH87Qu/l+9/6yb59S9+J2tvsJbsM+0L8uGPrhl9sFSdCF/cfMlxlW/943UaF2Nztb4HpKESvmtsruEnymsH3/rP6/T5nSaXfBtI9Z8L8/yy0e1Epd8DyA+3PMhhFV7Zu1B6/acimnuu5twD5M4uVc7E2Bd2MSiyBwQgkAOBbIZLV155pXz729+Wt956SxYtWrQEm6FDh8oRRxwhkydPzoFZK87QpsaSJ1W/lHrztTflrTfelpWGrShDRwz1E3usxhcPWJkt9a1/vE5jYCquVveANFTCd03FNfxkeezgW/95nDrfU+SWbwOh/nNjnm92dj9Z6fcA8qM9GYlX9l6VXv+piLYlV3PsAdrCLlXuhOwLuxB6aCEAgZwIZDFcuv766+XEE0+UNddcUyZNmiSbbbaZrLHGGjJ//nx57LHHZNasWfLSSy/JKaecIl/84hcb5TdhwgR59NFHZY899pAZM2YscZbbb79dfvSjH8mTTz4pL7zwQn09P/7xj5d53oceekjOPfdc6evrk2HDhsn2228v06ZNi/IurTY1ljypNprSywyOL3n64nIq3/rHaxeq/mvg6s/MRQHX3pR869+F+UBeQ77Zuw/zMOal3wPIj7D8sFTjlSXtTqzS6z8VUXJVTxZ2sNMTQAkBCJRCIIvh0u677y6vvfaa3HTTTbL66qsvxbb6uLyxY8fWH5NXDXCaesyZM0dOOukkef3117sOlyZOnCj//d//LZ/4xCfq4dLgwYOXOVx65JFH5IADDpANN9xQ9t57b6mu8fLLL5e1115bZs+eLdW7tUIebWosaUhCnE6nxZd0bFPv7Fv/eJ3GEbjCNQ0BhkuWXKljS9qdWDAPY+7bA4RFs1eTH/bMtRHxSktOryu9/vVkeivJVT1Z2MFOTwAlBCBQCoEshkubbLKJfOlLX5Jjjz12mVxPO+00+e53v1u/K6iJx8KFC2W33XarP5rvrLPO6jpc+s1vflO/W2mFFVaQatD0y1/+cpnDpWpY9vLLL8ttt90mw4cPry9p3rx5cvDBB8v06dNl//33D7rMNjWWNCRBVicT40sytMk39q1/vE5jCVzhmoZA711967+JM7YpJnVs7xbMw5iXfg8gP8Lyw1KNV5a0O7FKr/9URMlVPVnYwU5PACUEIFAKgSyGS7vssot89rOfrT8ab1mPk08+WR544AH54Q9/2Aj79+LfcsstsvHGG3cdLr3/YL2GS88//7zsuuuucuihh8rUqVOXuJ6dd95ZVlttNak+KjDk0abGkoYkxOl0WnxJxzb1zr71j9dpHIErXNMQ6L2rb/03ccY2xaSO7d2CeRjz0u8B5EdYfliq8cqSdidW6fWfiii5qicLO9jpCaCEAARKIZDFcOnqq6+WCy64oP44uHXXXXcptv/zP/8j48ePr4cx1XceWT+efvrp+qPrLr74Ytl2223rj7Lr9p1L7z9Xr+FSNaA66qij5NJLL5Vtttlmicupfl4N0B5//PH6HVDaR5saSxoSrctpdfiSlm/K3X3rH6/TuAFXuKYh0HtX3/pv4oxtikkd27sF8zDmpd8DyI+w/LBU45Ul7U6s0us/FVFyVU8WdrDTE0AJAQiUQiCL4dKjjz4q3/nOd6T691577SWbbbaZjBw5sv4eoscee0x+8IMfyFZbbSVf+cpXluL+yU9+MqkX7777ruy77771earhUvUIHS5ddtllcuaZZ8rNN99c7/X+R/Xz6vfVu7RGjRqlvraqsVy0aJGMGDFCvYeV8I033qhDDRs2zCokcRwI4IsDpERLRo8eHbSzb/3jdRDuZYrhClcNAev615xxIGmoY3u3BzLz0Pqv3PLtAewdDos4kPMjjJy9Gq/8mYfeA0qvf3+ibgpy1Y1Tt1Wwi8cutP71J0EJAQhAIIxAFsOl6iY6aNCgehhSPar//d7jvZ/1//l7v+/r63MiUO3z9ttvO60dPHiwDBkypF573XXXySmnnFJ/N9J6661X/yx0uFS9S+u8886TO++8U9Zff/0lznTuuefKhRdeKHfffbess846TufttqhNjSUNidrmpEJ8SYq35+ahjaVv/eN1Gq/hClcNAev615xxIGmoY3u3BzLz0Pqv3PLtAewdDos4kPMjjJy9Gq/8mYfeA0qvf3+ibgpy1Y1Tt1Wwi8cutP71J0EJAQhAIIxAFsOl888/f4mBks8l9f/OomVpn3322fqj7Fwe48aNk9NPP71+59Ruu+0m++23nxx++OGLpaHDJat3LlUH3mKLLVwuudE1vJW6UfzLDI4vefricirfj8TAaxeq/mvg6s/MRQHX3pR869+F+UBeQ77Zuw/zMOal3wPIj7D8sFTjlSXtTqzS6z8VUXJVTxZ2sNMTQAkBCJRCIIvhkgXMhQsXyty5c51CVe9Q2nLLLet3LFXfj3TNNdfI0KFDF2vHjBkjO+64o0yfPl1WX311WXnllZfaN+Q7l6p3ND3xxBN855KTWyxKRYBGMRXZ9Pv6vrDE6zSewBWuaQgwXLLkSh1b0u7EgnkYc98eICyavZr8sGeujYhXWnJ6Xen1ryfTW0mu6snCDnZ6AighAIFSCLR6uHTVVVfJ1VdfXX+EXIrHIYccsty9TzzxRJkwYcJS4XsNl5577rn6HVGHHnqo9H/n1c477yyrrrqqzJ49O+iS2tRY0pAEWZ1MjC/J0Cbf2Lf+8TqNJXCFaxoCvXf1rf8mztimmNSxvVswD2Ne+j2A/AjLD0s1XlnS7sQqvf5TESVX9WRhBzs9AZQQgEApBFo9XJo5c6ZU31/k+r1LvqZV7x566aWXlpJNmTJFttpqK5k8eXL9/Uvrrruu13CpWvyFL3xBqndTVd/lNHz48Fo/b948Ofjgg+WYY46RAw880Pe4S6xvU2NJQxJkdTIxviRDm3xj3/rH6zSWwBWuaQj03tW3/ps4Y5tiUsf2bsE8jHnp9wDyIyw/LNV4ZUm7E6v0+k9FlFzVk4Ud7PQEUEIAAqUQYLikcHJZ37n06KOPSvVP9ajeefTyyy/LQQcdVP/32muvLWPHjl0c7eGHH64HSNWX9o0fP17mz58vV1xxhay11lpyww03yLBhwxQn+6ukTY0lDUmQ1cnE+JIMbfKNfesfr9NYAle4piHQe1ff+m/ijG2KSR3buwXzMOal3wPIj7D8sFTjlSXtTqzS6z8VUXJVTxZ2sNMTQAkBCJRCgOGSwsllDZfOP/98qd5N1e1RvdNp1qxZS/zqwQcflHPPPbd+51U1TNpuu+1k2rRpMmrUKMWplpS0qbGkIQm2O8kG+JIEq8mmvvWP12lsgStc0xBguGTJlTq2pN2JBfMw5r49QFg0ezX5Yc9cGxGvtOT0utLrX0+mt5Jc1ZOFHez0BFBCAAKlEGC4VIqT/a6jTY0lDUmeSYgvefricirf+sdrF6r+a+Dqz8xFAVeGSy55EmsN+RaLpPs+MHdn1W2lbw8QFs1eTX7YM9dGxCstOb2u9PrXk2G4BLtUBPT7co/Us0MJAQjkRYDhUl5+RDtNmxpLnlSj2R51I3yJitN0M9/6x+s09sAVrmkIMFyy5EodW9LuxIJ5GHPfHiAsmr2a/LBnro2IV1pyel3p9a8nw3AJdqkI6PflHqlnhxICEMiLAMOlvPyIdpo2NZY8qUazPepG+BIVp+lmvvWP12nsgStc0xBguGTJlTq2pM1wKQZt3x4gRkzLPahJS9phsfAqjJ9GXXr9a5i4aMhVF0rd18AOdnoCKCEAgVIIMFwqxcl+19GmxpKGJM8kxJc8fXE5lW/947ULVf81cPVn5qKAK8MllzyJtYZ8i0XSfR+Yu7PqttK3BwiLZq8mP+yZayPilZacXld6/evJ9FaSq3qysIOdngBKCECgFAIMl0pxkuFSoU42d1k0is2xD43s+8ISr0OJd9fDFa5pCDBcsuRKHVvS7sSCeRhz3x4gLJq9mvywZ66NiFdacnpd6fWvJ8NwCXapCOj35R6pZ4cSAhDIi0Crh0tXXnmlXH311XLPPffkRTWD07SpseRJNYOE6XIEfMnTF5dT+dY/XrtQ9V8DV39mLgq49qbkW/8uzAfyGvLN3n2YhzEv/R5AfoTlh6Uaryxpd2KVXv+piJKrerKwg52eAEoIQKAUAq0eLpViQorraFNjSUOSIgPC98SXcIZN7eBb/3idxim4wjUNAYZLllypY0vanVgwD2Pu2wOERbNXkx/2zLUR8UpLTq8rvf71ZHoryVU9WdjBTk8AJQQgUAqBbIZLb775ptx1113ys5/9TF555RX5y1/+shTjQYMGyamnnloK+6TX0abGkoYkaSqoN8cXNbrGhb71j9dpLIMrXNMQYLhkyZU6tqTNcCkGbd8eIEZMyz2oSUvaYbHwKoyfRl16/WuYuGjIVRdK3dfADnZ6AighAIFSCGQxXHruuefkoIMOkt/+9reyaNGiZbKthkt9fX2lsE96HW1qLGlIkqaCenN8UaNrXOhb/3idxjK4wjUNAYZLllypY0vaDJdi0PbtAWLEtNyDmrSkHRYLr8L4adSl17+GiYuGXHWhxHBJTwl2sdmxHwQgkBeBLIZLBx54oDz00EPyta99TcaOHStrrrmmrLDCCnmRatlp2tRY0szlmVz4kqcvLqfyrX+8dqHqvwau/sxcFHBluOSSJ7HWkG+xSLrvA3N3Vt1W+vYAYdHs1eSHPXNtRLzSktPrSq9/PZneSnJVTxZ2sNMTQAkBCJRCIIvh0mabbSbbb7+9nHPOOaVwbfw62tRY0pA0ni5dD4Avefricirf+sdrF6r+a+Dqz8xFAVeGSy55EmsN+RaLpPs+MHdnxXApjBXqtASo5bR8B2L9pyJKrurJwg52egIoIQCBUghkMVzaeuut63csHXvssaVwbfw6fP+43OSBaUiapL/s2PiSpy8up/Ktf7x2oeq/Bq7+zFwUcGW45JInsdaQb7FIuu8Dc3dWA/GPy+RHWH5YqvHKknYnlu9rAPsT5hmRXNX7AjvY6QmghAAESiGQxXDp+OOPl5/97Gdy4403SvW9SjzCCbSpsaQhCfc7xQ74koKqzZ6+9Y/XaXyBK1zTEGC4ZMmVOrak3YkF8zDmvj1AWDR7Nflhz1wbEa+05PS60utfT6a3klzVk4Ud7PQEUEIAAqUQyGK49Morr8gBBxwgH/3oR2XatGmy1lprlcK3setoU2NJQ9JYmvQMjC95+uJyKt/6x2sXqv5r4OrPzEUBV4ZLLnkSaw35Fouk+z4wd2fVbaVvDxAWzV5Nftgz10bEKy05va70+teTYbgEu1QE9Ptyj9SzQwkBCORFIIvh0pgxY+TPf/6zvPTSSzWdVVddVVZeeeWlSFXvapo7d25eBDM9TZsaS55U80wifMnTF5dT+dY/XrtQ9V8DV39mLgq49qbkW/8uzAfyGvLN3n2YhzEv/R5AfoTlh6Uaryxpd2KVXv+piJKrerKwg52eAEoIQKAUAlkMl3bYYQdnnvfcc4/z2oG8sE2NJQ1JnpmKL3n64nIq3/rHaxeq/mvg6s/MRQFXhksueRJrDfkWi6T7PjB3Z9VtpW8PEBbNXk1+2DPXRsQrLTm9rvT615PprSRX9WRhBzs9AZQQgEApBLIYLpUCM6fraFNjSUOSU+b89Sz4kqcvLqfyrX+8dqHqvwau/sxcFHBluOSSJ7HWkG+xSLrvA3N3VgyXwlihTkuAWk7LdyDWfyqi5KqeLOxgpyeAEgIQKIUAw6VSnOx3Hb5/XG4SAw1Jk/SXHRtf8vTF5VS+9Y/XLlT918DVn5mLAq4Ml1zyJNYa8i0WSfd9YO7OaiD+cZn8CMsPSzVeWdLuxPJ9DWB/wjwjkqt6X2AHOz0BlBCAQCkEGC6V4mS/62hTY0lDkmcS4kuevricyrf+8dqFqv8auPozc1HAleGSS57EWkO+xSLpvg/M3VkxXApjhTotAWo5Ld+BWP+piJKrerKwg52eAEoIQKAUAo0Ml6ZPny6DBg2Sr3/96zJq1Cip/tvlUWlOPfVUl6UDfo3vH5ebBEZD0iT9ZcfGlzx9cTmVb/3jtQtV/zVw9WfmooBrb0q+9e/CfCCvId/s3Yd5GPPS7wHkR1h+WKrxypJ2J1bp9Z+KKLmqJws72OkJoIQABEoh0MhwafTo0fVw6fbbb5f1119fqv92eVSavr4+l6UDfk2bGksakjzTFV/y9MXlVL71j9cuVP3XwNWfmYsCrgyXXPIk1hryLRZJ931g7s6q20rfHiAsmr2a/LBnro2IV1pyel3p9a8n01tJrurJwg52egIoIQCBUgg0MlwqBV7O19GmxpKGJM9Mwpc8fXE5lW/947ULVf81cPVn5qKAK8MllzyJtYZ8i0XSfR+Yu7NiuBTGCnVaAtRyWr4Dsf5TESVX9WRhBzs9AZQQgEApBBguleJkv+vw/eNykxhoSJqkv+zY+JKnLy6n8q1/vHah6r8Grv7MXBRwZbjkkiex1pBvsUi67wNzd1YD8Y/L5EdYfliq8cqSdieW72sA+xPmGZFc1fsCO9jpCaCEAARKIcBwqRQn+11HmxpLGpI8kxBf8vTF5VS+9Y/XLlT918DVn5mLAq4Ml1zyJNYa8i0WSfd9YO7OiuFSGCvUaQlQy2n5DsT6T0WUXNWThR3s9ARQQgACpRDIarj061//Wh555BH5/e9/L2+//fZSjKvvXKC53p0AACAASURBVJoyZUop7JNeh+8fl5MeZjmb05A0SX/ZsfElT19cTuVb/3jtQtV/DVz9mbko4Nqbkm/9uzAfyGvIN3v3YR7GvPR7APkRlh+WaryypN2JVXr9pyJKrurJwg52egIoIQCBUghkMVxatGiRnHzyyXLdddfJu+++K9UQqfrZe4/3/rv6d19fXynsk15HmxpLGpKkqaDeHF/U6BoX+tY/XqexDK5wTUOA4ZIlV+rYknYnFszDmPv2AGHR7NXkhz1zbUS80pLT60qvfz2Z3kpyVU8WdrDTE0AJAQiUQiCL4dKVV14pp59+uuyzzz6y7777yl577SWTJ0+Wz33uc/LTn/5ULrnkEvnkJz8pRx99tKyzzjqlsE96HW1qLGlIkqaCenN8UaNrXOhb/3idxjK4wjUNAYZLllypY0vaDJdi0PbtAWLEtNyDmrSkHRYLr8L4adSl17+GiYuGXHWh1H0N7GCnJ4ASAhAohUAWw6U99thDhgwZIjfeeGPNdfTo0TJ16tT6n+rx/PPPyxe/+EU57LDD6qETj+UTaFNjSUOyfD+bWIEvTVCPE9O3/vE6Dvf+u8AVrmkI9N7Vt/6bOGObYlLH9m7BPIx56fcA8iMsPyzVeGVJuxOr9PpPRZRc1ZOFHez0BFBCAAKlEMhiuLTpppvW71o6/vjja64bbbSRHHzwwXLEEUcs5nzkkUfWH4l3++23l8I+6XW0qbGkIUmaCurN8UWNrnGhb/3jdRrL4ArXNAQYLllypY4taXdiwTyMuW8PEBbNXk1+2DPXRsQrLTm9rvT615PprSRX9WRhBzs9AZQQgEApBLIYLm211Vb1O5OOOeaYmmv13zvttJN885vfXMz5zDPPlGuuuUaeeOKJUtgnvY42NZY0JElTQb05vqjRNS70rX+8TmMZXOGahkDvXX3rv4kztikmdWzvFszDmJd+DyA/wvLDUo1XlrQ7sUqv/1REyVU9WdjBTk8AJQQgUAqBLIZL1XcsVd+ldN5559VcJ06cKL/+9a/ljjvukBVXXLH+WTV8Wrhwodx1112lsE96HW1qLGlIkqaCenN8UaNrXOhb/3idxjK4wjUNAYZLllypY0vanVgwD2Pu2wOERbNXkx/2zLUR8UpLTq8rvf71ZHoryVU9WdjBTk8AJQQgUAqBLIZLM2bMkO9///ty//3318OkW265RaZNmyYf//jHZeutt5bHH39cHnvsMfnqV79af+8Sj+UTaFNjSUOyfD+bWIEvTVCPE9O3/vE6Dvf+u8AVrmkI9N7Vt/6bOGObYlLH9m7BPIx56fcA8iMsPyzVeGVJuxOr9PpPRZRc1ZOFHez0BFBCAAKlEMhiuPTiiy/KfffdV38U3hprrFGzvfjii+U73/mOvPbaa7LSSivJ+PHj64/NGzJkSCnsk15HmxpLGpKkqaDeHF/U6BoX+tY/XqexDK5wTUOA4ZIlV+rYknYnFszDmPv2AGHR7NXkhz1zbUS80pLT60qvfz2Z3kpyVU8WdrDTE0AJAQiUQiCL4dKyYP7lL3+RP/7xjzJy5EgZPHhwKcxNrqNNjSUNiUlKeAfBF29k2Qh86x+v01gHV7imIcBwyZIrdWxJm+FSDNq+PUCMmJZ7UJOWtMNi4VUYP4269PrXMHHRkKsulLqvgR3s9ARQQgACpRDIYrg0ffp0GT16tEyePLkUro1fR5saSxqSxtOl6wHwJU9fXE7lW/947ULVfw1c/Zm5KODKcMklT2KtId9ikXTfB+burLqt9O0BwqLZq8kPe+baiHilJafXlV7/ejK9leSqnizsYKcngBICECiFQBbDpU033VQmTZokRx55ZClcG7+ONjWWNCSNpwvDpTwtUJ/Kt/6pQTXqnkK4wjUNAYZLllypY0vanVgwD2Pu2wOERbNXkx/2zLUR8UpLTq8rvf71ZBguwS4VAf2+3CP17FBCAAJ5EchiuLTXXnvJ+uuvL2eddVZedFp8mjY1ljyp5plo+JKnLy6n8q1/vHah6r8Grv7MXBRwZbjkkiex1pBvsUi67wNzd1bdVvr2AGHR7NXkhz1zbUS80pLT60qvfz0ZhkuwS0VAvy/3SD07lBCAQF4Eshgu3X777VJ9NN6sWbNkk002yYtQS0/TpsaSJ9U8kwxf8vTF5VS+9Y/XLlT918DVn5mLAq4Ml1zyJNYa8i0WSfd9YO7OiuFSGCvUaQlQy2n5DsT6T0WUXNWThR3s9ARQQgACpRDIYrg0Z84cueWWW+QnP/mJ7LLLLrLRRhvJGmusIYMGDVqK89ixY0thn/Q6fP+4nPQwy9mchqRJ+suOjS95+uJyKt/6x2sXqv5r4OrPzEUB196UfOvfhflAXkO+2bsP8zDmpd8DyI+w/LBU45Ul7U6s0us/FVFyVU8WdrDTE0AJAQiUQqCx4VL1TqUdd9xRxowZI6NHj64HSYsWLVqC6/uHS9Xvqv/u6+srhX3S62hTY0lDkjQV1Jvjixpd40Lf+sfrNJbBFa5pCDBcsuRKHVvS7sSCeRhz3x4gLJq9mvywZ66NiFdacnpd6fWvJ9NbSa7qycIOdnoCKCEAgVIINDZcqgZKU6dOrf+58cYbu75LqRvkcePGlcI+6XW0qbGkIUmaCurN8UWNrnGhb/3jdRrL4ArXNAQYLllypY4taTNcikHbtweIEdNyD2rSknZYLLwK46dRl17/GiYuGnLVhVL3NbCDnZ4ASghAoBQCWQyXSoGZ03W0qbGkIckpc/56FnzJ0xeXU/nWP167UPVfA1d/Zi4KuDJccsmTWGvIt1gk3feBuTurbit9e4CwaPZq8sOeuTYiXmnJ6XWl17+eTG8luaonCzvY6QmghAAESiHAcKkUJ/tdR5saSxqSPJMQX/L0xeVUvvWP1y5U/dfA1Z+ZiwKuDJdc8iTWGvItFkn3fWDuzorhUhgr1GkJUMtp+Q7E+k9FlFzVk4Ud7PQEUEIAAqUQYLhUipMMlwp1srnLolFsjn1oZIZLoQTj6KmhOBz77wJXhktpMqv7ruSbJe1OLJiHMfftAcKi2avJD3vm2oh4pSWn15Ve/3oyvZXkqp4s7GCnJ4ASAhAohUCjw6WPfOQjUv3j+hg0aJBcddVVrssH9Lo2NZY0JHmmKr7k6YvLqXzrH69dqPqvgas/MxcFXHtT8q1/F+YDeQ35Zu8+zMOYl34PID/C8sNSjVeWtDuxSq//VETJVT1Z2MFOTwAlBCBQCoFGh0u+EKvhUl9fn69sQK5vU2NJQ5JniuJLnr64nMq3/vHahar/Grj6M3NRwJXhkkuexFpDvsUi6b4PzN1ZdVvp2wOERbNXkx/2zLUR8UpLTq8rvf71ZHoryVU9WdjBTk8AJQQgUAqBRodLkydPlkmTJnmx9Hmnk9fGhS1uU2NJQ5Jn8uFLnr64nMq3/vHahar/Grj6M3NRwJXhkkuexFpDvsUi6b4PzN1ZMVwKY4U6LQFqOS3fgVj/qYiSq3qysIOdngBKCECgFAKNDpemTp0q1T884hPw/eNy/BO470hD4s7KciW+WNKOG8u3/vE6Lv/3doMrXNMQ6L2rb/03ccY2xaSO7d2CeRjz0u8B5EdYfliq8cqSdidW6fWfiii5qicLO9jpCaCEAARKIcBwqRQn+11HmxpLGpI8kxBf8vTF5VS+9Y/XLlT918DVn5mLAq4Ml1zyJNYa8i0WSfd9YO7OqttK3x4gLJq9mvywZ66NiFdacnpd6fWvJ9NbSa7qycIOdnoCKCEAgVIIMFwqxUmGS4U62dxl0Sg2xz40su8LS7wOJd5dD1e4piHAcMmSK3VsSbsTC+ZhzH17gLBo9mryw565NiJeacnpdaXXv54MwyXYpSKg35d7pJ4dSghAIC8CDJfy8iPaadrUWPKkGs32qBvhS1Scppv51j9ep7EHrnBNQ4DhkiVX6tiSNsOlGLR9e4AYMS33oCYtaYfFwqswfhp16fWvYeKiIVddKHVfAzvY6QmghAAESiHQ2HCpFIC5XkebGksakjyzCF/y9MXlVL71j9cuVP3XwNWfmYsCrgyXXPIk1hryLRZJ931g7s6q20rfHiAsmr2a/LBnro2IV1pyel3p9a8n01tJrurJwg52egIoIQCBUggwXCrFyX7X0abGkoYkzyTElzx9cTmVb/3jtQtV/zVw9WfmooArwyWXPIm1hnyLRdJ9H5i7s2K4FMYKdVoC1HJavgOx/lMRJVf1ZGEHOz0BlBCAQCkEGC55OjlhwgR59NFHZY899pAZM2YsVv/pT3+SH/zgB3LvvffKL37xC3n99ddl3XXXlc9//vMyefJkWWmllZaK9NBDD8m5554rfX19MmzYMNl+++1l2rRpMnLkSM9TLb3c94/LwQEDNqAhCYCXUIovCeEm3tq3/vE6jSFwhWsaAr139a3/Js7YppjUsb1bMA9jXvo9gPwIyw9LNV5Z0u7EKr3+UxElV/VkYQc7PQGUEIBAKQQYLnk4OWfOHDnppJPqwVH/4VI1VJoyZYp85jOfka233lpWXnnlegh16623yhZbbCGzZs2SFVZYYXG0Rx55RA444ADZcMMNZe+995YFCxbI5ZdfLmuvvbbMnj1bhg4d6nEyhktBsBB3JUCj2N7E8H1hiddpvIYrXNMQYLhkyZU6tqTdiQXzMOa+PUBYNHs1+WHPXBsRr7Tk9LrS619PpreSXNWThR3s9ARQQgACpRBguOTo5MKFC2W33Xar34V01llnLTVc+tWvflXvVL1b6f2P6p1JF154ocycOVN22mmnxb8aO3asvPzyy3LbbbfJ8OHD65/PmzdPDj74YJk+fbrsv//+jifrvqxNjSUNSZDVycT4kgxt8o196x+v01gCV7imIdB7V9/6b+KMbYpJHdu7BfMw5qXfA8iPsPywVOOVJe1OrNLrPxVRclVPFnaw0xNACQEIlEKA4ZKjkyeffLI88MADcsstt8jGG2+81HBpWds888wzsueee8phhx0mhxxySL3s+eefl1133VUOPfRQmTp16hLSnXfeWVZbbTW5/vrrHU/WfVmbGksakiCrk4nxJRna5Bv71j9ep7EErnBNQ6D3rr7138QZ2xSTOrZ3C+ZhzEu/B5AfYflhqcYrS9qdWKXXfyqi5KqeLOxgpyeAEgIQKIUAwyUHJ59++un6o+suvvhi2XbbbeuPsuv/sXjL2ua+++6Tr3zlK/XH6e277771smpAddRRR8mll14q22yzzRLS6uc//OEP5fHHH1/iY/QcjrnEkqqxXLRokYwYMcJXar7+jTfeqGNW3zvFIx8C+NKcF6NHjw4K7lv/eB2Ee5liuMJVQ8C6/jVnHEga6tje7YHMPLT+3/vjclteA2iyayDnh4ZXkxq88qcfeg/wfQ3gf8IyFeSq3lfYxWMXWv/6k6CEAAQgEEaA4dJy+L377rv1UGjkyJH1cKl6uA6XKm31MXpPPfWUzJ07V0aNGlXrL7vsMjnzzDPl5ptvrvd6/6P6efX76l1S763XWNymxpKGRONweg2+pGe8rAihjaVv/eN1Gq/hClcNAev615xxIGmoY3u3BzLz0Pqv3PLtAewdDos4kPMjjJy9Gq/8mYfeA0qvf3+ibgpy1Y1Tt1Wwi8cutP71J0EJAQhAIIzAgBkuVf8PvrffftuJ1uDBg2XIkCH12uuuu05OOeWU+ruR1ltvvfpnrsOls88+W/7jP/5DTjjhBJk4ceLi2BdccIGcd955cuedd8r666+/xJne+46mu+++W9ZZZx2n83Zb1Ka3xPNWarXNSYX4khRv0s196x+v09gBV7imIdB7V9/6b+KMbYpJHdu7BfMw5qXfA8iPsPywVOOVJe1OrNLrPxVRclVPFnaw0xNACQEIlEJgwAyXnn322fqj7Fwe48aNk9NPP10WLFggu+22m+y3335y+OGHL5a6DJf+8z//U/793/+9ftdT9ZF4739YvXOpirnFFlu4XHKja2hIGsW/zOD4kqcvLqfyfWGJ1y5U/dfA1Z+ZiwKuvSn51r8L84G8hnyzdx/mYcxLvweQH2H5YanGK0vanVil138qouSqnizsYKcngBICECiFwIAZLi1cuLD+aDqXR/UOpS233LJ+x1L1/UjXXHONDB06dLF0zJgxsuOOO8r06dNl9dVXl5VXXnmJbW+88UY57rjjZPfdd5cZM2ZI9U6o9z+W951L1TuannjiieDvXKpiMlxycZw13QjQKLY3L3xfWOJ1Gq/hCtc0BHrv6lv/TZyxTTGpY3u3YB7GvPR7APkRlh+WaryypN2JVXr9pyJKrurJwg52egIoIQCBUggMmOGSxrBDDjlEqo+n6/U48cQTZcKECYuXVB+fd9RRR8m2224rM2fOlA984ANLyZ977rn6HVGHHnqoTJ06dYnf77zzzrLqqqvK7NmzNUderGlTY0lDEmR1MjG+JEObfGPf+sfrNJbAFa5pCPTe1bf+mzhjm2JSx/ZuwTyMeen3APIjLD8s1XhlSbsTq/T6T0WUXNWThR3s9ARQQgACpRBguNTDyerdQy+99NJSK6ZMmSJbbbWVTJ48uf7+pXXXXbdeU70z6rDDDpNPfvKTcskll8iKK664zN2/8IUvSPVuqmoYNXz48HrdvHnz5OCDD5ZjjjlGDjzwwKAca1NjSUMSZHUyMb4kQ5t8Y9/6x+s0lsAVrmkI9N7Vt/6bOGObYlLH9m7BPIx56fcA8iMsPyzVeGVJuxOr9PpPRZRc1ZOFHez0BFBCAAKlEGC4pHCy23cuPfnkk/LlL39ZhgwZIkcffbQMGzZsiZ2rj9rbfPPNF//s4YcfrgdIo0ePlvHjx8v8+fPliiuukLXWWktuuOGGpfS+x2xTY0lD4uuuzXp8seGcIopv/eN1ChdE4ArXNAQYLllypY4taXdiwTyMuW8PEBbNXk1+2DPXRsQrLTm9rvT615PprSRX9WRhBzs9AZQQgEApBBguKZzsNlyqvmep+g6mZT3GjRsnp59++hK/fvDBB+Xcc8+Vvr6+epi03XbbybRp02TUqFGKUy0paVNjSUMSbHeSDfAlCVaTTX3rH6/T2AJXuKYh0HtX3/pv4oxtikkd27sF8zDmpd8DyI+w/LBU45Ul7U6s0us/FVFyVU8WdrDTE0AJAQiUQoDhUilO9ruONjWWNCR5JiG+5OmLy6l86x+vXaj6r4GrPzMXBVwZLrnkSaw15Fssku77wNydVbeVvj1AWDR7Nflhz1wbEa+05PS60utfT6a3klzVk4Ud7PQEUEIAAqUQYLhUipMMlwp1srnLolFsjn1oZN8XlngdSry7Hq5wTUOA4ZIlV+rYknYnFszDmPv2AGHR7NXkhz1zbUS80pLT60qvfz0ZhkuwS0VAvy/3SD07lBCAQF4EGC7l5Ue007SpseRJNZrtUTfCl6g4TTfzrX+8TmMPXOGahgDDJUuu1LElbYZLMWj79gAxYlruQU1a0g6LhVdh/DTq0utfw8RFQ666UOq+Bnaw0xNACQEIlEKA4VIpTva7jjY1ljQkeSYhvuTpi8upfOsfr12o+q+Bqz8zFwVcGS655EmsNeRbLJLu+8DcnVW3lb49QFg0ezX5Yc9cGxGvtOT0utLrX0+mt5Jc1ZOFHez0BFBCAAKlEGC4VIqTDJcKdbK5y6JRbI59aGTfF5Z4HUq8ux6ucE1DgOGSJVfq2JJ2JxbMw5j79gBh0ezV5Ic9c21EvNKS0+tKr389GYZLsEtFQL8v90g9O5QQgEBeBBgu5eVHtNO0qbHkSTWa7VE3wpeoOE03861/vE5jD1zhmoYAwyVLrtSxJW2GSzFo+/YAMWJa7kFNWtIOi4VXYfw06tLrX8PERUOuulDqvgZ2sNMTQAkBCJRCgOFSKU72u442NZY0JHkmIb7k6YvLqXzrH69dqPqvgas/MxcFXBkuueRJrDXkWyyS7vvA3J1Vt5W+PUBYNHs1+WHPXBsRr7Tk9LrS619PpreSXNWThR3s9ARQQgACpRBguFSKkwyXCnWyucuiUWyOfWhk3xeWeB1KvLsernBNQ4DhkiVX6tiSdicWzMOY+/YAYdHs1eSHPXNtRLzSktPrSq9/PRmGS7BLRUC/L/dIPTuUEIBAXgQYLuXlR7TTtKmx5Ek1mu1RN8KXqDhNN/Otf7xOYw9c4ZqGAMMlS67UsSVthksxaPv2ADFiWu5BTVrSDouFV2H8NOrS61/DxEVDrrpQ6r4GdrDTE0AJAQiUQoDhUilO9ruONjWWNCR5JiG+5OmLy6l86x+vXaj6r4GrPzMXBVwZLrnkSaw15Fssku77wNydVbeVvj1AWDR7Nflhz1wbEa+05PS60utfT6a3klzVk4Ud7PQEUEIAAqUQYLhUipMMlwp1srnLolFsjn1oZN8XlngdSry7Hq5wTUOA4ZIlV+rYknYnFszDmPv2AGHR7NXkhz1zbUS80pLT60qvfz0ZhkuwS0VAvy/3SD07lBCAQF4EGC7l5Ue007SpseRJNZrtUTfCl6g4TTfzrX+8TmMPXOGahgDDJUuu1LElbYZLMWj79gAxYlruQU1a0g6LhVdh/DTq0utfw8RFQ666UOq+Bnaw0xNACQEIlEKA4VIpTva7jjY1ljQkeSYhvuTpi8upfOsfr12o+q+Bqz8zFwVcGS655EmsNeRbLJLu+8DcnVW3lb49QFg0ezX5Yc9cGxGvtOT0utLrX0+mt5Jc1ZOFHez0BFBCAAKlEGC4VIqTDJcKdbK5y6JRbI59aGTfF5Z4HUq8ux6ucE1DgOGSJVfq2JJ2JxbMw5j79gBh0ezV5Ic9c21EvNKS0+tKr389GYZLsEtFQL8v90g9O5QQgEBeBBgu5eVHtNO0qbHkSTWa7VE3wpeoOE03861/vE5jD1zhmoYAwyVLrtSxJW2GSzFo+/YAMWJa7kFNWtIOi4VXYfw06tLrX8PERUOuulDqvgZ2sNMTQAkBCJRCgOFSKU72u442NZY0JHkmIb7k6YvLqXzrH69dqPqvgas/MxcFXBkuueRJrDXkWyyS7vvA3J1Vt5W+PUBYNHs1+WHPXBsRr7Tk9LrS619PpreSXNWThR3s9ARQQgACpRBguFSKkwyXCnWyucuiUWyOfWhk3xeWeB1KvLsernBNQ4DhkiVX6tiSdicWzMOY+/YAYdHs1eSHPXNtRLzSktPrSq9/PRmGS7BLRUC/L/dIPTuUEIBAXgQYLuXlR7TTtKmx5Ek1mu1RN8KXqDhNN/Otf7xOYw9c4ZqGAMMlS67UsSVthksxaPv2ADFiWu5BTVrSDouFV2H8NOrS61/DxEVDrrpQ6r4GdrDTE0AJAQiUQoDhUilO9ruONjWWNCR5JiG+5OmLy6l86x+vXaj6r4GrPzMXBVwZLrnkSaw15Fssku77wNydVbeVvj1AWDR7Nflhz1wbEa+05PS60utfT6a3klzVk4Ud7PQEUEIAAqUQYLhUipMMlwp1srnLolFsjn1oZN8XlngdSry7Hq5wTUOA4ZIlV+rYknYnFszDmPv2AGHR7NXkhz1zbUS80pLT60qvfz0ZhkuwS0VAvy/3SD07lBCAQF4EGC7l5Ue007SpseRJNZrtUTfCl6g4TTfzrX+8TmMPXOGahgDDJUuu1LElbYZLMWj79gAxYlruQU1a0g6LhVdh/DTq0utfw8RFQ666UOq+Bnaw0xNACQEIlEKA4VIpTva7jjY1ljQkeSYhvuTpi8upfOsfr12o+q+Bqz8zFwVcGS655EmsNeRbLJLu+8DcnVW3lb49QFg0ezX5Yc9cGxGvtOT0utLrX0+mt5Jc1ZOFHez0BFBCAAKlEGC4VIqTDJcKdbK5y6JRbI59aGTfF5Z4HUq8ux6ucE1DgOGSJVfq2JJ2JxbMw5j79gBh0ezV5Ic9c21EvNKS0+tKr389GYZLsEtFQL8v90g9O5QQgEBeBBgu5eVHtNO0qbHkSTWa7VE3wpeoOE03861/vE5jD1zhmoYAwyVLrtSxJW2GSzFo+/YAMWJa7kFNWtIOi4VXYfw06tLrX8PERUOuulDqvgZ2sNMTQAkBCJRCgOFSKU72u442NZY0JHkmIb7k6YvLqXzrH69dqPqvgas/MxcFXBkuueRJrDXkWyyS7vvA3J1Vt5W+PUBYNHs1+WHPXBsRr7Tk9LrS619PpreSXNWThR3s9ARQQgACpRBguFSKkwyXCnWyucuiUWyOfWhk3xeWeB1KvLsernBNQ4DhkiVX6tiSdicWzMOY+/YAYdHs1eSHPXNtRLzSktPrSq9/PRmGS7BLRUC/L/dIPTuUEIBAXgQYLuXlR7TTtKmx5Ek1mu1RN8KXqDhNN/Otf7xOYw9c4ZqGAMMlS67UsSVthksxaPv2ADFiWu5BTVrSDouFV2H8NOrS61/DxEVDrrpQ6r4GdrDTE0AJAQiUQoDhUilO9ruONjWWNCR5JiG+5OmLy6l86x+vXaj6r4GrPzMXBVwZLrnkSaw15Fssku77wNydVbeVvj1AWDR7Nflhz1wbEa+05PS60utfT6a3klzVk4Ud7PQEUEIAAqUQYLhUipMMlwp1srnLolFsjn1oZN8XlngdSry7Hq5wTUOA4ZIlV+rYknYnFszDmPv2AGHR7NXkhz1zbUS80pLT60qvfz0ZhkuwS0VAvy/3SD07lBCAQF4EGC7l5Ue007SpseRJNZrtUTfCl6g4TTfzrX+8TmMPXOGahgDDJUuu1LElbYZLMWj79gAxYlruQU1a0g6LhVdh/DTq0utfw8RFQ666UOq+Bnaw0xNACQEIlEKA4VIpTva7jjY1ljQkeSYhvuTpi8upfOsfr12o+q+Bqz8zFwVcGS655EmsNeRbLJLu+8DcnVW3lb49QFg0ezX5Yc9cGxGvtOT0utLrX0+mt5Jc1ZOFHez0BFBCAAKlEGC4VIqTDJcKdbK5y6JRbI59aGTfF5Z4HUq8ux6ucE1DgOGSJVfq2JJ2JxbMw5j79gBh0ezV5Ic9c21EvNKSHMwDtgAAIABJREFU0+tKr389GYZLsEtFQL8v90g9O5QQgEBeBBgu5eVHtNO0qbHkSTWa7VE3wpeoOE03861/vE5jD1zhmoYAwyVLrtSxJW2GSzFo+/YAMWJa7kFNWtIOi4VXYfw06tLrX8PERUOuulDqvgZ2sNMTQAkBCJRCgOFSKU72u442NZY0JHkmIb7k6YvLqXzrH69dqPqvgas/MxcFXBkuueRJrDXkWyyS7vvA3J1Vt5W+PUBYNHs1+WHPXBsRr7Tk9LrS619PpreSXNWThR3s9ARQQgACpRBguFSKkwyXCnWyucuiUWyOfWhk3xeWeB1KvLsernBNQ4DhkiVX6tiSdicWzMOY+/YAYdHs1eSHPXNtRLzSktPrSq9/PRmGS7BLRUC/L/dIPTuUEIBAXgQYLuXlR7TTtKmx5Ek1mu1RN8KXqDhNN/Otf7xOYw9c4ZqGAMMlS67UsSVthksxaPv2ADFiWu5BTVrSDouFV2H8NOrS61/DxEVDrrpQ6r4GdrDTE0AJAQiUQoDhUilO9ruONjWWNCR5JiG+5OmLy6l86x+vXaj6r4GrPzMXBVwZLrnkSaw15Fssku77wNydVbeVvj1AWDR7Nflhz1wbEa+05PS60utfT6a3klzVk4Ud7PQEUEIAAqUQYLhUipMMlwp1srnLolFsjn1oZN8XlngdSry7Hq5wTUOA4ZIlV+rYknYnFszDmPv2AGHR/n/27gX8qqpM4P/7g+QPQTowBCOIZj0pTw6ISiiVpaKFFxRIyAAxL6FycTBTkhzRRkkqbBQocUxD8hYXsQw0QbymRoyJjow4eBcNL48Yl0SQ/7P3+fP7C5zf+a31rr3W3nud73keH2dkv3vt83nfd7N+v7d9Tvho6iO8uXZFcqWV08fF3v96GYZL2PkS0J+Xe6TejkgEECiWAMOlYuUjs6sp08aSv1QzS3umJyIvmXIGPZlt/5NrP+nBFVc/AgyXQrrSxyG1GS5loW27B8hizZDnoCdDarutRa7c/DTRsfe/xsQkhlo1Uap+DHbY6QWIRACBWAQYLsWSyZ3eR5k2lmxIilmE5KWYeTG5Ktv+J9cmqvbH4GpvZhKBK8MlkzrJ6hjqLStJ8/Ngbm5V7UjbPYDbauGjqY/w5toVyZVWTh8Xe//rZWpHUqt6Weyw0wsQiQACsQgwXIolkwyXIs1kfm+LjWJ+9q4r2/5gSa5dxavH44qrHwGGSyFd6eOQ2pW1MHczt90DuK0WPpr6CG+uXZFcaeX0cbH3v16G4RJ2vgT05+UeqbcjEgEEiiXAcKlY+cjsasq0seQv1czSnumJyEumnEFPZtv/5NpPenDF1Y8Aw6WQrvRxSG2GS1lo2+4Bslgz5DnoyZDabmuRKzc/TXTs/a8xMYmhVk2Uqh+DHXZ6ASIRQCAWAYZLsWRyp/dRpo0lG5JiFiF5KWZeTK7Ktv/JtYmq/TG42puZRODKcMmkTrI6hnrLStL8PJibW1U70nYP4LZa+GjqI7y5dkVypZXTx8Xe/3qZ2pHUql4WO+z0AkQigEAsAgyXYskkw6VIM5nf22KjmJ+968q2P1iSa1fx6vG44upHgOFSSFf6OKR2ZS3M3cxt9wBuq4WPpj7Cm2tXJFdaOX1c7P2vl2G4hJ0vAf15uUfq7YhEAIFiCTBcKlY+MruaMm0s+Us1s7RneiLykiln0JPZ9j+59pMeXHH1I8BwKaQrfRxSm+FSFtq2e4As1gx5DnoypLbbWuTKzU8THXv/a0xMYqhVE6Xqx2CHnV6ASAQQiEWA4VIsmdzpfZRpY8mGpJhFSF6KmReTq7Ltf3Jtomp/DK72ZiYRuDJcMqmTrI6h3rKSND8P5uZW1Y603QO4rRY+mvoIb65dkVxp5fRxsfe/XqZ2JLWql8UOO70AkQggEIsAw6VYMslwKdJM5ve22CjmZ++6su0PluTaVbx6PK64+hFguBTSlT4OqV1ZC3M3c9s9gNtq4aOpj/Dm2hXJlVZOHxd7/+tlGC5h50tAf17ukXo7IhFAoFgCDJeKlY/MrqZMG0v+Us0s7ZmeiLxkyhn0ZLb9T679pAdXXP0IMFwK6Uofh9RmuJSFtu0eIIs1Q56Dngyp7bYWuXLz00TH3v8aE5MYatVEqfox2GGnFyASAQRiEWC4FEsmd3ofZdpYsiEpZhGSl2LmxeSqbPufXJuo2h+Dq72ZSQSuDJdM6iSrY6i3rCTNz4O5uVW1I233AG6rhY+mPsKba1ckV1o5fVzs/a+XqR1JreplscNOL0AkAgjEIsBwKZZMMlyKNJP5vS02ivnZu65s+4MluXYVrx6PK65+BBguhXSlj0NqV9bC3M3cdg/gtlr4aOojvLl2RXKlldPHxd7/ehmGS9j5EtCfl3uk3o5IBBAolgDDpWLlI7OrKdPGkr9UM0t7piciL5lyBj2Zbf+Taz/pwRVXPwIMl0K60schtRkuZaFtuwfIYs2Q56AnQ2q7rUWu3Pw00bH3v8bEJIZaNVGqfgx22OkFiEQAgVgEGC7Fksmd3keZNpZsSIpZhOSlmHkxuSrb/ifXJqr2x+Bqb2YSgSvDJZM6yeoY6i0rSfPzYG5uVe1I2z2A22rho6mP8ObaFcmVVk4fF3v/62VqR1KrelnssNMLEIkAArEIMFyyzOSIESNk2bJlMmDAAPnZz362Q/SUKVPSP3v11Vdl48aN0rlzZ+nTp4+MHj1a9tprr11Weuyxx+Saa66RlStXSps2beTII4+UCy+8UDp06GB5VbseXqaNJRsS53R7OQF58cIa5KS2/U+u/aQFV1z9CNQ+q23/53GNZVqTPg6fLczdzGO/B1AfbvURMppchdSurBV7//sSpVb1sthhpxcgEgEEYhFguGSRyQULFsjll1+eDo6qDZeGDx8u+++/v+y9997Stm3bdMg0d+5c2bJli8ybN0+6devWuNqf//xnOf3009PjTz75ZHn33XflxhtvlC5duqQxrVu3trgyhktOWARXFWCjWN7CsP3Bklz7yTWuuPoRYLgU0pU+DqldWQtzN3PbPYDbauGjqY/w5toVyZVWTh8Xe//rZWpHUqt6Weyw0wsQiQACsQgwXDLM5Pvvvy/HHnusnHbaaTJ16tSqw6Vqp3r66afT4dGoUaPkggsuaDxk4MCBsm7dOvnDH/4gn/zkJ9P//uCDD6bHXXzxxfKd73zH8MqqH1amjSUbEqdUewsmL95ovZ/Ytv/JtZ+U4IqrH4HaZ7Xt/zyusUxr0sfhs4W5m3ns9wDqw60+QkaTq5DalbVi739fotSqXhY77PQCRCKAQCwCDJcMM/mjH/1IHn30Ufn9738vPXr0MB4uJU8k9e3bV0455ZT0qafk9eKLL0r//v1l3LhxMnbs2B2u4Otf/7rsscceMmfOHMMrY7jkBEVwkwJsFMtbHLY/WJJrP7nGFVc/AgyXQrrSxyG1K2th7mZuuwdwWy18NPUR3ly7IrnSyunjYu9/vUztSGpVL4sddnoBIhFAIBYBhksGmfyf//mf9Omj6667Tr72ta+lH2VX7WPxklN99NFH8t5778nWrVvl9ddflxkzZshDDz2U/vvoo49OV0sGVN///vflhhtukMMPP3yHK0j++7333it//etfpWXLlgZXV/2QMm0s2ZCo0+w1kLx45fV6ctv+J9d+0oErrn4Eap/Vtv/zuMYyrUkfh88W5m7msd8DqA+3+ggZTa5CalfWir3/fYlSq3pZ7LDTCxCJAAKxCDBcaiaTybAoeeqoQ4cO6XApedUaLr322mvSr1+/xrO2b99ezjnnnB0+5u5Xv/qV/OQnP5Hf/e536bk+/kr+e/LnyVNSHTt2VNdZsrHctm1b+t1PRX9t2rQpvcQ2bdoU/VLr6vrIS37p7t69u9Pitv1Prp24mwzGFVeNQOj+11xjPcXQx+GzXc/mrv2//ZfLZfkZQFNd9VwfGq88Y8iVvb7rPcD2ZwD7K4wzglrV5xW77Oxc+19/JUQigAACbgJ1M1xKfsjavHmzkVaLFi1kt912S4+944475Iorrki/G2nvvfdO/1ut4dIHH3wgy5cvT9davXq13H333emwafTo0ZKcN3klTzFde+21cs8998i+++67wzVdc8018otf/EKWLFkie+21l9H1VjuoTBtLNiTqNHsNJC9eeWue3HVjadv/5NpPrnHFVSMQuv8111hPMfRx+GzXs7lr/yfZst0DhM+w24r1XB9ucuGjyZW9ues9IPb+txc1i6BWzZyqHYVddnau/a+/EiIRQAABN4G6GS6tWrUq/Sg7k9egQYPkqquukuT7ko499lj59re/LePHj28MrTVc2vn8a9askRNOOEFOPfVUOf/889M/DvXkUrLWwQcfbPKWcz2GR6lz5W9ycfJSzLyYXJXtR2KQaxNV+2NwtTczicC1tpJt/5uY1/Mx1Fv47GPuZh77PYD6cKuPkNHkKqR2Za3Y+9+XKLWql8UOO70AkQggEItA3QyX3n//fVm8eLFR3pInlHr37p0+sZR8P9Itt9wirVu3boxNnkRKvj/p4osvln/6p3+Sdu3a1TzvWWedJc8995w8/PDD6XHNfedS8kTTU089xXcuGWWLg3wJsFH0Jev/vLY/WJJrPznBFVc/AgyXQrrSxyG1K2th7mZuuwdwWy18NPUR3ly7IrnSyunjYu9/vUztSGpVL4sddnoBIhFAIBaBuhkuaRKWfJRd8vF0tV7//u//LiNGjKh5TPLU0ooVK9KBUfJ64YUX0ieixo0bJ2PHjt0h9utf/7rsvvvuMnfuXM0lN8YkH82XvBoaGpzOEyI4+cjCslxrCI+irEFe8stE0rcHHXSQ+gJs+59cq6lrBuKKq0YgdP9rrrGeYujj8NmuZ3PX/k+yZbsHCJ9htxXruT7c5MJHkyt7c9d7QOz9by9qFkGtmjlVOwq77Oxc+19/JUQigAACbgIMl2r4JcOgt956a5cjxowZI3369JHTTjst/f6lbt26SfJkVJs2bRq/q2l7UPK/5BgyZIj06tVLZs+e3Xiuk046KY1Jvsvpk5/8ZPrfH3zwQRk1apRMmDBBzjjjDKfMsrF04iMYgVwFXDeW9H+u6WNxBJwE6H8nPoIRKLWAa/8nb549QKlLgIuvcwHXewD9X+cFxNsvtYBr/5f6zXPxCCBQagGGS4r0VfvOpeQj9y677DLp37+/7LPPPulH2j3//POyYMGCdIVZs2ZJz549G1d7/PHH0wFS8qV9yfDpnXfekZtuukk6d+4s8+bNSwdVvBBAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQACBogkwXFJkpNpw6ZVXXpFf/vKX6Zdorl27Vj788EPp1KmTHHrooenTSPvuu+8uK/3pT3+Sa665RlauXJkOk4444gi58MILpWPHjoqrIgQBBBBAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQ8C/AcMm/MSsggAACCCCAAAIIIIAAAggggAACCCCAAAIIIIAAAtEIMFyKJpW8EQQQQAABBBBAAAEEEEAAAQQQQAABBBBAAAEEEEDAvwDDJf/GrIAAAggggAACCCCAAAIIIIAAAggggAACCCCAAAIIRCPAcCmaVPJGEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABBBBAAAH/AgyX/BuzAgIIIIAAAggggAACCCCAAAIIIIAAAggggAACCCAQjQDDpWhSyRtBAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABBPwLMFzyb8wKCCCAAAIIIIAAAggggAACCCCAAAIIIIAAAggggEA0AgyXokklbwQBBBBAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQ8C/AcMm/MSsggAACCCCAAAIIIIAAAggggAACCCCAAAIIIIAAAtEIMFyKJpW8EQQQQAABBBBAAAEEEEAAAQQQQAABBBBAAAEEEEDAvwDDJf/GrIAAAggggAACCCCAAAIIIIAAAggggAACCCCAAAIIRCPAcCmaVPJGEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABBBBAAAH/AgyX/BuzAgIIIIAAAggggAACCCCAAAIIIIAAAggggAACCCAQjQDDpWhSyRtBAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABBPwLMFzyb8wKCCCAAAIIIIAAAggggAACCCCAAAIIIIAAAggggEA0AgyXokklbwQBBBBAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQ8C/AcMm/MSsggAACCCCAAAIIIIAAAggggAACCCCAAAIIIIAAAtEIMFyKJpW8EQQQQAABBBBAAAEEEEAAAQQQQAABBBBAAAEEEEDAvwDDJf/GrIAAAggggAACCCCAAAIIIIAAAggggAACCCCAAAIIRCPAcCmaVO74Rp588sn0Pxx00EGRvkPeFgIINCVA/1MbCNSvAP1fv7nnnSOQCHAPoA4QqF8B+r9+c887RwABBBBAIC8Bhkt5yXte97//+7/TFQ4++GDPK7mf/n//93/Tk3Tv3t39ZJwhMwHykhll8BPZ9j+59pMiXHH1I1D7rLb9n8c1lmlN+jh8tjB3M4/9HkB9uNVHyGhyFVK7slbs/e9LlFrVy2KHnV6ASAQQiEWA4VIsmdzpfZRpY8mGpJhFSF6KmReTq7Ltf3Jtomp/DK72ZiYRuDJcMqmTrI6h3rKSND8P5uZW1Y603QO4rRY+mvoIb65dkVxp5fRxsfe/XqZ2JLWql8UOO70AkQggEIsAw6VYMslwKdJM5ve22CjmZ++6su0PluTaVbx6PK64+hFguBTSlT4OqV1ZC3M3c9s9gNtq4aOpj/Dm2hXJlVZOHxd7/+tlGC5h50tAf17ukXo7IhFAoFgCDJeKlY/MrqZMG0v+Us0s7ZmeiLxkyhn0ZLb9T679pAdXXP0IMFwK6Uofh9RmuJSFtu0eIIs1Q56Dngyp7bYWuXLz00TH3v8aE5MYatVEqfox2GGnFyASAQRiEWC4FEsmd3ofZdpYsiEpZhGSl2LmxeSqbPufXJuo2h+Dq72ZSQSuDJdM6iSrY6i3rCTNz4O5uVW1I233AG6rhY+mPsKba1ckV1o5fVzs/a+XqR1JreplscNOL0AkAgjEIsBwKZZMMlyKNJP5vS02ivnZu65s+4MluXYVrx6PK65+BBguhXSlj0NqV9bC3M3cdg/gtlr4aOojvLl2RXKlldPHxd7/ehmGS9j5EtCfl3uk3o5IBBAolgDDpWLlI7OrKdPGkr9UM0t7piciL5lyBj2Zbf+Taz/pwRVXPwIMl0K60schtRkuZaFtuwfIYs2Q56AnQ2q7rUWu3Pw00bH3v8bEJIZaNVGqfgx22OkFiEQAgVgEGC7Fksmd3keZNpZsSIpZhOSlmHkxuSrb/ifXJqr2x+Bqb2YSgSvDJZM6yeoY6i0rSfPzYG5uVe1I2z2A22rho6mP8ObaFcmVVk4fF3v/62VqR1KrelnssNMLEIkAArEIMFyKJZMMlyLNZH5vi41ifvauK9v+YEmuXcWrx+OKqx8BhkshXenjkNqVtTB3M7fdA7itFj6a+ghvrl2RXGnl9HGx979ehuESdr4E9OflHqm3IxIBBIolwHCpWPnI7GrKtLHkL9XM0p7pichLppxBT2bb/+TaT3pwxdWPAMOlkK70cUhthktZaNvuAbJYM+Q56MmQ2m5rkSs3P0107P2vMTGJoVZNlKofgx12egEiEUAgFgGGS7Fkcqf3UaaNJRuSYhYheSlmXkyuyrb/ybWJqv0xuNqbmUTgynDJpE6yOoZ6y0rS/DyYm1tVO9J2D+C2Wvho6iO8uXZFcqWV08fF3v96mdqR1KpeFjvs9AJEIoBALAIMl2LJJMOlSDOZ39tio5ifvevKtj9YkmtX8erxuOLqR4DhUkhX+jikdmUtzN3MbfcAbquFj6Y+wptrVyRXWjl9XOz9r5dhuISdLwH9eblH6u2IRACBYgkwXCpWPjK7mjJtLPlLNbO0Z3oi8pIpZ9CT2fY/ufaTHlxx9SPAcCmkK30cUpvhUhbatnuALNYMeQ56MqS221rkys1PEx17/2tMTGKoVROl6sdgh51egEgEEIhFgOFSLJnc6X2UaWPJhqSYRUheipkXk6uy7X9ybaJqfwyu9mYmEbgyXDKpk6yOod6ykjQ/D+bmVtWOtN0DuK0WPpr6CG+uXZFcaeX0cbH3v16mdiS1qpfFDju9AJEIIBCLAMOlWDLJcCnSTOb3ttgo5mfvurLtD5bk2lW8ejyuuPoRYLgU0pU+DqldWQtzN3PbPYDbauGjqY/w5toVyZVWTh8Xe//rZRguYedLQH9e7pF6OyIRQKBYAgyXipWPzK6mTBtL/lLNLO2Znoi8ZMoZ9GS2/U+u/aQHV1z9CDBcCulKH4fUZriUhbbtHiCLNUOeg54Mqe22Frly89NEx97/GhOTGGrVRKn6MdhhpxcgEgEEYhFguBRLJnd6H2XaWLIhKWYRkpdi5sXkqmz7n1ybqNofg6u9mUkErgyXTOokq2Oot6wkzc+DublVtSNt9wBuq4WPpj7Cm2tXJFdaOX1c7P2vl6kdSa3qZbHDTi9AJAIIxCLAcCmWTDJcijST+b0tNor52buubPuDJbl2Fa8ejyuufgQYLoV0pY9DalfWwtzN3HYP4LZa+GjqI7y5dkVypZXTx8Xe/3oZhkvY+RLQn5d7pN6OSAQQKJZAXQyXtmzZIjNnzpR58+bJW2+9JV27dpURI0bI8OHDpaGhoWZGFi5cKA888ICsWLFCXnrpJenUqZM89NBDzWZx8+bNMmDAgDTmnHPOkfPPP7/JmHfffVeOPfZYee+99+SKK66QIUOGNHv+5g4o08aSv1Sby2Y+f05e8nHPYlXb/ifXWajveg5ccfUjUPustv2fxzWWaU36OHy2MHczj/0eQH241UfIaHIVUruyVuz970uUWtXLYoedXoBIBBCIRaAuhkuXXHKJzJkzR4YOHSo9e/aURx55RO655x4ZN26cjB07tmYuTz31VHnmmWfkgAMOSAdFLVq0MBouzZgxQ2644QbZuHFjs8OliRMnyqJFi9JjGS7F0lrlfx9sFMubQ9sfLMm1n1zjiqsfAYZLIV3p45DalbUwdzO33QO4rRY+mvoIb65dkVxp5fRxsfe/XqZ2JLWql8UOO70AkQggEItA9MOllStXysCBA+WMM86QCRMmNOZt/PjxsmTJkvSf5Gmkpl5vvPFG+uctW7aUZND08ssvNztcevXVV+WEE06QMWPGyNSpU2sOl5INYPIEVXI9V199NcOlWDorgvfBRrG8SbT9wZJc+8k1rrj6EWC4FNKVPg6pzXApC23bPUAWa4Y8Bz0ZUtttLXLl5qeJjr3/NSYmMdSqiVL1Y7DDTi9AJAIIxCIQ/XApGdgkH4m3dOlS6dKlS2Peli9fLsOGDZNJkyal/zZ5mQ6Xzj77bNmwYYNcddVV0q9fvyaHS1u3bpXBgwdL9+7d03+PHDmS4ZJJIjgmiAAbxSDMXhax/cGSXHtJA//rez+suDbjatv/ntIUzWm5P4ZPJeZu5rHfA6gPt/oIGU2uQmpX1oq9/32JUqt6Weyw0wsQiQACsQhEP1xKnlhatWpV+lF4H38l34l04IEHpkOdK6+80iifJsOlxYsXy3nnnSd33nmntG3btuZwadasWXLttdemH9H3wgsvMFzq3t0oDxwURoCNYhhnH6vY/mBJrn1kgY928qOKa3Outv3f3Pnq/c+5P4avAMzdzGO/B1AfbvURMppchdRmuOSiTa3q9bDDTi9AJAIIxCIQ/XAp+Xi6Vq1ayfz583fJWd++fdPvUkq+G8nk1dxwadOmTXL88cfLUUcdJcn3PL322mtNDpfWrl0rxx57bPqdT6effro88cQTmQ+Xtm3blg64iv5K3JJXmzZtin6pdXV95CW/dCdPM7q8kl8s2fQ/uXbRbjoWV1w1AqH7X3ON9RRDH4fPdj2bu/Z/ki3bPUD4DLutWM/14SYXPppc2Zu73gNi7397UbMIatXMqdpR2GVn59r/+ishEgEEEHATiH64dPTRR0vHjh3l9ttv30XqiCOOkG7dusns2bONFJsbLiUfwTdnzhy59957Zffdd685XLrgggsk+T6ou+66S3bbbTeGSwyXjGow5EFsFENq77iW68bS9gdLcu0n17jiqhEI3f+aa6ynGPo4fLbr2dy1/xkuha9XVmxaoJ57WVsXrvcA258BtNcZWxy1qs8odtnZufa//kqIRAABBNwEoh8uhXpyKflYuxNPPDH9DqchQ4akWWnqyaXHH39cTjvtNPn1r38tydNTycvHk0vJeQ8++GC3CgkQzaPUAZAVS5AXBVpBQmw/Eodc+0kcrrj6Eah9Vtv+z+May7QmfRw+W5i7mcd+D6A+3OojZDS5CqldWSv2/vclSq3qZbHDTi9AJAIIxCIQ/XCpue9cGjRokEyePNkon7WeXDr33HPT701KPmKvoaEhPd+bb74pw4cPlxEjRqQffZc8QdW6det0CLXHHnvIj3/848Z1n3rqKfne974nyRNNxx13nHTu3Dl9okn7KtPGkg2JNst+48iLX1+fZ7ftf3LtJxu44upHgOFSSFf6OKR2ZS3M3cxt9wBuq4WPpj7Cm2tXJFdaOX1c7P2vl6kdSa3qZbHDTi9AJAIIxCIQ/XBp6tSpcv3118vSpUulS5cujXlbvny5DBs2TC699NJ0AGTyqjVcOumkkxp/GG7qXNddd50ceeSR0rt3b/n73/9ec8nf//73st9++5lcVtVjyrSxZEOiTrPXQPLildfryW37n1z7SQeuuPoRqH1W2/7P4xrLtCZ9HD5bmLuZx34PoD7c6iNkNLkKqV1ZK/b+9yVKreplscNOL0AkAgjEIhD9cOnZZ5+V5Omk5AmmCRMmNOZt/PjxsnjxYlmyZEn6lFDyWbFr1qyR9u3bS4cOHarmt9ZwKfl0C9zGAAAgAElEQVSou/Xr1+8Q984776TDq/79+8uAAQOkV69e6dNLDzzwgGzZsmWHY1etWiXXXHNNOuj60pe+JIcddpi0a9dOXWdl2liyIVGn2WsgefHK6/Xktv1Prv2kA1dc/QgwXArpSh+H1K6shbmbue0ewG218NHUR3hz7YrkSiunj4u9//UytSOpVb0sdtjpBYhEAIFYBKIfLiWJmjhxosyfP1+GDh0qPXr0kEcffVQWLVokY8eOlXHjxqW53P6dRx//b8l/X7ZsWfpP8po7d66sW7dOzjzzzPT/T56EGjhwYJO10NR3LlUL4DuXRPgCw2LdVtgoFisfNldj+4MlubbRNT8WV3MrmyNxra1l2/829vV4LPUWPuuYu5nHfg+gPtzqI2Q0uQqpXVkr9v73JUqt6mWxw04vQCQCCMQiUBfDpQ8//FBmzpyZDpjWrl0rXbt2TZ8QSp5E2v79SE0Nl6ZNmybTp0+vmu8+ffrI7NmzGS45dgMbEkdAT+HkxRNsgNPa/mBJrv0kBVdc/QgwXArpSh+H1K6shbmbue0ewG218NHUR3hz7YrkSiunj4u9//UytSOpVb0sdtjpBYhEAIFYBOpiuBRLsmzeR5k2lmxIbDIb7ljyEs4665Vs+59cZ50BfkHqRxRXE1fb/jc5Zz0fw/0xfPYxdzOP/R5AfbjVR8hochVSu7JW7P3vS5Ra1ctih51egEgEEIhFgOFSLJnc6X2UaWPJhqSYRUheipkXk6uy7X9ybaJqfwyu9mYmEbjWVrLtfxPzej6GeguffczdzGO/B1AfbvURMppchdRmuOSiTa3q9bDDTi9AJAIIxCLAcCmWTDJcijST+b0tNor52buubPuLJXLtKl49Hldc/QgwXArpSh+H1K6shbmbue0ewG218NHUR3hz7YrkSiunj4u9//UytSOpVb0sdtjpBYhEAIFYBBguxZJJhkuRZjK/t8VGMT9715Vtf7Ak167iDJf8COKqcbXtf80a9RTD/TF8tjF3M4/9HkB9uNVHyGhyFVK7slbs/e9LlFrVy2KHnV6ASAQQiEWA4VIsmWS4FGkm83tbbBTzs3dd2fYHS3LtKs4QxI8grhpX2/7XrFFPMdwfw2cbczfz2O8B1IdbfYSMJlchtRkuuWhTq3o97LDTCxCJAAKxCDBciiWTDJcizWR+b4uNYn72rivb/mKJXLuKMwTxI4irxtW2/zVr1FMM98fw2cbczTz2ewD14VYfIaPJVUhthksu2tSqXg877PQCRCKAQCwCDJdiySTDpUgzmd/bYqOYn73ryra/WCLXruIMQfwI4qpxte1/zRr1FMP9MXy2MXczj/0eQH241UfIaHIVUpvhkos2tarXww47vQCRCCAQiwDDpVgyyXAp0kzm97bYKOZn77qy7S+WyLWrOEMQP4K4alxt+1+zRj3FcH8Mn23M3cxjvwdQH271ETKaXIXUZrjkok2t6vWww04vQCQCCMQiwHAplkwyXIo0k/m9LTaK+dm7rmz7iyVy7SrOEMSPIK4aV9v+16xRTzHcH8NnG3M389jvAdSHW32EjCZXIbUZLrloU6t6Peyw0wsQiQACsQgwXIolkwyXIs1kfm+LjWJ+9q4r2/5iiVy7ijME8SOIq8bVtv81a9RTDPfH8NnG3M089nsA9eFWHyGjyVVIbYZLLtrUql4PO+z0AkQigEAsAgyXYskkw6VIM5nf22KjmJ+968q2v1gi167iDEH8COKqcbXtf80a9RTD/TF8tjF3M4/9HkB9uNVHyGhyFVKb4ZKLNrWq18MOO70AkQggEIsAw6VYMslwKdJM5ve22CjmZ++6su0vlsi1qzhDED+CuGpcbftfs0Y9xXB/DJ9tzN3MY78HUB9u9REymlyF1Ga45KJNrer1sMNOL0AkAgjEIsBwKZZMMlyKNJP5vS02ivnZu65s+4slcu0qzhDEjyCuGlfb/tesUU8x3B/DZxtzN/PY7wHUh1t9hIwmVyG1GS65aFOrej3ssNMLEIkAArEIMFyKJZMMlyLNZH5vi41ifvauK9v+Yolcu4ozBPEjiKvG1bb/NWvUUwz3x/DZxtzNPPZ7APXhVh8ho8lVSG2GSy7a1KpeDzvs9AJEIoBALAIMl2LJJMOlSDOZ39tio5ifvevKtr9YIteu4gxB/AjiqnG17X/NGvUUw/0xfLYxdzOP/R5AfbjVR8hochVSm+GSiza1qtfDDju9AJEIIBCLAMOlWDLJcCnSTOb3ttgo5mfvurLtL5bItas4QxA/grhqXG37X7NGPcVwfwyfbczdzGO/B1AfbvURMppchdRmuOSiTa3q9bDDTi9AJAIIxCLAcCmWTDJcijST+b0tNor52buubPuLJXLtKs4QxI8grhpX2/7XrFFPMdwfw2cbczfz2O8B1IdbfYSMJlchtRkuuWhTq3o97LDTCxCJAAKxCDBciiWTDJcizWR+b4uNYn72rivb/mKJXLuKMwTxI4irxtW2/zVr1FMM98fw2cbczTz2ewD14VYfIaPJVUhthksu2tSqXg877PQCRCKAQCwCDJdiySTDpUgzmd/bYqOYn73ryra/WCLXruIMQfwI4qpxte1/zRr1FMP9MXy2MXczj/0eQH241UfIaHIVUpvhkos2tarXww47vQCRCCAQiwDDpVgyyXAp0kzm97bYKOZn77qy7S+WyLWrOEMQP4K4alxt+1+zRj3FcH8Mn23M3cxjvwdQH271ETKaXIXUZrjkok2t6vWww04vQCQCCMQiwHAplkwyXIo0k/m9LTaK+dm7rmz7iyVy7SrOEMSPIK4aV9v+16xRTzHcH8NnG3M389jvAdSHW32EjCZXIbUZLrloU6t6Peyw0wsQiQACsQgwXIolkwyXIs1kfm+LjWJ+9q4r2/5iiVy7ijME8SOIq8bVtv81a9RTDPfH8NnG3M089nsA9eFWHyGjyVVIbYZLLtrUql4PO+z0AkQigEAsAgyXYskkw6VIM5nf22KjmJ+968q2v1gi167iDEH8COKqcbXtf80a9RTD/TF8tjF3M4/9HkB9uNVHyGhyFVKb4ZKLNrWq18MOO70AkQggEIsAw6VYMslwKdJM5ve22CjmZ++6su0vlsi1qzhDED+CuGpcbftfs0Y9xXB/DJ9tzN3MY78HUB9u9REymlyF1Ga45KJNrer1sMNOL0AkAgjEIsBwKZZMMlyKNJP5vS02ivnZu65s+4slcu0qzhDEjyCuGlfb/tesUU8x3B/DZxtzN/PY7wHUh1t9hIwmVyG1GS65aFOrej3ssNMLEIkAArEIMFyKJZMMlyLNZH5vi41ifvauK9v+Yolcu4ozBPEjiKvG1bb/NWvUUwz3x/DZxtzNPPZ7APXhVh8ho8lVSG2GSy7a1KpeDzvs9AJEIoBALAJ1MVzasmWLzJw5U+bNmydvvfWWdO3aVUaMGCHDhw+XhoaGmrlcuHChPPDAA7JixQp56aWXpFOnTvLQQw81m//NmzfLgAED0phzzjlHzj///MaYN998U+bPn5+e58UXX5Tk+vbZZx/51re+JSeffLK0bNmy2fM3d0CZfrBkQ9JcNvP5c/KSj3sWq9r2P7nOQn3Xc+CKqx+B2me17f88rrFMa9LH4bOFuZt57PcA6sOtPkJGk6uQ2gyXXLSpVb0edtjpBYhEAIFYBOpiuHTJJZfInDlzZOjQodKzZ0955JFH5J577pFx48bJ2LFja+by1FNPlWeeeUYOOOCAdFDUokULo+HSjBkz5IYbbpCNGzfuMlz6zW9+I1OmTJEjjzxSDjnkENltt93kwQcfTIdYJ5xwgkydOtW5vsr0gyUbEud0ezkBefHCGuSktv1Prv2kBVdc/QgwXArpSh+H1K6shbmbue0ewG218NHUR3hz7YrkSiunj4u9//UytSOpVb0sdtjpBYhEAIFYBKIfLq1cuVIGDhwoZ5xxhkyYMKExb+PHj5clS5ak/yRPIzX1euONN9I/T54mSgZNL7/8crPDpVdffTUdEo0ZMyYdFO385NKqVaukQ4cO0rFjxx2Wveiii+Suu+6SO++8U77whS841ViZNpZsSJxS7S2YvHij9X5i2/4n135SgiuufgRqn9W2//O4xjKtSR+Hzxbmbuax3wOoD7f6CBlNrkJqV9aKvf99iVKrelnssNMLEIkAArEIRD9cuvrqq9OPxFu6dKl06dKlMW/Lly+XYcOGyaRJk9J/m7xMh0tnn322bNiwQa666irp16/fLsOlpta6//775dxzz5Wf/OQnctJJJ5lcUpPHlGljyYbEKdXegsmLN1rvJ7btf3LtJyW44upHgOFSSFf6OKR2ZS3M3cxt9wBuq4WPpj7Cm2tXJFdaOX1c7P2vl6kdSa3qZbHDTi9AJAIIxCIQ/XApeWIpeVIo+Si8j7+S70Q68MADZfDgwXLllVca5dNkuLR48WI577zz0qeP2rZtazVcuv3229NhV/JxeocffrjRNTV1UJk2lmxInFLtLZi8eKP1fmLb/ifXflKCK65+BGqf1bb/87jGMq1JH4fPFuZu5rHfA6gPt/oIGU2uQmpX1oq9/32JUqt6Weyw0wsQiQACsQhEP1xKPp6uVatWMn/+/F1y1rdv3/S7lJJhjsmrueHSpk2b5Pjjj5ejjjpKku95eu2114yHS//4xz/Sp5WS72hKPqovuWaXV7Kx3LZtWzrgKvorcUtebdq0Kfql1tX1kZf80t29e3enxW37n1w7cTcZjCuuGoHQ/a+5xnqKoY/DZ7uezV37f/svl8vyM4Cmuuq5PjReecaQK3t913uA7c8A9lcYZwS1qs8rdtnZufa//kqIRAABBNwEoh8uHX300el3GyVPBe38OuKII6Rbt24ye/ZsI8XmhkvJR/DNmTNH7r33Xtl9992thksXXHCB3H333fKLX/wiHUi5vsq0sWRD4pptP/HkxY+ryVldN5a2/U+uTbJifwyu9mYmEbG7hu5/E/N6Pib2eitibuvZ3LX/GS4VsaLr95rquZe1WXe9B9j+DKC9ztjiqFV9RrHLzs61//VXQiQCCCDgJhD9cCnUk0svvPCCnHjiienH2g0ZMiTNiumTS1OmTJEbb7xRkgHTqFGj3DL6/0WX6ZF4HqXOJOWZn4S8ZE4a7IS2/U+u/aQGV1z9CNQ+q23/53GNZVqTPg6fLczdzGO/B1AfbvURMppchdSurBV7//sSpVb1sthhpxcgEgEEYhGIfrjU3HcuDRo0SCZPnmyUz1pPLp177rmSDJiSj9hraGhIz/fmm2/K8OHDZcSIEXL66aenT1C1bt16h7WmT58u06ZNkzPPPFMuuugio+swOahMG0s2JCYZDX8MeQlvntWKtv1PrrOS3/E8uOLqR6D2WW37P49rLNOa9HH4bGHuZh77PYD6cKuPkNHkKqR2Za3Y+9+XKLWql8UOO70AkQggEItA9MOlqVOnyvXXXy9Lly6VLl26NOZt+fLlMmzYMLn00kvTAZDJq9ZwKfm+pO1/sTZ1ruuuu06OPPLIxj9OBlE//elP5ZRTTpHLL7/c5BKMjynTxpINiXFagx5IXoJyZ7qYbf+T60z5G0+GK65+BGqf1bb/87jGMq1JH4fPFuZu5rHfA6gPt/oIGU2uQmpX1oq9/32JUqt6Weyw0wsQiQACsQhEP1x69tlnJXk6KXmCacKECY15Gz9+vCxevFiWLFkinTt3luSzYtesWSPt27eXDh06VM1vreHS448/LuvXr98h7p133kmHV/3795cBAwZIr1690qeXktctt9wiP/rRjyQZSiUfi7f9aaesCqtMG0s2JFllPdvzkJdsPUOezbb/ybWf7OCKqx8BhkshXenjkNqVtTB3M7fdA7itFj6a+ghvrl2RXGnl9HGx979epnYktaqXxQ47vQCRCCAQi0D0w6UkURMnTpT58+fL0KFDpUePHvLoo4/KokWLZOzYsTJu3Lg0l0888YSMHDlyh/+W/Pdly5al/ySvuXPnyrp169KPsEteyZNQAwcObLIWmvrOpWSolaydDJqS71lq0aLFDufYf//9xfXL/Mq0sWRDUszbCXkpZl5Mrsq2/8m1iar9Mbjam5lE4Fpbybb/Tczr+RjqLXz2MXczj/0eQH241UfIaHIVUruyVuz970uUWtXLYoedXoBIBBCIRaAuhksffvihzJw5Mx0wrV27Vrp27Zp+FF7yJNL2J4aaGi4l34eUfC9StVefPn1k9uzZ1sOlWudMTvbxoZe20Mq0sWRDos2y3zjy4tfX59lt+59c+8kGrrj6EWC4FNKVPg6pXVkLczdz2z2A22rho6mP8ObaFcmVVk4fF3v/62VqR1KrelnssNMLEIkAArEI1MVwKZZk2byPMm0s2ZDYZDbcseQlnHXWK9n2P7nOOgP8gtSPKK4mrrb9b3LOej6G+2P47GPuZh77PYD6cKuPkNHkKqR2Za3Y+9+XKLWql8UOO70AkQggEIsAw6VYMrnT+yjTxpINSTGLkLwUMy8mV2Xb/+TaRNX+GFztzUwicK2tZNv/Jub1fAz1Fj77mLuZx34PoD7c6iNkNLkKqc1wyUWbWtXrYYedXoBIBBCIRYDhUiyZZLgUaSbze1tsFPOzd13Z9hdL5NpVvHo8rrj6EWC4FNKVPg6pXVkLczdz2z2A22rho6mP8ObaFcmVVk4fF3v/62VqR1KrelnssNMLEIkAArEIMFyKJZMMlyLNZH5vi41ifvauK9v+YEmuXcUZLvkRxFXjatv/mjXqKYb7Y/hsY+5mHvs9gPpwq4+Q0eQqpHZlrdj735cotaqXxQ47vQCRCCAQiwDDpVgyyXAp0kzm97bYKOZn77qy7Q+W5NpVnCGIH0FcNa62/a9Zo55iuD+Gzzbmbuax3wOoD7f6CBlNrkJqM1xy0aZW9XrYYacXIBIBBGIRYLgUSyYZLkWayfzeFhvF/OxdV7b9xRK5dhVnCOJHEFeNq23/a9aopxjuj+Gzjbmbeez3AOrDrT5CRpOrkNoMl1y0qVW9HnbY6QWIRACBWASCDpdGjhypcmtoaJBZs2apYus1qEw/WLIhKWaVkpdi5sXkqmz7n1ybqNofg6u9mUkErrWVbPvfxLyej6Hewmcfczfz2O8B1IdbfYSMJlchtRkuuWhTq3o97LDTCxCJAAKxCAQdLh111FFqt/vvv18dW4+BZfrBkg1JMSuUvBQzLyZXZdv/5NpE1f4YXO3NTCJwZbhkUidZHUO9ZSVpfh7Mza2qHWm7B3BbLXw09RHeXLsiudLK6eNi73+9TO1IalUvix12egEiEUAgFoGgw6VY0MrwPsq0sWRDUsyKIi/FzIvJVdn2P7k2UbU/Bld7M5MIXBkumdRJVsdQb1lJmp8Hc3MrhktuVkT7FaCX/frWY//7EqVW9bLYYacXIBIBBGIRYLgUSyZ3eh+2v1zOk4ENSZ76Ta9NXoqZF5Orsu1/cm2ian8MrvZmJhG4MlwyqZOsjqHespI0Pw/m5lb1+Mtl6sOtPkJGk6uQ2pW1bH8GCH+FxVyRWtXnBTvs9AJEIoBALAKFGS49//zz8uKLL8rGjRtl4MCBsfjm9j7KtLFkQ5JbmdRcmLwUMy8mV2Xb/+TaRNX+GFztzUwicGW4ZFInWR1DvWUlaX4ezM2tGC65WRHtV4Be9utbj/3vS5Ra1ctih51egEgEEIhFIPfh0l/+8he57LLLZPXq1Y2mK1euTP/v5M/OPPNMmTp1qhx99NGxmAd5H7a/XA5yUU0swoYkT/2m1yYvxcyLyVXZ9j+5NlG1PwZXezOTCFxrK9n2v4l5PR9DvYXPPuZu5rHfA6gPt/oIGU2uQmpX1oq9/32JUqt6Weyw0wsQiQACsQjkOlxasWKFjBgxQlq3bi0nn3xyOmB66KGHZPtwKUFOhko9evSQn//857GYB3kfZdpYsiEJUhLWi5AXa7LCBNj2P7n2kzpccfUjwHAppCt9HFK7shbmbua2ewC31cJHUx/hzbUrkiutnD4u9v7Xy9SOpFb1sthhpxcgEgEEYhHIdbj03e9+V5IB04IFC2TPPfeU6dOny4wZM3YYLl1wwQXy9NNPyx//+MdYzIO8jzJtLNmQBCkJ60XIizVZYQJs+59c+0kdrrj6EWC4FNKVPg6pzXApC23bPUAWa4Y8Bz0ZUtttLXLl5qeJjr3/NSYmMdSqiVL1Y7DDTi9AJAIIxCKQ63Cpd+/e0r9/f7niiitSz2rDpZ/+9Kdyyy23yF//+tdYzIO8jzJtLNmQBCkJ60XIizVZYQJs+59c+0kdrrj6EWC4FNKVPg6pzXApC23bPUAWa4Y8Bz0ZUtttLXLl5qeJjr3/NSYmMdSqiRLDJb0SdlnbcT4EECiWQK7DpV69esm3vvUtufjii1OVasOlf//3f5dFixal37/Ey1ygTBtLNnPmeQ15JHkJqZ3tWrb9T66z9d9+Nlxx9SNQ+6y2/Z/HNZZpTfo4fLYwdzOP/R5AfbjVR8hochVSu7JW7P3vS5Ra1ctih51egEgEEIhFINfh0qBBg+QTn/iEzJkzp+pwaevWrXLcccdJhw4d5LbbbovFPMj7KNPGkg1JkJKwXoS8WJMVJsC2/8m1n9ThiqsfAYZLIV3p45DalbUwdzO33QO4rRY+mvoIb65dkVxp5fRxsfe/XqZ2JLWql8UOO70AkQggEItArsOlm2++WSZPniyjR4+WcePGpd+3tP07lzZv3ixXXXVVOlT6j//4Dzn55JNjMQ/yPsq0sWRDEqQkrBchL9ZkhQmw7X9y7Sd1uOLqR4DhUkhX+jikNsOlLLRt9wBZrBnyHPRkSG23tciVm58mOvb+15iYxFCrJkrVj8EOO70AkQggEItArsOl5MmksWPHytKlS6Vz587SunVreeWVV+RLX/qSPPfcc/L222/LMcccI9OmTYvFO9j7KNPGkg1JsLKwWoi8WHEV6mDb/ifXftKHK65+BBguhXSlj0NqM1zKQtt2D5DFmiHPQU+G1HZbi1y5+WmiY+9/jYlJDLVqosRwSa+EXdZ2nA8BBIolkOtwKaHYtm2b3HrrrekTSqtXr07//+T1mc98Rr797W/LyJEjpaGhoVhqJbiaMm0s2cwVs6DISzHzYnJVtv1Prk1U7Y/B1d7MJALX2kq2/W9iXs/HUG/hs4+5m3ns9wDqw60+QkaTq5DalbVi739fotSqXhY77PQCRCKAQCwCuQ+XPg65adMmef/996Vt27bSrl27WIxzeR9l2liyIcmlRJpdlLw0S1TYA2z7n1z7SSWuuPoRYLgU0pU+DqldWQtzN3PbPYDbauGjqY/w5toVyZVWTh8Xe//rZWpHUqt6Weyw0wsQiQACsQgUargUC2oR3keZNpZsSIpQMbteA3kpZl5Mrsq2/8m1iar9Mbjam5lE4MpwyaROsjqGestK0vw8mJtbVTvSdg/gtlr4aOojvLl2RXKlldPHxd7/ehmGS9j5EtCfl3uk3o5IBBAolkDQ4dKaNWvU775Lly7q2HoMLNPGkr9Ui1mh5KWYeTG5Ktv+J9cmqvbH4GpvZhKBa20l2/43Ma/nY6i38NnH3M089nsA9eFWHyGjyVVI7cpasfe/L1FqVS+LHXZ6ASIRQCAWgaDDpe7du6u/P2nlypWxmAd5H2XaWLIhCVIS1ouQF2uywgTY9j+59pM6XHH1I8BwKaQrfRxSu7IW5m7mtnsAt9XCR1Mf4c21K5IrrZw+Lvb+18vUjqRW9bLYYacXIBIBBGIRCDpcmjZt2i7DpSeffFIeffRR+exnPysHHXSQ/PM//7O888476f/q5sUXX5Qvf/nL6X8fO3ZsLOZB3keZNpZsSIKUhPUi5MWarDABtv1Prv2kDldc/QgwXArpSh+H1Ga4lIW27R4gizVDnoOeDKnttha5cvPTRMfe/xoTkxhq1USp+jHYYacXIBIBBGIRCDpc2hntsccek1GjRskVV1whJ5100i6mCxYskEsvvVSuv/56Oeyww2IxD/I+yrSxZEMSpCSsFyEv1mSFCbDtf3LtJ3W44upHgOFSSFf6OKQ2w6UstG33AFmsGfIc9GRIbbe1yJWbnyY69v7XmJjEUKsmSgyX9ErYZW3H+RBAoFgCuQ6XTjnlFEm+S+nqq69uUuX888+XN954Q26//Xa13JYtW2TmzJkyb948eeutt6Rr164yYsQIGT58eLMf07dw4UJ54IEHZMWKFfLSSy9Jp06d5KGHHmr2WjZv3iwDBgxIY8455xxJ3sfOr+TcyeBs9erVsscee8ixxx4r48ePl7Zt2zZ7/uYOKNPGks1cc9nM58/JSz7uWaxq2//kOgv1Xc+BK65+BGqf1bb/87jGMq1JH4fPFuZu5rHfA6gPt/oIGU2uQmpX1oq9/32JUqt6Weyw0wsQiQACsQjkOlzq1auXnHbaaVUHL9uBk8HT7NmzJfn4PO3rkksukTlz5sjQoUOlZ8+e8sgjj8g999wj48aNa/bj9k499VR55pln5IADDkgHRS1atDAaLs2YMUNuuOEG2bhxY9Xh0u9+9zu58MIL0yeyjj/++PTcN998s/Tu3VtuuummZodezVmUaWPJhqS5bObz5+QlH/csVrXtf3KdhTrDJT+KuNq62va/7fnr7Xjuj+Ezjrmbeez3AOrDrT5CRpOrkNoMl1y0qVW9HnbY6QWIRACBWARyHS595Stfkb333ltuvfXWJj2//e1vyyuvvJJ+L5PmtXLlShk4cKCcccYZMmHChMZTJE8ILVmyJP0neRqpqVfy1FTy5y1btpRk0PTyyy83O1x69dVX5YQTTpAxY8bI1KlTdxkuJU81HXnkkelTW8kTWcm5k1ficPnll8v06dPlmGOO0bzdxpgy/WDJhsQp1d6CyYs3Wu8ntu1/cu0nJbji6keg9llt+z+PayzTmvRx+Gxh7mYe+z2A+nCrj5DR5CqkdmWt2Pvflyi1qpfFDju9AJEIIBCLQK7DpRHU2IcAACAASURBVOS7lm655ZZ0+JM8RZQMW7a/1qxZI9OmTZPke5eSj69Lnj7SvJInn5KPxFu6dOkO51++fLkMGzZMJk2alP7b5GU6XDr77LNlw4YNctVVV0m/fv12GS4lT06deeaZMmXKlPS9b38lQ6dDDz1Uvva1r8l//ud/mlxSk8eUaWPJhsQp1d6CyYs3Wu8ntu1/cu0nJbji6keg9llt+z+PayzTmvRx+Gxh7mYe+z2A+nCrj5DR5CqkdmWt2Pvflyi1qpfFDju9AJEIIBCLQK7DpWQA893vfjfdBCVP73z605+WDh06yLvvvpt+N9LWrVvlkEMOkf/6r/+ST37ykyrz5ImlVatWpR+F9/FXMsg58MADZfDgwXLllVcandtkuLR48WI577zz5M4770y/O6nacOm6666Tn//857Jo0SL57Gc/u8PayZNab7/9ttx3331G19TUQYnptm3bMvn+JqcLMQjetGlTelSbNm0MjuaQUALkJZT0rut0797daXHb/ifXTtxNBuOKq0YgdP9rrrGeYujj8NmuZ3PX/k+yZbsHCJ9htxXruT7c5MJHkyt7c9d7QOz9by9qFkGtmjlVOwq77Oxc+19/JUQigAACbgK5DpeSS//oo49k7ty58vvf/z4dAq1fv17atWsn++23n5x44onyzW9+M/2eI+0r+Xi6Vq1ayfz583c5Rd++fdPvUkq+G8nk1dxwKfmLNfn+pKOOOip90uq1116rOlz60Y9+lD6xlTw9lbzXj7/+7d/+LX3KasWKFSaX1OQxZdpYsiFxSrW3YPLijbbZE7tuLG37n1w3mxLVAbiq2JoNit01dP83C17nB8Reb0VMbz2bu/Z/kk/bPUARa6DWNdVzfZCrsgnYX6/rPSD2/rcXNYvgvmLmVO0o7LKzc+1//ZUQiQACCLgJ5D5ccrv85qOPPvpo6dixY/rdRju/jjjiCOnWrZvMnj27+ROJNPudS8lH8M2ZM0fuvfde2X333ZscLk2cOFHmzZsnTz/9dDr4+vjroosukrvuukuSx4sbGhqMrqvaQWV6JJ5HqdVp9hpIXrzyej25bf+Taz/pwBVXPwK1z2rb/3lcY5nWpI/DZwtzN/PY7wHUh1t9hIwmVyG1K2vF3v++RKlVvSx22OkFiEQAgVgECjdcSj7KzWWosnNiQj259MILL6RPWiXf4TRkyJD0MvJ+cim5hoMPPrjwtcqGpJgpIi/FzIvJVdn+YEmuTVTtj8HV3swkAleGSyZ1ktUx1FtWkubnwdzcqtqRtnsAt9XCR1Mf4c21K5IrrZw+Lvb+18vUjqRW9bLYYacXIBIBBGIRyH24lAyT7rjjjsandf7xj39I69atJXkkdODAgTJ06FCnYVNz37k0aNAgmTx5slE+a30s3rnnnivJgCn5iL3tw7E333xThg8fLiNGjJDTTz89fYIqeW/NfedS8n1TyXc3ubzKtLFkQ+KSaX+x5MWfre8z2/Y/ufaTEVxx9SPAcCmkK30cUruyFuZu5rZ7ALfVwkdTH+HNtSuSK62cPi72/tfLMFzCzpeA/rzcI/V2RCKAQLEEch0ubd68Wc4++2x5/PHH04HMnnvumQ5g3n77bUkGM8n3MR122GEyc+bMXT4+zpRx6tSpcv3116ffY9SlS5fGsOT7joYNGyaXXnppOgAyedUaLp100kmNPww3da5kqHTkkUfKww8/LGeddZZMmTIlHaBtfyUehx56qHz1q1+Va665xuSSmjymTBtL/lJ1SrW3YPLijdb7iW37n1z7SQmuuPoRqH1W2/7P4xrLtCZ9HD5bmLuZx34PoD7c6iNkNLkKqV1ZK/b+9yVKreplscNOL0AkAgjEIpDrcGn69OmS/HP88cfL+eefL3vttVej65o1ayQZDC1cuFDGjh0rY8aMUZk/++yzkjydlDzBNGHChMZzjB8/Pn06aMmSJdK5c2dJvogwWbN9+/bSoUOHqmvVGi4lA7L169fvEPfOO++kw6v+/fvLgAEDpFevXunwLBkiJd/31LVr1/SprRYtWqRxt956q1x++eVy7bXXyje+8Q3V+90eVKaNJRsSp1R7CyYv3mi9n9i2/8m1n5Tgiqsfgdpnte3/PK6xTGvSx+Gzhbmbeez3AOrDrT5CRpOrkNqVtWLvf1+i1KpeFjvs9AJEIoBALAK5DpeSAcoee+whv/3tb5v0/Na3viXvvfee3HvvvWrziRMnyvz589OP2OvRo4c8+uijsmjRonRoNW7cuPS8TzzxhIwcOXKH/5b892XLlqX/JK+5c+fKunXr5Mwzz0z//+RJqI8/ebTzBTb1nUvJcQsWLEiHXX379pXjjjtOXn75ZZk1a5YcdNBBcvPNNzt9FGDZNpZsSNSl7TWQvHjl9Xpy2x8sybWfdOCKqx+B2me17f88rrFMa9LH4bOFuZt57PcA6sOtPkJGk6uQ2pW1Yu9/X6LUql4WO+z0AkQigEAsArkOl5JBT/JdRN/73vea9Lz66qvlpptukqefflpt/uGHH6YfrZcMmNauXZs+MZR8FF7yJNL270dqarg0bdq09Omqaq8+ffrI7Nmzm7yuWsOlJOgPf/hD+pF9yXc1JUO25Amn5Imqdu3aqd/r9sAybSzZkDin28sJyIsX1iAnte1/cu0nLbji6keg9llt+z+PayzTmvRx+Gxh7mYe+z2A+nCrj5DR5CqkdmWt2Pvflyi1qpfFDju9AJEIIBCLQK7DpS996Uvy5S9/WX7605826XnhhRemTxr96U9/isU8yPso08aSDUmQkrBehLxYkxUmwLb/ybWf1OGKqx8BhkshXenjkNqVtTB3M7fdA7itFj6a+ghvrl2RXGnl9HGx979epnYktaqXxQ47vQCRCCAQi0Cuw6Xke5buu+8+mTFjhnzta1/bxfTBBx9Mv2vp61//uiRPMPEyFyjTxpINiXleQx5JXkJqZ7uWbf+T62z9t58NV1z9CNQ+q23/53GNZVqTPg6fLczdzGO/B1AfbvURMppchdSurBV7//sSpVb1sthhpxcgEgEEYhHIdbj00ksvycknnywbNmyQ3r17y8EHHywdO3aUt99+O90Y/eUvf5FPfepT6XcyfeYzn4nFPMj7KNPGkg1JkJKwXoS8WJMVJsC2/8m1n9ThiqsfAYZLIV3p45DalbUwdzO33QO4rRY+mvoIb65dkVxp5fRxsfe/XqZ2JLWql8UOO70AkQggEItArsOlBPG5556Tyy67TJ588sldTA855BCZNGmS7LfffrF4B3sfZdpYsiEJVhZWC5EXK65CHWzb/+TaT/pwxdWPQO2z2vZ/HtdYpjXp4/DZwtzNPPZ7APXhVh8ho8lVSO3KWrH3vy9RalUvix12egEiEUAgFoHch0vbIV9//fV00LR+/Xpp166d7L///tK1a9dYnIO/jzJtLNmQBC8PowXJixFTIQ+y7X9y7SeNuOLqR4DhUkhX+jikdmUtzN3MbfcAbquFj6Y+wptrVyRXWjl9XOz9r5epHUmt6mWxw04vQCQCCMQiUJjhUiygRXkfZdpYsiEpStXseB3kpZh5Mbkq2/4n1yaq9sfgam9mEoErwyWTOsnqGOotK0nz82BublXtSNs9gNtq4aOpj/Dm2hXJlVZOHxd7/+tlGC5h50tAf17ukXo7IhFAoFgCDJeKlY/MrqZMG0v+Us0s7ZmeiLxkyhn0ZLb9T679pAdXXP0IMFwK6Uofh9SurIW5m7ntHsBttfDR1Ed4c+2K5Eorp4+Lvf/1MgyXsPMloD8v90i9HZEIIFAsgdyHS4899pjMmjVLVq1aJWvXrpWtW7fuItTQ0CDPPvtsseQKfjVl2ljyl2oxi4m8FDMvJldl2//k2kTV/hhc7c1MInCtrWTb/ybm9XwM9RY++5i7mcd+D6A+3OojZDS5CqldWSv2/vclSq3qZbHDTi9AJAIIxCKQ63Dp9ttvl8svv1y2bdsm++yzj3Ts2FFatGhR1Xb27NmxmAd5H2XaWLIhCVIS1ouQF2uywgTY9j+59pM6XHH1I8BwKaQrfRxSu7IW5m7mtnsAt9XCR1Mf4c21K5IrrZw+Lvb+18vUjqRW9bLYYacXIBIBBGIRyHW4dNRRR8kHH3wgv/rVr6R79+6xmBbifZRpY8mGpBAls8tFkJdi5sXkqmz7n1ybqNofg6u9mUkErgyXTOokq2Oot6wkzc+DublVtSNt9wBuq4WPpj7Cm2tXJFdaOX1c7P2vl2G4hJ0vAf15uUfq7YhEAIFiCeQ6XDrwwANl6NCh8sMf/rBYKhFcTZk2lvylWsyCIy/FzIvJVdn2P7k2UbU/Bld7M5MIXBkumdRJVsdQb1lJmp8Hc3MrhktuVkT7FaCX/frWY//7EqVW9bLYYacXIBIBBGIRyHW4NHjwYNlvv/3kqquuisWzMO/D9pfLeV44G5I89Ztem7wUMy8mV2Xb/+TaRNX+GFztzUwicGW4ZFInWR1DvWUlaX4ezM2t6vGXy9SHW32EjCZXIbUra9n+DBD+Cou5IrWqzwt22OkFiEQAgVgEch0u/fGPf5Qf/OAHcscdd8jnP//5WEwL8T7KtLFkQ1KIktnlIshLMfNiclW2/U+uTVTtj8HV3swkAleGSyZ1ktUx1FtWkubnwdzciuGSmxXRfgXoZb++9dj/vkSpVb0sdtjpBYhEAIFYBHIdLiWICxculMmTJ0vy/Uv777+/tG3btqrtwIEDYzEP8j5sf7kc5KKaWIQNSZ76Ta9NXoqZF5Orsu1/cm2ian8MrvZmJhG41lay7X8T83o+hnoLn33M3cxjvwdQH271ETKaXIXUrqwVe//7EqVW9bLYYacXIBIBBGIRyHW4tHHjxvTJpfvuu0+2bduWmjY0NOxgm/z35L+tXLkyFvMg76NMG0s2JEFKwnoR8mJNVpgA2/4n135ShyuufgQYLoV0pY9DalfWwtzN3HYP4LZa+GjqI7y5dkVypZXTx8Xe/3qZ2pHUql4WO+z0AkQigEAsArkOl5LB0oIFC+SAAw6QY445Rjp27CgtW7asajto0KBYzIO8jzJtLNmQBCkJ60XIizVZYQJs+59c+0kdrrj6EWC4FNKVPg6pzXApC23bPUAWa4Y8Bz0ZUtttLXLl5qeJjr3/NSYmMdSqiVL1Y7DDTi9AJAIIxCKQ63Dp0EMPlX333VduvfVWadGiRSymhXgfZdpYsiEpRMnschHkpZh5Mbkq2/4n1yaq9sfgam9mEoErwyWTOsnqGOotK0nz82BublXtSNs9gNtq4aOpj/Dm2hXJlVZOHxd7/+tlakdSq3pZ7LDTCxCJAAKxCOQ6XOrTp49885vflAkTJsTiWZj3UaaNJRuSwpTNDhdCXoqZF5Orsu1/cm2ian8MrvZmJhG4MlwyqZOsjqHespI0Pw/m5lYMl9ysiPYrQC/79a3H/vclSq3qZbHDTi9AJAIIxCKQ63DpvPPOk3Xr1smsWbNi8SzM+7D95XKeF86GJE/9ptcmL8XMi8lV2fY/uTZRtT8GV3szkwhcGS6Z1ElWx1BvWUmanwdzc6t6/OUy9eFWHyGjyVVI7cpatj8DhL/CYq5Irerzgh12egEiEUAgFoFch0tvvPGGDB8+XE488UQZPXq0tGrVKhbX3N9HmTaWbEhyL5eqF0BeipkXk6uy7X9ybaJqfwyu9mYmEbgyXDKpk6yOod6ykjQ/D+bmVgyX3KyI9itAL/v1rcf+9yVKreplscNOL0AkAgjEIpDrcGnkyJHy/vvvy3PPPSdt27aVffbZJ/33zq+GhgaebrKsONtfLluePtPD2ZBkypnZychLZpTBT2Tb/+TaT4pwxdWPAMOlkK70cUjtylqYu5nb7gHcVgsfTX2EN9euSK60cvq42PtfL1M7klrVy2KHnV6ASAQQiEUg1+FS9+7djRyT4dLKlSuNjuWgikCZNpZsSIpZteSlmHkxuSrb/ifXJqr2x+Bqb2YSgWttJdv+NzGv52Oot/DZx9zNPPZ7APXhVh8ho8lVSO3y/Q4gvE7TK1Kr+mxgh51egEgEEIhFINfhUiyIRXwfZfrBkg1JESuI/+VwMbNidlW2/U8PmrnaHoWrrZjZ8bgyXDKrlGyOot6ycbQ5C+Y2Wrsea7sHcFstfDT1Ed5cuyK50srp42Lvf71M7UhqVS+LHXZ6ASIRQCAWgVIOl5K/wJJ/Bg4cGEseMn8fZdpYsiHJPP2ZnJC8ZMKYy0ls+59c+0kTrrj6EWC4FNKVPg6pXVkLczdz2z2A22rho6mP8ObaFcmVVk4fF3v/62UYLmHnS0B/Xu6RejsiEUCgWAKlHC5Nnz5dZsyYwUfl1ailMm0s+Uu1WDeF7VdDXoqZF5Orsu1/cm2ian8MrvZmJhG4MlwyqZOsjqHespI0Pw/m5lbVjrTdA7itFj6a+ghvrl2RXGnl9HGx979ehuESdr4E9OflHqm3IxIBBIolwHCpWPnI7GrKtLHkL9XM0p7pichLppxBT2bb/+TaT3pwxdWPAMOlkK70cUjtylqYu5nb7gHcVgsfTX2EN9euSK60cvq42PtfL8NwCTtfAvrzco/U2xGJAALFEqiL4dKWLVtk5syZMm/ePHnrrbeka9euMmLECBk+fLg0NDTUzMjChQvlgQcekBUrVshLL70knTp1koceeqhqzJQpU2TZsmXy6quvysaNG6Vz587Sp08fGT16tOy11147xCTXdPvtt8ucOXPklVdekdatW8vnP/95Oeuss+SrX/2qc5WUaWPJX6rO6fZyAvLihTXISW37n1z7SQuuuPoRqH1W2/7P4xrLtCZ9HD5bmLuZx34PoD7c6iNkNLkKqV1ZK/b+9yVKreplscNOL0AkAgjEIlAXw6VLLrkkHeIMHTpUevbsKY888ojcc889Mm7cOBk7dmzNXJ566qnyzDPPyAEHHJAOl1q0aNHkcCkZVu2///6y9957S9u2bdMh09y5cyUZJCWDrW7dujWutf2aTjjhBPniF78oGzZsSI994YUX5Nprr5VvfOMbTjVWpo0lGxKnVHsLJi/eaL2f2Lb/ybWflOCKqx8BhkshXenjkNqVtTB3M7fdA7itFj6a+ghvrl2RXGnl9HGx979epnYktaqXxQ47vQCRCCAQi0D0w6WVK1fKwIED5YwzzpAJEyY05m38+PGyZMmS9J/kaaSmXm+88Ub65y1btpRk0PTyyy83OVyqdo6nn35aTj75ZBk1apRccMEF6SHr169Pn2jq16+fTJs2rTHs3XffTZ9a+vKXv5w+aeXyKtPGkg2JS6b9xZIXf7a+z2zb/+TaT0ZwxdWPQO2z2vZ/HtdYpjXp4/DZwtzNPPZ7APXhVh8ho8lVSO3KWrH3vy9RalUvix12egEiEUAgFoHoh0tXX311OqhZunSpdOnSpTFvy5cvl2HDhsmkSZPSf5u8NMOlZGDUt29fOeWUU+Tyyy9Pl0k+mu8rX/lK+rF8l156aePSH330kfTu3VsOP/xwueaaa0wuqcljyrSxZEPilGpvweTFG633E9v2P7n2kxJccfUjUPustv2fxzWWaU36OHy2MHczj/0eQH241UfIaHIVUruyVuz970uUWtXLYoedXoBIBBCIRSD64VLyxNKqVavSj8L7+Gvz5s1y4IEHyuDBg+XKK680yqfJcCkZEL333nuydetWef3112XGjBnpk07Jv48++ujGdY4//vj0zy+77LL0KabkY/F+9atfyb333iu//vWv02tzeSUby23btqUfz1f016ZNm9JLbNOmTdEvta6uj7zkl+7u3bs7LW7b/+TaibvJYFxx1QiE7n/NNdZTDH0cPtv1bO7a/9t/uVyWnwE01VXP9aHxyjOGXNnru94DbH8GsL/COCOoVX1escvOzrX/9VdCJAIIIOAmEP1wKflOo1atWsn8+fN3kUqeKEq+S+mGG24wUjQZLr322mvpx91tf7Vv317OOecc+c53vrPDGqtXr5bvf//78uyzzzb+909/+tPpEMp1sFS2HyzZkBiVX/CDyEtw8sYFXTeWtj9Ykms/ucYVV41A6P7XXGM9xdDH4bNdz+au/V+2nwE01VXP9aHxyjOGXNnru94DbH8GsL/COCOoVX1escvOzrX/9VdCJAIIIOAmEP1wKXlaqGPHjnL77bfvInXEEUdIt27dZPbs2UaKJsOlDz74QJKP3EuejEoGSHfffXc6bBo9erS0aNGicZ2//e1vknxkX7t27eTQQw9Nv4fpN7/5jSTDqWTY1bNnT6NrauqgMj0Sz6PUTqn2FkxevNF6P7Ft/5NrPynBFVc/ArXPatv/eVxjmdakj8NnC3M389jvAdSHW32EjCZXIbUra8Xe/75EqVW9LHbY6QWIRACBWARKOVxK/gJbuXKlDBo0qNk8hH5yaecLWrNmjSTXkAymzj///PSPk4/AS/5b8s8FF1zQGJL8rz6Sj8v71Kc+JXfddVez763WAWXaWLIhsUv1pg3/kM2bNkurNq2kTdvWdsEWR5MXC6yCHWrb/+TaTwJ9uYa6B/hRcT+rL1f3KyvGGWz7vxhXXdyrKFq91UP/F828uNVZ/cpivwfUc32Urf/rOVd53Tdi739frtSqXhY77PQCRCKAQCwCpRwu2eA3951LyYBq8uTJRqc0eXKp2onOOussee655+Thhx9O//jOO++UH/zgBzJnzpxdnlCaNGlS+pTVX/7yl3TIpH2VaWPJhsQsy39/d728+dJaueMnd8kbq9+UPT/3L/Kti06Sf/lMJ/lUh3ZmJ7E4irxYYBXsUNv+J9d+Epi1a+h7gB8V97Nm7ep+RcU6g23/F+vqi3c1Ram3eur/opgXrxrNrij2e0A91kdZ+78ec2XWpf6Oir3/fclRq3pZ7LDTCxCJAAKxCAQdLiWfIdrQ0GBtl8R8/LuJbE4wdepUuf7662Xp0qXSpUuXxtDko+uGDRsml156qQwfPtzolNrhUhK3YsUKeeqpp9J1Zs6cmX4kXjJEOuigg3ZY+5JLLkmHTo899ph06NDB6LqqHVSmjSUbkubTnPxQOefq38ttk3f97rBvTxwsQ743IPMBE3lpPi9FPcK2/8m1n0xm6ZrHPcCPivtZs3R1v5rincG2/4v3Dop1RUWot3rr/yKYF6sK7a4m9ntAvdVHmfu/3nJl16l+jo69//2oiVCrelnssNMLEIkAArEIBB0uJU/raIZLCfaPf/xjlXkylEqeTkqeYJowYULjOcaPHy+LFy+WJUuWSOfOnSX5SLrkI+zat2/f5FCn1nDp/ffflzZt2shuu+22w3Umf9kOGTJEevXq1fjdTvfdd5+MHTs2HW4lTyptf61bty79WLxWrVrJ/fffr3q/24PKtLFkQ9J8qp//7xdkdO//v353jvjFX6bI5w/+bPMnsjiCvFhgFexQ2/4n134SmKVrHvcAPyruZ83S1f1qincG2/4v3jso1hUVod7qrf+LYF6sKrS7mtjvAfVWH2Xu/3rLlV2n+jk69v73o8ZwycWVPtfrYae3IxIBBIolEHS4lNdbnzhxosyfP1+GDh0qPXr0kEcffVQWLVqUDnjGjRuXXtYTTzwhI0eO3OG/Jf992bJl6T/Ja+7cuZIMgM4888z0/0+ehBo4cGD6fyeDqssuu0z69+8v++yzj7Rs2VKef/55WbBgQfrns2bNavwIvA8//DAdOCXfG3XMMcdI37590+9huuOOO+S1116Tn/zkJ3LSSSc5cZVpY8lfqrVTnXy++tQzfykP/vZPTR54xLe+JBfccK60zvA7mMiLUwvmGmzb/+TaT7qycs3rHuBHxf2sWbm6X0kxz2Db/8V8F8W5qrzrrR77P2/z4lSf7kpivwfUU32Uvf/rKVe6bs0+Kvb+z16sckZqVS+LHXZ6ASIRQCAWgboYLiXDnOSj6JIB09q1a6Vr167pR+ElTyJtf5KqqeHStGnTZPr06VXz3adPn8ankV555RX55S9/KcmGLlkjWbNTp05y6KGHyqhRo2Tffffd4Rzr16+XG2+8Uf74xz/K66+/nv7ZF77whfQJq379+jnXV5k2lmxIaqd73dvvy8Rjr5RVy19o8sD9en9OJi+cKHt03N25drafgLxkRhn8RLb9T679pCgr17zuAX5U3M+alav7lRTzDLb9X8x3UZyryrve6rH/8zYvTvXpriT2e0A91UfZ+7+ecqXr1uyjYu//7MUqZ6RW9bLYYacXIBIBBGIRqIvhUizJsnkfZdpYsiGpndm8/leL5MWm44p1rG3/k2s/+cvKNa97gB8V97Nm5ep+JcU8g23/F/NdFOeq8q63euz/vM2LU326K4n9HlBP9VH2/q+nXOm6Nfuo2Ps/ezGGS66m9LleEDu9HZEIIFAsgaDDJe0TOcnTRcnHzvEyFyjTxpK/VJvPax6ft05ems9LUY+w7X9y7SeTWbrmcQ/wo+J+1ixd3a+meGew7f/ivYNiXVER6q3e+r8I5sWqQrurif0eUG/1Ueb+r7dc2XWqn6Nj738/ajy55OJKn+v1sNPbEYkAAsUSCDpcSj6GTvuaPXu2NrQu48q0seQv1eZL9O/vrpc5V/9ebps8f5eDh/1wsJx8/gD5VId2zZ/I4gjyYoFVsENt+59c+0lglq553AP8qLifNUtX96sp3hls+79476BYV1SEequ3/i+CebGq0O5qYr8H1Ft9lLn/6y1Xdp3q5+jY+9+PGsMlF1f6XK+Hnd6OSAQQKJZA0OFSsd563FdTpo0lf6ma1WLyw+WbL62V3/70Llmz+m/S5XOdZeiFJ8m/fKZT5oOl5IrIi1leiniUbf+Taz9ZzNo19D3Aj4r7WbN2db+iYp3Btv+LdfXFu5qi1Fs99X9RzItXjWZXFPs9oB7ro6z9X4+5MutSf0fF3v++5KhVvSx22OkFiEQAgVgEGC7Fksmd3keZNpZsSOyK8B8b/iEfbNos/0+bVtK6bWu7YIujyYsFVsEOte1/cu0n/2x1AgAAIABJREFUgb5cQ90D/Ki4n9WXq/uVFeMMtv1fjKsu7lUUrd7qof+LZl7c6qx+ZbHfA+q5PsrW//Wcq7zuG7H3vy9XalUvix12egEiEUAgFgGGS7FkkuFSpJnM722xUczP3nVl2x8sybWrePV4XHH1I1D7rLb9n8c1lmlN+jh8tjB3M4/9HkB9uNVHyGhyFVK7slbs/e9LlFrVy2KHnV6ASAQQiEUg1+HSyJEjjRwbGhpk1qxZRsdyUPk2lmxIilm15KWYeTG5KtsfLMm1iar9Mbjam5lE4MpwyaROsjqGestK0vw8mJtbVTvSdg/gtlr4aOojvLl2RXKlldPHxd7/epnakdSqXhY77PQCRCKAQCwCuQ6XjjrqqKqOGzZskHXr1kkyVOrYsaPstttucv/998diHuR9lGljyYYkSElYL0JerMkKE2Db/+TaT+pwxdWPAMOlkK70cUjtylqYu5nb7gHcVgsfTX2EN9euSK60cvq42PtfL8NwCTtfAvrzco/U2xGJAALFEsh1uFSL4vXXX5cpU6bI3/72N7nxxhulbdu2xZIr+NWUaWPJX6rFLCbyUsy8mFyVbf+TaxNV+2NwtTczicC1tpJt/5uY1/Mx1Fv47GPuZh77PYD6cKuPkNHkKqR2Za3Y+9+XKLWql8UOO70AkQggEItAYYdLCfCWLVtk0KBB0rt3b5k0aVIs5kHeR5k2lmxIgpSE9SLkxZqsMAG2/U+u/aQOV1z9CDBcCulKH4fUrqyFuZu57R7AbbXw0dRHeHPtiuRKK6ePi73/9TK1I6lVvSx22OkFiEQAgVgECj1cSpCvvPJK+cMf/iB/+tOfYjEP8j7KtLFkQxKkJKwXIS/WZIUJsO1/cu0ndbji6keA4VJIV/o4pDbDpSy0bfcAWawZ8hz0ZEhtt7XIlZufJjr2/teYmMRQqyZK1Y/BDju9AJEIIBCLQOGHSz/4wQ9k0aJF8tRTT8ViHuR9lGljyYYkSElYL0JerMkKE2Db/+TaT+pwxdWPAMOlkK70cUhthktZaNvuAbJYM+Q56MmQ2m5rkSs3P0107P2vMTGJoVZNlBgu6ZWwy9qO8yGAQLEECjtc2rZtm9x9991y8cUXy4EHHii33HJLseQKfjVl2liymStmMZGXYubF5Kps+59cm6jaH4OrvZlJBK4Ml0zqJKtjqLesJM3Pg7m5VbUjbfcAbquFj6Y+wptrVyRXWjl9XOz9r5epHUmt6mWxw04vQCQCCMQikOtwqV+/flUdt27dKu+88076nUtt27aVm266SXr06BGLeZD3UaaNJRuSICVhvQh5sSYrTIBt/5NrP6nDFVc/AgyXQrrSxyG1K2th7mZuuwdwWy18NPUR3ly7IrnSyunjYu9/vQzDJex8CejPyz1Sb0ckAggUSyDX4dKpp55aVaNFixay++67ywEHHCCDBw+WTp06FUutBFdTpo0lf6kWs6DISzHzYnJVtv1Prk1U7Y/B1d7MJALX2kq2/W9iXs/HUG/hs4+5m3ns9wDqw60+QkaTq5DalbVi739fotSqXhY77PQCRCKAQCwCuQ6XYkEs4vso08aSDUkRK4j/5XAxs2J2Vbb9Tw+audoehautmNnxuDJcMquUbI6i3rJxtDkL5jZaux5ruwdwWy18NPUR3ly7IrnSyunjYu9/vUztSGpVL4sddnoBIhFAIBYBhkuxZHKn91GmjSUbkmIWIXkpZl5Mrsq2/8m1iar9Mbjam5lE4MpwyaROsjqGestK0vw8mJtbVTvSdg/gtlr4aOojvLl2RXKlldPHxd7/ehmGS9j5EtCfl3uk3o5IBBAolgDDpWLlI7OrKdPGkr9UM0t7piciL5lyBj2Zbf+Taz/pwRVXPwIMl0K60schtStrYe5mbrsHcFstfDT1Ed5cuyK50srp42Lvf70MwyXsfAnoz8s9Um9HJAIIFEsg9+HS0qVLZfbs2fLss8/K3//+d/noo492EWpoaEj/nJe5QJk2lvylap7XkEeSl5Da2a5l2//kOlv/7WfDFVc/ArXPatv/eVxjmdakj8NnC3M389jvAdSHW32EjCZXIbUra8Xe/75EqVW9LHbY6QWIRACBWARyHS7NmzdPLrnkEmnZsqUccsgh0qlTJ/nEJz5R1fbHP/5xLOZB3keZNpZsSIKUhPUi5MWarDABtv1Prv2kDldc/QgwXArpSh+H1K6shbmbue0ewG218NHUR3hz7YrkSiunj4u9//UytSOpVb0sdtjpBYhEAIFYBHIdLn3jG9+QjRs3ym233SZ77bVXLKaFeB9l2liyISlEyexyEeSlmHkxuSrb/ifXJqr2x+Bqb2YSgSvDJZM6yeoY6i0rSfPzYG5uVe1I2z2A22rho6mP8ObaFcmVVk4fF3v/62UYLmHnS0B/Xu6RejsiEUCgWAK5Dpd69Oghp5xyivzwhz8slkoEV1OmjSV/qRaz4MhLMfNiclW2/U+uTVTtj8HV3swkAleGSyZ1ktUx1FtWkubnwdzciuGSmxXRfgXoZb++9dj/vkSpVb0sdtjpBYhEAIFYBHIdLh133HHSq1cvmTx5ciyehXkftr9czvPC2ZDkqd/02uSlmHkxuSrb/ifXJqr2x+Bqb2YSgSvDJZM6yeoY6i0rSfPzYG5uVY+/XKY+3OojZDS5CqldWcv2Z4DwV1jMFalVfV6ww04vQCQCCMQikOtw6be//a1MnTpVFixYIHvuuWcspoV4H2XaWLIhKUTJ7HIR5KWYeTG5Ktv+J9cmqvbH4GpvZhKBK8MlkzrJ6hjqLStJ8/Ngbm7FcMnNimi/AvSyX9967H9fotSqXhY77PQCRCKAQCwCuQ6Xli1bJjfddJP89a9/lZEjR8r+++8v7dq1q2r7xS9+MRbzIO/D9pfLQS6qiUXYkOSp3/Ta5KWYeTG5Ktv+J9cmqvbH4GpvZhKBa20l2/43Ma/nY6i38NnH3M089nsA9eFWHyGjyVVI7cpasfe/L1FqVS+LHXZ6ASIRQCAWgVyHS927d5eGhgbZtm1b+u9ar5UrV6rMt2zZIjNnzpR58+bJW2+9JV27dpURI0bI8OHDm11z4cKF8sADD8iKFSvkpZdekk6dOslDDz1U9TqmTJkiybDs1VdflY0bN0rnzp2lT58+Mnr0aNlrr712idm6davcdtttMnfuXHnxxRelVatW8rnPfU7GjBkjhx9+uOq9fjyoTBtLNiTO6fZyAvLihTXISW37n1z7SQuuuPoRqH1W2/7P4xrLtCZ9HD5bmLuZx34PoD7c6iNkNLkKqV1ZK/b+9yVKreplscNOL0AkAgjEIpDrcGnatGnNDni2Q48dO1Zlfskll8icOXNk6NCh0rNnT3nkkUfknnvukXHjxklz5zz11FPlmWeekQMOOCAdLrVo0aLJ4VIyrEqevNp7772lbdu26ZApGRwlw61ksNWtW7fG6//oo4/kvPPOkwcffFAGDRokPXr0kE2bNsn//d//pf/3kCFDVO+V4ZIzGyf4mAAbxfKWg+0PluTaT65xxdWPAMOlkK70cUjtylqYu5nb7gHcVgsfTX2EN9euSK60cvq42PtfL1M7klrVy2KHnV6ASAQQiEUg1+GSb8TkaaeBAwfKGWecIRMmTGhcbvz48bJkyZL0n+RppKZeb7zxRvrnLVu2lGTQ9PLLLzc5XKp2jqefflpOPvlkGTVqlFxwwQWNh9x8882SPOk0a9Ys6d27txeGMm0s2ZB4KQHnk5IXZ8LcTmDb/+TaT6pwxdWPQO2z2vZ/HtdYpjXp4/DZwtzNPPZ7APXhVh8ho8lVSO3KWrH3vy9RalUvix12egEiEUAgFoFSDpeSoUwyoEmGQ7VeV199dfqReEuXLpUuXbo0Hrp8+XIZNmyYTJo0Kf23yUszXHr33Xelb9++csopp8jll1+eLpM8tdSvX7/0CaVrr702/f+Tp5aSp52yfJVpY8mGJMvMZ3cu8pKdZegz2fY/ufaTIVxx9SNQ+6y2/Z/HNZZpTfo4fLYwdzOP/R5AfbjVR8hochVSu7JW7P3vS5Ra1ctih51egEgEEIhFoJTDpenTp8uMGTOkue9hSp5YWrVqVfpReB9/bd68WQ488EAZPHiwXHnllUa5NBkuJYOi9957T5LvU3r99dfTa0y+oyn599FHH52uk3z03fHHHy/nn3++rF27VubPn58Ol/bcc08555xz0kFUFq8ybSzZkGSR8ezPQV6yNw11Rtv+J9d+MoMrrn4Eap/Vtv/zuMYyrUkfh88W5m7msd8DqA+3+ggZTa5CalfWir3/fYlSq3pZ7LDTCxCJAAKxCEQ9XDrhhBOkVatW6QBn51fyRFHyXUo33HCDUS5NhkuvvfZa+lTS9lf79u3TgdF3vvOdxv+2ePFiGTNmjCR/1qZNm/TPk6eW7rjjDvnzn/8sl156qSTf3+T6SjaW27Zty/yJKNfrqhafDNeSV+LBqzgC5CW/XHTv3t1pcdv+J9dO3E0G44qrRiB0/2uusZ5i6OPw2a5nc9f+T7JluwcIn2G3Feu5PtzkwkeTK3tz13tA7P1vL2oWQa2aOfG7HL2TiZ1r/2d7dZwNAQQQMBeIeriUPC3UsWNHuf3223cROeKII6Rbt24ye/ZsIy2T4dIHH3wgyUfuJU9GrV69Wu6+++502DR69Ghp0aJFus5dd90lF110key2226yaNGi9BqS15YtWyQZhiVPPiVPWn3iE58wuq6mDirTxpLNnFOqvQWTF2+0zZ7YdWNp2//kutmUqA7AVcXWbFDsrqH7v1nwOj8g9norYnrr2dy1/5N82u4BilgDta6pnuuDXJVNwP56Xe8Bsfe/vahZBPcVM6dqR2GXnZ1r/+uvhEgEEEDATSDq4VLoJ5d2TsWaNWvSgVEymEo+Bi953XvvvXLeeedJnz59dhlsTZs2TZKP/Pvd734n+++/v1Nmy/RIPI9SO6XaWzB58Ubr/cS2/U+u/aQEV1z9CNQ+q23/53GNZVqTPg6fLczdzGO/B1AfbvURMppchdSurBV7//sSpVb1sthhpxcgEgEEYhGIerjU3HcuDRo0SCZPnmyUS5Mnl6qd6KyzzpLnnntOHn744fSPn3zyyfR7lY477jj5+c9/vkPIbbfdJpdddpn85je/kS9+8YtG19XUQWXaWLIhcUq1t2Dy4o3W+4lt+59c+0kJrrj6Eah9Vtv+z+May7QmfRw+W5i7mcd+D6A+3OojZDS5CqldWSv2/vclSq3qZbHDTi9AJAIIxCIQ9XBp6tSpcv3118vSpUulS5cujTlLPrpu2LBhVt9vpB0uJXErVqyQp556Kl1/w4YNcthhh8m//uu/SjJM+vgrGTZdd911snDhQvnc5z7nVGNl2liyIXFKtbdg8uKN1vuJbfufXPtJCa64+hFguBTSlT4OqV1ZC3M3c9s9gNtq4aOpj/Dm2hXJlVZOHxd7/+tlakdSq3pZ7LDTCxCJAAKxCEQ9XHr22WcleTopeYJpwoQJjTkbP368LF68WJYsWSKdO3eW5HNik4+wa9++vXTo0KFqbmsNl95//31p06ZN+j1KH38lf9EOGTJEevXqtcNH4CUfi3fffffJnXfeKds/VzW5huRppoaGhvS6kn+7vMq0sWRD4pJpf7HkxZ+t7zPb9j+59pMRXHH1I1D7rLb9n8c1lmlN+jh8tjB3M4/9HkB9uNVHyGhyFVK7slbs/e9LlFrVy2KHnV6ASAQQiEUg6uFSkqSJEyfK/PnzZejQodKjRw959NFHZdGiRTJ27FgZN25cmscnnnhCRo4cucN/S/77smXL0n+S19y5c2XdunVy5plnpv9/8iTUwIED0/87GVQlH2fXv39/2WeffaRly5by/PPPy4IFC9I/nzVrlvTs2bOxZl5++eV06JQMkJJ127Ztm15jEnPttdfKMccc41xfZdpYsiFxTreXE5AXL6xBTmrb/+TaT1pwxdWPQO2z2vZ/HtdYpjXp4/DZwtzNPPZ7APXhVh8ho8lVSO3KWrH3vy9RalUvix12egEiEUAgFoFSDpd+/etfy8033yz3339/s3n48MMPZebMmenwZu3atdK1a1cZPny4JE8ibX86qKnh0rRp02T69OlV1+jTp0/j00ivvPKK/PKXv0w3c8kayZqdOnWSQw89VEaNGiX77rvvLudYvXq1/OxnP0uHV5s3b5YvfOELMmbMGDn88MObfU8mB5RpY8mGxCSj4Y8hL+HNs1rRtv/JdVbyO54HV1z9CDBcCulKH4fUrqyFuZu57R7AbbXw0dRHeHPtiuRKK6ePi73/9TK1I6lVvSx22OkFiEQAgVgEch0uJU/tDB48uPEJoGqod911l8ybNy8dJvEyFyjTxpINiXleQx5JXkJqZ7uWbf+T62z9t58NV1z9CNQ+q23/53GN/297dwJtRXEufvsFFEHwKkPgH1AJooIxTMagxhhFUYkDAlGiMjkFjYJCVBCCiEMAI5AooMI1eoEE8YqCEREURFHjgBoBhYAhiuIQUCMEERm/VXW+cy7DPoeuqu7a3dW/XotlAl1dVc9b1bt6v7u7s1Qn89h/tDB3Mw/9HMD4cBsfPksTK5/aJXWFPv+TEmWs2stih529ACURQCAUgaIml9T7htTj6dSf8jZ1R5B6VNyyZctCMffSjywtLFmQeBkSxpUQF2Oy1BQwnf/EOpnQ4YprMgIkl3y6Mo99apfUhbmbuekawK02/6UZH/7NbWskVrZy9uVCn//2MhWXZKzay2KHnb0AJRFAIBSB1CeXRowYIQ8//LAsWrQoFHMv/cjSwpIFiZchYVwJcTEmS00B0/lPrJMJHa64JiNAcsmnK/PYpzbJpTi0TdcAcdTp8xjMSZ/abnURKzc/m9Khz38bkyhlGKtRlArvgx129gKURACBUAS8J5dmzJhRZnfTTTdJu3bt9J/dt+3bt8u//vUvUe9XatCggUyfPj0Ucy/9yNLCkgWJlyFhXAlxMSZLTQHT+U+skwkdrrgmI1DxUU3nfzHamKU6mcf+o4W5m3no5wDGh9v48FmaWPnULqkr9PmflChj1V4WO+zsBSiJAAKhCHhPLqlH4VWqVCmS344dO6Rq1aoyevToggmoSAfJ6U5ZWliyIEnnICUu6YxLlFaZzn9iHUXVfB9czc2ilMCV5FKUcRLXPoy3uCSjHwfz6FaF9jRdA7jV5r8048O/uW2NxMpWzr5c6PPfXqbikoxVe1nssLMXoCQCCIQi4D25VHoHkkocDRo0SCeNTjvttD08VQLqoIMOkpYtW0qtWrVC8fbWjywtLFmQeBsWRhURFyOuVO1sOv+JdTLhwxXXZARILvl0ZR771C6pC3M3c9M1gFtt/kszPvyb29ZIrGzl7MuFPv/tZUguYZeUgP1xOUfa21ESAQTSJeA9ubRz9wcOHFhucildTNlrTZYWlnyopnN8EZd0xiVKq0znP7GOomq+D67mZlFK4Fqxkun8j2Ke530Yb/6jj7mbeejnAMaH2/jwWZpY+dQuqSv0+Z+UKGPVXhY77OwFKIkAAqEIFDW5FApiGvuRpYUlC5I0jiB+OZzOqERrlen8Zw5GczXdC1dTsWj740pyKdpIiWcvxls8jiZHwdxEa899TdcAbrX5L8348G9uWyOxspWzLxf6/LeXqbgkY9VeFjvs7AUoiQACoQgUNbn0z3/+Uz/64qc//anUrFlTm27evFnuuecemT9/vn7f0iWXXCLnnXdeKN7e+pGlhSULEm/Dwqgi4mLElaqdTec/sU4mfLjimoxAxUc1nf/FaGOW6mQe+48W5m7moZ8DGB9u48NnaWLlU7ukrtDnf1KijFV7WeywsxegJAIIhCJQ1ORSv379ZOHChbJgwQKpXLmyNr3zzjvloYcekv3331+2bNkiW7dulQcffFBOOOGEUMy99CNLC0sWJF6GhHElxMWYLDUFTOc/sU4mdLjimowAySWfrsxjn9oldWHuZm66BnCrzX9pxod/c9saiZWtnH250Oe/vUzFJRmr9rLYYWcvQEkEEAhFoKjJpVNPPVWOOeYYGTlypPZUySSVRDr88MNl4sSJsn79eunUqZN8//vflwkTJoRi7qUfWVpYsiDxMiSMKyEuxmSpKWA6/4l1MqHDFddkBCo+qun8L0Ybs1Qn89h/tDB3Mw/9HMD4cBsfPksTK5/aJXWFPv+TEmWs2stih529ACURQCAUgaIml1q1aiXdu3eX66+/Xnu+8cYb0q1bN/nd734nHTp00H83dOhQee655/TdTWzRBbK0sGRBEj2uPvckLj61463LdP4T63j9S4+GK67JCFR8VNP5X4w2ZqlO5rH/aGHuZh76OYDx4TY+fJYmVj61S+oKff4nJcpYtZfFDjt7AUoigEAoAkVNLh133HE6ifSb3/xGe44dO1bGjRunE0nf+c539N+NGjVK38W0ePHiUMy99CNLC0sWJF6GhHElxMWYLDUFTOc/sU4mdLjimowAySWfrsxjn9oldWHuZm66BnCrzX9pxod/c9saiZWtnH250Oe/vUzFJRmr9rLYYWcvQEkEEAhFoKjJpV/84heybt06+ctf/iKVKlXSiaZ9991X///STd3VpBZJ8+fPD8XcSz+ytLBkQeJlSBhXQlyMyVJTwHT+E+tkQocrrskIVHxU0/lfjDZmqU7msf9oYe5mHvo5gPHhNj58liZWPrVL6gp9/iclyli1l8UOO3sBSiKAQCgCRU0uqSRS//795f/9v/+nk0qrV6+WIUOGyEUXXVTme/rpp0uTJk3k/vvvD8XcSz+ytLBkQeJlSBhXQlyMyVJTwHT+E+tkQocrrskIkFzy6co89qldUhfmbuamawC32vyXZnz4N7etkVjZytmXC33+28tUXJKxai+LHXb2ApREAIFQBIqaXFKIkyZNkhkzZmjP9u3bS69evcps1TuYfvWrX+l3Ml144YWhmHvpR5YWlixIvAwJ40qIizFZagqYzn9inUzocMU1GYGKj2o6/4vRxizVyTz2Hy3M3cxDPwcwPtzGh8/SxMqndkldoc//pEQZq/ay2GFnL0BJBBAIRaDoyaVQINPWjywtLFmQpG30lLSHuKQzLlFaZTr/iXUUVfN9cDU3i1ICV5JLUcZJXPsw3uKSjH4czKNbFdrTdA3gVpv/0owP/+a2NRIrWzn7cqHPf3uZiksyVu1lscPOXoCSCCAQigDJpVAiuVs/srSwZEGSzkFIXNIZlyitMp3/xDqKqvk+uJqbRSmBK8mlKOMkrn0Yb3FJRj8O5tGtSC65WVE6WQHmcrK+eZz/SYkyVu1lscPOXoCSCCAQikDRk0ubN2+WiRMnyuzZs+X999+XTZs2ydKlS7Wv+qCaOnWq9OjRQw477LBQzL30w/TLZS+NKqcSFiTF1C+/buKSzrhEaZXp/CfWUVTN98HV3CxKCVwrVjKd/1HM87wP481/9DF3Mw/9HMD4cBsfPksTK5/aJXWFPv+TEmWs2stih529ACURQCAUgaImlzZs2KATRyqZVKdOHalSpYqsXbtWli1bpn3Vv5900knStWtXueGGG0Ix99KPLC0sWZB4GRLGlRAXY7LUFDCd/8Q6mdDhimsyAiSXfLoyj31ql9SFuZu56RrArTb/pRkf/s1tayRWtnL25UKf//YyFZdkrNrLYoedvQAlEUAgFIGiJpeGDx+u71oaPHiwTiCNHTtW7r333rLkkkK+8sorZc2aNTJ9+vRQzL30I0sLSxYkXoaEcSXExZgsNQVM5z+xTiZ0uOKajADJJZ+uzGOf2iSX4tA2XQPEUafPYzAnfWq71UWs3PxsSoc+/21MopRhrEZRKrwPdtjZC1ASAQRCEShqcqlt27ZyxBFHyIQJE7SnSi6NGzdul+TSHXfcITNnzpRXX301FHMv/cjSwpIFiZchYVwJcTEmS00B0/lPrJMJHa64JiNAcsmnK/PYpzbJpTi0TdcAcdTp8xjMSZ/abnURKzc/m9Khz38bkyhlGKtRlEgu2SthF7cdx0MAgXQJFDW51Lx5c/1YvBtvvFGrFEoujRgxQqZMmSKLFy9Ol1zKW5OlhSWLuXQOJuKSzrhEaZXp/CfWUVTN98HV3CxKCVwrVjKd/1HM87wP481/9DF3Mw/9HMD4cBsfPksTK5/aJXWFPv+TEmWs2stih529ACURQCAUgaIml0455RRp1aqV/OEPfyg3uXT55ZfL6tWrZc6cOaGYe+lHlhaWLEi8DAnjSoiLMVlqCpjOf2KdTOhwxTUZAZJLPl2Zxz61S+rC3M3cdA3gVpv/0owP/+a2NRIrWzn7cqHPf3uZiksyVu1lscPOXoCSCCAQikBRk0vqXUtPPPGETJs2TZo2bbrHnUuvv/669OzZU9/dNHDgwFDMvfQjSwtLFiRehoRxJcTFmCw1BUznP7FOJnS44pqMAMkln67MY5/aJJfi0DZdA8RRp89jMCd9arvVRazc/GxKhz7/bUyilGGsRlEqvA922NkLUBIBBEIR8J5cUkmidu3ayWmnnSafffaZdOrUSTZt2iTdu3eXVatWyTPPPCMjR46Ut99+W6ZOnSoHHnigTkDVqVMnFHMv/cjSwpIFiZchYVwJcTEmS00B0/lPrJMJHa64JiNAcsmnK/PYpzbJpTi0TdcAcdTp8xjMSZ/abnURKzc/m9Khz38bkyhlGKtRlEgu2SthF7cdx0MAgXQJeE8uNWvWTHr37q3/qG3lypXSv39/effdd8tkKlWqJDt27JCjjjpKJ5qaNGnipLZ161YZP368PPbYY7J27Vpp2LChdOvWTbp27Sqqroq2WbNmyfPPP6/f+fTBBx9IvXr1ZMGCBQWL3HnnnbJw4UL56KOPZOPGjVK/fn1p06aNXH311XLwwQeXW83mzZvl3HPP1ce/6qqrpF+/fk79VYWztLBkMecc7kQOQFwSYfVyUNP5T6yTCQuuuCYjUPFRTed/MdqYpTqZx/6jhbnqIUmqAAAgAElEQVSbeejnAMaH2/jwWZpY+dQuqSv0+Z+UKGPVXhY77OwFKIkAAqEIFD25VAr5zjvv6ATO+vXrpUaNGtK8eXP9PqY4NvX4vUcffVS6dOkiLVq0kJdeeklmz54tffr0KUtylVePuqNKte3oo4/WyZ/KlSuXm1xSySr1eL9DDz1U90ElmdQj/1RySyW2DjnkkILVjBs3Th544AGdkCK5FEfEOUYcAiwU41AszjFMLyyJdTJxwhXXZARILvl0ZR771C6pC3M3c9M1gFtt/kszPvyb29ZIrGzl7MuFPv/tZSouyVi1l8UOO3sBSiKAQCgCqUkuJQW6bNky6dixo1x22WUyYMCAsmr69u0r8+bN03/U3UjlbZ9++qn+9ypVqpQ9uq+8O5cKHWPJkiVy/vnnS69eveT666/fYxeVgDrnnHPkmmuukVGjRpFcSmogcFxjARaKxmSpKWB6YUmskwkdrrgmI1DxUU3nfzHamKU6mcf+o4W5m3no5wDGh9v48FmaWPnULqkr9PmflChj1V4WO+zsBSiJAAKhCASfXBo9erR+JN78+fOlQYMGZXF788035eKLL5ZbbrlF/zfKVvpeKJPk0pdffiknnHCCXHjhhXLrrbfuUc2VV14pX3/9tYwYMUK/h4o7l6JEgn18CLBQ9KGcTB2mF5bEOpk44IprMgIkl3y6Mo99apfUhbmbuekawK02/6UZH/7NbWskVrZy9uVCn//2MhWXZKzay2KHnb0AJRFAIBSBoiSX1DuP1J+om3ov0sSJE6Puvst+6o6lFStW6Efh7byp9xy1bNlSOnfuLL/97W8jHTtKcmn79u3y1VdfybZt2+Tjjz8W9cg7lYxS/23Xrt0u9cydO1euvfZamT59un6MHsmlZpHiwE5+BFgo+nFOohbTC0tinUQU+II0GVVc9+ZqOv/3dry8/zvnR/8jAHM389DPAYwPt/HhszSx8qldUlfo8z8pUcaqvSx22NkLUBIBBEIRKEpyyRRPJZfU4+1sNvXIuapVq8rjjz++R3F1R5F6l5J631GULUpyafXq1TpJVLrVqlVL3410ySWX7FLFN998I2effbaceuqpot4JVVouzjuXduzYoZNWad+UhdqqV6+e9qbmqn3EpXjhbtbMLdGqLixN5j+xTibWuOJqI+B7/tu0MU9lmMf+o51nc9f5X/rlsskawH+E3WrM8/hwk/NfmliZm7ueA0yvAcxbGGYJxqp9XLGLz851/tu3hJIIIICAm0BRkks9e/aUHj16GLXc5E6nnQ+s7haqW7euTJ06dY/6TjnlFDnkkENk8uTJkdoSJbn07bffinrknrozauXKlTJz5kydbLr66qulcuXKZfWox/U9+uijMmfOHPmv//ovkksklyKNQZ87sVD0qb1rXa4LS9MLS2KdTKxxxdVGwPf8t2ljnsowj/1HO8/mrvOf5JL/8UqN5QvkeS7bjgvXc4DpNYBtO0Mrx1i1jyh28dm5zn/7llASAQQQcBMoSnKpd+/eov742HzfubR7nz755BNRbVCJqX79+ul//uc//ykdOnTQ73u64IIL9N8lceeSOu4xxxzjg9mpDm6lduJLrDBxSYw28QObPhKDWCcTElxxTUag4qOazv9itDFLdTKP/UcLczfz0M8BjA+38eGzNLHyqV1SV+jzPylRxqq9LHbY2QtQEgEEQhEIPrm0t3cuderUSYYNGxYpnlHuXCp0oCuuuEKWL18uL774ov7nX/3qVzrBpB7Hpx75p7bPPvtMunbtKt26dZNLL71U321VrVq1SO0qtFOWFpYsSKzDnGhB4pIob6IHN53/xDqZcOCKazICFR/VdP4Xo41ZqpN57D9amLuZh34OYHy4jQ+fpYmVT+2SukKf/0mJMlbtZbHDzl6AkgggEIpA8MmlUaNGyYQJE2T+/PnSoEGDsripR9ddfPHFMmTIEJ3UibLZJpdUucWLF8uiRYt0Needd56UfgiXV+/9998vbdu2jdKsgvtkaWHJgsQ6zIkWJC6J8iZ6cNP5T6yTCQeuuCYjUPFRTed/MdqYpTqZx/6jhbmbeejnAMaH2/jwWZpY+dQuqSv0+Z+UKGPVXhY77OwFKIkAAqEIBJ9cWrp0qai7k9QdTAMGDCiLW9++fWXu3Lkyb948qV+/vqhnxapH2NWqVUtq165dML4VJZfWr18v1atXl3333XeXsurDVj36rlWrVmXvdnr11Vdlw4YNu+z3xRdf6ERX+/bt5dxzz9X7q7uXbLcsLSxZkNhGOdlyxCVZ3ySPbjr/iXUy0cAV12QESC75dGUe+9QuqQtzN3PTNYBbbf5LMz78m9vWSKxs5ezLhT7/7WUqLslYtZfFDjt7AUoigEAoAt6TS8WAGzRokDz++OPSpUsXad68ubz88svy9NNP6/c+9enTRzfptddekx49euzyd+rvFy5cqP+obdq0abJu3Tq5/PLL9f9Xd0J17NhR/2+VqBo6dKhODjVq1EiqVKki7733nsyYMUP/+8SJE6VFixbldp93LonwAsNizI7y62ShmK54mLTG9MKSWJvoRt8X1+hWJnviWrGW6fw3sc/jvow3/1HH3M089HMA48NtfPgsTax8apfUFfr8T0qUsWovix129gKURACBUARykVzasmWLjB8/XieY1qxZIw0bNtSPwlN3IpW+86i85NKYMWNk7NixBePdpk2bsruRPvzwQ7nvvvv0gk7VoeqsV6+eHHfccdKrVy9p3LhxhWOG5BLJpbSdVFgopi0i0dtjemFJrKPbmuyJq4lW9H1xJbkUfbS478l4czc0PQLmpmK77m+6BnCrzX9pxod/c9saiZWtnH250Oe/vUzFJRmr9rLYYWcvQEkEEAhFIBfJpVCCZdKPLC0sWZCYRNbfvsTFn3XcNZnOf2IddwRKjocrrskIVHxU0/lfjDZmqU7msf9oYe5mHvo5gPHhNj58liZWPrVL6gp9/iclyli1l8UOO3sBSiKAQCgCJJdCieRu/cjSwpIFSToHIXFJZ1yitMp0/hPrKKrm++BqbhalBK4kl6KMk7j2YbzFJRn9OJhHtyq0p+kawK02/6UZH/7NbWskVrZy9uVCn//2MhWXZKzay2KHnb0AJRFAIBQBkkuhRJLkUqCRLF63WCgWz961ZtMLS2LtKl64PK64JiNAcsmnK/PYp3ZJXZi7mZuuAdxq81+a8eHf3LZGYmUrZ18u9PlvL0NyCbukBOyPyznS3o6SCCCQLgGSS+mKR2ytydLCkg/V2MIe64GIS6ycXg9mOv+JdTLhwRXXZARILvl0ZR771Ca5FIe26Rogjjp9HoM56VPbrS5i5eZnUzr0+W9jEqUMYzWKUuF9sMPOXoCSCCAQigDJpVAiuVs/srSwZEGSzkFIXNIZlyitMp3/xDqKqvk+uJqbRSmBK8mlKOMkrn0Yb3FJRj8O5tGtCu1pugZwq81/acaHf3PbGomVrZx9udDnv71MxSUZq/ay2GFnL0BJBBAIRYDkUiiRJLkUaCSL1y0WisWzd63Z9MKSWLuKFy6PK67JCJBc8unKPPapXVIX5m7mpmsAt9r8l2Z8+De3rZFY2crZlwt9/tvLkFzCLikB++NyjrS3oyQCCKRLgORSuuIRW2uytLDkQzW2sMd6IOISK6fXg5nOf2KdTHhwxTUZAZJLPl2Zxz61SS7FoW26BoijTp/HYE761Hari1i5+dmUDn3+25hEKcNYjaJUeB/ssLMXoCQCCIQiQHIplEju1o8sLSxZkKRzEBKXdMYlSqtM5z+xjqJqvg+u5mZRSuBKcinKOIlrH8ZbXJLRj4N5dKtCe5quAdxq81+a8eHf3LZGYmUrZ18u9PlvL1NxScaqvSx22NkLUBIBBEIRILkUSiRJLgUayeJ1i4Vi8exdaza9sCTWruKFy+OKazICJJd8ujKPfWqX1IW5m7npGsCtNv+lGR/+zW1rJFa2cvblQp//9jIkl7BLSsD+uJwj7e0oiQAC6RIguZSueMTWmiwtLPlQjS3ssR6IuMTK6fVgpvOfWCcTHlxxTUaA5JJPV+axT22SS3Fom64B4qjT5zGYkz613eoiVm5+NqVDn/82JlHKMFajKBXeBzvs7AUoiQACoQiQXAolkrv1I0sLSxYk6RyExCWdcYnSKtP5T6yjqJrvg6u5WZQSuJJcijJO4tqH8RaXZPTjYB7dqtCepmsAt9r8l2Z8+De3rZFY2crZlwt9/tvLVFySsWovix129gKURACBUARILoUSSZJLgUayeN1ioVg8e9eaTS8sibWreOHyuOKajADJJZ+uzGOf2iV1Ye5mbroGcKvNf2nGh39z2xqJla2cfbnQ57+9DMkl7JISsD8u50h7O0oigEC6BEgupSsesbUmSwtLPlRjC3usByIusXJ6PZjp/CfWyYQHV1yTESC55NOVeexTm+RSHNqma4A46vR5DOakT223uoiVm59N6dDnv41JlDKM1ShKhffBDjt7AUoigEAoAiSXQonkbv3I0sKSBUk6ByFxSWdcorTKdP4T6yiq5vvgam4WpQSuJJeijJO49mG8xSUZ/TiYR7cqtKfpGsCtNv+lGR/+zW1rJFa2cvblQp//9jIVl2Ss2stih529ACURQCAUAZJLoUSS5FKgkSxet1goFs/etWbTC0ti7SpeuDyuuCYjQHLJpyvz2Kd2SV2Yu5mbrgHcavNfmvHh39y2RmJlK2dfLvT5by9Dcgm7pATsj8s50t6OkgggkC4BkkvpikdsrcnSwpIP1djCHuuBiEusnF4PZjr/iXUy4cEV12QESC75dGUe+9QmuRSHtukaII46fR6DOelT260uYuXmZ1M69PlvYxKlDGM1ilLhfbDDzl6AkgggEIoAyaVQIrlbP7K0sGRBks5BSFzSGZcorTKd/8Q6iqr5Priam0UpgSvJpSjjJK59GG9xSUY/DubRrQrtaboGcKvNf2nGh39z2xqJla2cfbnQ57+9TMUlGav2sthhZy9ASQQQCEWA5FIokSS5FGgki9ctForFs3et2fTCkli7ihcujyuuyQiQXPLpyjz2qV1SF+Zu5qZrALfa/JdmfPg3t62RWNnK2ZcLff7by5Bcwi4pAfvjco60t6MkAgikS4DkUrriEVtrsrSw5EM1trDHeiDiEiun14OZzn9inUx4cMU1GQGSSz5dmcc+tUkuxaFtugaIo06fx2BO+tR2q4tYufnZlA59/tuYRCnDWI2iVHgf7LCzF6AkAgiEIkByKZRI7taPLC0sWZCkcxASl3TGJUqrTOc/sY6iar4PruZmUUrgSnIpyjiJax/GW1yS0Y+DeXSrQnuargHcavNfmvHh39y2RmJlK2dfLvT5by9TcUnGqr0sdtjZC1ASAQRCESC5FEokSS4FGsnidYuFYvHsXWs2vbAk1q7ihcvjimsyAiSXfLoyj31ql9SFuZu56RrArTb/pRkf/s1tayRWtnL25UKf//YyJJewS0rA/ricI+3tKIkAAukSILmUrnjE1posLSz5UI0t7LEeiLjEyun1YKbzn1gnEx5ccU1GgOSST1fmsU9tkktxaJuuAeKo0+cxmJM+td3qIlZufjalQ5//NiZRyjBWoygV3gc77OwFKIkAAqEIkFwKJZK79SNLC0sWJOkchMQlnXGJ0irT+U+so6ia74OruVmUEriSXIoyTuLah/EWl2T042Ae3arQnqZrALfa/JdmfPg3t62RWNnK2ZcLff7by1RckrFqL4sddvYClEQAgVAESC6FEkmSS4FGsnjdYqFYPHvXmk0vLIm1q3jh8rjimowAySWfrsxjn9oldWHuZm66BnCrzX9pxod/c9saiZWtnH250Oe/vQzJJeySErA/LudIeztKIoBAugRILqUrHrG1JksLSz5UYwt7rAciLrFyej2Y6fwn1smEB1dckxEgueTTlXnsU5vkUhzapmuAOOr0eQzmpE9tt7qIlZufTenQ57+NSZQyjNUoSoX3wQ47ewFKIoBAKAIkl0KJ5G79yNLCkgVJOgchcUlnXKK0ynT+E+soqub74GpuFqUEriSXooyTuPZhvMUlGf04mEe3KrSn6RrArTb/pRkf/s1tayRWtnL25UKf//YyFZdkrNrLYoedvQAlEUAgFIHgk0tbt26V8ePHy2OPPSZr166Vhg0bSrdu3aRr165SqVKlCuM4a9Ysef7552Xx4sXywQcfSL169WTBggUFy9x5552ycOFC+eijj2Tjxo1Sv359adOmjVx99dVy8MEHl5X56quvZPr06TJ//nxZuXKl3veQQw6Rc845R3r27Cn77bdfLGMrSwtLFiSxhDz2gxCX2Em9HdB0/hPrZEKDK67JCFR8VNP5X4w2ZqlO5rH/aGHuZh76OYDx4TY+fJYmVj61S+oKff4nJcpYtZfFDjt7AUoigEAoAsEnlwYPHiyPPvqodOnSRVq0aCEvvfSSzJ49W/r06SO9e/euMI7du3eXd955R44++midXKpcuXK5ySWVrGra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C8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhAnBh3MAACAASURBVCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSJJcCjSSxesWC8Xi2bvWbPrFErF2FScJkowgrjaupvPfpo48leH86D/amLuZh34OYHy4jQ+fpYmVT22SSy7ajFV7PeywsxegJAIIhCJAcimUSO7WjzfffFP/TaVKlVLfwx07dmSmranHjLGBxCVGTMNDqXnbunVrw1L/t7vp/CfW1tQVFsQVVxsB3/Pfpo15KsM89h/tPJu7zn8VLdM1gP8Iu9WY5/HhJue/NLEyN3c9B4Q+/81Fo5VgrEZzKrQXdvHZuc5/+5ZQEgEEEHATILnk5pfa0iwsUxsaGobAXgVcF5bM/70SswMCqRVg/qc2NDQMgcQFXOd/HpJLiQeBChAoooDrOYBrgCIGj6oRcBRwnf+O1VMcAQQQsBYguWRNR0EEEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABBBBAIH8CJJfyF3N6jAACCCCAAAIIIIAAAggggAACCCCAAAIIIIAAAghYC5BcsqajIAIIIIAAAggggAACCCCAAAIIIIAAAggggAACCCCQPwGSS/mLOT1GAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABBKwFSC5Z01EQAQQQQAABBBBAAAEEEEAAAQQQQAABBBBAAAEEEMifAMml/MWcHiOAAAIIIIAAAggggAACCCCAAAIIIIAAAggggAAC1gIkl6zpKIgAAggggAACCCCAAAIIIIAAAggggAACCCCAAAII5E+A5FL+Yk6PEUAAAQQQQAABBBBAAAEEEEAAAQQQQAABBBBAAAFrAZJL1nQURAABBBBAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQTyJ0ByKX8xp8cIIIAAAggggAACCCCAAAIIIIAAAggggAACCCCAgLUAySVrOgoigAACCCCAAAIIIIAAAggggAACCCCAAAIIIIAAAvkTILmUv5jTYwQQQAABBBBAAAEEEEAAAQQQQAABBBBAAAEEEEDAWoDkkjUdBRFAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQACB/AmQXMpfzOkxAggggAACCCCAAAIIIIAAAggggAACCCCAAAIIIGAtQHLJmo6CCCCAAAIIIIAAAggggAACCCCAAAIIIIAAAggggED+BEgu5S/mqehxt27dZOHChXLuuefKyJEjy23T5s2b9T4ffPCBXHXVVdKvX7899p01a5ZMmDBBVq5cKQceeKD87Gc/k759+0qNGjVS0dcsNaK8uHz11Vcyffp0mT9/vnbeuHGjHHLIIXLOOedIz549Zb/99tujm6+88orcfffdsmzZMqlevbq0bdtWbrzxRqldu3aWSDLR1q1bt8r48ePlsccek7Vr10rDhg1FxbJr165SqVKlcvtgE9fSg0WZm5nAq6CRvly7d+8ur7/++h4tqVKliixdujTrjHu035erqnjbtm3y8MMPy7Rp0+T999+XqlWrSpMmTeSaa66Rk046KThbOrSngK/xduedd+p1zUcffaQ/I+vXry9t2rSRq6++Wg4++OBchcaX+c6oefhMytUgKtBZ1qjZGAF7u8YzuT7gGi8bMffZStvPF9VGNZ6ef/55Wbx4sf5uoV69erJgwYKCzTf5TM/KOt6XnelncxbmeRrtTj31VPn444/3GL+HHnqoPPvssz6nJXUhgAACQnKJQeBdYMaMGXLrrbfqL1/2llwaN26cPPDAA3rfQsmlv/zlLzphcfzxx8vZZ5+tF4qTJk2SY489Vh566KEKv1j33vGUV1hRXFRSSX0Ze+KJJ2rrmjVr6i/RZs6cKcccc4xMnjxZ1BfhpZv6ovzSSy+Vpk2byvnnny9ffvmlPPjgg9KgQQP9JW+1atVSrpGt5g0ePFgeffRR6dKli7Ro0UJeeuklmT17tvTp00d69+5dbmdM47rzgfY2N7MlWLi1vlzVReny5cvlN7/5zS4NqVy5sj5Hhrb5ct2+fbtce+218sILL0inTp2kefPm8s0338g//vEP/b8vuOCC0GjpTwEBX+NNJfPVZ566qFc/blFJJvV5p76QUIl/9YOMvGy+zPP2mZSX8VOon6xRsxH9vV3jmVwfcI2XjZj7bqXt54tqp1pvv/POO3L00Ufr7wzUOru85JLJZ3pW1vG+7Ew+m7Myz9Nop5JL++67r/4R086bWoO2a9fO99SkPgQQyLkAyaWcDwDf3V+/fr2+s0jd7TJq1KgKk0vqixl1Z4xKaqh9d08uqV+pqrthVMJi6tSpZcmNKVOm6OTV2LFj5fTTT/fdxUzWt7e4qFiobfcvx9SdSffee+8e1h07dpR169bJU089Jfvvv78uq77g7dWrlwwcOFAuueSSTDqlsdHqzjDlfdlll8mAAQPKmqju3ps3b57+o36ZV2gzjWvpMfY2N9PoZNomn67qonTVqlXlXuCatj3N+/t0VT80UL88nThxov7BAVv+BHyOt0K6S5Ys0T+wUJ99119/fS4CUAzzPHwm5WLwlNNJ1qjZiP7e4qR6EfX6gGu8bMTcdytdPl9UWz/99FN9TaR+EGmz9i7vM93mWHmw29tnc1bmeTHG3d7s1PhRySV1l7x6QgMbAgggUGwBkkvFjkDO6r/tttvk5ZdflieffFL/cryiO5euvPJK+frrr2XEiBFy2mmn7ZFcUndnXH755frLQ3WxUrqphcpxxx0nJ598svzhD3/ImbBdd03isnMN6m6LDh06yHXXXVf2qxn16Kn27dsXvGvmjDPO0I8uVHfZsMUjMHr0aP1IPHUXkkq0lm5vvvmmXHzxxXLLLbfo/5psheK6c/m9zU2TutK6r0/X0otSFUN1Z436xVlFjzNMq1mUdvlyVXctqc8N9Tlzzz33iPr/pbZR2sk+YQj4Gm/laam7dk844QS58MIL9Y9e8rAVwzwPn0l5GDvl9ZE1ajaiv7c4mVwfcI2XjZj7bmWcny82CaHyPtOzsI4vht3ePpuzMs/TaKfmXmlySf2Y7ttvv9VPlmFDAAEEiiVAcqlY8jms991339W/4L3//vt14kc9Pqa85NLcuXP144zUe37UF62FkkvqOL///e/l6aeflsMOO2wX0Ysuukg+//xznjcbYZyZxGX3w7344otyxRVX6C/N1JdnalOJwxtuuEE/znD3d5qov58zZ468/fbbuzxGL0Iz2aUcAXXH0ooVK/Sj8HbeVJK1ZcuW0rlzZ/ntb39r5FcorqUHiDI3jSpL6c4+XdVF6VtvvSX77LOPbNq0SV8cqESsmi916tRJqZBds3y5qkffqUelqvf0rVmzRh5//HGdXPrud7+rf6hQer6y6wWlsiLga7yVeqgkpnqXnXrXl3oOvnp8qHrkjvpvXh5R4ts8L59JWZlzcbeTNWrcoskcL0qcTK4PuMZLJk5ZP2qcny9RkktRP9OzsI73bRflszkr8zyNdmouq+SSetfyjh07ZMuWLVKrVi393Zq69il9ckzW5zztRwCB7AiQXMpOrDLdUrU4U1/m1a5dWyeX1FZeckl9Aai+FFQfmOr5tqtXry6YXFK/kPvzn/8s6g6N3X+poe6kUXcBqBd2spUvYBKX3Y+iyqrHG6pHBKgFZN26dfUuf/zjH+V3v/udqGcoqxjvvKm/V/+u7l4r3Z/4uAmoR0dWrVpVf3m++6Z+Ma+eK64SfVG38uKqykedm1HrSvN+Pl3VoyLVYw3UfFEXCOpl1+pdLeoxlOq///Vf/5VmKqO2+XJV5yT1SFV1oVW9enWdUFI/VHjkkUdEvfNhyJAhop6nzxa2gK/xVqpYul4p/f9q/Kmxl6dHwfo0z9NnUtgztXDvWKNmI+pR42RyfcA1XjZi77uVcX6+REkuRf1Mz8I63qdd1M/mrMzzNNqpuafuDGvdurU0adJEP+3n+eef1z+6Vn+n3oet3sfEhgACCPgSILnkSzqQetQXn+qOiCibeklm6Yea+kLvjjvu0O/gUS+7Vlt5ySV167F6bJq6w0V9qVpecmnQoEH6JdkquaG+XN9569+/vzzxxBPy97//PdjHS+3cXx9x2T3mpbeIqwSgWqCXbuoX2uoxVLNnz5bGjRvvUqz0HU3qPUAHH3xwlGHEPnsRUL+GV4k69d6x3bdTTjlFJyjUAjPqVl5cVfmoczNqXWnez6drIQd1DlRzq3fv3voRk6FsvlzV+V99DqjPIHWhVfq+uK1bt+p3+am7S9TdfupuMbZwBXyNt1JB9VgS9YMXtU5auXKlzJw5U/84Rr1sWa2J8rD5NM/TZ1KWxw5r1GxEL+k4mVwfcI2XjTHju5Vxfr5ESS65fKanbR3v0y7qZ3NW5nka7cqbe+qpPuqH3MOHD9dPL2FDAAEEfAmQXPIlHUg96vFb6nbbKFunTp30+5LU84l/9rOfiXpUXd++fcuKFkou/fOf/9Tv8FHvibngggv0vty5tHftpOOyewv+9Kc/ye23317wPRImv0zce8/YY28Ccf6aqqK4mszNvbU5C//uy7UiC/XuOJWgLZQ4zIJhoTb6clU/TlCPVm3Tps0eydUxY8bI2LFjC95dmVVX2l1YwNd4K8//k08+0clM9SWWekxJHjZf5nn7TMry2GGNmo3oJR0nk+uDrNzRkI3IhtPKOD9foiSXdpcz/UxP0zrel53JZ3NW5nka7cqb1Rs2bJAf/vCHeu05atSocCY/PUEAgdQLkFxKfYjS1cD169frR6BF2dQdSscee6y+Y0k9Z1s9wq5atWplRdWvedUvQdSt5AcddJB+tN2vfvUrUYsS9Riv0hfaf/bZZ/rxRd26dZNLL71U36WhjrO35/SqZ9BGbWuU/qR5n6TjsnPf1ePX1C+NzjrrLBk5cuQev8be2zPV1R1NixYt4p1LMQ2ovT0HWiV5hw0bttfa9hZXk7m518oysIMv14ooVOw2btyo7+IMZfPl+re//U0nv9V5Sv2Kb+ft4YcflqFDh4pKpv7oRz8KhZZ+FBDwNd4qwlfvJVy+fLmod9nlYfNlnrfPpCyPHdao2Yhe0nEyuT7gGi8bY8Z3K+P6fFHttkkuqXImn+lpWsf7sjP5bM7KPE+jXUVzT/2w7gc/+IE8+OCDvqco9SGAQI4FSC7lOPi+uq4eB6Meg1bRdvPNN+vk0XnnnacfZVfRphYibdu21V/UqAXenXfeKR07diwroh5Ho34p9NOf/lTUY9jYCguYxKX0COqxhjfccIOcfPLJ+pf/hR4ppZKD6k419Sgv9UivnbczzjhDP+pQvUeGLR4B9aukCRMm6HeMNWjQoOyg6tFMF198caR3y0SJq8ncjKdnxT2KL9fyeqneYaDOY4cffrioZEgomy9X9ezx448/Xl9c7e5X+siIWbNm6eeUs4Ur4Gu8VSSovsBS739UP6rIw+bLPG+fSXkYOzv3kTVqNiJuEieT6wOu8bIRf9+tjOPzpbTNtsmlqJ/paVvH+7Iz+WzOyjxPo115c0899ltdP6onDakfAbMhgAACvgRILvmSznE96gsVdRfR7pt62br6ZUXPnj31+5fUOzFeffVVUbfz7rx98cUX+gvy9u3b6w/KVq1a6buXVBJJvVOmYcOG+iXtpe8zmDJlitx66636vT9nnnlmjuUr7rpJXNSR1F1g1113nf6lv0pm7P6eq51rUwtL9QtIlbTYf//99T+98MIL0qtXLxkwYICoXwCxxSOwdOlSUb+MU6bKtnRTj6BUMVOJ3fr164t6uap6nIN6wXzt2rXL9osaV5O5GU/PinsUX67qfKfeC7Tffvvt0mF19+Zdd90lv/71r/ULW0PZfLkqL/VYvGeffVamT58uzZo104RqHqi7mdSdsWpulN4hG4ov/dhVwNd4U5931atX3+PlyerHMuoRv2rdYvLuuyzH0Zd53j6TsjwmbNrOGtVGzX8Z0zhFvT7gGs9/LLNQo+vny859rCi5ZPKZnpV1vC87k8/mrMzzNNqpJNIBBxywx5NgSh81qBJi6tF4bAgggIAvAZJLvqSpZw+BQu9cKsRU3juX1L4zZszQX6ifcMIJ+gvDVatWycSJE6V169YyadIkvji0GHeF4qJ+da0eTai+BO/fv7/+Em3nTT0CUZmXbmphqZId6gtd9cWaShA+9NBDOsnx2GOP7VHeopkU2UlAPaZQPdauS5cu0rx5c3n55Zfl6aef1neOqTvI1Pbaa69Jjx49dvk707jujl7R3AwhQD5cVVzUu1jU+UvNI5XsUH+nkiJq/qi7bkoTtCGYqj74cFX1qM8Ddf5Rpmrs16hRQ8+T9957T//44PTTTw+FlH5UIOBjvKkkvXrUovoRTKNGjfTFvhpnao2iNrUuadGiRW7i5MPcdL2YG/zAO8oaNRsBLu8az+T6gGu8bMTadyttP19UOxcuXKj/qE09RWPdunVy+eWX6/+vnv5Q+iQUk8/0LK3jfdiZfjZnZZ6nzU5dz6inyKgfUh988MH6R9fPP/+8/qH2T37yE/1DYLUWZUMAAQR8CZBc8iVNPXsIxJFcUgdVd8eoD1D1uIUDDzxQf7mj7tpQ73BiMxcoFBe1gFHvxipvU3fOjBgxYpd//utf/6ofS7hs2TKdTFJ3md144436rjO2eAW2bNki48eP11+cr1mzRt/Np5KB6ld5pXdmFEou2cR155aHnlzy4aoM1WML3nnnHfn8889l27Zt+iJBPUJS3emnEiKhbT5cS81WrlypfdWXCerC6/vf/76ou2ZPOumk0FjpTzkCPsbbhx9+KPfdd5+89dZb+hys6qxXr55+NImax40bN85VfHyYm36BlasABNxZ1qjZCG5F13gm1wdc42Uj3j5bafv5oto4ZswY/YV8oU09TaX0DmOTz/QsreN92Nl8NmdhnqfN7t1335Vx48aJuqvqyy+/1Nf73/ve9/RTftRTgdQPgtkQQAABnwIkl3xqUxcCCCCAAAIIIIAAAggggAACCCCAAAIIIIAAAgggkHEBkksZDyDNRwABBBBAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQR8CpBc8qlNXQgggAACCCCAAAIIIIAAAggggAACCCCAAAIIIIBAxgVILmU8gDQfAQQQQAABBBBAAAEEEEAAAQQQQAABBBBAAAEEEPApQHLJpzZ1IYAAAggggAACCCCAAAIIIIAAAggggAACCCCAAAIZFyC5lPEA0nwEEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABBBBAwKcAySWf2tSFAAIIIIAAAggggAACCCCAAAIIIIAAAggggAACCGRcgORSxgNI8xFAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABnwIkl3xqUxcCCCCAAAIIIIAAAggggAACCCCAAAIIIIAAAgggkHEBkksZDyDNRwABBBBAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQR8CpBc8qlNXQgggAACCCCAAAIIIIAAAggggAACCCCAAAIIIIBAxgVILmU8gDQfAQQQQAABBBBAAAEEEEAAAQQQQAABBBBAAAEEEPApQHLJpzZ1ISAi3377rWzYsEEOOuggqVKlSpnJggULZP78+VK1alXp0qWLNGnSBC8EEAhAYMaMGda96Nixo3VZCiKAAAIIIIBAegS4BkhPLGgJAggggAACCCCAQDwCJJficeQoCEQWGDp0qDzxxBPy0ksvSY0aNXS5Rx99VIYMGSI7duzQ/1/9/bRp06Rx48aRj8uOCCCQToFmzZpJpUqVyuZ31FaqMsuWLYu6O/shgEAKBXr06GHVKjX/J06caFWWQgggkE4BrgHSGRdahUBSAqwBkpLluAgggAACaRIguZSmaNCWXAj87Gc/00mje++9t6y/bdu21V8+33XXXfLFF19I//79Re03fPjwXJjQSQRCFnj99detu9emTRvrshREAIHiC5x66qnWjXjuueesy1IQAQTSJ8A1QPpiQosQSFKANUCSuhwbAQQQQCAtAiSX0hIJ2pEbgWOPPVYuuOACGTBggO7zihUrpEOHDjJw4EDp2bOn/rsbbrhB3n77bZk7d25uXOgoAggggAACCCCAAAKhCnANEGpk6RcCCCCAAAIIIJBfAZJL+Y09PS+SQOvWreWiiy7SdyepbfLkyTJs2DCZOXNm2XuWRo8erR+Hs2jRoiK1kmoRQAABBBBAAAEEEEAgLgGuAeKS5DgIIIAAAggggAACaREguZSWSNCO3Aicc845UqtWLZ1UUpt6FvOqVavkhRdeKDMYNGiQLFiwQL+XiQ0BBMIT2Lx5s04gz549W95//33ZtGmTLF26VHf073//u0ydOlWfGw477LDwOk+PEEBAC7z33nt6/m/cuFE6duyICgIIBC7ANUDgAaZ7CBgIsAYwwGJXBBBAAIFUC5BcSnV4aFyIAvfdd5/cfffdcsYZZ0i1atXkySeflMsuu0xuvPHGsu6ef/75st9++8mf//znEAnoEwK5FtiwYYNOHKlkUp06daRKlSqydu1aWbZsmXZR/37SSSdJ165d9SMy2RBAICyBN954Q4YOHSorV64s61jp/Ff/dvnll8uoUaOkXbt2YXWc3iCQcwGuAXI+AOg+AiLCGoBhgAACCCAQmgDJpdAiSn9SL6DuUFBfGM+bN0927NghJ554oowZM0b2339/3Xb1ZdPZZ58tvXv31n/YEEAgLIHhw4fru5YGDx6sE0hjx46Ve++9tyy5pHp75ZVXypo1a2T69OlhdZ7eIJBzgcWLF0u3bt30j0vUD0nUZ766U7k0uaR4VFKpefPm8vvf/z7nWnQfgbAEuAYIK570BgFTAdYApmLsjwACCCCQBQGSS1mIEm0MUuA///mPVKpUSWrWrLlL/7788kv9pXLDhg3lgAMOCLLvdAqBPAu0bdtWjjjiCJkwYYJmUMmlcePG7fLl8h133KHfw/bqq6/mmYq+IxCcwC9/+UtRXy7NmDFDvvvd7xac/9dff70sWbJEnnnmmeD6T4cQQECEawBGAQL5FGANkM+402sEEEAgdAGSS6FHmP6lTmDgwIHSrFkz6dmzZ+raRoMQQCB5AXVHgnosXumjMAsll0aMGCFTpkzRX0KzIYBAOALHHnustG/fXlQCubzk8l133aUfi/v222+H03F6ggACwjUAgwCBfAuwBsh3/Ok9AgggEKoAyaVQI0u/UivQsmVL/cWy+mUyGwII5E/glFNOkVatWskf/vCHcr9cVu9cWb16tcyZMyd/QPQYgYAF1Nz/xS9+ob9kLi+5dPPNN8vTTz+t38vAhgAC4QhwDRBOLOkJAjYCrAFs1CiDAAIIIJB2AZJLaY8Q7QtOoHPnztK4cWP9sm42BBDIn4B619ITTzwh06ZNk6ZNm+7xWKzXX39d39moktClX0DnT4keIxCmQKdOnWSfffaRRx99tGByadu2bXLWWWdJ7dq15eGHHw4TgV4hkFMBrgFyGni6jcD/L8AagKGAAAIIIBCiAMmlEKNKn1ItMGvWLP2F8eTJk6VFixapbiuNQwCB+AU+++wzUReX6sXe3bt3l1WrVul3q4wcOVI/Bmvq1Kly4IEH6gRUnTp14m8AR0QAgaIJTJo0SYYNGyZXX3219OnTR79vrfSda5s3bxb1SEyVVLr99tvl/PPPL1o7qRgBBOIX4BogflOOiECWBFgDZClatBUBBBBAIKoAyaWoUuyHQEwC6iXeTz75pLz22mty5plnI1xl4QAAE4lJREFUylFHHaW/QK5UqdIeNXTs2DGmWjkMAgikSWDlypXSv39/effdd8uapc4BO3bs0OcElWhq0qRJmppMWxBAIAYBdWdS7969Zf78+VK/fn2pVq2afPjhh/LjH/9Yli9fLp9//rmcfvrpMmbMmBhq4xAIIJAmAa4B0hQN2oKAfwHWAP7NqREBBBBAIHkBkkvJG1MDArsINGvWTCeS1JfIO287J5fUv6n/v2zZMvQQQCBggXfeeUcWL14s69evlxo1akjz5s31+5jYEEAgXAH1GT9lyhR9h5JKNJeuB773ve/JRRddpB+JWegHJ+GK0DME8iHANUA+4kwvEahIgDUA4wMBBBBAIDQBkkuhRZT+pF5g+vTpkduoHp3FhgACCCCAAAJhCnzzzTdlyeWaNWuG2Ul6hQACWoBrAAYCAgjsLMAagPGAAAIIIBCCAMmlEKJIHxBAAAEEEEAAAQQQQAABBBBAAAEEEEAAAQQQQAABTwIklzxBUw0CCCCAQD4F1COubDb1WKyJEyfaFKUMAgikROCTTz6xbkmDBg2sy1IQAQQQQAABBIorwBqguP7UjgACCCDgR4Dkkh9nakFgDwH1Mu+ZM2fKihUrZMOGDaIeh3PkkUfKueeeK6eccgpiCCAQiMCpp566R0+2bNkia9eu1X+/zz77yEEHHSRfffWVbN26Vb9rpW7durLvvvvKc889F4gC3UAgnwKl71ix6T3vXbRRowwC6RfgGiD9MaKFCMQhwBogDkWOgQACCCCQdgGSS2mPEO0LTmDz5s3Sr18//aWxeqFnoS+W1ZfRv//976Vq1arB9Z8OIZB3gX//+99y6aWXSp06deTaa6+VFi1a6ISSOh8sWrRI7rnnHlH7PPTQQzrpxIYAAtkVGDNmjJ7fO29/+9vf5OWXX5bDDjtMWrdurc8FX3zxhbz11lvy/vvvy4knnqj/vnfv3tntOC1HAIE9BLgGYFAgkC8B1gD5ije9RQABBPIqQHIpr5Gn30UTuOuuu+SPf/yjnHTSSdKnTx9p3rx52RfLS5YsEbUIfemll+SKK66Q66+/vmjtpGIEEEhGYODAgaLuSFAv9t79S2dV4/bt26Vz585y1FFHyfDhw5NpBEdFAIGiCLzyyivSq1cvueOOO+S8887bow0zZsyQIUOGyIQJE+T4448vShupFAEEkhHgGiAZV46KQFYEWANkJVK0EwEEEEDARIDkkokW+yIQg4BKKn3nO9+Rxx9/vODR1N0LP//5z/Ujs1588cUYauQQCCCQJgH1hXGXLl3k17/+dbnNGjVqlEybNk3URSgbAgiEI3DhhReKepfS6NGjy+2Uurv5008/lalTp4bTcXqCAAL6h2VcAzAQEMivAGuA/MaeniOAAAIhC5BcCjm69C2VAq1atZIePXpU+MWy+tJp0qRJ8vbbb6eyDzQKAQTsBdTjrtq3b1/hXUk33XSTzJkzR9Tjs9gQQCAcAbUG6Nmzp348bnmbWgNMnjyZ+R9O2OkJAlqAawAGAgL5FmANkO/403sEEEAgVAGSS6FGln6lVkDdsXDooYfKyJEjy22jehze6tWr5ZFHHkltP2gYAgjYCajkskocP/DAA9KmTZs9DvLaa6/px2L+8Ic/lP/5n/+xq4RSCCCQSoGf/OQneg0wZcqUctt30UUXyYcffqjfy8SGAALhCHANEE4s6QkCNgKsAWzUKIMAAgggkHYBkktpjxDtC07gr3/9q1x11VUyYsQIOeuss/bo31NPPSWDBg2S+++/X0444YTg+k+HEMi7wDvvvCPdu3eXTZs26XeqqF8x1q5dW7788kt9p4JKLlWrVk3+9Kc/ydFHH513LvqPQFAC6l1Lf/7zn6Vjx476vYvqEXml2yeffKLfu6jeu9S1a1cZPHhwUH2nMwjkXYBrgLyPAPqfdwHWAHkfAfQfAQQQCFOA5FKYcaVXKRYYO3as/gJZXWAefvjhu3yxrO5m+Mc//iE//vGPRT06a+etUqVKcs0116S4ZzQNAQSiCixdulRuu+22go++VHN/yJAhctRRR0U9HPshgEBGBL7++mv55S9/KW+99ZZUqVJFv3+lNLms3rW4bds2fdfif//3f8v++++fkV7RTAQQiCLANUAUJfZBIFwB1gDhxpaeIYAAAnkWILmU5+jT96IINGvWzKpelVxatmyZVVkKIYBAOgU+/vhjWb58uWzYsEFq1qwpTZs2lYYNG6azsbQKAQRiEdi+fbtMmzZNnnzySVmxYkXZ/D/yyCOlQ4cO8vOf/1wqV64cS10cBAEE0iPANUB6YkFLECiWAGuAYslTLwIIIIBAUgIkl5KS5bgIlCPw+uuvW9sUej+L9cEoiAACCCCAAAIIIIAAAl4EuAbwwkwlCCCAAAIIIIAAAh4FSC55xKYqBFwE1J0N69ev3+X9DC7HoywCCBRXYN26dfLMM8/I3//+97I7F9Svms844ww58MADi9s4akcAAW8CO3bsEHV3MhsCCCBQSIBrAMYFAuEKsAYIN7b0DAEEEMiLAMmlvESafmZeQD2nfdy4cTwaL/ORpAMIiDz++ONy++23y6ZNm0RdVO68Va9eXW6++Wbp3LkzVAggEKCAmvOPPPKIPPHEEzq5rM4D1apVE5Vc7tixo3Tp0oVkU4Bxp0sI2ApwDWArRzkE0ifAGiB9MaFFCCCAAAJuAiSX3PwojYA3AS4svVFTEQKJCixYsECuvPJKfXdS9+7dRT3usk6dOvLFF1/IwoULZdKkSfouxfHjx8tJJ52UaFs4OAII+BXYvHmznv+vvvqqTiB997vflbp168rnn38un332mah3MRx//PF6/letWtVv46gNAQRSKcA1QCrDQqMQMBZgDWBMRgEEEEAAgQwIkFzKQJBoIgJKgAtLxgECYQiohNLKlStl+vTpUr9+/T069a9//UvfvXD44YfL5MmTw+g0vUAAAS2gPsvVn7PPPlv69esnBx98cJnMJ598IqNGjZJZs2ZJ79695ZprrkENAQQQ4BqAMYBAIAKsAQIJJN1AAAEEENhFgOQSAwKBjAiQXMpIoGgmAnsR+OEPfyidOnWSwYMHl7vnHXfcoZNPb775Jp4IIBCQwJlnnqnvWvzf//3fcnv1i1/8Qr766iuZM2dOQD2nKwggYCvANYCtHOUQSJcAa4B0xYPWIIAAAgjEI0ByKR5HjoJA4gJcWCZOTAUIeBFo3bq1qC+Pb7rppnLrGzFihH4ny9/+9jcvbaISBBDwI9C8eXO59NJL5de//nW5FY4ePVoeeughWbJkiZ9GUQsCCKRagGuAVIeHxiEQWYA1QGQqdkQAAQQQyJAAyaUMBYum5luAC8t8x5/ehyPQpUsXWbt2rTz55JNSs2bNPTq2YcMGOffcc6VevXo6wcSGAALhCPz4xz+WE088Ue66665yO3XjjTfKyy+/LH/961/D6Tg9QQABawGuAazpKIhAqgRYA6QqHDQGAQQQQCAmAZJLMUFyGASSFuDCMmlhjo+AHwGVVFJfHjdq1EiuvPJKUY/Jq1OnjnzxxRf6MXgTJkyQVatW6S+fzznnHD+NohYEEPAioN6z9Oyzz8q4cePk5JNP3qPOF154Qb9r6YwzzhB1BxMbAgggwDUAYwCBMARYA4QRR3qBAAIIILCrAMklRgQCGRHgwjIjgaKZCEQQUF8sqz87duzYY+9KlSrpL5fVHzYEEAhL4IMPPpDzzz9fvv76azn22GPlmGOOkbp168rnn38ub731lrzxxhtywAEH6Hcyfe973wur8/QGAQSsBLgGsGKjEAKpE2ANkLqQ0CAEEEAAgRgESC7FgMghEPAhwIWlD2XqQMCfwPvvvy8zZ86UFStWiHoUnnpEXtOmTfXdSnyp7C8O1ISAb4Hly5fL0KFDC75TTd3JeMstt8iRRx7pu1nUhwACKRXgGiClgaFZCFgIsAawQKMIAggggECqBUgupTo8NA6B/xOYO3euzJs3T4YPHw4LAggggAACCGRc4OOPPxb1JdPOyeWGDRtmvFc0HwEE4hbgGiBuUY6HQPEFWAMUPwa0AAEEEEAgHgGSS/E4chQEjAW+/fZbeeWVV0TdvbBx48ayR2Cpv1dfNNWqVUsqV65sfFwKIIAAAggggAACCCCAQDoFuAZIZ1xoFQIIIIAAAggggIC5AMklczNKIOAsMGvWLLnttttk3bp1+p0r6h0ry5Yt08ddsmSJdOnSRd+h1LFjR+e6OAACCKRPYPv27fLss8/Ke++9J2vWrJEtW7bs0Uh1Xhg2bFj6Gk+LEEAAAQQQQMBKgGsAKzYKIYAAAggggAACCKRUgORSSgNDs8IVeOmll+SXv/ylHHLIIdKzZ0/9zoWnnnqqLLmkeq7euaL+/b777gsXgp4hkFOBf/zjH3LVVVeJehyGSi6Xt+2cdM4pFd1GIEgBddfyxIkT9fvWVHJ527Zte/RTzf+lS5cG2X86hUBeBbgGyGvk6TcC/yfAGoDRgAACCCAQmgDJpdAiSn9SL9C1a1f56KOPdELpgAMOkEIv6e3fv7+8+eab+h1LbAggEJZAt27d5I033pAePXpI+/bt5Tvf+U65j8Dk/SthxZ7eIDB16lS59dZbdWK5UaNGUrdu3XLn/+TJkwFDAIGABLgGCCiYdAUBCwHWABZoFEEAAQQQSL0AyaXUh4gGhibQunVr/bi7W265RXetUHJp1KhRMmnSJFm0aFFo3ac/COReoEWLFtK2bVu5++67c28BAAJ5Ezj11FNFvW/lj3/8ozRr1ixv3ae/CORagGuAXIefziMgrAEYBAgggAACIQqQXAoxqvQp1QLHHHOMdO7cWQYPHlxucmngwIHy3HPPyWuvvZbqvtA4BBAwFzj55JPlzDPPlEGDBpkXpgQCCGRaoGXLlvq9ir/5zW8y3Q8ajwAC5gJcA5ibUQKBkARYA4QUTfqCAAIIIFAqQHKJsYCAZ4ELL7xQ/vOf/8iTTz6pH4Wz+51LmzZt0l88N2nSRB588EHPraM6BBBIWuDOO++UBQsWyPTp06Vq1apJV8fxEUAgRQLqxyVHHnmkjBgxIkWtoikIIOBDgGsAH8rUgUB6BVgDpDc2tAwBBBBAwF6A5JK9HSURsBKYMWOG3HTTTdKpUyd999JDDz0k48aNk2XLlsm///1vGTJkiMydO1fuueceOf30063qoBACCKRXYPPmzdK7d29R/+3bt6/+onn//fdPb4NpGQIIxCbwzDPP6DXAI488IkcccURsx+VACCCQfgGuAdIfI1qIQJICrAGS1OXYCCCAAALFEiC5VCx56s21wG233SZTpkyRfffdVw444ACdVDr00EPl448/lq1bt0r37t15ZE6uRwidD11APfZywIABsmHDhnK7WqlSJVm6dGnoFPQPgdwJzJo1S4YNG6bfvdC0aVOpUaNGQQP1fkY2BBAIS4BrgLDiSW8QMBVgDWAqxv4IIIAAAmkXILmU9gjRvmAFnn/+eXn44YdlyZIlsn79ev3l0g9+8AO56KKLpF27dsH2m44hkHeBxx57TN+1uGPHDmnUqJHUrVtXPyKz0DZ58uS8c9F/BIIS2Lhxo75z6dlnn9XnALWpRPLOm/p79XfqjmY2BBAIT4BrgPBiSo8QiCLAGiCKEvsggAACCGRNgORS1iJGezMvsHDhQqlZs6YcddRRme8LHUAAAXOBM844Q9+xpN6p1qxZM/MDUAIBBDIroBJL6tFYRx99tH70rUouV6lSpWB/1ONz2RBAIBwBrgHCiSU9QcBGgDWAjRplEEAAAQTSLkByKe0Ron3BCXz/+9/XdyfdfPPNwfWNDiGAwN4FWrZsKRdccIG+e4kNAQTyJXDcccdJ48aN9aNxy7tjMV8i9BaB/AhwDZCfWNNTBAoJsAZgXCCAAAIIhChAcinEqNKnVAucfPLJ+rF3JJdSHSYah0BiAupuBHXH0vDhwxOrgwMjgEA6Bdq0aSM///nP9TvX2BBAIF8CXAPkK970FoHdBVgDMCYQQAABBEIUILkUYlTpU6oF1BfK8+fPl7/85S9SrVq1VLeVxiGAQPwCc+bMkUGDBsnUqVPliCOOiL8CjogAAqkVuPbaa2XdunUyceLE1LaRhiGAQDICXAMk48pREciKAGuArESKdiKAAAIImAiQXDLRYl8EYhD45ptvpHfv3vL111/r/6pHZNSuXTuGI3MIBBDIgoB638rs2bPllVdekfPOO0+aNm0qNWrUKNj0jh07ZqFLtBEBBCIKfPrpp9K1a1fp0KGDXH311VK1atWIJdkNAQSyLsA1QNYjSPsRcBNgDeDmR2kEEEAAgXQKkFxKZ1xoVcACRx11lO7djh07pFKlSuX2VP3b0qVLA5agawjkU0A9Ek/Nb3UOKN12PxeUnh+WLVuWTyR6jUCgAj169JD169fL8uXLdVK5UaNGBZPL6pzA3U2BDgK6lVsBrgFyG3o6joAWYA3AQEAAAQQQCFGA5FKIUaVPqRbo3r175PZNnjw58r7siAAC2RCYPn165Iaq9zOxIYBAOAIquRxlU8klkstRpNgHgewIcA2QnVjRUgSSEGANkIQqx0QAAQQQKLYAyaViR4D6EUAAAQQQQAABBBBAAAEEEEAAAQQQQAABBBBAAIEMCZT/TK4MdYKmIoAAAggggAACCCCAAAIIIIAAAggggAACCCCAAAII+BH4/wAbie1QH5vWlwAAAABJRU5ErkJggg==" width="1499.5555555555557">


