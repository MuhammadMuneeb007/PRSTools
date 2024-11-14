Step1-DownloadSampleDataandSetupEnvironment
===========================================

Download Anaconda
-----------------

To run the PRS Notebooks, set up a conda environment. Some tools require
a dedicated (separate) conda environment, and for such tools,
instructions are provided in their corresponding Notebooks.

Windows:
~~~~~~~~

1. Visit the `Anaconda website <https://www.anaconda.com/download>`__.
2. Choose the appropriate installer (64-bit is recommended for most
   systems).
3. Run the installer and follow the on-screen instructions.

Linux:
~~~~~~

1. Open a terminal window and run these commands:
   ``bash     wget https://repo.anaconda.com/archive/Anaconda3-2024.01-Linux-x86_64.sh     chmod +x Anaconda3-2024.01-Linux-x86_64.sh     bash Anaconda3-2024.01-Linux-x86_64.sh     source ~/.bashrc``

Create a Dedicated Environment for Your Project:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Open a terminal window (or Anaconda Prompt on Windows).

2. Create a new environment named “prstools” with Python 3.10:
   ``bash     conda create -n prstools python=3.10``

3. Activate the environment: ``bash     conda activate prstools``

Install R within the Environment:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install R version 4.3.2 from the conda-forge channel:

.. code:: bash

   conda install -c conda-forge r-base=4.3.2
   conda install -c conda-forge r-essentials

Data Preparation
----------------

In this section, we will download the sample data.

Download Sample Data
~~~~~~~~~~~~~~~~~~~~

The first step is to download the sample data. You can find the base
data
`here <https://drive.google.com/file/d/1RWjk49QNZj9zvJHc9X_wyZ51fdy6xQjv/view>`__
and the target data
`here <https://drive.google.com/file/d/1uhJR_3sn7RA8U5iYQbcmTp6vFdQiF4F2/view>`__.

Extracting Target Data
~~~~~~~~~~~~~~~~~~~~~~

Extract the target data and organize it with the following files in a
directory named ``SampleData1``:

-  ``SampleData1.bed``: Plink Files (bed, bim, fam) containing genotype
   information.
-  ``SampleData1.bim``
-  ``SampleData1.fam``
-  ``SampleData1.cov``: Contains covariate information.
-  ``SampleData1.height``: Contains the actual height phenotype.
-  ``SampleData1.gz``: GWAS summary statistic file.

**Ensure to rename all files to the ``SampleData1`` prefix.**

::

   .
   ├── SampleData1
   │   ├── SampleData1.bed
   │   ├── SampleData1.bim
   │   ├── SampleData1.cov
   │   ├── SampleData1.fam
   │   ├── SampleData1.gz
   │   └── SampleData1.height

View Bed Files
~~~~~~~~~~~~~~

To view Bed files, visit `Plink Binary Files
Documentation <https://zzz.bwh.harvard.edu/plink/binary.shtml>`__.

Analyzing Quantitative Trait (Height)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We will begin by analyzing the quantitative trait (height). The dataset
used for the continuous trait analysis is sourced from the tutorial by
`choishingwan <https://choishingwan.github.io/PRS-Tutorial/base/>`__,
and we express our gratitude for making it publicly available.

Polygenic Risk Score (PRS) Tools Overview
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PRS tools can be categorized into three types:

1. **Type 1: Summary Statistics File Only**

   -  These tools use a summary statistic file (e.g., gwas.txt) and
      estimate betas or posterior probabilities. They may require
      population-specific linkage disequilibrium calculations.

2. **Type 2: Summary Statistic File and Individual Data Set**

   -  This category of PRS tools requires both a summary statistic file
      and an individual dataset. They optimize various hyperparameters
      for beta estimation or snp weight estimation.

3. **Type 3: Individual Data Set Only**

   -  These tools focus on polygenic risk score calculation for
      individuals without relying on a summary statistic file. They
      incorporate linkage equilibrium for betas or snp weight
      estimation.

4. **Type 4: Multi-ancestry/Disease Tools**

   -  These tools use GWAS from multiple populations or multiple
      diseases to calculate posterior probabilities or new Betas.

Data Splitting and Hyperparameter Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The target data is split into a 75% training or hyperparameter
optimization set and a 25% testing set for final performance evaluation.

For binary phenotypes, we used a stratified cross-validation technique
to split the data.

.. code:: python

   skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

For continuous phenotypes, we used cross-validation to split the data
into 5 folds.

.. code:: python

   kf = KFold(n_splits=5, shuffle=True, random_state=42)

Hyperparameter optimization can be categorized into two types:

1. **Genotype data Quality controls related Hyperparameters**

   -  This includes considerations such as pruning and clumping during
      quality controls on the target data.

2. **Tools-Specific Hyperparameters**

   -  When using the GWAS File and an individual dataset, specific
      hyperparameters, such as lambdas for lasso, may be considered.

Full Credit to Shing Wan Choi
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  GitHub: `Shing Wan Choi <https://github.com/choishingwan>`__

Tutorial Link from where the code has been inferred: `Shing Wan Choi’s
PRS Tutorial <https://choishingwan.github.io/PRS-Tutorial/base/>`__

.. code:: ipython3

    # The following code is for Linux
    !conda config --add channels defaults
    !conda config --add channels bioconda
    !conda config --add channels conda-forge
    !conda install bioconda::bedtools
    !conda install bioconda::samtools
    # You can execute the code from terminal after activitating the conda environment.
    # conda activate genetics.
    

View Data
---------

In this section, we will view the sample data.

.. code:: ipython3

    import os
    import pandas as pd
    import subprocess
    
    # Set the directory where the files are located
    filedirec = "SampleData1"
    
    # Define file paths for different data files
    BED = filedirec + os.sep + filedirec
    BIM = filedirec + os.sep + filedirec+".bim"
    FAM = filedirec + os.sep + filedirec+".fam"
    COV = filedirec + os.sep + filedirec+".cov"
    Height = filedirec + os.sep + filedirec+".height"
    GWAS = filedirec + os.sep + filedirec+".gz"
    
    # Read the first 10 rows of the BED file (matrix, cannot be viewed directly)
    # This is a placeholder comment since viewing a BED file is not applicable in this context
    #https://en.wikipedia.org/wiki/BED_(file_format)
    
    # Read the first 10 rows of the BIM file (PLINK map file)
    temp = pd.read_csv(BIM, sep="\s+", header=None, nrows=10)
    print("Columns of BIM file:")
    print(temp.columns)
    print("First 10 rows of BIM file:")
    print(temp.head())
    
    # Read the first 10 rows of the FAM file (PLINK pedigree file)
    temp = pd.read_csv(FAM, sep="\s+", header=None, nrows=10)
    print("Columns of FAM file:")
    print(temp.columns)
    print("First 10 rows of FAM file:")
    print(temp.head())
    
    # Read the first 10 rows of the COV file (covariate information file)
    temp = pd.read_csv(COV, sep="\s+", nrows=10)
    print("Columns of COV file:")
    print(temp.columns)
    print("First 10 rows of COV file:")
    print(temp.head())
    
    # Read the first 10 rows of the Height file (file containing actual height phenotypes)
    temp = pd.read_csv(Height, sep="\s+", nrows=10)
    print("Columns of Height file:")
    print(temp.columns)
    print("First 10 rows of Height file:")
    print(temp.head())
    
    # Read the first 10 rows of the GWAS file (GWAS summary statistic file)
    temp = pd.read_csv(GWAS, sep="\s+", nrows=10)
    print("Columns of GWAS file:")
    print(temp.columns)
    print("First 10 rows of GWAS file:")
    print(temp.head())
    


.. parsed-literal::

    Columns of BIM file:
    Index([0, 1, 2, 3, 4, 5], dtype='int64')
    First 10 rows of BIM file:
       0           1         2       3  4  5
    0  1   rs3131962  0.490722  756604  A  G
    1  1  rs12562034  0.495714  768448  0  0
    2  1   rs4040617  0.500708  779322  G  A
    3  1  rs79373928  0.587220  801536  G  T
    4  1  rs11240779  0.620827  808631  G  A
    Columns of FAM file:
    Index([0, 1, 2, 3, 4, 5], dtype='int64')
    First 10 rows of FAM file:
             0        1  2  3  4  5
    0  HG00096  HG00096  0  0  1 -9
    1  HG00097  HG00097  0  0  2 -9
    2  HG00099  HG00099  0  0  2 -9
    3  HG00100  HG00100  0  0  2 -9
    4  HG00101  HG00101  0  0  1 -9
    Columns of COV file:
    Index(['FID', 'IID', 'Sex'], dtype='object')
    First 10 rows of COV file:
           FID      IID  Sex
    0  HG00096  HG00096    1
    1  HG00097  HG00097    2
    2  HG00099  HG00099    2
    3  HG00100  HG00100    2
    4  HG00101  HG00101    1
    Columns of Height file:
    Index(['FID', 'IID', 'Height'], dtype='object')
    First 10 rows of Height file:
           FID      IID      Height
    0  HG00096  HG00096  169.132169
    1  HG00097  HG00097  171.256259
    2  HG00099  HG00099  171.534380
    3  HG00101  HG00101  169.850176
    4  HG00102  HG00102  172.788361
    Columns of GWAS file:
    Index(['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'OR', 'INFO', 'MAF'], dtype='object')
    First 10 rows of GWAS file:
       CHR      BP         SNP A1 A2       N        SE         P        OR  \
    0    1  756604   rs3131962  A  G  388028  0.003017  0.483171  0.997887   
    1    1  768448  rs12562034  A  G  388028  0.003295  0.834808  1.000687   
    2    1  779322   rs4040617  G  A  388028  0.003033  0.428970  0.997604   
    3    1  801536  rs79373928  G  T  388028  0.008413  0.808999  1.002036   
    4    1  808631  rs11240779  G  A  388028  0.002428  0.590265  1.001308   
    
           INFO       MAF  
    0  0.890558  0.369390  
    1  0.895894  0.336846  
    2  0.897508  0.377368  
    3  0.908963  0.483212  
    4  0.893213  0.450410  
    

Important Note
--------------

1. This notebook needs to be executed only once to download the sample
   data.
2. If someone wants to use the data in a different format, ensure you
   follow the same directory structure.
3. Please ensure that your data has the same number of columns and the
   same headers. Cheers!
4. For continuous phenotype GWAS, the ``SampleData1/SampleData1.gz``
   file should have BETAs, and for binary phenotypes, it should have OR
   instead of BETAs. If BETAs are not available, we convert OR to BETAs
   using ``BETA = np.log(OR)`` and convert BETAs to OR using
   ``OR = np.exp(BETA)``.
5. ``import numpy as np``; ``np`` is the NumPy module.

