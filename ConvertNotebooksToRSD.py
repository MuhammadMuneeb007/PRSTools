import os
#os.system("jupyter nbconvert --to rst *.ipynb")
#os.system("make html")



import subprocess
import os

# List of Jupyter notebook filenames (you can modify this list as per your requirements)
notebooks = [
    "Introduction.ipynb",
    "Step1-DownloadSampleDataandSetupEnvironment.ipynb",
    "Step2-GWASAndIndividualGenotypeDataQualityControls.ipynb",
    "NotebooksWorkflow.ipynb",
    "Plink.ipynb",
    "PRSice-2.ipynb",
    "DBSLMM.ipynb",
    "GCTA.ipynb",
    "Lassosum.ipynb",
    "LDpred-2-Inf.ipynb",
    "LDpred-2-Grid.ipynb",
    "LDpred-2-Auto.ipynb",
    "LDpred-2-Lassosum.ipynb",
    "LDpred-funct-code.ipynb",
    "GCTB-SBayesR.ipynb",
    "GCTB-SBayesRC.ipynb",
    "LDAK-GenotypeData.ipynb",
    "LDAK-GWAS.ipynb",
    "PRScs.ipynb",
    "PRScsx.ipynb",
    "tlpSum.ipynb",
    "PRSbils.ipynb",
    "CTPR.ipynb",
    "NPS.ipynb",
    "SDPR.ipynb",
    "JAMPred.ipynb",
    "EBPRS.ipynb",
    "PANPRS.ipynb",
    "BOLT-LMM.ipynb",
    "RapidoPGS-single.ipynb",
    "LDpred-gibbs.ipynb",
    "LDpred-fast.ipynb",
    "LDpred-p+t.ipynb",
    "LDpred-inf.ipynb",
    "AnnoPredCode.ipynb",
    "PleioPredCode.ipynb",
    "smtpredwMtOLS.ipynb",
    "smtpredwMtSBLUP.ipynb",
    "C+T.ipynb",
    "viprs-simple.ipynb",
    "viprs-grid.ipynb",
    "HAIL.ipynb",
    "GEMMA-LM.ipynb",
    "GEMMA-LMM.ipynb",
    "GEMMA-BSLMM.ipynb",
    "MTG2.ipynb",
    "SCT.ipynb",
    "XPBLUP.ipynb",
    "CTSLEB.ipynb",
    "PolyPred.ipynb"
]

# Build the Jupyter book
for notebook in notebooks:
    result = subprocess.run(['jupyter-book', 'build', notebook], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error building {notebook}: {result.stderr}")
    else:
        print(f"Successfully built {notebook}")
