import pandas as pd
import os
import numpy as np
import pandas as pd
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# Suppress the SettingWithCopyWarning
pd.options.mode.chained_assignment = None

phenotypes =  [
    "asthma",
    "asthma_19",
    "SampleData1",
    "body_mass_index_bmi",  
    
    "depression",
    "gastro_oesophageal_reflux_gord_gastric_reflux",
    "high_cholesterol",
    "hypothyroidism_myxoedema",
    "irritable_bowel_syndrome",
    "migraine", 
]

phenotypes = [
     
    "body_mass_index_bmi_13",
    "body_mass_index_bmi_14",
    "body_mass_index_bmi_15",
    "body_mass_index_bmi_16",
    "body_mass_index_bmi_17",
    "body_mass_index_bmi_19",
    "body_mass_index_bmi_22",
    "body_mass_index_bmi_23",
    "body_mass_index_bmi_24",
    "body_mass_index_bmi_25",
    "body_mass_index_bmi_29",
    "body_mass_index_bmi_2",
    "body_mass_index_bmi_30",
    "body_mass_index_bmi_31",
    "body_mass_index_bmi_33",
    "body_mass_index_bmi_34",
    "body_mass_index_bmi_35",
    "body_mass_index_bmi_36",
    "body_mass_index_bmi_37",
    "body_mass_index_bmi_41",
    "body_mass_index_bmi_4",
    "body_mass_index_bmi_9",
    "body_mass_index_bmi_0",
    "body_mass_index_bmi_11",
    "body_mass_index_bmi",
]
"""
   

"""


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
        "lasso_parameters_count",
        
        'gctb_ld_model',
        'BayesModel', 
   
        "ldaksubmodel", 
        "ldakmodel", 
        "ldakpower",    
        
        
        "PRScs_ref_dir" ,
        "PRScs_phi" ,
        "PRScs_va"  ,
        "PRScs_vb" ,
        "PRScs_number_of_iteration" ,
        "PRScs_burning_iteration" ,
        "PRScs_thining_iteration" ,                                  
                                      
        "tlpsum_lambda"  ,
        "tlpsum_tau" ,
        "tlpsum_s" ,
        "convergencethreshold",
        "iteration", 

        "PRSbils_ref_dir",
         

        "PRSbils_thres",
        "PRSbils_n_iter",
        "PRSbils_beta_std",

        "PRSbils_va",
        "PRSbils_vb",
        "PRSbils_vc",

        "PRSbils_vd",
        "PRSbils_ve",
        "PRSbils_vf",

        "PRSbils_flip",
        "PRSbils_fixtau",   

        "CTPR_f",
        "CTPR_lamda",
        "CTPR_penalty",
        "CTPR_withoutsum", 

        "window_shift",



        "SDPR_M_val",
        "SDPR_opt_llk_val",
        "SDPR_iter_val_val",
        "SDPR_burn_val",
        "SDPR_thin_val",
        "SDPR_r2_val",
        "SDPR_a_val",
        "SDPR_c_val",
        "SDPR_a0k_val",
        "SDPR_b0k_val",
        "SDPR_referenceset",   

        "JAMPred_Iteration",
        "JAMPred_Beta_Lambda",
        "JAMPred_Beta_Binomial_Prior",
        "JAMPred_Min_Effect",
        "JAMPred_Max_Effect",
        "JAMPred_Min_Var",
        "JAMPred_Max_Var",

        "PlinkLDtype",
        "panprs_n_iter",
        "panprs_z_scale",
        "panprs_len_lim_lambda",
        "panprs_sub_tuning",
        "panprs_len_lambda",
        "panprs_sparse_beta",
        "panprs_parameters_count",
          
        'BOLTmodel',

        'genomic_built',
        
        "ldradius",
        "ldfilename",
        "colname",
 

        "gibsfraction",
        "gibsburn",
        "gibsiterations",

        "clumpingandpruningtime",


        "Tier",
        "pvalue_AnnoPred",
        "datafile",
        
        "hyp_search",
        "method",
        "epsilon_steps",
        "pi_steps",

        "gemmamodel",
        "relatedmatrixname",
        "lmmmodel",  

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
    print(common_rows.head().to_markdown())    
    for i in range(1, len(allfoldsframe_dropped)):
        # Get the next DataFrame
        next_df = allfoldsframe_dropped[i]

        # Count unique rows in the current DataFrame and the next DataFrame
        unique_in_common = common_rows.shape[0]
        unique_in_next = next_df.shape[0]

        print(next_df.head().to_markdown())
        #continue
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



method_results = [
    "Plink",
    "PRSice",
    "GCTA",
    "lassosum",
    "DBSLMM",
    "ldpred2_inf",
    "ldpred2_grid",
    "ldpred2_auto",
    "ldpred2_lassosum2",
    "LDpred-funct",
    "SBayesR",
    "SBayesRC",
    "ldak-genotype",
    "ldak-gwas",
    "PRScs",
    "PRScsx",
    "PRSbils",
    #"ctpr",
    #"NPS",

    #'PleioPred'
    "tlpSum",
    "SDPR",
    "JAMPred",
    "EBPRS",
    "PANPRS",
    "Bolt-lmm",
    "AnnoPred",
    "RapidoPGS-single",
    "RapidoPGS-multi",
    "LDpred-gibbs",
    "LDpred-p+t",
    "LDpred-inf",
    "LDpred-fast",
    "smtpred-wMtOLS",
    "smtpred-wMtSBLUP",
    "viprs_simple",
    "viprs_grid",
    "CandT",
    "HAIL",
    "GEMMA_LM",
    "GEMMA_LMM",
    "GEMMA_BSLMM",
    "MTG2",
    "sct",
    "XP-BLUP",
    "CTSLEB",
    
 
]
"""

    "NPS",

    "SDPR",
    "JAMPred",
    "EBPRS",  
    "PANPRS",
    "BOLT-LMM",

    "RapidoPGS-single",
    "RapidoPGS-multi",
    "LDpred-gibbs",
    "LDpred-p+t",
    "LDpred-inf",
    "LDpred-fast",

    "AnnoPred",
    "smtpred-wMtOLS",
    "smtpred-wMtSBLUP",
    "C+T",

    "HAIL",



    "Plink",
"PRSice-2",
"GCTA",
"DBSLMM",
"lassosum",
"ldpred2_inf",
"ldpred2_grid",
"ldpred2_auto",  
"ldpred2_lassosum2",  
"LDpred-funct",

"GCTB-SBayesR",
"GCTB-SBayesRC",
"LDAK-GWAS",
"LDAK-GenotypeData",
"PRScs",
"PRScsx",
"tlpSum",
"PRSbils",
"ctpr",
"NPS",

"SDPR",
"JAMPred",
"EBPRS",  
"PANPRS",
"BOLT-LMM",
"RapidoPGS-single",
"LDpred-gibbs",
"LDpred-p+t",
"LDpred-inf",
"LDpred-fast",

"AnnoPred",
"smtpred-wMtOLS",
"smtpred-wMtSBLUP",
"CandT",   
"viprs_simple",
"viprs_grid",
"HAIL",
"GEMMA-LM",
"GEMMA-LMM",
"GEMMA-BSLMM",

"lassosum",
"ldpred2_inf",
"ldpred2_grid",
"ldpred2_auto",  
"ldpred2_lassosum2",  
"LDpred-funct",

"GCTB-SBayesR",
"GCTB-SBayesRC",
"LDAK-GWAS",
"LDAK-GenotypeData",
"PRScs",
"PRScsx",
"tlpSum",
"PRSbils",
"ctpr",
"NPS",

"SDPR",
"JAMPred",
"EBPRS",  
"PANPRS",
"BOLT-LMM",
"RapidoPGS-single",
"LDpred-gibbs",
"LDpred-p+t",
"LDpred-inf",
"LDpred-fast",

"AnnoPred",
"smtpred-wMtOLS",
"smtpred-wMtSBLUP",
"CandT",   
"viprs_simple",
"viprs_grid",
"HAIL",
"GEMMA-LM",
"GEMMA-LMM",
"GEMMA-BSLMM",


"MTG2",
"sct",
"CTSLEB",
"XPBLUP",
"PleioPred",
"PolyPred",








"GCTB-SBayesR",
"GCTB-SBayesRC",
"LDAK-GWAS",
"LDAK-GenotypeData",
"PRScs",
"PRScsx",
"tlpSum",
"PRSbils",
"ctpr",
"NPS",

"SDPR",
"JAMPred",
"EBPRS",  
"PANPRS",
"BOLT-LMM",
"RapidoPGS-single",
"LDpred-gibbs",
"LDpred-p+t",
"LDpred-inf",
"LDpred-fast",

"AnnoPred",
"smtpred-wMtOLS",
"smtpred-wMtSBLUP",
"CandT",   
"viprs_simple",
"viprs_grid",
"HAIL",
"GEMMA-LM",
"GEMMA-LMM",
"GEMMA-BSLMM",
"MTG2",
"sct",
"CTSLEB",
"XPBLUP",
"PleioPred",
"PolyPred",


"""  




method_results = [
   "Plink",
"PRSice-2",
"GCTA",
"DBSLMM",
"lassosum",
"ldpred2_inf",
"ldpred2_grid",
#"ldpred2_auto",  
"ldpred2_lassosum2",  
"LDpred-funct",

   
]


for method1results in method_results:
    for phenotype in phenotypes:
        allframes = []
        countnumberofframes = 0
        for fold in range(0,5):
            try:
                data = pd.read_csv(phenotype+os.sep+"Fold_"+str(fold)+os.sep+method1results+os.sep+"Results.csv")
                print(len(data),method1results,phenotype)
                if method1results=="GCTA":
                    if 'lambda' in data.columns:
                        data.rename(columns={'lambda': 'heritability_lambda'}, inplace=True)
                        data.rename(columns={'model': 'temp_lambda'}, inplace=True)
                if method1results=="ldpred2_lassosum2":
                    if 'lambda' in data.columns:
                        data.rename(columns={'lambda': 'lasso_lambda'}, inplace=True)
                        data.rename(columns={'delta': 'lasso_delta'}, inplace=True)
                        data.rename(columns={'delta': 'lasso_delta'}, inplace=True)
 
                try:
                    data['pvalue'] = data['pvalue'].str.replace("Pt_", "")
                except:
                    pass
                data = data.replace({True: 1, False: 0})
                #data['pvalue'] = data['pvalue'].astype(float)
                allframes.append(data)
                #print("YES!",method1results,phenotype+os.sep+"Fold_"+str(fold))
                countnumberofframes = countnumberofframes+1
            except:
                 
                pass
                #print("NO!",method1results,phenotype+os.sep+"Fold_"+str(fold))
                #exit(0)
        if countnumberofframes>0:
            extracted_common_rows_list = find_common_rows(allframes)
            averaged_df = sum_and_average_columns(extracted_common_rows_list)
            averaged_df["Method"] = method1results
            averaged_df["Phenotype"] = phenotype
            averaged_df.to_csv(phenotype+os.sep+method1results+"_heritability.csv")
        else:
            averaged_df = pd.DataFrame()
            averaged_df["Method"] = method1results
            averaged_df["Phenotype"] = phenotype
            averaged_df.to_csv(phenotype+os.sep+method1results+"_heritability.csv")            
            print(method1results,phenotype)
            #exit(0)
        
    

allframes = []
for method1results in method_results:
    for phenotype in phenotypes:
        allframes.append(pd.read_csv(phenotype+os.sep+method1results+"_heritability.csv"))
        pass

print(len(allframes))
combined_df = pd.concat(allframes, ignore_index=True)
print(len(combined_df))
del combined_df["Unnamed: 0"]

df = combined_df.copy()



def basedonbesttrainingmethod1(df):
    sorted_df = df.sort_values(by=['Train_best_model'], ascending=[False])
    sorted_df = sorted_df.drop_duplicates(subset=['Phenotype'],keep='first')

    sorted_df = sorted_df.sort_values(by=['Phenotype'], ascending=[True])
    df = sorted_df.copy()
    df = sorted_df[["Phenotype","Train_best_model","Test_best_model"]]

    #df['Phenotype'] = df['Phenotype'].str.replace(r'_\d+', '', regex=True)
    #sorted_df = df.sort_values(by=['Train_best_model'], ascending=[False])
    #sorted_df = sorted_df.drop_duplicates(subset=['Phenotype'],keep='first')



    
    print(sorted_df[["Phenotype","Train_best_model","Test_best_model"]].to_markdown())
    return sorted_df[["Phenotype","Train_best_model","Test_best_model"]]
def basedonsumanddifferencemethod1(df):
    # Extracting unique phenotypes and assigning colors
    df['Difference'] = abs(df['Train_best_model'] - df['Test_best_model'])
    df['Sum'] = df['Train_best_model'] + df['Test_best_model']

    # Step 2: Sort the DataFrame
    sorted_df = df.sort_values(by=['Sum', 'Difference'], ascending=[False,True])
    sorted_df = sorted_df.drop_duplicates(subset=['Phenotype'],keep='first')
    
    #sorted_df.to_csv("sorted_data.csv", index=False)
    sorted_df = sorted_df.sort_values(by=['Phenotype'], ascending=[True])
    df = sorted_df.copy()
    #del sorted_df["Method"]
    #sorted_df = sorted_df.rename(columns={'Phenotype': 'Phenotype', 'Train_best_model': 'GCTA_Train', 'Test_best_model': 'GCTA_Test','Tool':'Method'})

    df = sorted_df[["Phenotype","Train_best_model","Test_best_model"]]
 
    #print(sorted_df[["Phenotype","Train_best_model","Test_best_model"]].to_markdown())

    df = sorted_df[["Phenotype","Train_best_model","Test_best_model"]]
    return df
 
#basedonbesttrainingmethod1(df)

#basedonsumanddifferencemethod(df)
results = {}

# Iterate over each unique combination of Method and Phenotype
for method in df["Method"].unique():
    for phenotype in df["Phenotype"].unique():
        # Filter the DataFrame for the current combination of Method and Phenotype
        t = df[(df["Method"] == method) & (df["Phenotype"] == phenotype)]
        
        if not t.empty:
            # Perform the desired operation (replace basedonsumanddifferencemethod1 with your actual function)
            x = basedonsumanddifferencemethod1(t)["Test_best_model"].values[0]
           
            
            # Store the result in the dictionary
            if method not in results:
                results[method] = {}
            results[method][phenotype] = x

# Convert the dictionary to a DataFrame
results_df = pd.DataFrame(results).T

# Print the resulting DataFrame
print(results_df.T.to_markdown())
print(results_df.to_markdown())

df = results_df
# Replace NaN and negative values with 0
df = df.fillna(0)  # Replace NaN with 0
df[df < 0] = 0    # Replace negative values with 0

# Create a single heatmap for the entire DataFrame with min 0 and max 1
plt.figure(figsize=(20, 20))
sns.heatmap(df, annot=True, cmap='viridis', cbar=True, vmin=0, vmax=1)
plt.title('Heatmap of All Variables (Min: 0, Max: 1)')
plt.xlabel('Variables')
plt.ylabel('Models')
 
plt.savefig("Analysis1-BetaAnalysis/Results.png")
