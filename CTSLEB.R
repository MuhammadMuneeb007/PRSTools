args <- commandArgs(trailingOnly = TRUE)
print(args)
library(data.table)
library(CTSLEB)
library(data.table)
library(dplyr)
library(devtools)
library(caret)
library(SuperLearner)
library(ranger)
library(glmnet)


# Argument Descriptions

# Argument Descriptions
# Argument one is the directory for gwas1. Example: `SampleData1`
# Argument two is the directory for gwas2. Example: `SampleData1`
# Argument three is the file name for EUR data. Example: `SampleData1\\SampleData1_CTSLEB.txt`
# Argument four is the file name for AFR data. Example: `SampleData1\\SampleData1_CTSLEB.txt`
# Argument five is Plink path
# Argument six is Plink2 path

# Argument seven is train genotype data for population 1
# Argument eight is train genotype data for population 2


# For our implementatkion it is the same as the train data, as we got the new weights and then used it for validation. 
# Argument nine is test genotype data for population 1
# Argument ten is test genotype data for population 2


# Argument eleven is outdirectory
# Argument twelve is outprefix









sum_target <- fread(paste0(args[3]),header=T)
sum_ref <- fread(paste0(args[4]),header=T)
head(sum_target)
head(sum_ref)

PRS_farm <- SetParamsFarm(plink19_exec = args[5],
                          plink2_exec = args[6])

print(PRS_farm)

prs_mat <- dimCT(results_dir = args[11],
                 sum_target = sum_target,
                 sum_ref = sum_ref,
                 ref_plink = args[7],
                 target_plink = args[8],
                 test_target_plink = args[10],
                 out_prefix = args[12],
                 params_farm = PRS_farm)


 

total_samples <- nrow(prs_mat)

# Calculate the number of samples for training (75% of total)
train_size <- floor(0.75 * total_samples)


 
# Assign the remaining 25% to the testing set
#prs_test <- prs_mat[(train_size + 1):total_samples, ]

# Assign the first 75% to the training set
prs_tune <- prs_mat[1:train_size, ]


n.total.prs <- length(pthres)^2*length(r2_vec)*length(wc_base_vec)
prs_r2_vec_test <- rep(0,n.total.prs)

y_tune <- fread(paste0(args[13]))

for(p_ind in 1:n.total.prs){
  model <- lm(y_tune$V1~prs_tune[,(2+p_ind)])
  prs_r2_vec_test[p_ind] <- summary(model)$r.square
}

max_ind <- which.max(prs_r2_vec_test)

 




best_snps <- colnames(prs_tune)[max_ind+2]
#calculate eb effect using EB coefficients
prs_mat_eb <- CalculateEBEffectSize(bfile = args[10] ,
                                    snp_ind = best_snps,
                                    plink_list = plink_list,
                                    out_prefix = args[12],
                                    results_dir = args[11],
                                    params_farm = PRS_farm)

 



prs_tune <- prs_mat_eb[1:train_size,]
prs_validation <- prs_mat_eb[(train_size + 1):total_samples, ]
 
 


y_tune <- fread(paste0(args[13]))

Cleaned_Data <- PRS_Clean(Tune_PRS = prs_tune,Tune_Y = y_tune,Validation_PRS = prs_validation)

y_vad_file <- paste0(args[14])
y_vad <- fread(y_vad_file)

prs_tune_sl <- Cleaned_Data$Cleaned_Tune_PRS
prs_valid_sl <- Cleaned_Data$Cleaned_Validation_PRS
 




SL.library <- c(
  "SL.glmnet",
  "SL.ridge"
)

sl <- SuperLearner(Y = y_tune$V1, 
                   X = prs_tune_sl[,-c(1,2)], 
                   family = gaussian(),
                   SL.library = SL.library)

y_pred_valid <- predict(sl, prs_valid_sl[,-c(1,2)], onlySL = TRUE)

 
model <- lm(y_vad$V1~y_pred_valid$pred)
r2_ctsleb <- summary(model)$r.square
print(r2_ctsleb)
y_pred_tune <- predict(sl, prs_tune_sl[,-c(1,2)], onlySL = TRUE)
Predicted_Tune_Y <- y_pred_tune$pred
Tune_PRS <- prs_tune_sl[,-c(1,2)]

Final_Betas <- ExtractFinalBetas(Tune_PRS = prs_tune_sl[,-c(1,2)],Predicted_Tune_Y = y_pred_tune$pred,prs_mat_eb = prs_mat_eb,unique_infor_post = unique_infor_post,pthres = pthres)


write.table(Final_Betas,file = paste0(args[11],"Final_PRS_Coefficients"),col.names = T,row.names = F,quote=F)
system(paste0(args[6]," --threads 2 --score ",args[11],"Final_PRS_Coefficients cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile ",args[10]," --out ",args[11],"Final_Betas_Evaluation"))

Final_Betas_Evaluation <- read.delim(paste0(args[11],"Final_Betas_Evaluation.sscore"), header=FALSE, comment.char="#")
cor(Final_Betas_Evaluation[,5],c(y_pred_tune$pred,y_pred_valid$pred))














