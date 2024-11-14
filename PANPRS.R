library(PANPRSnext)
args <- commandArgs(trailingOnly = TRUE) 
print(args)

#1. Argument one is the directory. Example: `SampleData1`
#2. Argument two is the file name. Example: `SampleData1\\Fold_0`
#3. Argument three is the output file name. Example: `train_data`
#4. Argument four is the specific function to be called. Example: `train_data.QC.clumped.pruned`
#5. Argument five is the GWAS. Example: `PANPRSGWAS
 

#6. Argument six is iteration. Example: n_iter
#7. Argument seven is z_scale. Example: 2.5
#8. Argument eight is len_lim_lambda. Example: 10
#9. Argument nine is sub_tuning. Example: 0.1
#10. Argument ten is len_lambda. Example: 100
#panprs_n_iter = 1000,
#panprs_z_scale = 1,
#panprs_len_lim_lambda = 20,
#panprs_sub_tuning = 20,
#panprs_len_lambda = 20



result <-paste(".",args[2],args[5],sep="//")
GWAS <- read.table(result, header=TRUE,sep = "\t")

result <-paste(".",args[2],"train_data_LDFile.ld",sep="//")
LDfile <- read.table(result, header = TRUE, sep = "")

summary_z <- data.frame(
  Zobs1 = GWAS$z_statistic
  
)

output <- gsPEN_R(
  summary_z = summary_z,
  n_vec = GWAS$N[1],
  debug_output=TRUE,
  plinkLD = LDfile,
  n_iter = as.numeric(args[6]),
  z_scale = as.numeric(args[7]),
  tau_factor = c(1/25,1),
  len_lim_lambda = as.numeric(args[8]),
  sub_tuning = as.numeric(args[9]),
  lim_lambda = c(0.5, 0.9),
  len_lambda = as.numeric(args[10]),
    sparse_beta = as.logical(args[11])
)

result <-paste(".",args[2],"PANPRSGWASBETAS",sep="//")


if (file.exists(result)) {
  file.remove(result)
  print(paste("File", result, "deleted."))
}
write.table(output$beta_matrix, file = result)



result <-paste(".",args[2],"PANPRSGWASARGUMENTS",sep="//")

if (file.exists(result)) {
  file.remove(result)
  print(paste("File", result, "deleted."))
}
write.table(output$all_tuning_matrix, file = result)




