library(lassosum)
library(penRegSum)
library(data.table)
library(pROC)
package_directory <- system.file(package = "lassosum")
args <- commandArgs(trailingOnly = TRUE)
print(args)

#1. Argument one is the directory. Example: `SampleData1`
#2. Argument two is the file name. Example: `SampleData1\\Fold_0`
#3. Argument three is the output file name. Example: `train_data`
#4. Argument four is the specific function to be called. Example: `train_data.QC.clumped.pruned`
#5. Argument five is the GWAS. Example: `EBPRSGWAS.txt`
# Argument Descriptions:
# 1. Directory
# 2. File name
# 3. Output file name
# 4. Specific function to be called
# 5. Number of PCA
# 6. convergence
# 7. iteration

#args <- c(
#  "SampleData1",
#  "SampleData1\\Fold_0",
#  "train_data",
#  "train_data.QC.clumped.pruned",
#  "6",
#  "convergence",
#  "iteration"
#)

#args <- c("SampleData1","SampleData1\\Fold_0","train_data","train_data.QC.clumped.pruned","EBPRSGWAS.txt"  )


ld.file <- "EUR.hg19.bed"
result <- paste(".", args[1], paste(args[1], toString(".txt"), sep = ""), sep = "//")
bimfile <- paste(".", args[2], paste(args[4], toString(".bim"), sep = ""), sep = "//")
bimfile <-fread(bimfile)

result <- paste(".", args[1], paste(args[1], toString(".txt"), sep = ""), sep = "//")
ss <- fread(result)

ss <- ss[ss$SNP %in% bimfile$V2, ]

ss <- ss[!P == 0]
cor <- p2cor(p = ss$P,
             n = ss$N,
             sign = ss$BETA)

lambdas <- c(0.001,0.01) 
taus <- c(0.01,0.1)
s <- c(0.5)
bfile <- paste(".", args[2], args[4], sep = "//")


out <-tlpSum(
  cors= cor,
  bfile = bfile,
  lambdas = lambdas,
  taus = taus,
  s =  s,
  #s = 0.5,
  thr = as.numeric(args[6]),
  #init = NULL,
  maxIter = as.numeric(args[7]),
  #extract = NULL,
  ldBlocks = ld.file,
  #corBim = NULL
)

 
selected_columns <- data.frame(
  lambdas = out$lambdas,
  taus = out$taus,
  s = out$s
)

#result <-paste(".",args[2],paste(args[3],toString(".lassosum_betas"), sep = ""),sep="//")
result <-paste(".",args[2],paste(args[3],toString(".tlpSum_parameters"), sep = ""),sep="//")
write.table(selected_columns, file = result, row.names = FALSE,sep=",", quote = FALSE)

result <-paste(".",args[2],paste(args[3],toString(".tlpSum_betas"), sep = ""),sep="//")
write.table(out$beta, file = result, row.names = FALSE,sep=",", quote = FALSE)
print("IterationDOne")
