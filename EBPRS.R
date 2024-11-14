args <- commandArgs(trailingOnly = TRUE)
print(args)
# Argument Descriptions

#1. Argument one is the directory. Example: `SampleData1`
#2. Argument two is the file name. Example: `SampleData1\\Fold_0`
#3. Argument three is the output file name. Example: `train_data`
#4. Argument four is the specific function to be called. Example: `train_data.QC.clumped.pruned`
#5. Argument five is the GWAS. Example: `EBPRSGWAS.txt`
#6. Argument six is cases
#7. Argument six is controls

#args <- c("SampleData1","SampleData1\\Fold_0","train_data","train_data.QC.clumped.pruned","EBPRSGWAS.txt"  )

 


library("EBPRS")
result <-paste(".",args[1],args[5],sep="//")

train <- fread(result)  
#N1 <-20791
#N0 <- 323124

result <-paste(".",args[2],args[4],sep="//")
test <- read_plink(result)

 
result2 <- EBPRS(train = train, test = test, N1 = as.integer(args[6]), N0 = as.integer(args[7]))

result <-paste(".",args[2],"newEBPRS.txt",sep="//")
result2$result
write.table(result2$result, file = result, row.names = FALSE,sep="\t", quote = FALSE)



