# Load required libraries
library(lassosum)
library(data.table)
library(bigsnpr)
library(magrittr)
 
 
# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)

# Argument Descriptions:
# 1. Directory
# 2. File name
# 3. Output file name
# 4. Specific function to be called
# 5. Number of PCA
# 6. Use SNPs in GWAS that are also in HapMap
# 7. Reference panel name

#args <- c(
#  "SampleData1",
#  "SampleData1\\Fold_0",
#  "train_data",
#  "train_data.QC.clumped.pruned",
#  "6",
#  "hapmap",
#  "same"
#)












# Load phenotype, covariate, and principal component data
result <- paste(".", args[2], paste(args[3], toString(".PHENO"), sep = ""), sep = "//")
phenotype <- fread(result)

result <- paste(".", args[2], paste(args[3], toString(".cov"), sep = ""), sep = "//")
covariate <- fread(result)

result <- paste(".", args[2], paste(args[3], toString(".eigenvec"), sep = ""), sep = "//")
pcs <- fread(result)

# Rename columns
result <- paste(".", args[1], paste(args[1], toString(".txt"), sep = ""), sep = "//")
sum.stat <- result
ss <- fread(sum.stat)
 
# Load SNP information
info <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/25503788",
  fname = "map_hm3_ldpred2.rds"
))

# Filter SNPs based on HapMap
if (args[6] == "hapmap") {
  ss <- ss[ss$SNP %in% info$rsid, ]
  print("Filtered SNPs based on HapMap")
}

# Remove rows where P == 0
ss <- ss[!P == 0]

# Calculate correlation
cor <- p2cor(p = ss$P,
             n = ss$N,
             sign = ss$BETA)

# Define file paths
result <- paste(".", args[2], args[4], sep = "//")
bfile <- result
ld.file <- "EUR.hg19"
delta <- c(0.1, 0.2)
nlambda <- c(0.001, 0.1, 1)

# Create combinations for lassosum parameters
combinations <- expand.grid(nlambda, delta)
colnames(combinations) <- c("lambda", "delta")


# Write all lassosum betas results to file
result <- paste(".", args[2], paste(args[3], toString(".lassosum_parameters"), sep = ""), sep = "//")
write.table(combinations, file = result, row.names = FALSE, sep = ",", quote = FALSE)



# Run lassosum.pipeline based on the condition of args[7]
if (args[7] == "same") {
  out <- lassosum.pipeline(
    cor = cor,
    chr = ss$CHR,
    pos = ss$BP,
    A1 = ss$A1,
    A2 = ss$A2,
    lambda = nlambda,
    s = delta,
    ref.bfile = bfile,
    test.bfile = bfile,
    LDblocks = ld.file
  )
  print("Same referrence panel as input file.")
} else {
  out <- lassosum.pipeline(
    cor = cor,
    chr = ss$CHR,
    pos = ss$BP,
    A1 = ss$A1,
    A2 = ss$A2,
    lambda = nlambda,
    s = delta,
    ref.bfile = args[7],
    test.bfile = bfile,
    LDblocks = ld.file
  )
}

# Combine lassosum beta results into a data frame
m <- matrix(0, ncol = 0, nrow = length(out$beta[[1]][, 1]))
m <- data.frame(m)

# Get all the betas from the lassosum output.
# Save all the betas from the lassosum.
for (i in 1:length(delta)) {
  x <- out$beta[[i]]
  x <- as.data.frame(x)
  for (loop2 in 1:dim(x)[2]) {
    m[, paste("V", loop2, i, sep = "")] <- x[[paste("V", loop2, sep = "")]]
  }
}

# Write all lassosum betas results to file
result <- paste(".", args[2], paste(args[3], toString(".lassosum_betas"), sep = ""), sep = "//")
write.table(m, file = result, row.names = FALSE, sep = ",", quote = FALSE)

 
t <- out$sumstats
selected_rows <- ss[t$order, ]

newgwas <- selected_rows[, c("SNP", "A1", "BETA")]  

# Write lassosum GWAS results to file
result <- paste(".", args[2], paste(args[3], toString(".lassosum_gwas"), sep = ""), sep = "//")
write.table(newgwas, file = result, row.names = FALSE, sep = ",", quote = FALSE)
