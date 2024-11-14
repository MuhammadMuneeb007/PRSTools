library(R2BGLiMS)
#data("JAMPred_Example")
#help(JAMPred)
#1. Argument one is the directory. Example: `SampleData1`
#2. Argument two is the file name. Example: `SampleData1\\Fold_0`
#3. Argument three is the output file name. Example: `train_data`
#4. Argument four is the specific function to be called. Example: `train_data.QC.clumped.pruned`
#5. Argument five is the GWAS. Example: `EBPRSGWAS.txt`
args <- commandArgs(trailingOnly = TRUE)
#result <-paste(".",args[1],args[4],sep="//")
#data_without_headers <- read.table("SampleData1/Fold_0/train_data.QC.clumped.pruned.JAMPred.JAMPred", header = FALSE, sep = "\t")
data_without_headers <- read.table(args[4], header = FALSE, sep = "\t")
#data_without_headers <- read.table(result, header = FALSE, sep = "\t")

#data_without_headers <- read.csv("SampleData4\Fold_0\train_data.QC.clump.prune.raw", header = FALSE)
#data_without_headers <- read.table("SampleData4/Fold_0/train_data.QC.clump.prune.JAMPred", header = FALSE, sep = "\t")


your_data_frame <- as.data.frame(lapply(data_without_headers, as.numeric))
your_data_frame <- your_data_frame[-1, ]
your_data_frame <- replace(your_data_frame, is.na(your_data_frame), 0)


#your_data_frame <- as.numeric(your_data_frame)
#data_without_headers <- read.table("a1.dose", header = FALSE, sep = "\t")

# Get the number of columns in the data
num_columns <- ncol(your_data_frame)
num_columns <- ncol(data_without_headers)

# Generate headers based on the number of columns
column_names <- paste("SNP", 1:num_columns, sep = "")
#your_data_frame[,"SNP1"]

# Assign headers to the data
colnames(your_data_frame) <- column_names
colnames(data_without_headers) <- column_names


#read GWAS
result <-paste(".",args[2],paste(args[1], toString(".JAMPred"), sep = ""),sep="//")
#result <-paste(".","SampleData1",paste("SampleData1", toString(".JAMPred"), sep = ""),sep="//")
#print(result)
#exit(0)

GWAS <- read.table(result, header=TRUE,sep = "\t")
print(GWAS)
#phenofile <- read.table("SampleData1/Fold_0/train_data.QC.clumped.pruned.fam", header = FALSE, sep = "")
# Replace 0 with 1 in specific column
#phenofile$V6[phenofile$V6 == 1] <- 0

# Replace 2 with 1 in specific column
#phenofile$V6[phenofile$V6 == 2] <- 1

# Count the number of 0 and 1 in the specific column
#controls <- sum(phenofile$V6 == 0)
#cases <- sum(phenofile$V6 == 1)
#gwas_cases <- 20791 
#gwas_controls <-  323124
#gwas_numbers <- gwas_cases+gwas_controls

betas <-GWAS$BETA
#betas<- exp(betas)

num_columns <- length(betas)

# Generate headers based on the number of columns
column_names <- paste("SNP", 1:num_columns, sep = "")

# Assign headers to the data
names(betas) <- column_names
Se <-GWAS$SE

num_columns <- length(Se)
#num_columns
# Generate headers based on the number of columns
column_names <- paste("SNP", 1:num_columns, sep = "")
#column_names
#snps <- chromosome.snps[[1]]
#snps
#data.validation[,snps]

your_data_frame[,column_names]
print("xx1")
data.validation <- your_data_frame
print("xx2")
your_data_frame <- as.matrix(your_data_frame)
print("xx3")
 


# Assign headers to the data
names(Se) <- column_names
#your_data_frame[,column_names]
#snps <- data_without_headers
n <-GWAS$N[1]
 


iterationsinmillions <- c( 0.01)
betaslambda <- c(0.1)
betabinomialprior <- c(1)
effectsduniformpriormin <- c(0.05,0.01)
effectsduniformpriormax <- c(2,3)
residualvarinvgammapriormin <- c(0.01)
residualvarinvgammapriormax <- c(0.01)
param_combinations <- expand.grid(
  Iteration = iterationsinmillions,
  Beta_Lambda = betaslambda,
  Beta_Binomial_Prior = betabinomialprior,
  Min_Effect = effectsduniformpriormin,
  Max_Effect = effectsduniformpriormax,
  Min_Var = residualvarinvgammapriormin,
  Max_Var = residualvarinvgammapriormax
)

result <-paste(".",args[2],paste(args[1], toString(".JAMPred_Parameters"), sep = ""),sep="//")
write.table(param_combinations, file = result, row.names = FALSE)


myDataFrame <- data.frame(matrix(runif(num_columns * nrow(param_combinations)), nrow = num_columns, ncol = nrow(param_combinations)))

numColumns <- nrow(param_combinations)
columnNames <- paste("Column", 1:numColumns, sep="_")


for (i in 1:nrow(param_combinations)) {
  row_values <- param_combinations[i, ]
  
  iter <- row_values$Iteration
  beta <- row_values$Beta_Lambda
  prior <- row_values$Beta_Binomial_Prior
  min_effect <- row_values$Min_Effect
  max_effect <- row_values$Max_Effect
  min_var <- row_values$Min_Var
  max_var <- row_values$Max_Var
  
  
  jampred.res.lin <- JAMPred(
    marginal.betas = betas,
    n.training =n,
    ref.geno = your_data_frame[,column_names],
    total.snps.genome.wide = num_columns, # Total SNPs across all chromosomes
    
    n.mil = iter,
    beta.binom.b.lambda=beta,
    beta.binom.a = prior,
    effect.sd.uniform.prior = c(min_effect,max_effect),
    residual.var.invgamma.prior = c(min_var, max_var),
    seed = 1 # For re-producibility. If not set a random seed is used
  )
  
  # Print values for each iteration
  print(paste(
    "Iteration:", iter,
    "Beta Lambda:", beta,
    "Beta Binomial Prior:", prior,
    "Min Effect:", min_effect,
    "Max Effect:", max_effect,
    "Min Var:", min_var,
    "Max Var:", max_var
  ))
  #myDataFrame[[i]] <- jampred.res.lin$step2.posterior.mean.snp.weights
  myDataFrame[[paste("X", toString(i), sep = "")]] <- jampred.res.lin$step2.posterior.mean.snp.weights

}
result <-paste(".",args[2],paste(args[1], toString(".JAMPred_Effects"), sep = ""),sep="//")

write.csv(myDataFrame, file = result, row.names = FALSE)


