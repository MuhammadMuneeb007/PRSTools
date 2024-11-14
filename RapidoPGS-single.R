library(RapidoPGS)



#1. Argument one is the Directory Example: `SampleData1`
#2. Argument two is the GWASfile Example: `SampleData1/SampleData1.txt`
#3. Argument three is the built. Example: `hg38` or `hg19`
#4. Argument 4 is the traittype Example: `cc` for binary  or `quant` for continous.
#5. Argument 5 is the niumber of people for required for continous phenotypes.


# Read the GWAS data file
args <- commandArgs(trailingOnly = TRUE)
directory <- args[1] 
gwasfile <- args[2]
build <- args[3]
traittype <- args[4]
n <- args[5]
print(directory)
print(gwasfile)
print(build)


ds <- read.table(gwasfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ds <- as.data.frame(ds)
print(names(ds))
print(data.frame(Column = 1:length(names(ds)), Name = names(ds)))
print(head(ds))
print(class(ds))


# Ensure required columns are present
required_columns <- c("CHR", "BP", "SNPID", "REF", "ALT","N", "ALT_FREQ", "BETA", "SE", "P")
 
# Subset relevant columns and remove non-autosomal chromosomes
ds <- ds[, c("SNPID", "CHR", "BP", "REF", "ALT", "N","ALT_FREQ", "BETA", "SE", "P")]
print(nrow(ds))
# Remove non-autosomal chromosomes
ds <- ds[ds$CHR != "X", ]
print(nrow(ds))
# Convert CHR to numeric
ds$CHR <- as.numeric(ds$CHR)
print(nrow(ds))
# Order by CHR and BP
ds <- ds[order(ds$CHR, ds$BP), ]
print(nrow(ds)) 
# Remove rows with NA in BETA or ALT_FREQ
ds<- ds[!is.na(ds$BETA) & !is.na(ds$ALT_FREQ), ]
print(nrow(ds))

# Compute the polygenic score
full_PGS <- rapidopgs_single(ds, trait = traittype , build = build,N=as.numeric(n))

# Print the result
print(head(full_PGS))

# Save the result to a file
 
result <-paste(".",args[1],"RapidoPGS-single-gwas.txt",sep="//")

if (file.exists(result)) {
file.remove(result)
print(paste("File", result, "deleted."))
}
 
write.table(full_PGS, file = result, sep = "\t", row.names = FALSE, quote = FALSE)





