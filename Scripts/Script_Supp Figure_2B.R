# Load packages
library(gplots)
library(RColorBrewer)

# Set working directory
setwd("/Users/WIMM/Documents/BRCA_2018/Github/")

# Read file
df <- read.table("Data/TCGA_final.txt", header=TRUE, sep="\t", na.strings="NA", check.names=FALSE)

# CIN25 gene list
genes <- c("TPX2", "PRC1", "FOXM1", "CDC2", "C20orf24", "MCM2", "H2AFZ", "TOP2A", "PCNA", "UBE2C", "MELK", "TRIP13", "CNAP1", "MCM7", "RNASEH2A", "RAD51AP1", "KIF20A", "CDC45L", "MAD2L1", "ESPL1", "CCNB2", "FEN1", "TTK", "CCT5", "RFC4")
genes <- intersect(names(df), genes)

# Subset relevant columns
df <- df[,c("BRCA_status", genes)]

# Remove missing values
df <- na.omit(df)

# Compute FC and p-values
fc.brca1 <- NULL
fc.brca2 <- NULL
p.brca1 <- NULL
p.brca2 <- NULL

for(i in 1:length(genes)) {

    # Subset gene
    df.small <- df[, c("BRCA_status", genes[i])]
    
    # Retrieve expression values
    proficient <- df.small[which(df.small$BRCA_status=="BRCA-proficient"), 2]
    brca1 <- df.small[which(df.small$BRCA_status=="BRCA1-deficient"), 2]
    brca2 <- df.small[which(df.small$BRCA_status=="BRCA2-deficient"), 2]
    
    # Compute FC
    fc.brca1[i] <- log2(median((2^brca1 - 1), na.rm=TRUE) / median((2^proficient - 1), na.rm=TRUE))
    fc.brca2[i] <- log2(median((2^brca2 - 1), na.rm=TRUE) / median((2^proficient - 1), na.rm=TRUE))
       
    
    # Compute p.values
    p.brca1[i] <- wilcox.test(brca1, proficient, na.rm=TRUE)$p.value
    p.brca2[i] <- wilcox.test(brca2, proficient, na.rm=TRUE)$p.value
    
}

# Tabulate
data.frame(genes, fc.brca1, p.brca1, fc.brca2, p.brca2)

# Heatmap
  # Set factor levels
    df$BRCA_status <- as.factor(df$BRCA_status)
    df$BRCA_status <- factor(df$BRCA_status, levels=c("BRCA1-deficient", "BRCA2-deficient", "BRCA-proficient"))

    # Sort rows according to factor levels
    df <- df[order(df$BRCA_status),]
    
    # Transpose table
    n <- df$BRCA_status
    df_transposed <- as.data.frame(t(df[,-1]))
    colnames(df_transposed) <- n

    # Convert table to matrix
    df_transposed <- as.matrix(df_transposed)
    
    # Definition
    data <- df_transposed
    color.labels <- c(rep("purple", length(which(df$BRCA_status=="BRCA1-deficient"))), rep("pink", length(which(df$BRCA_status=="BRCA2-deficient"))), rep("yellow", length(which(df$BRCA_status=="BRCA-proficient"))))
    
    # Plot and save
    pdf(file="Figures/Supp Figure_2B.pdf")

    heatmap.2(df_transposed, dendrogram='none', Rowv=FALSE, Colv=FALSE, labCol=FALSE, trace='none', scale="row", cexRow=1, margins=c(12,8), key=FALSE, lwid=c(0.1,4), lhei=c(0.1,4), col=bluered, ColSideColors=color.labels)

    dev.off()
