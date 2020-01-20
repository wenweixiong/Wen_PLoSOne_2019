setwd("/Users/WIMM/Documents/BRCA_2018/Github/")

df <- read.table("Data/WSI_final.txt", header=TRUE, sep="\t", na.strings="NA", check.names=FALSE)

df <- df[, !duplicated(colnames(df))]

write.table(df, "/Users/WIMM/Documents/BRCA_2018/Github/Data/WSI_final.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


df <- read.table("Data/TCGA_final.txt", header=TRUE, sep="\t", na.strings="NA", check.names=FALSE)

df <- df[, !duplicated(colnames(df))]

write.table(df, "/Users/WIMM/Documents/BRCA_2018/Github/Data/TCGA_final.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
