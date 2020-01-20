# Set working directory
setwd("/Users/WIMM/Documents/BRCA_2018/Github/")

#####################################################################
######################### WSI (BRCA1) ###############################
#####################################################################

# Read file
df <- read.table("Data/WSI_final.txt", header=TRUE, sep="\t", na.strings="NA", stringsAsFactors=FALSE, check.names=FALSE)

# Subset relevant genotypes
df <- df[which(df$BRCA_status %in% c("BRCA1-deficient", "BRCA-proficient")), ]

# Subset samples with transcriptomic data available
df <- df[which(!is.na(df$T.Cell.Inflammed.Signature.Score)),]

# Define variables
    # Categorical
    variables.cat <- c("Tumour.Grade", "TNBC", "PAM50")

    # Continuous
    variables.cont <- "Age"

# Set factor levels
    # Grade
    df$Tumour.Grade <- factor(df$Tumour.Grade, levels=c("I", "II", "III"))
    
    # TNBC
    df$TNBC <- factor(df$TNBC, levels=c("non-TNBC", "TNBC"), labels=c("No", "Yes"))

    # PAM50
    df$PAM50 <- factor(df$PAM50, levels=c("LumA", "LumB", "Basal", "Her2", "Normal"), labels=c("Luminal A", "Luminal B", "Basal-like", "HER2-enriched", "Normal-like"))

# Split data frame
df.brca <- df[which(df$BRCA_status=="BRCA1-deficient"), ]
df.proficient <- df[which(df$BRCA_status=="BRCA-proficient"), ]
   
# Build model for categorical variables
results.list <- list()

for(i in 1:length(variables.cat)) {
    
    # Subset data frames
    df.brca.small <- df.brca[, variables.cat[i], drop=FALSE]
    df.brca.small$Cohort <- "BRCA1-deficient"
    df.proficient.small <- df.proficient[, variables.cat[i], drop=FALSE]
    df.proficient.small$Cohort <- "BRCA-proficient"
    
    df.small <- rbind.data.frame(df.brca.small, df.proficient.small)
    
    # Create factor level for cohort
    df.small$Cohort <- factor(df.small$Cohort, level=c("BRCA-proficient", "BRCA1-deficient"))
    
    # Create frequency table
    tbl.freq <- table(df.small[,1], df.small[,2])
    
    # Create proprotion table
    tbl.prop <- round(prop.table(table(df.small[,1], df.small[,2]), 2)*100, digits=1)
    
    # Fisher exact test
    p.value <- rep(fisher.test(as.matrix(tbl.freq))$p.value, nrow(tbl.freq))
    
    # Tabulate
    results <- cbind.data.frame(as.data.frame.matrix(tbl.freq), as.data.frame.matrix(tbl.prop), p.value)
    
    # Create level column
    results <- data.frame("level"=row.names(results), results, stringsAsFactors=FALSE)
    
    # Create variable column
    results <- data.frame("variable"=variables.cat[i], results, stringsAsFactors=FALSE)
    
    # Save into list
    results.list[[i]] <- results

    
}

results <- do.call(rbind.data.frame, results.list)

results.cat <- as.data.frame(apply(results, 2, function(x) as.character(as.factor(x))))

# Build model for continuous variable
    # Compute median (BRCA)
    middle <- fivenum(df.brca$Age)[3]
    q1 <- fivenum(df.brca$Age)[1]
    q3 <- fivenum(df.brca$Age)[4]
    cont.brca <- paste(middle, " (", q1, "-", q3, ")", sep="")

    # Compute median (WT)
    middle <- fivenum(df.proficient$Age)[3]
    q1 <- fivenum(df.proficient$Age)[1]
    q3 <- fivenum(df.proficient$Age)[4]
    cont.proficient <- paste(middle, " (", q1, "-", q3, ")", sep="")

    # Compute p.value
    p.value <- wilcox.test(df.brca$Age, df.proficient$Age, na.rm=TRUE)$p.value

    # Tabulate (match categorical variable data frame)
    results.cont <- data.frame("variable"="Age", "level"="-", "BRCA.proficient"=cont.proficient, "BRCA1.deficient"=cont.brca, "BRCA.proficient.1"="-", "BRCA1.deficient.1"="-", p.value)

# Merge
    # Merge
    results <- rbind.data.frame(results.cont, results.cat)

    # Rename columns
    names(results) <- c("variable", "level", "WSI_BRCA.proficient.freq", "WSI_BRCA1.deficient.freq", "WSI_BRCA.proficient.prop", "WSI_BRCA1.deficient.prop", "p.value")

    # Save as new object
    results.WSI.BRCA1 <- results

#####################################################################
######################### WSI (BRCA2) ###############################
#####################################################################

# Read file
df <- read.table("Data/WSI_final.txt", header=TRUE, sep="\t", na.strings="NA", stringsAsFactors=FALSE, check.names=FALSE)

# Subset relevant genotypes
df <- df[which(df$BRCA_status %in% c("BRCA2-deficient", "BRCA-proficient")), ]

# Subset samples with transcriptomic data available
df <- df[which(!is.na(df$T.Cell.Inflammed.Signature.Score)),]

# Define variables
    # Categorical
    variables.cat <- c("Tumour.Grade", "TNBC", "PAM50")

    # Continuous
    variables.cont <- "Age"

# Set factor levels
    # Grade
    df$Tumour.Grade <- factor(df$Tumour.Grade, levels=c("I", "II", "III"))
    
    # TNBC
    df$TNBC <- factor(df$TNBC, levels=c("non-TNBC", "TNBC"), labels=c("No", "Yes"))

    # PAM50
    df$PAM50 <- factor(df$PAM50, levels=c("LumA", "LumB", "Basal", "Her2", "Normal"), labels=c("Luminal A", "Luminal B", "Basal-like", "HER2-enriched", "Normal-like"))

# Split data frame
df.brca <- df[which(df$BRCA_status=="BRCA2-deficient"), ]
df.proficient <- df[which(df$BRCA_status=="BRCA-proficient"), ]
   
# Build model for categorical variables
results.list <- list()

for(i in 1:length(variables.cat)) {
    
    # Subset data frames
    df.brca.small <- df.brca[, variables.cat[i], drop=FALSE]
    df.brca.small$Cohort <- "BRCA2-deficient"
    df.proficient.small <- df.proficient[, variables.cat[i], drop=FALSE]
    df.proficient.small$Cohort <- "BRCA-proficient"
    
    df.small <- rbind.data.frame(df.brca.small, df.proficient.small)
    
    # Create factor level for cohort
    df.small$Cohort <- factor(df.small$Cohort, level=c("BRCA-proficient", "BRCA2-deficient"))
    
    # Create frequency table
    tbl.freq <- table(df.small[,1], df.small[,2])
    
    # Create proprotion table
    tbl.prop <- round(prop.table(table(df.small[,1], df.small[,2]), 2)*100, digits=1)
    
    # Fisher exact test
    p.value <- rep(fisher.test(as.matrix(tbl.freq))$p.value, nrow(tbl.freq))
    
    # Tabulate
    results <- cbind.data.frame(as.data.frame.matrix(tbl.freq), as.data.frame.matrix(tbl.prop), p.value)
    
    # Create level column
    results <- data.frame("level"=row.names(results), results, stringsAsFactors=FALSE)
    
    # Create variable column
    results <- data.frame("variable"=variables.cat[i], results, stringsAsFactors=FALSE)
    
    # Save into list
    results.list[[i]] <- results

    
}

results <- do.call(rbind.data.frame, results.list)

results.cat <- as.data.frame(apply(results, 2, function(x) as.character(as.factor(x))))

# Build model for continuous variable
    # Compute median (BRCA)
    middle <- fivenum(df.brca$Age)[3]
    q1 <- fivenum(df.brca$Age)[1]
    q3 <- fivenum(df.brca$Age)[4]
    cont.brca <- paste(middle, " (", q1, "-", q3, ")", sep="")

    # Compute median (WT)
    middle <- fivenum(df.proficient$Age)[3]
    q1 <- fivenum(df.proficient$Age)[1]
    q3 <- fivenum(df.proficient$Age)[4]
    cont.proficient <- paste(middle, " (", q1, "-", q3, ")", sep="")

    # Compute p.value
    p.value <- wilcox.test(df.brca$Age, df.proficient$Age, na.rm=TRUE)$p.value

    # Tabulate (match categorical variable data frame)
    results.cont <- data.frame("variable"="Age", "level"="-", "BRCA.proficient"=cont.proficient, "BRCA2.deficient"=cont.brca, "BRCA.proficient.1"="-", "BRCA2.deficient.1"="-", p.value)

# Merge
    # Merge
    results <- rbind.data.frame(results.cont, results.cat)

    # Rename columns
    names(results) <- c("variable", "level", "WSI_BRCA.proficient.freq", "WSI_BRCA2.deficient.freq", "WSI_BRCA.proficient.prop", "WSI_BRCA2.deficient.prop", "p.value")

    # Save as new object
    results.WSI.BRCA2 <- results

#####################################################################
######################### TCGA (BRCA1) ##############################
#####################################################################

# Read file
df <- read.table("Data/TCGA_final.txt", header=TRUE, sep="\t", na.strings="NA", stringsAsFactors=FALSE, check.names=FALSE)

# Subset relevant genotypes
df <- df[which(df$BRCA_status %in% c("BRCA1-deficient", "BRCA-proficient")), ]

# Subset samples with transcriptomic data available
df <- df[which(!is.na(df$T.Cell.Inflammed.Signature.Score)),]

# Define variables
    # Categorical
    variables.cat <- c("Tumour.Grade", "TNBC", "PAM50")

    # Continuous
    variables.cont <- "Age"

# Set factor levels
    # Grade
    df$Tumour.Grade <- factor(df$Tumour.Grade, levels=c("I", "II", "III"))
    
    # TNBC
    df$TNBC <- factor(df$TNBC, levels=c("non-TNBC", "TNBC"), labels=c("No", "Yes"))

    # PAM50
    df$PAM50 <- factor(df$PAM50, levels=c("LumA", "LumB", "Basal", "Her2", "Normal"), labels=c("Luminal A", "Luminal B", "Basal-like", "HER2-enriched", "Normal-like"))

# Split data frame
df.brca <- df[which(df$BRCA_status=="BRCA1-deficient"), ]
df.proficient <- df[which(df$BRCA_status=="BRCA-proficient"), ]
   
# Build model for categorical variables
results.list <- list()

for(i in 1:length(variables.cat)) {
    
    # Subset data frames
    df.brca.small <- df.brca[, variables.cat[i], drop=FALSE]
    df.brca.small$Cohort <- "BRCA1-deficient"
    df.proficient.small <- df.proficient[, variables.cat[i], drop=FALSE]
    df.proficient.small$Cohort <- "BRCA-proficient"
    
    df.small <- rbind.data.frame(df.brca.small, df.proficient.small)
    
    # Create factor level for cohort
    df.small$Cohort <- factor(df.small$Cohort, level=c("BRCA-proficient", "BRCA1-deficient"))
    
    # Create frequency table
    tbl.freq <- table(df.small[,1], df.small[,2])
    
    # Create proprotion table
    tbl.prop <- round(prop.table(table(df.small[,1], df.small[,2]), 2)*100, digits=1)
    
    # Fisher exact test
    p.value <- rep(fisher.test(as.matrix(tbl.freq))$p.value, nrow(tbl.freq))
    
    # Tabulate
    results <- cbind.data.frame(as.data.frame.matrix(tbl.freq), as.data.frame.matrix(tbl.prop), p.value)
    
    # Create level column
    results <- data.frame("level"=row.names(results), results, stringsAsFactors=FALSE)
    
    # Create variable column
    results <- data.frame("variable"=variables.cat[i], results, stringsAsFactors=FALSE)
    
    # Save into list
    results.list[[i]] <- results

    
}

results <- do.call(rbind.data.frame, results.list)

results.cat <- as.data.frame(apply(results, 2, function(x) as.character(as.factor(x))))

# Build model for continuous variable
    # Compute median (BRCA)
    middle <- fivenum(df.brca$Age)[3]
    q1 <- fivenum(df.brca$Age)[1]
    q3 <- fivenum(df.brca$Age)[4]
    cont.brca <- paste(middle, " (", q1, "-", q3, ")", sep="")

    # Compute median (WT)
    middle <- fivenum(df.proficient$Age)[3]
    q1 <- fivenum(df.proficient$Age)[1]
    q3 <- fivenum(df.proficient$Age)[4]
    cont.proficient <- paste(middle, " (", q1, "-", q3, ")", sep="")

    # Compute p.value
    p.value <- wilcox.test(df.brca$Age, df.proficient$Age, na.rm=TRUE)$p.value

    # Tabulate (match categorical variable data frame)
    results.cont <- data.frame("variable"="Age", "level"="-", "BRCA.proficient"=cont.proficient, "BRCA1.deficient"=cont.brca, "BRCA.proficient.1"="-", "BRCA1.deficient.1"="-", p.value)

# Merge
    # Merge
    results <- rbind.data.frame(results.cont, results.cat)

    # Rename columns
    names(results) <- c("variable", "level", "TCGA_BRCA.proficient.freq", "TCGA_BRCA1.deficient.freq", "TCGA_BRCA.proficient.prop", "TCGA_BRCA1.deficient.prop", "p.value")

    # Save as new object
    results.TCGA.BRCA1 <- results

#####################################################################
######################### TCGA (BRCA2) ##############################
#####################################################################

# Read file
df <- read.table("Data/TCGA_final.txt", header=TRUE, sep="\t", na.strings="NA", stringsAsFactors=FALSE, check.names=FALSE)

# Subset relevant genotypes
df <- df[which(df$BRCA_status %in% c("BRCA2-deficient", "BRCA-proficient")), ]

# Subset samples with transcriptomic data available
df <- df[which(!is.na(df$T.Cell.Inflammed.Signature.Score)),]

# Define variables
    # Categorical
    variables.cat <- c("Tumour.Grade", "TNBC", "PAM50")

    # Continuous
    variables.cont <- "Age"

# Set factor levels
    # Grade
    df$Tumour.Grade <- factor(df$Tumour.Grade, levels=c("I", "II", "III"))
    
    # TNBC
    df$TNBC <- factor(df$TNBC, levels=c("non-TNBC", "TNBC"), labels=c("No", "Yes"))

    # PAM50
    df$PAM50 <- factor(df$PAM50, levels=c("LumA", "LumB", "Basal", "Her2", "Normal"), labels=c("Luminal A", "Luminal B", "Basal-like", "HER2-enriched", "Normal-like"))

# Split data frame
df.brca <- df[which(df$BRCA_status=="BRCA2-deficient"), ]
df.proficient <- df[which(df$BRCA_status=="BRCA-proficient"), ]
   
# Build model for categorical variables
results.list <- list()

for(i in 1:length(variables.cat)) {
    
    # Subset data frames
    df.brca.small <- df.brca[, variables.cat[i], drop=FALSE]
    df.brca.small$Cohort <- "BRCA2-deficient"
    df.proficient.small <- df.proficient[, variables.cat[i], drop=FALSE]
    df.proficient.small$Cohort <- "BRCA-proficient"
    
    df.small <- rbind.data.frame(df.brca.small, df.proficient.small)
    
    # Create factor level for cohort
    df.small$Cohort <- factor(df.small$Cohort, level=c("BRCA-proficient", "BRCA2-deficient"))
    
    # Create frequency table
    tbl.freq <- table(df.small[,1], df.small[,2])
    
    # Create proprotion table
    tbl.prop <- round(prop.table(table(df.small[,1], df.small[,2]), 2)*100, digits=1)
    
    # Fisher exact test
    p.value <- rep(fisher.test(as.matrix(tbl.freq))$p.value, nrow(tbl.freq))
    
    # Tabulate
    results <- cbind.data.frame(as.data.frame.matrix(tbl.freq), as.data.frame.matrix(tbl.prop), p.value)
    
    # Create level column
    results <- data.frame("level"=row.names(results), results, stringsAsFactors=FALSE)
    
    # Create variable column
    results <- data.frame("variable"=variables.cat[i], results, stringsAsFactors=FALSE)
    
    # Save into list
    results.list[[i]] <- results

    
}

results <- do.call(rbind.data.frame, results.list)

results.cat <- as.data.frame(apply(results, 2, function(x) as.character(as.factor(x))))

# Build model for continuous variable
    # Compute median (BRCA)
    middle <- fivenum(df.brca$Age)[3]
    q1 <- fivenum(df.brca$Age)[1]
    q3 <- fivenum(df.brca$Age)[4]
    cont.brca <- paste(middle, " (", q1, "-", q3, ")", sep="")

    # Compute median (WT)
    middle <- fivenum(df.proficient$Age)[3]
    q1 <- fivenum(df.proficient$Age)[1]
    q3 <- fivenum(df.proficient$Age)[4]
    cont.proficient <- paste(middle, " (", q1, "-", q3, ")", sep="")

    # Compute p.value
    p.value <- wilcox.test(df.brca$Age, df.proficient$Age, na.rm=TRUE)$p.value

    # Tabulate (match categorical variable data frame)
    results.cont <- data.frame("variable"="Age", "level"="-", "BRCA.proficient"=cont.proficient, "BRCA2.deficient"=cont.brca, "BRCA.proficient.1"="-", "BRCA2.deficient.1"="-", p.value)

# Merge
    # Merge
    results <- rbind.data.frame(results.cont, results.cat)

    # Rename columns
    names(results) <- c("variable", "level", "TCGA_BRCA.proficient.freq", "TCGA_BRCA2.deficient.freq", "TCGA_BRCA.proficient.prop", "TCGA_BRCA2.deficient.prop", "p.value")

    # Save as new object
    results.TCGA.BRCA2 <- results

#####################################################################
######################### TCGA (BRCA2) ##############################
#####################################################################

# Merge
header <- results.WSI.BRCA1[,c(1:2)]
results.WSI.BRCA1 <- results.WSI.BRCA1[,-c(1:2)]
results.WSI.BRCA2 <- results.WSI.BRCA2[,-c(1:2)]
results.TCGA.BRCA1 <- results.TCGA.BRCA1[,-c(1:2)]
results.TCGA.BRCA2 <- results.TCGA.BRCA2[,-c(1:2)]
results <- cbind.data.frame(header, results.WSI.BRCA1, results.WSI.BRCA2, results.TCGA.BRCA1, results.TCGA.BRCA2)

# Write table
write.table(results, "Tables/Supp Table_2.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

