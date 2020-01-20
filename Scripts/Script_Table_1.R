# Set working directory
setwd("/Users/WIMM/Documents/BRCA_2018/Github/")

###################################################################
############################### WSI ###############################
###################################################################

# Read file
df <- read.table("Data/WSI_final.txt", header=TRUE, sep="\t", na.strings="NA", check.names=FALSE)

# Subset relevant columns
df <- df[,c("Sample", "BRCA_status", "Age", "Tumour.Grade", "T.Stage", "ER", "PR", "HER2", "TNBC", "PAM50", "T.Cell.Inflammed.Signature.Score")]

# Set factor levels
    # BRCA
    df$BRCA_status <- factor(df$BRCA_status, c("BRCA-proficient", "BRCA1-deficient", "BRCA2-deficient"))
    
    # Age
    df$Age <- ifelse(df$Age <= median(df$Age, na.rm=TRUE), "=<median", ">median")
    df$Age <- factor(df$Age, levels=c(">median", "=<median"))

    # Grade
    df$Tumour.Grade <- factor(df$Tumour.Grade, levels=c("I", "II", "III"))
    
    # Stage
    df$T.Stage <- factor(df$T.Stage, levels=c("I", "II", "III", "IV"))
    
    # ER
    df$ER <- factor(df$ER, levels=c("positive", "negative"))
    
    # PR
    df$PR <- factor(df$PR, levels=c("positive", "negative"))
    
    # HER2
    df$HER2 <- factor(df$HER2, levels=c("positive", "negative"))
    
    # TNBA
    df$TNBC <- factor(df$TNBC, levels=c("non-TNBC", "TNBC"), labels=c("No", "Yes"))
    
    # PAM50
    df$PAM50 <- factor(df$PAM50, levels=c("LumA", "LumB", "Basal", "Her2", "Normal"))

# Build model
variables <- c("BRCA_status", "Age", "Tumour.Grade", "T.Stage", "ER", "PR", "HER2", "TNBC", "PAM50")

model.list <- list()

for(i in 1:length(variables)) {

    model.list[[i]] <- lm(as.formula(paste("T.Cell.Inflammed.Signature.Score ~ ", variables[i], sep="")), df)

}

# Retrieve coefficients
coef <- round(as.data.frame(unlist(lapply(model.list, function(x) x$coefficient[-1]))), digits=2)

# Retreive SE
se <- round(as.data.frame(unlist(lapply(model.list, function(x) summary(x)$coefficient[-1,2]))), digits=2)

# Retrieve p-values
p.value <- as.data.frame(unlist(lapply(model.list, function(x) summary(x)$coefficient[-1,4])))

# Tabulate
results <- cbind.data.frame(coef, se, p.value)
names(results) <- paste("WSI_", c("coef", "se", "p.value"), sep="")

# Save as new object
results.WSI <- results

###################################################################
############################### TCGA ###############################
###################################################################

# Read file
df <- read.table("Data/TCGA_final.txt", header=TRUE, sep="\t", na.strings="NA", check.names=FALSE)

# Subset relevant columns
df <- df[,c("Sample", "BRCA_status", "Age", "Tumour.Grade", "T.Stage", "ER", "PR", "HER2", "TNBC", "PAM50", "T.Cell.Inflammed.Signature.Score")]

# Set factor levels
    # BRCA
    df$BRCA_status <- factor(df$BRCA_status, c("BRCA-proficient", "BRCA1-deficient", "BRCA2-deficient"))
    
    # Age
    df$Age <- ifelse(df$Age <= median(df$Age, na.rm=TRUE), "=<median", ">median")
    df$Age <- factor(df$Age, levels=c(">median", "=<median"))

    # Grade
    df$Tumour.Grade <- factor(df$Tumour.Grade, levels=c("I", "II", "III"))
    
    # Stage
    df$T.Stage <- factor(df$T.Stage, levels=c("I", "II", "III", "IV"))
    
    # ER
    df$ER <- factor(df$ER, levels=c("positive", "negative"))
    
    # PR
    df$PR <- factor(df$PR, levels=c("positive", "negative"))
    
    # HER2
    df$HER2 <- factor(df$HER2, levels=c("positive", "negative"))
    
    # TNBC
    df$TNBC <- factor(df$TNBC, levels=c("non-TNBC", "TNBC"), labels=c("No", "Yes"))
    
    # PAM50
    df$PAM50 <- factor(df$PAM50, levels=c("LumA", "LumB", "Basal", "Her2", "Normal"))

# Build model
variables <- c("BRCA_status", "Age", "Tumour.Grade", "T.Stage", "ER", "PR", "HER2", "TNBC", "PAM50")

model.list <- list()

for(i in 1:length(variables)) {

    model.list[[i]] <- lm(as.formula(paste("T.Cell.Inflammed.Signature.Score ~ ", variables[i], sep="")), df)

}

# Retrieve coefficients
coef <- round(as.data.frame(unlist(lapply(model.list, function(x) x$coefficient[-1]))), digits=2)

# Retreive SE
se <- round(as.data.frame(unlist(lapply(model.list, function(x) summary(x)$coefficient[-1,2]))), digits=2)

# Retrieve p-values
p.value <- as.data.frame(unlist(lapply(model.list, function(x) summary(x)$coefficient[-1,4])))

# Tabulate
results <- cbind.data.frame(coef, se, p.value)
names(results) <- paste("TCGA_", c("coef", "se", "p.value"), sep="")

# Save as new object
results.TCGA <- results

###################################################################
############################### MERGE #############################
###################################################################

# Merge
results <- cbind.data.frame(results.WSI, results.TCGA)

# Add level column
results <- data.frame("level"=row.names(results), results, stringsAsFactors=FALSE)

# Write table
write.table(results, "Tables/Table_1.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
