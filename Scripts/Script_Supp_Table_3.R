# Set working directory
setwd("/Users/WIMM/Documents/BRCA_2018/Github/")

###################################################################
############################### WSI ###############################
###################################################################

# Read file
df <- read.table("Data/WSI_final.txt", header=TRUE, sep="\t", na.strings="NA", check.names=FALSE)

# Subset relevant columns
df <- df[,c("Sample", "BRCA_status", "Age", "Tumour.Grade", "TNBC", "PAM50", "T.Cell.Inflammed.Signature.Score")]

# Set factor levels
    # BRCA
    df$BRCA_status <- factor(df$BRCA_status, c("BRCA-proficient", "BRCA1-deficient", "BRCA2-deficient"))
    
    # Age
    df$Age <- ifelse(df$Age <= median(df$Age, na.rm=TRUE), "=<median", ">median")
    df$Age <- factor(df$Age, levels=c(">median", "=<median"))

    # Grade
    df$Tumour.Grade <- factor(df$Tumour.Grade, levels=c("I", "II", "III"))
    
    # TNBC
    df$TNBC <- factor(df$TNBC, levels=c("non-TNBC", "TNBC"), labels=c("No", "Yes"))
    
    # PAM50
    df$PAM50 <- factor(df$PAM50, levels=c("LumA", "LumB", "Basal", "Her2", "Normal"))

# Adjusted for age, grade, TNBC
    # Build model
    variables <- paste(c("BRCA_status", "Age", "Tumour.Grade", "TNBC"), collapse=" + ")
    formula <- as.formula(paste("T.Cell.Inflammed.Signature.Score ~ ", variables, sep=""))
    model <- lm(formula, df)

    # Retrieve coefficients
    coef <- round(model$coefficient[-1], digits=2)

    # Retreive SE
    se <- round(summary(model)$coefficient[-1,2], digits=2)

    # Retrieve p-values
    p.value <- summary(model)$coefficient[-1,4]

    # Tabulate
    results <- cbind.data.frame(coef, se, p.value)
    names(results) <- paste("WSI_", c("coef", "se", "p.value"), sep="")

    # Save as new object
    results.WSI.model.1 <- results[c(1:2),]

# Adjusted for age, grade, PAM50
    # Build model
    variables <- paste(c("BRCA_status", "Age", "Tumour.Grade", "PAM50"), collapse=" + ")
    formula <- as.formula(paste("T.Cell.Inflammed.Signature.Score ~ ", variables, sep=""))
    model <- lm(formula, df)

    # Retrieve coefficients
    coef <- round(model$coefficient[-1], digits=2)

    # Retreive SE
    se <- round(summary(model)$coefficient[-1,2], digits=2)

    # Retrieve p-values
    p.value <- summary(model)$coefficient[-1,4]

    # Tabulate
    results <- cbind.data.frame(coef, se, p.value)
    names(results) <- paste("WSI_", c("coef", "se", "p.value"), sep="")

    # Save as new object
    results.WSI.model.2 <- results[c(1:2),]

# Merge
results.WSI <- rbind.data.frame(results.WSI.model.1, results.WSI.model.2)

###################################################################
############################### TCGA ###############################
###################################################################

# Read file
df <- read.table("Data/TCGA_final.txt", header=TRUE, sep="\t", na.strings="NA", check.names=FALSE)

# Subset relevant columns
df <- df[,c("Sample", "BRCA_status", "Age", "Tumour.Grade", "TNBC", "PAM50", "T.Cell.Inflammed.Signature.Score")]

# Set factor levels
    # BRCA
    df$BRCA_status <- factor(df$BRCA_status, c("BRCA-proficient", "BRCA1-deficient", "BRCA2-deficient"))
    
    # Age
    df$Age <- ifelse(df$Age <= median(df$Age, na.rm=TRUE), "=<median", ">median")
    df$Age <- factor(df$Age, levels=c(">median", "=<median"))

    # Grade
    df$Tumour.Grade <- factor(df$Tumour.Grade, levels=c("I", "II", "III"))
    
    # TNBC
    df$TNBC <- factor(df$TNBC, levels=c("non-TNBC", "TNBC"), labels=c("No", "Yes"))
    
    # PAM50
    df$PAM50 <- factor(df$PAM50, levels=c("LumA", "LumB", "Basal", "Her2", "Normal"))

# Adjusted for age, grade, TNBC
    # Build model
    variables <- paste(c("BRCA_status", "Age", "Tumour.Grade", "TNBC"), collapse=" + ")
    formula <- as.formula(paste("T.Cell.Inflammed.Signature.Score ~ ", variables, sep=""))
    model <- lm(formula, df)

    # Retrieve coefficients
    coef <- round(model$coefficient[-1], digits=2)

    # Retreive SE
    se <- round(summary(model)$coefficient[-1,2], digits=2)

    # Retrieve p-values
    p.value <- summary(model)$coefficient[-1,4]

    # Tabulate
    results <- cbind.data.frame(coef, se, p.value)
    names(results) <- paste("TCGA_", c("coef", "se", "p.value"), sep="")

    # Save as new object
    results.TCGA.model.1 <- results[c(1:2),]

# Adjusted for age, grade, PAM50
    # Build model
    variables <- paste(c("BRCA_status", "Age", "Tumour.Grade", "PAM50"), collapse=" + ")
    formula <- as.formula(paste("T.Cell.Inflammed.Signature.Score ~ ", variables, sep=""))
    model <- lm(formula, df)

    # Retrieve coefficients
    coef <- round(model$coefficient[-1], digits=2)

    # Retreive SE
    se <- round(summary(model)$coefficient[-1,2], digits=2)

    # Retrieve p-values
    p.value <- summary(model)$coefficient[-1,4]

    # Tabulate
    results <- cbind.data.frame(coef, se, p.value)
    names(results) <- paste("TCGA_", c("coef", "se", "p.value"), sep="")

    # Save as new object
    results.TCGA.model.2 <- results[c(1:2),]

# Merge
results.TCGA <- rbind.data.frame(results.TCGA.model.1, results.TCGA.model.2)

###################################################################
############################### MERGE #############################
###################################################################

# Merge
results <- cbind.data.frame(results.WSI, results.TCGA)

# Add level column
results <- data.frame("level"=row.names(results), results, stringsAsFactors=FALSE)
results$level <- gsub("1$", "", results$level)
results$level <- gsub("BRCA_status", "", results$level)

# Write table
write.table(results, "Tables/Supp Table_3.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
