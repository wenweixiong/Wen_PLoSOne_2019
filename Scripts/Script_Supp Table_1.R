# Set working directory
setwd("/Users/WIMM/Documents/BRCA_2018/Github/")

# Read file
    # WSI
    df.wsi <- read.table("Data/WSI_final.txt", header=TRUE, sep="\t", na.strings="NA", stringsAsFactors=FALSE, check.names=FALSE)

    # TCGA
    df.tcga <- read.table("Data/TCGA_final.txt", header=TRUE, sep="\t", na.strings="NA", stringsAsFactors=FALSE, check.names=FALSE)

# Create new variable for germline BRCA status
    # WSI
    df.wsi$BRCA_status_germline <- df.wsi$BRCA_status_allele
    df.wsi$BRCA_status_germline[which(df.wsi$BRCA_status_germline=="BRCA1 germline pathogenic (biallelic)")] <- "BRCA1-deficient"
    df.wsi$BRCA_status_germline[which(df.wsi$BRCA_status_germline=="BRCA2 germline pathogenic (biallelic)")] <- "BRCA2-deficient"
    df.wsi$BRCA_status_germline[which(df.wsi$BRCA_status_germline %in% c("BRCA1 somatic pathogenic (biallelic)", "BRCA1 promoter hypermethylation", "BRCA2 somatic pathogenic (biallelic)", "BRCA2 homozygous deletion"))] <- NA
    df.wsi$BRCA_status_germline[which(df.wsi$BRCA_status_germline %in% c("BRCA1 germline pathogenic (monoallelic)", "BRCA1 somatic pathogenic (monoallelic)", "BRCA2 germline pathogenic (monoallelic)", "BRCA2 somatic pathogenic (monoallelic)"))] <- NA
    df.wsi$BRCA_status_germline[-which(df.wsi$BRCA_status_germline=="BRCA1-deficient" | df.wsi$BRCA_status_germline=="BRCA2-deficient" | is.na(df.wsi$BRCA_status_germline))] <-  "BRCA-proficient"

    # TCGA
    df.tcga$BRCA_status_germline <- df.tcga$BRCA_status_allele
    df.tcga$BRCA_status_germline[which(df.tcga$BRCA_status_germline=="BRCA1 germline pathogenic (biallelic)")] <- "BRCA1-deficient"
    df.tcga$BRCA_status_germline[which(df.tcga$BRCA_status_germline=="BRCA2 germline pathogenic (biallelic)")] <- "BRCA2-deficient"
    df.tcga$BRCA_status_germline[which(df.tcga$BRCA_status_germline %in% c("BRCA1 somatic pathogenic (biallelic)", "BRCA1 homozygous deletion", "BRCA1 promoter hypermethylation", "BRCA2 somatic pathogenic (biallelic)", "BRCA2 homozygous deletion"))] <- NA
    df.tcga$BRCA_status_germline[which(df.tcga$BRCA_status_germline %in% c("BRCA1 germline pathogenic (monoallelic)", "BRCA1 somatic pathogenic (monoallelic)", "BRCA2 germline pathogenic (monoallelic)", "BRCA2 somatic pathogenic (monoallelic)"))] <- NA
    df.tcga$BRCA_status_germline[-which(df.tcga$BRCA_status_germline=="BRCA1-deficient" | df.tcga$BRCA_status_germline=="BRCA2-deficient" | is.na(df.tcga$BRCA_status_germline))] <-  "BRCA-proficient"

# Create new variable for somatic BRCA status
    # WSI
    df.wsi$BRCA_status_somatic <- df.wsi$BRCA_status_allele
    df.wsi$BRCA_status_somatic[which(df.wsi$BRCA_status_somatic %in% c("BRCA1 somatic pathogenic (biallelic)", "BRCA1 promoter hypermethylation"))] <- "BRCA1-deficient"
    df.wsi$BRCA_status_somatic[which(df.wsi$BRCA_status_somatic %in% c("BRCA2 somatic pathogenic (biallelic)", "BRCA2 homozygous deletion"))] <- "BRCA2-deficient"
    df.wsi$BRCA_status_somatic[which(df.wsi$BRCA_status_somatic %in% c("BRCA1 germline pathogenic (biallelic)", "BRCA2 germline pathogenic (biallelic)"))] <- NA
    df.wsi$BRCA_status_somatic[which(df.wsi$BRCA_status_somatic %in% c("BRCA1 germline pathogenic (monoallelic)", "BRCA1 somatic pathogenic (monoallelic)", "BRCA2 germline pathogenic (monoallelic)", "BRCA2 somatic pathogenic (monoallelic)"))] <- NA
    df.wsi$BRCA_status_somatic[-which(df.wsi$BRCA_status_somatic=="BRCA1-deficient" | df.wsi$BRCA_status_somatic=="BRCA2-deficient" | is.na(df.wsi$BRCA_status_somatic))] <-  "BRCA-proficient"

    # TCGA
    df.tcga$BRCA_status_somatic <- df.tcga$BRCA_status_allele
    df.tcga$BRCA_status_somatic[which(df.tcga$BRCA_status_somatic %in% c("BRCA1 somatic pathogenic (biallelic)", "BRCA1 promoter hypermethylation", "BRCA1 homozygous deletion"))] <- "BRCA1-deficient"
    df.tcga$BRCA_status_somatic[which(df.tcga$BRCA_status_somatic %in% c("BRCA2 somatic pathogenic (biallelic)", "BRCA2 homozygous deletion"))] <- "BRCA2-deficient"
    df.tcga$BRCA_status_somatic[which(df.tcga$BRCA_status_somatic %in% c("BRCA1 germline pathogenic (biallelic)", "BRCA2 germline pathogenic (biallelic)"))] <- NA
    df.tcga$BRCA_status_somatic[which(df.tcga$BRCA_status_somatic %in% c("BRCA1 germline pathogenic (monoallelic)", "BRCA1 somatic pathogenic (monoallelic)", "BRCA2 germline pathogenic (monoallelic)", "BRCA2 somatic pathogenic (monoallelic)"))] <- NA
    df.tcga$BRCA_status_somatic[-which(df.tcga$BRCA_status_somatic=="BRCA1-deficient" | df.tcga$BRCA_status_somatic=="BRCA2-deficient" | is.na(df.tcga$BRCA_status_somatic))] <-  "BRCA-proficient"

# Define variables
    # Categorical
    variables.cat <- c("Menopause.Status", "Tumour.Grade", "ER", "PR", "HER2", "TNBC", "BRCA_status", "BRCA_status_germline", "BRCA_status_somatic")

    # Continuous
    variables.cont <- "Age"

# Set factor levels
    # Menopause Status
    df.wsi$Menopause.Status <- factor(df.wsi$Menopause.Status, c("pre-menopausal", "post-menopausal"))
    df.tcga$Menopause.Status <- factor(df.tcga$Menopause.Status, c("pre-menopausal", "post-menopausal"))
    
    # Grade
    df.wsi$Tumour.Grade <- factor(df.wsi$Tumour.Grade, levels=c("I", "II", "III"))
    df.tcga$Tumour.Grade <- factor(df.tcga$Tumour.Grade, levels=c("I", "II", "III"))
    
    # ER
    df.wsi$ER <- factor(df.wsi$ER, levels=c("positive", "negative"))
    df.tcga$ER <- factor(df.tcga$ER, levels=c("positive", "negative"))
    
    # PR
    df.wsi$PR <- factor(df.wsi$PR, levels=c("positive", "negative"))
    df.tcga$PR <- factor(df.tcga$PR, levels=c("positive", "negative"))
    
    # HER2
    df.wsi$HER2 <- factor(df.wsi$HER2, levels=c("positive", "negative"))
    df.tcga$HER2 <- factor(df.tcga$HER2, levels=c("positive", "negative"))
    
    # TNBC
    df.wsi$TNBC <- factor(df.wsi$TNBC, levels=c("non-TNBC", "TNBC"), labels=c("No", "Yes"))
    df.tcga$TNBC <- factor(df.tcga$TNBC, levels=c("non-TNBC", "TNBC"), labels=c("No", "Yes"))
    
    # BRCA status
    df.wsi$BRCA_status <- factor(df.wsi$BRCA_status, c("BRCA-proficient", "BRCA1-deficient", "BRCA2-deficient"))
    df.tcga$BRCA_status <- factor(df.tcga$BRCA_status, c("BRCA-proficient", "BRCA1-deficient", "BRCA2-deficient"))

    # BRCA status (germline)
    df.wsi$BRCA_status_germline <- factor(df.wsi$BRCA_status_germline, c("BRCA-proficient", "BRCA1-deficient", "BRCA2-deficient"))
    df.tcga$BRCA_status_germline <- factor(df.tcga$BRCA_status_germline, c("BRCA-proficient", "BRCA1-deficient", "BRCA2-deficient"))

    # BRCA status (somatic)
    df.wsi$BRCA_status_somatic <- factor(df.wsi$BRCA_status_somatic, c("BRCA-proficient", "BRCA1-deficient", "BRCA2-deficient"))
    df.tcga$BRCA_status_somatic <- factor(df.tcga$BRCA_status_somatic, c("BRCA-proficient", "BRCA1-deficient", "BRCA2-deficient"))

# Build model for categorical variable
results.list <- list()

for(i in 1:length(variables.cat)) {
    
    # Subset data frames
    df.wsi.small <- df.wsi[, variables.cat[i], drop=FALSE]
    df.wsi.small$Cohort <- "WSI"
    df.tcga.small <- df.tcga[, variables.cat[i], drop=FALSE]
    df.tcga.small$Cohort <- "TCGA"
    df.small <- rbind.data.frame(df.wsi.small, df.tcga.small)
    
    # Create factor level for cohort
    df.small$Cohort <- factor(df.small$Cohort, level=c("WSI", "TCGA"))
    
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
    # Compute median (WSI)
    middle <- fivenum(df.wsi$Age)[3]
    q1 <- fivenum(df.wsi$Age)[1]
    q3 <- fivenum(df.wsi$Age)[4]
    cont.wsi <- paste(middle, " (", q1, "-", q3, ")", sep="")

    # Compute median (TCGA)
    middle <- fivenum(df.tcga$Age)[3]
    q1 <- fivenum(df.tcga$Age)[1]
    q3 <- fivenum(df.tcga$Age)[4]
    cont.tcga <- paste(middle, " (", q1, "-", q3, ")", sep="")

    # Compute p.value
    p.value <- wilcox.test(df.wsi$Age, df.tcga$Age, na.rm=TRUE)$p.value

    # Tabulate (match categorical variable data frame)
    results.cont <- data.frame("variable"="Age", "level"="-", "WSI"=cont.wsi, "TCGA"=cont.tcga, "WSI.1"="-", "TCGA.1"="-", p.value)

# Merge
    # Merge
    results <- rbind.data.frame(results.cont, results.cat)

    # Rename columns
    names(results) <- c("variable", "level", "WSI.freq", "TCGA.freq", "WSI.prop", "TCGA.prop", "p.value")

    # Write table
    write.table(results, "Tables/Supp Table_1.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
