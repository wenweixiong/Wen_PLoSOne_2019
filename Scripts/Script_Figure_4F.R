# Load packages
library(ggplot2)

# Set working directory
setwd("/Users/WIMM/Documents/BRCA_2018/Github/")

# Read file
df <- read.table("Data/TCGA_final.txt", header=TRUE, sep="\t", na.strings="NA", check.names=FALSE)

# Subset relevant columns
df <- df[, c("BRCA_status", "PTEN_Status", "PTEN_RPPA", "T.Cell.Inflammed.Signature.Score")]

# Remove missing values
df <- na.omit(df)

# Define PTEN MT, WT
df$PTEN_Status <- ifelse(df$PTEN_Status=="No CNA or point mutations" | df$PTEN_Status=="Amplification", "Wildtype", "Mutant")

# Set factor levels
df$PTEN_Status <- as.factor(df$PTEN_Status)
df$PTEN_Status <- factor(df$PTEN_Status, levels=c("Wildtype", "Mutant"))

# Subset BRCA1-deficient samples
df <- df[which(df$BRCA_status %in% c("BRCA1-deficient", "BRCA-proficient")),]

# Merge BRCA, PTEN status
df$BRCA_PTEN_status <- paste(df$BRCA_status, "_PTEN-", df$PTEN_Status, sep="")
df$BRCA_PTEN_status <- gsub("BRCA\\-proficient.*", "BRCA\\-proficient", df$BRCA_PTEN_status)

# Set factor levels
table(df$BRCA_PTEN_status)
df$BRCA_PTEN_status <- as.factor(df$BRCA_PTEN_status)
df$BRCA_PTEN_status <- factor(df$BRCA_PTEN_status, levels=c("BRCA-proficient", "BRCA1-deficient_PTEN-Wildtype", "BRCA1-deficient_PTEN-Mutant"))

# Compute p-values
    # BRCA1-deficient_PTEN-Wildtype vs. BRCA-proficient
    wilcox.test(T.Cell.Inflammed.Signature.Score ~ BRCA_PTEN_status, data=subset(df, BRCA_PTEN_status %in% c("BRCA1-deficient_PTEN-Wildtype", "BRCA-proficient")))$p.value

    # BRCA1-deficient_PTEN-Mutant vs. BRCA-proficient
    wilcox.test(T.Cell.Inflammed.Signature.Score ~ BRCA_PTEN_status, data=subset(df, BRCA_PTEN_status %in% c("BRCA1-deficient_PTEN-Mutant", "BRCA-proficient")))$p.value

# Tabulate counts
n.proficient <- length(df$BRCA_PTEN_status[which(df$BRCA_PTEN_status=="BRCA-proficient")])
n.brca1.pten.wt <- length(df$BRCA_PTEN_status[which(df$BRCA_PTEN_status=="BRCA1-deficient_PTEN-Wildtype")])
n.brca1.pten.mt <- length(df$BRCA_PTEN_status[which(df$BRCA_PTEN_status=="BRCA1-deficient_PTEN-Mutant")])

# Density plot
    # Definition
    data <- df
    x <- data$BRCA_PTEN_status
    y <- data$T.Cell.Inflammed.Signature.Score
    maintitle <- ""
    xtitle <- ""
    ytitle="T cell-inflamed signature score"
    xticks=c(paste("BRCA-proficent\nbreast cancers\n(n=", n.proficient, ")", sep=""),
             paste("BRCA1-deficient\nbreast cancers\n(PTEN wildtype)\n(n=", n.brca1.pten.wt, ")", sep=""),
             paste("BRCA1-deficient\nbreast cancers\n(PTEN mutant)\n(n=", n.brca1.pten.mt, ")", sep=""))

    # Plot
    plot <- ggplot(data, aes(x=x, y=y, fill=x)) +
            geom_boxplot(outlier.colour = "white") +
            geom_jitter(height = 0, width = 0.05, color="black") +
            scale_fill_manual(values=c("gray", "blue", "red")) +
            scale_x_discrete(labels=xticks) +
            labs(title=NULL, x=NULL, y=ytitle) +
            theme(panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.background=element_blank(),
                  axis.line = element_line(colour="black"),
                  axis.text=element_text(size=18),
                  axis.title=element_text(size=18),
                  legend.title=element_text(size=18),
                  legend.text=element_text(size=18),
                  axis.text.x=element_text(colour="black"),
                  legend.position="none")
    # Save plot
    ggsave("Figures/Figure_4G.pdf", plot, width=9, height=7)
