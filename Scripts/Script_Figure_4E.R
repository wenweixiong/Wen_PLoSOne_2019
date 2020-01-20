# Load packages
library(ggplot2)
library(Rtsne)

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

# Density plot
    # Definition
    data <- df
    x <- data$T.Cell.Inflammed.Signature.Score
    z <- data$BRCA_PTEN_status
    maintitle <- ""
    xtitle <- "T cell-inflamed signature score"
    ytitle <- "Density"
    legend.title <- "BRCA and PTEN status"
    legend.label <- c("BRCA-proficient", "BRCA1-deficient (PTEN wildtype)", "BRCA1-deficient (PTEN mutant)")

    # Plot
    plot <- ggplot(data, aes(x=x, fill=z)) +
            geom_density(alpha=0.5) +
            scale_fill_manual(values=c("grey", "blue", "red"), name=legend.title, labels=legend.label) +
            labs(title=maintitle, x=xtitle, y=ytitle) +
            theme(panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.background=element_blank(),
                  axis.line=element_line(colour="black"),
                  axis.text=element_text(size=13),
                  axis.title=element_text(size=18),
                  axis.text.x=element_text(colour="black"),
                  plot.title=element_text(hjust=0.5, size=18))

    # Save plot
    ggsave("Figures/Figure_4F.pdf", plot=plot, width=9, height=7)
