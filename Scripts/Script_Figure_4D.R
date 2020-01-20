# Load packages
library(ggplot2)
library(Rtsne)

# Set working directory
setwd("/Users/WIMM/Documents/BRCA_2018/Github/")

# Read file
df <- read.table("Data/TCGA_final.txt", header=TRUE, sep="\t", na.strings="NA", check.names=FALSE)

# Subset relevant columns
genes <- c("CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK", "CD2", "HLA-DRA", "CXCL13", "IL2RG", "NKG7", "HLA-E", "CXCR6", "LAG3", "TAGAP", "CXCL10", "STAT1", "GZMB")
df <- df[, c("BRCA_status", "PTEN_Status", "PTEN_RPPA", "T.Cell.Inflammed.Signature.Score", genes)]

# Remove missing values
df <- na.omit(df)

# Define PTEN MT, WT
df$PTEN_Status <- ifelse(df$PTEN_Status=="No CNA or point mutations" | df$PTEN_Status=="Amplification", "Wildtype", "Mutant")

# Set factor levels
df$PTEN_Status <- as.factor(df$PTEN_Status)
df$PTEN_Status <- factor(df$PTEN_Status, levels=c("Wildtype", "Mutant"))

# Subset BRCA-deficient samples
df <- df[which(df$BRCA_status %in% c("BRCA1-deficient", "BRCA2-deficient")),]

# Statify inflammed signature
df$T.Cell.Inflammed.Signature.Score <- ifelse(df$T.Cell.Inflammed.Signature.Score > median(df$T.Cell.Inflammed.Signature.Score), "> median", "=< median")

# Reduce dimensions
    # Set seed
    set.seed(123)

    # Reduce dimension
    tsne <- Rtsne(df[,genes], dims=2, perplexity=10)

    # Retrieve dimension
    tsne <- data.frame(tsne$Y)

# Scatterplot (label by T.Cell.Inflamed.Signature.Score)
    # Definitions
    data <- df
    x <- tsne[,1]
    y <- tsne[,2]
    z <- data$PTEN_Status
    xtitle <- "tSNE1"
    ytitle <- "tSNE2"
    legend.title <- "PTEN mutation\nstatus"

    # Plot
    plot <- ggplot() +
        geom_point(data, mapping=aes(x=x, y=y, color=z), size=4) +
        scale_color_manual(name=legend.title, breaks=unique(z), values=c("lightblue", "pink")) +
        labs(x=xtitle, y=ytitle) +
        theme(legend.title=element_text(size=10),
              legend.text=element_text(size=10),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
              axis.line=element_line(colour = "black"),
              axis.text=element_text(size=13),
              axis.title=element_text(size=18),
              axis.text.x=element_text(colour="black"))

    # Save plot
    ggsave("Figures/Figure_4D.pdf", plot=plot, width=8, height=7)



