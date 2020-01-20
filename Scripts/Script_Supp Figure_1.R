# Load packages
library(ggplot2)

# Set working directory
setwd("/Users/WIMM/Documents/BRCA_2018/Github/")

# Read file
df <- read.table("Data/WSI_final.txt", header=TRUE, sep="\t", na.strings="NA", check.names=FALSE)

# Subset relevant columns
df <- df[ , c("BRCA_status", "Signature.3.Whole.Exome")]

# Remove missing rows
df <- na.omit(df)

# Compute p-values
    # BRCA1-deficient vs. BRCA-proficient
    wilcox.test(Signature.3.Whole.Exome ~ BRCA_status, data=subset(df, BRCA_status %in% c("BRCA-proficient", "BRCA1-deficient")))$p.value
    
    # BRCA2-deficient vs. BRCA-proficient
    wilcox.test(Signature.3.Whole.Exome ~ BRCA_status, data=subset(df, BRCA_status %in% c("BRCA-proficient", "BRCA2-deficient")))$p.value
    
# Tabulate sample size per genotype
n.proficient <- length(df$BRCA_status[which(df$BRCA_status=="BRCA-proficient")])
n.brca1 <- length(df$BRCA_status[which(df$BRCA_status=="BRCA1-deficient")])
n.brca2 <- length(df$BRCA_status[which(df$BRCA_status=="BRCA2-deficient")])

# Boxplot
    # Definition
    data <- df
    x <- data$BRCA_status
    y <- data$Signature.3.Whole.Exome
    maintitle <- "WSI"
    xtitle <- ""
    ytitle <- "Proportion of mutational signature 3"
    xlabels <- c(paste("BRCA-proficient\nbreast cancers\n(n=", n.proficient, ")", sep=""),
     paste("BRCA1-deficient\nbreast cancers\n(n=", n.brca1, ")", sep=""),
     paste("BRCA2-deficient\nbreast cancers\n(n=", n.brca2, ")", sep=""))
    
    
    # Plot
    plot <- ggplot(data, aes(x=x, y=y)) +
        geom_boxplot(outlier.colour="white") +
        geom_jitter(height=0, width=0.05, color="blue") +
        ylim(0, 1.0) +
        labs(title=maintitle, x=xtitle, y=ytitle) +
        scale_x_discrete(labels=xlabels) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
              axis.line = element_line(colour="black"),
              axis.text=element_text(size=16),
              axis.title=element_text(size=18),
              axis.text.x=element_text(colour="black"),
              plot.title=element_text(size=18))

    # Sava plot
    ggsave("Figures/Supp Figure_1.pdf", plot, width=8, height=8)
