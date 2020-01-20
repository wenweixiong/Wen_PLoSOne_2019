# Load packages
library(ggplot2)

# Set working directory
setwd("/Users/WIMM/Documents/BRCA_2018/Github/")

# Read file
df <- read.table("Data/TCGA_final.txt", header=TRUE, sep="\t", na.strings="NA", check.names=FALSE)

# Subset relevant columns
df <- df[ , c("BRCA_status", "Mutation.Load.Whole.Exome")]

# Remove missing rows
df <- na.omit(df)

# Compute p-values
    # BRCA1-deficient vs. BRCA-proficient
    wilcox.test(Mutation.Load.Whole.Exome ~ BRCA_status, data=subset(df, BRCA_status %in% c("BRCA-proficient", "BRCA1-deficient")))$p.value
    
    # BRCA2-deficient vs. BRCA-proficient
    wilcox.test(Mutation.Load.Whole.Exome ~ BRCA_status, data=subset(df, BRCA_status %in% c("BRCA-proficient", "BRCA2-deficient")))$p.value
    
# Tabulate sample size per genotype
n.proficient <- length(df$BRCA_status[which(df$BRCA_status=="BRCA-proficient")])
n.brca1 <- length(df$BRCA_status[which(df$BRCA_status=="BRCA1-deficient")])
n.brca2 <- length(df$BRCA_status[which(df$BRCA_status=="BRCA2-deficient")])

# Boxplot
    # Definition
    data <- df
    x <- data$BRCA_status
    y <- log10(data$Mutation.Load.Whole.Exome)
    maintitle <- "TCGA"
    xtitle <- ""
    ytitle <- "log10(Number of mutations)"
    xlabels <- c(paste("BRCA-proficient\nbreast cancers\n(n=", n.proficient, ")", sep=""),
     paste("BRCA1-deficient\nbreast cancers\n(n=", n.brca1, ")", sep=""),
     paste("BRCA2-deficient\nbreast cancers\n(n=", n.brca2, ")", sep=""))
    ymin <- floor(min(y))
    ymax <- ceiling(max(y))
    
    # Plot
    plot <- ggplot(data, aes(x=x, y=y)) +
        geom_boxplot(outlier.colour="white") +
        geom_jitter(height=0, width=0.05, color="blue") +
        ylim(ymin, ymax) +
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
    ggsave("Figures/Figure_1D.pdf", plot, width=8, height=8)
