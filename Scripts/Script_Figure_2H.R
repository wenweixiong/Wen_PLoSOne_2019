# Load packages
library(ggplot2)
library(ggrepel)

# Set working directory
setwd("/Users/WIMM/Documents/BRCA_2018/Github/")

# Read file
df <- read.table("Data/TCGA_final.txt", header=TRUE, sep="\t", na.strings="NA", check.names=FALSE)

# CIN25 gene list
genes <- c("CTLA4", "PDCD1", "LAG3", "BTLA", "CD160", "IDO1", "IL10", "IL10RB", "TGFB1", "TGFBR1", "VTCN1", "CD244", "LGALS9", "HAVCR2", "ADORA2A", "TIGIT", "CSF1R", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KDR", "CD96", "PVRL2")

genes <- intersect(names(df), genes)

# Subset relevant columns
df <- df[,c("BRCA_status", genes)]

# Remove missing values
df <- na.omit(df)

# Compute FC and p-values
fc.brca <- NULL
p.brca <- NULL

for(i in 1:length(genes)) {

    # Subset gene
    df.small <- df[, c("BRCA_status", genes[i])]
    
    # Retrieve expression values
    proficient <- df.small[which(df.small$BRCA_status=="BRCA-proficient"), 2]
    brca <- df.small[which(df.small$BRCA_status=="BRCA2-deficient"), 2]
    
    # Compute FC
    fc.brca[i] <- log2(median((2^brca - 1), na.rm=TRUE) / median((2^proficient - 1), na.rm=TRUE))
    
    # Compute p.values
    p.brca[i] <- wilcox.test(brca, proficient, na.rm=TRUE)$p.value
    
}

# Tabulate
results <- data.frame(genes, fc.brca, p.brca, stringsAsFactors=FALSE)

# Indicate direction
results$direction <- ifelse(results$fc.brca > 0, "Up", "Down")

# Volcano plot
    # Definition
    data <- results
    x <- data$fc.brca
    y <- -log10(data$p.brca)
    z <- data$direction
    label <- data$genes
    maintitle <- "TCGA"
    subtitle <- "BRCA2-deficient vs. BRCA-proficient breast cancers"
    xtitle <- "log2(FC)"
    ytitle <- "-log10(P value)"
    ymin <- floor(min(y))
    ymax <- ceiling(max(y))
    yinterval <- 2
    xmin <- max(abs(floor(min(x))), ceiling(max(x))) * -1
    xmax <- max(abs(floor(min(x))), ceiling(max(x)))
    xinterval <- 1
    
    # Plot
    plot <- ggplot(data, aes(x=x, y=y, label=label)) +
               geom_point(shape=20, size=4, mapping=aes(color=z), show.legend=FALSE) +
               geom_text_repel(aes(label=label)) +
               geom_vline(xintercept = 0) +
               geom_hline(yintercept = 1.3) +
               scale_colour_manual(values=c("blue", "red")) +
               scale_x_continuous(limits=c(xmin, xmax), breaks=seq(xmin, xmax, by=xinterval)) +
               scale_y_continuous(limits=c(ymin, ymax), breaks=seq(ymin, ymax, by=yinterval)) +
               labs(title=maintitle, subtitle=subtitle, x=xtitle, y=ytitle) +
               theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     panel.background=element_blank(),
                     axis.line = element_line(colour="black"),
                     axis.text.x=element_text(colour="black"),
                     axis.text=element_text(size=12),
                     axis.title=element_text(size=16),
                     plot.title=element_text(size=16),
                     plot.subtitle=element_text(size=12))
                     
    # Sava plot
    ggsave("Figures/Figure_2H.pdf", plot, width=8, height=8)
