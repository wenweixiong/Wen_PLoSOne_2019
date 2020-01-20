# Load packages
library(ggplot2)
library(ggrepel)

# Set working directory
setwd("/Users/WIMM/Documents/BRCA_2018/Github/")

# Read file
df <- read.table("Data/WSI_final.txt", header=TRUE, sep="\t", na.strings="NA", check.names=FALSE)

# CIN25 gene list
cells <-  c("ADC", "MDSC", "TEM_CD4", "TEM_CD8", "CD56_DIM", "TREG", "PDC", "MAST", "TFH", "MEM_B_CELL", "IMM_DC", "EOS", "CD56_BRIGHT", "MON", "MAC", "TGD", "ACT_CD4", "NKT", "ACT_B_CELL", "TCM_CD4", "TH1", "NK", "TH17", "IMM_B_CELL", "TH2", "TCM_CD8", "ACT_CD8", "NEU")

cells <- intersect(names(df), cells)

# Subset relevant columns
df <- df[,c("BRCA_status", cells)]

# Remove missing values
df <- na.omit(df)

# Compute FC and p-values
df.brca <- NULL
p.brca <- NULL

for(i in 1:length(cells)) {

    # Subset gene
    df.small <- df[, c("BRCA_status", cells[i])]
    
    # Retrieve expression values
    proficient <- df.small[which(df.small$BRCA_status=="BRCA-proficient"), 2]
    brca <- df.small[which(df.small$BRCA_status=="BRCA2-deficient"), 2]
    
    # Compute FC
    df.brca[i] <- median(brca, na.rm=TRUE) - median(proficient, na.rm=TRUE)
    
    # Compute p.values
    p.brca[i] <- wilcox.test(brca, proficient, na.rm=TRUE)$p.value
    
}

# Tabulate
results <- data.frame(cells, df.brca, p.brca, stringsAsFactors=FALSE)

# Indicate direction
results$direction <- ifelse(results$df.brca > 0, "Up", "Down")

# Volcano plot
    # Definition
    data <- results
    x <- data$df.brca
    y <- -log10(data$p.brca)
    z <- data$direction
    label <- data$cells
    maintitle <- "WSI"
    subtitle <- "BRCA2-deficient vs. BRCA-proficient breast cancers"
    xtitle <- "NES difference"
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
    ggsave("Figures/Figure_3C.pdf", plot, width=8, height=8)

## Latest NES values here differ slightly from published ones because GSEA algorithm uses permutation. Therefore, slightly difference NES values will be returned each round.
