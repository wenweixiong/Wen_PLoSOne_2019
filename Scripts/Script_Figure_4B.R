# Load packages
library(ggplot2)
library(scales)

# Set working directory
setwd("/Users/WIMM/Documents/BRCA_2018/Github/")

# Read file
df <- read.table("Data/TCGA_final.txt", header=TRUE, sep="\t", na.strings="NA", stringsAsFactors=FALSE, check.names=FALSE)

# Subset relevant columns
df <- df[, c("BRCA_status", "PTEN_Status", "PTEN_RPPA", "T.Cell.Inflammed.Signature.Score")]

# Remove missing values
df <- na.omit(df)

# Set factor levels
df$PTEN_Status <- as.factor(df$PTEN_Status)
df$PTEN_Status <- factor(df$PTEN_Status, levels=c("No CNA or point mutations", "Amplification", "Deep Deletion", "Truncating mutation (putative driver)", "Missense Mutation (putative driver)", "Missense Mutation (putative passenger)", "Inframe Mutation (putative passenger)"))

# Define PTEN MT, WT
df$PTEN_Status <- ifelse(df$PTEN_Status=="No CNA or point mutations" | df$PTEN_Status=="Amplification", "Wildtype", "Mutant")

# Set factor levels
df$PTEN_Status <- as.factor(df$PTEN_Status)
df$PTEN_Status <- factor(df$PTEN_Status, levels=c("Wildtype", "Mutant"))

# Fisher exact test
    # BRCA1
    df.small <- df[which(df$BRCA_status %in% c("BRCA-proficient", "BRCA1-deficient")),]
    tbl <- table(df.small$PTEN_Status, df.small$BRCA_status)
    fisher.test(as.matrix(tbl))$p.value

    # BRCA2
    df.small <- df[which(df$BRCA_status %in% c("BRCA-proficient", "BRCA2-deficient")),]
    tbl <- table(df.small$PTEN_Status, df.small$BRCA_status)
    fisher.test(as.matrix(tbl))$p.value

# Tabulate sample size per genotype
n.proficient <- length(df$BRCA_status[which(df$BRCA_status=="BRCA-proficient")])
n.brca1 <- length(df$BRCA_status[which(df$BRCA_status=="BRCA1-deficient")])
n.brca2 <- length(df$BRCA_status[which(df$BRCA_status=="BRCA2-deficient")])

# Barplot
    # Definition
    data <- as.data.frame(table(df$BRCA_status, df$PTEN_Status))
    x <- data[,1]
    y <- data$Freq
    z <- data[,2]
    maintitle <- ""
    xtitle <- ""
    ytitle <- "Samples (%)"
    xticks <- c(paste("BRCA-\nproficient\nbreast\ncancers\n(n=", n.proficient, ")", sep=""),
                paste("BRCA1-\ndeficient\nbreast\ncancers\n(n=", n.brca1, ")", sep=""),
                paste("BRCA2-\ndeficient\nbreast\ncancers\n(n=", n.brca2, ")", sep=""))
    legend.title <- "PTEN mutation status"
    legend.label <- c("Mutant", "Wildtype")

    # Plot
    plot <- ggplot(data, aes(x=x, y=y, fill=z)) +
        geom_bar(position="fill", stat="identity") +
        scale_y_continuous(labels=percent_format()) +
        scale_fill_manual(values=c("blue", "red"), name=legend.title, breaks=legend.label) +
        scale_x_discrete(labels=xticks) +
        labs(title=NULL, x=NULL, y=ytitle) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
              axis.line=element_line(colour="black"),
              axis.text=element_text(size=13),
              axis.title=element_text(size=18),
              axis.text.x=element_text(colour="black"),
              plot.title=element_text(size=18))

    # Save plot
    ggsave("Figures/Figure_4B.pdf", plot, width=8, height=7)
    

