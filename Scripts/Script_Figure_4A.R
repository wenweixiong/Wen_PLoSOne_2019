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

# Set factor levels
df$PTEN_Status <- as.factor(df$PTEN_Status)
df$PTEN_Status <- factor(df$PTEN_Status, levels=c("No CNA or point mutations", "Amplification", "Deep Deletion", "Truncating mutation (putative driver)", "Missense Mutation (putative driver)", "Missense Mutation (putative passenger)", "Inframe Mutation (putative passenger)"))

# Compute p-values
    # Define groups
    nonref <- c("Amplification", "Deep Deletion", "Truncating mutation (putative driver)", "Missense Mutation (putative driver)", "Missense Mutation (putative passenger)", "Inframe Mutation (putative passenger)")
    ref <- "No CNA or point mutations"

    # Compute
    p.value <- NULL

    for(i in 1:length(nonref)) {

        x <- df$PTEN_RPPA[which(df$PTEN_Status==ref)]
        y <- df$PTEN_RPPA[which(df$PTEN_Status==nonref[i])]
        p.value[i] <- wilcox.test(x, y, alternative="two.sided", paired=FALSE)$p.value

    }

    # Tabulate
    data.frame("Comparison"=paste(nonref, " vs. No CNA or point mutation", sep=""), p.value)

# Tabulate sample size per group
freq <- as.data.frame(table(df$PTEN_Status))$Freq

# Violin plot
    # Definition
    data <- df
    x <- data$PTEN_Status
    y <- data$PTEN_RPPA
    maintitle <- ""
    xtitle <- ""
    ytitle <- "PTEN expression"
    xticks <- c(paste("No CNA or\npoint mutations\n(n=", freq[1], ")", sep=""),
                paste("Amplification\n(n=", freq[2], ")", sep=""),
                paste("Deep deletion\n(n=", freq[3], ")", sep=""),
                paste("Truncating mutation\n(putative driver)\n(n=", freq[4], ")", sep=""),
                paste("Missense mutation\n(putative driver)\n(n=", freq[5], ")", sep=""),
                paste("Missense mutation\n(putative passenger)\n(n=", freq[6], ")", sep=""),
                paste("Inframe mutation\n(putative passenger)\n(n=", freq[7], ")", sep=""))

    # Plot
    plot <- ggplot(data, aes(x=x, y=y)) +
        geom_violin(aes(fill=x), show.legend=FALSE, color="transparent", alpha=0.5) +
        geom_jitter(height=0, width=0.05, color="black", size=1) +
        scale_x_discrete(labels=xticks) +
        stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median,
        geom="crossbar", width=0.5, fatten=1, color="red") +
        labs(title=NULL, x=NULL, y=ytitle) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
              axis.line = element_line(colour="black"),
              axis.text=element_text(size=13),
              axis.title=element_text(size=18),
              axis.text.x=element_text(colour="black"))

    # Save plot
    ggsave("Figures/Figure_4A.pdf", plot, width=14, height=7)
