#  Remember to read in the file of rna counts into rna_norm_counts!!!!!!!!

rna_norm_counts <- read.csv("/your_directory/rna_norm_counts.tsv", sep = "\t")

control <- rna_norm_counts[,2:7]
rownames(control) <- rna_norm_counts[,1]

experiment <- rna_norm_counts[,8:ncol(rna_norm_counts)]
rownames(experiment) <- rna_norm_counts[,1]

data_ec <- data.frame(experiment, control)

group_ec <- c(rep("experiment", ncol(experiment)), rep("control", ncol(control)))

library(limma)
design_ec <- model.matrix(~0+factor(group_ec))
colnames(design_ec) <- c("experiment","control")
rownames(design_ec) <- colnames(data_ec)

# Create expression and design matrices to construct linear regression models
contrast_ec <- makeContrasts(experiment - control, levels = design_ec)
# Compare the linear regression model to the design matrix
fit_ec <- lmFit(data_ec, design_ec)
# Compare the linear regression model to the contrast matrix
fit2_ec <- contrasts.fit(fit_ec, contrast_ec)
# Bayesian test
fit3_ec <- eBayes(fit2_ec)
# Extracting a table named DEG_X_1_10 of top-ranked genes from a linear model fit
DEG_ec <- topTable(fit3_ec, coef = "experiment - control", n = Inf)

# Up-regulation and down-regulation of genes
library(dplyr)
dif_ec <- DEG_ec %>% mutate(change = case_when( # Adding a new column, and naming it "change"
  logFC >= 1 & adj.P.Val <= 0.05 ~"UP",
  logFC <= -1 & adj.P.Val <= 0.05 ~"DOWN",
  TRUE ~"UNCHANGED"
))

sig_genes <- rownames(dif_ec[dif_ec$adj.P.Val < 0.05, ])

write.csv(dif_ec, "differential_treatment_vs_control.csv", row.names = TRUE)
