rna_norm_counts <- read.csv("your own directory")
sample_sheet <- read.csv("your own directory")


library(limma)

data_X_1_10 <- rna_norm_counts[,-c(1:7)]
rownames(data_X_1_10) <- rna_norm_counts[,1]

group_X_1_10 <- c(rep("X1",72), rep("X10",72))
design_X_1_10 <- model.matrix(~0+factor(group_X_1_10))
colnames(design_X_1_10) <- c("Tumour","Normal")
rownames(design_X_1_10) <- colnames(data_X_1_10)

# Create expression and design matrices to construct linear regression models
contrast_X_1_10 <- makeContrasts(Tumour - Normal, levels = design_X_1_10)
# Compare the linear regression model to the design matrix
fit_X_1_10 <- lmFit(data_X_1_10, design_X_1_10)
# Compare the linear regression model to the contrast matrix
fit2_X_1_10 <- contrasts.fit(fit_X_1_10, contrast_X_1_10)
# Bayesian test
fit3_X_1_10 <- eBayes(fit2_X_1_10)
# Extracting a table named DEG_X_1_10 of top-ranked genes from a linear model fit
DEG_X_1_10 <- topTable(fit3_X_1_10, coef = "Tumour - Normal", n = Inf)

# Up-regulation and down-regulation of genes
library(dplyr)
dif_X_1_10 <- DEG_X_1_10 %>% mutate(change = case_when( # Adding a new column, and naming it "change"
  logFC >= 1 & adj.P.Val <= 0.05 ~"UP",
  logFC <= -1 & adj.P.Val <= 0.05 ~"DOWN",
  TRUE ~"UNCHANGED"
))
