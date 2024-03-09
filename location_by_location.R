
rna_norm_counts <- read.csv("your own directory/rna_norm_counts.csv")
sample_sheet <- read.csv("your own directory/sample_sheet.csv")

control <- rna_norm_counts[,2:7]
rownames(control) <- rna_norm_counts[,1]

L01X1 <- rna_norm_counts[,8:13]
rownames(L01X1) <- rna_norm_counts[,1]

L02X1 <- rna_norm_counts[,14:19]
rownames(L02X1) <- rna_norm_counts[,1]

L03X1 <- rna_norm_counts[,20:25]
rownames(L03X1) <- rna_norm_counts[,1]

L04X1 <- rna_norm_counts[,26:31]
rownames(L04X1) <- rna_norm_counts[,1]

L05X1 <- rna_norm_counts[,32:37]
rownames(L05X1) <- rna_norm_counts[,1]

L06X1 <- rna_norm_counts[,38:43]
rownames(L06X1) <- rna_norm_counts[,1]

L07X1 <- rna_norm_counts[,44:49]
rownames(L07X1) <- rna_norm_counts[,1]

L08X1 <- rna_norm_counts[,50:55]
rownames(L08X1) <- rna_norm_counts[,1]

L09X1 <- rna_norm_counts[,56:61]
rownames(L09X1) <- rna_norm_counts[,1]

L10X1 <- rna_norm_counts[,62:67]
rownames(L10X1) <- rna_norm_counts[,1]

L11X1 <- rna_norm_counts[,68:73]
rownames(L11X1) <- rna_norm_counts[,1]

L12X1 <- rna_norm_counts[,74:79]
rownames(L12X1) <- rna_norm_counts[,1]



my_limma <- function(a, b) {
  data <- data.frame(a, b)
  group <- c(rep("control", ncol(a)), rep("treatment", ncol(b)))
  
  library(limma)
  design <- model.matrix(~0+factor(group))
  colnames(design) <- c("control","treatment")
  rownames(design) <- colnames(data)
  
  # Create expression and design matrices to construct linear regression models
  contrast <- makeContrasts(control - treatment, levels = design)
  # Compare the linear regression model to the design matrix
  fit <- lmFit(data, design)
  # Compare the linear regression model to the contrast matrix
  fit2 <- contrasts.fit(fit, contrast)
  # Bayesian test
  fit3 <- eBayes(fit2)
  # Extracting a table named DEG_X_1_10 of top-ranked genes from a linear model fit
  DEG <- topTable(fit3, coef = "control - treatment", n = Inf)
  
  # Up-regulation and down-regulation of genes
  library(dplyr)
  dif <- DEG %>% mutate(change = case_when( # Adding a new column, and naming it "change"
    logFC >= 1 & adj.P.Val <= 0.05 ~"UP",
    logFC <= -1 & adj.P.Val <= 0.05 ~"DOWN",
    TRUE ~"UNCHANGED"
  ))
  
  return(dif)
}


dif_control_L01X1 <- my_limma(control, L01X1)

dif_control_L02X1 <- my_limma(control, L02X1)

dif_control_L03X1 <- my_limma(control, L03X1)

dif_control_L04X1 <- my_limma(control, L04X1)

dif_control_L05X1 <- my_limma(control, L05X1)

dif_control_L06X1 <- my_limma(control, L06X1)

dif_control_L07X1 <- my_limma(control, L07X1)

dif_control_L08X1 <- my_limma(control, L08X1)

dif_control_L09X1 <- my_limma(control, L09X1)

dif_control_L10X1 <- my_limma(control, L10X1)

dif_control_L11X1 <- my_limma(control, L11X1)

dif_control_L12X1 <- my_limma(control, L12X1)




L01X10 <- rna_norm_counts[,80:85]
rownames(L01X10) <- rna_norm_counts[,1]

L02X10 <- rna_norm_counts[,86:91]
rownames(L02X10) <- rna_norm_counts[,1]

L03X10 <- rna_norm_counts[,92:97]
rownames(L03X10) <- rna_norm_counts[,1]

L04X10 <- rna_norm_counts[,98:103]
rownames(L04X10) <- rna_norm_counts[,1]

L05X10 <- rna_norm_counts[,104:109]
rownames(L05X10) <- rna_norm_counts[,1]

L06X10 <- rna_norm_counts[,110:115]
rownames(L06X10) <- rna_norm_counts[,1]

L07X10 <- rna_norm_counts[,116:121]
rownames(L07X10) <- rna_norm_counts[,1]

L08X10 <- rna_norm_counts[,122:127]
rownames(L08X10) <- rna_norm_counts[,1]

L09X10 <- rna_norm_counts[,128:133]
rownames(L09X10) <- rna_norm_counts[,1]

L10X10 <- rna_norm_counts[,134:139]
rownames(L10X10) <- rna_norm_counts[,1]

L11X10 <- rna_norm_counts[,140:145]
rownames(L11X10) <- rna_norm_counts[,1]

L12X10 <- rna_norm_counts[,146:151]
rownames(L12X10) <- rna_norm_counts[,1]


dif_control_L01X10 <- my_limma(control, L01X10)

dif_control_L02X10 <- my_limma(control, L02X10)

dif_control_L03X10 <- my_limma(control, L03X10)

dif_control_L04X10 <- my_limma(control, L04X10)

dif_control_L05X10 <- my_limma(control, L05X10)

dif_control_L06X10 <- my_limma(control, L06X10)

dif_control_L07X10 <- my_limma(control, L07X10)

dif_control_L08X10 <- my_limma(control, L08X10)

dif_control_L09X10 <- my_limma(control, L09X10)

dif_control_L10X10 <- my_limma(control, L10X10)

dif_control_L11X10 <- my_limma(control, L11X10)

dif_control_L12X10 <- my_limma(control, L12X10)


