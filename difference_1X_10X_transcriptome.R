rna_norm_counts <- read.csv("your own directory/rna_norm_counts.csv")
sample_sheet <- read.csv("your own directory/sample_sheet.csv")


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


# In order to keep the consistence above, I set the threshold of logFC between -1 and +1
right_X_1_10 <- 1
left_X_1_10 <- (-right_X_1_10)

# To make the matrices more suitable for plotting, by changing some elements
dif1_X_1_10 <- dif_X_1_10
colnames(dif1_X_1_10)[1] = "log2FC"
Changes <- ifelse((dif1_X_1_10$P.Value < 0.05 & abs(dif1_X_1_10$log2FC)> abs(left_X_1_10)), ifelse(dif1_X_1_10$log2FC > right_X_1_10,"UP","DOWN"), "UNCHANGED")

# Do a minus log10 process with the p-values
dif1_X_1_10[,4] <- (-log10(dif1_X_1_10[,4]))
colnames(dif1_X_1_10)[4] = "minuslog10PValue"

# Prepare the labels for the genes which are highly concerned
a_X_1_10<-subset(dif1_X_1_10,change=="UP")
a1_X_1_10<-head(a_X_1_10[order(abs(a_X_1_10[,4]),decreasing=TRUE),],5)
b_X_1_10<-subset(dif1_X_1_10,change=="DOWN")
b1_X_1_10<-head(b_X_1_10[order(abs(b_X_1_10[,4]),decreasing=TRUE),],5)
data1_X_1_10<-rbind(a1_X_1_10,b1_X_1_10)


# Load the package of "ggplot"
library(ggplot2)

# In case that ggrepel has not been installed:
# install.packages("ggrepel")
# Load the package of ggrepel to label the points of concerned genes
library(ggrepel)
# Plotting the volcano map with labels
pdf("/rds/homes/t/txj322/00000M5_Group/plots/gene_volcano.pdf")
ggplot(dif1_X_1_10, aes(log2FC, minuslog10PValue))+
  geom_point(aes(col=Changes), size=0.1)+
  scale_color_manual(values=c("blue","grey","red"))+
  geom_vline(xintercept=c(left_X_1_10,right_X_1_10), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  theme(plot.width = 100, plot.height = 100)+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
  theme_bw()+
  coord_cartesian(xlim = c(-6000, 6000), ylim = c(-1, 25))+
  labs(title = "Volcano Map (with Labels) -- The Appropriate (minus) log10 P-value")+
  labs(x="log2(FoldChange)",y="-log10(P-value)")+
  str(dif1_X_1_10, max.level = c(-1, 1))+
geom_label_repel(
  data = data1_X_1_10, 
  aes(log2FC, minuslog10PValue, label = rownames(data1_X_1_10)),size = 2, max.overlaps = 30)
dev.off()

