rna_norm_counts <- read.csv("your own directory/rna_norm_counts.csv")
sample_sheet <- read.csv("your own directory/sample_sheet.csv")


labels_gene_raw <- sample_sheet$REF
labels_gene <- rep("Control",length(labels_gene_raw))

for (i in 1:length(labels_gene_raw)) {
  if (labels_gene_raw[i] == "1x") {
    labels_gene[i] <- "X1"
  }
  if (labels_gene_raw[i] == "10x") {
    labels_gene[i] <- "X10"
  }
  if (labels_gene_raw[i] == "Control") {
    labels_gene[i] <- "Control"
  }
}

data_gene <- t(rna_norm_counts[,-1])
colnames(data_gene) <- NULL
rownames(data_gene) <- NULL
#rownames(data_gene) <- labels_gene


pca_result_gene <- prcomp(data_gene, scale. = TRUE)


pca_data_gene <- as.data.frame(pca_result_gene$x)
pca_data_gene$labels <- labels_gene_raw



library(ggplot2)
library(ggsci)
pdf("your own directory/gene_pca_2D.pdf")
ggplot(pca_data_gene, aes(PC1, PC2, color = labels)) + 
  geom_point() +
  stat_ellipse(aes(fill = labels), alpha = 0.2,
               geom ="polygon",type = "norm")+
  scale_fill_aaas()+
  scale_color_aaas()+
  theme_bw()+
  labs(title = "2D PCA with Labels")
dev.off()





library(scatterplot3d)
pdf("your own directory/gene_pca_3D.pdf")
scatterplot3d(pca_data_gene[,1:3], 
              color = c(rep("#00AFBB", 6), rep("#E7B800", 72), rep("#FC4E07", 72)),
              pch = 15,
              lty.hide = 2
)
legend("topleft",c('Control','x1','x2'),
       fill=c("#00AFBB", "#E7B800", "#FC4E07"),box.col=NA)
dev.off()

