
raw_water<-read.csv("/your_own_directory/water_chemicals.tsv", sep = "\t", header = TRUE)

water <- t(raw_water[,-c(1,2)])
colnames(water) <- raw_water[,1]

hc_water<-hclust(dist(water), method = "average")

# install.packages("ggdendro")
library(ggdendro)
library(ggplot2)

df_water<-dendro_data(hc_water,type = "rectangle")
df1_water<-df_water$segments
df2_water<-df_water$labels


pdf("/your_own_directory/water_chemicals.pdf")
ggplot()+
  geom_segment(data=df1_water,aes(x=x,y=y,
                                  xend=xend,
                                  yend=yend))+
  geom_text(data=df2_water,aes(x=x,y=y-20,label=label),size=3)+
  geom_point(data=df2_water,aes(x=x,y=y))+
  scale_y_continuous(expand = c(0.1,0))+
  theme_void()
dev.off()
