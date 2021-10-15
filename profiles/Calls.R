df<-read.csv("splatt_256procs_med5tns.csv",header=TRUE,skip=5)
library(ggplot2)

#g<-ggplot(data=df, aes(x=dose, y=len, fill=supp)) +
#  geom_bar(stat="identity")

df <- df[,colSums(df !=0) > 0]
df <- df[df$Comm != "NULL", ]

library(tidyverse)
df <- df %>% select(ends_with(c("_Calls", "Comm")))


#str(df)

library(reshape)

library(RColorBrewer)
#library(viridis)

df2 <- melt(df)
str(df2)
colnames(df2) <- c("Communicator", "Call", "NCalls")
#heatmap(as.matrix(df2[,-1]))
library(RColorBrewer)
library(svglite)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
g <- ggplot(df2, aes(x = Communicator, y = Call, fill = NCalls)) +
  geom_tile() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_distiller(palette = 18,direction = 1)
  #scale_fill_viridis(direction=1)
#scale_fill_gradient(low="yellow",high="black",space="Lab")
#scale_fill_gradient(low = "#ff2D00",high = "#ffffff",guide = "colorbar")
ggsave(file="calls_splat.svg", plot=g, width=10, height=8)
print(g)
