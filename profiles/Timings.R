df<-read.csv("lammps_64procs.csv",header=TRUE,skip=5)
library(ggplot2)

#g<-ggplot(data=df, aes(x=dose, y=len, fill=supp)) +
#  geom_bar(stat="identity")
df <- df[,colSums(df !=0) > 0]
df <- df[df$Comm != "NULL", ]

library(tidyverse)
df <- df %>% select(ends_with(c("_Time", "Comm")))
#df <- df %>% select()

#str(df)

library(reshape)

library(RColorBrewer)
#library(viridis)

df2 <- melt(df)
str(df2)
colnames(df2) <- c("Communicator", "Call", "Time")
#heatmap(as.matrix(df2[,-1]))
#library(RColorBrewer)
#coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
g <- ggplot(df2, aes(x = Communicator, y = Call, fill = Time)) +
  geom_tile() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_distiller(palette = 18,direction = 1)
  #scale_fill_viridis()
#scale_fill_gradient(low="yellow",high="black",space="Lab")
#scale_fill_gradient(low = "#ff2D00",high = "#ffffff",guide = "colorbar")
library(svglite)
ggsave(file="splatt_128procs.svg", plot=g, width=10, height=8)
print(g)
