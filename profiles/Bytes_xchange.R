df<-read.csv("splatt_128procs_med5tns.csv",header=TRUE,skip=5)
library(ggplot2)
p<-ggplot(data=df[df$Send_Calls > 0, ], aes(x=Comm, y=Send_Calls, fill=Comm)) +
  geom_bar(stat="identity")
#g<-ggplot(data=df, aes(x=dose, y=len, fill=supp)) +
#  geom_bar(stat="identity")

df <- df[df$Comm != "NULL", ]

library(tidyverse)
df <- df %>% select(ends_with(c("_Bytes", "Comm")))


str(df)

library(reshape)

library(RColorBrewer)

df2 <- melt(df)
#str(as.matrix(df2[-1,]))
colnames(df2) <- c("Communicator", "Call", "Bytes")
g <- ggplot(df2, aes(x = Communicator, y = Call, fill = Bytes)) +
geom_tile() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  scale_fill_gradient(low = "#ff0000",
                    high = "#ffffff",
                    guide = "colorbar")
print(g)
