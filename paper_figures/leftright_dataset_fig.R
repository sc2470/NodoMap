#Bar plot of left and right cells from Buchanan and inhouse
library(ggplot2)
df <- data.frame(table(data@meta.data$Position))

df$Var1 <- c("Left","Right")

pdf("figure_leftright_barplot.pdf", width = 3, height = 3)
ggplot(df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = c("yellowgreen","plum3"), width = 0.5) + theme_classic() +
  xlab(NULL) +
  ylab("Single Cell / Nucleus") +
  theme(axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 20, vjust = -4)) 
dev.off()






