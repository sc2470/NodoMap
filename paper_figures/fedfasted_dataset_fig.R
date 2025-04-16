#Bar plot of fed and fasted cells from inhouse
##Already subset inhouse
library(ggplot2)

df <- data.frame(table(data@meta.data$Nutr.cond))
df$Var1 <- c("Adlib", "Fast")

pdf("Figure_fedfasted_barplot.pdf", height = 3, width = 3)
ggplot(df, aes(x = Var1, y = Freq)) +
  geom_bar(stat="identity", fill = c('#117733', '#AA4499'), width = 0.5) + theme_classic() +
  xlab(NULL) +
  ylab("Single Cell / Nucleus") +
  theme(axis.text = element_text(size = 20, color = "black"),
        axis.title = element_text(size = 20, vjust = -4)) 

dev.off()

