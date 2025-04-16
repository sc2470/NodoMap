#QC1
#Aim: draw a line plot with the average of asw at each resolution from each PCs
#Calculate the overall average and SD of the asw at different PCs and resolutions

setwd("/home/rstudio/2023_April_Integration/asw_2023_Jun")

name <- c("Run", "ASW", "PC_Res")
data <- list.files(path = ".", pattern = "[0-9].asw.txt")
##convert the list to data frame with rows and column combined + the file names
data <- setNames(do.call(rbind, Map("cbind", lapply(data, read.delim), V3 = data)), name)
##calculate the average of asw
p <- aggregate(data$ASW, list(data$PC_Res), FUN=mean)
colnames(p)[2] <- "mean"
##calculate the standard deviation of asw avgerage
q <- aggregate(data$ASW, list(data$PC_Res), FUN=sd)
p$sd <- q$x[match(p$Group.1, q$Group.1)]

#string pattern example: "nodose-integrated-pcs-30-resolution-0.6.asw.txt"
#remove everything up to and including pcs-, and then everything from the subsequent -
p$pc <- gsub(".*pcs-|-.*", "", p$Group.1)
#remove everything up to and including resolution-, and then everything from the subsequent .asw
p$res <- gsub(".*resolution-|.asw.*", "", p$Group.1)

write.table(p, "asw_avg_PC30-60_Res0.6-0.8.txt", sep = "\t", quote = F, col.names = NA)


#Line plot
data <- read.delim("asw_avg_PC30-60_Res0.6-0.8.txt")

library(ggplot2)

data$res <- as.factor(data$res)
data$pc <- as.factor(data$pc)

pdf("overall_asw_line_plot.pdf", width = 4, height = 4)
ggplot(data, aes(x=res, y=mean, group=pc)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color=pc), width=.1, size=1) +
  geom_line(aes(color=pc), size=1) +
  geom_point(aes(color=pc, shape=pc), size=2) +
  theme_classic() +
  xlab("Resolutions") +
  ylab("Average Silhouette Width") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black", family = "Helvetica"),
        axis.title = element_text(face = "bold", size = 10, family = "Helvetica"),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(face = "bold", size =10),
        legend.title = element_text(face = "bold", size = 10))
dev.off()








