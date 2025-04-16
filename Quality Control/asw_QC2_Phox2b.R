#QC2: Use the nodose ganglia cell cluster specific asw to draw line graph to identify the best pc&res model
#Phox2b is the marker of nodose ganglia cells
#Identify the clusters where nodose ganglia cells are clustered ("Phox2b_cluster")
#Extract the pergroup_asw of Phox2b expressed clusters and average their asw from run 1-5

setwd("/home/rstudio/2023_April_Integration/asw_2023_Jun")

#pc30, res0.6
phox2b <- read.delim("/home/rstudio/2023_April_Integration/Phox2b/Phox2b_pc30_res0.6.txt")

asw <- list.files(path = ".", pattern = "nodose-integrated-pcs-30-resolution-0.6run[1-5].sw.txt_pergroup_asw.txt")
asw <- lapply(asw, read.delim)

#Extract the asw of phox2b cluster in five runs
phox2b_asw <- list()
for (i in 1:5) {
  phox2b_asw[[i]] <-asw[[i]]$x[match(phox2b$x, asw[[i]]$Group.1)]
}
#The average phox2b asw of each run
phox2b_avg_asw <- sapply(phox2b_asw, mean)
#The overall average phox2b sw at pc30 res0.6
phox2b_overall_asw <- c(mean(phox2b_avg_asw), sd(phox2b_avg_asw))




#pc30
res.list<-seq(0.6, 2.0, 0.2)
phox2b_asw <- list()
phox2b_asw_pc30 <- NULL

for (i in res.list) {
  phox2b <- read.delim(paste0("/home/rstudio/2023_April_Integration/Phox2b/Phox2b_pc30_res",i,".txt"))
  asw <- list.files(path = ".", pattern = paste0("nodose-integrated-pcs-30-resolution-",i,"run[1-5].sw.txt_pergroup_asw.txt"))
  asw <- lapply(asw, read.delim)
  for (j in 1:5) {
    phox2b_asw[[j]] <-asw[[j]]$x[match(phox2b$x, asw[[j]]$Group.1)]
  }
  avg_asw <- sapply(phox2b_asw, mean)
  overall_asw <- data.frame(avg=mean(avg_asw), sd=sd(avg_asw), row.names = i)
  phox2b_asw_pc30 <- rbind(phox2b_asw_pc30, overall_asw)
}

phox2b_asw_pc30$pc <- 30
write.table(phox2b_asw_pc30, file = "phox2b_asw_pc30.txt", sep = "\t", quote = F)



#pc40
res.list<-seq(0.6, 2.0, 0.2)
phox2b_asw <- list()
phox2b_asw_pc40 <- NULL


for (i in res.list) {
  phox2b <- read.delim(paste0("/home/rstudio/2023_April_Integration/Phox2b/Phox2b_pc40_res",i,".txt"))
  asw <- list.files(path = ".", pattern = paste0("nodose-integrated-pcs-40-resolution-",i,"run[1-5].sw.txt_pergroup_asw.txt"))
  asw <- lapply(asw, read.delim)
  for (j in 1:5) {
    phox2b_asw[[j]] <-asw[[j]]$x[match(phox2b$x, asw[[j]]$Group.1)]
  }
  avg_asw <- sapply(phox2b_asw, mean)
  overall_asw <- data.frame(avg=mean(avg_asw), sd=sd(avg_asw), row.names = i)
  phox2b_asw_pc40 <- rbind(phox2b_asw_pc40, overall_asw)
}

phox2b_asw_pc40$pc <- 40
write.table(phox2b_asw_pc40, file = "phox2b_asw_pc40.txt", sep = "\t", quote = F)



#pc50
res.list<-seq(0.6, 2.0, 0.2)
phox2b_asw <- list()
phox2b_asw_pc50 <- NULL


for (i in res.list) {
  phox2b <- read.delim(paste0("/home/rstudio/2023_April_Integration/Phox2b/Phox2b_pc50_res",i,".txt"))
  asw <- list.files(path = ".", pattern = paste0("nodose-integrated-pcs-50-resolution-",i,"run[1-5].sw.txt_pergroup_asw.txt"))
  asw <- lapply(asw, read.delim)
  for (j in 1:5) {
    phox2b_asw[[j]] <-asw[[j]]$x[match(phox2b$x, asw[[j]]$Group.1)]
  }
  avg_asw <- sapply(phox2b_asw, mean)
  overall_asw <- data.frame(avg=mean(avg_asw), sd=sd(avg_asw), row.names = i)
  phox2b_asw_pc50 <- rbind(phox2b_asw_pc50, overall_asw)
}

phox2b_asw_pc50$pc <- 50
write.table(phox2b_asw_pc50, file = "phox2b_asw_pc50.txt", sep = "\t", quote = F)



#pc60
res.list<-seq(0.6, 2.0, 0.2)
phox2b_asw <- list()
phox2b_asw_pc60 <- NULL


for (i in res.list) {
  phox2b <- read.delim(paste0("/home/rstudio/2023_April_Integration/Phox2b/Phox2b_pc60_res",i,".txt"))
  asw <- list.files(path = ".", pattern = paste0("nodose-integrated-pcs-60-resolution-",i,"run[1-5].sw.txt_pergroup_asw.txt"))
  asw <- lapply(asw, read.delim)
  for (j in 1:5) {
    phox2b_asw[[j]] <-asw[[j]]$x[match(phox2b$x, asw[[j]]$Group.1)]
  }
  avg_asw <- sapply(phox2b_asw, mean)
  overall_asw <- data.frame(avg=mean(avg_asw), sd=sd(avg_asw), row.names = i)
  phox2b_asw_pc60 <- rbind(phox2b_asw_pc60, overall_asw)
}

phox2b_asw_pc60$pc <- 60
write.table(phox2b_asw_pc60, file = "phox2b_asw_pc60.txt", sep = "\t", quote = F)


#Line plot
temp <- list.files(path = ".", pattern = "phox2b_asw_pc")

phox2b <- do.call(rbind, lapply(temp, read.delim))
phox2b$res <- seq(0.6, 2.0, 0.2)

rownames(phox2b) <- unlist(lapply(temp, row.names))

library(ggplot2)

phox2b$res <- as.factor(phox2b$res)
phox2b$pc <- as.factor(phox2b$pc)

pdf("phox2b_asw_line_plot.pdf", width = 4, height = 4)
ggplot(phox2b, aes(x=res, y=avg, group=pc)) +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd, color=pc), width=.1, size=1) +
  geom_line(aes(color=pc), size=1) +
  geom_point(aes(color=pc, shape=pc), size=2) +
  theme_classic() +
  xlab("Resolutions") +
  ylab("Average Silhouette Width of Phox2b") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black", family = "Helvetica"),
        axis.title = element_text(face = "bold", size = 10, family = "Helvetica"),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 10))
dev.off()


#Line plot with the average of QC1 and QC2

#Calculate the average asw at resolution from 0.6 to 2
p <- aggregate(data$mean, list(data$res), FUN=mean)
q <- aggregate(data$mean, list(data$res), FUN=sd)

p$sd <- q$x[match(p$Group.1, q$Group.1)]
p$asw <- "Overall"


p1 <- aggregate(phox2b$avg, list(phox2b$res), FUN=mean)
q1 <- aggregate(phox2b$avg, list(phox2b$res), FUN=sd)

p1$sd <- q$x[match(p1$Group.1, q1$Group.1)]
p1$asw <- "Phox2b"

#Combine two dataframe
avg_res <- rbind(p,p1)

colnames(avg_res)[1] <- "res"
colnames(avg_res)[2] <- "avg"

#Line plot

pdf("avg_asw_line_plot.pdf", width = 25, height = 25)

ggplot(avg_res, aes(x=res, y=avg, group=asw)) +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd, color=asw), width=.1, size=1) +
  geom_line(aes(color=asw), size=3) +
  geom_point(aes(color=asw, shape=asw), size=5) +
  theme_classic() +
  xlab("Resolutions") +
  ylab("The avgerage of average silhouette width") +
  theme(axis.text = element_text(face = "bold", size = 40, color = "black"),
        axis.title = element_text(face = "bold", size = 40),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2),
        legend.text = element_text(face = "bold", size = 40),
        legend.title = element_text(face = "bold", size = 40),
        legend.position = c(0.8,0.8),
        legend.spacing = unit(2.0, 'cm'),
        legend.key.size = unit(1.4, "cm")) 

dev.off()


