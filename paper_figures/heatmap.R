library(dplyr)
library(pheatmap)
library(viridis)
library(RColorBrewer)

#Plotting the top10 -log(B-H p-value) of each cluster

fedfast <- c("RHO GTPase cycle",
             "Estrogen receptor signalling",
             "mTOR signalling",
             "Epithelial adherens junction signalling",
             "D-myo-inositol-5-phosphate metabolism",
             "Sirtuin signalling",
             "3-phosphoinositide degradation",
             "3-phosphoinositide biosynthesis",
             "TP53 Regulates Metabolic Genes",
             "NOTCH1 signalling",
             "ER translocation",
             "Translation elongation",
             "Translation initiation",
             "Translation termination",
             "Selenoamino acid metabolism",
             "Nonsense-mediated decay (NMD)",
             "Response to amino acid deficiency",
             "EIF2 signalling",
             "Nucleolar and cytosolic rRNA processing",
             #"Sirtuin Signaling Pathway",
             #"Eukaryotic Translation Initiation",
             #"SRP-dependent cotranslational protein targeting to membrane",
             #"Sirtuin Signaling Pathway"
             "Thermogenesis by uncoupling proteins",
             #"EIF2 Signaling",
             #"Nonsense-Mediated Decay (NMD)",
             #"Major pathway of rRNA processing in the nucleolus and cytosol",
             "Protein ubiquitination",
             "Pre-mRNA processing",
             "Deubiquitination",
             "MHC class I antigen activities",
             "Autophagy",
             "Cardiac hypertrophy signalling",
             #"Epithelial Adherens Junction Signaling",
             "NGF Signalling",
             "AMPK Signalling",
             "Mitotic telophase/cytokinesis",
             #"Protein Ubiquitination Pathway",
             "Semaphorin neuronal repulsive signalling",
             "Renin-angiotensin signalling")

NGN2 <- c(9.67,3.61,3.43,3.23,2.98,2.91,2.87,2.69,2.68,2.68,
          0.687,1.48,0.949,0.788,1.21,1.26,1.07,0.945,0.542,
          2.25,1.43,1.42,1.69,
          1.21,2.63,1.69,0.766,1.66,0.424,0.736,0.152)

NGN6 <- c(1.33,2.57,6.1,0.689,0.148,9.49,0.157,0.123,1.45,0.341,
          24.1,23.9,23.6,23.1,22.3,22.2,21.9,16.6,16.3,
          4.48,8.97,1.77,5.6,
          2.11,0.617,0.155,0.81,0.367,NA,0,0.249)

NGN9 <- c(0.389,1.54,1.16,0.0837,0.748,5.41,0.789,1.24,0.101,0.146,
          5.41,3.28,5.85,3.3,2.73,4.88,2.92,5.26,4.76,
          5.32,3.8,3.38,3.33,
          1.38,0.376,0,0.118,1.2,NA,0.161,0.118)

NGN11 <- c(1.93,1.45,0.765,2.44,0.472,0.665,0.5,0.665,1.18,0.614,
           0.191,0.166,0,0,0,0,0,0.656,0.402,
           0.871,2.15,1.74,1.17,
           3.51,2.76,2.76,2.34,2.28,2.15,2.15,2.15)

data <- data.frame(
  row.names=fedfast,
  NGN2=NGN2,
  NGN6=NGN6,
  NGN9=NGN9,
  NGN11=NGN11)

data <- t(data)

pdf("IPA_heatmap_fedfast.pdf", width = 15, height = 7)
pheatmap(data, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlGn")))(100),
         cellwidth = 25,
         cellheight = 25,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         scale = "column",
         angle_col = 315,
         fontfamily = "Helvetica",
         fontsize = 20)
dev.off()


leftright <- c("Thermogenesis by uncoupling proteins",
               "Mitochondrial dysfunction",
               "Oxidative phosphorylation",
               "Synaptogenesis signalling",
               "RHO GTPase cycle",
               "Estrogen receptor signalling",
               "Clathrin-mediated endocytosis",
               "Glutaminergic Receptor signalling",
               "MHC class I antigen activities",
               "Clathrin-mediated endocytosis recognition",
               #"Synaptogenesis Signaling Pathway",
               #"Glutaminergic Receptor Signaling Pathway (Enhanced)",
               "Opioid signalling",
               "L1CAM interactions",
               "Translation termination",
               "NFAT in cardiac hypertrophy",
               #"Estrogen Receptor Signaling",
               "Neurexins and neuroligins",
               "Nonsense-mediated decay (NMD)",
               "mTOR signalling",
               #"L1CAM interactions",
               #"RHO GTPase cycle",
               "Protein ubiquitination",
               "CLEAR signalling",
               #"Synaptogenesis Signaling Pathway",
               "MHC class II antigen presentation",
               #"Glutaminergic Receptor Signaling Pathway (Enhanced)",
               "COPI-mediated anterograde transport",
               "Golgi transport and retrograde ER traffic",
               #"Role of NFAT in Cardiac Hypertrophy",
               #"Class I MHC mediated antigen processing and presentation",
               #"L1CAM interactions",
               #"Protein Ubiquitination Pathway"
               #"Synaptogenesis Signaling Pathway",
               #"Glutaminergic Receptor Signaling Pathway (Enhanced)",
               #"RHO GTPase cycle",
               "Netrin signalling",
               "Huntington's Disease Signalling",
               #"Intra-Golgi and retrograde Golgi-to-ER traffic",
               "PTEN Regulation")

NGN2 <- c(6.91,6.91,6.91,6.91,5.32,5.14,5.06,5,4.95,4.36,
          2.56,1.47,1.36,2.65,2.28,2.62,1.44,
          3.17,1.2,0.757,1.77,1.9,
          3.1,3.83,1.02)

NGN5 <- c(1.04,5.71,1.66,18.6,4.01,9.81,9,12.8,0.608,5.27,
          12.8,10.6,10.3,10.2,9.25,9.19,9.12,
          1.81,5.39,4.9,3.23,4.9,
          7.7,7.7,2.61)

NGN10 <- c(2.11,9.5,2.73,12.1,13.6,7.41,5.05,11.6,5.88,2.24,
           8.13,18,0,10.6,3.77,0.32,4.07,
           13.6,12.2,12.1,11.6,11.3,
           5.93,6.72,1.21)

NGN17 <- c(2.43,5,1.98,8.49,7.92,3.41,6.08,8.49,11.5,3.77,
           4.16,9.09,0.419,4.3,4.07,0.994,4.12,
           8.93,4.49,4.81,5.99,6.97,
           7.13,7.07,6.82)

data <- data.frame(
  row.names=leftright,
  NGN2=NGN2,
  NGN5=NGN5,
  NGN10=NGN10,
  NGN17=NGN17)

pdf("IPA_heatmap_leftright.pdf", width = 10, height = 12)
pheatmap(data,
         cellwidth = 25,
         cellheight = 25,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",
         fontfamily = "Helvetica",
         fontsize = 20)

dev.off()







