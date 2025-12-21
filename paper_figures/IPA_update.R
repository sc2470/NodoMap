library(dplyr)
library(pheatmap)
library(viridis)
library(RColorBrewer)

#Plotting the top10 -log(p-value) of NGN3, NGN9, NGN16, NGN18
leftright <- c("Oxidative Phosphorylation",
               "Mitochondrial Dysfunction",
               "Respiratory electron transport",
               "RHO GTPase cycle",
               "Synaptogenesis signalling Pathway", 
               "Clathrin-mediated endocytosis",
               "Estrogen Receptor signalling",
               "Class I MHC mediated antigen processing and presentation",
               "Signalling by ROBO receptors",
               "Glutaminergic receptor signalling pathway (Enhanced)", 
               "L1CAM interactions",
               #"RHO GTPase cycle", 
               "Protein ubiquitination pathway", 
               #"Synaptogenesis signalling Pathway",
               "CLEAR signalling Pathway", 
               #"Glutaminergic Receptor signalling Pathway (Enhanced)", 
               "MHC class II antigen presentation", 
               "Intra-Golgi and retrograde Golgi-to-ER traffic", 
               "COPI-mediated anterograde transport", 
               "Role of NFAT in Cardiac Hypertrophy",
               #"Class I MHC mediated antigen processing and presentation",
               #"Synaptogenesis signalling Pathway",
               #"RHO GTPase cycle",
               #"Protein Ubiquitination Pathway",
               #"Glutaminergic Receptor signalling Pathway (Enhanced)",
               #"L1CAM interactions",
               "Insulin secretion signalling pathway",
               "Huntington's Disease signalling",
               "Netrin signalling",
               #"Intra-Golgi and retrograde Golgi-to-ER traffic",
               #"RHO GTPase cycle",
               #"Synaptogenesis signalling Pathway",
               "Response of EIF2AK4 (GCN2) to amino acid deficiency",
               #"Glutaminergic Receptor signalling Pathway (Enhanced)",
               "Serotonin receptor signalling",
               #"signalling by ROBO receptors",
               #"Class I MHC mediated antigen processing and presentation",
               #"Clathrin-mediated endocytosis",
               #"Netrin signalling",
               "Nonsense-Mediated Decay (NMD)")

NGN3 <- c(12.7,12.3,11.9,11.5,11.4,8.83,
          8.74,8.54,8.48,8.41,3.29,5.66,
          2.82,1.87,3.51,3.18,4.77,6.33,
          6.67,5.51,2.46,3.46,3.8)

NGN9 <- c(3.98,14.4,4.4,21.2,17.3,7.76,
          10.8,9.27,2.63,16,24.2,18.7,
          16.8,15.6,15.1,15.1,14.8,9.28,
          9.6,8.91,0.329,7.52,0.634)

NGN16 <- c(3.71,8.86,4.03,13.1,13.2,9.07,
           5.78,15.8,6.26,12.7,12.5,12.9,
           7.35,6.65,10.3,8.04,6.71,10.8,
           10.7,10.5,2.85,6.07,2.45)

NGN18 <- c(1.3,5.04,1.68,12.4,10.2,8.34,
           6.65,8.56,8.7,9.26,2.36,4.44,
           5.94,3.09,5.58,5.41,6.6,3.51,
           2.44,8.08,9.39,9.12,8)

data <- data.frame(
  row.names=leftright,
  NGN3=NGN3,
  NGN9=NGN9,
  NGN16=NGN16,
  NGN18=NGN18)

pdf("IPA_heatmap_leftright.pdf", width = 12, height = 12)
pheatmap(data,
         cellwidth = 25,
         cellheight = 25,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",
         fontfamily = "Helvetica",
         fontsize = 15)

dev.off()


# ------------------------------------------------------------------------


# #Plotting the top10 -log(p-value) of NGN2, NGN5, NGN8, NGN10--------------
fedfast <- c("RHO GTPase cycle",
             "Estrogen Receptor Signalling",
             "mTOR Signalling",
             "Epithelial Adherens Junction Signalling",
             "D-myo-inositol-5-phosphate Metabolism",
             "Sirtuin Signalling Pathway",
             "3-phosphoinositide Degradation",
             "3-phosphoinositide Biosynthesis",
             "TP53 Regulates Metabolic Genes",
             "Signalling by NOTCH1",
             "SRP-dependent cotranslational protein targeting to membrane",
             "Eukaryotic Translation Elongation",
             "Eukaryotic Translation Initiation",
             "Eukaryotic Translation Termination",
             "Selenoamino acid metabolism",
             "Nonsense-Mediated Decay (NMD)",
             "Response of EIF2AK4 (GCN2) to amino acid deficiency",
             "EIF2 Signalling",
             "Major pathway of rRNA processing in the nucleolus and cytosol",
             #"Sirtuin Signalling Pathway",
             #"SRP-dependent cotranslational protein targeting to membrane"
             #"Sirtuin Signalling Pathway",
             "Electron transport, ATP synthesis, and heat production by uncoupling proteins",
             #"EIF2 Signalling",
             #"Nonsense-Mediated Decay (NMD)",
             #"Major pathway of rRNA processing in the nucleolus and cytosol"
             "Protein Ubiquitination Pathway",
             "Processing of Capped Intron-Containing Pre-mRNA",
             "Deubiquitination",
             #"Eukaryotic Translation Termination",
             "Class I MHC mediated antigen processing and presentation",
             "Autophagy",
             "Cardiac Hypertrophy Signalling",
             #"Epithelial Adherens Junction Signalling",
             "NGF Signalling",
             "AMPK Signalling",
             "Mitotic Telophase/Cytokinesis",
             #"Protein Ubiquitination Pathway",
             "Semaphorin Neuronal Repulsive Signalling Pathway",
             "Renin-Angiotensin Signalling")

NGN2 <- c(9.67,3.61,3.43,3.23,2.98,2.91,2.87,2.69,
          2.68,2.68,0.687,1.48,0.949,0.788,1.21,1.26,
          1.07,0.945,0.542,2.25,1.43,1.42,1.69,1.21,
          2.63,1.69,0.766,1.66,0.424,0.736,0.152)

NGN5 <- c(1.33,2.57,6.1,0.689,0.148,9.49,0.157,0.123,
          1.45,0.341,24.1,23.9,23.6,23.1,22.3,22.2,
          21.9,16.6,16.3,4.48,8.97,1.77,5.6,2.11,
          0.617,0.155,0.81,0.367,NA,0,0.249)

NGN8 <- c(0.389,1.54,1.16,0.0837,0.748,5.41,0.789,1.24,
          0.101,0.146,5.41,3.28,5.85,3.3,2.73,4.88,
          2.92,5.26,4.76,5.32,3.8,3.38,3.33,1.38,
          0.376,0,0.118,1.2,NA,0.161,0.118)

NGN10 <- c(1.93,1.45,0.765,2.44,0.472,0.665,0.5,0.665,
           1.18,0.614,0.191,0.166,0,0,0,0,
           0,0.656,0.402,0.871,2.15,1.74,1.17,3.51,
           2.76,2.76,2.34,2.28,2.15,2.15,2.15)

data <- data.frame(
  row.names=fedfast,
  NGN2=NGN2,
  NGN5=NGN5,
  NGN8=NGN8,
  NGN10=NGN10)

data <- t(data)

pdf("IPA_heatmap_fedfast.pdf", width = 20, height = 10)
pheatmap(data, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlGn")))(100),
         cellwidth = 25,
         cellheight = 25,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         scale = "column",
         angle_col = 315,
         fontfamily = "Helvetica",
         fontsize = 15)
dev.off()
















