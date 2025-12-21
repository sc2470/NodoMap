
# fedfast_marker_supp_table -----------------------------------------------

df <- fedfastmarkers
map <- c("0"=	"EC1	Cldn5/Ly6c1/Ly6a",
         "1"=	"SGC1	Mmd2/Ednrb/Fabp7",
         "2"=	"SGC2	Cxcl10/Acsbg1/Phgdh",
         "3"=	"MGC1	Ncmap/Fam178b/Cldn19",
         "4"=	"SGC3	Bcan/Fbln2/Fbln5",
         "5"=	"NGN1	C1ql3/Adcyap1/St18",
         "6"=	"HC1	C1qa/C1qc/C1qb",
         "7"=	"SGC4	Atf3/Emp1/Cebpd",
         "8"=	"SGC5	Pou3f1/Col1a2/Gas7",
         "9"=	"FB1	Lum/Ccl11/Smoc2",
         "10"=	"SGC6	Ifit3/Cdh19/Isg15",
         "11"=	"MGC2	Prx/Pllp/Bcas1",
         "12"=	"FB2	Acta2/Myh11/Ndufa4l2",
         "13"=	"NGN2	Htr3b/Htr3a/Gpr65",
         "14"=	"GC1	Kcna1/Fxyd3/Mbp",
         "15"=	"NGN3	Cysltr2/Chrnb3/Kcnip4",
         "16"=	"MGC3	Slc36a2/Nr4a2/Ugt8a",
         "17"=	"NGN4	Uts2b/Cckar/Vip",
         "18"=	"SGC7	G0s2/Sdc4/Ccn1",
         "19"=	"MGC4	Mt2/Tcim/Drp2",
         "21"=	"NGN5	Kcng1/Vmn1r85/Trpv1",
         "22"=	"SGC8	Ttyh1/Ptprz1/Sfrp5",
         "23"=	"JGN1	Mrgprd/Tmem45b/Grik1",
         "24"=	"NGN6	Lypd1/Ddc/Tafa2",
         "25"=	"GC2	Sostdc1/Plekhb1/Tsc22d4",
         "26"=	"JGN2	Gfra3/Tac1/Tmem255a",
         "27"=	"FB3	Igfbp6/Thbs4/Islr",
         "28"=	"NGN7	Rbp4/Gda/Kcnk9",
         "29"=	"NGN8	Trpa1/Efcab6/Nos1",
         "30"=	"EC2	Selp/Fabp4/Vwf",
         "31"=	"NGN9	Miat/Snhg11/Meg3",
         "32"=	"NGN10	Olfm3/Tmem233/P2ry1",
         "33"=	"NGN11	Sprr1a/Ecel1/Ucn",
         "34"=	"NGN12	Slc18a3/Pappa2/Runx3",
         "35"=	"EC3	Podxl/Depp1/Ptprb",
         "36"=	"NGN13	Lox/Pkib/F2r",
         "37"=	"NGN14	Slc17a7/Hapln4/Lgi3",
         "38"=	"JGN3	Trappc3l/Nptx1/Tafa1",
         "39"=	"HC2	Cd52/Lsp1/Rac2",
         "40"=	"NGN15	Bmp3/Gal/Pcdh9",
         "41"=	"JGN4	Wfdc2/C1ql4/Rarres1",
         "42"=	"NGN16	Gabra1/Chodl/Gabrb2",
         "43"=	"NGN17	Thsd7b/Cacng5/Gata3",
         "44"=	"FB4	Hbb-bs/Bnc2/Ebf2",
         "45"=	"GC3	Tyrobp/Lyz2/Fcer1g",
         "46"=	"NGN18	Chrm2/Col24a1/Lrp1b",
         "47"=	"JGN5	Trpm8/Foxp2/Cdh8",
         "48"=	"MGC5	Smoc2/Ccl11/Lum",
         "49"=	"MGC6	Ifitm1/Kcnj8/Rgs5",
         "50"=	"NGN19	Slc6a2/Npy/Hand2",
         "51"=	"NGN20	Glp1r/Amigo2/Npy2r",
         "52"=	"NGN21	Olfr78/Runx3/Avpr1a")

df$col3 <- map[as.character(fedfastmarkers$cluster)]
df$cluster_class <- sub("\\s.*","",df$col3)
df$cluster_label <- sub("^\\S+\\s+","",df$col3)
df <- df[,-9]

write.table(df, print(paste0("fedfasted_marker_merged", ".txt")), sep='\t', quote=F)



# leftright_marker_supp_table ---------------------------------------------------------------

df1 <- leftrightmarkers
df1$col3 <- map[as.character(leftrightmarkers$cluster)]
df1$cluster_class <- sub("\\s.*","",df1$col3)
df1$cluster_label <- sub("^\\S+\\s+","",df1$col3)
df1 <- df1[,-9]

write.table(df1, print(paste0("leftright_marker_merged", ".txt")), sep='\t', quote=F)









