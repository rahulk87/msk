source("ccf.R")
suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("RColorBrewer"));
suppressPackageStartupMessages(library("plyr"));
suppressPackageStartupMessages(library("dplyr"));
suppressPackageStartupMessages(library("tidyr"));
#suppressPackageStartupMessages(library("stringr"));
suppressPackageStartupMessages(library("magrittr"));
suppressPackageStartupMessages(library("foreach"));
suppressPackageStartupMessages(library("VariantAnnotation"));

file<-read.csv("test.txt", sep="\t", header=T)

a=1
for(a in 1:length(file$ID)){
	#print(file$lcn[a])

	tcn <- file$tcn[a]
	lcn <- file$lcn[a]
	vaf <- file$TUMOR_MAF1
	purity <- file$Purities

	ccfFit <- computeCCF(vaf = vaf, tcn, lcn, purity = purity)
	ccf <- ccfFit$ccf
	print (ccf)
	}

