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
sink("output.txt")
for(a in 1:length(file$ID)){
	#print(file$lcn[a])
	tcn <- file$tcn[a]
	lcn <- file$lcn[a]
	vaf <- file$TUMOR_MAF[a]
	vaf1 <- file$TUMOR_MAF1[a]
	purity <- file$Purities[a]
	ref <- file$REF[a]
	alt <- file$ALT[a]

	ccfFit <- computeCCF(vaf = vaf, tcn, lcn, purity = purity)
	ccf <- ccfFit$ccf
	cat (ccf)
	cat("\t")
	conf <- confCCF(alt = alt, ref = ref, tcn, lcn, purity = purity,multiplicity = ccfFit$multiplicity)
	ccfLower <- conf$lower
	ccfUpper <- conf$upper
	clonalStatus <- ifelse(round(ccfLower, 2) >= 0.75, "clonal",
                           ifelse(round(ccfLower, 2) < 0.75 & ccfFit$ccf >= 0.8, 'likely_clonal',
                                  "subclonal"))
	#cat (ccfLower)
        #cat("\t")
	#cat (clonalStatus)
	#cat("\t")

	## Wild type
	ccfFit1 <- computeCCF(vaf = vaf1, tcn, lcn, purity = purity)
        ccf1 <- ccfFit1$ccf
        cat (ccf1)
        #cat("\t")
        conf1 <- confCCF(alt = alt, ref = ref, tcn, lcn, purity = purity,multiplicity = ccfFit1$multiplicity)
        ccfLower1 <- conf1$lower
        ccfUpper1 <- conf1$upper
        clonalStatus1 <- ifelse(round(ccfLower1, 2) >= 0.75, "clonal",
                           ifelse(round(ccfLower1, 2) < 0.75 & ccfFit1$ccf >= 0.8, 'likely_clonal',
                                  "subclonal"))
        #cat (ccfLower1)
        #cat("\t")
	#cat (clonalStatus1)
        cat("\n")
	}
sink()
