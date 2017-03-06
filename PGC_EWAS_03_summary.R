########################################################################################
# Forest Plots and Results Summary
########################################################################################

rm(list=ls())
library(rmeta)
library(ChAMP)
library(forestplot)

########################################################################################
# Step 1: Load and subset the data
########################################################################################

rm(list=ls())

# Checking results for all the sites

load("PGC_EWAS_nonSmoke_inverseNorm_intersectSites.Rdata")
load("PGC_EWAS_nonSmoke_inverseNorm_remainingSites.Rdata")

all<-pval.combined.marot
rm(pval.combined.marot)
colnames(all)<-c("CpG", "p")
Studies<-rep(c("DNHS, GTP, MRS, VA, PRISMO, WTC, DUKE, AS, MIR"), nrow(all))
all<-cbind(all, Studies)
all<-all[, colnames(results)]
all(colnames(all)==colnames(results))
results<-rbind(all, results)

pvalue <- as.numeric(results[, "p"])
chisq <- qchisq(1-pvalue,1)
lambda = median(chisq)/qchisq(0.5,1) 
lambda # 1.15
FDR<-p.adjust(pvalue, method="fdr", n=nrow(results))
sum(FDR<=0.05) # 5

sum(pvalue<=5*10^-5) # 47
sum(pvalue<=5*10^-6) # 12
sum(pvalue<=5*10^-7) # 4

results<-cbind(results, FDR)
save(results, file="PGC_EWAS_inverseNorm_allResults.Rdata")

nStudies<-unlist(lapply(strsplit(results[, "Studies"], ", "), function(x) length(x)))
sum(nStudies==9) # 442966
sum(nStudies==8) # 9273
sum(nStudies==7) # 1020
sum(nStudies==6) # 1862
sum(nStudies==5) # 568
sum(nStudies==4) # 235
sum(nStudies==3) # 137
sum(nStudies==2) # 60
sum(nStudies==1) # 73

rm(list=ls()[-match("results", ls())])

load("PGC_EWAS_DataPrep_nonSmoke.Rdata")

results<-results[as.numeric(results[, "p"])<=5*10^-5, ] # subsetting only most significant results
results<-results[order(as.numeric(results[, "p"])), ]

rownames(results)<-results[, "CpG"]
data(probe.features)
probe.features<-probe.features[results[, "CpG"], c("gene", "feature")]
all(rownames(results)==rownames(probe.features))
probe.features$pvalue<-signif(as.numeric(results[, "p"]), 2)
probe.features$FDR<-signif(as.numeric(results[, "FDR"]), 3)

write.csv(probe.features, "PGC_EWAS_TopResults_09.15.16.csv")

# Subsetting individual site data
results<-results[results[, "Studies"]=="DNHS, GTP, MRS, VA, PRISMO, WTC, DUKE, AS, MIR, INTRUST", ]
results<-results[as.numeric(results[, "FDR"])<=0.05, ]

colnames(DNHS.coef) 
DNHS.coef<-DNHS.coef[rownames(results),c("PTSDpm", "N.subjects")] # just the PTSD coefficient and number of subjects needed
DNHS.results<-DNHS.results[rownames(results),]
colnames(DNHS.coef)<-c("PTSD", "N.subjects") # need consistent column headings
all(rownames(DNHS.coef)==rownames(DNHS.results))
DNHS<-data.frame(DNHS.coef, DNHS.results)
str(DNHS) # should all be numbers nut factors
DNHS$s.e.<-DNHS$PTSD/DNHS$t # calculate the standard error
DNHS$weight<-1/(DNHS$s.e.^2) # calculate weight
str(DNHS) # should all be numbers not factors before converting to matrix
DNHS.m<-as.matrix(DNHS ) # matrices speed up calculations over dataframes
all(DNHS.m==DNHS) 
DNHS<-DNHS.m
rm(DNHS.results, DNHS.m)

colnames(GTP.coef)
GTP.coef<-GTP.coef[rownames(results),c("Current_PTSD_01_sum", "N.subjects")]
GTP.results<-GTP.results[rownames(results),]
colnames(GTP.coef)<-c("PTSD", "N.subjects")
all(rownames(GTP.coef)==rownames(GTP.results))
GTP<-data.frame(GTP.coef, GTP.results)
str(GTP)
GTP$s.e.<-GTP$PTSD/GTP$t # calculate the standard error
GTP$weight<-1/(GTP$s.e.^2) # calculate weight
str(GTP)
GTP.m<-as.matrix(GTP)
all(GTP.m==GTP)
GTP<-GTP.m
rm(GTP.results, GTP.m)

colnames(MRS.coef)
MRS.coef<-MRS.coef[rownames(results),c("PTSDbroad", "N.subjects")]
MRS.results<-MRS.results[rownames(results),]
colnames(MRS.coef)<-c("PTSD", "N.subjects")
all(rownames(MRS.coef)==rownames(MRS.results))
MRS<-data.frame(MRS.coef, MRS.results)
str(MRS)
MRS$s.e.<-MRS$PTSD/MRS$t # calculate standard error
MRS$weight<-1/(MRS$s.e.^2) # calculate weight
str(MRS)
MRS.m<-as.matrix(MRS)
all(MRS.m==MRS)
MRS<-MRS.m
rm(MRS.results, MRS.m)

colnames(VA.coef)
VA.coef<-VA.coef[rownames(results),c("PTSD_C", "N.subjects")]
VA.results<-VA.results[rownames(results),]
colnames(VA.coef)<-c("PTSD", "N.subjects")
all(rownames(VA.coef)==rownames(VA.results))
VA<-data.frame(VA.coef, VA.results)
str(VA)
VA$s.e.<-VA$PTSD/VA$t # calculate the standard error
VA$weight<-1/(VA$s.e.^2) # calculate weight
str(VA)
VA.m<-as.matrix(VA)
all(VA.m==VA)
VA<-VA.m
rm(VA.results, VA.m)

colnames(PRISMO.coef)
PRISMO.coef<-PRISMO.coef[rownames(results),c("ptssTRUE", "N.subjects")]
PRISMO.results<-PRISMO.results[rownames(results),]
colnames(PRISMO.coef)<-c("PTSD", "N.subjects")
all(rownames(PRISMO.coef)==rownames(PRISMO.results))
PRISMO<-data.frame(PRISMO.coef, PRISMO.results)
str(PRISMO)
PRISMO$s.e.<-PRISMO$PTSD/PRISMO$t # calculate the standard error
PRISMO$weight<-1/(PRISMO$s.e.^2) # calculate weight
str(PRISMO)
PRISMO.m<-as.matrix(PRISMO)
all(PRISMO.m==PRISMO)
PRISMO<-PRISMO.m
rm(PRISMO.results, PRISMO.m)

colnames(WTC.coef)
WTC.coef<-WTC.coef[rownames(results),c("PTSD", "N.subjects")]
WTC.results<-WTC.results[rownames(results),]
colnames(WTC.coef)<-c("PTSD", "N.subjects")
all(rownames(WTC.coef)==rownames(WTC.results))
WTC<-data.frame(WTC.coef, WTC.results)
str(WTC)
WTC$s.e.<-WTC$PTSD/WTC$t # calculate the standard error
WTC$weight<-1/(WTC$s.e.^2) # calculate weight
str(WTC)
WTC.m<-as.matrix(WTC)
all(WTC.m==WTC)
WTC<-WTC.m
rm(WTC.results, WTC.m)

colnames(DUKE.coef)
DUKE.coef<-DUKE.coef[rownames(results),c("currptsd", "N.subjects")]
DUKE.results<-DUKE.results[rownames(results),]
colnames(DUKE.coef)<-c("PTSD", "N.subjects")
all(rownames(DUKE.coef)==rownames(DUKE.results))
DUKE<-data.frame(DUKE.coef, DUKE.results)
str(DUKE)
DUKE$s.e.<-DUKE$PTSD/DUKE$t # calculate the standard error
DUKE$weight<-1/(DUKE$s.e.^2) # calculate weight
str(DUKE)
DUKE.m<-as.matrix(DUKE)
all(DUKE.m==DUKE)
DUKE<-DUKE.m
rm(DUKE.results, DUKE.m)

colnames(AS.coef)
AS.coef<-AS.coef[rownames(results),c("d_pts30_t2", "N.subjects")]
AS.results<-AS.results[rownames(results),]
colnames(AS.coef)<-c("PTSD", "N.subjects")
all(rownames(AS.coef)==rownames(AS.results))
AS<-data.frame(AS.coef, AS.results)
str(AS)
AS$s.e.<-AS$PTSD/AS$t # calculate the standard error
AS$weight<-1/(AS$s.e.^2) # calculate weight
str(AS)
AS.m<-as.matrix(AS)
all(AS.m==AS)
AS<-AS.m
rm(AS.coef,AS.results, AS.m)

colnames(MIR.coef)
MIR.coef<-MIR.coef[rownames(results),c("curr_ptsd", "N.subjects")]
MIR.results<-MIR.results[rownames(results),]
colnames(MIR.coef)<-c("PTSD", "N.subjects")
all(rownames(MIR.coef)==rownames(MIR.results))
MIR<-data.frame(MIR.coef, MIR.results)
str(MIR)
MIR$s.e.<-MIR$PTSD/MIR$t # calculate the standard error
MIR$weight<-1/(MIR$s.e.^2) # calculate weight
str(MIR)
MIR.m<-as.matrix(MIR)
all(MIR.m==MIR)
MIR<-MIR.m
rm(MIR.coef,MIR.results, MIR.m)

colnames(TRUST.coef)
TRUST.coef<-TRUST.coef[rownames(results),c("PTSD", "N.subjects")]
TRUST.results<-TRUST.results[rownames(results),]
colnames(TRUST.coef)<-c("PTSD", "N.subjects")
all(rownames(TRUST.coef)==rownames(TRUST.results))
TRUST<-data.frame(TRUST.coef, TRUST.results)
str(TRUST)
TRUST$s.e.<-TRUST$PTSD/TRUST$t # calculate the standard error
TRUST$weight<-1/(TRUST$s.e.^2) # calculate weight
str(TRUST)
TRUST.m<-as.matrix(TRUST)
all(TRUST.m==TRUST)
TRUST<-TRUST.m
rm(TRUST.coef,TRUST.results, TRUST.m)

all(rownames(results)==rownames(DNHS))
all(rownames(results)==rownames(GTP))
all(rownames(results)==rownames(MRS))
all(rownames(results)==rownames(VA))
all(rownames(results)==rownames(PRISMO))
all(rownames(results)==rownames(WTC))
all(rownames(results)==rownames(DUKE))
all(rownames(results)==rownames(AS))
all(rownames(results)==rownames(MIR))
all(rownames(results)==rownames(TRUST))

rm(list=ls()[grep("ebayes", ls())])

save.image("PGC_Meta-Analysis_Data_Limma_DGMBPWDAMI_04_09.15.16.Rdata")

rm(list=ls())

########################################################################################
# Step 2: Generate Plots
########################################################################################

load("PGC_Meta-Analysis_Data_Limma_DGMBPWDAMI_04_09.15.16.Rdata")

FEM<-matrix(nrow=nrow(results), ncol=8)
rownames(FEM)<-rownames(results)
colnames(FEM)<-c("beta_avg", "var_avg", "lower", "upper", "sig", "Q", "Q_pval", "I2")

all(rownames(results)==rownames(DNHS))
all(rownames(results)==rownames(GTP))
all(rownames(results)==rownames(MRS))
all(rownames(results)==rownames(VA))
all(rownames(results)==rownames(WTC))
all(rownames(results)==rownames(PRISMO))
all(rownames(results)==rownames(DUKE))
all(rownames(results)==rownames(AS))
all(rownames(results)==rownames(MIR))
all(rownames(results)==rownames(TRUST))
all(rownames(results)==rownames(FEM))


for(ii in 1:nrow(FEM)){
  weights<-c(DNHS[ii, "weight"], GTP[ii, "weight"], MRS[ii, "weight"], VA[ii, "weight"],
             PRISMO[ii, "weight"], WTC[ii, "weight"], DUKE[ii, "weight"], AS[ii, "weight"],
             MIR[ii, "weight"], TRUST[ii, "weight"])
  betas<-c(DNHS[ii, "PTSD"], GTP[ii, "PTSD"], MRS[ii, "PTSD"], VA[ii, "PTSD"],
           PRISMO[ii, "PTSD"], WTC[ii, "PTSD"], DUKE[ii, "PTSD"], AS[ii, "PTSD"],
           MIR[ii, "PTSD"], TRUST[ii, "PTSD"])
  avg<-sum(betas*weights)/sum(weights)
  FEM[ii, "beta_avg"]<-avg
  var<-1/sum(weights)
  FEM[ii, "var_avg"]<-var
  lower<-avg-(1.96*sqrt(var))
  FEM[ii, "lower"]<-lower
  upper<-avg+(1.96*sqrt(var))
  FEM[ii, "upper"]<-upper
  if(upper > 0 && lower > 0){
    FEM[ii, "sig"]<-"significant"
  } else if(upper< 0 && lower < 0){
    FEM[ii, "sig"]<-"significant"
  } else {
    FEM[ii, "sig"]<-"non-significant"
  }
  Q<-sum(weights*((betas-avg)^2))
  FEM[ii, "Q"]<-Q
  FEM[ii, "Q_pval"]<-pchisq(Q, df=10-1, lower.tail=FALSE) # 10 studies - 1 for degrees of freedom
  I2<-max(0, (Q-(10-1))/Q)
  FEM[ii, "I2"]<-I2
  if(ii%%10000==0){
    print(ii)
  }
}

pdf("PGC_MA_ForestPlots_Limma_DGMBPWDAMI_09.15.16.pdf")
for(ii in 1:nrow(DNHS)){
  betas<-c(DNHS[ii, "PTSD"], GTP[ii, "PTSD"], MRS[ii, "PTSD"], 
           PRISMO[ii, "PTSD"], DUKE[ii, "PTSD"], VA[ii, "PTSD"], WTC[ii, "PTSD"], 
           AS[ii, "PTSD"], MIR[ii, "PTSD"], TRUST[ii, "PTSD"])
  stderr<-c(DNHS[ii, "s.e."], GTP[ii, "s.e."], MRS[ii, "s.e."], 
            PRISMO[ii, "s.e."], DUKE[ii, "s.e."], VA[ii, "s.e."], WTC[ii, "s.e."],
            AS[ii, "s.e."], MIR[ii, "s.e."], TRUST[ii, "s.e."])
  lower<-betas-1.96*stderr
  upper<-betas+1.96*stderr
  weights<-c(DNHS[ii, "weight"], GTP[ii, "weight"], MRS[ii, "weight"], 
             PRISMO[ii, "weight"], DUKE[ii, "weight"], VA[ii, "weight"], WTC[ii, "weight"],
             AS[ii, "weight"], MIR[ii, "weight"], TRUST[ii, "weight"])
  betaC<-sum(betas*weights)/sum(weights)
  stderrC<-1/sum(weights)
  lowerC<-betaC-(1.96*sqrt(stderrC))
  upperC<-betaC+(1.96*sqrt(stderrC))
  beta.table<-c("",
                "Beta", round(betas,3), NA, round(betaC,3))
  betas<-c(NA, NA, round(betas,3), NA, round(betaC,3))
  lower<-c(NA, NA, lower, NA, lowerC)
  upper<-c(NA, NA, upper, NA, upperC)
  studies<-c(rownames(DNHS)[ii], "Study", "DNHS", "GTP", "MRS", "PRISMO", "VA-M", "VA-NCP", "WTC", "AS", 
             "MIR", "TRUST", NA, "Summary")
  summary<-cbind(studies, beta.table)
  # Call forestplot
  cochrane<-data.frame(betas, lower, upper)
  forestplot(summary, cochrane, 
             #boxsize=0.5, lwd.ci=1.5, 
             new_page = TRUE, 
             is.summary=c(TRUE,TRUE,rep(FALSE,10),TRUE),
             xlog=FALSE, #lineheight=unit(1,"cm"),
             col=fpColors(box="black",line="black", summary="royalblue"),
             xticks=(c(-0.4, -0.2, 0.1)),
             txt_gp = fpTxtGp( label = list(gpar(fontfamily="", cex=2)),
                               ticks = gpar(fontfamily = "", cex=1)))
}
dev.off()

########################################################################################
# Step 3: QQ Plots of Combined Results
########################################################################################

rm(list=ls())

load("PGC_EWAS_DGMBPWDAMI_allResults_09.15.16.Rdata")

ggd.qqplot = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=2, cex.lab=4.5, cex.axis=4, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}

pvalue <- as.numeric(results[, "p"])
chisq <- qchisq(1-pvalue,1)
lambda = median(chisq)/qchisq(0.5,1)
lambda # 1.472844

jpeg("PGC_MA_ResultsMarot_QQplot_Limma_DGMBPWDAMI_09.15.16.jpg", width=1279, height=960, quality=100, units="px")
par(mar=c(10,10,4,2), mgp=c(6,2,0))
ggd.qqplot(pvalue,
           main = paste(""))
dev.off()

rm(list=ls())


########################################################################################
# Step 4: Rank Results
########################################################################################

load("PGC_EWAS_DGMBPWDAMI_allResults_09.15.16.Rdata")
load("PGC_EWAS_DataPrep_09.15.16.Rdata")

rank<-matrix(nrow=nrow(results), ncol=10)
rownames(rank)<-rownames(results)
colnames(rank)<-c("AS", "DNHS", "GTP", "MIR", "MRS", "PRISMO", "DUKE", "TRUST", "VA", "WTC")

all(rownames(DNHS.results)==rownames(AS.results))
all(rownames(DNHS.results)==rownames(GTP.results))
all(rownames(DNHS.results)==rownames(DUKE.results))
all(rownames(DNHS.results)==rownames(MIR.results))
all(rownames(DNHS.results)==rownames(MRS.results))
all(rownames(DNHS.results)==rownames(VA.results))
all(rownames(DNHS.results)==rownames(PRISMO.results))
all(rownames(DNHS.results)==rownames(WTC.results))
all(rownames(DNHS.results)==rownames(TRUST.results))

sum(is.na(match(rownames(results), rownames(AS.results))))
sum(is.na(match(rownames(results), rownames(DNHS.results))))
sum(is.na(match(rownames(results), rownames(DUKE.results))))
sum(is.na(match(rownames(results), rownames(GTP.results))))
sum(is.na(match(rownames(results), rownames(MIR.results))))
sum(is.na(match(rownames(results), rownames(MRS.results))))
sum(is.na(match(rownames(results), rownames(VA.results))))
sum(is.na(match(rownames(results), rownames(PRISMO.results))))
sum(is.na(match(rownames(results), rownames(WTC.results))))
sum(is.na(match(rownames(results), rownames(TRUST.results))))

AS.results<-AS.results[rownames(results), ]
DNHS.results<-DNHS.results[rownames(results), ]
DUKE.results<-DUKE.results[rownames(results), ]
GTP.results<-GTP.results[rownames(results), ]
MIR.results<-MIR.results[rownames(results), ]
MRS.results<-MRS.results[rownames(results), ]
VA.results<-VA.results[rownames(results), ]
PRISMO.results<-PRISMO.results[rownames(results), ]
WTC.results<-WTC.results[rownames(results), ]

all(rownames(rank)==rownames(AS.results))
all(rownames(rank)==rownames(DNHS.results))
all(rownames(rank)==rownames(DUKE.results))
all(rownames(rank)==rownames(GTP.results))
all(rownames(rank)==rownames(MIR.results))
all(rownames(rank)==rownames(MRS.results))
all(rownames(rank)==rownames(WTC.results))
all(rownames(rank)==rownames(PRISMO.results))
all(rownames(rank)==rownames(VA.results))

rank[, "AS"]<-AS.results[, "rank"]
rank[, "DNHS"]<-DNHS.results[, "rank"]
rank[, "DUKE"]<-DUKE.results[, "rank"]
rank[, "GTP"]<-GTP.results[, "rank"]
rank[, "MIR"]<-MIR.results[, "rank"]
rank[, "MRS"]<-MRS.results[, "rank"]
rank[, "WTC"]<-WTC.results[, "rank"]
rank[, "PRISMO"]<-PRISMO.results[, "rank"]
rank[, "VA"]<-VA.results[, "rank"]

write.csv(rank, "PGC_MA_resultsMarot_Limma_rank_DGMBPWDAM_04.16.16.csv")

rm(list=ls())

########################################################################################
# Step 6: t-statistic tests for MA results
########################################################################################

# All sites

load("PGC_Meta-Analysis_Data_DGMBPWDAM_Combined_04.16.16.Rdata")

all(rownames(DNHS.results)==rownames(GTP.results))
all(rownames(DNHS.results)==rownames(DUKE.results))
all(rownames(DNHS.results)==rownames(MIR.results))
all(rownames(DNHS.results)==rownames(MRS.results))
all(rownames(DNHS.results)==rownames(WTC.results))
all(rownames(DNHS.results)==rownames(PRISMO.results))
all(rownames(DNHS.results)==rownames(VA.results))
all(rownames(DNHS.results)==rownames(AS.results))

tstats<-cbind(AS.results[, "t"], DNHS.results[, "t"], GTP.results[, "t"], MIR.results[, "t"], MRS.results[, "t"], 
              PRISMO.results[, "t"], DUKE.results[, "t"], VA.results[, "t"],  WTC.results[,"t"])

colnames(tstats)<-c("AS", "DNHS", "GTP", "MIR", "MRS","PRISMO", "VA-D", "VA-NCP", "WTC")
tstatsCor<-cor(tstats)

write.csv(tstatsCor, "PGC_Meta_t-stat_study_Limma_correlations_allSites_DGMBPWDAM_04.16.16.csv")

# Top 50,000 sites
load("PGC_MA_ResultsMarot_DGMBPWDAM_04.16.16.Rdata")
top<-rownames(results.marot)[1:50000]

sum(is.na(match(top, rownames(tstats))))
tstatsTop<-tstats[top, ]
tstatsTopCor<-cor(tstatsTop)

write.csv(tstatsTopCor, "PGC_Meta_t-stat_study_Limma_correlations_top50k_DGMBPWDAM_04.16.16.csv")

# Top results

sum(results.marot[, "pval.combined.marot"]<5*10^-5) # 45
results<-results.marot[results.marot[, "pval.combined.marot"]<5*10^-5,]

tstatsSig<-tstats[rownames(results), ]
tstatsCorSig<-cor(tstatsSig)
tstatsCorSig

write.csv(tstatsCorSig, "PGC_Meta_t-stat_study_Limma_correlations_top20_DGMBPWDAM_04.16.16.csv")
rm(list=ls())


########################################################################################
# Step 7: t-statistic correlation heatmaps
########################################################################################

library(reshape2)
library(ggplot2)
rm(list=ls())
t<-read.csv("PGC_Meta_t-stat_study_Limma_correlations_allSites_DGMBPWDAM_04.16.16.csv", row.names=1)
colnames(t)[match("VA.D", colnames(t))]<-"VA-M"
colnames(t)[match("VA.NCP", colnames(t))]<-"VA-NCP"
m<-melt(t)
m$variable<-as.character(m$variable)
m$var2<-rep(colnames(t), 9)

jpeg("PGC_Limma_heatmap_allSites_04.16.16.jpg", width=640, height=480, units="px")
ggplot(m, aes(x=variable, y= var2))+geom_tile(aes(fill=abs(value)))+
  scale_fill_gradient(low="white",high ="red", name="", 
                      limits=c(-0.1,1), breaks=c(0, 1), guide=FALSE)+
  xlab("")+ylab("")+
  theme(panel.border=element_blank())+theme_bw()+
  theme(text=element_text(size=28), axis.text.x = element_text(angle=90, vjust=1),
        legend.position="bottom")
dev.off()

p<-ggplot(m, aes(x=variable, y= var2))+geom_tile(aes(fill=value))+
  scale_fill_gradient(low="white",high ="red", name="", 
                      limits=c(-0.1,1), breaks=c(0, 0.5, 1))+
  xlab("")+ylab("")+
  theme(panel.border=element_blank())+theme_bw()+
  theme(text=element_text(size=28), axis.text.x = element_text(angle=90, vjust=1),
        legend.position="bottom")

library(gridExtra)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

legend <- g_legend(p)

png("PGC_Limma_heatmap_allSites_04.16.16_legend.jpg", width=480, height=480, units="px", bg="transparent")
grid.arrange(legend)
dev.off()

rm(list=ls())
t<-read.csv("PGC_Meta_t-stat_study_Limma_correlations_top50k_DGMBPWDAM_04.16.16.csv", row.names=1)
colnames(t)[match("VA.D", colnames(t))]<-"VA-M"
colnames(t)[match("VA.NCP", colnames(t))]<-"VA-NCP"
m<-melt(t)
m$variable<-as.character(m$variable)
m$var2<-rep(colnames(t), 9)

jpeg("PGC_Limma_heatmap_50kSites_04.16.16.jpg", width=640, height=480, units="px")
ggplot(m, aes(x=variable, y= var2))+geom_tile(aes(fill=value))+
  scale_fill_gradient(low="white",high ="red", name="Correlations", 
                      limits=c(-0.1,1), breaks=c(0, 1), guide=FALSE)+
  xlab("")+ylab("")+theme(panel.border=element_blank())+theme_bw()+
  theme(text=element_text(size=28), axis.text.x = element_text(angle=90, vjust=1),
        legend.position="bottom")
dev.off()

rm(list=ls())
t<-read.csv("PGC_Meta_t-stat_study_Limma_correlations_top20_DGMBPWDAM_04.16.16.csv", row.names=1)
colnames(t)[match("VA.D", colnames(t))]<-"VA-D"
colnames(t)[match("VA.NCP", colnames(t))]<-"VA-NCP"
m<-melt(t)
m$variable<-as.character(m$variable)
m$var2<-rep(colnames(t), 9)
dm<-m
dm[dm$value==1,"value"]<-0
jpeg("PGC_Limma_heatmap_topSites_04.16.16.jpg", width=480, height=480, units="px")
ggplot(m, aes(x=variable, y= var2))+geom_tile(aes(fill=value))+
  scale_fill_gradient(low="white",high ="red", name="Correlations", 
                      limits=c(-0.1,1), breaks=c(0,  1), guide=F)+
  xlab("")+ylab("")+theme(panel.border=element_blank())+theme_bw()+
  theme(text=element_text(size=28), axis.text.x = element_text(angle=90, vjust=1),
        legend.position="bottom")
dev.off()

rm(list=ls())

########################################################################################
# Step 8: Subjects
########################################################################################

load("PGC_Meta-Analysis_Data_DGMBPWDA_Combined_03.09.16.Rdata")

subjs<-matrix(nrow=nrow(DNHS.results), ncol=8)
rownames(subjs)<-rownames(DNHS.results)
all(rownames(subjs)==rownames(AS.results))
all(rownames(subjs)==rownames(DNHS.results))
all(rownames(subjs)==rownames(GTP.results))
all(rownames(subjs)==rownames(MRS.results))
all(rownames(subjs)==rownames(WTC.results))
all(rownames(subjs)==rownames(PRISMO.results))
all(rownames(subjs)==rownames(VA.results))
all(rownames(subjs)==rownames(DUKE.results))
colnames(subjs)<-c("AS","DNHS", "GTP", "MRS", "WTC", "PRISMO", "VA", "DUKE")

subjs[, "AS"]<-AS.coef[, "N.subjects"]
subjs[, "DNHS"]<-DNHS.coef[, "N.subjects"]
subjs[, "GTP"]<-GTP.coef[, "N.subjects"]
subjs[, "MRS"]<-MRS.coef[, "N.subjects"]
subjs[, "WTC"]<-WTC.coef[, "N.subjects"]
subjs[, "PRISMO"]<-PRISMO.coef[, "N.subjects"]
subjs[, "VA"]<-VA.coef[, "N.subjects"]
subjs[, "DUKE"]<-DUKE.coef[, "N.subjects"]

min(apply(subjs, 1, sum))  # 1155
max(apply(subjs, 1, sum)) # 1229

rm(list=ls())

########################################################################################
# Step 9: Effect Size Scatterplot
########################################################################################

load("/Users/ar3054/Documents/R/PGC_EWAS/Final_04.19.16/Military/PGC_Meta-Analysis_Military_ES_04.16.16.Rdata")
load("/Users/ar3054/Documents/R/PGC_EWAS/Final_02.28.16/Civilian/PGC_Meta-Analysis_Data_DGW_Civilian_ES.Rdata")

cRes<-cbind(names(civRes$TestStatistic), civRes$TestStatistic)
colnames(cRes)<-c("civProbe", "civZ")

mRes<-cbind(names(milRes$TestStatistic), milRes$TestStatistic)
colnames(mRes)<-c("milProbe", "milZ")

sum(is.na(match(rownames(cRes), rownames(mRes)))) # 9,344

sites<-intersect(rownames(cRes), rownames(mRes))

cRes<-cRes[sites,]
mRes<-mRes[sites,]
all(rownames(cRes)==rownames(mRes))
res<-cbind(cRes, mRes)
res<-res[, c("civZ", "milZ")]
class(res)<-"numeric"

df<-data.frame(res)
str(df)

png("PGC_Meta-Analysis_ES_plot_04.16.16.png", width=480, height=480, units="px")
ggplot(df, aes(x=civZ, y=milZ))+geom_point(alpha=0.2)+theme_bw()+
  xlab("Civilian Effect Size")+ylab("Military Effect Size")
dev.off()

cor(df$civZ, df$milZ) #-0.0404349

# Top 50K sites
load("PGC_MA_ResultsMarot_DGMBPWDAM_04.16.16.Rdata")
top<-rownames(results.marot)[1:50000]
df<-df[top,]

png("PGC_Meta-Analysis_ES_top50k_plot_04.16.16.png", width=480, height=480, units="px")
ggplot(df, aes(x=civZ, y=milZ))+geom_point(alpha=0.2)+theme_bw()+
  xlab("Civilian Effect Size")+ylab("Military Effect Size")
dev.off()

cor(df$civZ, df$milZ) # 0.5314942

# Top
top<-rownames(results.marot)[1:5000]
df<-df[top,]

png("PGC_Meta-Analysis_ES_top5k_plot_04.16.16.png", width=480, height=480, units="px")
ggplot(df, aes(x=civZ, y=milZ))+geom_point(alpha=0.2)+theme_bw()+
  xlab("Civilian Effect Size")+ylab("Military Effect Size")
dev.off()

cor(df$civZ, df$milZ) # 0.6875985

# 5x10^-5 sites
head(results.marot)
top<-rownames(results.marot[results.marot[, "pval.combined.marot"]<=5*10^-5, ]) # 45 sites
df<-df[top,]

png("PGC_Meta-Analysis_ES_top45_plot_04.16.16.png", width=480, height=480, units="px")
ggplot(df, aes(x=civZ, y=milZ))+geom_point(alpha=0.2)+theme_bw()+
  xlab("Civilian Effect Size")+ylab("Military Effect Size")
dev.off()

cor(df$civZ, df$milZ) # 0.7358501

