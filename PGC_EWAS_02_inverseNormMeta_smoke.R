########################################################################################
# PGC Meta-Analysis: CpG Sites in All Studies
########################################################################################

rm(list=ls())

########################################################################################
# Step 1: Load and subset data
########################################################################################

load("PGC_EWAS_combinedData_smoke.Rdata")

# The minimum adjusted p-values are:
min(DNHS.results[,"adj.P.Val"]) #  
min(GTP.results[,"adj.P.Val"]) # 
min(MRS.results[,"adj.P.Val"]) # 
min(VA.results[,"adj.P.Val"]) # 
min(PRISMO.results[,"adj.P.Val"]) # 
min(WTC.results[,"adj.P.Val"]) #  
min(DUKE.results[, "adj.P.Val"]) # 
min(AS.results[, "adj.P.Val"]) # 
min(MIR.results[, "adj.P.Val"]) #

# Step 1A: Convert data.frames to matrices. This speeds up the later loops
str(DNHS.results) 
DNHS.results.m<-as.matrix(DNHS.results)
all(DNHS.results==DNHS.results.m)
DNHS.results<-DNHS.results.m
rm(DNHS.results.m)

str(GTP.results)
GTP.results.m<-as.matrix(GTP.results)
all(GTP.results==GTP.results.m)
GTP.results<-GTP.results.m
rm(GTP.results.m)

str(MRS.results)
MRS.results.m<-as.matrix(MRS.results)
all(MRS.results==MRS.results.m)
MRS.results<-MRS.results.m
rm(MRS.results.m)

str(VA.results)
VA.results.m<-as.matrix(VA.results)
all(VA.results==VA.results.m)
VA.results<-VA.results.m
rm(VA.results.m)

str(PRISMO.results)
PRISMO.results.m<-as.matrix(PRISMO.results)
all(PRISMO.results==PRISMO.results.m)
PRISMO.results<-PRISMO.results.m
rm(PRISMO.results.m)

str(WTC.results)
WTC.results.m<-as.matrix(WTC.results)
all(WTC.results==WTC.results.m)
WTC.results<-WTC.results.m
rm(WTC.results.m)

str(DUKE.results)
DUKE.results.m<-as.matrix(DUKE.results)
all(DUKE.results==DUKE.results.m)
DUKE.results<-DUKE.results.m
rm(DUKE.results.m)

str(AS.results)
AS.results.m<-as.matrix(AS.results)
all(AS.results==AS.results.m)
AS.results<-AS.results.m
rm(AS.results.m)

str(MIR.results)
MIR.results.m<-as.matrix(MIR.results)
all(MIR.results==MIR.results.m)
MIR.results<-MIR.results.m
rm(MIR.results.m)

# Step 1B: Check that all the rownames 
sum(is.na(match(rownames(DNHS.coef), rownames(DNHS.results))))
sum(is.na(match(rownames(DNHS.coef), rownames(DNHS.ebayes))))
head(DNHS.results)
DNHS.results<-DNHS.results[rownames(DNHS.ebayes),]
head(DNHS.results)
head(DNHS.coef)
DNHS.coef<-DNHS.coef[rownames(DNHS.ebayes),]
head(DNHS.coef)
all(rownames(DNHS.results)==rownames(DNHS.coef))
all(rownames(DNHS.results)==rownames(DNHS.ebayes))

sum(is.na(match(rownames(GTP.coef), rownames(GTP.results))))
sum(is.na(match(rownames(GTP.coef), rownames(GTP.ebayes))))
head(GTP.results)
GTP.results<-GTP.results[rownames(GTP.ebayes),]
head(GTP.results)
head(GTP.coef)
GTP.coef<-GTP.coef[rownames(GTP.ebayes),]
head(GTP.coef)
all(rownames(GTP.results)==rownames(GTP.coef))
all(rownames(GTP.results)==rownames(GTP.ebayes))

sum(is.na(match(rownames(MRS.coef), rownames(MRS.results))))
sum(is.na(match(rownames(MRS.coef), rownames(MRS.ebayes))))
head(MRS.results)
MRS.results<-MRS.results[rownames(MRS.ebayes),]
head(MRS.results)
head(MRS.coef)
MRS.coef<-MRS.coef[rownames(MRS.ebayes),]
head(MRS.coef)
all(rownames(MRS.results)==rownames(MRS.coef))
all(rownames(MRS.results)==rownames(MRS.ebayes))

sum(is.na(match(rownames(VA.coef), rownames(VA.results))))
sum(is.na(match(rownames(VA.coef), rownames(VA.ebayes))))
head(VA.results)
VA.results<-VA.results[rownames(VA.ebayes),]
head(VA.results)
head(VA.coef)
VA.coef<-VA.coef[rownames(VA.ebayes),]
head(VA.coef)
all(rownames(VA.results)==rownames(VA.coef))
all(rownames(VA.results)==rownames(VA.ebayes))

sum(is.na(match(rownames(PRISMO.coef), rownames(PRISMO.results))))
sum(is.na(match(rownames(PRISMO.coef), rownames(PRISMO.ebayes))))
head(PRISMO.results)
PRISMO.results<-PRISMO.results[rownames(PRISMO.ebayes),]
head(PRISMO.results)
head(PRISMO.coef)
PRISMO.coef<-PRISMO.coef[rownames(PRISMO.ebayes),]
head(PRISMO.coef)
all(rownames(PRISMO.results)==rownames(PRISMO.coef))
all(rownames(PRISMO.results)==rownames(PRISMO.ebayes))

sum(is.na(match(rownames(WTC.coef), rownames(WTC.results))))
sum(is.na(match(rownames(WTC.coef), rownames(WTC.ebayes))))
head(WTC.results)
WTC.results<-WTC.results[rownames(WTC.ebayes),]
head(WTC.results)
head(WTC.coef)
WTC.coef<-WTC.coef[rownames(WTC.ebayes),]
head(WTC.coef)
all(rownames(WTC.results)==rownames(WTC.coef))
all(rownames(WTC.results)==rownames(WTC.ebayes))

sum(is.na(match(rownames(DUKE.coef), rownames(DUKE.results))))
sum(is.na(match(rownames(DUKE.coef), rownames(DUKE.ebayes))))
head(DUKE.results)
DUKE.results<-DUKE.results[rownames(DUKE.ebayes),]
head(DUKE.results)
head(DUKE.coef)
DUKE.coef<-DUKE.coef[rownames(DUKE.ebayes),]
head(DUKE.coef)
all(rownames(DUKE.results)==rownames(DUKE.coef))
all(rownames(DUKE.results)==rownames(DUKE.ebayes))

sum(is.na(match(rownames(AS.coef), rownames(AS.results))))
sum(is.na(match(rownames(AS.coef), rownames(AS.ebayes))))
head(AS.results)
AS.results<-AS.results[rownames(AS.ebayes),]
head(AS.results)
head(AS.coef)
AS.coef<-AS.coef[rownames(AS.ebayes),]
head(AS.coef)
all(rownames(AS.results)==rownames(AS.coef))
all(rownames(AS.results)==rownames(AS.ebayes))

sum(is.na(match(rownames(MIR.coef), rownames(MIR.results))))
sum(is.na(match(rownames(MIR.coef), rownames(MIR.ebayes))))
head(MIR.results)
MIR.results<-MIR.results[rownames(MIR.ebayes),]
head(MIR.results)
head(MIR.coef)
MIR.coef<-MIR.coef[rownames(MIR.ebayes),]
head(MIR.coef)
all(rownames(MIR.results)==rownames(MIR.coef))
all(rownames(MIR.results)==rownames(MIR.ebayes))

# Step 1C: Determine the number of probes available in all studies
dnhs.sites<-rownames(DNHS.coef)
gtp.sites<-rownames(GTP.coef)
mrs.sites<-rownames(MRS.coef)
va.sites<-rownames(VA.coef)
prismo.sites<-rownames(PRISMO.coef)
wtc.sites<-rownames(WTC.coef)
duke.sites<-rownames(DUKE.coef)
as.sites<-rownames(AS.coef)
mir.sites<-rownames(MIR.coef)

sites<-append(dnhs.sites, gtp.sites)
sites<-append(sites, mrs.sites)
sites<-append(sites, va.sites)
sites<-append(sites, prismo.sites)
sites<-append(sites, wtc.sites)
sites<-append(sites, duke.sites)
sites<-append(sites, as.sites)
sites<-append(sites, mir.sites)
sites<-unique(sites)
length(sites) 

all<-intersect(dnhs.sites, gtp.sites)
all<-intersect(all, mrs.sites)
all<-intersect(all, va.sites)
all<-intersect(all, prismo.sites)
all<-intersect(all, wtc.sites)
all<-intersect(all, duke.sites)
all<-intersect(all, as.sites)
all<-intersect(all, mir.sites)

sum(is.na(match(all, rownames(DNHS.coef)))) # these should all be 0 
sum(is.na(match(all, rownames(GTP.coef))))
sum(is.na(match(all, rownames(MRS.coef))))
sum(is.na(match(all, rownames(VA.coef))))
sum(is.na(match(all, rownames(PRISMO.coef))))
sum(is.na(match(all, rownames(WTC.coef))))
sum(is.na(match(all, rownames(DUKE.coef))))
sum(is.na(match(all, rownames(AS.coef))))
sum(is.na(match(all, rownames(MIR.coef))))
length(all) 

rm(dnhs.sites, gtp.sites, mrs.sites, va.sites, wtc.sites, prismo.sites, duke.sites, 
   as.sites, mir.sites)

########################################################################################
# Step 2: Data Prep
########################################################################################

# Step 2A: Calculate one-sided p-values for each CpG site's t-statistic
DNHS.oneSided<-pt(DNHS.ebayes$t, df = (DNHS.ebayes$df.prior+DNHS.ebayes$df.residual))
all(rownames(DNHS.oneSided)==rownames(DNHS.results))
all(rownames(DNHS.oneSided)==rownames(DNHS.coef))
all(rownames(DNHS.oneSided)==rownames(DNHS.ebayes))

GTP.oneSided<-pt(GTP.ebayes$t, df = (GTP.ebayes$df.prior+GTP.ebayes$df.residual))
all(rownames(GTP.oneSided)==rownames(GTP.results))
all(rownames(GTP.oneSided)==rownames(GTP.coef))
all(rownames(GTP.oneSided)==rownames(GTP.ebayes))

MRS.oneSided<-pt(MRS.ebayes$t, df = (MRS.ebayes$df.prior+MRS.ebayes$df.residual))
all(rownames(MRS.oneSided)==rownames(MRS.results))
all(rownames(MRS.oneSided)==rownames(MRS.coef))
all(rownames(MRS.oneSided)==rownames(MRS.ebayes))

VA.oneSided<-pt(VA.ebayes$t, df = (VA.ebayes$df.prior+VA.ebayes$df.residual))
all(rownames(VA.oneSided)==rownames(VA.results))
all(rownames(VA.oneSided)==rownames(VA.coef))
all(rownames(VA.oneSided)==rownames(VA.ebayes))

PRISMO.oneSided<-pt(PRISMO.ebayes$t, df = (PRISMO.ebayes$df.prior+PRISMO.ebayes$df.residual))
all(rownames(PRISMO.oneSided)==rownames(PRISMO.results))
all(rownames(PRISMO.oneSided)==rownames(PRISMO.coef))
all(rownames(PRISMO.oneSided)==rownames(PRISMO.ebayes))

WTC.oneSided<-pt(WTC.ebayes$t, df = (WTC.ebayes$df.prior+WTC.ebayes$df.residual))
all(rownames(WTC.oneSided)==rownames(WTC.results))
all(rownames(WTC.oneSided)==rownames(WTC.coef))
all(rownames(WTC.oneSided)==rownames(WTC.ebayes))

DUKE.oneSided<-pt(DUKE.ebayes$t, df = (DUKE.ebayes$df.prior+DUKE.ebayes$df.residual))
all(rownames(DUKE.oneSided)==rownames(DUKE.results))
all(rownames(DUKE.oneSided)==rownames(DUKE.coef))
all(rownames(DUKE.oneSided)==rownames(DUKE.ebayes))

AS.oneSided<-pt(AS.ebayes$t, df = (AS.ebayes$df.prior+AS.ebayes$df.residual))
all(rownames(AS.oneSided)==rownames(AS.results))
all(rownames(AS.oneSided)==rownames(AS.coef))
all(rownames(AS.oneSided)==rownames(AS.ebayes))

MIR.oneSided<-pt(MIR.ebayes$t, df = (MIR.ebayes$df.prior+MIR.ebayes$df.residual))
all(rownames(MIR.oneSided)==rownames(MIR.results))
all(rownames(MIR.oneSided)==rownames(MIR.coef))
all(rownames(MIR.oneSided)==rownames(MIR.ebayes))

save.image("PGC_EWAS_DataPrep_smoke.Rdata")

# Step 2B: Subset all and non-all sites

studies<-c("DNHS", "GTP", "MRS", "VA", "PRISMO", "WTC", "DUKE", "AS", "MIR")

for(ii in 1:length(studies)){
  coef<-get(paste(studies[ii], ".coef", sep=""))
  coef<-coef[all, ]
  assign(paste(studies[ii], ".coef", sep=""), coef)
  rm(coef)
  results<-get(paste(studies[ii], ".results", sep=""))
  results<-results[all, ]
  assign(paste(studies[ii], ".results", sep=""), results)
  rm(results)
  ebayes<-get(paste(studies[ii], ".ebayes", sep=""))
  ebayes<-ebayes[all, ]
  assign(paste(studies[ii], ".ebayes", sep=""), ebayes)
  rm(ebayes)
  oneSided<-get(paste(studies[ii], ".oneSided", sep=""))
  oneSided<-oneSided[all, ]
  assign(paste(studies[ii], ".oneSided", sep=""), oneSided)
  rm(oneSided)
}

save.image("PGC_EWAS_DataPrep_smoke_intersectSites.Rdata")

rm(list=ls())

load("PGC_EWAS_DataPrep_smoke.Rdata")

studies<-c("DNHS", "GTP", "MRS", "VA", "PRISMO", "WTC", "DUKE", "AS", "MIR")

for(ii in 1:length(studies)){
  coef<-get(paste(studies[ii], ".coef", sep=""))
  coef<-coef[!rownames(coef)%in%all,]
  assign(paste(studies[ii], ".coef", sep=""), coef)
  rm(coef)
  results<-get(paste(studies[ii], ".results", sep=""))
  results<-results[!rownames(results)%in%all,]
  assign(paste(studies[ii], ".results", sep=""), results)
  rm(results)
  ebayes<-get(paste(studies[ii], ".ebayes", sep=""))
  ebayes<-ebayes[!rownames(ebayes)%in%all,]
  assign(paste(studies[ii], ".ebayes", sep=""), ebayes)
  rm(ebayes)
  oneSided<-get(paste(studies[ii], ".oneSided", sep=""))
  oneSided<-oneSided[!rownames(oneSided)%in%all, ]
  assign(paste(studies[ii], ".oneSided", sep=""), oneSided)
  rm(oneSided)
}

save.image("PGC_EWAS_DataPrep_smoke_remainingSites.Rdata")

rm(list=ls())

########################################################################################
# Step 3: Meta-Analysis of Sites in All Studies
########################################################################################

load("PGC_EWAS_DataPrep_smoke_intersectSites.Rdata")

# Step 3A: Check that all rownames line up
all(rownames(DNHS.coef)==rownames(DNHS.results))
all(rownames(DNHS.coef)==rownames(DNHS.ebayes))
all(rownames(DNHS.coef)==names(DNHS.oneSided))

all(rownames(DNHS.coef)==rownames(GTP.coef))
all(rownames(DNHS.coef)==rownames(GTP.results))
all(rownames(DNHS.coef)==rownames(GTP.ebayes))
all(rownames(DNHS.coef)==names(GTP.oneSided))

all(rownames(DNHS.coef)==rownames(MRS.coef))
all(rownames(DNHS.coef)==rownames(MRS.results))
all(rownames(DNHS.coef)==rownames(MRS.ebayes))
all(rownames(DNHS.coef)==names(MRS.oneSided))

all(rownames(DNHS.coef)==rownames(VA.coef))
all(rownames(DNHS.coef)==rownames(VA.results))
all(rownames(DNHS.coef)==rownames(VA.ebayes))
all(rownames(DNHS.coef)==names(VA.oneSided))

all(rownames(DNHS.coef)==rownames(PRISMO.coef))
all(rownames(DNHS.coef)==rownames(PRISMO.results))
all(rownames(DNHS.coef)==rownames(PRISMO.ebayes))
all(rownames(DNHS.coef)==names(PRISMO.oneSided))

all(rownames(DNHS.coef)==rownames(WTC.coef))
all(rownames(DNHS.coef)==rownames(WTC.results))
all(rownames(DNHS.coef)==rownames(WTC.ebayes))
all(rownames(DNHS.coef)==names(WTC.oneSided))

all(rownames(DNHS.coef)==rownames(DUKE.coef))
all(rownames(DNHS.coef)==rownames(DUKE.results))
all(rownames(DNHS.coef)==rownames(DUKE.ebayes))
all(rownames(DNHS.coef)==names(DUKE.oneSided))

all(rownames(DNHS.coef)==rownames(AS.coef))
all(rownames(DNHS.coef)==rownames(AS.results))
all(rownames(DNHS.coef)==rownames(AS.ebayes))
all(rownames(DNHS.coef)==names(AS.oneSided))

all(rownames(DNHS.coef)==rownames(MIR.coef))
all(rownames(DNHS.coef)==rownames(MIR.results))
all(rownames(DNHS.coef)==rownames(MIR.ebayes))
all(rownames(DNHS.coef)==names(MIR.oneSided))

pval.combined.marot<-NULL

for(ii in 1:nrow(DNHS.coef)){
  probe.subjs<-c(DNHS.coef[ii, "N.subjects"], GTP.coef[ii, "N.subjects"], 
                 MRS.coef[ii, "N.subjects"], VA.coef[ii, "N.subjects"],
                 PRISMO.coef[ii, "N.subjects"], WTC.coef[ii, "N.subjects"],
                 DUKE.coef[ii, "N.subjects"], AS.coef[ii, "N.subjects"],
                 MIR.coef[ii, "N.subjects"]) 
  weight<-sqrt(probe.subjs/sum(probe.subjs))
  DNHS.Zi<-qnorm(1-as.numeric(DNHS.oneSided[ii]));DNHS.Zi
  DNHS.Z<-DNHS.Zi*weight[1]; DNHS.Z # multiply by weight 1
  GTP.Zi<-qnorm(1-as.numeric(GTP.oneSided[ii]))
  GTP.Z<-GTP.Zi*weight[2] # multiply by weight 2
  MRS.Zi<-qnorm(1-as.numeric(MRS.oneSided[ii]))
  MRS.Z<-MRS.Zi*weight[3] # multiply by weight 3
  VA.Zi<-qnorm(1-as.numeric(VA.oneSided[ii]))
  VA.Z<-VA.Zi*weight[4] # multiply by weight 4
  PRISMO.Zi<-qnorm(1-as.numeric(PRISMO.oneSided[ii]))
  PRISMO.Z<-PRISMO.Zi*weight[5] # multiply by weight 5
  WTC.Zi<-qnorm(1-as.numeric(WTC.oneSided[ii]))
  WTC.Z<-WTC.Zi*weight[6] # multiply by weight 6
  DUKE.Zi<-qnorm(1-as.numeric(DUKE.oneSided[ii]))
  DUKE.Z<-DUKE.Zi*weight[7] # multiply by weight 7
  AS.Zi<-qnorm(1-as.numeric(AS.oneSided[ii]))
  AS.Z<-AS.Zi*weight[8] # multiply by weight 8
  MIR.Zi<-qnorm(1-as.numeric(MIR.oneSided[ii]))
  MIR.Z<-MIR.Zi*weight[9] # multiply by weight 9
  Sg<-DNHS.Z+GTP.Z+MRS.Z+VA.Z+PRISMO.Z+WTC.Z+DUKE.Z+AS.Z+MIR.Z # add studies together
  p<-2*(1-pnorm(abs(Sg))) # calculated two-tailed p-value
  pval.combined.marot<-append(pval.combined.marot, p)
  if(ii%%10000==0){
    print(ii)
  }
}

pval.combined.marot<-cbind(rownames(DNHS.results), pval.combined.marot)

save(pval.combined.marot, file="PGC_EWAS_smoke_inverseNorm_intersectSites.Rdata")

rm(list=ls())

########################################################################################
# Step 4: Meta-analysis in remaining sties
########################################################################################

load("PGC_EWAS_DataPrep_smoke_remainingSites.Rdata")

sites<-sites[!sites%in%all]

# Step 2B: Add to list of sites the Study IDs
studies<-matrix(FALSE, nrow=length(sites), ncol=9)
colnames(studies)<-c("DNHS", "GTP", "MRS", "VA", "PRISMO", "WTC", "DUKE", "AS", "MIR")
rownames(studies)<-sites

# DNHS
studies[rownames(studies)%in%rownames(DNHS.results), "DNHS"]<-TRUE
miss<-names(which(studies[, "DNHS"]==FALSE)) 
unmat<-rownames(studies)[is.na(match(rownames(studies), rownames(DNHS.results)))]
sum(is.na(match(miss, unmat))) # should be 0
rm(miss, unmat)

# GTP
studies[rownames(studies)%in%rownames(GTP.results), "GTP"]<-TRUE
miss<-names(which(studies[, "GTP"]==FALSE)) 
unmat<-rownames(studies)[is.na(match(rownames(studies), rownames(GTP.results)))]
sum(is.na(match(miss, unmat))) # should be 0
rm(miss, unmat)

# MRS
studies[rownames(studies)%in%rownames(MRS.results), "MRS"]<-TRUE
miss<-names(which(studies[, "MRS"]==FALSE)) 
unmat<-rownames(studies)[is.na(match(rownames(studies), rownames(MRS.results)))]
sum(is.na(match(miss, unmat))) # should be 0
rm(miss, unmat)

# VA
studies[rownames(studies)%in%rownames(VA.results), "VA"]<-TRUE
miss<-names(which(studies[, "VA"]==FALSE)) 
unmat<-rownames(studies)[is.na(match(rownames(studies), rownames(VA.results)))]
sum(is.na(match(miss, unmat))) # should be 0
rm(miss, unmat)

# PRISMO
studies[rownames(studies)%in%rownames(PRISMO.results), "PRISMO"]<-TRUE
miss<-names(which(studies[, "PRISMO"]==FALSE)) 
unmat<-rownames(studies)[is.na(match(rownames(studies), rownames(PRISMO.results)))]
sum(is.na(match(miss, unmat))) # should be 0
rm(miss, unmat)

# WTC
studies[rownames(studies)%in%rownames(WTC.results), "WTC"]<-TRUE
miss<-names(which(studies[, "WTC"]==FALSE)) 
unmat<-rownames(studies)[is.na(match(rownames(studies), rownames(WTC.results)))]
sum(is.na(match(miss, unmat))) # should be 0
rm(miss, unmat)

# DUKE
studies[rownames(studies)%in%rownames(DUKE.results), "DUKE"]<-TRUE
miss<-names(which(studies[, "DUKE"]==FALSE)) 
unmat<-rownames(studies)[is.na(match(rownames(studies), rownames(DUKE.results)))]
sum(is.na(match(miss, unmat))) # should be 0
rm(miss, unmat)

# AS
studies[rownames(studies)%in%rownames(AS.results), "AS"]<-TRUE
miss<-names(which(studies[, "AS"]==FALSE)) 
unmat<-rownames(studies)[is.na(match(rownames(studies), rownames(AS.results)))]
sum(is.na(match(miss, unmat))) # should be 0
rm(miss, unmat)

# MIRECC
studies[rownames(studies)%in%rownames(MIR.results), "MIR"]<-TRUE
miss<-names(which(studies[, "MIR"]==FALSE)) 
unmat<-rownames(studies)[is.na(match(rownames(studies), rownames(MIR.results)))]
sum(is.na(match(miss, unmat))) # should be 0
rm(miss, unmat)

head(studies)
all(apply(studies, 1, sum)>0) # should be true

(dim(studies)[1]*dim(studies)[2])-sum(studies) # this is the number of missing values overall

sum(apply(studies, 1, function(x) 9-sum(x))>0) # this is the number of CpG sites with less than 9 studies
sum(apply(studies, 1, function(x) 8-sum(x))>0) # this is the number of CpG sites with less than 8 studies
sum(apply(studies, 1, function(x) 7-sum(x))>0) # this is the number of CpG sites with less than 7 studies
sum(apply(studies, 1, function(x) 6-sum(x))>0) # this is the number of CpG sites with less than 6 studies
sum(apply(studies, 1, function(x) 5-sum(x))>0) # this is the number of CpG sites with less than 5 studies
sum(apply(studies, 1, function(x) 4-sum(x))>0) # this is the number of CpG sites with less than 4 studies
sum(apply(studies, 1, function(x) 3-sum(x))>0) # this is the number of CpG sites with less than 3 studies
sum(apply(studies, 1, function(x) 2-sum(x))>0) # this is the number of CpG sites with less than 2 studies
sum(apply(studies, 1, function(x) 1-sum(x))>0) # this is the number of CpG sites with less than 1 studies, should be 0

# Step 2C: Loop to calculate p-values
results<-matrix(nrow=nrow(studies), ncol=3)
colnames(results)<-c("CpG", "Studies", "p")

start<-proc.time()[3]
for(ii in 1:nrow(studies)){
  cpg<-rownames(studies)[ii]  # CpG site of analysis
  inc<-colnames(studies)[studies[ii,]] # studies included in analysis
  results[ii, "CpG"]<-cpg
  results[ii, "Studies"]<-paste(inc, collapse=", ")
  
  # Table of parameters for the inverse normal method
  tab<-matrix(nrow=length(inc), ncol=4)
  rownames(tab)<-inc
  colnames(tab)<-c("N", "Weight", "Zi", "Z")
  
  # Getting the number of subjects and calculating weights
  for(jj in 1:length(inc)){
    temp<-get(paste(inc[jj], ".coef", sep=""))
    tab[inc[jj], "N"]<-temp[cpg, "N.subjects"]
    rm(temp)
  }
  tab[, "Weight"]<-sqrt(tab[, "N"]/sum(tab[, "N"]))
  
  # Getting the Z and weighted Z values for each included study
  
  for(ll in 1:length(inc)){
    temp<-get(paste(inc[ll],".oneSided", sep=""))
    tab[inc[ll], "Zi"]<-qnorm(1-as.numeric(temp[cpg]))
    rm(temp)
  }
  tab[, "Z"]<-tab[, "Zi"]*tab[, "Weight"]
  
  # Calculating effect and two-sided p-value
  Sg<-sum(tab[, "Z"]) # add studies together
  results[ii, "p"]<-2*(1-pnorm(abs(Sg))) # calculated two-tailed p-value
  rm(tab)
  
  if(ii%%1000==0){
    print(paste(ii, proc.time()[3]-start, sep=": "))
  }
}

head(results) # check that things look good

save(results, file="PGC_EWAS_smoke_inverseNorm_remainingSites.Rdata")
rm(list=ls())



