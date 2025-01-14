##################Brain Structures####################
setwd("D:/GWAS/GWAS1/Lifestyle Factors")
lo <-list.files(pattern = ".txt$")
lo<-lo[2]
library(openxlsx)
library(TwoSampleMR)
library(data.table)
f<-fread("SmkInit.txt")
head(f)
lo1<-c("Brainstem","left.accumbens","left.amygdala",
       "left.caudate","left.hippocampus","left.insula",
       "left.pallidum","left.putamen","left.thalamus",
       "right.accumbens","right.amygdala","right.caudate",
       "right.hippocampus","right.insula","right.pallidum",
       "right.putamen","right.thalamus")
lo<-c("SmkInit.txt","DrnkWk.txt","CigDay.txt")
lo1<-c("SmkInit","DrnkWk","CigDay")
for (i in 1:3) {
  f<-fread(lo[i])
  format_MR_data(
    dat = f,
    type_name = lo1[i],
    type = "exposure",
    type_binary = F,
    SNP = "SNP",
    chr = "CHR",
    pos = "POS",
    beta = "b",
    se = "se",
    eaf = "Freq",
    effect_allele = "A1",
    other_allele = "A2",
    pval = "p",
    N = 805431 ,
    gene = "gene",
    info = "info",
    ncase = NULL,
    clump_local = FALSE,
    save_path = "D:/GWAS/GWAS1/Mental Health copy"
  )
  
}


###############Biochemical markers###############
library(openxlsx)
library(TwoSampleMR)
library(data.table)

setwd("D:/GWAS/GWAS1/Mental Health")
f<-fread("ADHD.meta")
head(f)
f$beta<-log(f$OR)
f<-subset(f, grepl("rs", SNP))
format_MR_data(
  dat = f,
  type_name = "ADHD",
  type = "exposure",
  type_binary = F,
  SNP = "SNP",
  chr = "CHR",
  pos = "BP",
  beta = "beta",
  se = "se",
  eaf = "Freq",
  effect_allele = "A1",
  other_allele = "A2",
  pval = "p",
  N = 225534,
  gene = "gene",
  info = "INFO",
  ncase = 38691,
  clump_kb = 10000,
  clump_r2 = 0.001,
  pval_threshold = 5e-08,
  clump_local = F,
  save_path =  "D:/GWAS/GWAS1/Mental Health copy"
)

####################CM11###################
setwd("D:/GWAS/GWAS1")
library(data.table)
f<-fread("SA.80")
head(f)

format_MR_data(
  dat = f,
  type_name = "AS",
  type = "exposure",
  type_binary = F,
  SNP = "SNP",
  chr = "CHR.y",
  pos = "BP.y",
  beta = "beta",
  se = "se",
  eaf = "Freq",
  effect_allele = "A1",
  other_allele = "A2",
  pval = "p",
  N = 518233,
  gene = "gene",
  info = "INFO",
  ncase = 26568 ,
  clump_kb = 10000,
  clump_r2 = 0.001,
  pval_threshold = 5e-08,
  clump_local = FALSE,
  save_path = "D:/"
)

exposure_to_outcome_data(
  exp_full_file = "D:/CM1_full_exposure.rds",
  exp_full_data = NULL,
  out_name ="CM1",
  save_path = "D:"
)


###############CM_SA##############
setwd("D:/CM_AS")
library(TwoSampleMR)
f<-readRDS("AS_full_exposure.rds")
dat<-subset(f,f$pval.exposure<1e-06)
dat_clumped<-TwoSampleMR::clump_data(dat,clump_kb = 10000,clump_r2 = 0.001,)
g<-readRDS("CM_full_outcome.rds")
t2d_out <- TwoSampleMR::format_data(
  dat=g,
  type = "outcome",
  snps = dat_clumped$SNP,
  header = TRUE,
  phenotype_col = "outcome",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "pval.outcome",
  chr_col = "chr.outcome",
  samplesize_col = "samplesize.outcome",
  eaf_col = "eaf.outcome",
  pos_col = "pos.outcome",
  ncase_col = "ncase.outcome",
  ncontrol_col = "ncontrol.outcome"
) 
mydata <- harmonise_data(
  exposure_dat=dat_clumped,
  outcome_dat=t2d_out,
  action= 2)
R2<-NULL
R2<-mydata$beta.exposure*mydata$beta.exposure/(mydata$beta.exposure*mydata$beta.exposure+mydata$se.exposure*mydata$se.exposure*mydata$samplesize.exposure)
F_statistics<-R2*mydata$samplesize.exposure/(1-R2)
mydata<-cbind(mydata,R2,F_statistics)
writexl::write_xlsx(x = mydata,path = paste0("AS_to_MC","_mydata.xlsx"))
het=mr_heterogeneity(mydata)
writexl::write_xlsx(x = het,path = paste0("AS_to_MC","_het.xlsx"))
pleio=mr_pleiotropy_test(mydata)
writexl::write_xlsx(x = pleio,path = paste0("AS_to_MC","_pleio.xlsx"))
res_mi=generate_odds_ratios(mr(mydata))
writexl::write_xlsx(x = res_mi,path = paste0("AS_to_MC","_res_mi.xlsx"))

presso<-run_mr_presso(mydata, NbDistribution = 1000, SignifThreshold = 0.05)
presso1<-c(presso[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs,presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue)
presso1<-t(presso1)
presso1<-cbind(presso1,attributes(presso)[["exposure"]],attributes(presso)[["outcome"]])
presso1<-as.data.frame(presso1)
names(presso1)<-c("Global Test`$RSSobs","Test`$Pvalue","exposure","outcome")
writexl::write_xlsx(x = presso1,path = paste0("AS_to_MC","_PRESSO_TEST.xlsx"))
presso4<-presso[[1]][["Main MR results"]]
OR = exp(presso4$`Causal Estimate`)
OR_up = exp(presso4$`Causal Estimate`+1.96*presso4$Sd)
OR_down = exp(presso4$`Causal Estimate`-1.96*presso4$Sd)
press<-NULL
press<-cbind(attributes(presso)[["exposure"]],attributes(presso)[["outcome"]],presso4,OR,OR_down,OR_up)
press<-as.data.frame(press)
writexl::write_xlsx(x = press,path = paste0("AS_to_MC","_PRESSO_MR.xlsx"))

###LDSC#####
ldsc(
  exp_name="AS",
  out_name="CM",
  expfilename = "D:/CM_AS/AS_full_exposure.rds",
  outfilename = "D:/CM_AS孟/CM_full_outcome.rds",
  exp_dat = "AS_full_exposure.rds",
  out_dat = "CM_full_outcome.rds",
  exp_type ="con",
  out_type ="con",
  pop_rate = 0.46,
  remove_HLA = FALSE,
  run_MRlap = T,
  MAF = 0.01,
  pval = 1e-06,
  kb = 10000,
  r2 = 0.001,
  MR_reverse = 5e-06,
  save_logfiles = TRUE,
  save_path ="D:/CM_AS"
)

###placo####
https://github.com/RayDebashree/PLACO?tab=readme-ov-file

require(devtools)#devtools
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.1.1.R?raw=TRUE")

###------------- Code for formally testing pleiotropic association of two phenotypes with SNPs using GWAS summary statistics-----------------
#
message("============================================================")
message("                  PLACO v0.1.1 is loaded")
message("============================================================")
message("If you use this software, please cite:")
message("Ray et al.(2020) A powerful method for pleiotropic analysis")
message("    under composite null hypothesis identifies novel shared")
message("    loci between type 2 diabetes and prostate cancer.")
message("    BioRxiv https://doi.org/10.1101/2020.04.11.037630")
message("------------------------------------------------------------")
message("")

############################################
#---------------- Function for normal product based tail probability calculation 
# (Using modified Bessel function of the 2nd kind with order 0)
.pdfx<-function(x) besselK(x=abs(x),nu=0)/pi
.p.bessel<-function(z, varz, AbsTol=1e-13){
  p1<-2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]/sqrt(varz[1])),Inf, abs.tol=AbsTol)$value)
  p2<-2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]/sqrt(varz[2])),Inf, abs.tol=AbsTol)$value)
  p0<-2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]),Inf, abs.tol=AbsTol)$value)
  pval.compnull<-p1+p2-p0
  return(pval.compnull)
}

#---------------- Function for estimating the variances for PLACO
var.placo<-function(Z.matrix, P.matrix, p.threshold=1e-4){
  # Here Z.matrix is the pxk matrix of Z-scores where p is total no. of variants in the dataset, and k is the no. of traits 
  # Similarly, P.matrix is the corresponding pxk matrix of p-values where p is total no. of variants in the dataset, and k is the no. of traits
  # p.threshold determines which variants are deemed to have no association marginally (default: 1e-4)
  # checks
  k<-ncol(Z.matrix)
  if(k!=2) stop("This method is meant for 2 traits only. Columns correspond to traits.")
  ZP<-cbind(Z.matrix,P.matrix)
  ZP<-na.omit(ZP)

  rows.alt<-which(ZP[,3]<p.threshold & ZP[,4]<p.threshold)
  if(length(rows.alt)>0){
    ZP<-ZP[-rows.alt,]
    if(nrow(ZP)==0) stop(paste("No 'null' variant left at p-value threshold",p.threshold))
    if(nrow(ZP)<30) warning(paste("Too few 'null' variants at p-value threshold",p.threshold))
  }
  varz<-diag(var(ZP[,c(1,2)]))
  return(varz)
}

#---------------- Function for estimating correlation matrix of the Z's
cor.pearson<-function(Z.matrix, P.matrix, p.threshold=1e-4){
  # Here Z.matrix is the pxk matrix of Z-scores where p is total no. of variants in the dataset, and k is the no. of traits 
  # Similarly, P.matrix is the corresponding pxk matrix of p-values where p is total no. of variants in the dataset, and k is the no. of traits
  # p.threshold determines which variants are deemed to have no association marginally (default: 1e-4)
  # checks
  k<-ncol(Z.matrix)
  if(k!=2) stop("This method is meant for 2 traits only.")
  # estimating correlation
    row.exclude<-which( apply(P.matrix, MARGIN = 1, function(x) any(x < p.threshold)) == TRUE )
    if(length(row.exclude)>0) Z.matrix<-Z.matrix[-row.exclude,]
  R<-cor(Z.matrix)
  return(R)
}

############################################
placo<-function(Z, VarZ, AbsTol=.Machine$double.eps^0.8){
        # Z: vector of Z-scores of size k=2 (i.e., collection of Z-scores of a particular SNP for k=2 traits)
        # VarZ: vector of variances of Z-scores (covariance assumed 0; so need to be independent traits)
        # AbsTol: absolute tolerance (accuracy paramater) for numerical integration.
   # checks        
   k<-length(Z)    
   if(k!=2) stop("This method is meant for 2 traits only.")
   if(length(VarZ)!=k) stop("Provide variance estimates for 2 traits as obtained using var.placo() function.")
   
   # test of pleiotropy: PLACO
   pvalue.b=.p.bessel(z=Z, varz=VarZ, AbsTol=AbsTol)
   return(list(T.placo=prod(Z), p.placo=pvalue.b))
}


library(data.table)
f<-fread("CM_AS.txt")
dat<-subset(f,f$p.placo<5e-06)
writexl::write_xlsx(x = dat,path = paste0("MC_to_AS","PLACO.xlsx"))
data<-select(f,c("SNP","chr.exposure","pos.exposure","p.placo","effect_allele.exposure"))
write.table(data,"cm_as_placo.txt",sep = "\t",quote = F,row.names = F)

##########cm1#########
setwd("D:/GWAS/GWAS1/Lifestyle Factors copy")
lo <-list.files(pattern = "outcome.rds$")
CM<-subset(CM_full_exposure,pval.exposure<5e-08)
CM_clump_data<-TwoSampleMR::clump_data(CM,clump_kb = 10000,clump_r2 = 0.001)
CM_clump_data<-clump_data_local(
  data = CM,
  SNP = "SNP",
  pop = "EUR",
  clump_kb = 10000,
  clump_r2 = 0.001,
  p1 = 5e-08,
  p2 = 1,
  core = 10,
  save_path = "D:/"
)
lo1<-strsplit(lo,"f")
for (i in 1:16){
  g<-readRDS(lo[i])
  t2d_out <- TwoSampleMR::format_data(
    dat=g,
    type = "outcome",
    snps = CM_clump_data$SNP,
    header = TRUE,
    phenotype_col = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    chr_col = "chr.outcome",
    samplesize_col = "samplesize.outcome",
    eaf_col = "eaf.outcome",
    pos_col = "pos.outcome",
    ncase_col = "ncase.outcome",
    ncontrol_col = "ncontrol.outcome"
  )
  mydata <- harmonise_data(
    exposure_dat=CM_clump_data,
    outcome_dat=t2d_out,
    action= 2)
  R2<-NULL
  R2<-mydata$beta.exposure*mydata$beta.exposure/(mydata$beta.exposure*mydata$beta.exposure+mydata$se.exposure*mydata$se.exposure*mydata$samplesize.exposure)
  F_statistics<-R2*mydata$samplesize.exposure/(1-R2)
  mydata<-cbind(mydata,R2,F_statistics)
  writexl::write_xlsx(x = mydata,path = paste0("D:/GWAS/GWAS1/Lifestyle Factors copy","/",lo1[[i]][1],"CM_mydata.xlsx"))
  directionality_test <- directionality_test(mydata)
  writexl::write_xlsx(x = directionality_test,path = paste0("D:/GWAS/GWAS1/Lifestyle Factors copy","/",lo1[[i]][1],"CM_directionality_test.xlsx"))
  het=mr_heterogeneity(mydata)
  writexl::write_xlsx(x = het,path = paste0("D:/GWAS/GWAS1/Lifestyle Factors copy","/",lo1[[i]][1],"CM_het.xlsx"))
  pleio=mr_pleiotropy_test(mydata)
  writexl::write_xlsx(x = pleio,path = paste0("D:/GWAS/GWAS1/Lifestyle Factors copy","/",lo1[[i]][1],"CM_pleio.xlsx"))
  res_mi=generate_odds_ratios(mr(mydata))
  writexl::write_xlsx(x = res_mi,path = paste0("D:/GWAS/GWAS1/Lifestyle Factors copy","/",lo1[[i]][1],"CM_res_mi.xlsx"))
}

lo1 <-list.files(pattern = "mydata.xlsx$")
library(openxlsx)
for(j in 1:16){
  f<-read.xlsx(lo1[j])
  presso<-run_mr_presso(f, NbDistribution = 1000, SignifThreshold = 0.05)
  presso1<-c(presso[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs,presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue)
  presso1<-t(presso1)
  presso1<-cbind(presso1,attributes(presso)[["exposure"]],attributes(presso)[["outcome"]])
  presso1<-as.data.frame(presso1)
  names(presso1)<-c("Global Test`$RSSobs","Test`$Pvalue","exposure","outcome")
  writexl::write_xlsx(x = presso1,path = paste0("D:/GWAS/GWAS1/Lifestyle Factors copy/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],".PRESSO_TEST.xlsx"))
  presso4<-presso[[1]][["Main MR results"]]
  OR = exp(presso4$`Causal Estimate`)
  OR_up = exp(presso4$`Causal Estimate`+1.96*presso4$Sd)
  OR_down = exp(presso4$`Causal Estimate`-1.96*presso4$Sd)
  press<-NULL
  press<-cbind(attributes(presso)[["exposure"]],attributes(presso)[["outcome"]],presso4,OR,OR_down,OR_up)
  press<-as.data.frame(press)
  writexl::write_xlsx(x = press,path = paste0("D:/GWAS/GWAS1/Lifestyle Factors copy","/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],".PRESSO_MR.xlsx"))
}
#############cm2##############
setwd("D:/GWAS/GWAS1/Mental Health copy")
lo <-list.files(pattern = "outcome.rds$")
CM<-subset(CM_full_exposure,pval.exposure<5e-08)
CM_clump_data<-TwoSampleMR::clump_data(CM,clump_kb = 10000,clump_r2 = 0.001)
CM_clump_data<-clump_data_local(
  data = CM,
  SNP = "SNP",
  pop = "EUR",
  clump_kb = 10000,
  clump_r2 = 0.001,
  p1 = 5e-08,
  p2 = 1,
  core = 10,
  save_path = "D:/"
)
lo1<-strsplit(lo,"f")
for (i in 1:14){
  g<-readRDS(lo[i])
  t2d_out <- TwoSampleMR::format_data(
    dat=g,
    type = "outcome",
    snps = CM_clump_data$SNP,
    header = TRUE,
    phenotype_col = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    chr_col = "chr.outcome",
    samplesize_col = "samplesize.outcome",
    eaf_col = "eaf.outcome",
    pos_col = "pos.outcome",
    ncase_col = "ncase.outcome",
    ncontrol_col = "ncontrol.outcome"
  )
  mydata <- harmonise_data(
    exposure_dat=CM_clump_data,
    outcome_dat=t2d_out,
    action= 2)
  R2<-NULL
  R2<-mydata$beta.exposure*mydata$beta.exposure/(mydata$beta.exposure*mydata$beta.exposure+mydata$se.exposure*mydata$se.exposure*mydata$samplesize.exposure)
  F_statistics<-R2*mydata$samplesize.exposure/(1-R2)
  mydata<-cbind(mydata,R2,F_statistics)
  writexl::write_xlsx(x = mydata,path = paste0("D:/GWAS/GWAS1/Mental Health copy","/",lo1[[i]][1],"CM_mydata.xlsx"))
  directionality_test <- directionality_test(mydata)
  writexl::write_xlsx(x = directionality_test,path = paste0("D:/GWAS/GWAS1/Mental Health copy","/",lo1[[i]][1],"CM_directionality_test.xlsx"))
  het=mr_heterogeneity(mydata)
  writexl::write_xlsx(x = het,path = paste0("D:/GWAS/GWAS1/Mental Health copy","/",lo1[[i]][1],"CM_het.xlsx"))
  pleio=mr_pleiotropy_test(mydata)
  writexl::write_xlsx(x = pleio,path = paste0("D:/GWAS/GWAS1/Mental Health copy","/",lo1[[i]][1],"CM_pleio.xlsx"))
  res_mi=generate_odds_ratios(mr(mydata))
  writexl::write_xlsx(x = res_mi,path = paste0("D:/GWAS/GWAS1/Mental Health copy","/",lo1[[i]][1],"CM_res_mi.xlsx"))
}

lo1 <-list.files(pattern = "mydata.xlsx$")
library(openxlsx)
for(j in 1:14){
  f<-read.xlsx(lo1[j])
  presso<-run_mr_presso(f, NbDistribution = 1000, SignifThreshold = 0.05)
  presso1<-c(presso[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs,presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue)
  presso1<-t(presso1)
  presso1<-cbind(presso1,attributes(presso)[["exposure"]],attributes(presso)[["outcome"]])
  presso1<-as.data.frame(presso1)
  names(presso1)<-c("Global Test`$RSSobs","Test`$Pvalue","exposure","outcome")
  writexl::write_xlsx(x = presso1,path = paste0("D:/GWAS/GWAS1/Mental Health copy/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],".PRESSO_TEST.xlsx"))
  presso4<-presso[[1]][["Main MR results"]]
  OR = exp(presso4$`Causal Estimate`)
  OR_up = exp(presso4$`Causal Estimate`+1.96*presso4$Sd)
  OR_down = exp(presso4$`Causal Estimate`-1.96*presso4$Sd)
  press<-NULL
  press<-cbind(attributes(presso)[["exposure"]],attributes(presso)[["outcome"]],presso4,OR,OR_down,OR_up)
  press<-as.data.frame(press)
  writexl::write_xlsx(x = press,path = paste0("D:/GWAS/GWAS1/Mental Health copy","/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],".PRESSO_MR.xlsx"))
}

##########cm3##########
setwd("D:/GWAS/GWAS1/Physical Health copy")
lo <-list.files(pattern = "outcome.rds$")
CM<-subset(CM_full_exposure,pval.exposure<5e-08)
CM_clump_data<-TwoSampleMR::clump_data(CM,clump_kb = 10000,clump_r2 = 0.001)
CM_clump_data<-clump_data_local(
  data = CM,
  SNP = "SNP",
  pop = "EUR",
  clump_kb = 10000,
  clump_r2 = 0.001,
  p1 = 5e-08,
  p2 = 1,
  core = 10,
  save_path = "D:/"
)
lo1<-strsplit(lo,"f")
for (i in 5:14){
  g<-readRDS(lo[4])
  t2d_out <- TwoSampleMR::format_data(
    dat=g,
    type = "outcome",
    snps = CM_clump_data$SNP,
    header = TRUE,
    phenotype_col = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    chr_col = "chr.outcome",
    samplesize_col = "samplesize.outcome",
    eaf_col = "eaf.outcome",
    pos_col = "pos.outcome",
    ncase_col = "ncase.outcome",
    ncontrol_col = "ncontrol.outcome"
  )
  mydata <- harmonise_data(
    exposure_dat=CM_clump_data,
    outcome_dat=t2d_out,
    action= 2)
  R2<-NULL
  R2<-mydata$beta.exposure*mydata$beta.exposure/(mydata$beta.exposure*mydata$beta.exposure+mydata$se.exposure*mydata$se.exposure*mydata$samplesize.exposure)
  F_statistics<-R2*mydata$samplesize.exposure/(1-R2)
  mydata<-cbind(mydata,R2,F_statistics)
  writexl::write_xlsx(x = mydata,path = paste0("D:/GWAS/GWAS1/Physical Health copy","/",lo1[[i]][1],"CM_mydata.xlsx"))
  directionality_test <- directionality_test(mydata)
  writexl::write_xlsx(x = directionality_test,path = paste0("D:/GWAS/GWAS1/Physical Health copy","/",lo1[[i]][1],"CM_directionality_test.xlsx"))
  het=mr_heterogeneity(mydata)
  writexl::write_xlsx(x = het,path = paste0("D:/GWAS/GWAS1/Physical Health copy","/",lo1[[i]][1],"CM_het.xlsx"))
  pleio=mr_pleiotropy_test(mydata)
  writexl::write_xlsx(x = pleio,path = paste0("D:/GWAS/GWAS1/Physical Health copy","/",lo1[[i]][1],"CM_pleio.xlsx"))
  res_mi=generate_odds_ratios(mr(mydata))
  writexl::write_xlsx(x = res_mi,path = paste0("D:/GWAS/GWAS1/Physical Health copy","/",lo1[[i]][1],"CM_res_mi.xlsx"))
}

lo1 <-list.files(pattern = "mydata.xlsx$")
library(openxlsx)
for(j in 2:13){
  f<-read.xlsx(lo1[j])
  presso<-run_mr_presso(f, NbDistribution = 1000, SignifThreshold = 0.05)
  presso1<-c(presso[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs,presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue)
  presso1<-t(presso1)
  presso1<-cbind(presso1,attributes(presso)[["exposure"]],attributes(presso)[["outcome"]])
  presso1<-as.data.frame(presso1)
  names(presso1)<-c("Global Test`$RSSobs","Test`$Pvalue","exposure","outcome")
  writexl::write_xlsx(x = presso1,path = paste0("D:/GWAS/GWAS1/Physical Health copy/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],".PRESSO_TEST.xlsx"))
  presso4<-presso[[1]][["Main MR results"]]
  OR = exp(presso4$`Causal Estimate`)
  OR_up = exp(presso4$`Causal Estimate`+1.96*presso4$Sd)
  OR_down = exp(presso4$`Causal Estimate`-1.96*presso4$Sd)
  press<-NULL
  press<-cbind(attributes(presso)[["exposure"]],attributes(presso)[["outcome"]],presso4,OR,OR_down,OR_up)
  press<-as.data.frame(press)
  writexl::write_xlsx(x = press,path = paste0("D:/GWAS/GWAS1/Physical Health copy","/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],".PRESSO_MR.xlsx"))
}

##########cm4##########
setwd("D:/GWAS/GWAS2/Brain Structures copy")
lo <-list.files(pattern = "outcome.rds$")
CM<-subset(CM_full_exposure,pval.exposure<5e-08)
CM_clump_data<-TwoSampleMR::clump_data(CM,clump_kb = 10000,clump_r2 = 0.001)
CM_clump_data<-clump_data_local(
  data = CM,
  SNP = "SNP",
  pop = "EUR",
  clump_kb = 10000,
  clump_r2 = 0.001,
  p1 = 5e-08,
  p2 = 1,
  core = 10,
  save_path = "D:/"
)
lo1<-strsplit(lo,"_")
for (i in 1:17){
  g<-readRDS(lo[i])
  t2d_out <- TwoSampleMR::format_data(
    dat=g,
    type = "outcome",
    snps = CM_clump_data$SNP,
    header = TRUE,
    phenotype_col = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    chr_col = "chr.outcome",
    samplesize_col = "samplesize.outcome",
    eaf_col = "eaf.outcome",
    pos_col = "pos.outcome",
    ncase_col = "ncase.outcome",
    ncontrol_col = "ncontrol.outcome"
  )
  mydata <- harmonise_data(
    exposure_dat=CM_clump_data,
    outcome_dat=t2d_out,
    action= 2)
  R2<-NULL
  R2<-mydata$beta.exposure*mydata$beta.exposure/(mydata$beta.exposure*mydata$beta.exposure+mydata$se.exposure*mydata$se.exposure*mydata$samplesize.exposure)
  F_statistics<-R2*mydata$samplesize.exposure/(1-R2)
  mydata<-cbind(mydata,R2,F_statistics)
  writexl::write_xlsx(x = mydata,path = paste0("D:/GWAS/GWAS2/Brain Structures copy","/",lo1[[i]][1],"CM_mydata.xlsx"))
  directionality_test <- directionality_test(mydata)
  writexl::write_xlsx(x = directionality_test,path = paste0("D:/GWAS/GWAS2/Brain Structures copy","/",lo1[[i]][1],"CM_directionality_test.xlsx"))
  het=mr_heterogeneity(mydata)
  writexl::write_xlsx(x = het,path = paste0("D:/GWAS/GWAS2/Brain Structures copy","/",lo1[[i]][1],"CM_het.xlsx"))
  pleio=mr_pleiotropy_test(mydata)
  writexl::write_xlsx(x = pleio,path = paste0("D:/GWAS/GWAS2/Brain Structures copy","/",lo1[[i]][1],"CM_pleio.xlsx"))
  res_mi=generate_odds_ratios(mr(mydata))
  writexl::write_xlsx(x = res_mi,path = paste0("D:/GWAS/GWAS2/Brain Structures copy","/",lo1[[i]][1],"CM_res_mi.xlsx"))
}

lo1 <-list.files(pattern = "mydata.xlsx$")
library(openxlsx)
for(j in 12:17){
  f<-read.xlsx(lo1[j])
  presso<-run_mr_presso(f, NbDistribution = 1000, SignifThreshold = 0.05)
  presso1<-c(presso[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs,presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue)
  presso1<-t(presso1)
  presso1<-cbind(presso1,attributes(presso)[["exposure"]],attributes(presso)[["outcome"]])
  presso1<-as.data.frame(presso1)
  names(presso1)<-c("Global Test`$RSSobs","Test`$Pvalue","exposure","outcome")
  writexl::write_xlsx(x = presso1,path = paste0("D:/GWAS/GWAS2/Brain Structures copy/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],".PRESSO_TEST.xlsx"))
  presso4<-presso[[1]][["Main MR results"]]
  OR = exp(presso4$`Causal Estimate`)
  OR_up = exp(presso4$`Causal Estimate`+1.96*presso4$Sd)
  OR_down = exp(presso4$`Causal Estimate`-1.96*presso4$Sd)
  press<-NULL
  press<-cbind(attributes(presso)[["exposure"]],attributes(presso)[["outcome"]],presso4,OR,OR_down,OR_up)
  press<-as.data.frame(press)
  writexl::write_xlsx(x = press,path = paste0("D:/GWAS/GWAS2/Brain Structures copy","/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],".PRESSO_MR.xlsx"))
}


setwd("D:/GWAS/GWAS1/Physical Health copy/")
a1<-"data.txt"
lo1 <-list.files(pattern = "directionality_test.xlsx$")
for (i in 1:13) {
  a<-readxl::read_excel(lo1[i])
  a1<-rbind(a1,a)
}
writexl::write_xlsx(x = a1,"MR_Brain Structures_directionality_test.xlsx")

a2<-"data.txt"
lo2 <-list.files(pattern = "res_mi.xlsx$")
for (i in 1:13) {
  a<-readxl::read_excel(lo2[i])
  a2<-rbind(a2,a)
}
writexl::write_xlsx(x = a2,"MR_Brain Structures.xlsx")

a3<-"data.txt"
lo3 <-list.files(pattern = "pleio.xlsx$")
for (i in 1:13) {
  a<-readxl::read_excel(lo3[i])
  a3<-rbind(a3,a)
}
writexl::write_xlsx(x = a3,"MR_Brain Structures_pleio.xlsx")

a4<-"data.txt"
lo4 <-list.files(pattern = "het.xlsx$")
for (i in 1:13) {
  a<-readxl::read_excel(lo4[i])
  a4<-rbind(a4,a)
}
writexl::write_xlsx(x = a4,"MR_Brain Structures_het.xlsx")


a5<-"data.txt"
lo5 <-list.files(pattern = "PRESSO_MR.xlsx$")
for (i in 1:13) {
  a<-readxl::read_excel(lo5[i])
  a5<-rbind(a5,a)
}
writexl::write_xlsx(x = a5,"MR_Brain Structures_PRESSO_MR.xlsx")

a6<-"data.txt"
lo6 <-list.files(pattern = "PRESSO_TEST.xlsx$")
for (i in 1:13) {
  a<-readxl::read_excel(lo6[i])
  a6<-rbind(a6,a)
}
writexl::write_xlsx(x = a6,"MR_Brain Structures_PRESSO_TEST.xlsx")

##########AS1########
setwd("D:/GWAS/GWAS1/Lifestyle Factors copy")
lo <-list.files(pattern = "exposure.rds$")
CM_clump_data<-clump_data_local(
  data = CM,
  SNP = "SNP",
  pop = "EUR",
  clump_kb = 10000,
  clump_r2 = 0.001,
  p1 = 5e-08,
  p2 = 1,
  core = 10,
  save_path = "D:/"
)
lo1<-strsplit(lo,"f")
for (i in 1:16){
  g<-readRDS(lo[i])
  g<-subset(g,pval.exposure<5e-08)
  g<-TwoSampleMR::clump_data(g,clump_kb = 10000,clump_r2 = 0.001)
  t2d_out <- TwoSampleMR::format_data(
    dat=AS_full_outcome,
    type = "outcome",
    snps = g$SNP,
    header = TRUE,
    phenotype_col = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    chr_col = "chr.outcome",
    samplesize_col = "samplesize.outcome",
    eaf_col = "eaf.outcome",
    pos_col = "pos.outcome",
    ncase_col = "ncase.outcome",
    ncontrol_col = "ncontrol.outcome"
  )
  mydata <- harmonise_data(
    exposure_dat=g,
    outcome_dat=t2d_out,
    action= 2)
  R2<-NULL
  R2<-mydata$beta.exposure*mydata$beta.exposure/(mydata$beta.exposure*mydata$beta.exposure+mydata$se.exposure*mydata$se.exposure*mydata$samplesize.exposure)
  F_statistics<-R2*mydata$samplesize.exposure/(1-R2)
  mydata<-cbind(mydata,R2,F_statistics)
  writexl::write_xlsx(x = mydata,path = paste0("D:/GWAS/GWAS1/Lifestyle Factors copy","/",lo1[[i]][1],"AS_mydata.xlsx"))
  directionality_test <- directionality_test(mydata)
  writexl::write_xlsx(x = directionality_test,path = paste0("D:/GWAS/GWAS1/Lifestyle Factors copy","/",lo1[[i]][1],"AS_directionality_test.xlsx"))
  het=mr_heterogeneity(mydata)
  writexl::write_xlsx(x = het,path = paste0("D:/GWAS/GWAS1/Lifestyle Factors copy","/",lo1[[i]][1],"AS_het.xlsx"))
  pleio=mr_pleiotropy_test(mydata)
  writexl::write_xlsx(x = pleio,path = paste0("D:/GWAS/GWAS1/Lifestyle Factors copy","/",lo1[[i]][1],"AS_pleio.xlsx"))
  res_mi=generate_odds_ratios(mr(mydata))
  writexl::write_xlsx(x = res_mi,path = paste0("D:/GWAS/GWAS1/Lifestyle Factors copy","/",lo1[[i]][1],"AS_res_mi.xlsx"))
}

lo1 <-list.files(pattern = "mydata.xlsx$")
library(openxlsx)
for(j in 13:16){
  f<-read.xlsx(lo1[j])
  presso<-run_mr_presso(f, NbDistribution = 1000, SignifThreshold = 0.05)
  presso1<-c(presso[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs,presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue)
  presso1<-t(presso1)
  presso1<-cbind(presso1,attributes(presso)[["exposure"]],attributes(presso)[["outcome"]])
  presso1<-as.data.frame(presso1)
  names(presso1)<-c("Global Test`$RSSobs","Test`$Pvalue","exposure","outcome")
  writexl::write_xlsx(x = presso1,path = paste0("D:/GWAS/GWAS1/Lifestyle Factors copy/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],"AS.PRESSO_TEST.xlsx"))
  presso4<-presso[[1]][["Main MR results"]]
  OR = exp(presso4$`Causal Estimate`)
  OR_up = exp(presso4$`Causal Estimate`+1.96*presso4$Sd)
  OR_down = exp(presso4$`Causal Estimate`-1.96*presso4$Sd)
  press<-NULL
  press<-cbind(attributes(presso)[["exposure"]],attributes(presso)[["outcome"]],presso4,OR,OR_down,OR_up)
  press<-as.data.frame(press)
  writexl::write_xlsx(x = press,path = paste0("D:/GWAS/GWAS1/Lifestyle Factors copy","/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],"AS.PRESSO_MR.xlsx"))
}


##########AS2##########
setwd("D:/GWAS/GWAS1/Mental Health copy")
lo <-list.files(pattern = "exposure.rds$")
CM_clump_data<-clump_data_local(
  data = CM,
  SNP = "SNP",
  pop = "EUR",
  clump_kb = 10000,
  clump_r2 = 0.001,
  p1 = 5e-08,
  p2 = 1,
  core = 10,
  save_path = "D:/"
)
lo1<-strsplit(lo,"f")
for (i in 9:9){
  g<-readRDS(lo[i])
  g<-subset(g,pval.exposure<1e-05)
  g<-TwoSampleMR::clump_data(g,clump_kb = 10000,clump_r2 = 0.001)
  t2d_out <- TwoSampleMR::format_data(
    dat=AS_full_outcome,
    type = "outcome",
    snps = g$SNP,
    header = TRUE,
    phenotype_col = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    chr_col = "chr.outcome",
    samplesize_col = "samplesize.outcome",
    eaf_col = "eaf.outcome",
    pos_col = "pos.outcome",
    ncase_col = "ncase.outcome",
    ncontrol_col = "ncontrol.outcome"
  )
  mydata <- harmonise_data(
    exposure_dat=g,
    outcome_dat=t2d_out,
    action= 2)
  R2<-NULL
  R2<-mydata$beta.exposure*mydata$beta.exposure/(mydata$beta.exposure*mydata$beta.exposure+mydata$se.exposure*mydata$se.exposure*mydata$samplesize.exposure)
  F_statistics<-R2*mydata$samplesize.exposure/(1-R2)
  mydata<-cbind(mydata,R2,F_statistics)
  writexl::write_xlsx(x = mydata,path = paste0("D:/GWAS/GWAS1/Mental Health copy","/",lo1[[i]][1],"AS_mydata.xlsx"))
  directionality_test <- directionality_test(mydata)
  writexl::write_xlsx(x = directionality_test,path = paste0("D:/GWAS/GWAS1/Mental Health copy","/",lo1[[i]][1],"AS_directionality_test.xlsx"))
  het=mr_heterogeneity(mydata)
  writexl::write_xlsx(x = het,path = paste0("D:/GWAS/GWAS1/Mental Health copy","/",lo1[[i]][1],"AS_het.xlsx"))
  pleio=mr_pleiotropy_test(mydata)
  writexl::write_xlsx(x = pleio,path = paste0("D:/GWAS/GWAS1/Mental Health copy","/",lo1[[i]][1],"AS_pleio.xlsx"))
  res_mi=generate_odds_ratios(mr(mydata))
  writexl::write_xlsx(x = res_mi,path = paste0("D:/GWAS/GWAS1/Mental Health copy","/",lo1[[i]][1],"AS_res_mi.xlsx"))
}

lo1 <-list.files(pattern = "mydata.xlsx$")
library(openxlsx)
for(j in 9:16){
  f<-read.xlsx(lo1[j])
  presso<-run_mr_presso(f, NbDistribution = 1000, SignifThreshold = 0.05)
  presso1<-c(presso[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs,presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue)
  presso1<-t(presso1)
  presso1<-cbind(presso1,attributes(presso)[["exposure"]],attributes(presso)[["outcome"]])
  presso1<-as.data.frame(presso1)
  names(presso1)<-c("Global Test`$RSSobs","Test`$Pvalue","exposure","outcome")
  writexl::write_xlsx(x = presso1,path = paste0("D:/GWAS/GWAS1/Mental Health copy/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],"AS.PRESSO_TEST.xlsx"))
  presso4<-presso[[1]][["Main MR results"]]
  OR = exp(presso4$`Causal Estimate`)
  OR_up = exp(presso4$`Causal Estimate`+1.96*presso4$Sd)
  OR_down = exp(presso4$`Causal Estimate`-1.96*presso4$Sd)
  press<-NULL
  press<-cbind(attributes(presso)[["exposure"]],attributes(presso)[["outcome"]],presso4,OR,OR_down,OR_up)
  press<-as.data.frame(press)
  writexl::write_xlsx(x = press,path = paste0("D:/GWAS/GWAS1/Mental Health copy","/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],"AS.PRESSO_MR.xlsx"))
}

##########AS3########
setwd("D:/GWAS/GWAS1/Physical Health copy")
lo <-list.files(pattern = "exposure.rds$")
CM_clump_data<-clump_data_local(
  data = CM,
  SNP = "SNP",
  pop = "EUR",
  clump_kb = 10000,
  clump_r2 = 0.001,
  p1 = 5e-08,
  p2 = 1,
  core = 10,
  save_path = "D:/"
)
lo1<-strsplit(lo,"f")
for (i in 4:14){
  g<-readRDS(lo[i])
  g<-subset(g,pval.exposure<5e-08)
  g<-TwoSampleMR::clump_data(g,clump_kb = 10000,clump_r2 = 0.001)
  t2d_out <- TwoSampleMR::format_data(
    dat=AS_full_outcome,
    type = "outcome",
    snps = g$SNP,
    header = TRUE,
    phenotype_col = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    chr_col = "chr.outcome",
    samplesize_col = "samplesize.outcome",
    eaf_col = "eaf.outcome",
    pos_col = "pos.outcome",
    ncase_col = "ncase.outcome",
    ncontrol_col = "ncontrol.outcome"
  )
  mydata <- harmonise_data(
    exposure_dat=g,
    outcome_dat=t2d_out,
    action= 2)
  R2<-NULL
  R2<-mydata$beta.exposure*mydata$beta.exposure/(mydata$beta.exposure*mydata$beta.exposure+mydata$se.exposure*mydata$se.exposure*mydata$samplesize.exposure)
  F_statistics<-R2*mydata$samplesize.exposure/(1-R2)
  mydata<-cbind(mydata,R2,F_statistics)
  writexl::write_xlsx(x = mydata,path = paste0("D:/GWAS/GWAS1/Physical Health copy","/",lo1[[i]][1],"AS_mydata.xlsx"))
  directionality_test <- directionality_test(mydata)
  writexl::write_xlsx(x = directionality_test,path = paste0("D:/GWAS/GWAS1/Physical Health copy","/",lo1[[i]][1],"AS_directionality_test.xlsx"))
  het=mr_heterogeneity(mydata)
  writexl::write_xlsx(x = het,path = paste0("D:/GWAS/GWAS1/Physical Health copy","/",lo1[[i]][1],"AS_het.xlsx"))
  pleio=mr_pleiotropy_test(mydata)
  writexl::write_xlsx(x = pleio,path = paste0("D:/GWAS/GWAS1/Physical Health copy","/",lo1[[i]][1],"AS_pleio.xlsx"))
  res_mi=generate_odds_ratios(mr(mydata))
  writexl::write_xlsx(x = res_mi,path = paste0("D:/GWAS/GWAS1/Physical Health copy","/",lo1[[i]][1],"AS_res_mi.xlsx"))
}

lo1 <-list.files(pattern = "AS_mydata.xlsx$")
library(openxlsx)
for(j in 1:14){
  f<-read.xlsx(lo1[j])
  presso<-run_mr_presso(f, NbDistribution = 1000, SignifThreshold = 0.05)
  presso1<-c(presso[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs,presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue)
  presso1<-t(presso1)
  presso1<-cbind(presso1,attributes(presso)[["exposure"]],attributes(presso)[["outcome"]])
  presso1<-as.data.frame(presso1)
  names(presso1)<-c("Global Test`$RSSobs","Test`$Pvalue","exposure","outcome")
  writexl::write_xlsx(x = presso1,path = paste0("D:/GWAS/GWAS1/Physical Health copy/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],"AS.PRESSO_TEST.xlsx"))
  presso4<-presso[[1]][["Main MR results"]]
  OR = exp(presso4$`Causal Estimate`)
  OR_up = exp(presso4$`Causal Estimate`+1.96*presso4$Sd)
  OR_down = exp(presso4$`Causal Estimate`-1.96*presso4$Sd)
  press<-NULL
  press<-cbind(attributes(presso)[["exposure"]],attributes(presso)[["outcome"]],presso4,OR,OR_down,OR_up)
  press<-as.data.frame(press)
  writexl::write_xlsx(x = press,path = paste0("D:/GWAS/GWAS1/Physical Health copy","/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],"AS.PRESSO_MR.xlsx"))
}


##########AS4###########
setwd("D:/GWAS/GWAS2/Brain Structures copy")
lo <-list.files(pattern = "exposure.rds$")
CM_clump_data<-clump_data_local(
  data = CM,
  SNP = "SNP",
  pop = "EUR",
  clump_kb = 10000,
  clump_r2 = 0.001,
  p1 = 5e-08,
  p2 = 1,
  core = 10,
  save_path = "D:/"
)
lo1<-strsplit(lo,"_")
for (i in 1:17){
  g<-readRDS(lo[i])
  g<-subset(g,pval.exposure<5e-08)
  g<-TwoSampleMR::clump_data(g,clump_kb = 10000,clump_r2 = 0.001)
  t2d_out <- TwoSampleMR::format_data(
    dat=AS_full_outcome,
    type = "outcome",
    snps = g$SNP,
    header = TRUE,
    phenotype_col = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    chr_col = "chr.outcome",
    samplesize_col = "samplesize.outcome",
    eaf_col = "eaf.outcome",
    pos_col = "pos.outcome",
    ncase_col = "ncase.outcome",
    ncontrol_col = "ncontrol.outcome"
  )
  mydata <- harmonise_data(
    exposure_dat=g,
    outcome_dat=t2d_out,
    action= 2)
  R2<-NULL
  R2<-mydata$beta.exposure*mydata$beta.exposure/(mydata$beta.exposure*mydata$beta.exposure+mydata$se.exposure*mydata$se.exposure*mydata$samplesize.exposure)
  F_statistics<-R2*mydata$samplesize.exposure/(1-R2)
  mydata<-cbind(mydata,R2,F_statistics)
  writexl::write_xlsx(x = mydata,path = paste0("D:/GWAS/GWAS2/Brain Structures copy","/",lo1[[i]][1],"AS_mydata.xlsx"))
  directionality_test <- directionality_test(mydata)
  writexl::write_xlsx(x = directionality_test,path = paste0("D:/GWAS/GWAS2/Brain Structures copy","/",lo1[[i]][1],"AS_directionality_test.xlsx"))
  het=mr_heterogeneity(mydata)
  writexl::write_xlsx(x = het,path = paste0("D:/GWAS/GWAS2/Brain Structures copy","/",lo1[[i]][1],"AS_het.xlsx"))
  pleio=mr_pleiotropy_test(mydata)
  writexl::write_xlsx(x = pleio,path = paste0("D:/GWAS/GWAS2/Brain Structures copy","/",lo1[[i]][1],"AS_pleio.xlsx"))
  res_mi=generate_odds_ratios(mr(mydata))
  writexl::write_xlsx(x = res_mi,path = paste0("D:/GWAS/GWAS2/Brain Structures copy","/",lo1[[i]][1],"AS_res_mi.xlsx"))
}

lo1 <-list.files(pattern = "AS_mydata.xlsx$")
library(openxlsx)
for(j in 12:17){
  f<-read.xlsx(lo1[j])
  presso<-run_mr_presso(f, NbDistribution = 1000, SignifThreshold = 0.05)
  presso1<-c(presso[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs,presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue)
  presso1<-t(presso1)
  presso1<-cbind(presso1,attributes(presso)[["exposure"]],attributes(presso)[["outcome"]])
  presso1<-as.data.frame(presso1)
  names(presso1)<-c("Global Test`$RSSobs","Test`$Pvalue","exposure","outcome")
  writexl::write_xlsx(x = presso1,path = paste0("D:/GWAS/GWAS2/Brain Structures copy/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],"AS.PRESSO_TEST.xlsx"))
  presso4<-presso[[1]][["Main MR results"]]
  OR = exp(presso4$`Causal Estimate`)
  OR_up = exp(presso4$`Causal Estimate`+1.96*presso4$Sd)
  OR_down = exp(presso4$`Causal Estimate`-1.96*presso4$Sd)
  press<-NULL
  press<-cbind(attributes(presso)[["exposure"]],attributes(presso)[["outcome"]],presso4,OR,OR_down,OR_up)
  press<-as.data.frame(press)
  writexl::write_xlsx(x = press,path = paste0("D:/GWAS/GWAS2/Brain Structures copy","/",attributes(presso)[["exposure"]],".",attributes(presso)[["outcome"]],"AS.PRESSO_MR.xlsx"))
}




#########cbind#########
setwd("D:/GWAS/GWAS2/Brain Structures copy/)
a1<-"data.txt"
lo1 <-list.files(pattern = "directionality_test.xlsx$")
for (i in 1:17) {
  a<-readxl::read_excel(lo1[i])
  a1<-rbind(a1,a)
}
writexl::write_xlsx(x = a1,"MR_Brain Structures_directionality_test.xlsx")

a2<-"data.txt"
lo2 <-list.files(pattern = "res_mi.xlsx$")
for (i in 1:17) {
  a<-readxl::read_excel(lo2[i])
  a2<-rbind(a2,a)
}
writexl::write_xlsx(x = a2,"MR_Brain Structures.xlsx")

a3<-"data.txt"
lo3 <-list.files(pattern = "pleio.xlsx$")
for (i in 1:17) {
  a<-readxl::read_excel(lo3[i])
  a3<-rbind(a3,a)
}
writexl::write_xlsx(x = a3,"MR_Brain Structures_pleio.xlsx")

a4<-"data.txt"
lo4 <-list.files(pattern = "het.xlsx$")
for (i in 1:17) {
  a<-readxl::read_excel(lo4[i])
  a4<-rbind(a4,a)
}
writexl::write_xlsx(x = a4,"MR_Brain Structures_het.xlsx")

a5<-"data.txt"
lo5 <-list.files(pattern = "PRESSO_MR.xlsx$")
for (i in 1:12) {
  a<-readxl::read_excel(lo5[i])
  a5<-rbind(a5,a)
}
writexl::write_xlsx(x = a5,"MR_Brain Structures_PRESSO_MR.xlsx")

a6<-"data.txt"
lo6 <-list.files(pattern = "PRESSO_TEST.xlsx$")
for (i in 1:12) {
  a<-readxl::read_excel(lo6[i])
  a6<-rbind(a6,a)
}
writexl::write_xlsx(x = a6,"MR_Brain Structures_PRESSO_TEST.xlsx")


#########中介分析#########
library(openxlsx)
a<-read.xlsx("D:/GWAS/GWAS1/Physical Health copy/MR_Brain Structures_.xlsx")
b<-read.xlsx("D:/GWAS/GWAS1/Physical Health copy/MR_Physical Health.xlsx")

a<-subset(a,method=="Inverse variance weighted"|
            method=="Wald ratio") 
b<-subset(b,method=="Inverse variance weighted"|
            method=="Wald ratio") 

mediated<-as.numeric(a$b)*as.numeric(b$b)

S=sqrt(as.numeric(a$b)^2*as.numeric(b$se)^2+as.numeric(b$b)^2*as.numeric(a$se)^2)

Z=mediated/S
P=2*pnorm(q=abs(Z), lower.tail=FALSE) #Z=1.96  P=0.05

mean_indirect <- as.numeric(a$b)*as.numeric(b$b)

c<-as.numeric(a$b)
d<-as.numeric(b$b)
e<-as.numeric(a$se)
f<-as.numeric(b$se)

indirect_and_prop <- function(exposure_mediator,
                              exposure_mediator_se,
                              mediator_outcome,
                              mediator_outcome_se,
                              exposure_outcome,
                              exposure_outcome_se) {
  
  mean_indirect <- exposure_mediator*mediator_outcome
  
  m1 <- eval(stats::D(expression(exposure_mediator*mediator_outcome), "exposure_mediator"))
  
  m2 <- eval(stats::D(expression(exposure_mediator*mediator_outcome), "mediator_outcome"))
  
  indirect_se <- sqrt((m1^2)*exposure_mediator_se^2 + (m2^2)*mediator_outcome_se^2)
  
  mean_prop_mediated <- mean_indirect / exposure_outcome
  
  m3 <- eval(stats::D(expression(mean_indirect / exposure_outcome), "mean_indirect"))
  
  m4 <- eval(stats::D(expression(mean_indirect / exposure_outcome), "exposure_outcome"))
  
  prop_mediated_se <- sqrt((m3^2)*indirect_se^2 + (m4^2)*exposure_outcome_se^2)
  
  indirect_Z <- as.numeric(mean_indirect/indirect_se)
  
  indirect_pval <- as.numeric(2*stats::pnorm(-abs(indirect_Z)))
  
  res <- list(indirect = mean_indirect,
              indirect_se = indirect_se,
              prop_mediated = mean_prop_mediated,
              prop_mediated_se = prop_mediated_se,
              lci_prop_mediated = mean_prop_mediated - 1.96*prop_mediated_se,
              uci_prop_mediated = mean_prop_mediated + 1.96*prop_mediated_se,
              indirect_Z = indirect_Z,
              indirect_pval = indirect_pval)
  
  return (dplyr::bind_rows(res))
}

dd<-indirect_and_prop(c,
                      e,
                      d,
                      f,
                      1.223,
                      0.19)
cc<-cbind(a$outcome,dd)
writexl::write_xlsx(x = cc,"D:/GWAS/GWAS1/Physical Health copy/m_Physical Health.xlsx")


format_MR_data <- function(dat=NULL,
                           filepath,
                           type_name,
                           type,
                           type_binary=TRUE,
                           SNP="SNP",
                           chr="chr",
                           pos="pos",
                           beta="beta",
                           se="se",
                           eaf="eaf",
                           effect_allele="effect_allele",
                           other_allele="other_allele",
                           pval="pval",
                           N="N",
                           gene="gene",
                           info="info",
                           ncase,
                           filter_maf=0,
                           clump_kb=10000,
                           clump_r2=0.001,
                           pval_threshold=5e-8,
                           clump_local=TRUE,
                           save_path){
  
  if (.checkaccess()){
    dat <-  .format_MR_data(dat,
                            filepath,
                            type_name,
                            type,
                            type_binary,
                            SNP,
                            chr,
                            pos,
                            beta,
                            se,
                            eaf,
                            effect_allele,
                            other_allele,
                            pval,
                            N,
                            gene,
                            info,
                            ncase,
                            filter_maf,
                            clump_kb,
                            clump_r2,
                            pval_threshold,
                            clump_local,
                            save_path)
  }
  return(dat)
}


.format_MR_data <- function(dat,
                            filepath,
                            type_name,
                            type,
                            type_binary,
                            SNP,
                            chr,
                            pos,
                            beta,
                            se,
                            eaf,
                            effect_allele,
                            other_allele,
                            pval,
                            N,
                            gene,
                            info,
                            ncase,
                            filter_maf,
                            clump_kb,
                            clump_r2,
                            pval_threshold,
                            clump_local,
                            save_path){

  
  if(!is.null(dat)){
    
    dat <- dat%>%
      data.frame()
    
  }else if(!is.null(filepath)){
    
    dat <-data.table::fread(filepath)%>%
      data.frame()
    
  }
  
  if(!is.null(chr) && all(c(chr,pos) %in% colnames(dat))){
    
    dat <- dat%>%
      dplyr::mutate(chr=as.integer(chr),
                    pos=as.integer(pos))
  }
  
  if(!is.null(effect_allele) && !is.null(other_allele) && all(c(effect_allele,other_allele) %in% colnames(dat))){
    
    dat <- dat%>%
      dplyr::mutate_at(c(effect_allele,other_allele),toupper)
    
  }
  
  if(type=="exposure"){
    
    if(class(N)=="character"){
      
      dat <- dat%>%
        TwoSampleMR::format_data(.,
                                 type=type,
                                 snp_col = SNP,
                                 chr_col = chr,
                                 pos_col = pos,
                                 beta_col = beta,
                                 se_col =se ,
                                 eaf_col = eaf,
                                 effect_allele_col = effect_allele,
                                 other_allele_col = other_allele,
                                 pval_col = pval,
                                 samplesize_col=N,
                                 gene_col = gene,
                                 info_col = info)
    }else{
      
      dat <- dat%>%
        TwoSampleMR::format_data(.,
                                 type=type,
                                 snp_col = SNP,
                                 chr_col = chr,
                                 pos_col = pos,
                                 beta_col = beta,
                                 se_col =se ,
                                 eaf_col = eaf,
                                 effect_allele_col = effect_allele,
                                 other_allele_col = other_allele,
                                 pval_col = pval,
                                 gene_col = gene,
                                 info_col = info)
    }
    
    dat <- dat%>%
      dplyr::mutate(exposure=type_name)%>%
      dplyr::mutate(id.exposure=type_name)
    
    
    if(any(is.na(dat$eaf.exposure))){
      
      
      dat$eaf.exposure <- as.numeric(dat$eaf.exposure)
      
      dat[is.na(dat$eaf.exposure),]$eaf.exposure <- 0.5
      
    }
    
    if((is.na(dat$se.exposure[1])) && ("beta.exposure" %in% colnames(dat))  && ("pval.exposure" %in% colnames(dat))){
      
      
      dat$se.exposure <- sqrt(((dat$beta.exposure)^2)/qchisq(dat$pval.exposure, 1, lower.tail=FALSE))
      
    }
    
    if(!("samplesize.exposure" %in% colnames(dat))){
      
      if(class(N)=="numeric"){
        
        dat$samplesize.exposure <- as.numeric(N)
        
      }
      
    }
    
    if(type_binary==TRUE){
      
      dat$ncase.exposure <- ncase
      
      dat$ncontrol.exposure <- dat$samplesize.exposure - dat$ncase.exposure
      
    }
    
    if(!is.null(eaf)){
      
      if("eaf.exposure" %in% colnames(dat) && all(!is.na(dat$eaf.exposure))){
        
        dat <- dat%>%
          dplyr::filter(eaf.exposure > filter_maf)%>%
          dplyr::filter(eaf.exposure < (1-filter_maf))
        
      }
      
    }
    
    saveRDS(dat,file = paste0(save_path,"/",type_name,"_full_exposure.rds"))
    
    dat <- dat%>%
      dplyr::filter(pval.exposure < pval_threshold)%>%
      clump_data(.,clump_local = clump_local,
                                   save_path = save_path,
                                   clump_kb = clump_kb,
                                   clump_r2 = clump_r2,
                                   clump_p1 = pval_threshold)%>%
      dplyr::arrange(.,SNP)
    
    utils::write.table(x = dat,file = paste0(save_path,"/",type_name,".txt"),row.names =FALSE)
    
    
  }
  
  if(type=="outcome"){
    
    if(class(N)=="character"){
      
      dat <- dat%>%
        TwoSampleMR::format_data(.,
                                 type=type,
                                 snp_col = SNP,
                                 chr_col = chr,
                                 pos_col = pos,
                                 beta_col = beta,
                                 se_col =se ,
                                 eaf_col = eaf,
                                 effect_allele_col = effect_allele,
                                 other_allele_col = other_allele,
                                 pval_col = pval,
                                 samplesize_col=N)%>%
        dplyr::mutate(outcome=type_name)%>%
        dplyr::mutate(id.outcome=type_name)
      
    }else{
      
      dat <- dat%>%
        TwoSampleMR::format_data(.,
                                 type=type,
                                 snp_col = SNP,
                                 chr_col = chr,
                                 pos_col = pos,
                                 beta_col = beta,
                                 se_col =se ,
                                 eaf_col = eaf,
                                 effect_allele_col = effect_allele,
                                 other_allele_col = other_allele,
                                 pval_col = pval)%>%
        dplyr::mutate(outcome=type_name)%>%
        dplyr::mutate(id.outcome=type_name)
      
    }
    
    if(any(is.na(dat$eaf.outcome))){
      
      dat$eaf.outcome <- as.numeric(dat$eaf.outcome)
      
      dat[is.na(dat$eaf.outcome),]$eaf.outcome <- 0.5
      
    }
    
    if((is.na(dat$se.outcome[1])) && ("beta.outcome" %in% colnames(dat)) && ("pval.outcome" %in% colnames(dat))){
      
      
      dat$se.outcome <- sqrt(((dat$beta.outcome)^2)/qchisq(dat$pval.outcome, 1, lower.tail=FALSE))
      
    }
    
    if(!"samplesize.outcome" %in% colnames(dat)){
      
      if(class(N)=="numeric"){
        
        dat$samplesize.outcome <- as.numeric(N)
        
      }
      
      
    }
    
    if(type_binary==TRUE){
      
      dat$ncase.outcome <- ncase
      
      dat$ncontrol.outcome <- dat$samplesize.outcome - dat$ncase.outcome
      
    }
    
    if(!is.null(eaf)){
      
      if("eaf.outcome" %in% colnames(dat) && all(!is.na(dat$eaf.outcome))){
        
        dat <- dat%>%
          dplyr::filter(eaf.outcome > filter_maf)%>%
          dplyr::filter(eaf.outcome < (1-filter_maf))
        
      }
      
    }
    
    saveRDS(dat,file = paste0(save_path,"/",type_name,"_full_outcome.rds"))
    
  }
  
  return(dat)
  
  
}

