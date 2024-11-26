# Association analysis in DMP: 991 Shortbred genes ~  (filter zero/or not and then log2-transformed)~phenotype+age+sex+read.counts+DNA concentration
# Habrok
setwd("/scratch/p303998/SV_MWAS/Rdata_1216/")
# setwd("~/Downloads/2023_09_02/")
# setwd("/groups/umcg-fu/tmp01/users/umcg-yzh/R/2023_09_13/")
# functions
source(paste0("/scratch/p303998/SV_MWAS/Rdata_1216/Step1/", "Part1_functions.R"))
#source(paste0("~/Downloads/2023_09_02/", "Part1_functions.R"))
#source(paste0("/groups/umcg-fu/tmp01/users/umcg-yzh/R/2023_09_13/", "Part1_functions.R"))

#### 0. import data ####
setwd("/scratch/p303998/SV_MWAS/Rdata_0828/")
#setwd("~/Downloads/2023_09_02/")
#setwd("/groups/umcg-fu/tmp01/users/umcg-yzh/R/2023_09_13/")
# SV data
load("./SV_rawdata/dsgv_full.RData")
load("./SV_rawdata/vsgv_full.RData")
# abundance data
load("./abundance_rawdata/s_abun.RData")

# sample number data
para <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/running_list_higher_30.rds")
# para <- readRDS("~/Downloads/2023_09_02/running_list_higher_30.rds")
# para <- readRDS("/groups/umcg-fu/tmp01/users/umcg-yzh/R/2023_09_13/running_list_higher_30.rds")
sample_number <- read.delim("./Running_list/sampleCohorts.id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
sample_number$Cohort_2=NA
sample_number$Cohort_2[sample_number$Cohort==c("300OB")]=c("300ob")
sample_number$Cohort_2[sample_number$Cohort==c("300TZFG")]=c("300tzfg")
sample_number$Cohort_2[sample_number$Cohort==c("500FG_FSK")]=c("500fg_fsk")
sample_number$Cohort_2[sample_number$Cohort==c("LLD1")]=c("lld1")
sample_number$Cohort_2[sample_number$Cohort==c("DMP")]=c("dmp")
# physical score
#physical_score <- read.delim("~/Downloads/2023_09_02/202312_physical_activity_DAG3_LLD_To-Yue.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
physical_score <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/202312_physical_activity_DAG3_LLD_To-Yue.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
physical_score$ID=physical_score$DAG3_sampleID
physical_score$ID[is.na(physical_score$ID)]=physical_score$LLDEEP_ID[is.na(physical_score$ID)]
row.names(physical_score)=physical_score$ID
physical_score=physical_score[row.names(physical_score)%in%c(sample_number$New_SV_ID),]
# Frailty index
#Frailty_index <- read.delim("~/Downloads/2023_09_02/Frailty_index.csv",sep=",",row.names = NULL,header = T,check.names = F,fill = F)
Frailty_index <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/Frailty_index.csv",sep=",",row.names = NULL,header = T,check.names = F,fill = F)
#trans_ID <- read.delim("~/Downloads/2023_09_02/key_DAG3_pseudo_id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
trans_ID <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/key_DAG3_pseudo_id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
Frailty_index=merge(Frailty_index,trans_ID[,1:3],by.x = "PROJECT_PSEUDO_ID",by.y = "PROJECT_PSEUDO_ID",all.x=T)
Frailty_index=Frailty_index[!is.na(Frailty_index$DAG3_sampleID),c("FI64","FI60","FI41_B","FI41_FU","FI39_B","FI39_FU","DAG3_sampleID")]%>% `rownames<-`(.[,'DAG3_sampleID']) %>% dplyr::select(-'DAG3_sampleID')
# Phenotype data
load("./Phenotype_rawdata/full_phen.RData")
full_phen=removeColsAllNa(full_phen)
full_phen=merge(full_phen,physical_score[,c("total_scor_VAL","MVPA_scor_VAL")],by.x = "row.names",by.y = "row.names",all.x = T)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
full_phen=merge(full_phen,Frailty_index[,c("FI64","FI60","FI41_B","FI41_FU","FI39_B","FI39_FU")],by.x = "row.names",by.y = "row.names",all.x = T)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
Pheno_info <- read.delim("./Phenotype_info/phenInfoSumm.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
Pheno_info=Pheno_info[!Pheno_info$Class==c(""),]
Pheno_info=Pheno_info[-which(Pheno_info$UnifiedName==c("MeatFreqPerWeek")&Pheno_info$Unit==c("4 point scale")),]
# SV_info
dsgv_info_anno_ld <- read.delim("./SV_info/20230827_full_dsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
vsgv_info_anno_ld <- read.delim("./SV_info/20230827_full_vsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
all_sgv_info_anno_ld = rbind(dsgv_info_anno_ld,vsgv_info_anno_ld)
load("./SV_info/info.RData")
all_sv_info_anno <- readRDS("./SV_annotation/all_sv_info_anno.rds")
load("/scratch/p303998/SV_MWAS/Rdata_0828/SV_annotation/svAnnoDb.RData")
# pheno_cov
# read counts data
load("./readcounts_rawdata/full_read.RData")
full_read$log10_counts=log10(full_read$V2)
pheno_cov=merge(full_phen[,c("DNA.Concentration","Age","Sex","BristolType")],full_read[,c("V1","log10_counts"),drop=F],by.x = "row.names",by.y = "V1",all=F)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
# pheno_cov=pheno_cov[pheno_cov$Age>18,]
# vSV and dSV together
sgv_full=merge(dsgv_full,vsgv_full,by.x = "row.names",by.y = "row.names",all=F)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
# 109 SVs protein
SVs109__protein <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/SVs109__protein.rds")
# shortbred result
shortbred_all <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/shortbred_all.rds")

#Shortbred data
# all
shortbred_all<- do.call(rbind, list(Shortbred1, Shortbred2, Shortbred3, Shortbred4, Shortbred5, Shortbred6, Shortbred7, Shortbred8))
shortbred_all=as.data.frame(shortbred_all)
shortbred_all=shortbred_all[,-which(colSums(shortbred_all)==0)] # 8398 1005
saveRDS(shortbred_all,paste0("/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/","shortbred_all",".rds"))
shortbred_binary=shortbred_all
shortbred_binary[shortbred_binary>0]=1
saveRDS(shortbred_binary,paste0("/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/","shortbred_binary",".rds"))
#### 1. Shortbred data: NO transfer zero to NAs! Covariants not including SV and relative abundance ####
#### 1.1 shortbred result association ####
# Age
Shortbred_all <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/shortbred_all.rds")
cohort = c("dmp")
pheno_cohort = full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],c("Age","BMI")]
pheno_cov_abundance =NULL
pheno_cov_abundance = merge(pheno_cov, get(paste0("sp_re_", cohort, "_clr")), by.x = "row.names", by.y = "row.names", all = FALSE) %>%
  `rownames<-`(.[, 'Row.names']) %>%
  dplyr::select(-'Row.names')
covar=c("DNA.Concentration", "Sex", "log10_counts")
sample_name = Reduce(intersect, list(row.names(pheno_cohort), row.names(Shortbred_all), row.names(pheno_cov_abundance)))
result_age = shortbred_correct_noSV(Shortbred_all[sample_name,], pheno_cohort[sample_name,],sgv_full[sample_name,],
                                    pheno_cov_abundance[sample_name,], covar,SVs109__protein) 
for (i in c("Beta","SE","p","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate","fdr.p")){`result_age`[,i]=as.numeric(`result_age`[,i])}
result_age=result_age[!result_age$X==c("BMI"),]
# Sex
covar=c("DNA.Concentration", "Age", "log10_counts")
pheno_cohort = full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],c("Sex","BMI")]
result_sex = shortbred_correct_noSV(Shortbred_all[sample_name,], pheno_cohort[sample_name,],sgv_full[sample_name,],
                                    pheno_cov_abundance[sample_name,], covar,SVs109__protein) 
for (i in c("Beta","SE","p","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate","fdr.p")){`result_sex`[,i]=as.numeric(`result_sex`[,i])}
result_sex=result_sex[!result_sex$X==c("BMI"),]
result_sex_age=rbind(result_age,result_sex)
result_sex_age=merge(result_sex_age,SVs109__protein[,c("GeneName","Product.mixed","KoNumber","SV_Name","QueryId")],by.x = "Y", by.y = "GeneName", all.x = T)
saveRDS(result_sex_age, paste0("/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/", "result_sex_age","_no_SVabun_corrected_with_zero", ".rds"))

#### 1.2 Shortbred cut-phenotypes running ####
cohort = c("dmp")
pheno_cohort = full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],setdiff(unique(para$V2[para$V1==c("dmp")]),c("Age","Sex"))]
#### 1.3 Template for the script (Habrak) ####
script_template <- '
# Habrok
setwd("/scratch/p303998/SV_MWAS/Rdata_1216/")
# functions
source(paste0("/scratch/p303998/SV_MWAS/Rdata_1216/Step1/", "Part1_functions.R"))
#### 1. import data ####
setwd("/scratch/p303998/SV_MWAS/Rdata_0828/")
# SV data
load("./SV_rawdata/dsgv_full.RData")
load("./SV_rawdata/vsgv_full.RData")
# abundance data
load("./abundance_rawdata/s_abun.RData")
# sample number data
para <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/running_list_higher_30.rds")
sample_number <- read.delim("./Running_list/sampleCohorts.id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
sample_number$Cohort_2=NA
sample_number$Cohort_2[sample_number$Cohort==c("300OB")]=c("300ob")
sample_number$Cohort_2[sample_number$Cohort==c("300TZFG")]=c("300tzfg")
sample_number$Cohort_2[sample_number$Cohort==c("500FG_FSK")]=c("500fg_fsk")
sample_number$Cohort_2[sample_number$Cohort==c("LLD1")]=c("lld1")
sample_number$Cohort_2[sample_number$Cohort==c("DMP")]=c("dmp")
# Phenotype data
load("./Phenotype_rawdata/full_phen.RData")
full_phen=removeColsAllNa(full_phen)
Pheno_info <- read.delim("./Phenotype_info/phenInfoSumm.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
Pheno_info=Pheno_info[!Pheno_info$Class==c(""),]
Pheno_info=Pheno_info[-which(Pheno_info$UnifiedName==c("MeatFreqPerWeek")&Pheno_info$Unit==c("4 point scale")),]
# physical score
physical_score <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/202312_physical_activity_DAG3_LLD_To-Yue.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
physical_score$ID=physical_score$DAG3_sampleID
physical_score$ID[is.na(physical_score$ID)]=physical_score$LLDEEP_ID[is.na(physical_score$ID)]
row.names(physical_score)=physical_score$ID
physical_score=physical_score[row.names(physical_score)%in%c(sample_number$New_SV_ID),]
# Frailty index
Frailty_index <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/Frailty_index.csv",sep=",",row.names = NULL,header = T,check.names = F,fill = F)
trans_ID <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/key_DAG3_pseudo_id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
Frailty_index=merge(Frailty_index,trans_ID[,1:3],by.x = "PROJECT_PSEUDO_ID",by.y = "PROJECT_PSEUDO_ID",all.x=T)
Frailty_index=Frailty_index[!is.na(Frailty_index$DAG3_sampleID),c("FI64","FI60","FI41_B","FI41_FU","FI39_B","FI39_FU","DAG3_sampleID")]%>% `rownames<-`(.[,"DAG3_sampleID"]) %>% dplyr::select(-"DAG3_sampleID")
# Phenotype data
load("./Phenotype_rawdata/full_phen.RData")
full_phen=removeColsAllNa(full_phen)
full_phen=merge(full_phen,physical_score[,c("total_scor_VAL","MVPA_scor_VAL")],by.x = "row.names",by.y = "row.names",all.x = T)%>% `rownames<-`(.[,"Row.names"]) %>% dplyr::select(-"Row.names")
full_phen=merge(full_phen,Frailty_index[,c("FI64","FI60","FI41_B","FI41_FU","FI39_B","FI39_FU")],by.x = "row.names",by.y = "row.names",all.x = T)%>% `rownames<-`(.[,"Row.names"]) %>% dplyr::select(-"Row.names")
# read counts data
load("./readcounts_rawdata/full_read.RData")
# SV_info
dsgv_info_anno_ld <- read.delim("./SV_info/20230827_full_dsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
vsgv_info_anno_ld <- read.delim("./SV_info/20230827_full_vsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
all_sgv_info_anno_ld = rbind(dsgv_info_anno_ld,vsgv_info_anno_ld)
load("./SV_info/info.RData")
# pheno_cov
full_read$log10_counts=log10(full_read$V2)
pheno_cov=merge(full_phen[,c("DNA.Concentration", "Age", "Sex","BristolType")], full_read[,c("V1", "log10_counts"), drop=FALSE], by.x = "row.names", by.y = "V1", all=FALSE) %>% `rownames<-`(.[, "Row.names"]) %>% dplyr::select(-"Row.names")
# vSV and dSV together
sgv_full=merge(dsgv_full, vsgv_full, by.x = "row.names", by.y = "row.names", all=FALSE) %>% `rownames<-`(.[, "Row.names"]) %>% dplyr::select(-"Row.names")
# shortbred
Shortbred_all <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/shortbred_all.rds")
SVs109__protein <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/SVs109__protein.rds")

#### 1. relative abundance ####
#### DMP ####
sp_re_dmp <- s_abun[sample_number$New_SV_ID[sample_number$Cohort==c("DMP")],]
sp_re_dmp_clr=do_clr_externalWeighting(t(sp_re_dmp),t(sp_re_dmp)) %>% data.frame(.)  %>% t(.); sp_re_dmp_clr=as.data.frame(sp_re_dmp_clr)

#### dmp ####
cohort = c("dmp")
pheno_cohort = full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],setdiff(unique(para$V2[para$V1==c("dmp")]),c("Age","Sex"))]
pheno_cov_abundance =NULL
pheno_cov_abundance = merge(pheno_cov, get(paste0("sp_re_", cohort, "_clr")), by.x = "row.names", by.y = "row.names", all = FALSE) %>%
    `rownames<-`(.[, "Row.names"]) %>%
    dplyr::select(-"Row.names")
sample_name = Reduce(intersect, list(row.names(pheno_cohort), row.names(Shortbred_all), row.names(pheno_cov_abundance),row.names(sgv_full)))
covar = c("DNA.Concentration", "Sex", "log10_counts", "Age")
result = shortbred_correct_noSV(Shortbred_all[sample_name,], pheno_cohort[sample_name,{NUM_RANGE}],sgv_full[sample_name,],
                          pheno_cov_abundance[sample_name,], covar,SVs109__protein) 
for (i in c("Beta","SE","p","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate","fdr.p")){`result`[,i]=as.numeric(`result`[,i])}
saveRDS(result,paste0("/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/batch/no_SVabun_corrected/with_zeros/","_range_{NUM_RANGE}.rds"))
'
#### 1.4 generate scripts ####
generate_scripts <- function(start, end, step) {
  for (i in seq(start, end, step)) {
    num_range <- paste(i, min(i + step - 1, end), sep = ":")
    script_content <- gsub("\\{NUM_RANGE\\}", num_range, script_template)
    script_name <- paste0("analysis_script_", num_range, ".R")
    writeLines(script_content, con = script_name)
  }
}
setwd("/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/scripts/no_SVabun_corrected/with_zeros/")
generate_scripts(1, ncol(pheno_cohort), 3)
#### 1.5 linux submit ####
# cd /scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/scripts/no_SVabun_corrected/with_zeros/
# vim submit_jobs.sh
#!/bin/bash
# SCRIPT_DIR="/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/scripts/no_SVabun_corrected/with_zeros/" #####!!!!! Don't forget to change the folder name.
# module load RPlus
# for script in $SCRIPT_DIR/*.R
# do
# ls $script
# sbatch --job-name=myRjob --output=output_%j.txt --time=01:00:00 --cpus-per-task=4 --mem=10G --wrap="module load RPlus; Rscript $script"
# done
# chmod +x submit_jobs.sh

#### 1.6 rbind all data ####
result_sex_age_no_SVabun_corrected_with_zero <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/result_sex_age_no_SVabun_corrected_with_zero.rds")
setpath <- "/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/batch/no_SVabun_corrected/with_zeros/"
list = list.files(setpath,pattern = '.rds',recursive = TRUE,full.names = TRUE)%>%unique()
dataframes_list <- lapply(list, readRDS)
combined_dataframe <- do.call(rbind, dataframes_list)
combined_dataframe=merge(combined_dataframe,SVs109__protein[,c("GeneName","Product.mixed","KoNumber","SV_Name","QueryId")],by.x = "Y", by.y = "GeneName", all.x = T)
combined_dataframe=rbind(combined_dataframe,result_sex_age_no_SVabun_corrected_with_zero)
saveRDS(combined_dataframe,paste0("/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/","shortbred_result_dmp_no_SVabun_corrected_with_zero",".rds"))
# filter zero and rates
shortbred_result_dmp_no_SVabun_corrected_with_zero <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/shortbred_result_dmp_no_SVabun_corrected_with_zero.rds")
shortbred_result_with_zero=shortbred_result_dmp_no_SVabun_corrected_with_zero[!is.na(shortbred_result_dmp_no_SVabun_corrected_with_zero$p),]
shortbred_result_with_zero=shortbred_result_with_zero[!shortbred_result_with_zero$X%in%c("FI64","FI60","FI39_B","FI39_FU",
                                                                                         "BristolType","BristolFreq","CleanReadCount","Latitude","Longitude","FattyLiverIndexT1Class","DNA.Concentration"),]
shortbred_result_with_zero=shortbred_result_with_zero[shortbred_result_with_zero$N>80,]
shortbred_result_with_zero=shortbred_result_with_zero[-which(shortbred_result_with_zero$x_uniq_N==2&shortbred_result_with_zero$x_non_zero_rate<0.1),]
shortbred_result_with_zero=shortbred_result_with_zero[-which(shortbred_result_with_zero$x_uniq_N==2&shortbred_result_with_zero$x_non_zero_rate>0.9),]
shortbred_result_with_zero$fdr.p=p.adjust(shortbred_result_with_zero$p,method = "fdr")
shortbred_result_with_zero = merge(shortbred_result_with_zero, Cazy_file, by.x = "Y", by.y = "PROTEIN_ID", all.x = T)
shortbred_result_with_zero = merge(shortbred_result_with_zero, GO_file, by.x = "QueryId", by.y = "PROTEIN_ID", all.x = T)
# final results
saveRDS(shortbred_result_with_zero,paste0("/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/","shortbred_result_with_zero_filtered",".rds"))
write.table(shortbred_result_with_zero, paste0("/scratch/p303998/SV_MWAS/Rdata_0828/","shortbred_no_SVabun",".tsv"),sep = "\t", quote = F,col.names = T,row.names = F)
#### 1.7 summary: Figure 3 ####
shortbred_with_zero <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/shortbred_result_with_zero_filtered.rds")
for (i in c("Y","X")){`shortbred_with_zero`[,i]=as.character(`shortbred_with_zero`[,i])}
length(unique(shortbred_with_zero$Y))#1004
length(unique(shortbred_with_zero$X))#113 phenotypes
#### 1.8 Volcano plot for Shortbred genes associated with age ####
#shortbred_with_zero <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/shortbred_result/shortbred_result_with_zero_filtered.rds")
shortbred_with_zero <- readRDS("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/shortbred_result_with_zero_filtered.rds")
for (i in c("Y","X")){`shortbred_with_zero`[,i]=as.character(`shortbred_with_zero`[,i])}
length(unique(shortbred_with_zero$Y))#1004
length(unique(shortbred_with_zero$X))#113 phenotypes
shortbred_fdr=0.05/(1004*113)
shortbred_with_zero=shortbred_with_zero[shortbred_with_zero$X==c("Age"),]
shortbred_with_zero$replication=NA
shortbred_with_zero$replication[shortbred_with_zero$p<shortbred_fdr&shortbred_with_zero$Beta>0&shortbred_with_zero$X==c("Age")]=c("significant-positeive")
shortbred_with_zero$replication[shortbred_with_zero$p<shortbred_fdr&shortbred_with_zero$Beta<0&shortbred_with_zero$X==c("Age")]=c("significant-negative")
shortbred_with_zero$replication[shortbred_with_zero$p>shortbred_fdr]=c("non-significant")
shortbred_with_zero$replication=as.factor(shortbred_with_zero$replication)
shortbred_with_zero$replication=factor(shortbred_with_zero$replication,levels = c("significant-positeive","significant-negative","non-significant"))
library(ggplot2) ##绘图使用
library(ggprism) ##设置主题私用
library(ggrepel) ##给火山图加标签使用
ggplot(shortbred_with_zero, aes(x =Beta, y= -log10(p), color=replication)) +
  geom_point(alpha=1, size=2.5) +
  scale_color_manual(values=c('brown','steelblue','gray')) +
  xlim(c(min(shortbred_with_zero$Beta), max(shortbred_with_zero$Beta))) +  ##调整x轴的取值范围，可以根据max(abs(shortbred_with_zero$Beta))，获得差异基因最大值是多少，再进行取舍
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05/(1004*113)), lty=4,col="black",lwd=0.8) + 
  labs(x="Beta effect", y="-log10(P)") +
  ggtitle("Replication in LLD cohort") + 
  theme(plot.title = element_text(hjust = 0.5),legend.position="right",legend.title = element_blank())+
  theme_prism(border = T)
