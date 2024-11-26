# Association analysis in 5 cohorts (model:SV~age+sex+DNA.concentrate+readcounts+species abundance)
# Habrok
# setwd("/scratch/p303998/SV_MWAS/Rdata_1216/")
# setwd("~/Downloads/2023_09_02/")
# setwd("/groups/umcg-fu/tmp01/users/umcg-yzh/R/2023_09_13/")
# functions
#source(paste0("/scratch/p303998/SV_MWAS/Rdata_1216/Step1/", "Part1_functions.R"))
#source(paste0("~/Downloads/2023_09_02/", "Part1_functions.R"))
#source(paste0("/groups/umcg-fu/tmp01/users/umcg-yzh/R/2023_09_13/", "Part1_functions.R"))
#### 1.0 import data ####
setwd("/scratch/p303998/SV_MWAS/Rdata_0828/")
#setwd("~/Downloads/2023_09_02/")
#setwd("/groups/umcg-fu/tmp01/users/umcg-yzh/R/2023_09_13/")
# SV data
load("./SV_rawdata/dsgv_full.RData")
load("./SV_rawdata/vsgv_full.RData")
# abundance data
load("./abundance_rawdata/s_abun.RData")
# sample number data
para <- readRDS("./running_list_higher_30.rds")
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
# read counts data
load("./readcounts_rawdata/full_read.RData")
# SV_info
dsgv_info_anno_ld <- read.delim("./SV_info/20230827_full_dsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
vsgv_info_anno_ld <- read.delim("./SV_info/20230827_full_vsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
all_sgv_info_anno_ld = rbind(dsgv_info_anno_ld,vsgv_info_anno_ld)
load("./SV_info/info.RData")
all_sv_info_anno <- readRDS("./SV_annotation/all_sv_info_anno.rds")
#load("/scratch/p303998/SV_MWAS/Rdata_0828/SV_annotation/svAnnoDb.RData")
load("./SV_annotation/svAnnoDb.RData")
# pheno_cov
full_read$log10_counts=log10(full_read$V2)
pheno_cov=merge(full_phen[,c("DNA.Concentration","Age","Sex","BristolType")],full_read[,c("V1","log10_counts"),drop=F],by.x = "row.names",by.y = "V1",all=F)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
# pheno_cov=pheno_cov[pheno_cov$Age>18,]
# vSV and dSV together
sgv_full=merge(dsgv_full,vsgv_full,by.x = "row.names",by.y = "row.names",all=F)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
#### 2. relative abundance ####
#### DMP ####
sp_re_dmp <- s_abun[sample_number$New_SV_ID[sample_number$Cohort==c("DMP")],]
sp_re_dmp_clr=do_clr_externalWeighting(t(sp_re_dmp),t(sp_re_dmp)) %>% data.frame(.)  %>% t(.); sp_re_dmp_clr=as.data.frame(sp_re_dmp_clr)
#### LLD-baseline ####
sp_re_lld_baseline <- na.omit(s_abun[sample_number$New_SV_ID[sample_number$Cohort==c("LLD1")],])
sp_re_lld1_clr=do_clr_externalWeighting(t(sp_re_lld_baseline),t(sp_re_lld_baseline)) %>% data.frame(.)  %>% t(.); sp_re_lld1_clr=as.data.frame(sp_re_lld1_clr)
#### LLD-follow up-FSK ####
sp_re_lld2_FSK <- na.omit(s_abun[sample_number$New_SV_ID[sample_number$Cohort==c("LLD2_FSK")],])
sp_re_lld2_FSK_clr=do_clr_externalWeighting(t(sp_re_lld2_FSK),t(sp_re_lld2_FSK)) %>% data.frame(.)  %>% t(.); sp_re_lld2_FSK_clr=as.data.frame(sp_re_lld2_FSK_clr)
#### 500FG-FSK ####
sp_re_500FG_FSK <- na.omit(s_abun[sample_number$New_SV_ID[sample_number$Cohort==c("500FG_FSK")],])
sp_re_500fg_fsk_clr=do_clr_externalWeighting(t(sp_re_500FG_FSK),t(sp_re_500FG_FSK)) %>% data.frame(.)  %>% t(.); sp_re_500fg_fsk_clr=as.data.frame(
  
)
#### 300OB ####
sp_re_300OB <- na.omit(s_abun[sample_number$New_SV_ID[sample_number$Cohort==c("300OB")],])
sp_re_300ob_clr=do_clr_externalWeighting(t(sp_re_300OB),t(sp_re_300OB)) %>% data.frame(.)  %>% t(.); sp_re_300ob_clr=as.data.frame(sp_re_300ob_clr)
#### 300TZFG ####
sp_re_300TZFG <- na.omit(s_abun[sample_number$New_SV_ID[sample_number$Cohort==c("300TZFG")],])
sp_re_300tzfg_clr=do_clr_externalWeighting(t(sp_re_300TZFG),t(sp_re_300TZFG)) %>% data.frame(.)  %>% t(.); sp_re_300tzfg_clr=as.data.frame(sp_re_300tzfg_clr)
#### 3.0 run association with age for all 5 cohorts in parallel ####
no_cores <- 6  # Leave one core free for system stability
registerDoParallel(cores = no_cores)
foreach(cohort = c("dmp","300ob","300tzfg","500fg_fsk","lld1"), .combine = 'rbind', .packages = c('dplyr'),.verbose = TRUE) %dopar% {
  running_info_cohort = para[para$V1 == cohort,]
  pheno_cohort = full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],]
  pheno_cov_abundance = merge(pheno_cov, get(paste0("sp_re_", cohort, "_clr")), by.x = "row.names", by.y = "row.names", all = FALSE) %>%
    `rownames<-`(.[, 'Row.names']) %>%
    dplyr::select(-'Row.names')
  sample_name = Reduce(intersect, list(row.names(pheno_cohort), row.names(sgv_full), row.names(pheno_cov_abundance)))
  
  # Determine covariates based on cohort
  covar = switch(cohort,
                 "dmp" = c("DNA.Concentration", "Sex", "log10_counts"),
                 "300ob" = c("Sex", "log10_counts"),
                 "300tzfg" = c("Sex", "log10_counts"),
                 "500fg_fsk" = c("Sex", "log10_counts"),
                 "lld1" = c("Sex", "log10_counts"))
  
  result = SV_lm_glm(pheno_cohort[sample_name, "Age",drop=F], sgv_full[sample_name,],
                     pheno_cov_abundance[sample_name,], covar, running_info_cohort,cohort)
  saveRDS(result, paste0("/scratch/p303998/SV_MWAS/Rdata_1216/cohort_result/Age/",cohort, ".rds"))
}
#### 5. Filter the results according to the zero rate ####
# dmp
dmp$p=as.numeric(dmp$p)
dmp$p[which(dmp$p==0)]=NA
for (i in c("X","Y")){dmp[,i]=as.character(dmp[,i])}
for (i in c("Beta","SE","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate")){dmp[,i]=as.numeric(dmp[,i])}
dmp$species=sapply(strsplit(as.character(dmp$Y), "\\:"), "[", 1)
dmp$species[grep("Phascolarctobacterium sp. CAG:207",dmp$Y)]=c("Phascolarctobacterium sp. CAG:207")
dmp$species[grep("Phascolarctobacterium sp. CAG:266",dmp$Y)]=c("Phascolarctobacterium sp. CAG:266")
# one level ---> 0
dmp$Beta[which(dmp$Beta==c("one level"))]=NA
dmp$SE[which(dmp$SE==c("one level"))]=NA
# dSV level ---> 0.1<zero number<0.9
dmp$p[dmp$N<sum(sample_number$Cohort_2==c("dmp"),na.rm = T)*0.01]=NA
dmp$p[dmp$y_uniq_N==2&dmp$y_non_zero_rate>0.9]=NA
dmp$p[dmp$y_uniq_N==2&dmp$y_non_zero_rate<0.1]=NA
dmp$fdr.p=p.adjust(dmp$p,method="fdr")
dmp=dmp[dmp$Y%in%dmp$Y[dmp$X==c("Age")],]

# lld1
lld1$p=as.numeric(lld1$p)
lld1$p[which(lld1$p==0)]=NA
lld1$species=sapply(strsplit(as.character(lld1$Y), "\\:"), "[", 1)
lld1$species[grep("Phascolarctobacterium sp. CAG:207",lld1$Y)]=c("Phascolarctobacterium sp. CAG:207")
lld1$species[grep("Phascolarctobacterium sp. CAG:266",lld1$Y)]=c("Phascolarctobacterium sp. CAG:266")
for (i in c("X","Y")){lld1[,i]=as.character(lld1[,i])}
for (i in c("Beta","SE","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate")){lld1[,i]=as.numeric(lld1[,i])}
# one level ---> 0
lld1$Beta[which(lld1$Beta==c("one level"))]=NA
lld1$SE[which(lld1$SE==c("one level"))]=NA
# dSV level ---> 0.1<zero number<0.9
lld1$p[lld1$N<sum(sample_number$Cohort_2==c("lld1"),na.rm = T)*0.05]=NA
lld1$p[lld1$y_uniq_N==2&lld1$y_non_zero_rate>0.9]=NA
lld1$p[lld1$y_uniq_N==2&lld1$y_non_zero_rate<0.1]=NA
lld1$fdr.p=p.adjust(lld1$p,method="fdr")

# 300ob
`300ob`$p=as.numeric(`300ob`$p)
`300ob`$p[which(`300ob`$p==0)]=NA
`300ob`$species=sapply(strsplit(as.character(`300ob`$Y), "\\:"), "[", 1)
`300ob`$species[grep("Phascolarctobacterium sp. CAG:207",`300ob`$Y)]=c("Phascolarctobacterium sp. CAG:207")
`300ob`$species[grep("Phascolarctobacterium sp. CAG:266",`300ob`$Y)]=c("Phascolarctobacterium sp. CAG:266")
for (i in c("X","Y")){`300ob`[,i]=as.character(`300ob`[,i])}
for (i in c("Beta","SE","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate")){`300ob`[,i]=as.numeric(`300ob`[,i])}
# one level ---> 0
`300ob`$Beta[which(`300ob`$Beta==c("one level"))]=NA
`300ob`$SE[which(`300ob`$SE==c("one level"))]=NA
# dSV level ---> 0.1<zero number<0.9
`300ob`$p[`300ob`$N<sum(sample_number$Cohort_2==c("300ob"),na.rm = T)*0.1]=NA
`300ob`$p[`300ob`$y_uniq_N==2&`300ob`$y_non_zero_rate>0.9]=NA
`300ob`$p[`300ob`$y_uniq_N==2&`300ob`$y_non_zero_rate<0.1]=NA
`300ob`$fdr.p=p.adjust(`300ob`$p,method="fdr")

# 300tzfg
`300tzfg`$p=as.numeric(`300tzfg`$p)
`300tzfg`$p[which(`300tzfg`$p==0)]=NA
`300tzfg`$species=sapply(strsplit(as.character(`300tzfg`$Y), "\\:"), "[", 1)
`300tzfg`$species[grep("Phascolarctobacterium sp. CAG:207",`300tzfg`$Y)]=c("Phascolarctobacterium sp. CAG:207")
`300tzfg`$species[grep("Phascolarctobacterium sp. CAG:266",`300tzfg`$Y)]=c("Phascolarctobacterium sp. CAG:266")
for (i in c("X","Y")){`300tzfg`[,i]=as.character(`300tzfg`[,i])}
for (i in c("Beta","SE","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate")){`300tzfg`[,i]=as.numeric(`300tzfg`[,i])}
# one level ---> 0
`300tzfg`$Beta[which(`300tzfg`$Beta==c("one level"))]=NA
`300tzfg`$SE[which(`300tzfg`$SE==c("one level"))]=NA
# dSV level ---> 0.1<zero number<0.9
`300tzfg`$p[`300tzfg`$N<sum(sample_number$Cohort_2==c("300tzfg"),na.rm = T)*0.1]=NA
`300tzfg`$p[`300tzfg`$y_uniq_N==2&`300tzfg`$y_non_zero_rate>0.9]=NA
`300tzfg`$p[`300tzfg`$y_uniq_N==2&`300tzfg`$y_non_zero_rate<0.1]=NA
`300tzfg`$fdr.p=p.adjust(`300tzfg`$p,method="fdr")


# 500fg_fsk
`500fg_fsk`$p=as.numeric(`500fg_fsk`$p)
`500fg_fsk`$p[which(`500fg_fsk`$p==0)]=NA
`500fg_fsk`$species=sapply(strsplit(as.character(`500fg_fsk`$Y), "\\:"), "[", 1)
`500fg_fsk`$species[grep("Phascolarctobacterium sp. CAG:207",`500fg_fsk`$Y)]=c("Phascolarctobacterium sp. CAG:207")
`500fg_fsk`$species[grep("Phascolarctobacterium sp. CAG:266",`500fg_fsk`$Y)]=c("Phascolarctobacterium sp. CAG:266")
for (i in c("X","Y")){`500fg_fsk`[,i]=as.character(`500fg_fsk`[,i])}
for (i in c("Beta","SE","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate")){`500fg_fsk`[,i]=as.numeric(`500fg_fsk`[,i])}
# one level ---> 0
`500fg_fsk`$Beta[which(`500fg_fsk`$Beta==c("one level"))]=NA
`500fg_fsk`$SE[which(`500fg_fsk`$SE==c("one level"))]=NA
# dSV level ---> 0.1<zero number<0.9
`500fg_fsk`$p[`500fg_fsk`$N<sum(sample_number$Cohort_2==c("500fg_fsk"),na.rm = T)*0.1]=NA
`500fg_fsk`$p[`500fg_fsk`$y_uniq_N==2&`500fg_fsk`$y_non_zero_rate>0.9]=NA
`500fg_fsk`$p[`500fg_fsk`$y_uniq_N==2&`500fg_fsk`$y_non_zero_rate<0.1]=NA
`500fg_fsk`$fdr.p=p.adjust(`500fg_fsk`$p,method="fdr")

#### 6. save the final files ####
for( i in c("dmp","300ob","300tzfg","500fg_fsk","lld1")){
  print(i)
  saveRDS(get(i),paste0("/scratch/p303998/SV_MWAS/Rdata_1216/Step1/Results/All/Filter_zero/","filter_zero_",i,".rds"))
}
#### 7. Read in files ####
object_names <- c("dmp", "300ob", "300tzfg", "500fg_fsk", "lld1")
directory <- "/scratch/p303998/SV_MWAS/Rdata_1216/Step1/Results/All/filter_zero/"
data_list <- list()
for (name in object_names) {
  file_path <- paste0(directory, "filter_zero_", name, ".rds")
  data_list[[name]] <- readRDS(file_path)
}
all_sv_info_anno <- readRDS("./SV_annotation/all_sv_info_anno.rds")
dsgv_info_anno_ld <- read.delim("./SV_info/20230827_full_dsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
vsgv_info_anno_ld <- read.delim("./SV_info/20230827_full_vsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
all_sgv_info_anno_ld = rbind(dsgv_info_anno_ld,vsgv_info_anno_ld)
#### 8. All SVs: Age-SVs FDR threshold ####
filter_zero_dmp <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/Step1/Results/All/filter_zero/filter_zero_dmp.rds")
dmp_age=filter_zero_dmp%>%na.omit(.)
unique_ld_counts_dSV_age <- dsgv_info_anno_ld[dsgv_info_anno_ld$SV_Name%in%dmp_age$Y,] %>% dplyr::group_by(organism) %>% dplyr::summarise(unique_ld_count = dplyr::n_distinct(LD_block))
unique_ld_counts_vSV_age <- vsgv_info_anno_ld[vsgv_info_anno_ld$SV_Name%in%dmp_age$Y,] %>% dplyr::group_by(organism) %>% dplyr::summarise(unique_ld_count = dplyr::n_distinct(LD_block))
unique_ld_counts_age=merge(unique_ld_counts_dSV_age,unique_ld_counts_vSV_age,by.x = "organism",by.y = "organism",all=T)
colnames(unique_ld_counts_age)[2:3]=c("dSV","vSV")
unique_ld_counts_age$ld_sum=unique_ld_counts_age$dSV+unique_ld_counts_age$vSV
unique_ld_counts_age$organism%in%unique(dmp_age$species)
dSV_ld_sum_age=sum(unique_ld_counts_age$dSV)
vSV_ld_sum_age=sum(unique_ld_counts_age$vSV)
SV_ld_sum_age=sum(unique_ld_counts_age$ld_sum)
fdr_threshold_age=0.05/SV_ld_sum_age
#### 9. replication in other 4 cohorts ####
#### 9.1 merge all files ####
# Reorder the list to make sure 'dmp' is the first element
dmp_age_sig=dmp_age[dmp_age$p<fdr_threshold_age,]
data_list[["dmp"]]=dmp_age_sig
data_list <- data_list[c("dmp", "300ob", "300tzfg", "500fg_fsk", "lld1")]
for(i in seq_along(data_list)) {
  # Rename columns except for the first two (assuming they are "X" and "Y")
  colnames(data_list[[i]]) <- ifelse(1:ncol(data_list[[i]]) <= 2, 
                                     colnames(data_list[[i]]), 
                                     paste(colnames(data_list[[i]]), names(data_list)[i], sep="_"))
}
merged_data <- Reduce(function(x, y) merge(x, y, by = c("X", "Y"), all.x = TRUE), data_list)
merged_data$Beta_300ob[which(is.na(merged_data$p_300ob))]=NA
merged_data$SE_300ob[which(is.na(merged_data$p_300ob))]=NA
merged_data$Beta_300tzfg[which(is.na(merged_data$p_300tzfg))]=NA
merged_data$SE_300tzfg[which(is.na(merged_data$p_300tzfg))]=NA
merged_data$Beta_500fg_fsk[which(is.na(merged_data$p_500fg_fsk))]=NA
merged_data$SE_500fg_fsk[which(is.na(merged_data$p_500fg_fsk))]=NA
merged_data$Beta_300tzfg[which(is.na(merged_data$p_300tzfg))]=NA
merged_data$SE_300tzfg[which(is.na(merged_data$p_300tzfg))]=NA
merged_data$Beta_lld1[which(is.na(merged_data$p_lld1))]=NA
merged_data$SE_lld1[which(is.na(merged_data$p_lld1))]=NA

#### 9.2 replication ####
replication_result=data.frame(SV=NA,p_other_cohorts=NA,"300ob"=NA,"300tzfg"=NA,"500fg_fsk"=NA,"lld1"=NA)
for (i in 1:nrow(merged_data)){
  print(i)
  replication_result[i,"SV"]=merged_data$Y[i]
  for (j in c("300ob","300tzfg","500fg_fsk","lld1")){
    print(j)
    # j=c("300ob")
    Beta_name=paste0("Beta_",j)
    P_name=paste0("p_",j)
    if (!j==c("lld1")){j=paste0("X",j)}
    if (is.na(merged_data[i,P_name])){merged_data[i,P_name]=1}
    if (merged_data[i,Beta_name]*merged_data[i,"Beta_dmp"]>0&merged_data[i,P_name]<0.05){replication_result[i,j]=1}else{replication_result[i,j]=0}
  }
}
replication_result$p_other_cohorts=replication_result$X300ob+replication_result$X300tzfg+replication_result$X500fg_fsk+replication_result$lld1
dmp_age_sig=dmp_age_sig[dmp_age_sig$Y%in%replication_result$SV[replication_result$p_other_cohorts>0],]
saveRDS(dmp_age_sig, paste0("/scratch/p303998/SV_MWAS/Rdata_1216/","meta_age_significant",".rds"))
saveRDS(replication_result, paste0("/scratch/p303998/SV_MWAS/Rdata_1216/","replication_result",".rds"))

#### 6. Meta-analysis for 300tzfg and 3000OB ####
i=c("300tzfg")
SV_name=replication_result$SV[replication_result[,grep(i,colnames(replication_result))]==1]
test=merged_data[merged_data$Y%in%SV_name,c("Beta_dmp","SE_dmp","Beta_300tzfg","SE_300tzfg","Y","X")]%>%na.omit(.)
meta_300TZFG=my_batch_meta_lm(test,c("DMP", "300TZFG"),c(1,3),c(2,4),row_var_col = 5,col_var_col = 6)
meta_300TZFG=meta_300TZFG[["table"]]

i=c("300ob")
SV_name=replication_result$SV[replication_result[,grep(i,colnames(replication_result))]==1]
test=merged_data[merged_data$Y%in%SV_name,c("Beta_dmp","SE_dmp","Beta_300ob","SE_300ob","Y","X")]%>%na.omit(.)
meta_300OB=my_batch_meta_lm(test,c("DMP", "300OB"),c(1,3),c(2,4),row_var_col = 5,col_var_col = 6)
meta_300OB=meta_300TZFG[["table"]]





