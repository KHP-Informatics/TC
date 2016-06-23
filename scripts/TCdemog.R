mmseTS <- vector("list",5)

camdTS <- mmse_both_cpath

mmseTS[[1]] <- mmse_aac
mmseTS[[2]] <- mmseADNIts
mmseTS[[3]] <- mmseAIBLts
mmseTS[[4]] <- mmse1013ts
mmseTS[[5]] <- mmse1014ts

mmse_demog <- data.frame(combined=NA,adni=NA,aibl=NA,c1013=NA,c1014=NA)

rownames(mmse_demog) <- "Num ind any data"

for (i in 1:5){
  
  mmse_demog["Num ind > 1 TP",i] <- mmseTS[[i]]@numTS

  tmp_lengths <- unlist(lapply(mmseTS[[i]]@tps,length))
  
  mmse_demog["Median num TP",i] <- median(tmp_lengths)

  mmse_demog["LQ num TP",i] <- as.numeric(quantile(unlist(lapply(mmseTS[[i]]@tps,length)),0.25))

  mmse_demog["UQ TP",i] <- as.numeric(quantile(unlist(lapply(mmseTS[[i]]@tps,length)),0.75))
  
  last_tp <- c()
  
  for (j in 1:length(tmp_lengths)){
   
    last_tp[j] <- mmseTS[[i]]@tps[[j]][tmp_lengths[j]]
    
  }
  
  mmse_demog["Median years",i] <- signif(median(last_tp) / 365.25,2)

  mmse_demog["LQ years",i] <- signif(as.numeric(quantile(last_tp,0.25))/ 365.25,2)

  mmse_demog["UQ years",i] <- signif(as.numeric(quantile(last_tp,0.75))/ 365.25,2)
  
  
}

setwd("../data")

# ADNI

setwd("adni_sep_2015/")

diag <- read.csv("DXSUM_PDXCONV_ADNIALL.csv",header=T,as.is=T)

adni_bl_diag <- array(dim=c(length(mmseADNIts@subj_id),1))
adni_lv_diag <- array(dim=c(length(mmseADNIts@subj_id),1))
adni_bl_date <- list()

for (i in 1:length(mmseADNIts@subj_id)){
  
  ind <- which(diag[,"RID"] == mmseADNIts@subj_id[i])
  
  ind2 <- which(diag[ind,"VISCODE"]=="bl")
  
  adni_bl_diag[i] <- diag[ind[ind2],"DXCURREN"]
  
  adni_bl_date[[i]] <- as.Date(diag[ind[ind2],"EXAMDATE"])
  
  ind3 <- order(as.Date(diag[ind[!is.na(diag[ind,"DXCURREN"])],"EXAMDATE"]),decreasing = T)[1]
  
  adni_lv_diag[i] <- diag[ind[!is.na(diag[ind,"DXCURREN"])][ind3],"DXCURREN"]
  
}

mmse_demog[c("Num CTL","Num MCI","Num AD"),"adni"] <- table(adni_bl_diag)

demog <- read.csv("PTDEMOG.csv",header=T,as.is=T)

adni_df <- data.frame(age=NA,education=NA,gender=NA)

adni_years_educ <- numeric(length(mmseADNIts@subj_id))

for (i in 1:length(mmseADNIts@subj_id)){
  
  ind <- which(demog[,"RID"] == mmseADNIts@subj_id[i])
  
  dob <- as.Date(paste(demog[ind[1],"PTDOBYY"],demog[ind[1],"PTDOBMM"],1,sep="-"))
  
  adni_years_educ[i] <- demog[ind[1],"PTEDUCAT"]
  
  adni_df[i,] <- c(signif(as.numeric(adni_bl_date[[1]] - dob)/365.25,2),demog[ind[1],"PTEDUCAT"],demog[ind[1],"PTGENDER"])
  
}

mmse_demog[c("Num Male","Num Female"),"adni"] <- table(adni_df[,"gender"])

mmse_demog["Median bl age","adni"] <- median(adni_df[,"age"])

mmse_demog["LQ bl age","adni"] <- as.numeric(quantile(adni_df[,"age"],0.25))

mmse_demog["UQ bl age","adni"] <- as.numeric(quantile(adni_df[,"age"],0.75))

apoe <- read.csv("APOERES.csv",header=T,as.is=T)

adni_num_apoe <- array(dim=c(length(mmseADNIts@subj_id),1))

for (i in 1:length(mmseADNIts@subj_id)){
  
  ind <- which(apoe[,"RID"] == mmseADNIts@subj_id[i])
  
  adni_num_apoe[i] <- (apoe[ind[1],"APGEN1"] == 4) + (apoe[ind[1],"APGEN2"] == 4)
  
}

mmse_demog[c("Num 0 E4","Num 1 E4","Num 2 E4"),"adni"] <- table(adni_num_apoe)

setwd("..")



# AIBL

setwd("aibl_28apr_2015/")

#aibl_bl_date[[i]] <- as.Date(diag[ind[ind2],"EXAMDATE"])

visits <- read.csv("aibl_mmse_28-Apr-2015.csv",header=T,as.is=T)

aibl_bl_date <- list()

for (i in 1:length(mmseAIBLts@subj_id)){
  
  ind <- which(visits[,"RID"] == mmseAIBLts@subj_id[i])
  
  ind2 <- which(visits[ind,"VISCODE"]=="bl")
  
  aibl_bl_date[[i]] <- as.Date(visits[ind[ind2],"EXAMDATE"],"%m/%d/%Y")
  
}


diag <- read.csv("aibl_pdxconv_28-Apr-2015.csv",header=T,as.is=T)

aibl_bl_diag <- array(dim=c(length(mmseAIBLts@subj_id),1))
aibl_lv_diag <- array(dim=c(length(mmseAIBLts@subj_id),1))

for (i in 1:length(mmseAIBLts@subj_id)){
  
  ind <- which(diag[,"RID"] == mmseAIBLts@subj_id[i])
  
  ind2 <- which(diag[ind,"VISCODE"]=="bl")
  
  aibl_bl_diag[i] <- diag[ind[ind2],"DXCURREN"]
  
  inds <- which(diag[ind,"DXCURREN"]!= -4)
  
  ind3 <- order(diag[ind[inds],"VISCODE"],decreasing = T)[1]
  
  aibl_lv_diag[i] <- diag[ind[inds[ind3]],"DXCURREN"]
    
}

mmse_demog[c("Num CTL","Num MCI","Num AD"),"aibl"] <- table(aibl_bl_diag)

demog <- read.csv("aibl_ptdemog_28-Apr-2015.csv",header=T,as.is=T)

aibl_df <- data.frame(age=NA,gender=NA)

for (i in 1:length(mmseAIBLts@subj_id)){
  
  ind <- which(demog[,"RID"] == mmseAIBLts@subj_id[i])
  
  dob <- as.Date(paste(substr(demog[ind[1],"PTDOB"],2,5),"07",1,sep="-"))
  
  aibl_df[i,] <- c(signif(as.numeric(aibl_bl_date[[1]] - dob)/365.25,2),demog[ind[1],"PTGENDER"])
  
}

mmse_demog[c("Num Male","Num Female"),"aibl"] <- table(aibl_df[,"gender"])

mmse_demog["Median bl age","aibl"] <- median(aibl_df[,"age"])

mmse_demog["LQ bl age","aibl"] <- as.numeric(quantile(aibl_df[,"age"],0.25))

mmse_demog["UQ bl age","aibl"] <- as.numeric(quantile(aibl_df[,"age"],0.75))



apoe <- read.csv("aibl_apoeres_28-Apr-2015.csv",header=T,as.is=T)

aibl_num_apoe <- array(dim=c(length(mmseAIBLts@subj_id),1))

for (i in 1:length(mmseAIBLts@subj_id)){
  
  ind <- which(apoe[,"RID"] == mmseAIBLts@subj_id[i])
  
  aibl_num_apoe[i] <- (apoe[ind[1],"APGEN1"] == 4) + (apoe[ind[1],"APGEN2"] == 4)
  
}

mmse_demog[c("Num 0 E4","Num 1 E4","Num 2 E4"),"aibl"] <- table(aibl_num_apoe)

setwd("..")





# CAMD-1013  (Tombaugh and McIntyre (1992) for MMSE cutoffs, except no better than MCI due to recruitment of MCI and AD)

setwd("camd_june_2015/")

c1013_bl_mmse <- array(dim=c(length(mmse1013ts@subj_id),1))
c1013_lv_mmse <- array(dim=c(length(mmse1013ts@subj_id),1))

for (i in 1:length(mmse1013ts@subj_id)){
  
  c1013_bl_mmse[i] <- mmse1013ts@data[[i]][1]
  c1013_lv_mmse[i] <- mmse1013ts@data[[i]][length(mmse1013ts@data[[i]])]
  
}

c1013_diag <- (c1013_bl_mmse<18)+2
c1013_lv_diag <- (c1013_lv_mmse<18)+2

mmse_demog[c("Num MCI","Num AD"),"c1013"] <- table(c1013_diag)


demog <- read.csv("dm.csv",header=T,as.is=T)

c1013_df <- data.frame(age=NA,gender=NA)

for (i in 1:length(mmse1013ts@subj_id)){
  
  ind <- which(demog[,"USUBJID"] == mmse1013ts@subj_id[i])
  
  c1013_df[i,] <- c(demog[ind[1],"AGE"],as.numeric(demog[ind[1],"SEX"]=="F")+1)
  
}

mmse_demog[c("Num Male","Num Female"),"c1013"] <- table(c1013_df[,"gender"])

mmse_demog["Median bl age","c1013"] <- median(c1013_df[,"age"])

mmse_demog["LQ bl age","c1013"] <- as.numeric(quantile(c1013_df[,"age"],0.25))

mmse_demog["UQ bl age","c1013"] <- as.numeric(quantile(c1013_df[,"age"],0.75))


# CAMD-1014  (Tombaugh and McIntyre (1992) for MMSE cutoffs, except no better than MCI due to recruitment of MCI and AD)


c1014_bl_mmse <- array(dim=c(length(mmse1014ts@subj_id),1))
c1014_lv_mmse <- array(dim=c(length(mmse1014ts@subj_id),1))

for (i in 1:length(mmse1014ts@subj_id)){
  
  c1014_bl_mmse[i] <- mmse1014ts@data[[i]][1]
  c1014_lv_mmse[i] <- mmse1014ts@data[[i]][length(mmse1014ts@data[[i]])]
  
}

c1014_diag <- (c1014_bl_mmse<18)+2
c1014_lv_diag <- (c1014_lv_mmse<18)+2

mmse_demog[c("Num MCI","Num AD"),"c1014"] <- table(c1014_diag)


demog <- read.csv("dm.csv",header=T,as.is=T)

c1014_df <- data.frame(age=NA,gender=NA)

for (i in 1:length(mmse1014ts@subj_id)){
  
  ind <- which(demog[,"USUBJID"] == mmse1014ts@subj_id[i])
  
  c1014_df[i,] <- c(demog[ind[1],"AGE"],as.numeric(demog[ind[1],"SEX"]=="F")+1)
  
}


mmse_demog[c("Num Male","Num Female"),"c1014"] <- table(c1014_df[,"gender"])

mmse_demog["Median bl age","c1014"] <- median(c1014_df[,"age"])

mmse_demog["LQ bl age","c1014"] <- as.numeric(quantile(c1014_df[,"age"],0.25))

mmse_demog["UQ bl age","c1014"] <- as.numeric(quantile(c1014_df[,"age"],0.75))


mmse_demog[c("Num CTL","Num MCI","Num AD"),"combined"] <- c(sum(mmse_demog["Num CTL",2:5],na.rm=T),sum(mmse_demog["Num MCI",2:5],na.rm=T),sum(mmse_demog["Num AD",2:5],na.rm=T))

mmse_demog[c("Num Male","Num Female"),"combined"] <- c(sum(mmse_demog["Num Male",2:5],na.rm=T),sum(mmse_demog["Num Female",2:5],na.rm=T))

comb_age <- c(adni_df[,"age"],aibl_df[,"age"],c1013_df[,"age"],c1014_df[,"age"])

comb_gender <- c(adni_df[,"gender"],aibl_df[,"gender"],c1013_df[,"gender"],c1014_df[,"gender"])

mmse_demog["Median bl age","combined"] <- median(comb_age)

mmse_demog["LQ bl age","combined"] <- as.numeric(quantile(comb_age,0.25))

mmse_demog["UQ bl age","combined"] <- as.numeric(quantile(comb_age,0.75))


setwd("../..")

setwd("results/")

write.csv(mmse_demog[-1,],file="mmse_demog.csv")

bl_diags <- c(camd_diag,adni_bl_diag,aibl_bl_diag)
lv_diags <- c(camd_lv_diag,adni_lv_diag,aibl_lv_diag)


print("Table 1 (demographics) stored in results/")

