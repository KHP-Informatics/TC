# script to convert multiTS to LCMM format

install.packages("NormPsy")
library(NormPsy)

install.packages("lcmm")
library(lcmm)

aac_lcmm <- data.frame(array(dim=c(length(which(!is.na(mmse_aac@mat_data))),4)))

colnames(aac_lcmm) <- c("ID","MMSE","TIME","AGE")

k <- 1

for (i in 1:mmse_aac@numTS){

  for (j in 1:length(which(!is.na(mmse_aac@mat_data[i,])))){
    
    aac_lcmm[k,] <- c(i,mmse_aac@mat_data[i,j],mmse_aac@mat_tps[i,j],aac_demog_999[i,"bl_age"]*365.25 + mmse_aac@mat_tps[i,j])
    
    k <- k + 1
    
  }
  
}

aac_lcmm[,"norm"] <- normMMSE(aac_lcmm[,2])

aac_lcmm[,"AGE50"] <- aac_lcmm[,"AGE"] - 50*365.25

lcmm_cl <- array(dim=c(2412,4))

mspl_raw_zero <- lcmm(MMSE ~ TIME,random=~ TIME,mixture=~ TIME, subject="ID",ng=2,data=aac_lcmm,link = "splines")
lcmm_cl[,1] <- mspl_raw_zero$pprob[,"class"]

mspl_norm_zero <- lcmm(norm ~ TIME,random=~ TIME,mixture=~ TIME, subject="ID",ng=2,data=aac_lcmm,link = "splines")
lcmm_cl[,2] <- mspl_norm_zero$pprob[,"class"]



aac_lcmm <- aac_lcmm[complete.cases(aac_lcmm),] # THIS EXCLUDES 13 CAMD SUBJECTS WITH ANONYMISED AGES, WHOSE AGES HAVE BEEN SET TO NA



comp.case <- unique(aac_lcmm[,1])

mspl_raw_age <- lcmm(MMSE ~ AGE50,random=~ AGE50,mixture=~ AGE50, subject="ID",ng=2,data=aac_lcmm,link = "splines")
lcmm_cl[comp.case,3] <- mspl_raw_age$pprob[,"class"]

mspl_norm_age <- lcmm(norm ~ AGE50,random=~ AGE50,mixture=~ AGE50, subject="ID",ng=2,data=aac_lcmm,link = "splines")
lcmm_cl[comp.case,4] <- mspl_norm_age$pprob[,"class"]

colnames(lcmm_cl) <- c("raw_zero","raw_age","norm_zero","norm_age")

# now compare to AD risk factors

tc_lcmm_df <- data.frame(cbind(lcmm_cl,aac_999_df))


for (i in 1:5){tc_lcmm_df[,i] <- factor(tc_lcmm_df[,i])}

TCvsLCMM <- data.frame(method=
    rep(c("Temporal Clustering","LCMM study time","LCMM age","normalised LCMM study time","normalised LCMM age"),each=4),
  variable = rep(c("MMSE first visit","1 APOE e4","2 APOE e4","CSF tau"),5),N=NA,lod=NA,se=NA,abs_z=NA,pvalue=NA)

# AAC MMSE at first visit
k_mmse_full <- summary(glm(k ~ bl_age + bl_mmse + gender + cohort, data = tc_lcmm_df, family = "binomial"))
raw_zero_mmse_full <- summary(glm(raw_zero ~ bl_age + bl_mmse + gender + cohort, data = tc_lcmm_df, family = "binomial"))
raw_age_mmse_full <- summary(glm(raw_age ~ bl_age + bl_mmse + gender + cohort, data = tc_lcmm_df, family = "binomial"))
norm_zero_mmse_full <- summary(glm(norm_zero ~ bl_age + bl_mmse + gender + cohort, data = tc_lcmm_df, family = "binomial"))
norm_age_mmse_full <- summary(glm(norm_age ~ bl_age + bl_mmse + gender + cohort, data = tc_lcmm_df, family = "binomial"))

TCvsLCMM[1,3:7] <- c(k_mmse_full$df.null,k_mmse_full$coef[3,])
TCvsLCMM[5,3:7] <- c(raw_zero_mmse_full$df.null,raw_zero_mmse_full$coef[3,])
TCvsLCMM[9,3:7] <- c(raw_age_mmse_full$df.null,raw_age_mmse_full$coef[3,])
TCvsLCMM[13,3:7] <- c(norm_zero_mmse_full$df.null,norm_zero_mmse_full$coef[3,])
TCvsLCMM[17,3:7] <- c(norm_age_mmse_full$df.null,norm_age_mmse_full$coef[3,])

# AAC APOE
k_apoe_full <- summary(glm(k ~ bl_age + bl_mmse + apoe + gender + cohort, data = tc_lcmm_df, family = "binomial"))
raw_zero_apoe_full <- summary(glm(raw_zero ~ bl_age + bl_mmse + apoe+  gender+cohort, data = tc_lcmm_df, family = "binomial"))
raw_age_apoe_full <- summary(glm(raw_age ~ bl_age + bl_mmse + apoe+ gender+cohort, data = tc_lcmm_df, family = "binomial"))
norm_zero_apoe_full <- summary(glm(norm_zero ~ bl_age + bl_mmse + apoe + gender+cohort, data = tc_lcmm_df, family = "binomial"))
norm_age_apoe_full <- summary(glm(norm_age ~ bl_age + bl_mmse + apoe  + gender+cohort, data = tc_lcmm_df, family = "binomial"))

TCvsLCMM[2:3,3:7] <- cbind(k_apoe_full$df.null,k_apoe_full$coef[4:5,])
TCvsLCMM[6:7,3:7] <- cbind(raw_zero_apoe_full$df.null,raw_zero_apoe_full$coef[4:5,])
TCvsLCMM[10:11,3:7] <- cbind(raw_age_apoe_full$df.null,raw_age_apoe_full$coef[4:5,])
TCvsLCMM[14:15,3:7] <- cbind(norm_zero_apoe_full$df.null,norm_zero_apoe_full$coef[4:5,])
TCvsLCMM[18:19,3:7] <- cbind(norm_age_apoe_full$df.null,norm_age_apoe_full$coef[4:5,])

# AAC tau
k_tau_full <- summary(glm(k ~ bl_age + bl_mmse + apoe + gender + tau, data = tc_lcmm_df, family = "binomial"))
raw_zero_tau_full <- summary(glm(raw_zero ~ bl_age + bl_mmse + apoe + gender + tau, data = tc_lcmm_df, family = "binomial"))
raw_age_tau_full <- summary(glm(raw_age ~ bl_age + bl_mmse + apoe + gender + tau, data = tc_lcmm_df, family = "binomial"))
norm_zero_tau_full <- summary(glm(norm_zero ~ bl_age + bl_mmse + apoe + gender + tau, data = tc_lcmm_df, family = "binomial"))
norm_age_tau_full <- summary(glm(norm_age ~ bl_age + bl_mmse + apoe  + gender + tau, data = tc_lcmm_df, family = "binomial"))

TCvsLCMM[4,3:7] <- c(k_tau_full$df.null,k_tau_full$coef[7,])
TCvsLCMM[8,3:7] <- c(raw_zero_tau_full$df.null,raw_zero_tau_full$coef[7,])
TCvsLCMM[12,3:7] <- c(raw_age_tau_full$df.null,raw_age_tau_full$coef[7,])
TCvsLCMM[16,3:7] <- c(norm_zero_tau_full$df.null,norm_zero_tau_full$coef[7,])
TCvsLCMM[20,3:7] <- c(norm_age_tau_full$df.null,norm_age_tau_full$coef[7,])



TCvsLCMM[,"N"] <- TCvsLCMM[,"N"] + 1

TCvsLCMM[,4:7] <- signif(TCvsLCMM[,4:7],2)

setwd("../results/")

write.csv(TCvsLCMM,file="TCvsLCMM_main.csv",row.names = F)





#plotting
#pred <- predictY(lin_raw_age,newdata=datnew,var.time="AGE50",draws=T)
#plot(pred$times[["AGE50"]],pred$pred[,"Ypred_50_class2"],ylim=c(0,100))
#points(pred$times[["TIME"]],pred$pred[,"Ypred_50_class1"],col="red")