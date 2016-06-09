# Script to collect and summarise Disease Progression Score 
# simulation performed by TCsimDPS.R on a computing cluster

#install.packages('ggplot2')
library(ggplot2)

setwd("../results/")
setwd("DPSsim")

# as above but for AAC-style simulations


for (theta_num in 1:5){

  setwd("aac_sim")
  
  k <- 1
  
# Data frame to store data in format for plotting
sim_df <- data.frame(bias=NA,range=NA,theta1=NA,theta2=NA)




files <- list.files()

align_sim <- array(dim=c(length(files),5,3,3,2))

for (i in 1:length(files)){
  
  load(files[i])
  
  align_sim[i,,,,] <- sim_thetas
  
}



for (iter in 1:50){
  
  for (bias in 1:3){
    
    for (range in 1:3){
      
        sim_df[k,1:4] <- c(bias,range,align_sim[iter,theta_num,bias,range,])
        
        k <- k + 1 
      
      }
    
  }
  
}

sim_df[,"range_bias"] <- paste(sim_df[,"range"],sim_df[,"bias"])

setwd("../")

sim_df[,"range"] <- factor(sim_df[,"range"])
levels(sim_df[,"range"]) <- c("10000","20000","50000")

real_theta2 <- (seq(1,10,2)  * 10^-4)[theta_num] # theta_2 values used in the 5 simulations

# generate plots for Supplementary Table 1
file_name <- paste("align_",theta_num,"_1.pdf",sep="")

p <- ggplot(sim_df, aes(factor(range),theta1))
p + geom_boxplot(aes(fill = factor(bias))) + xlab("Range (+/-)") + ylab("Theta 1") + scale_fill_discrete(name="Bias")+ geom_hline(yintercept = 29,color="red",linetype="dashed",size=1)
ggsave(file_name,width=5,height=2.8)

# generate plots for Supplementary Table 2
file_name <- paste("align_",theta_num,"_2.pdf",sep="")


p <- ggplot(sim_df, aes(factor(range),theta2))
p + geom_boxplot(aes(fill = factor(bias))) + xlab("Range (+/-)") + ylab("Theta 2") + scale_fill_discrete(name="Bias")+ geom_hline(yintercept =  real_theta2,color="red",linetype="dashed",size=1)
ggsave(file_name,width=5,height=2.8)

}

print("Suppl Figure 1&2 plots created in results/DPSsim")
setwd("..")