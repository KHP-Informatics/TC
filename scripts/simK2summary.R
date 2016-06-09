thetas <- array(c(29,29,29,29,29,0.0005,0.0003,0.0005,0.0006,0.0008,29,29,29,29,29,0.0005,0.0001,0.0001,0.001,0.001),dim=c(5,4))


setwd("simK2scripts")

source("new_summary.R")
summary_df <- data.frame(theta_num=NA,k=NA,fix=NA,label_swap=NA,theta1=NA,theta2=NA)

# Load in the extra 30 iterations for normal
for (theta_num in 1:5){
  
  file_name <- paste("new_",theta_num,".RData",sep="")
  
  load(file_name)
  
  tmp_k1 <- cbind(theta_num,1,F,sim_df[,c(1:3)])
    
  tmp_k2 <- cbind(theta_num,2,F,sim_df[,c(1,4:5)])
  
  names(tmp_k2)[c(2:3,5:6)]<-c("k","fix","theta1","theta2")
  names(tmp_k1)[c(2:3,5:6)]<-c("k","fix","theta1","theta2")
  
  summary_df <- rbind(summary_df,tmp_k1,tmp_k2)
  
}

setwd("../../../scripts/simK2scripts/")

source("fix_summary.R")



# Load in the 50 fixed baseline iterations
for (theta_num in 1:5){
  
  file_name <- paste("fix_",theta_num,".RData",sep="")
  
  load(file_name)
  
  tmp_k1 <- cbind(theta_num,1,T,sim_df[,1:3])
    
  tmp_k2 <- cbind(theta_num,2,T,sim_df[,c(1,4:5)])
  
  names(tmp_k2)[c(2:3,5:6)]<-c("k","fix","theta1","theta2")
  names(tmp_k1)[c(2:3,5:6)]<-c("k","fix","theta1","theta2")
  
  summary_df <- rbind(summary_df,tmp_k1,tmp_k2)
  
}


summary_df[,"theta_num_k"] <- paste(summary_df[,1]," (",summary_df[,2],")",sep="")

summary_df <- summary_df[-1,]

setwd("../..")

p <- ggplot(summary_df, aes(factor(theta_num_k), theta2))

p + geom_boxplot(aes(fill = factor(fix))) +
  geom_segment(aes(x=0.5, xend=1.45, y=0.0005, yend=0.0005), linetype=1,col="red") +
  geom_segment(aes(x=1.55, xend=2.5, y=0.0005, yend=0.0005), linetype=1,col="red")+
  geom_segment(aes(x=2.5, xend=3.5, y=0.0003, yend=0.0003), linetype=1,col="red") +
  geom_segment(aes(x=3.5, xend=4.5, y=0.0001, yend=0.0001), linetype=1,col="red") + 
  geom_segment(aes(x=4.5, xend=5.5, y=0.0005, yend=0.0005), linetype=1,col="red") + 
  geom_segment(aes(x=5.5, xend=6.5, y=0.0001, yend=0.0001), linetype=1,col="red")+ 
  geom_segment(aes(x=6.5, xend=7.5, y=0.0006, yend=0.0006), linetype=1,col="red") + 
  geom_segment(aes(x=7.5, xend=8.5, y=0.001, yend=0.001), linetype=1,col="red")+ 
  geom_segment(aes(x=8.5, xend=9.5, y=0.0008, yend=0.0008), linetype=1,col="red") + 
  geom_segment(aes(x=9.5, xend=10.5, y=0.001, yend=0.001), linetype=1,col="red")+ 
  geom_segment(aes(x=2.5, xend=2.5, y=0.00015, yend=0.0012), linetype=2)+ 
  geom_segment(aes(x=4.5, xend=4.5, y=0.00015, yend=0.0012), linetype=2)+ 
  geom_segment(aes(x=6.5, xend=6.5, y=0.00015, yend=0.0012), linetype=2)+ 
  geom_segment(aes(x=8.5, xend=8.5, y=0.00015, yend=0.0012), linetype=2)+ 
  xlab("Simulations (Cluster)") + ylab("Theta 2") +
    scale_fill_discrete(name="Sim type & method",labels=c("standard", "fixed baseline"))

ggsave("simk2_supp_theta2.pdf",width=7,height=2.8)

p <- ggplot(summary_df, aes(factor(theta_num_k), theta1))

p + geom_boxplot(aes(fill = factor(fix))) +
  geom_segment(aes(x=0.5, xend=10.5, y=29, yend=29), linetype=1,col="red") +
  geom_segment(aes(x=2.5, xend=2.5, y=28, yend=45), linetype=2)+ 
  geom_segment(aes(x=4.5, xend=4.5, y=28, yend=45), linetype=2)+ 
  geom_segment(aes(x=6.5, xend=6.5, y=28, yend=45), linetype=2)+ 
  geom_segment(aes(x=8.5, xend=8.5, y=28, yend=45), linetype=2)+ 
  xlab("Simulations (Cluster)") + ylab("Theta 1") +
  scale_fill_discrete(name="Sim type & method",labels=c("separate baseline", "fixed baseline"))

ggsave("simk2_supp_theta1.pdf",width=7,height=2.8)




setwd("K2sim/fix/aac_sim/")

sim_names <- "sim"
tmp_file_name <- "aac_discrim.RData"

tmp_path <- "../../results/K2sim/fix/aac_sim/"

setwd("../../../../scripts/simK2scripts/")
source("K2discrimCollect.R")

install.packages('mclust')
library(mclust)

load("aac_discrim.RData")

thres <- c(0,2,5,10)

rand_df <- data.frame(thres=NA,sim_num=NA,adjRand=NA,num_left=NA)

k <- 1

for (i in 1:4){

  for (j in 1:200){
    
    pos_iter_ind <- which(pos_ctls[,1]==j)
    
    pos_thres_ind <- which(pos_ctls[pos_iter_ind,4] > thres[i])
    
    rand_df[k,1:2] <- c(thres[i],j) 
    
    rand_df[k,3] <- adjustedRandIndex(pos_ctls[pos_iter_ind[pos_thres_ind],2],pos_ctls[pos_iter_ind[pos_thres_ind],3])
      
    rand_df[k,4] <- length(pos_thres_ind)
    
    k <- k + 1
    
  } 
  
}

rand_df[,"thres"] <- factor(rand_df[,"thres"],thres)

pdf("new_discrim.pdf",2,4)
par(mfrow = c(2, 1))
par(cex=0.6)
par(mar=c(0.5,0.5,0.5,0.5),oma=c(4,4,0.5,0.5))
par(tcl=-0.25)
par(mgp=c(2,0.6,0))
boxplot(adjRand~thres,rand_df,ylim=c(0,1),xlab="Discrimination score threshold",ylab="Adjusted Rand Index",xaxt="n")

boxplot((num_left/24.12)~thres,rand_df,ylim=c(0,100),xlab="Discrimination score threshold",ylab="Number of individuals left")
mtext("Discrimination score threshold",side=1,outer=TRUE,cex=0.7,line=2.2)
mtext("Percentage individuals left              Clustering accuracy (ARI)",side=2,outer=TRUE,cex=0.7,line=2.2)
dev.off()




median(rand_df[rand_df[,"thres"]==0,"adjRand"],na.rm=T)
median(rand_df[rand_df[,"thres"]==2,"adjRand"],na.rm=T)
median(rand_df[rand_df[,"thres"]==5,"adjRand"],na.rm=T)
median(rand_df[rand_df[,"thres"]==10,"adjRand"],na.rm=T)

median(rand_df[rand_df[,"thres"]==0,"num_left"],na.rm=T)
median(rand_df[rand_df[,"thres"]==2,"num_left"],na.rm=T)
median(rand_df[rand_df[,"thres"]==5,"num_left"],na.rm=T)
median(rand_df[rand_df[,"thres"]==10,"num_left"],na.rm=T)
