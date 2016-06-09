setwd("../../results/K2sim/standard")

rb <- c("1_1","1_2","1_3","2_1","2_2","2_3","3_1","3_2","3_3")

thetas <- array(c(29,29,29,29,29,0.0005,0.0003,0.0005,0.0006,0.0008,29,29,29,29,29,0.0005,0.0001,0.0001,0.001,0.001),dim=c(5,4))


for (theta_num in 1:5){

sim_df <- 
data.frame(label_swap=logical(100),theta1_1=NA,theta1_2=NA,theta2_1=NA,theta2_2=NA)

k <- 1

setwd("aac_sim")

setwd(paste("theta",theta_num,sep=""))

files <- list.files()

for (i in 1:length(files)){
  
  load(files[i])

    
    
    #if (sum(table(results[[j]]@clusters)[array(c(1,2,1,2),dim=c(2,2))]) < sum(table(results[[j]]@clusters)[array(c(1,2,2,1),dim=c(2,2))])){
    if (dist(rbind(thetas[theta_num,c(2,4)],simRes@params[c(2,4)])) < dist(rbind(thetas[theta_num,c(2,4)],simRes@params[c(4,2)]))) { 
     
      sim_df[k,1] = T
      
      sim_df[k,2:5] <- simRes@params
      
    } else {
    
      sim_df[k,c(4,5,2,3)] <- simRes@params
        
    }
    
    
    k <- k + 1
     
  
  
}


setwd("..")
setwd("..")

file_name <- paste("new_",theta_num,".RData",sep="")

save(sim_df,file=file_name)

}

