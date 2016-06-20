setwd(tmp_path)

pos_ctls <- data.frame(iter=NA,real_clusters=NA,est_clusters=NA,discrim_score=NA,num_tp=NA,length=NA,real_shift=NA,mmse_at_first_visit=NA)

for (sim_num in 2:5){
  
  print(sim_num)
  
  setwd(paste("theta",sim_num,sep=""))
  
  for (iter in 1:100){
  
    file_name <- paste(sim_names,iter,".RData",sep="")
  
    load(file_name)

    tmp_table <- data.frame(iter=iter + (sim_num-2)*100,real_clusters=simRes@clusters[,"real_clusters"],est_clusters=simRes@clusters[,"est_clusters"],discrim_score=simRes@loss[,3],num_tp=NA,length=NA,real_shift=simRes@shifts[,1],mmse_at_first_visit=simRes@sim_data[,1])
  
    sim_length <- c()

    sim_num_tp <- c()

    for (i in 1:2412){
  
      sim_num_tp[i] <- length(which(!is.na(simRes@sim_tps[i,])))
  
      sim_length[i] <- simRes@sim_tps[i,sim_num_tp[i]]
    
    }

    tmp_table[,"num_tp"] <- sim_num_tp
    tmp_table[,"length"] <- sim_length
    
    pos_ctls <- rbind(pos_ctls,tmp_table)
  
  }
  
  setwd("..")
  
}

setwd("..")
setwd("..")
setwd("..")

save(pos_ctls,file=tmp_file_name)