nruns<-100 
methods<-3

censoringvalues_all<-c(0.1,0.2,0.3,0.5)
n1_all<-c(200,300,500,1000,2000)
q_all<-c(10)
q_true_all<-c(2)
rho_all<-c(0.3)


for (cen_i in 1:length(censoringvalues_all)) {
  for (n_j in 1:length(n1_all)) {
    for (q_i in 1:length(q_all)) {
      for (q_true_i in 1:length(q_true_all)) {
        for (rho_k in 1:length(rho_all)) { 
          n1<-n1_all[n_j]
          mfunctions<-q_all[q_i]
          mfunctions_true<-q_true_all[q_true_i]
          q<-q_all[q_i]
          q_true<-q_true_all[q_true_i]
          censoringvalues<-censoringvalues_all[cen_i]
          rho<-rho_all[rho_k]
          
          tp_scalar<-matrix(0,nruns,methods)
          fp_scalar<-matrix(0,nruns,methods)
          tn_scalar<-matrix(0,nruns,methods)
          fn_scalar<-matrix(0,nruns,methods)
          
          tp_dev<-matrix(0,nruns,methods)
          fp_dev<-matrix(0,nruns,methods)
          tn_dev<-matrix(0,nruns,methods)
          fn_dev<-matrix(0,nruns,methods)
          
          tp_fun<-matrix(0,nruns,methods)
          fp_fun<-matrix(0,nruns,methods)
          tn_fun<-matrix(0,nruns,methods)
          fn_fun<-matrix(0,nruns,methods)
          
          censoring<-matrix(0,nruns,1)
          for(runs_i in 1:nruns){
            filename=paste(runs_i,"_rho",rho,"_censoringvalues",censoringvalues,"_n",n1,"_mfunctions",mfunctions,"_mfunctions_true",mfunctions_true,"_q",q,"_q_true",q_true,".csv",sep="")
            
            result<-read.csv(filename)
            censoring[runs_i]<-result[4,2]
            result<-result[1:3,-1]
            tp_scalar[runs_i,]<-result[,1]
            fp_scalar[runs_i,]<-result[,2]
            tn_scalar[runs_i,]<-result[,3]
            fn_scalar[runs_i,]<-result[,4]
            tp_dev[runs_i,]<-result[,5]
            fp_dev[runs_i,]<-result[,6]
            tn_dev[runs_i,]<-result[,7]
            fn_dev[runs_i,]<-result[,8]
            tp_fun[runs_i,]<-result[,9] 
            fp_fun[runs_i,]<-result[,10]
            tn_fun[runs_i,]<-result[,11]
            fn_fun[runs_i,]<-result[,12]
          }
          tpr_scalar<-tp_scalar/(tp_scalar+fn_scalar)
          fpr_scalar<-fp_scalar/(fp_scalar+tn_scalar)
          tpr_fun<-tp_fun/(tp_fun+fn_fun)
          fpr_fun<-fp_fun/(fp_fun+tn_fun)
          tpr_dev<-tp_dev/(tp_dev+fn_dev)
          fpr_dev<-fp_dev/(fp_dev+tn_dev)
          
          
          
          fdr_scalar<-fp_scalar/(tp_scalar+fp_scalar)
          fnr_scalar<-fn_scalar/(tp_scalar+fn_scalar)
          sensitivity_scalar<-1-fnr_scalar
          specificity_scalar<-1-fp_scalar/(fp_scalar+tn_scalar)
          
          fdr_dev<-fp_dev/(tp_dev+fp_dev)
          fnr_dev<-fn_dev/(tp_dev+fn_dev)
          sensitivity_dev<-1-fnr_dev
          specificity_dev<-1-fp_dev/(fp_dev+tn_dev)
          
          
          fdr_fun<-fp_fun/(tp_fun+fp_fun)
          fnr_fun<-fn_fun/(tp_fun+fn_fun)
          sensitivity_fun<-1-fnr_fun
          specificity_fun<-1-fp_fun/(fp_fun+tn_fun)
          
          fdr_all<-(fp_fun+fp_scalar+fp_dev)/(tp_fun+fp_fun+tp_dev+fp_dev+tp_scalar+fp_scalar)
          fnr_all<-(fn_fun+fn_scalar+fn_dev)/(tp_fun+fn_fun+tp_dev+fn_dev+tp_scalar+fn_scalar)
          sensitivity_all<-1-fnr_all
          specificity_all<-1-(fp_fun+fp_scalar+fp_dev)/(fp_fun+tn_fun+fp_dev+tn_dev+fp_scalar+tn_scalar)
          
          
          fp_all<-fp_fun+fp_scalar+fp_dev
          tp_all<-tp_fun +tp_scalar+tp_dev
          
          
          fn_all<-fn_fun+fn_scalar+fn_dev
          tn_all<-tn_fun +tn_scalar +tn_dev
          
          
          apply(fdr_all,2,mean)
          apply(fnr_all,2,mean)
          apply(sensitivity_all,2,mean)
          apply(specificity_all,2,mean)
          
          apply(fp_all,2,mean)
          apply(tp_all,2,mean)
          
          apply(fn_all,2,mean)
          apply(tn_all,2,mean)
          
          apply(tpr_scalar,2,mean)
          apply(fpr_scalar,2,mean)
          apply(tpr_fun,2,mean)
          apply(fpr_fun,2,mean)
          apply(tpr_dev,2,mean)
          apply(fpr_dev,2,mean)
          
          summary(censoring)
          plot(censoring)
          
          
          save_file<-paste("q",q,"true",q_true,"n",n1,"c",censoringvalues,"rho",rho,".RData",sep="")
          
          save.image(save_file)
        }
      }
      
    }
    
  }
}

censoringvalues_all<-c(0.1,0.2,0.3,0.5)
n1_all<-c(200,300,500,1000,2000)
q_all<-c(10)
q_true_all<-c(2)
rho_all<-c(0.3)


for (cen_i in 1:length(censoringvalues_all)) {
  for (n_j in 1:length(n1_all)) {
    for (q_i in 1:length(q_all)) {
      for (q_true_i in 1:length(q_true_all)) {
        for (rho_k in 1:length(rho_all)) { 
          n1<-n1_all[n_j]
          mfunctions<-q_all[q_i]
          mfunctions_true<-q_true_all[q_true_i]
          q<-q_all[q_i]
          q_true<-q_true_all[q_true_i]
          censoringvalues<-censoringvalues_all[cen_i]
          rho<-rho_all[rho_k]
          
          tp_scalar<-matrix(0,nruns,methods)
          fp_scalar<-matrix(0,nruns,methods)
          tn_scalar<-matrix(0,nruns,methods)
          fn_scalar<-matrix(0,nruns,methods)
          
          tp_dev<-matrix(0,nruns,methods)
          fp_dev<-matrix(0,nruns,methods)
          tn_dev<-matrix(0,nruns,methods)
          fn_dev<-matrix(0,nruns,methods)
          
          tp_fun<-matrix(0,nruns,methods)
          fp_fun<-matrix(0,nruns,methods)
          tn_fun<-matrix(0,nruns,methods)
          fn_fun<-matrix(0,nruns,methods)
          
          censoring<-matrix(0,nruns,1)
          for(runs_i in 1:nruns){
            filename=paste(runs_i,"_rho",rho,"_censoringvalues",censoringvalues,"_n",n1,"_mfunctions",mfunctions,"_mfunctions_true",mfunctions_true,"_q",q,"_q_true",q_true,".csv",sep="")
            
            result<-read.csv(filename)
            censoring[runs_i]<-result[4,2]
            result<-result[1:3,-1]
            tp_scalar[runs_i,]<-result[,1]
            fp_scalar[runs_i,]<-result[,2]
            tn_scalar[runs_i,]<-result[,3]
            fn_scalar[runs_i,]<-result[,4]
            tp_dev[runs_i,]<-result[,5]
            fp_dev[runs_i,]<-result[,6]
            tn_dev[runs_i,]<-result[,7]
            fn_dev[runs_i,]<-result[,8]
            tp_fun[runs_i,]<-result[,9] 
            fp_fun[runs_i,]<-result[,10]
            tn_fun[runs_i,]<-result[,11]
            fn_fun[runs_i,]<-result[,12]
          }
          tpr_scalar<-tp_scalar/(tp_scalar+fn_scalar)
          fpr_scalar<-fp_scalar/(fp_scalar+tn_scalar)
          tpr_fun<-tp_fun/(tp_fun+fn_fun)
          fpr_fun<-fp_fun/(fp_fun+tn_fun)
          tpr_dev<-tp_dev/(tp_dev+fn_dev)
          fpr_dev<-fp_dev/(fp_dev+tn_dev)
          
          
          
          fdr_scalar<-fp_scalar/(tp_scalar+fp_scalar)
          fnr_scalar<-fn_scalar/(tp_scalar+fn_scalar)
          sensitivity_scalar<-1-fnr_scalar
          specificity_scalar<-1-fp_scalar/(fp_scalar+tn_scalar)
          
          fdr_dev<-fp_dev/(tp_dev+fp_dev)
          fnr_dev<-fn_dev/(tp_dev+fn_dev)
          sensitivity_dev<-1-fnr_dev
          specificity_dev<-1-fp_dev/(fp_dev+tn_dev)
          
          
          fdr_fun<-fp_fun/(tp_fun+fp_fun)
          fnr_fun<-fn_fun/(tp_fun+fn_fun)
          sensitivity_fun<-1-fnr_fun
          specificity_fun<-1-fp_fun/(fp_fun+tn_fun)
          
          fdr_all<-(fp_fun+fp_scalar+fp_dev)/(tp_fun+fp_fun+tp_dev+fp_dev+tp_scalar+fp_scalar)
          fnr_all<-(fn_fun+fn_scalar+fn_dev)/(tp_fun+fn_fun+tp_dev+fn_dev+tp_scalar+fn_scalar)
          sensitivity_all<-1-fnr_all
          specificity_all<-1-(fp_fun+fp_scalar+fp_dev)/(fp_fun+tn_fun+fp_dev+tn_dev+fp_scalar+tn_scalar)
          
          
          fp_all<-fp_fun+fp_scalar+fp_dev
          tp_all<-tp_fun +tp_scalar+tp_dev
          
          
          fn_all<-fn_fun+fn_scalar+fn_dev
          tn_all<-tn_fun +tn_scalar +tn_dev
          
          
          apply(fdr_all,2,mean)
          apply(fnr_all,2,mean)
          apply(sensitivity_all,2,mean)
          apply(specificity_all,2,mean)
          
          apply(fp_all,2,mean)
          apply(tp_all,2,mean)
          
          apply(fn_all,2,mean)
          apply(tn_all,2,mean)
          
          apply(tpr_scalar,2,mean)
          apply(fpr_scalar,2,mean)
          apply(tpr_fun,2,mean)
          apply(fpr_fun,2,mean)
          apply(tpr_dev,2,mean)
          apply(fpr_dev,2,mean)
          
          summary(censoring)
          plot(censoring)
          
          
          save_file<-paste("q",q,"true",q_true,"n",n1,"c",censoringvalues,"rho",rho,".RData",sep="")
          
          save.image(save_file)
        }
      }
      
    }
    
  }
}

censoringvalues_all<-c(0.1,0.2,0.3,0.5)
n1_all<-c(200,300,500,1000,2000)
q_all<-c(10)
q_true_all<-c(2)
rho_all<-c(0.8)


for (cen_i in 1:length(censoringvalues_all)) {
  for (n_j in 1:length(n1_all)) {
    for (q_i in 1:length(q_all)) {
      for (q_true_i in 1:length(q_true_all)) {
        for (rho_k in 1:length(rho_all)) { 
          n1<-n1_all[n_j]
          mfunctions<-q_all[q_i]
          mfunctions_true<-q_true_all[q_true_i]
          q<-q_all[q_i]
          q_true<-q_true_all[q_true_i]
          censoringvalues<-censoringvalues_all[cen_i]
          rho<-rho_all[rho_k]
          
          tp_scalar<-matrix(0,nruns,methods)
          fp_scalar<-matrix(0,nruns,methods)
          tn_scalar<-matrix(0,nruns,methods)
          fn_scalar<-matrix(0,nruns,methods)
          
          tp_dev<-matrix(0,nruns,methods)
          fp_dev<-matrix(0,nruns,methods)
          tn_dev<-matrix(0,nruns,methods)
          fn_dev<-matrix(0,nruns,methods)
          
          tp_fun<-matrix(0,nruns,methods)
          fp_fun<-matrix(0,nruns,methods)
          tn_fun<-matrix(0,nruns,methods)
          fn_fun<-matrix(0,nruns,methods)
          
          censoring<-matrix(0,nruns,1)
          for(runs_i in 1:nruns){
            filename=paste(runs_i,"_rho",rho,"_censoringvalues",censoringvalues,"_n",n1,"_mfunctions",mfunctions,"_mfunctions_true",mfunctions_true,"_q",q,"_q_true",q_true,".csv",sep="")
            
            result<-read.csv(filename)
            censoring[runs_i]<-result[4,2]
            result<-result[1:3,-1]
            tp_scalar[runs_i,]<-result[,1]
            fp_scalar[runs_i,]<-result[,2]
            tn_scalar[runs_i,]<-result[,3]
            fn_scalar[runs_i,]<-result[,4]
            tp_dev[runs_i,]<-result[,5]
            fp_dev[runs_i,]<-result[,6]
            tn_dev[runs_i,]<-result[,7]
            fn_dev[runs_i,]<-result[,8]
            tp_fun[runs_i,]<-result[,9] 
            fp_fun[runs_i,]<-result[,10]
            tn_fun[runs_i,]<-result[,11]
            fn_fun[runs_i,]<-result[,12]
          }
          tpr_scalar<-tp_scalar/(tp_scalar+fn_scalar)
          fpr_scalar<-fp_scalar/(fp_scalar+tn_scalar)
          tpr_fun<-tp_fun/(tp_fun+fn_fun)
          fpr_fun<-fp_fun/(fp_fun+tn_fun)
          tpr_dev<-tp_dev/(tp_dev+fn_dev)
          fpr_dev<-fp_dev/(fp_dev+tn_dev)
          
          
          
          fdr_scalar<-fp_scalar/(tp_scalar+fp_scalar)
          fnr_scalar<-fn_scalar/(tp_scalar+fn_scalar)
          sensitivity_scalar<-1-fnr_scalar
          specificity_scalar<-1-fp_scalar/(fp_scalar+tn_scalar)
          
          fdr_dev<-fp_dev/(tp_dev+fp_dev)
          fnr_dev<-fn_dev/(tp_dev+fn_dev)
          sensitivity_dev<-1-fnr_dev
          specificity_dev<-1-fp_dev/(fp_dev+tn_dev)
          
          
          fdr_fun<-fp_fun/(tp_fun+fp_fun)
          fnr_fun<-fn_fun/(tp_fun+fn_fun)
          sensitivity_fun<-1-fnr_fun
          specificity_fun<-1-fp_fun/(fp_fun+tn_fun)
          
          fdr_all<-(fp_fun+fp_scalar+fp_dev)/(tp_fun+fp_fun+tp_dev+fp_dev+tp_scalar+fp_scalar)
          fnr_all<-(fn_fun+fn_scalar+fn_dev)/(tp_fun+fn_fun+tp_dev+fn_dev+tp_scalar+fn_scalar)
          sensitivity_all<-1-fnr_all
          specificity_all<-1-(fp_fun+fp_scalar+fp_dev)/(fp_fun+tn_fun+fp_dev+tn_dev+fp_scalar+tn_scalar)
          
          
          fp_all<-fp_fun+fp_scalar+fp_dev
          tp_all<-tp_fun +tp_scalar+tp_dev
          
          
          fn_all<-fn_fun+fn_scalar+fn_dev
          tn_all<-tn_fun +tn_scalar +tn_dev
          
          
          apply(fdr_all,2,mean)
          apply(fnr_all,2,mean)
          apply(sensitivity_all,2,mean)
          apply(specificity_all,2,mean)
          
          apply(fp_all,2,mean)
          apply(tp_all,2,mean)
          
          apply(fn_all,2,mean)
          apply(tn_all,2,mean)
          
          apply(tpr_scalar,2,mean)
          apply(fpr_scalar,2,mean)
          apply(tpr_fun,2,mean)
          apply(fpr_fun,2,mean)
          apply(tpr_dev,2,mean)
          apply(fpr_dev,2,mean)
          
          summary(censoring)
          plot(censoring)
          
          
          save_file<-paste("q",q,"true",q_true,"n",n1,"c",censoringvalues,"rho",rho,".RData",sep="")
          
          save.image(save_file)
        }
      }
      
    }
    
  }
}
censoringvalues_all<-c(0.1,0.2,0.3,0.5)
n1_all<-c(200,300,500,1000,2000)
q_all<-c(30)
q_true_all<-c(2)
rho_all<-c(0.5)


for (cen_i in 1:length(censoringvalues_all)) {
  for (n_j in 1:length(n1_all)) {
    for (q_i in 1:length(q_all)) {
      for (q_true_i in 1:length(q_true_all)) {
        for (rho_k in 1:length(rho_all)) { 
          n1<-n1_all[n_j]
          mfunctions<-q_all[q_i]
          mfunctions_true<-q_true_all[q_true_i]
          q<-q_all[q_i]
          q_true<-q_true_all[q_true_i]
          censoringvalues<-censoringvalues_all[cen_i]
          rho<-rho_all[rho_k]
          
          tp_scalar<-matrix(0,nruns,methods)
          fp_scalar<-matrix(0,nruns,methods)
          tn_scalar<-matrix(0,nruns,methods)
          fn_scalar<-matrix(0,nruns,methods)
          
          tp_dev<-matrix(0,nruns,methods)
          fp_dev<-matrix(0,nruns,methods)
          tn_dev<-matrix(0,nruns,methods)
          fn_dev<-matrix(0,nruns,methods)
          
          tp_fun<-matrix(0,nruns,methods)
          fp_fun<-matrix(0,nruns,methods)
          tn_fun<-matrix(0,nruns,methods)
          fn_fun<-matrix(0,nruns,methods)
          
          censoring<-matrix(0,nruns,1)
          for(runs_i in 1:nruns){
            filename=paste(runs_i,"_rho",rho,"_censoringvalues",censoringvalues,"_n",n1,"_mfunctions",mfunctions,"_mfunctions_true",mfunctions_true,"_q",q,"_q_true",q_true,".csv",sep="")
            
            result<-read.csv(filename)
            censoring[runs_i]<-result[4,2]
            result<-result[1:3,-1]
            tp_scalar[runs_i,]<-result[,1]
            fp_scalar[runs_i,]<-result[,2]
            tn_scalar[runs_i,]<-result[,3]
            fn_scalar[runs_i,]<-result[,4]
            tp_dev[runs_i,]<-result[,5]
            fp_dev[runs_i,]<-result[,6]
            tn_dev[runs_i,]<-result[,7]
            fn_dev[runs_i,]<-result[,8]
            tp_fun[runs_i,]<-result[,9] 
            fp_fun[runs_i,]<-result[,10]
            tn_fun[runs_i,]<-result[,11]
            fn_fun[runs_i,]<-result[,12]
          }
          tpr_scalar<-tp_scalar/(tp_scalar+fn_scalar)
          fpr_scalar<-fp_scalar/(fp_scalar+tn_scalar)
          tpr_fun<-tp_fun/(tp_fun+fn_fun)
          fpr_fun<-fp_fun/(fp_fun+tn_fun)
          tpr_dev<-tp_dev/(tp_dev+fn_dev)
          fpr_dev<-fp_dev/(fp_dev+tn_dev)
          
          
          
          fdr_scalar<-fp_scalar/(tp_scalar+fp_scalar)
          fnr_scalar<-fn_scalar/(tp_scalar+fn_scalar)
          sensitivity_scalar<-1-fnr_scalar
          specificity_scalar<-1-fp_scalar/(fp_scalar+tn_scalar)
          
          fdr_dev<-fp_dev/(tp_dev+fp_dev)
          fnr_dev<-fn_dev/(tp_dev+fn_dev)
          sensitivity_dev<-1-fnr_dev
          specificity_dev<-1-fp_dev/(fp_dev+tn_dev)
          
          
          fdr_fun<-fp_fun/(tp_fun+fp_fun)
          fnr_fun<-fn_fun/(tp_fun+fn_fun)
          sensitivity_fun<-1-fnr_fun
          specificity_fun<-1-fp_fun/(fp_fun+tn_fun)
          
          fdr_all<-(fp_fun+fp_scalar+fp_dev)/(tp_fun+fp_fun+tp_dev+fp_dev+tp_scalar+fp_scalar)
          fnr_all<-(fn_fun+fn_scalar+fn_dev)/(tp_fun+fn_fun+tp_dev+fn_dev+tp_scalar+fn_scalar)
          sensitivity_all<-1-fnr_all
          specificity_all<-1-(fp_fun+fp_scalar+fp_dev)/(fp_fun+tn_fun+fp_dev+tn_dev+fp_scalar+tn_scalar)
          
          
          fp_all<-fp_fun+fp_scalar+fp_dev
          tp_all<-tp_fun +tp_scalar+tp_dev
          
          
          fn_all<-fn_fun+fn_scalar+fn_dev
          tn_all<-tn_fun +tn_scalar +tn_dev
          
          
          apply(fdr_all,2,mean)
          apply(fnr_all,2,mean)
          apply(sensitivity_all,2,mean)
          apply(specificity_all,2,mean)
          
          apply(fp_all,2,mean)
          apply(tp_all,2,mean)
          
          apply(fn_all,2,mean)
          apply(tn_all,2,mean)
          
          apply(tpr_scalar,2,mean)
          apply(fpr_scalar,2,mean)
          apply(tpr_fun,2,mean)
          apply(fpr_fun,2,mean)
          apply(tpr_dev,2,mean)
          apply(fpr_dev,2,mean)
          
          summary(censoring)
          plot(censoring)
          
          
          save_file<-paste("q",q,"true",q_true,"n",n1,"c",censoringvalues,"rho",rho,".RData",sep="")
          
          save.image(save_file)
        }
      }
      
    }
    
  }
}

censoringvalues_all<-c(0.1,0.2,0.3,0.5)
n1_all<-c(200,300,500,1000,2000)
q_all<-c(50)
q_true_all<-c(2)
rho_all<-c(0.5)


for (cen_i in 1:length(censoringvalues_all)) {
  for (n_j in 1:length(n1_all)) {
    for (q_i in 1:length(q_all)) {
      for (q_true_i in 1:length(q_true_all)) {
        for (rho_k in 1:length(rho_all)) { 
          n1<-n1_all[n_j]
          mfunctions<-q_all[q_i]
          mfunctions_true<-q_true_all[q_true_i]
          q<-q_all[q_i]
          q_true<-q_true_all[q_true_i]
          censoringvalues<-censoringvalues_all[cen_i]
          rho<-rho_all[rho_k]
          
          tp_scalar<-matrix(0,nruns,methods)
          fp_scalar<-matrix(0,nruns,methods)
          tn_scalar<-matrix(0,nruns,methods)
          fn_scalar<-matrix(0,nruns,methods)
          
          tp_dev<-matrix(0,nruns,methods)
          fp_dev<-matrix(0,nruns,methods)
          tn_dev<-matrix(0,nruns,methods)
          fn_dev<-matrix(0,nruns,methods)
          
          tp_fun<-matrix(0,nruns,methods)
          fp_fun<-matrix(0,nruns,methods)
          tn_fun<-matrix(0,nruns,methods)
          fn_fun<-matrix(0,nruns,methods)
          
          censoring<-matrix(0,nruns,1)
          for(runs_i in 1:nruns){
            filename=paste(runs_i,"_rho",rho,"_censoringvalues",censoringvalues,"_n",n1,"_mfunctions",mfunctions,"_mfunctions_true",mfunctions_true,"_q",q,"_q_true",q_true,".csv",sep="")
            
            result<-read.csv(filename)
            censoring[runs_i]<-result[4,2]
            result<-result[1:3,-1]
            tp_scalar[runs_i,]<-result[,1]
            fp_scalar[runs_i,]<-result[,2]
            tn_scalar[runs_i,]<-result[,3]
            fn_scalar[runs_i,]<-result[,4]
            tp_dev[runs_i,]<-result[,5]
            fp_dev[runs_i,]<-result[,6]
            tn_dev[runs_i,]<-result[,7]
            fn_dev[runs_i,]<-result[,8]
            tp_fun[runs_i,]<-result[,9] 
            fp_fun[runs_i,]<-result[,10]
            tn_fun[runs_i,]<-result[,11]
            fn_fun[runs_i,]<-result[,12]
          }
          tpr_scalar<-tp_scalar/(tp_scalar+fn_scalar)
          fpr_scalar<-fp_scalar/(fp_scalar+tn_scalar)
          tpr_fun<-tp_fun/(tp_fun+fn_fun)
          fpr_fun<-fp_fun/(fp_fun+tn_fun)
          tpr_dev<-tp_dev/(tp_dev+fn_dev)
          fpr_dev<-fp_dev/(fp_dev+tn_dev)
          
          
          
          fdr_scalar<-fp_scalar/(tp_scalar+fp_scalar)
          fnr_scalar<-fn_scalar/(tp_scalar+fn_scalar)
          sensitivity_scalar<-1-fnr_scalar
          specificity_scalar<-1-fp_scalar/(fp_scalar+tn_scalar)
          
          fdr_dev<-fp_dev/(tp_dev+fp_dev)
          fnr_dev<-fn_dev/(tp_dev+fn_dev)
          sensitivity_dev<-1-fnr_dev
          specificity_dev<-1-fp_dev/(fp_dev+tn_dev)
          
          
          fdr_fun<-fp_fun/(tp_fun+fp_fun)
          fnr_fun<-fn_fun/(tp_fun+fn_fun)
          sensitivity_fun<-1-fnr_fun
          specificity_fun<-1-fp_fun/(fp_fun+tn_fun)
          
          fdr_all<-(fp_fun+fp_scalar+fp_dev)/(tp_fun+fp_fun+tp_dev+fp_dev+tp_scalar+fp_scalar)
          fnr_all<-(fn_fun+fn_scalar+fn_dev)/(tp_fun+fn_fun+tp_dev+fn_dev+tp_scalar+fn_scalar)
          sensitivity_all<-1-fnr_all
          specificity_all<-1-(fp_fun+fp_scalar+fp_dev)/(fp_fun+tn_fun+fp_dev+tn_dev+fp_scalar+tn_scalar)
          
          
          fp_all<-fp_fun+fp_scalar+fp_dev
          tp_all<-tp_fun +tp_scalar+tp_dev
          
          
          fn_all<-fn_fun+fn_scalar+fn_dev
          tn_all<-tn_fun +tn_scalar +tn_dev
          
          
          apply(fdr_all,2,mean)
          apply(fnr_all,2,mean)
          apply(sensitivity_all,2,mean)
          apply(specificity_all,2,mean)
          
          apply(fp_all,2,mean)
          apply(tp_all,2,mean)
          
          apply(fn_all,2,mean)
          apply(tn_all,2,mean)
          
          apply(tpr_scalar,2,mean)
          apply(fpr_scalar,2,mean)
          apply(tpr_fun,2,mean)
          apply(fpr_fun,2,mean)
          apply(tpr_dev,2,mean)
          apply(fpr_dev,2,mean)
          
          summary(censoring)
          plot(censoring)
          
          
          save_file<-paste("q",q,"true",q_true,"n",n1,"c",censoringvalues,"rho",rho,".RData",sep="")
          
          save.image(save_file)
        }
      }
      
    }
    
  }
}
censoringvalues_all<-c(0.1,0.2,0.3,0.5)
rho_all<-c(0.5 )
n1_all<-c(200,300,500,1000,2000)
nruns<-100  
methods<-6
q<-50
q_true<-2
mfunctions<-q
mfunctions_true<-q_true

for (cen_i in 1:length(censoringvalues_all)) {
  for (n_j in 1:length(n1_all)) {
    for (rho_k in 1:length(rho_all)) {
      
      n1<-n1_all[n_j]
      
      censoringvalues<-censoringvalues_all[cen_i]
      rho<-rho_all[rho_k]
      tp_scalar<-matrix(0,nruns,methods)
      fp_scalar<-matrix(0,nruns,methods)
      tn_scalar<-matrix(0,nruns,methods)
      fn_scalar<-matrix(0,nruns,methods)
      
      tp_dev<-matrix(0,nruns,methods)
      fp_dev<-matrix(0,nruns,methods)
      tn_dev<-matrix(0,nruns,methods)
      fn_dev<-matrix(0,nruns,methods)
      
      tp_fun<-matrix(0,nruns,methods)
      fp_fun<-matrix(0,nruns,methods)
      tn_fun<-matrix(0,nruns,methods)
      fn_fun<-matrix(0,nruns,methods)
      
      censoring<-matrix(0,nruns,1)
      for(runs_i in 1:nruns){
        filename=paste(runs_i,"_rho",rho,"_censoringvalues",censoringvalues,"_n",n1,"_mfunctions",mfunctions,"_mfunctions_true",mfunctions_true,"_q",q,"_q_true",q_true,".csv",sep="")
        
        result<-read.csv(filename)
        censoring[runs_i]<-result[4,2]
        result<-result[-4,-1]
        tp_scalar[runs_i,]<-result[,1]
        fp_scalar[runs_i,]<-result[,2]
        tn_scalar[runs_i,]<-result[,3]
        fn_scalar[runs_i,]<-result[,4]
        tp_dev[runs_i,]<-result[,5]
        fp_dev[runs_i,]<-result[,6]
        tn_dev[runs_i,]<-result[,7]
        fn_dev[runs_i,]<-result[,8]
        tp_fun[runs_i,]<-result[,9] 
        fp_fun[runs_i,]<-result[,10]
        tn_fun[runs_i,]<-result[,11]
        fn_fun[runs_i,]<-result[,12]
      }
      tpr_scalar<-tp_scalar/(tp_scalar+fn_scalar)
      fpr_scalar<-fp_scalar/(fp_scalar+tn_scalar)
      tpr_fun<-tp_fun/(tp_fun+fn_fun)
      fpr_fun<-fp_fun/(fp_fun+tn_fun)
      tpr_dev<-tp_dev/(tp_dev+fn_dev)
      fpr_dev<-fp_dev/(fp_dev+tn_dev)
      
      
      
      fdr_scalar<-fp_scalar/(tp_scalar+fp_scalar)
      fnr_scalar<-fn_scalar/(tp_scalar+fn_scalar)
      sensitivity_scalar<-1-fnr_scalar
      specificity_scalar<-1-fp_scalar/(fp_scalar+tn_scalar)
      
      fdr_dev<-fp_dev/(tp_dev+fp_dev)
      fnr_dev<-fn_dev/(tp_dev+fn_dev)
      sensitivity_dev<-1-fnr_dev
      specificity_dev<-1-fp_dev/(fp_dev+tn_dev)
      
      
      fdr_fun<-fp_fun/(tp_fun+fp_fun)
      fnr_fun<-fn_fun/(tp_fun+fn_fun)
      sensitivity_fun<-1-fnr_fun
      specificity_fun<-1-fp_fun/(fp_fun+tn_fun)
      
      fdr_all<-(fp_fun+fp_scalar+fp_dev)/(tp_fun+fp_fun+tp_dev+fp_dev+tp_scalar+fp_scalar)
      fnr_all<-(fn_fun+fn_scalar+fn_dev)/(tp_fun+fn_fun+tp_dev+fn_dev+tp_scalar+fn_scalar)
      sensitivity_all<-1-fnr_all
      specificity_all<-1-(fp_fun+fp_scalar+fp_dev)/(fp_fun+tn_fun+fp_dev+tn_dev+fp_scalar+tn_scalar)
      
      
      fp_all<-fp_fun+fp_scalar+fp_dev
      tp_all<-tp_fun +tp_scalar+tp_dev
      
      
      fn_all<-fn_fun+fn_scalar+fn_dev
      tn_all<-tn_fun +tn_scalar +tn_dev
      
      
      save_file<-paste("q",q,"true",q_true,"n",n1,"c",censoringvalues,"rho",rho,"_stability_selection.RData",sep="")
      
      save.image(save_file)
    }
  }
}

censoringvalues_all<-c(0.1,0.2,0.3,0.5)
rho_all<-c(0.5 )
n1_all<-c(200,300,500,1000,2000)
nruns<-100  
methods<-6
q<-100
q_true<-2
mfunctions<-q
mfunctions_true<-q_true

for (cen_i in 1:length(censoringvalues_all)) {
  for (n_j in 1:length(n1_all)) {
    for (rho_k in 1:length(rho_all)) {
      
      n1<-n1_all[n_j]
      
      censoringvalues<-censoringvalues_all[cen_i]
      rho<-rho_all[rho_k]
      tp_scalar<-matrix(0,nruns,methods)
      fp_scalar<-matrix(0,nruns,methods)
      tn_scalar<-matrix(0,nruns,methods)
      fn_scalar<-matrix(0,nruns,methods)
      
      tp_dev<-matrix(0,nruns,methods)
      fp_dev<-matrix(0,nruns,methods)
      tn_dev<-matrix(0,nruns,methods)
      fn_dev<-matrix(0,nruns,methods)
      
      tp_fun<-matrix(0,nruns,methods)
      fp_fun<-matrix(0,nruns,methods)
      tn_fun<-matrix(0,nruns,methods)
      fn_fun<-matrix(0,nruns,methods)
      
      censoring<-matrix(0,nruns,1)
      for(runs_i in 1:nruns){
        filename=paste(runs_i,"_rho",rho,"_censoringvalues",censoringvalues,"_n",n1,"_mfunctions",mfunctions,"_mfunctions_true",mfunctions_true,"_q",q,"_q_true",q_true,".csv",sep="")
        
        result<-read.csv(filename)
        censoring[runs_i]<-result[4,2]
        result<-result[-4,-1]
        tp_scalar[runs_i,]<-result[,1]
        fp_scalar[runs_i,]<-result[,2]
        tn_scalar[runs_i,]<-result[,3]
        fn_scalar[runs_i,]<-result[,4]
        tp_dev[runs_i,]<-result[,5]
        fp_dev[runs_i,]<-result[,6]
        tn_dev[runs_i,]<-result[,7]
        fn_dev[runs_i,]<-result[,8]
        tp_fun[runs_i,]<-result[,9] 
        fp_fun[runs_i,]<-result[,10]
        tn_fun[runs_i,]<-result[,11]
        fn_fun[runs_i,]<-result[,12]
      }
      tpr_scalar<-tp_scalar/(tp_scalar+fn_scalar)
      fpr_scalar<-fp_scalar/(fp_scalar+tn_scalar)
      tpr_fun<-tp_fun/(tp_fun+fn_fun)
      fpr_fun<-fp_fun/(fp_fun+tn_fun)
      tpr_dev<-tp_dev/(tp_dev+fn_dev)
      fpr_dev<-fp_dev/(fp_dev+tn_dev)
      
      
      
      fdr_scalar<-fp_scalar/(tp_scalar+fp_scalar)
      fnr_scalar<-fn_scalar/(tp_scalar+fn_scalar)
      sensitivity_scalar<-1-fnr_scalar
      specificity_scalar<-1-fp_scalar/(fp_scalar+tn_scalar)
      
      fdr_dev<-fp_dev/(tp_dev+fp_dev)
      fnr_dev<-fn_dev/(tp_dev+fn_dev)
      sensitivity_dev<-1-fnr_dev
      specificity_dev<-1-fp_dev/(fp_dev+tn_dev)
      
      
      fdr_fun<-fp_fun/(tp_fun+fp_fun)
      fnr_fun<-fn_fun/(tp_fun+fn_fun)
      sensitivity_fun<-1-fnr_fun
      specificity_fun<-1-fp_fun/(fp_fun+tn_fun)
      
      fdr_all<-(fp_fun+fp_scalar+fp_dev)/(tp_fun+fp_fun+tp_dev+fp_dev+tp_scalar+fp_scalar)
      fnr_all<-(fn_fun+fn_scalar+fn_dev)/(tp_fun+fn_fun+tp_dev+fn_dev+tp_scalar+fn_scalar)
      sensitivity_all<-1-fnr_all
      specificity_all<-1-(fp_fun+fp_scalar+fp_dev)/(fp_fun+tn_fun+fp_dev+tn_dev+fp_scalar+tn_scalar)
      
      
      fp_all<-fp_fun+fp_scalar+fp_dev
      tp_all<-tp_fun +tp_scalar+tp_dev
      
      
      fn_all<-fn_fun+fn_scalar+fn_dev
      tn_all<-tn_fun +tn_scalar +tn_dev
      
      
      save_file<-paste("q",q,"true",q_true,"n",n1,"c",censoringvalues,"rho",rho,"_stability_selection.RData",sep="")
      
      save.image(save_file)
    }
  }
}
censoringvalues_all<-c(0.1,0.2,0.3,0.5)
rho_all<-c(0.5 )
n1_all<-c(200,300,500,1000,2000)
nruns<-100  
methods<-6
q<-200
q_true<-2
mfunctions<-q
mfunctions_true<-q_true

for (cen_i in 1:length(censoringvalues_all)) {
  for (n_j in 1:length(n1_all)) {
    for (rho_k in 1:length(rho_all)) {
      
      n1<-n1_all[n_j]
      
      censoringvalues<-censoringvalues_all[cen_i]
      rho<-rho_all[rho_k]
      tp_scalar<-matrix(0,nruns,methods)
      fp_scalar<-matrix(0,nruns,methods)
      tn_scalar<-matrix(0,nruns,methods)
      fn_scalar<-matrix(0,nruns,methods)
      
      tp_dev<-matrix(0,nruns,methods)
      fp_dev<-matrix(0,nruns,methods)
      tn_dev<-matrix(0,nruns,methods)
      fn_dev<-matrix(0,nruns,methods)
      
      tp_fun<-matrix(0,nruns,methods)
      fp_fun<-matrix(0,nruns,methods)
      tn_fun<-matrix(0,nruns,methods)
      fn_fun<-matrix(0,nruns,methods)
      
      censoring<-matrix(0,nruns,1)
      for(runs_i in 1:nruns){
        filename=paste(runs_i,"_rho",rho,"_censoringvalues",censoringvalues,"_n",n1,"_mfunctions",mfunctions,"_mfunctions_true",mfunctions_true,"_q",q,"_q_true",q_true,".csv",sep="")
        
        result<-read.csv(filename)
        censoring[runs_i]<-result[4,2]
        result<-result[-4,-1]
        tp_scalar[runs_i,]<-result[,1]
        fp_scalar[runs_i,]<-result[,2]
        tn_scalar[runs_i,]<-result[,3]
        fn_scalar[runs_i,]<-result[,4]
        tp_dev[runs_i,]<-result[,5]
        fp_dev[runs_i,]<-result[,6]
        tn_dev[runs_i,]<-result[,7]
        fn_dev[runs_i,]<-result[,8]
        tp_fun[runs_i,]<-result[,9] 
        fp_fun[runs_i,]<-result[,10]
        tn_fun[runs_i,]<-result[,11]
        fn_fun[runs_i,]<-result[,12]
      }
      tpr_scalar<-tp_scalar/(tp_scalar+fn_scalar)
      fpr_scalar<-fp_scalar/(fp_scalar+tn_scalar)
      tpr_fun<-tp_fun/(tp_fun+fn_fun)
      fpr_fun<-fp_fun/(fp_fun+tn_fun)
      tpr_dev<-tp_dev/(tp_dev+fn_dev)
      fpr_dev<-fp_dev/(fp_dev+tn_dev)
      
      
      
      fdr_scalar<-fp_scalar/(tp_scalar+fp_scalar)
      fnr_scalar<-fn_scalar/(tp_scalar+fn_scalar)
      sensitivity_scalar<-1-fnr_scalar
      specificity_scalar<-1-fp_scalar/(fp_scalar+tn_scalar)
      
      fdr_dev<-fp_dev/(tp_dev+fp_dev)
      fnr_dev<-fn_dev/(tp_dev+fn_dev)
      sensitivity_dev<-1-fnr_dev
      specificity_dev<-1-fp_dev/(fp_dev+tn_dev)
      
      
      fdr_fun<-fp_fun/(tp_fun+fp_fun)
      fnr_fun<-fn_fun/(tp_fun+fn_fun)
      sensitivity_fun<-1-fnr_fun
      specificity_fun<-1-fp_fun/(fp_fun+tn_fun)
      
      fdr_all<-(fp_fun+fp_scalar+fp_dev)/(tp_fun+fp_fun+tp_dev+fp_dev+tp_scalar+fp_scalar)
      fnr_all<-(fn_fun+fn_scalar+fn_dev)/(tp_fun+fn_fun+tp_dev+fn_dev+tp_scalar+fn_scalar)
      sensitivity_all<-1-fnr_all
      specificity_all<-1-(fp_fun+fp_scalar+fp_dev)/(fp_fun+tn_fun+fp_dev+tn_dev+fp_scalar+tn_scalar)
      
      
      fp_all<-fp_fun+fp_scalar+fp_dev
      tp_all<-tp_fun +tp_scalar+tp_dev
      
      
      fn_all<-fn_fun+fn_scalar+fn_dev
      tn_all<-tn_fun +tn_scalar +tn_dev
      
      
      save_file<-paste("q",q,"true",q_true,"n",n1,"c",censoringvalues,"rho",rho,"_stability_selection.RData",sep="")
      
      save.image(save_file)
    }
  }
}
censoringvalues_all<-c(0.1,0.2,0.3,0.5)
rho_all<-c(0.5 )
n1_all<-c(200,300,500,1000,2000)
nruns<-100  
methods<-6
q<-300
q_true<-2
mfunctions<-q
mfunctions_true<-q_true

for (cen_i in 1:length(censoringvalues_all)) {
  for (n_j in 1:length(n1_all)) {
    for (rho_k in 1:length(rho_all)) {
      
      n1<-n1_all[n_j]
      
      censoringvalues<-censoringvalues_all[cen_i]
      rho<-rho_all[rho_k]
      tp_scalar<-matrix(0,nruns,methods)
      fp_scalar<-matrix(0,nruns,methods)
      tn_scalar<-matrix(0,nruns,methods)
      fn_scalar<-matrix(0,nruns,methods)
      
      tp_dev<-matrix(0,nruns,methods)
      fp_dev<-matrix(0,nruns,methods)
      tn_dev<-matrix(0,nruns,methods)
      fn_dev<-matrix(0,nruns,methods)
      
      tp_fun<-matrix(0,nruns,methods)
      fp_fun<-matrix(0,nruns,methods)
      tn_fun<-matrix(0,nruns,methods)
      fn_fun<-matrix(0,nruns,methods)
      
      censoring<-matrix(0,nruns,1)
      for(runs_i in 1:nruns){
        filename=paste(runs_i,"_rho",rho,"_censoringvalues",censoringvalues,"_n",n1,"_mfunctions",mfunctions,"_mfunctions_true",mfunctions_true,"_q",q,"_q_true",q_true,".csv",sep="")
        
        result<-read.csv(filename)
        censoring[runs_i]<-result[4,2]
        result<-result[-4,-1]
        tp_scalar[runs_i,]<-result[,1]
        fp_scalar[runs_i,]<-result[,2]
        tn_scalar[runs_i,]<-result[,3]
        fn_scalar[runs_i,]<-result[,4]
        tp_dev[runs_i,]<-result[,5]
        fp_dev[runs_i,]<-result[,6]
        tn_dev[runs_i,]<-result[,7]
        fn_dev[runs_i,]<-result[,8]
        tp_fun[runs_i,]<-result[,9] 
        fp_fun[runs_i,]<-result[,10]
        tn_fun[runs_i,]<-result[,11]
        fn_fun[runs_i,]<-result[,12]
      }
      tpr_scalar<-tp_scalar/(tp_scalar+fn_scalar)
      fpr_scalar<-fp_scalar/(fp_scalar+tn_scalar)
      tpr_fun<-tp_fun/(tp_fun+fn_fun)
      fpr_fun<-fp_fun/(fp_fun+tn_fun)
      tpr_dev<-tp_dev/(tp_dev+fn_dev)
      fpr_dev<-fp_dev/(fp_dev+tn_dev)
      
      
      
      fdr_scalar<-fp_scalar/(tp_scalar+fp_scalar)
      fnr_scalar<-fn_scalar/(tp_scalar+fn_scalar)
      sensitivity_scalar<-1-fnr_scalar
      specificity_scalar<-1-fp_scalar/(fp_scalar+tn_scalar)
      
      fdr_dev<-fp_dev/(tp_dev+fp_dev)
      fnr_dev<-fn_dev/(tp_dev+fn_dev)
      sensitivity_dev<-1-fnr_dev
      specificity_dev<-1-fp_dev/(fp_dev+tn_dev)
      
      
      fdr_fun<-fp_fun/(tp_fun+fp_fun)
      fnr_fun<-fn_fun/(tp_fun+fn_fun)
      sensitivity_fun<-1-fnr_fun
      specificity_fun<-1-fp_fun/(fp_fun+tn_fun)
      
      fdr_all<-(fp_fun+fp_scalar+fp_dev)/(tp_fun+fp_fun+tp_dev+fp_dev+tp_scalar+fp_scalar)
      fnr_all<-(fn_fun+fn_scalar+fn_dev)/(tp_fun+fn_fun+tp_dev+fn_dev+tp_scalar+fn_scalar)
      sensitivity_all<-1-fnr_all
      specificity_all<-1-(fp_fun+fp_scalar+fp_dev)/(fp_fun+tn_fun+fp_dev+tn_dev+fp_scalar+tn_scalar)
      
      
      fp_all<-fp_fun+fp_scalar+fp_dev
      tp_all<-tp_fun +tp_scalar+tp_dev
      
      
      fn_all<-fn_fun+fn_scalar+fn_dev
      tn_all<-tn_fun +tn_scalar +tn_dev
      
      
      save_file<-paste("q",q,"true",q_true,"n",n1,"c",censoringvalues,"rho",rho,"_stability_selection.RData",sep="")
      
      save.image(save_file)
    }
  }
}
censoringvalues_all<-c(0.1,0.2,0.3,0.5)
rho_all<-c(0.3 )
n1_all<-c(200,300,500,1000,2000)
nruns<-100  
methods<-6
q<-300
q_true<-2
mfunctions<-q
mfunctions_true<-q_true

for (cen_i in 1:length(censoringvalues_all)) {
  for (n_j in 1:length(n1_all)) {
    for (rho_k in 1:length(rho_all)) {
      
      n1<-n1_all[n_j]
      
      censoringvalues<-censoringvalues_all[cen_i]
      rho<-rho_all[rho_k]
      tp_scalar<-matrix(0,nruns,methods)
      fp_scalar<-matrix(0,nruns,methods)
      tn_scalar<-matrix(0,nruns,methods)
      fn_scalar<-matrix(0,nruns,methods)
      
      tp_dev<-matrix(0,nruns,methods)
      fp_dev<-matrix(0,nruns,methods)
      tn_dev<-matrix(0,nruns,methods)
      fn_dev<-matrix(0,nruns,methods)
      
      tp_fun<-matrix(0,nruns,methods)
      fp_fun<-matrix(0,nruns,methods)
      tn_fun<-matrix(0,nruns,methods)
      fn_fun<-matrix(0,nruns,methods)
      
      censoring<-matrix(0,nruns,1)
      for(runs_i in 1:nruns){
        filename=paste(runs_i,"_rho",rho,"_censoringvalues",censoringvalues,"_n",n1,"_mfunctions",mfunctions,"_mfunctions_true",mfunctions_true,"_q",q,"_q_true",q_true,".csv",sep="")
        
        result<-read.csv(filename)
        censoring[runs_i]<-result[4,2]
        result<-result[-4,-1]
        tp_scalar[runs_i,]<-result[,1]
        fp_scalar[runs_i,]<-result[,2]
        tn_scalar[runs_i,]<-result[,3]
        fn_scalar[runs_i,]<-result[,4]
        tp_dev[runs_i,]<-result[,5]
        fp_dev[runs_i,]<-result[,6]
        tn_dev[runs_i,]<-result[,7]
        fn_dev[runs_i,]<-result[,8]
        tp_fun[runs_i,]<-result[,9] 
        fp_fun[runs_i,]<-result[,10]
        tn_fun[runs_i,]<-result[,11]
        fn_fun[runs_i,]<-result[,12]
      }
      tpr_scalar<-tp_scalar/(tp_scalar+fn_scalar)
      fpr_scalar<-fp_scalar/(fp_scalar+tn_scalar)
      tpr_fun<-tp_fun/(tp_fun+fn_fun)
      fpr_fun<-fp_fun/(fp_fun+tn_fun)
      tpr_dev<-tp_dev/(tp_dev+fn_dev)
      fpr_dev<-fp_dev/(fp_dev+tn_dev)
      
      
      
      fdr_scalar<-fp_scalar/(tp_scalar+fp_scalar)
      fnr_scalar<-fn_scalar/(tp_scalar+fn_scalar)
      sensitivity_scalar<-1-fnr_scalar
      specificity_scalar<-1-fp_scalar/(fp_scalar+tn_scalar)
      
      fdr_dev<-fp_dev/(tp_dev+fp_dev)
      fnr_dev<-fn_dev/(tp_dev+fn_dev)
      sensitivity_dev<-1-fnr_dev
      specificity_dev<-1-fp_dev/(fp_dev+tn_dev)
      
      
      fdr_fun<-fp_fun/(tp_fun+fp_fun)
      fnr_fun<-fn_fun/(tp_fun+fn_fun)
      sensitivity_fun<-1-fnr_fun
      specificity_fun<-1-fp_fun/(fp_fun+tn_fun)
      
      fdr_all<-(fp_fun+fp_scalar+fp_dev)/(tp_fun+fp_fun+tp_dev+fp_dev+tp_scalar+fp_scalar)
      fnr_all<-(fn_fun+fn_scalar+fn_dev)/(tp_fun+fn_fun+tp_dev+fn_dev+tp_scalar+fn_scalar)
      sensitivity_all<-1-fnr_all
      specificity_all<-1-(fp_fun+fp_scalar+fp_dev)/(fp_fun+tn_fun+fp_dev+tn_dev+fp_scalar+tn_scalar)
      
      
      fp_all<-fp_fun+fp_scalar+fp_dev
      tp_all<-tp_fun +tp_scalar+tp_dev
      
      
      fn_all<-fn_fun+fn_scalar+fn_dev
      tn_all<-tn_fun +tn_scalar +tn_dev
      
      
      save_file<-paste("q",q,"true",q_true,"n",n1,"c",censoringvalues,"rho",rho,"_stability_selection.RData",sep="")
      
      save.image(save_file)
    }
  }
}
censoringvalues_all<-c(0.1,0.2,0.3,0.5)
rho_all<-c(0.8 )
n1_all<-c(200,300,500,1000,2000)
nruns<-100  
methods<-6
q<-300
q_true<-2
mfunctions<-q
mfunctions_true<-q_true

for (cen_i in 1:length(censoringvalues_all)) {
  for (n_j in 1:length(n1_all)) {
    for (rho_k in 1:length(rho_all)) {
      
      n1<-n1_all[n_j]
      
      censoringvalues<-censoringvalues_all[cen_i]
      rho<-rho_all[rho_k]
      tp_scalar<-matrix(0,nruns,methods)
      fp_scalar<-matrix(0,nruns,methods)
      tn_scalar<-matrix(0,nruns,methods)
      fn_scalar<-matrix(0,nruns,methods)
      
      tp_dev<-matrix(0,nruns,methods)
      fp_dev<-matrix(0,nruns,methods)
      tn_dev<-matrix(0,nruns,methods)
      fn_dev<-matrix(0,nruns,methods)
      
      tp_fun<-matrix(0,nruns,methods)
      fp_fun<-matrix(0,nruns,methods)
      tn_fun<-matrix(0,nruns,methods)
      fn_fun<-matrix(0,nruns,methods)
      
      censoring<-matrix(0,nruns,1)
      for(runs_i in 1:nruns){
        filename=paste(runs_i,"_rho",rho,"_censoringvalues",censoringvalues,"_n",n1,"_mfunctions",mfunctions,"_mfunctions_true",mfunctions_true,"_q",q,"_q_true",q_true,".csv",sep="")
        
        result<-read.csv(filename)
        censoring[runs_i]<-result[4,2]
        result<-result[-4,-1]
        tp_scalar[runs_i,]<-result[,1]
        fp_scalar[runs_i,]<-result[,2]
        tn_scalar[runs_i,]<-result[,3]
        fn_scalar[runs_i,]<-result[,4]
        tp_dev[runs_i,]<-result[,5]
        fp_dev[runs_i,]<-result[,6]
        tn_dev[runs_i,]<-result[,7]
        fn_dev[runs_i,]<-result[,8]
        tp_fun[runs_i,]<-result[,9] 
        fp_fun[runs_i,]<-result[,10]
        tn_fun[runs_i,]<-result[,11]
        fn_fun[runs_i,]<-result[,12]
      }
      tpr_scalar<-tp_scalar/(tp_scalar+fn_scalar)
      fpr_scalar<-fp_scalar/(fp_scalar+tn_scalar)
      tpr_fun<-tp_fun/(tp_fun+fn_fun)
      fpr_fun<-fp_fun/(fp_fun+tn_fun)
      tpr_dev<-tp_dev/(tp_dev+fn_dev)
      fpr_dev<-fp_dev/(fp_dev+tn_dev)
      
      
      
      fdr_scalar<-fp_scalar/(tp_scalar+fp_scalar)
      fnr_scalar<-fn_scalar/(tp_scalar+fn_scalar)
      sensitivity_scalar<-1-fnr_scalar
      specificity_scalar<-1-fp_scalar/(fp_scalar+tn_scalar)
      
      fdr_dev<-fp_dev/(tp_dev+fp_dev)
      fnr_dev<-fn_dev/(tp_dev+fn_dev)
      sensitivity_dev<-1-fnr_dev
      specificity_dev<-1-fp_dev/(fp_dev+tn_dev)
      
      
      fdr_fun<-fp_fun/(tp_fun+fp_fun)
      fnr_fun<-fn_fun/(tp_fun+fn_fun)
      sensitivity_fun<-1-fnr_fun
      specificity_fun<-1-fp_fun/(fp_fun+tn_fun)
      
      fdr_all<-(fp_fun+fp_scalar+fp_dev)/(tp_fun+fp_fun+tp_dev+fp_dev+tp_scalar+fp_scalar)
      fnr_all<-(fn_fun+fn_scalar+fn_dev)/(tp_fun+fn_fun+tp_dev+fn_dev+tp_scalar+fn_scalar)
      sensitivity_all<-1-fnr_all
      specificity_all<-1-(fp_fun+fp_scalar+fp_dev)/(fp_fun+tn_fun+fp_dev+tn_dev+fp_scalar+tn_scalar)
      
      
      fp_all<-fp_fun+fp_scalar+fp_dev
      tp_all<-tp_fun +tp_scalar+tp_dev
      
      
      fn_all<-fn_fun+fn_scalar+fn_dev
      tn_all<-tn_fun +tn_scalar +tn_dev
      
      
      save_file<-paste("q",q,"true",q_true,"n",n1,"c",censoringvalues,"rho",rho,"_stability_selection.RData",sep="")
      
      save.image(save_file)
    }
  }
}
censoringvalues_all<-c(0.1,0.2,0.3,0.5)
rho_all<-c(0.5 )
n1_all<-c(200,300,500,1000,2000)
nruns<-100  
methods<-6
q<-300
q_true<-2
mfunctions<-q
mfunctions_true<-q_true

for (cen_i in 1:length(censoringvalues_all)) {
  for (n_j in 1:length(n1_all)) {
    for (rho_k in 1:length(rho_all)) {
      
      n1<-n1_all[n_j]
      
      censoringvalues<-censoringvalues_all[cen_i]
      rho<-rho_all[rho_k]
      tp_scalar<-matrix(0,nruns,methods)
      fp_scalar<-matrix(0,nruns,methods)
      tn_scalar<-matrix(0,nruns,methods)
      fn_scalar<-matrix(0,nruns,methods)
      
      tp_dev<-matrix(0,nruns,methods)
      fp_dev<-matrix(0,nruns,methods)
      tn_dev<-matrix(0,nruns,methods)
      fn_dev<-matrix(0,nruns,methods)
      
      tp_fun<-matrix(0,nruns,methods)
      fp_fun<-matrix(0,nruns,methods)
      tn_fun<-matrix(0,nruns,methods)
      fn_fun<-matrix(0,nruns,methods)
      
      censoring<-matrix(0,nruns,1)
      for(runs_i in 1:nruns){
        filename=paste(runs_i,"_rho",rho,"_censoringvalues",censoringvalues,"_n",n1,"_mfunctions",mfunctions,"_mfunctions_true",mfunctions_true,"_q",q,"_q_true",q_true,".csv",sep="")
        
        result<-read.csv(filename)
        censoring[runs_i]<-result[4,2]
        result<-result[-4,-1]
        tp_scalar[runs_i,]<-result[,1]
        fp_scalar[runs_i,]<-result[,2]
        tn_scalar[runs_i,]<-result[,3]
        fn_scalar[runs_i,]<-result[,4]
        tp_dev[runs_i,]<-result[,5]
        fp_dev[runs_i,]<-result[,6]
        tn_dev[runs_i,]<-result[,7]
        fn_dev[runs_i,]<-result[,8]
        tp_fun[runs_i,]<-result[,9] 
        fp_fun[runs_i,]<-result[,10]
        tn_fun[runs_i,]<-result[,11]
        fn_fun[runs_i,]<-result[,12]
      }
      tpr_scalar<-tp_scalar/(tp_scalar+fn_scalar)
      fpr_scalar<-fp_scalar/(fp_scalar+tn_scalar)
      tpr_fun<-tp_fun/(tp_fun+fn_fun)
      fpr_fun<-fp_fun/(fp_fun+tn_fun)
      tpr_dev<-tp_dev/(tp_dev+fn_dev)
      fpr_dev<-fp_dev/(fp_dev+tn_dev)
      
      
      
      fdr_scalar<-fp_scalar/(tp_scalar+fp_scalar)
      fnr_scalar<-fn_scalar/(tp_scalar+fn_scalar)
      sensitivity_scalar<-1-fnr_scalar
      specificity_scalar<-1-fp_scalar/(fp_scalar+tn_scalar)
      
      fdr_dev<-fp_dev/(tp_dev+fp_dev)
      fnr_dev<-fn_dev/(tp_dev+fn_dev)
      sensitivity_dev<-1-fnr_dev
      specificity_dev<-1-fp_dev/(fp_dev+tn_dev)
      
      
      fdr_fun<-fp_fun/(tp_fun+fp_fun)
      fnr_fun<-fn_fun/(tp_fun+fn_fun)
      sensitivity_fun<-1-fnr_fun
      specificity_fun<-1-fp_fun/(fp_fun+tn_fun)
      
      fdr_all<-(fp_fun+fp_scalar+fp_dev)/(tp_fun+fp_fun+tp_dev+fp_dev+tp_scalar+fp_scalar)
      fnr_all<-(fn_fun+fn_scalar+fn_dev)/(tp_fun+fn_fun+tp_dev+fn_dev+tp_scalar+fn_scalar)
      sensitivity_all<-1-fnr_all
      specificity_all<-1-(fp_fun+fp_scalar+fp_dev)/(fp_fun+tn_fun+fp_dev+tn_dev+fp_scalar+tn_scalar)
      
      
      fp_all<-fp_fun+fp_scalar+fp_dev
      tp_all<-tp_fun +tp_scalar+tp_dev
      
      
      fn_all<-fn_fun+fn_scalar+fn_dev
      tn_all<-tn_fun +tn_scalar +tn_dev
      
      
      save_file<-paste("q",q,"true",q_true,"n",n1,"c",censoringvalues,"rho",rho,"_stability_selection.RData",sep="")
      
      save.image(save_file)
    }
  }
}
