simulation<-function(runs_i){
  
  #runs_i<-1
  set.seed(runs_i)
  sample_i<-sample(1:1e5,2,replace = FALSE)
  #library(locpol)
  library(MASS)
  #library(grpCox)
  library(grpreg)
  library(survcomp)
  library(fda)
  #time1<-Sys.time() 
  
  n1<-2000
  n2<-0+n1
  n<-0+2000
  
  use_n_all<-c(200,300,500,1000,2000)
  
  mfunctions<-300
  mfunctions_true<-2
  
  mfunctions_true_dev<-2
  q<-300
  q_true<-2 
  
  
  rho<-0.5 
  
  h<-calh<-0.01
  designpoints<-seq(0,1,calh)
  
  t<-1/h+1
  
  numberofbasis<-10
  fourierbasis10 <- create.fourier.basis(rangeval=c(0,1),nbasis=numberofbasis+1,period=1,dropind = 1)#
  plot(fourierbasis10)
  
  
  fourierbasis11 <- create.fourier.basis(rangeval=c(0,1),nbasis=numberofbasis+1,period=1)#
  plot(fourierbasis11)
  
  
  little_Sigma<-matrix(0,10,10)
  
  for (j in 1:5){
    little_Sigma[j,j]<- (j)^(-2)#16*
  }
  for (j in 6:10){
    little_Sigma[j,j]<-(j-5)^(-2)	 #16*
  }
  
  
  
  Sigma<-matrix(0,nrow(little_Sigma)*mfunctions,nrow(little_Sigma)*mfunctions)
  
  for(mfunctions_i in 1:mfunctions){
    Sigma[((mfunctions_i-1)*nrow(little_Sigma)+1):(mfunctions_i*nrow(little_Sigma)),((mfunctions_i-1)*nrow(little_Sigma)+1):(mfunctions_i*nrow(little_Sigma))]<-little_Sigma
    
  }
  
  
  
  
  truebeta_list<-list()
  
  truebeta_list[[1]]<-fd(c(0.5, rep(0,9)),fourierbasis10)
  truebeta_list[[2]]<-fd(c(0,-0.3,rep(0,8)),fourierbasis10)
  
  for (turebetalist_i in (mfunctions_true+1):mfunctions) {
    truebeta_list[[turebetalist_i]]<-fd(  rep(0,numberofbasis) ,fourierbasis10)
    
  }
  
  
  
  #summary(rawX[[yita_i]]%*%truebeta[,yita_i]*calh)
  
  
  true_beta_flag<-c(rep(1,mfunctions_true), rep(0,mfunctions-mfunctions_true))
  
  gamma<-c(c(0.5, -0.3  ),rep(0,q-q_true))
  
  eta_dev<-c(0,0.001, -0.003 ,rep(0,mfunctions-1-mfunctions_true_dev))
  
  correlationmatrix<-matrix(0,nrow(Sigma),q)
  
  
  SigmaZ<-rho^t(sapply(1:q, function(i, j) abs(i-j), 1:q))
  
  Bigsigma<-rbind(cbind(Sigma, correlationmatrix), cbind(t(correlationmatrix), SigmaZ))
  
  mu<-rep(0,nrow(Bigsigma))
  set.seed(sample_i[1])
  data<-mvrnorm(n,mu,Bigsigma)
  
  
  ximatrix<-list()
  for(ximatrix_i in 1:mfunctions){
    ximatrix[[ximatrix_i]]<-data[,((ximatrix_i-1)*nrow(little_Sigma)+1):(ximatrix_i*nrow(little_Sigma))]
  }
  
  Vxmatrix<-ximatrix
  for(vmatrix_i in 1:mfunctions){
    
    if((vmatrix_i>0)&(vmatrix_i<=(mfunctions-2))){
      Vxmatrix[[vmatrix_i]]<-ximatrix[[vmatrix_i]] +rho*(ximatrix[[vmatrix_i+1]]+ximatrix[[vmatrix_i+2]])
    }
    
    if(vmatrix_i ==(mfunctions-1)){
      Vxmatrix[[vmatrix_i]]<-ximatrix[[vmatrix_i]] +rho*(ximatrix[[vmatrix_i+1]]+ximatrix[[1]])
    }
    if(vmatrix_i ==mfunctions){
      Vxmatrix[[vmatrix_i]]<-ximatrix[[vmatrix_i]] +rho*(ximatrix[[1]]+ximatrix[[2]])
    }
  }
  
  
  
  
  rawX<-list()
  observeX<-list()
  X<-list()
  dev1<-list()
  dev1_int_true<-matrix(0,n,mfunctions)
  dev1_int_obs<-matrix(0,n,mfunctions)
  
  for(mfunctions_i in 1:mfunctions){
    rawX[[mfunctions_i]]<- fd(t(Vxmatrix[[mfunctions_i]]),fourierbasis10)
    dev1[[mfunctions_i]]<-deriv.fd( rawX[[mfunctions_i]],2)
    
    
    for(dev1_i in 1:n){
      dev1fd_i<-fd(dev1[[mfunctions_i]]$coefs[,dev1_i],fourierbasis10)
      dev1_int_true[dev1_i,mfunctions_i] <-sqrt(inprod(dev1fd_i,dev1fd_i))
    }
    
    
    observeX[[mfunctions_i]]<- t(eval.fd(designpoints, rawX[[mfunctions_i]]))
    
    fit_fd <-smooth.basis(argvals=designpoints,y=t(observeX[[mfunctions_i]]) ,fourierbasis11) $fd
    X[[mfunctions_i]] <-fit_fd
    
    dev1_obs_fd<-deriv.fd( X[[mfunctions_i]],2)
    
    for(dev1_i in 1:n){
      dev1fd_i<-fd(dev1_obs_fd$coefs[,dev1_i],dev1_obs_fd$basis)
      dev1_int_obs[dev1_i,mfunctions_i] <-sqrt(inprod(dev1fd_i,dev1fd_i))
    }
    
    
  }
  
  
  
  
  
  
  Z<-data[,(nrow(Sigma)+1):(nrow(Sigma)+q)]
  
  
  ###### generate the time-to-event data
  yita<-rep(0,n)
  for(yita_i in 1:mfunctions_true){
    yita<-yita+ inprod(rawX[[yita_i]],truebeta_list[[yita_i]])  #observeX[[yita_i]]%*%truebeta[,yita_i]*calh
  }
  parameter<-exp(yita+Z%*%gamma+dev1_int_true%*%eta_dev)
  
  summary(yita)
  
  summary(Z%*%gamma)
  
  summary(dev1_int_true%*%eta_dev)
  
  
  
  failuretime<-rexp(n, rate = parameter)
  
  
  ##############################################################
  ##############################################################
  ######################censoringvalues all  ######################
  ##############################################################
  ##############################################################
  censoringvalues_all<-c( 0.1,0.2,0.3,0.5)
  tau_all<-c( 80,30,18,7)
  
  tau<-7  
  set.seed(sample_i[2])
  censoringtime<-runif(n,0,tau) 
  event<-rep(0,n)
  event<-as.numeric(failuretime<censoringtime)
  sum(event)/n  
  
  for (use_n in use_n_all) {
    ############################################################
    #####################    get fpca       ####################
    ############################################################
    
    
    p_group<-c()
    fun_Eigval<-list()
    fun_Eigfun<-list()
    eigenscore<-list()
    sn_chose<-rep(0,mfunctions)
    
    all_designmatrix<-cbind(Z[1:use_n,],dev1_int_obs[1:use_n,])
    for(pca_i in 1:mfunctions){
      use_x_n<-fd(coef =( X[[pca_i]]$coefs)[,1:use_n] ,basisobj= X[[pca_i]]$basis)
      pca_result<-pca.fd(use_x_n, 10) 
      sn_chose[pca_i]<-which(cumsum(pca_result$varprop)>0.95)[1]
      pca_result<-pca.fd(use_x_n, sn_chose[pca_i]) 
      fun_Eigfun[[pca_i]]=pca_result$harmonics
      eigenscore[[pca_i]]<-pca_result$scores
      fun_Eigval[[pca_i]]<-pca_result$values[1:(sn_chose[pca_i])]
      all_designmatrix<-cbind(all_designmatrix,eigenscore[[pca_i]])
    }
    
    p_group <-c(1:q,(1:mfunctions)+q,rep((1:mfunctions)+mfunctions +q, sn_chose ))
    m_fun<-c()
    for(m_i in 1:mfunctions){
      m_fun[m_i]<-sqrt(sum(fun_Eigval[[m_i]]))
    }
    
    m <- c(rep(1,q),rep(1,mfunctions),m_fun)
     
    designmatrix<-all_designmatrix
    for (tau_j in 1:length(tau_all)) {
      ##############################################################
      ##############################################################
      ######################censoringvalues  ######################
      ##############################################################
      ##############################################################
      
      censoringvalues<-censoringvalues_all[tau_j]
      tau<- tau_all[tau_j]
      
      
      set.seed(sample_i[2])
      censoringtime<-runif(n,0,tau)
      
      
      event<-rep(0,n)
      event<-as.numeric(failuretime<censoringtime)
      sum(event)/n  
      
      time<-failuretime*event+censoringtime*(rep(1,n)-event)
      
      
      
      
      
      ############################################################
      ############################################################
      #####################       fit models   ###################
      ############################################################
      ############################################################
      
      
      all_y <- data.frame(illt=time[1:n2], ills=event[1:n2])
      names(all_y) <- c("time", "status")
      
      
      y <- data.frame(illt=time[1:use_n], ills=event[1:use_n])
      names(y) <- c("time", "status")
      set.seed(2023)
      myfold<-sample(1:5,use_n,replace = TRUE)
      ############################################################################
      
      use_method<-c("grLasso","grSCAD","grMCP")
      result<-matrix(0, nrow = 4, ncol =  4+4+4) 
      result_sta<-matrix(0, nrow = 3 , ncol =  4+4+4+2) 
      for(method_use_i in 1:length(use_method)){
        fit <- cv.grpsurv( X=designmatrix , y=y,group=p_group ,group.multiplier=m,  penalty=use_method[method_use_i] , fold = myfold  )
        
        
        
        fit_coef<- fit$fit$beta[,fit$min] 
        
        fit_coef_scalar<-fit_coef[1:q]
        tp_scalar<-sum(fit_coef_scalar[which(gamma!=0)]!=0)
        fp_scalar<-sum(fit_coef_scalar!=0)-tp_scalar
        tn_scalar<-sum(fit_coef_scalar[which(gamma==0)]==0)
        fn_scalar<-sum(fit_coef_scalar==0)-tn_scalar
        
        fit_coef_dev<-fit_coef[(q+1):(q+mfunctions)]
        tp_dev<-sum(fit_coef_dev[which(eta_dev!=0)]!=0)
        fp_dev<-sum(fit_coef_dev!=0)-tp_dev
        tn_dev<-sum(fit_coef_dev[which(eta_dev==0)]==0)
        fn_dev<-sum(fit_coef_dev==0)-tn_dev
        
        
        
        
        fit_coef_fun<-fit_coef[-c(1:(q+mfunctions))]
        fit_coef_fun_flag<-rep(0,mfunctions)
        for(subfit_i in 1:mfunctions){
          if(sum(abs(fit_coef_fun[which(as.numeric(p_group[-c(1:(q+mfunctions))]-q-mfunctions)==subfit_i)]))>0) fit_coef_fun_flag[subfit_i]<-1
        }
        tp_fun<-sum(fit_coef_fun_flag[which(true_beta_flag!=0)]!=0)
        fp_fun<-sum(fit_coef_fun_flag!=0)-tp_fun
        tn_fun<-sum(fit_coef_fun_flag[which(true_beta_flag==0)]==0)
        fn_fun<-sum(fit_coef_fun_flag==0)-tn_fun
        
        
        
        result[method_use_i,]<-c(tp_scalar,fp_scalar,tn_scalar,fn_scalar,
                                 tp_dev,fp_dev,tn_dev,fn_dev,
                                 tp_fun,fp_fun,tn_fun,fn_fun  )
        
        #####################    Stability  grLasso  ###################
        control_ev<-1
        thresholdvalues_lambda<-sqrt(0.8*(mfunctions+mfunctions+q) *control_ev)
        set.seed(9)
        #################
        
        
        fit <- cv.grpsurv(designmatrix,y,p_group,group.multiplier=m,penalty=use_method[method_use_i], fold = myfold)
        
        
        lam<-matrix(0,dim(fit$fit$beta)[2], q+mfunctions+mfunctions)
        for (lambdai in 1:(dim(fit$fit$beta)[2])) {
          fit_coef_fun<-fit$fit$beta[-c(1:(q+mfunctions)),lambdai]
          fit_coef_fun_flag<-rep(0,mfunctions)
          for(subfit_i in 1:mfunctions){
            if(sum(abs(fit_coef_fun[which(as.numeric(p_group)==(subfit_i+(q+mfunctions)))-(q+mfunctions)]))>0) fit_coef_fun_flag[subfit_i]<-1
          }
          lam[lambdai,]<- c(as.numeric(fit$fit$beta[c(1:(q+mfunctions)),lambdai]!=0),fit_coef_fun_flag)
        }
        
        sel_num<-apply(lam!=0, 1, sum)
        
        start_id<-1
        if(sel_num[1]==0) {
          start_id<-min(which(sel_num==0))
        }
        
        for (ss in (start_id+1):length(sel_num)) {
          if (sum(apply(lam[start_id:ss,]!=0, 2, sum)!=0) >=thresholdvalues_lambda)
            break
        }
        
        ss<-ss-1
        
        
        
        newlam<-fit$fit$lambda[start_id:ss]
        numLambda = length(newlam)
        
        numSample = 100
        nz = array(NA,c(numLambda,length(p_group),numSample))
        for(ii in 1:numSample) {
          cat(ii,fill=T)
          sampleInds = sample(1:use_n,size=use_n/2)
          
          aaa=grpsurv(designmatrix[sampleInds,],y[sampleInds,],p_group,group.multiplier=m,penalty=use_method[method_use_i],lambda =newlam )
          usecoe<-aaa$beta
          while(dim(usecoe)[2]<numLambda){
            cat(ii,fill=T)
            sampleInds = sample(1:use_n,size=use_n/2)
            
            aaa=grpsurv(designmatrix[sampleInds,],y[sampleInds,],p_group,group.multiplier=m,penalty=use_method[method_use_i],lambda =newlam )
            usecoe<-aaa$beta
            while(dim(usecoe)[2]<numLambda){
              cat(ii,fill=T)
              sampleInds = sample(1:use_n,size=use_n/2)
              
              aaa=grpsurv(designmatrix[sampleInds,],y[sampleInds,],p_group,group.multiplier=m,penalty=use_method[method_use_i],lambda =newlam )
              usecoe<-aaa$beta
              while(dim(usecoe)[2]<numLambda){
                cat(ii,fill=T)
                sampleInds = sample(1:use_n,size=use_n/2)
                
                aaa=grpsurv(designmatrix[sampleInds,],y[sampleInds,],p_group,group.multiplier=m,penalty=use_method[method_use_i],lambda =newlam )
                usecoe<-aaa$beta
                while(dim(usecoe)[2]<numLambda){
                  cat(ii,fill=T)
                  sampleInds = sample(1:use_n,size=use_n/2)
                  
                  aaa=grpsurv(designmatrix[sampleInds,],y[sampleInds,],p_group,group.multiplier=m,penalty=use_method[method_use_i],lambda =newlam )
                  usecoe<-aaa$beta
                  while(dim(usecoe)[2]<numLambda){
                    cat(ii,fill=T)
                    sampleInds = sample(1:use_n,size=use_n/2)
                    
                    aaa=grpsurv(designmatrix[sampleInds,],y[sampleInds,],p_group,group.multiplier=m,penalty=use_method[method_use_i],lambda =newlam )
                    usecoe<-aaa$beta
                    while(dim(usecoe)[2]<numLambda){
                      cat(ii,fill=T)
                      sampleInds = sample(1:use_n,size=use_n/2)
                      
                      aaa=grpsurv(designmatrix[sampleInds,],y[sampleInds,],p_group,group.multiplier=m,penalty=use_method[method_use_i],lambda =newlam )
                      usecoe<-aaa$beta
                      while(dim(usecoe)[2]<numLambda){
                        cat(ii,fill=T)
                        sampleInds = sample(1:use_n,size=use_n/2)
                        
                        aaa=grpsurv(designmatrix[sampleInds,],y[sampleInds,],p_group,group.multiplier=m,penalty=use_method[method_use_i],lambda =newlam )
                        usecoe<-aaa$beta
                        while(dim(usecoe)[2]<numLambda){
                          cat(ii,fill=T)
                          sampleInds = sample(1:use_n,size=use_n/2)
                          
                          aaa=grpsurv(designmatrix[sampleInds,],y[sampleInds,],p_group,group.multiplier=m,penalty=use_method[method_use_i],lambda =newlam )
                          usecoe<-aaa$beta
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          
          
          nz[,,ii]=t(1*(usecoe!=0))
        }
        
        pr_glasso=apply(nz,c(1,2),sum)/numSample
        
        # compute power and FDR
        maxPr = apply(pr_glasso,2,max);
        STAL<-c(0.6)
        selVar = maxPr >= STAL;
        
        
        fit_coef_scalar<-as.numeric(selVar)[1:q]
        tp_scalar<-sum(fit_coef_scalar[which(gamma!=0)]!=0)
        fp_scalar<-sum(fit_coef_scalar!=0)-tp_scalar
        tn_scalar<-sum(fit_coef_scalar[which(gamma==0)]==0)
        fn_scalar<-sum(fit_coef_scalar==0)-tn_scalar
        
        fit_coef_dev<-as.numeric(selVar)[(q+1):(q+mfunctions)]
        tp_dev<-sum(fit_coef_dev[which(eta_dev!=0)]!=0)
        fp_dev<-sum(fit_coef_dev!=0)-tp_dev
        tn_dev<-sum(fit_coef_dev[which(eta_dev==0)]==0)
        fn_dev<-sum(fit_coef_dev==0)-tn_dev
        
        fit_coef_fun<-as.numeric(selVar)[-c(1:(q+mfunctions))]
        fit_coef_fun_flag<-rep(0,mfunctions)
        for(subfit_i in 1:mfunctions){
          if(sum(abs(fit_coef_fun[which(as.numeric(p_group[-c(1:(q+mfunctions))]-q-mfunctions)==subfit_i)]))>0) fit_coef_fun_flag[subfit_i]<-1
          
        }
        tp_fun<-sum(fit_coef_fun_flag[which(true_beta_flag!=0)]!=0)
        fp_fun<-sum(fit_coef_fun_flag!=0)-tp_fun
        tn_fun<-sum(fit_coef_fun_flag[which(true_beta_flag==0)]==0)
        fn_fun<-sum(fit_coef_fun_flag==0)-tn_fun
        
        
        Stability_grLasso<-c(tp_scalar,fp_scalar,tn_scalar,fn_scalar,
                             tp_dev,fp_dev,tn_dev,fn_dev,
                             tp_fun,fp_fun,tn_fun,fn_fun ,STAL,control_ev)
        
        result_sta[method_use_i,]<-Stability_grLasso
        
        
        
      }   
      
      
      
      
      
      
      ############################################################################
      
      result[4,1]<-1-sum(event)/n
      
      result_all<-cbind(result,matrix(0,nrow = 4,ncol = 2))
      
      result_all<-rbind(result_all,result_sta)
      
      rownames(result_all)<-c( "result_grouplasso","result_groupSCAD","result_groupMCP","censoringvalues","s_glasso","s_glSCAD","s_gMCP")
      colnames(result_all)<-c( "tp_scalar","fp_scalar","tn_scalar","fn_scalar","tp_dev","fp_dev","tn_dev","fn_dev","tp_fun","fp_fun","tn_fun","fn_fun", "STAL","control_ev" )
      result_all
      
      filename=paste(runs_i,"_rho",rho,"_censoringvalues",censoringvalues,"_n",use_n,"_mfunctions",mfunctions,"_mfunctions_true",mfunctions_true,"_q",q,"_q_true",q_true,".csv",sep="")
      write.csv(result_all,filename)
      
      
      
    }
    
    
    
  }
  
  
  
}
