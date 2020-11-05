
library(survival)
library(boot)
library(ph2bayes)
library(reshape2)
library(ggplot2)
library(tictoc)
####################################################################
simTrial_lPFS1=function(pt=9, n=500, ppm=20, lPFS=0.6, fol_cut = "Yes", mfol =10,
                        npts=200, drop ="Yes",droprate =0.1,time_drop=12,
                        m=200, interim ="Yes", seed1=123, seed2=145, seed3=243){
  #simulating the enrolling time
  set.seed(seed1)
  enrol.time = -log(runif(n))/ppm
  for(i in 2:length(enrol.time))
  {  
    enrol.time[i] = enrol.time[i]+enrol.time[i-1]
  }
  #simulating the event time
  lambda0 = -log(lPFS)/pt
  set.seed(seed2)
  event.time = -log(runif(n))/lambda0
  #generating censoring status based on drop-out for no reason
  if (drop=="Yes"){
    lambda1 = -log(1-droprate)/time_drop
    set.seed(seed3)
    cens.time = -log(runif(n))/lambda1
    obs.time = ifelse(event.time <cens.time, event.time, cens.time)
    status = ifelse(event.time <cens.time, 1,0)
  }else{
    obs.time = event.time
  }
  #Calculating the study time
  study.time = obs.time+enrol.time
  #stopping the trial based on the mfol (follow-up time)
  if(fol_cut == "Yes"){
    cut = enrol.time[n]+mfol
    obs.time = ifelse(study.time <= cut ,obs.time, cut - enrol.time)
    if(drop=="Yes"){
      status = ifelse(study.time>cut, 0, status)
    }else{
      status = ifelse(study.time>cut, 0, 1)}
  }else{
    #stopping the trial based on the enrolling npts with PFS information
    if(drop=="Yes"){
      cut = event.time[which(event.time>=pt & status!=0)[npts]]
      status = ifelse(study.time>cut, 0, status)
      obs.time = ifelse(study.time <= cut ,obs.time, cut - enrol.time)
    }else{
      cut = event.time[which(event.time>=pt)[npts]]
      status = ifelse(study.time>cut, 0, 1)
      obs.time = ifelse(study.time <= cut ,obs.time, cut - enrol.time)
    }
  }
  sim_dat = data.frame(enrol.time =enrol.time, event.time=event.time,
                       study.time=study.time, obs.time=obs.time, status=status, ID = seq(1,n))
  if(interim=="Yes"){
    interim_dat=sim_dat[which(sim_dat$obs.time>=pt |sim_dat$obs.time<pt & sim_dat$study.time<=cut)[1:m],]
    return(list(full_trial=sim_dat, interim=interim_dat, cut=cut))
  }else{
    return(list(full_trial=sim_dat, cut=cut))
  }
}

###################################################################
#now we define a function to calculate the predictive probability and PET
#along with other relevent statistics
predProb_interim_PET = function(data, pt, prior.data, n_max, P_T, theta, P_L, P_U,index){
  data = data[index,]
  n_1 = nrow(data)
  prior.survival =  survfit(Surv(prior.data$obs.time, prior.data$status) ~ 1,
                            error="greenwood" , conf.type= "plain")
  smr0 = summary(prior.survival, times=pt)
  prior.mean=smr0$surv
  prior.var = (smr0$std.err)^2
  interim.survival= survfit(Surv(data$obs.time,data$status) ~ 1,
                            error="greenwood", conf.type= "plain")
  smr = summary(interim.survival, times=pt)
  interim.var = (smr$std.err)^2
  s_hat=smr$surv
  #bootstrap mean of the predictive posterior distribution
  predPost.mean = ifelse(!is.null(s_hat)&!is.null(interim.var),
                         tryCatch(((prior.var*s_hat)+(interim.var*prior.mean))/
                                    (interim.var+prior.var), error=function(e){NA}), NA)
  #bootstrap variance of the predictive posterior distribution
  predPost.var = ifelse(!is.null(interim.var),
                        tryCatch((n_1*interim.var)/n_max+(1/prior.var+1/interim.var)^(-1),
                                 error=function(e) NA), NA)
  max.var = ifelse(!is.null(interim.var),
                   tryCatch(n_1*interim.var/n_max, error=function(e) NA), NA)
  sigma_star = sqrt((1/prior.var+1/max.var)^(-1))
  delta = ((max.var+prior.var)*(theta-sigma_star*qnorm(1-P_T))-
             (max.var*prior.mean))/prior.var
  PredProb = 1-pnorm((delta-predPost.mean)/
                       sqrt(predPost.var), mean=0, sd =1)
  
  s_L = ifelse(!is.null(interim.var), tryCatch((delta-(interim.var*prior.mean/(interim.var+prior.var))-
                                                  sqrt(predPost.var)*qnorm(1-P_L))/(prior.var/(interim.var+prior.var)),error=function(e) NA) ,0) 
  
  
  s_U = ifelse(!is.null(interim.var), tryCatch((delta-(interim.var*prior.mean/(interim.var+prior.var))-
                                                  sqrt(predPost.var)*qnorm(1-P_U))/(prior.var/(interim.var+prior.var)),
                                               error=function(e) NA) ,1) 
  
  PET = ifelse(s_L==0, 0, pnorm((s_L-s_hat)/sqrt(interim.var), mean=0, sd=1))
  predProb= vector(length=7)
  predProb[1]=predPost.mean
  predProb[2]=predPost.var
  predProb[3]=delta
  predProb[4]=PredProb
  predProb[5]=s_L
  predProb[6]=s_U
  predProb[7]=PET
  return(predProb)
}

#######################################################################
#Function to perform bootstarping for predictive probability function
#here we dont need a confidence intervals
Predictive_probability_PET = function(data, statistic, R, pt, prior.data, n_max, P_T, theta, P_L, P_U, seed){
  set.seed(seed)
  results = boot(data =data, statistic=statistic, R=R, pt=pt, prior.data=prior.data,
                 n_max=n_max, P_T=P_T, theta=theta, P_L=P_L, P_U=P_U)
  bs_stat = apply(results$t,2,mean, na.rm=T)
  
  return(bs_stat)
}
###############################################################################################
#Operating characteristic function that output results for
#both beta-binomial and normal approximation
op_Char1_sim= function(simNum=1000, pt=9, n=500, n.prior=300, ppm=20, fol_cut= "Yes", mfol=10,
                       npts=200, drop ="Yes",droprate =0.1,time_drop=12, m=200, P_T, P_L, P_U,
                       theta, target, Dl, R){
  
  
  simDat0 =simTrial_lPFS1(pt=pt, n=n.prior, ppm=ppm, lPFS=theta, fol_cut=fol_cut ,mfol=mfol+1,
                          npts=npts, drop=drop ,droprate=droprate ,time_drop=time_drop, m=m,
                          interim ="No", seed1=231, seed2=9876, seed3 = 123)
  
  historical = simDat0$full_trial
  
  idx = seq(theta-Dl, target+Dl, by=0.03)
  idx = as.numeric(as.character(idx))
  petMat = matrix(ncol=5,nrow=length(idx))
  for (i in 1:length(idx)){
    PP= matrix(ncol=4, nrow=simNum)
    for (j in 1:simNum){
      ranSeed = sample(1:1000, 4, replace=F)
      simDat = simTrial_lPFS1(pt=pt, n=n, ppm=ppm, lPFS=idx[i], fol_cut=fol_cut , mfol=mfol,
                              npts=npts, drop=drop ,droprate=droprate ,time_drop=time_drop,
                              m=m, interim ="Yes", seed1=ranSeed[1], seed2=ranSeed[2], seed3=ranSeed[3])
      
      full_trial=simDat$full_trial
      interim=simDat$interim
      Predprob_norm= Predictive_probability_PET(data=interim, statistic =predProb_interim_PET, R=R,
                                                pt=pt, prior.data=historical,n_max=n, P_T=P_T,
                                                theta=theta, P_L=P_L, P_U=P_U, seed=123)
      
      y = nrow(subset(interim, obs.time==pt|obs.time>pt))
      PFS_hist = subset(historical, obs.time==pt|obs.time>pt)
      alpha_e = nrow(PFS_hist)/nrow(historical)
      beta_e = 1-alpha_e
      Predprob_betbin = predprob(y=y, n=m, nmax=n, alpha_e=alpha_e, beta_e= beta_e, p_s=theta, theta_t=P_T)
      
      bound = stopbound_pred(theta = P_L, type="futility", nmax=n,alpha_e=alpha_e ,
                             beta_e= beta_e, p_s=theta, theta_t=P_T)
      
      L_index =ifelse(length(which(bound[1]==m))==0,length(which(bound[1]<m)),which(bound[1]==m))
      L_n= bound[2][L_index,]
      PET_betbin = pbinom(L_n, m, y/nrow(interim))
      
      PP[j,1]=Predprob_norm[4]
      PP[j,2]=Predprob_norm[7]
      PP[j,3]=Predprob_betbin
      PP[j,4]= PET_betbin
    }
    PP_mean= apply(PP,2, mean)
    
    petMat[i,1]=idx[i]
    petMat[i,2]=PP_mean[1]
    petMat[i,3]=PP_mean[2]
    petMat[i,4]=PP_mean[3]
    petMat[i,5]=PP_mean[4]
  }
  return(petMat) 
}  

####################################################################################

pred1=op_Char1_sim(simNum=500, pt=9, n=50, n.prior=50, ppm=10, fol_cut= "Yes", mfol=10, npts=25, drop ="Yes",
                   droprate =0.1,time_drop=12, m=30, P_T=0.9, P_L=0.1, P_U=0.9, theta=0.4, target=0.55,
                   Dl=0.2, R=500)

pred2=op_Char1_sim(simNum=500, pt=9, n=100, n.prior=100, ppm=20, fol_cut= "Yes", mfol=20, npts=50, drop ="Yes",
                   droprate =0.1,time_drop=12, m=50, P_T=0.9, P_L=0.1, P_U=0.9, theta=0.4, target=0.55,
                   Dl=0.2, R=500)

pred4=op_Char1_sim(simNum=500, pt=9, n=300, n.prior=300, ppm=20, fol_cut= "Yes", mfol=20, npts=150, drop ="Yes",
                   droprate =0.1,time_drop=12, m=150, P_T=0.9, P_L=0.1, P_U=0.9, theta=0.4, target=0.55,
                   Dl=0.2, R=500)

##########################################################################################################

library(reshape2)
library(ggplot2)
pp1= as.data.frame(pred1[-19,])
colnames(pp1) = c("lPFS", "PPN", "PETN", "PPB", "PETB")
pp1_new = pp1[,c(1,2,4)]
pp1_reshaped = melt(pp1_new, id.vars ="lPFS")
colnames(pp1_reshaped)=c("lPFS", "Method", "predprob")


p1_PP = ggplot(data=pp1_reshaped, aes(x=lPFS, y=predprob, color=Method))+geom_point()+geom_line()+
  scale_x_continuous(name ="Landmark PFS at 9 months", breaks=seq(0.2,0.74, 0.06))+
  theme(plot.margin= unit(c(.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=10))+
  scale_y_continuous(name ="Predictive probability")+
  scale_color_manual(labels = c("Normal approximation", "Beta-binomial"), values = c("red", "blue")) +
  ggtitle("Normal approximation vs beta-binomial")+
  geom_vline(xintercept=0.4)+geom_text(aes(x=0.393, label="benchmark PFS", y=0.3),
                                       color ="black",angle=90, size=3.2)+
  geom_text(aes(x=0.53, label="Max sample size = 50", y=0.70),
            color ="black",angle=0, size=3.2)+
  geom_text(aes(x=0.53, label="Interim sample size = 30", y=0.75),
            color ="black",angle=0, size=3.2)
#-----------------------------------------------------------------------------------------
pet1 = pp1[,c(1,3,5)]
pet1_reshaped = melt(pet1, id.vars ="lPFS")
colnames(pet1_reshaped)=c("lPFS", "Method", "PET")


p1_PET = ggplot(data=pp02_reshaped, aes(x=lPFS, y=PET, color=Method))+geom_line()+
  scale_x_continuous(name ="Landmark PFS at 9 months", breaks=seq(0.2,0.74, 0.06))+
  theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=10))+
  scale_y_continuous(name ="Probability of early termination (PET)")+
  scale_color_manual(labels = c("Normal approximation", "Beta-binomial"), values = c("red", "blue")) +
  ggtitle("Normal approximation vs beta-binomial")+
  geom_vline(xintercept=0.4)+geom_text(aes(x=0.39, label="benchmark PFS", y=0.25),
                                       color ="black",angle=90, size=3.2)+
  geom_text(aes(x=0.61, label="Max sample size = 50", y=0.98),
            color ="black",angle=0, size=3.2)+
  geom_text(aes(x=0.61, label="Interim sample size = 30", y=0.93),
            color ="black",angle=0, size=3.2)

#we repeated the above plot code for the other two plots with (Max sample size of 100 and 300)
############################################################################
#Now with max sample size of 100 and different interim size
pred20=op_Char1_sim(simNum=500, pt=9, n=100, n.prior=100, ppm=20, fol_cut= "Yes", mfol=20, npts=20, drop ="Yes",
                    droprate =0.1,time_drop=12, m=20, P_T=0.9, P_L=0.1, P_U=0.9, theta=0.4, target=0.55,
                    Dl=0.2, R=500)

pred40=op_Char1_sim(simNum=500, pt=9, n=100, n.prior=100, ppm=20, fol_cut= "Yes", mfol=20, npts=40, drop ="Yes",
                    droprate =0.1,time_drop=12, m=40, P_T=0.9, P_L=0.1, P_U=0.9, theta=0.4, target=0.55,
                    Dl=0.2, R=500)


pred60=op_Char1_sim(simNum=500, pt=9, n=100, n.prior=100, ppm=20, fol_cut= "Yes", mfol=20, npts=60, drop ="Yes",
                    droprate =0.1,time_drop=12, m=60, P_T=0.9, P_L=0.1, P_U=0.9, theta=0.4, target=0.55,
                    Dl=0.2, R=500)

pred80=op_Char1_sim(simNum=500, pt=9, n=100, n.prior=100, ppm=20, fol_cut= "Yes", mfol=20, npts=80, drop ="Yes",
                    droprate =0.1,time_drop=12, m=80, P_T=0.9, P_L=0.1, P_U=0.9, theta=0.4, target=0.55,
                    Dl=0.2, R=500)

###########################################################################
pp20= as.data.frame(pet20[-19,-1])
colnames(pp20) = c("lPFS", "PPN", "PETN", "PPB", "PETB")
pp020 = pp20[,c(1,3,5)]
pp020_reshaped = melt(pp020, id.vars ="lPFS")
colnames(pp020_reshaped)=c("lPFS", "Method", "PET")


p20_PET = ggplot(data=pp020_reshaped, aes(x=lPFS, y=PET, color=Method))+geom_line()+
  scale_x_continuous(name ="Landmark PFS at 9 months", breaks=seq(0.2,0.74, 0.06))+
  theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=10))+
  scale_y_continuous(name ="Probability of early termination (PET)")+
  scale_color_manual(labels = c("Normal approximation", "Beta-binomial"), values = c("red", "blue")) +
  ggtitle("Normal approximation vs beta-binomial")+
  geom_vline(xintercept=0.4)+geom_text(aes(x=0.39, label="benchmark PFS", y=0.25),
                                       color ="black",angle=90, size=3.2)+
  geom_text(aes(x=0.61, label="Max sample size = 100", y=0.98),
            color ="black",angle=0, size=3.2)+
  geom_text(aes(x=0.61, label="Interim sample size = 20", y=0.93),
            color ="black",angle=0, size=3.2)

#----------------------------------------------
pp40= as.data.frame(pet40[-19,-1])
colnames(pp40) = c("lPFS", "PPN", "PETN", "PPB", "PETB")
pp040 = pp40[,c(1,3,5)]
pp040_reshaped = melt(pp040, id.vars ="lPFS")
colnames(pp040_reshaped)=c("lPFS", "Method", "PET")


p40_PET = ggplot(data=pp040_reshaped, aes(x=lPFS, y=PET, color=Method))+geom_line()+
  scale_x_continuous(name ="Landmark PFS at 9 months", breaks=seq(0.2,0.74, 0.06))+
  theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=10))+
  scale_y_continuous(name ="Probability of early termination (PET)")+
  scale_color_manual(labels = c("Normal approximation", "Beta-binomial"), values = c("red", "blue")) +
  ggtitle("Normal approximation vs beta-binomial")+
  geom_vline(xintercept=0.4)+geom_text(aes(x=0.39, label="benchmark PFS", y=0.25),
                                       color ="black",angle=90, size=3.2)+
  geom_text(aes(x=0.61, label="Max sample size = 100", y=0.98),
            color ="black",angle=0, size=3.2)+
  geom_text(aes(x=0.61, label="Interim sample size = 40", y=0.93),
            color ="black",angle=0, size=3.2)

#-----------------------------------------
pp60= as.data.frame(pet60[-19,-1])
colnames(pp60) = c("lPFS", "PPN", "PETN", "PPB", "PETB")
pp060 = pp60[,c(1,3,5)]
pp060_reshaped = melt(pp060, id.vars ="lPFS")
colnames(pp060_reshaped)=c("lPFS", "Method", "PET")


p60_PET = ggplot(data=pp060_reshaped, aes(x=lPFS, y=PET, color=Method))+geom_line()+
  scale_x_continuous(name ="Landmark PFS at 9 months", breaks=seq(0.2,0.74, 0.06))+
  theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=10))+
  scale_y_continuous(name ="Probability of early termination (PET)")+
  scale_color_manual(labels = c("Normal approximation", "Beta-binomial"), values = c("red", "blue")) +
  ggtitle("Normal approximation vs beta-binomial")+
  geom_vline(xintercept=0.4)+geom_text(aes(x=0.39, label="benchmark PFS", y=0.25),
                                       color ="black",angle=90, size=3.2)+
  geom_text(aes(x=0.61, label="Max sample size = 100", y=0.98),
            color ="black",angle=0, size=3.2)+
  geom_text(aes(x=0.61, label="Interim sample size = 60", y=0.93),
            color ="black",angle=0, size=3.2)

#-------------------------------------------------------

pp80= as.data.frame(pet80[-19,-1])
colnames(pp80) = c("lPFS", "PPN", "PETN", "PPB", "PETB")
pp080 = pp80[,c(1,3,5)]
pp080_reshaped = melt(pp080, id.vars ="lPFS")
colnames(pp080_reshaped)=c("lPFS", "Method", "PET")


p80_PET = ggplot(data=pp080_reshaped, aes(x=lPFS, y=PET, color=Method))+geom_line()+
  scale_x_continuous(name ="Landmark PFS at 9 months", breaks=seq(0.2,0.74, 0.06))+
  theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=10))+
  scale_y_continuous(name ="Probability of early termination (PET)")+
  scale_color_manual(labels = c("Normal approximation", "Beta-binomial"), values = c("red", "blue")) +
  ggtitle("Normal approximation vs beta-binomial")+
  geom_vline(xintercept=0.4)+geom_text(aes(x=0.39, label="benchmark PFS", y=0.25),
                                       color ="black",angle=90, size=3.2)+
  geom_text(aes(x=0.61, label="Max sample size = 100", y=0.98),
            color ="black",angle=0, size=3.2)+
  geom_text(aes(x=0.61, label="Interim sample size = 80", y=0.93),
            color ="black",angle=0, size=3.2)


#################################################################
###################################################################
#now we define a function to calculate the predictive probability and PET
#along with other relevent statistics without taking into account the historical data
predProb_interim_PET_no_hist = function(data, pt, n_max, P_T, theta, P_L, P_U,index){
  data = data[index,]
  n_1 = nrow(data)
  prior.mean=0
  prior.var = 100
  interim.survival= survfit(Surv(data$obs.time,data$status) ~ 1,
                            error="greenwood", conf.type= "plain")
  smr = summary(interim.survival, times=pt)
  interim.var = (smr$std.err)^2
  s_hat=smr$surv
  #bootstrap mean of the predictive posterior distribution
  predPost.mean = ifelse(!is.null(s_hat)&!is.null(interim.var),
                         tryCatch(((prior.var*s_hat)+(interim.var*prior.mean))/
                                    (interim.var+prior.var), error=function(e){NA}), NA)
  #bootstrap variance of the predictive posterior distribution
  predPost.var = ifelse(!is.null(interim.var),
                        tryCatch((n_1*interim.var)/n_max+(1/prior.var+1/interim.var)^(-1),
                                 error=function(e) NA), NA)
  max.var = ifelse(!is.null(interim.var),
                   tryCatch(n_1*interim.var/n_max, error=function(e) NA), NA)
  sigma_star = sqrt((1/prior.var+1/max.var)^(-1))
  delta = ((max.var+prior.var)*(theta-sigma_star*qnorm(1-P_T))-
             (max.var*prior.mean))/prior.var
  PredProb = 1-pnorm((delta-predPost.mean)/
                       sqrt(predPost.var), mean=0, sd =1)
  
  s_L = ifelse(!is.null(interim.var), tryCatch((delta-(interim.var*prior.mean/(interim.var+prior.var))-
                                                  sqrt(predPost.var)*qnorm(1-P_L))/(prior.var/(interim.var+prior.var)),error=function(e) NA) ,0) 
  
  
  s_U = ifelse(!is.null(interim.var), tryCatch((delta-(interim.var*prior.mean/(interim.var+prior.var))-
                                                  sqrt(predPost.var)*qnorm(1-P_U))/(prior.var/(interim.var+prior.var)),
                                               error=function(e) NA) ,1) 
  
  PET = ifelse(s_L==0, 0, pnorm((s_L-s_hat)/sqrt(interim.var), mean=0, sd=1))
  predProb= vector(length=7)
  predProb[1]=predPost.mean
  predProb[2]=predPost.var
  predProb[3]=delta
  predProb[4]=PredProb
  predProb[5]=s_L
  predProb[6]=s_U
  predProb[7]=PET
  return(predProb)
}

#######################################################################
#Function to perform bootstarping for predictive probability function
#here we dont need a confidence intervals
Predictive_probability_PET_no_hist = function(data, statistic, R, pt, n_max, P_T, theta, P_L, P_U, seed){
  set.seed(seed)
  results = boot(data =data, statistic=statistic, R=R, pt=pt, 
                 n_max=n_max, P_T=P_T, theta=theta, P_L=P_L, P_U=P_U)
  bs_stat = apply(results$t,2,mean, na.rm=T)
  
  return(bs_stat)
}
###############################################################################################
#Operating characteristic function that output results for
#both beta-binomial and normal approximation without historical data
op_Char1_sim_no_hist=function(simNum=1000, pt=9, n=500, ppm=20, fol_cut= "Yes", mfol=10,
                              npts=200, drop ="Yes",droprate =0.1,time_drop=12, m=200, P_T, P_L, P_U,
                              theta, target, Dl, R){
  idx = seq(theta-Dl, target+Dl, by=0.03)
  idx = as.numeric(as.character(idx))
  petMat = matrix(ncol=5,nrow=length(idx))
  for (i in 1:length(idx)){
    PP= matrix(ncol=4, nrow=simNum)
    for (j in 1:simNum){
      ranSeed = sample(1:1000, 4, replace=F)
      simDat = simTrial_lPFS1(pt=pt, n=n, ppm=ppm, lPFS=idx[i], fol_cut=fol_cut , mfol=mfol,
                              npts=npts, drop=drop ,droprate=droprate ,time_drop=time_drop,
                              m=m, interim ="Yes", seed1=ranSeed[1], seed2=ranSeed[2], seed3=ranSeed[3])
      
      full_trial=simDat$full_trial
      interim=simDat$interim
      Predprob_norm= Predictive_probability_PET_no_hist(data=interim, statistic =predProb_interim_PET_no_hist, R=R,
                                                        pt=pt,n_max=n, P_T=P_T, theta=theta, P_L=P_L, P_U=P_U, seed=123)
      
      y = nrow(subset(interim, obs.time==pt|obs.time>pt))
      alpha_e = 0.5
      beta_e = 0.5
      Predprob_betbin = predprob(y=y, n=m, nmax=n, alpha_e=alpha_e, beta_e= beta_e, p_s=theta, theta_t=P_T)
      
      bound = stopbound_pred(theta = P_L, type="futility", nmax=n,alpha_e=alpha_e ,
                             beta_e= beta_e, p_s=theta, theta_t=P_T)
      
      L_index =ifelse(length(which(bound[1]==m))==0,length(which(bound[1]<m)),which(bound[1]==m))
      L_n= bound[2][L_index,]
      PET_betbin = pbinom(L_n, m, y/nrow(interim))
      
      PP[j,1]=Predprob_norm[4]
      PP[j,2]=Predprob_norm[7]
      PP[j,3]=Predprob_betbin
      PP[j,4]= PET_betbin
    }
    PP_mean= apply(PP,2, mean)
    
    petMat[i,1]=idx[i]
    petMat[i,2]=PP_mean[1]
    petMat[i,3]=PP_mean[2]
    petMat[i,4]=PP_mean[3]
    petMat[i,5]=PP_mean[4]
  }
  return(petMat) 
}  


###########################
pred20_no_hist=op_Char1_sim_no_hist(simNum=500, pt=9, n=100, ppm=20, fol_cut= "Yes", mfol=20, npts=20, drop ="Yes",
                                    droprate =0.1,time_drop=12, m=20, P_T=0.9, P_L=0.1, P_U=0.9, theta=0.4, target=0.55,
                                    Dl=0.2, R=500)
pred20_no_hist_dat= as.data.frame(pred20_no_hist)
pred20_no_hist_dat = write.csv(pred20_no_hist_dat,"~/Desktop//pred20_no_hist_dat.csv")
pet20_no_hist= read.csv("~/Desktop//pred20_no_hist_dat.csv")
#----------------------------------------------------------------------------------------
pred40_no_hist=op_Char1_sim_no_hist(simNum=500, pt=9, n=100, ppm=20, fol_cut= "Yes", mfol=20, npts=40, drop ="Yes",
                                    droprate =0.1,time_drop=12, m=40, P_T=0.9, P_L=0.1, P_U=0.9, theta=0.4, target=0.55,
                                    Dl=0.2, R=500)
pred40_no_hist_dat= as.data.frame(pred40_no_hist)
pred40_no_hist_dat = write.csv(pred40_no_hist_dat,"~/Desktop//pred40_no_hist_dat.csv")
pet40_no_hist= read.csv("~/Desktop//pred40_no_hist_dat.csv")
#---------------------------------------------------------------------------------------------
pred60_no_hist=op_Char1_sim_no_hist(simNum=500, pt=9, n=100, ppm=20, fol_cut= "Yes", mfol=20, npts=60, drop ="Yes",
                                    droprate =0.1,time_drop=12, m=60, P_T=0.9, P_L=0.1, P_U=0.9, theta=0.4, target=0.55,
                                    Dl=0.2, R=500)
pred60_no_hist_dat= as.data.frame(pred60_no_hist)
pred60_no_hist_dat = write.csv(pred60_no_hist_dat,"~/Desktop//pred60_no_hist_dat.csv")
pet60_no_hist= read.csv("~/Desktop//pred60_no_hist_dat.csv")
#-------------------------------------------------------------------------------------------
pred80_no_hist=op_Char1_sim_no_hist(simNum=500, pt=9, n=100, ppm=20, fol_cut= "Yes", mfol=20, npts=80, drop ="Yes",
                                    droprate =0.1,time_drop=12, m=80, P_T=0.9, P_L=0.1, P_U=0.9, theta=0.4, target=0.55,
                                    Dl=0.2, R=500)
pred80_no_hist_dat= as.data.frame(pred80_no_hist)
pred80_no_hist_dat = write.csv(pred80_no_hist_dat,"~/Desktop//pred80_no_hist_dat.csv")
pet80_no_hist= read.csv("~/Desktop//pred80_no_hist_dat.csv")
##########################################################################
pp20_no_hist= as.data.frame(pet20_no_hist[-19,-1])
colnames(pp20_no_hist) = c("lPFS", "PPN", "PETN", "PPB", "PETB")
pp020_no_hist = pp20_no_hist[,c(1,3,5)]
pp020_no_hist_reshaped = melt(pp020_no_hist, id.vars ="lPFS")
colnames(pp020_no_hist_reshaped)=c("lPFS", "Method", "PET")


p20_no_hist_PET = ggplot(data=pp020_no_hist_reshaped, aes(x=lPFS, y=PET, color=Method))+geom_line()+
  scale_x_continuous(name ="Landmark PFS at 9 months", breaks=seq(0.2,0.74, 0.06))+
  theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=10))+
  scale_y_continuous(name ="Probability of early termination (PET)")+
  scale_color_manual(labels = c("Normal approximation", "Beta-binomial"), values = c("red", "blue")) +
  ggtitle("Normal approximation vs beta-binomial")+
  geom_vline(xintercept=0.4)+geom_text(aes(x=0.39, label="benchmark PFS", y=0.25),
                                       color ="black",angle=90, size=3.2)+
  geom_text(aes(x=0.61, label="Max sample size = 100", y=0.98),
            color ="black",angle=0, size=3.2)+
  geom_text(aes(x=0.61, label="Interim sample size = 20", y=0.93),
            color ="black",angle=0, size=3.2)

#-----------------------------------------
pp40_no_hist= as.data.frame(pet40_no_hist[-19,-1])
colnames(pp40_no_hist) = c("lPFS", "PPN", "PETN", "PPB", "PETB")
pp040_no_hist = pp40_no_hist[,c(1,3,5)]
pp040_no_hist_reshaped = melt(pp040_no_hist, id.vars ="lPFS")
colnames(pp040_no_hist_reshaped)=c("lPFS", "Method", "PET")


p40_no_hist_PET = ggplot(data=pp040_no_hist_reshaped, aes(x=lPFS, y=PET, color=Method))+geom_line()+
  scale_x_continuous(name ="Landmark PFS at 9 months", breaks=seq(0.2,0.74, 0.06))+
  theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=10))+
  scale_y_continuous(name ="Probability of early termination (PET)")+
  scale_color_manual(labels = c("Normal approximation", "Beta-binomial"), values = c("red", "blue")) +
  ggtitle("Normal approximation vs beta-binomial")+
  geom_vline(xintercept=0.4)+geom_text(aes(x=0.39, label="benchmark PFS", y=0.25),
                                       color ="black",angle=90, size=3.2)+
  geom_text(aes(x=0.61, label="Max sample size = 100", y=0.98),
            color ="black",angle=0, size=3.2)+
  geom_text(aes(x=0.61, label="Interim sample size = 40", y=0.93),
            color ="black",angle=0, size=3.2)

#--------------------------------------------------------
pp60_no_hist= as.data.frame(pet60_no_hist[-19,-1])
colnames(pp60_no_hist) = c("lPFS", "PPN", "PETN", "PPB", "PETB")
pp060_no_hist = pp60_no_hist[,c(1,3,5)]
pp060_no_hist_reshaped = melt(pp060_no_hist, id.vars ="lPFS")
colnames(pp060_no_hist_reshaped)=c("lPFS", "Method", "PET")


p60_no_hist_PET = ggplot(data=pp060_no_hist_reshaped, aes(x=lPFS, y=PET, color=Method))+geom_line()+
  scale_x_continuous(name ="Landmark PFS at 9 months", breaks=seq(0.2,0.74, 0.06))+
  theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=10))+
  scale_y_continuous(name ="Probability of early termination (PET)")+
  scale_color_manual(labels = c("Normal approximation", "Beta-binomial"), values = c("red", "blue")) +
  ggtitle("Normal approximation vs beta-binomial")+
  geom_vline(xintercept=0.4)+geom_text(aes(x=0.39, label="benchmark PFS", y=0.25),
                                       color ="black",angle=90, size=3.2)+
  geom_text(aes(x=0.61, label="Max sample size = 100", y=0.98),
            color ="black",angle=0, size=3.2)+
  geom_text(aes(x=0.61, label="Interim sample size = 60", y=0.93),
            color ="black",angle=0, size=3.2)

######################################################
pp80_no_hist= as.data.frame(pet80_no_hist[-19,-1])
colnames(pp80_no_hist) = c("lPFS", "PPN", "PETN", "PPB", "PETB")
pp080_no_hist = pp80_no_hist[,c(1,3,5)]
pp080_no_hist_reshaped = melt(pp080_no_hist, id.vars ="lPFS")
colnames(pp080_no_hist_reshaped)=c("lPFS", "Method", "PET")


p80_no_hist_PET = ggplot(data=pp080_no_hist_reshaped, aes(x=lPFS, y=PET, color=Method))+geom_line()+
  scale_x_continuous(name ="Landmark PFS at 9 months", breaks=seq(0.2,0.74, 0.06))+
  theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=10))+
  scale_y_continuous(name ="Probability of early termination (PET)")+
  scale_color_manual(labels = c("Normal approximation", "Beta-binomial"), values = c("red", "blue")) +
  ggtitle("Normal approximation vs beta-binomial")+
  geom_vline(xintercept=0.4)+geom_text(aes(x=0.39, label="benchmark PFS", y=0.25),
                                       color ="black",angle=90, size=3.2)+
  geom_text(aes(x=0.61, label="Max sample size = 100", y=0.98),
            color ="black",angle=0, size=3.2)+
  geom_text(aes(x=0.61, label="Interim sample size = 80", y=0.93),
            color ="black",angle=0, size=3.2)




