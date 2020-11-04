library(survival)
library(boot)
library(survminer)
library(ph2bayes)
library(ggplot2)
library(gridExtra)
ggsurv = function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                  cens.col = 'purple', lty.est = 1, lty.ci = 3,
                  cens.shape = 3, back.white = F, xlab = 'Time in months',
                  ylab = 'Proportion not progressing', main = ""){
  
  library(ggplot2)
  strata <- ifelse(is.null(s$strata) ==T, 1, length(s$strata))
  stopifnot(length(surv.col) == 1 | length(surv.col) == strata)
  stopifnot(length(lty.est) == 1 | length(lty.est) == strata)
  
  ggsurv.s <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'purple', lty.est = 1, lty.ci = 3,
                       cens.shape = 3, back.white = F, xlab = 'Time in months',
                       ylab = 'Proportion not progressing', main = ''){
    
    dat <- data.frame(time = c(0, s$time),
                      surv = c(1, s$surv),
                      up = c(1, s$upper),
                      low = c(1, s$lower),
                      cens = c(0, s$n.censor))
    dat.cens <- subset(dat, cens != 0)
    
    col <- ifelse(surv.col == 'gg.def', 'black', surv.col)
    
    pl <- ggplot(dat, aes(x = time, y = surv)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(col = col, lty = lty.est)
    
    pl <- if(CI == T | CI == 'def') {
      pl + geom_step(aes(y = up), color = col, lty = lty.ci) +
        geom_step(aes(y = low), color = col, lty = lty.ci)
    } else (pl)
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  
  ggsurv.m <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 3,
                       cens.shape = 3, back.white = F, xlab = 'Time in months',
                       ylab = 'Proportion not progressing', main = '') {
    n <- s$strata
    
    groups <- factor(unlist(strsplit(names
                                     (s$strata), '='))[seq(2, 2*strata, by = 2)])
    gr.name <-  unlist(strsplit(names(s$strata), '='))[1]
    gr.df <- vector('list', strata)
    ind <- vector('list', strata)
    n.ind <- c(0,n); n.ind <- cumsum(n.ind)
    for(i in 1:strata) ind[[i]] <- (n.ind[i]+1):n.ind[i+1]
    
    for(i in 1:strata){
      gr.df[[i]] <- data.frame(
        time = c(0, s$time[ ind[[i]] ]),
        surv = c(1, s$surv[ ind[[i]] ]),
        up = c(1, s$upper[ ind[[i]] ]),
        low = c(1, s$lower[ ind[[i]] ]),
        cens = c(0, s$n.censor[ ind[[i]] ]),
        group = rep(groups[i], n[i] + 1))
    }
    
    dat <- do.call(rbind, gr.df)
    dat.cens <- subset(dat, cens != 0)
    
    pl <- ggplot(dat, aes(x = time, y = surv, group = group)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(aes(col = group, lty = group))
    
    col <- if(length(surv.col == 1)){
      scale_colour_manual(name = gr.name, values = rep(surv.col, strata))
    } else{
      scale_colour_manual(name = gr.name, values = surv.col)
    }
    
    pl <- if(surv.col[1] != 'gg.def'){
      pl + col
    } else {pl + scale_colour_discrete(name = gr.name)}
    
    line <- if(length(lty.est) == 1){
      scale_linetype_manual(name = gr.name, values = rep(lty.est, strata))
    } else {scale_linetype_manual(name = gr.name, values = lty.est)}
    
    pl <- pl + line
    
    pl <- if(CI == T) {
      if(length(surv.col) > 1 && length(lty.est) > 1){
        stop('Either surv.col or lty.est should be of length 1 in order
             to plot 95% CI with multiple strata')
      }else if((length(surv.col) > 1 | surv.col == 'gg.def')[1]){
        pl + geom_step(aes(y = up, color = group), lty = lty.ci) +
          geom_step(aes(y = low, color = group), lty = lty.ci)
      } else{pl +  geom_step(aes(y = up, lty = group), col = surv.col) +
          geom_step(aes(y = low,lty = group), col = surv.col)}
    } else {pl}
    
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  pl <- if(strata == 1) {ggsurv.s(s, CI , plot.cens, surv.col ,
                                  cens.col, lty.est, lty.ci,
                                  cens.shape, back.white, xlab,
                                  ylab, main)
  } else {ggsurv.m(s, CI, plot.cens, surv.col ,
                   cens.col, lty.est, lty.ci,
                   cens.shape, back.white, xlab,
                   ylab, main)}
  pl
}


#function to simulate a clinical trial
simTrial=function(pt=9, n=500, ppm=20, mPFS=6, fol_cut = "Yes", mfol =20, npts=200, 
                  drop ="Yes",droprate =0.1,time_drop=12, m=200, interim ="Yes", 
                  seed1=123, seed2=145, seed3=243){
  #simulating the enrolling time
  set.seed(seed1)
  enrol.time = -log(runif(n))/ppm
  for(i in 2:length(enrol.time))
  {   
    enrol.time[i] = enrol.time[i]+enrol.time[i-1]
  }
  #simulating the event time
  lambda0 = log(2)/mPFS
  set.seed(seed2)
  event.time = -log(runif(n))/lambda0
  #generating censoring status based on drop-out for no reason
  if (drop=="Yes"){
    set.seed(seed3)
    lambda1 = -log(1-droprate)/time_drop
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
    interim_dat=sim_dat[1:m,]
    cut2=interim_dat$enrol.time[m]+mfol
    rest = sim_dat[m+1:n,]
    rest_interim = subset(rest, enrol.time<cut2| enrol.time==cut2)
    interim = rbind(interim_dat, rest_interim)
    interim$status = ifelse(interim$study.time>cut2, 0, interim$status)
    
    return(list(full_trial=sim_dat, interim=interim, cut1=cut, cut2=cut2))
  }else{
    return(list(full_trial=sim_dat, cut=cut)) 
  }
}

#try the normal approximation and beta-binomial models with 200 pts and 60
#pts as well with half data used for interim analysis
#simulate full trial with 200 pts and use 100 pts as interim data
simDat = simTrial(pt=6, n=200,ppm=10, mPFS=12, fol_cut = "Yes", mfol =6, npts=75, 
                  drop ="Yes",droprate =0.2, time_drop=12, m=75, interim="Yes", 
                  seed1=123, seed2=145, seed3 = 243)

full_trial=simDat$full_trial
cut1 = simDat$cut1
cut2 = simDat$cut2
interim=simDat$interim
nrow(full_trial)
nrow(interim)
table(interim$status)
tail(full_trial)
head(interim)

pt=6
#K-M plot for full data
full_trial.survival= survfit(Surv(full_trial$obs.time,full_trial$status) ~ 1, 
                             error="greenwood", conf.type= "plain")
smr = summary(full_trial.survival, times=pt) 
survminer::surv_median(full_trial.survival)
table(full_trial$status)
p1 = ggsurv(full_trial.survival)
pl = p1+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,30, 2), expand = c(0, 0))+scale_y_continuous(breaks=seq(0,1, 0.1),expand = c(0, 0))+ggtitle("Full data")+theme(plot.margin = unit(c(1.1,1.1,1.1,1.1),"cm"))+
  geom_segment(aes(x=6, y=0, xend=6, yend=0.69), linetype="dashed")+geom_segment(aes(x=0, y=0.69, xend=6, yend=0.69), linetype="dashed")+
  geom_text(aes(x=12.2, y=0.69, label="PFS at 6 months: 0.69 \n 95% CI: (0.62, 0.75)"), size=3.2)

pl
#----------------------------------------------------------
#K-M plot for interim data
interim.survival= survfit(Surv(interim$obs.time,interim$status) ~ 1, 
                          error="greenwood", conf.type= "plain")

smr1 = summary(interim.survival, times=pt) 
survminer::surv_median(interim.survival)
table(interim$status)
p11 = ggsurv(interim.survival)
pl1 = p11+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,40, 2), expand = c(0, 0))+scale_y_continuous(breaks=seq(0,1, 0.1),expand = c(0, 0))+ggtitle("Interim data")+theme(plot.margin = unit(c(1.1,1.1,1.1,1.1),"cm"))+
  geom_segment(aes(x=6, y=0, xend=6, yend=0.75), linetype="dashed")+geom_segment(aes(x=0, y=0.75, xend=6, yend=0.75), linetype="dashed")+
  geom_text(aes(x=12.2, y=0.75, label="PFS at 6 months: 0.75 \n 95% CI: (0.68, 0.83)"), size=3.2)
pl1

#---------------------------------------------------------------
simDat2 = simTrial(pt=6, n=200, ppm=10, mPFS=9, fol_cut = "Yes", mfol =6, npts=75, 
                   drop ="Yes",droprate =0.2,time_drop=12, m=75, interim="No", 
                   seed1=231, seed2=9876, seed3 = 123)

historical = simDat2$full_trial
cut1_historical = simDat2$cut
table(historical$status)

#K-M plot for historical data
historical.survival= survfit(Surv(historical$obs.time,historical$status) ~ 1, 
                             error="greenwood", conf.type= "plain")
smr2 = summary(historical.survival, times=pt) 
survminer::surv_median(historical.survival)
table(historical$status)

p12 = ggsurv(historical.survival)
pl2 = p12+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,40, 2), expand = c(0, 0))+scale_y_continuous(breaks=seq(0,1, 0.1),expand = c(0, 0))+ggtitle("Historical data")+theme(plot.margin = unit(c(1.1,1.1,1.1,1.1),"cm"))+
  geom_segment(aes(x=6, y=0, xend=6, yend=0.65), linetype="dashed")+geom_segment(aes(x=0, y=0.65, xend=6, yend=0.65), linetype="dashed")+
  geom_text(aes(x=12.5, y=0.65, label="PFS at 6 months: 0.65 \n 95% CI: (0.59, 0.72)"), size=3.2)

pl2

grid.arrange(pl,pl1, pl2, nrow = 1)

########################################################################
#proportion of responses in all simulated data sets
#1) proportion of responses in the historical data
pt=6
PFS9_hist = subset(historical, obs.time==pt|obs.time>pt)
#Tproportion of responses
alpha_e = nrow(PFS9_hist)/nrow(historical)
#2) proportion of responses in the interim data
PFS9_interim = subset(interim, obs.time==pt|obs.time>pt)
#Tproportion of responses
alpha_e = nrow(PFS9_interim)/nrow(interim)
#2) proportion of responses in the full data
PFS9_full = subset(full_trial, obs.time==pt|obs.time>pt)
#Tproportion of responses
alpha_e = nrow(PFS9_full)/nrow(full_trial)
########################################################################
#calculating PPOS using the proposed method (without bootstraping)
predProb_interim = function(data, pt, prior.data, n_max, P_T, theta, index){
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
  predProb= vector(length=4)
  predProb[1]=predPost.mean
  predProb[2]=predPost.var
  predProb[3]=delta
  predProb[4]=PredProb
  
  return(predProb) 
}

#final function to calculate PPOS with bootstraping
Predictive_probability = function(data, statistic, R, pt, prior.data, n_max, P_T, theta, seed){
  set.seed(seed)
  results = boot(data =data, statistic=statistic, R=R, pt=pt, prior.data=prior.data, 
                 n_max=n_max, P_T=P_T, theta=theta)
  bs_stat = apply(results$t,2,mean, na.rm=T)
  bootCI_percent = matrix(ncol=2, nrow=4)
  for(i in 1:4){
    bootCI_percent[i,1] =boot.ci(results, index=i, type=c("perc"))$percent[4] 
    bootCI_percent[i,2] =boot.ci(results, index=i, type=c("perc"))$percent[5] 
  }
  str=paste0("    Statistics                    95% Percentile CI\nPosterior mean         =",
             round(bs_stat[1],3),"     (",
             round(bootCI_percent[1,1],3),", ",round(bootCI_percent[1,2],3),")\nPosterior variance     =",
             round(bs_stat[2],3), "     (", round(bootCI_percent[2,1],3), ","
             ,round(bootCI_percent[2,2],3),")\nSuccess Threshold      =", round(bs_stat[3],3),
             "     (" ,round(bootCI_percent[3,1],3), ","  ,round(bootCI_percent[3,2],3),
             ")\nPredictive Probability =", round(bs_stat[4],3), "     (", 
             round(bootCI_percent[4,1],3),"," , round(bootCI_percent[4,2],3),")")
  
  cat(str)
  
}

#######################################################

Predictive_probability(data=interim, statistic =predProb_interim, R=1000, pt=6, prior.data=historical,
                       n_max=200, P_T=0.9, theta=0.65, seed=123)


Predictive_probability(data=interim, statistic =predProb_interim, R=1000, pt=6, prior.data=historical,
                       n_max=200, P_T=0.8, theta=0.65, seed=1)

Predictive_probability(data=interim, statistic =predProb_interim, R=1000, pt=6, prior.data=historical,
                       n_max=200, P_T=0.7, theta=0.65, seed=12356)

#P_T=0.6, theta =0.56
Predictive_probability(data=interim, statistic =predProb_interim, R=1000, pt=6, prior.data=historical,
                       n_max=200, P_T=0.6, theta=0.65, seed=123)

Predictive_probability(data=interim, statistic =predProb_interim, R=1000, pt=6, prior.data=historical,
                       n_max=200, P_T=0.5, theta=0.65, seed=12365)
#------------------------------------------------------------------------
library(ph2bayes)
PFS9 = subset(interim, obs.time==pt|obs.time>pt)
y = nrow(PFS9)
#number of patient at interim
n = nrow(interim)
#maximum number of patient (interim+future)
nmax= nrow(full_trial)
PFS9_hist = subset(historical, obs.time==pt|obs.time>pt)
#The hyperparameter(shape1) of the beta prior: proportion of responses in historical data
alpha_e = nrow(PFS9_hist)/nrow(historical)
#The hyperparameter(shape2) of the beta prior: proportion of non-responders in historical data
beta_e = 1-alpha_e
#p_s:the efficacy threshold from the standard treatment  
#theta_t= The prespecified target probability

predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.9) 
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.8) 
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.7) 
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.6)
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.5)
####################################################################################################
simDat = simTrial(pt=6, n=200,ppm=10, mPFS=12, fol_cut = "Yes", mfol =6, npts=75, 
                  drop ="Yes",droprate =0.1, time_drop=12, m=75, interim="Yes", 
                  seed1=123, seed2=145, seed3 = 243)

full_trial=simDat$full_trial
cut1 = simDat$cut1
cut2 = simDat$cut2
interim=simDat$interim
nrow(full_trial)
nrow(interim)
table(interim$status)
tail(full_trial)
head(interim)

pt=6
#K-M plot for full data
full_trial.survival= survfit(Surv(full_trial$obs.time,full_trial$status) ~ 1, 
                             error="greenwood", conf.type= "plain")
smr = summary(full_trial.survival, times=pt) 
survminer::surv_median(full_trial.survival)
table(full_trial$status)

p1 = ggsurv(full_trial.survival)
pl = p1+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,30, 2), expand = c(0, 0))+scale_y_continuous(breaks=seq(0,1, 0.1),expand = c(0, 0))+ggtitle("Full data")+theme(plot.margin = unit(c(1.1,1.1,1.1,1.1),"cm"))+
  geom_segment(aes(x=6, y=0, xend=6, yend=0.69), linetype="dashed")+geom_segment(aes(x=0, y=0.69, xend=6, yend=0.69), linetype="dashed")+
  geom_text(aes(x=12.2, y=0.69, label="PFS at 6 months: 0.69 \n 95% CI: (0.62, 0.75)"), size=3.2)

pl
#----------------------------------------------------------
#K-M plot for interim data
interim.survival= survfit(Surv(interim$obs.time,interim$status) ~ 1, 
                          error="greenwood", conf.type= "plain")

smr1 = summary(interim.survival, times=pt) 

survminer::surv_median(interim.survival)
table(interim$status)
p11 = ggsurv(interim.survival)
pl1 = p11+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,40, 2), expand = c(0, 0))+scale_y_continuous(breaks=seq(0,1, 0.1),expand = c(0, 0))+ggtitle("Interim data")+theme(plot.margin = unit(c(1.1,1.1,1.1,1.1),"cm"))+
  geom_segment(aes(x=6, y=0, xend=6, yend=0.75), linetype="dashed")+geom_segment(aes(x=0, y=0.75, xend=6, yend=0.75), linetype="dashed")+
  geom_text(aes(x=12.2, y=0.75, label="PFS at 6 months: 0.75 \n 95% CI: (0.68, 0.83)"), size=3.2)
pl1

#---------------------------------------------------------------

simDat2 = simTrial(pt=6, n=200, ppm=10, mPFS=9, fol_cut = "Yes", mfol =6, npts=75, 
                   drop ="Yes",droprate =0.1,time_drop=12, m=75, interim="No", 
                   seed1=8642, seed2=9753, seed3 = 123745)

historical = simDat2$full_trial
cut1_historical = simDat2$cut

table(historical$status)

#K-M plot for historical data
historical.survival= survfit(Surv(historical$obs.time,historical$status) ~ 1, 
                             error="greenwood", conf.type= "plain")
smr2 = summary(historical.survival, times=pt) 

survminer::surv_median(historical.survival)

table(historical$status)

p12 = ggsurv(historical.survival)
pl2 = p12+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,40, 2), expand = c(0, 0))+scale_y_continuous(breaks=seq(0,1, 0.1),expand = c(0, 0))+ggtitle("Historical data")+theme(plot.margin = unit(c(1.1,1.1,1.1,1.1),"cm"))+
  geom_segment(aes(x=6, y=0, xend=6, yend=0.65), linetype="dashed")+geom_segment(aes(x=0, y=0.65, xend=6, yend=0.65), linetype="dashed")+
  geom_text(aes(x=12.5, y=0.65, label="PFS at 6 months: 0.65 \n 95% CI: (0.58, 0.72)"), size=3.2)

pl2

grid.arrange(pl,pl1, pl2, nrow = 1)

########################################################################
#proportion of responses in all simulated data sets
#1) proportion of responses in the historical data
pt=6
PFS9_hist = subset(historical, obs.time==pt|obs.time>pt)
#Tproportion of responses
alpha_e = nrow(PFS9_hist)/nrow(historical)
#2) proportion of responses in the interim data
PFS9_interim = subset(interim, obs.time==pt|obs.time>pt)
#Tproportion of responses
alpha_e = nrow(PFS9_interim)/nrow(interim)
#2) proportion of responses in the full data
PFS9_full = subset(full_trial, obs.time==pt|obs.time>pt)
#Tproportion of responses
alpha_e = nrow(PFS9_full)/nrow(full_trial)
########################################################################

Predictive_probability(data=interim, statistic =predProb_interim, R=1000, pt=6, prior.data=historical,
                       n_max=200, P_T=0.9, theta=0.65, seed=123)
Predictive_probability(data=interim, statistic =predProb_interim, R=1000, pt=6, prior.data=historical,
                       n_max=200, P_T=0.8, theta=0.65, seed=123)

Predictive_probability(data=interim, statistic =predProb_interim, R=1000, pt=6, prior.data=historical,
                       n_max=200, P_T=0.7, theta=0.65, seed=12356)

Predictive_probability(data=interim, statistic =predProb_interim, R=1000, pt=6, prior.data=historical,
                       n_max=200, P_T=0.6, theta=0.65, seed=123)

Predictive_probability(data=interim, statistic =predProb_interim, R=1000, pt=6, prior.data=historical,
                       n_max=200, P_T=0.5, theta=0.65, seed=12365)

#------------------------------------------------------------------------
library(ph2bayes)
PFS9 = subset(interim, obs.time==pt|obs.time>pt)
y = nrow(PFS9)
#number of patient at interim
n = nrow(interim)
#maximum number of patient (interim+future)
nmax= nrow(full_trial)
PFS9_hist = subset(historical, obs.time==pt|obs.time>pt)
#The hyperparameter(shape1) of the beta prior: proportion of responses in historical data
alpha_e = nrow(PFS9_hist)/nrow(historical)
#The hyperparameter(shape2) of the beta prior: proportion of non-responders in historical data
beta_e = 1-alpha_e
#p_s:the efficacy threshold from the standard treatment  
#theta_t= The prespecified target probability

predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.9) 
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.8) 
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.7) 
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.6)
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.5)

#######################################################################

simDat = simTrial(pt=6, n=100, ppm=5, mPFS=12,fol_cut = "Yes", mfol =6, npts=40,
                  drop ="Yes",droprate =0.1, time_drop=12, m=40, interim="Yes",
                  seed1=12348, seed2=145982, seed3 = 2436)

full_trial=simDat$full_trial
cut1 = simDat$cut1
cut2 = simDat$cut2
interim=simDat$interim
nrow(full_trial)
nrow(interim)
table(interim$status)
tail(full_trial)
head(interim)

pt=6
#K-M plot for full data
full_trial.survival= survfit(Surv(full_trial$obs.time,full_trial$status) ~ 1, 
                             error="greenwood", conf.type= "plain")
smr = summary(full_trial.survival, times=pt) 
survminer::surv_median(full_trial.survival)
table(full_trial$status)

p1 = ggsurv(full_trial.survival)
pl = p1+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,30, 2), expand = c(0, 0))+scale_y_continuous(breaks=seq(0,1, 0.1),expand = c(0, 0))+ggtitle("Full data")+theme(plot.margin = unit(c(1.1,1.1,1.1,1.1),"cm"))+
  geom_segment(aes(x=6, y=0, xend=6, yend=0.71), linetype="dashed")+geom_segment(aes(x=0, y=0.71, xend=6, yend=0.71), linetype="dashed")+
  geom_text(aes(x=12.2, y=0.78, label="PFS at 6 months: 0.71 \n 95% CI: (0.62, 0.81)"), size=3.2)

pl
#----------------------------------------------------------
#K-M plot for interim data
interim.survival= survfit(Surv(interim$obs.time,interim$status) ~ 1, 
                          error="greenwood", conf.type= "plain")

smr1 = summary(interim.survival, times=pt) 

survminer::surv_median(interim.survival)

table(interim$status)

p11 = ggsurv(interim.survival)
pl1 = p11+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,40, 2), expand = c(0, 0))+scale_y_continuous(breaks=seq(0,1, 0.1),expand = c(0, 0))+ggtitle("Interim data")+theme(plot.margin = unit(c(1.1,1.1,1.1,1.1),"cm"))+
  geom_segment(aes(x=6, y=0, xend=6, yend=0.72), linetype="dashed")+geom_segment(aes(x=0, y=0.72, xend=6, yend=0.72), linetype="dashed")+
  geom_text(aes(x=11, y=0.8, label="PFS at 6 months: 0.72 \n 95% CI: (0.60, 0.84)"), size=3.2)
pl1

#---------------------------------------------------------------

simDat2 = simTrial(pt=6, n=100, ppm=5, mPFS=10,fol_cut = "Yes", mfol =6, npts=40,
                  drop ="Yes",droprate =0.1, time_drop=12, m=40, interim="Yes",
                  seed1=231, seed2=987, seed3 = 1238)


historical = simDat2$full_trial
cut1_historical = simDat2$cut
table(historical$status)


#K-M plot for historical data
historical.survival= survfit(Surv(historical$obs.time,historical$status) ~ 1, 
                             error="greenwood", conf.type= "plain")
smr2 = summary(historical.survival, times=pt) 

survminer::surv_median(historical.survival)

table(historical$status)

p12 = ggsurv(historical.survival)
pl2 = p12+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,40, 2), expand = c(0, 0))+scale_y_continuous(breaks=seq(0,1, 0.1),expand = c(0, 0))+ggtitle("Historical data")+theme(plot.margin = unit(c(1.1,1.1,1.1,1.1),"cm"))+
  geom_segment(aes(x=6, y=0, xend=6, yend=0.68), linetype="dashed")+geom_segment(aes(x=0, y=0.68, xend=6, yend=0.68), linetype="dashed")+
  geom_text(aes(x=12.5, y=0.68, label="PFS at 6 months: 0.68 \n 95% CI: (0.58, 0.77)"), size=3.2)

pl2

grid.arrange(pl,pl1, pl2, nrow = 1)

########################################################################
#proportion of responses in all simulated data sets
#1) proportion of responses in the historical data
pt=6
PFS9_hist = subset(historical, obs.time==pt|obs.time>pt)
#Tproportion of responses
alpha_e = nrow(PFS9_hist)/nrow(historical)
#2) proportion of responses in the interim data
PFS9_interim = subset(interim, obs.time==pt|obs.time>pt)
#Tproportion of responses
alpha_e = nrow(PFS9_interim)/nrow(interim)
#2) proportion of responses in the full data
PFS9_full = subset(full_trial, obs.time==pt|obs.time>pt)
#Tproportion of responses
alpha_e = nrow(PFS9_full)/nrow(full_trial)

########################################################################
Predictive_probability(data=interim, statistic =predProb_interim, R=1000, pt=6, prior.data=historical,
                       n_max=100, P_T=0.9, theta=0.65, seed=123)

Predictive_probability(data=interim, statistic =predProb_interim, R=1000, pt=6, prior.data=historical,
                       n_max=100, P_T=0.8, theta=0.65, seed=123)


Predictive_probability(data=interim, statistic =predProb_interim, R=1000, pt=6, prior.data=historical,
                       n_max=100, P_T=0.7, theta=0.65, seed=123)

Predictive_probability(data=interim, statistic =predProb_interim, R=1000, pt=6, prior.data=historical,
                       n_max=100, P_T=0.6, theta=0.65, seed=123)

Predictive_probability(data=interim, statistic =predProb_interim, R=1000, pt=6, prior.data=historical,
                       n_max=100, P_T=0.5, theta=0.65, seed=123)

#------------------------------------------------------------------------
library(ph2bayes)
PFS9 = subset(interim, obs.time==pt|obs.time>pt)
y = nrow(PFS9)
#number of patient at interim
n = nrow(interim)
#maximum number of patient (interim+future)
nmax= nrow(full_trial)

PFS9_hist = subset(historical, obs.time==pt|obs.time>pt)
#The hyperparameter(shape1) of the beta prior: proportion of responses in historical data
alpha_e = nrow(PFS9_hist)/nrow(historical)

#The hyperparameter(shape2) of the beta prior: proportion of non-responders in historical data
beta_e = 1-alpha_e

#p_s:the efficacy threshold from the standard treatment  
#theta_t= The prespecified target probability

predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.9) 
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.8) 
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.7) 
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.6)
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.5)

##############################################################################
#with non-informative prior
library(survival)
library(boot)
library(survminer)
library(ph2bayes)
library(ggplot2)
library(gridExtra)
ggsurv = function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                  cens.col = 'purple', lty.est = 1, lty.ci = 3,
                  cens.shape = 3, back.white = F, xlab = 'Time in months',
                  ylab = 'Proportion surviving', main = ""){
  
  library(ggplot2)
  strata <- ifelse(is.null(s$strata) ==T, 1, length(s$strata))
  stopifnot(length(surv.col) == 1 | length(surv.col) == strata)
  stopifnot(length(lty.est) == 1 | length(lty.est) == strata)
  
  ggsurv.s <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'purple', lty.est = 1, lty.ci = 3,
                       cens.shape = 3, back.white = F, xlab = 'Time in months',
                       ylab = 'Proportion surviving', main = ''){
    
    dat <- data.frame(time = c(0, s$time),
                      surv = c(1, s$surv),
                      up = c(1, s$upper),
                      low = c(1, s$lower),
                      cens = c(0, s$n.censor))
    dat.cens <- subset(dat, cens != 0)
    
    col <- ifelse(surv.col == 'gg.def', 'black', surv.col)
    
    pl <- ggplot(dat, aes(x = time, y = surv)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(col = col, lty = lty.est)
    
    pl <- if(CI == T | CI == 'def') {
      pl + geom_step(aes(y = up), color = col, lty = lty.ci) +
        geom_step(aes(y = low), color = col, lty = lty.ci)
    } else (pl)
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  
  ggsurv.m <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 3,
                       cens.shape = 3, back.white = F, xlab = 'Time in months',
                       ylab = 'Proportion surviving', main = '') {
    n <- s$strata
    
    groups <- factor(unlist(strsplit(names
                                     (s$strata), '='))[seq(2, 2*strata, by = 2)])
    gr.name <-  unlist(strsplit(names(s$strata), '='))[1]
    gr.df <- vector('list', strata)
    ind <- vector('list', strata)
    n.ind <- c(0,n); n.ind <- cumsum(n.ind)
    for(i in 1:strata) ind[[i]] <- (n.ind[i]+1):n.ind[i+1]
    
    for(i in 1:strata){
      gr.df[[i]] <- data.frame(
        time = c(0, s$time[ ind[[i]] ]),
        surv = c(1, s$surv[ ind[[i]] ]),
        up = c(1, s$upper[ ind[[i]] ]),
        low = c(1, s$lower[ ind[[i]] ]),
        cens = c(0, s$n.censor[ ind[[i]] ]),
        group = rep(groups[i], n[i] + 1))
    }
    
    dat <- do.call(rbind, gr.df)
    dat.cens <- subset(dat, cens != 0)
    
    pl <- ggplot(dat, aes(x = time, y = surv, group = group)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(aes(col = group, lty = group))
    
    col <- if(length(surv.col == 1)){
      scale_colour_manual(name = gr.name, values = rep(surv.col, strata))
    } else{
      scale_colour_manual(name = gr.name, values = surv.col)
    }
    
    pl <- if(surv.col[1] != 'gg.def'){
      pl + col
    } else {pl + scale_colour_discrete(name = gr.name)}
    
    line <- if(length(lty.est) == 1){
      scale_linetype_manual(name = gr.name, values = rep(lty.est, strata))
    } else {scale_linetype_manual(name = gr.name, values = lty.est)}
    
    pl <- pl + line
    
    pl <- if(CI == T) {
      if(length(surv.col) > 1 && length(lty.est) > 1){
        stop('Either surv.col or lty.est should be of length 1 in order
             to plot 95% CI with multiple strata')
      }else if((length(surv.col) > 1 | surv.col == 'gg.def')[1]){
        pl + geom_step(aes(y = up, color = group), lty = lty.ci) +
          geom_step(aes(y = low, color = group), lty = lty.ci)
      } else{pl +  geom_step(aes(y = up, lty = group), col = surv.col) +
          geom_step(aes(y = low,lty = group), col = surv.col)}
    } else {pl}
    
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  pl <- if(strata == 1) {ggsurv.s(s, CI , plot.cens, surv.col ,
                                  cens.col, lty.est, lty.ci,
                                  cens.shape, back.white, xlab,
                                  ylab, main)
  } else {ggsurv.m(s, CI, plot.cens, surv.col ,
                   cens.col, lty.est, lty.ci,
                   cens.shape, back.white, xlab,
                   ylab, main)}
  pl
}


#function to simulate a clinical trial
simTrial=function(pt=9, n=500, ppm=20, mPFS=6, fol_cut = "Yes", mfol =20, npts=200, 
                  drop ="Yes",droprate =0.1,time_drop=12, m=200, interim ="Yes", 
                  seed1=123, seed2=145, seed3=243){
  #simulating the enrolling time
  set.seed(seed1)
  enrol.time = -log(runif(n))/ppm
  for(i in 2:length(enrol.time))
  {   
    enrol.time[i] = enrol.time[i]+enrol.time[i-1]
  }
  #simulating the event time
  lambda0 = log(2)/mPFS
  set.seed(seed2)
  event.time = -log(runif(n))/lambda0
  #generating censoring status based on drop-out for no reason
  if (drop=="Yes"){
    set.seed(seed3)
    lambda1 = -log(1-droprate)/time_drop
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
    interim_dat=sim_dat[1:m,]
    cut2=interim_dat$enrol.time[m]+mfol
    rest = sim_dat[m+1:n,]
    rest_interim = subset(rest, enrol.time<cut2| enrol.time==cut2)
    interim = rbind(interim_dat, rest_interim)
    interim$status = ifelse(interim$study.time>cut2, 0, interim$status)
    
    return(list(full_trial=sim_dat, interim=interim, cut1=cut, cut2=cut2))
  }else{
    return(list(full_trial=sim_dat, cut=cut)) 
  }
}


simDat = simTrial(pt=6, n=100, ppm=5, mPFS=12,fol_cut = "Yes", mfol =6, npts=40,
                  drop ="Yes",droprate =0.1, time_drop=12, m=40, interim="Yes",
                  seed1=12348, seed2=145982, seed3 = 2436)

full_trial=simDat$full_trial
cut = simDat$cut
interim=simDat$interim
nrow(full_trial)

tail(full_trial)
head(interim)

pt=6

#K-M plot for full data
full_trial.survival= survfit(Surv(full_trial$obs.time,full_trial$status) ~ 1, 
                             error="greenwood", conf.type= "plain")

smr = summary(full_trial.survival, times=pt) 

survminer::surv_median(full_trial.survival)

table(full_trial$status)

#----------------------------------------------------------
#K-M plot for interim data
interim.survival= survfit(Surv(interim$obs.time,interim$status) ~ 1, 
                          error="greenwood", conf.type= "plain")

grid.arrange(pl,pl1,nrow = 1)

smr1 = summary(interim.survival, times=pt) 

survminer::surv_median(interim.survival)

table(interim$status)

########################################################################
#calculating PPOS using the proposed method (without bootstraping)
predProb_interim_no_hist = function(data, pt, n_max, P_T, theta, index){
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
  predProb= vector(length=4)
  predProb[1]=predPost.mean
  predProb[2]=predPost.var
  predProb[3]=delta
  predProb[4]=PredProb
  
  return(predProb) 
}

#final function to calculate PPOS with bootstraping
Predictive_probability_no_hist = function(data, statistic, R, pt, n_max, P_T, theta, seed){
  set.seed(seed)
  results = boot(data =data, statistic=statistic, R=R, pt=pt, 
                 n_max=n_max, P_T=P_T, theta=theta)
  bs_stat = apply(results$t,2,mean, na.rm=T)
  bootCI_percent = matrix(ncol=2, nrow=4)
  for(i in 1:4){
    bootCI_percent[i,1] =boot.ci(results, index=i, type=c("perc"))$percent[4] 
    bootCI_percent[i,2] =boot.ci(results, index=i, type=c("perc"))$percent[5] 
  }
  str=paste0("    Statistics                    95% Percentile CI\nPosterior mean         =",
             round(bs_stat[1],3),"     (",
             round(bootCI_percent[1,1],3),", ",round(bootCI_percent[1,2],3),")\nPosterior variance     =",
             round(bs_stat[2],3), "     (", round(bootCI_percent[2,1],3), ","
             ,round(bootCI_percent[2,2],3),")\nSuccess Threshold      =", round(bs_stat[3],3),
             "     (" ,round(bootCI_percent[3,1],3), ","  ,round(bootCI_percent[3,2],3),
             ")\nPredictive Probability =", round(bs_stat[4],3), "     (", 
             round(bootCI_percent[4,1],3),"," , round(bootCI_percent[4,2],3),")")
  
  cat(str)
  
}

#######################################################
Predictive_probability_no_hist(data=interim, statistic =predProb_interim_no_hist, R=1000, pt=6,
                               n_max=100, P_T=0.9, theta=0.65, seed=254)

Predictive_probability_no_hist(data=interim, statistic =predProb_interim_no_hist, R=1000, pt=6, 
                               n_max=100, P_T=0.8, theta=0.65, seed=123)


Predictive_probability_no_hist(data=interim, statistic =predProb_interim_no_hist, R=1000, pt=6, 
                               n_max=100, P_T=0.7, theta=0.65, seed=123)

Predictive_probability_no_hist(data=interim, statistic =predProb_interim_no_hist, R=1000, pt=6,
                               n_max=100, P_T=0.6, theta=0.65, seed=123)

Predictive_probability_no_hist(data=interim, statistic =predProb_interim_no_hist, R=1000, pt=6,
                               n_max=100, P_T=0.5, theta=0.65, seed=123)

#------------------------------------------------------------------------
library(ph2bayes)
PFS9 = subset(interim, obs.time==pt|obs.time>pt)
y = nrow(PFS9)
#number of patient at interim
n = nrow(interim)

#maximum number of patient (interim+future)
nmax= nrow(full_trial)
alpha_e = 0.5

#The hyperparameter(shape2) of the beta prior: proportion of non-responders in historical data
beta_e = 0.5

#p_s:the efficacy threshold from the standard treatment  
#theta_t= The prespecified target probability

predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.9) 
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.8) 
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.7) 
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.6)
predprob(y, n, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.5)







