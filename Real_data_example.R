library(survival)
library(boot)
library(survminer)
library(ph2bayes)
library(ggplot2)

ggsurv = function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                  cens.col = 'purple', lty.est = 1, lty.ci = 2,
                  cens.shape = 3, back.white = F, xlab = 'Time in months',
                  ylab = 'proportion surviving', main = ""){
  
  library(ggplot2)
  strata <- ifelse(is.null(s$strata) ==T, 1, length(s$strata))
  stopifnot(length(surv.col) == 1 | length(surv.col) == strata)
  stopifnot(length(lty.est) == 1 | length(lty.est) == strata)
  
  ggsurv.s <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'purple', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time in months',
                       ylab = 'proportion surviving', main = ''){
    
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
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time in months',
                       ylab = 'proportion surviving', main = '') {
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
pt=9
dat0= read.csv("~/Desktop/dat2.csv")
colnames(dat0)[1]="USUBJID"
dat1= read.csv("~/Desktop/adsl.csv")
names(dat0)
names(dat1)
dat2 = subset(dat1, ARMN==1)
dat3 = dat2[, c(1,2,32,33,35,36,38,39,40,87,138,139,152,153,155,163,165,166,168,169,171,206, 207,208,209)]
dat = merge(dat0, dat3, by="USUBJID")
dat = dat[with(dat, order(TRTSDTC, TRTSTM)),]
names(dat)[28] = "obs.time"
names(dat)[29] = "Censor"
dat$status = ifelse(dat$Censor==1,0, 1)
table(dat$status)
table(dat$status)/sum(table(dat$status))
#table(dat$DSFUREAS)
censored_dat = subset(dat, status==0)
head(censored_dat)
#table(censored_dat$DSTREAS)
table(dat$DSFUREAS)
(table(dat$DSTREAS))*100/nrow(dat)
names(dat)
range(dat$obs.time)
surv_dat= survfit(Surv(dat$obs.time,dat$status) ~ 1,
                  error="greenwood", conf.type= "plain")
pt=9
survminer::surv_median(surv_dat)
smr =summary(surv_dat, times=pt)

p0 = ggsurv(surv_dat)
p0 = p0+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,40, 2))+theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=11))+
  scale_y_continuous(name ="Proportion surviving")+
  ggtitle("")+geom_text(x = 28, y=0.95, label="PFS at 9 months: 0.76  \n 95% CI: (0.73, 0.79)", size=3.2)

# + geom_text(aes(x=0.61, label="Max sample size = 871", y=0.98),
#                                                   color ="black",angle=0, size=3.2)+
#   geom_text(aes(x=0.61, label="Interim sample size = 150", y=0.93),
#             color ="black",angle=0, size=3.2)
p0
#------------------------------------------------------
library(lubridate)
# 4 months after 100th patient enrollment
dat$study.time = as.Date(dat$TRTSDTC)+ dat$obs.time
enrol_100 = dat[100,"TRTSDTC"]
cut_off1 = as.Date(enrol_100) %m+% months(4)
interim1 = subset(dat,as.Date(TRTSDTC) <= cut_off1)
dim(interim1)
table(interim1$status)
#interim1$EOSDTC= as.Date(interim1$EOSDTC, format = "%Y-%m-%d")
interim1$status = ifelse(as.Date((interim1$study.time))>as.Date(cut_off1), 0, interim1$status)
table(interim1$status)
interim.surv1= survfit(Surv(interim1$obs.time,interim1$status) ~ 1,
                       error="greenwood", conf.type= "plain")
smr1 = summary(interim.surv1, times=pt)

survminer::surv_median(interim.surv1)

p1 = ggsurv(interim.surv1)
p1 = p1+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,40, 2))+theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=11))+
  scale_y_continuous(name ="Proportion surviving")+
  ggtitle("Interim analysis 4 months after the 100th patient's enrollment")+geom_text(aes(x=26, y=0.98, label="PFS at 9 months: 0.77 \n 95% CI: (0.72, 0.82)"), size=3.2)
p1

######################################################################
# 4 months after 200th patient enrollment
enrol_200 = dat[200,"TRTSDTC"]
cut_off2 = as.Date(enrol_200) %m+% months(4)
interim2 = subset(dat,as.Date(TRTSDTC) <= cut_off2)
dim(interim2)
table(interim2$status)
interim2$status = ifelse(as.Date((interim2$study.time))>as.Date(cut_off2), 0, interim2$status)
table(interim2$status)
interim.surv2= survfit(Surv(interim2$obs.time,interim2$status) ~ 1,
                       error="greenwood", conf.type= "plain")
smr2 = summary(interim.surv2, times=pt)

survminer::surv_median(interim.surv2)

p2 = ggsurv(interim.surv2)
p2 = p2+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,40, 2))+theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=11))+
  scale_y_continuous(name ="Proportion surviving")+
  ggtitle("Interim analysis 4 months after the 200th patient's enrollment")+geom_text(aes(x=26, y=0.98, label="PFS at 9 months: 0.77 \n 95% CI: (0.73, 0.82)"), size=3.2)
p2

###################################################################
# 4 months after 300th patient enrollment
enrol_300 = dat[300,"TRTSDTC"]
cut_off3 = as.Date(enrol_300) %m+% months(6)
interim3 = subset(dat,as.Date(TRTSDTC) <= cut_off3)
dim(interim3)
table(interim3$status)
interim3$status = ifelse(as.Date((interim3$study.time))>as.Date(cut_off3), 0, interim3$status)
table(interim3$status)

interim.surv3= survfit(Surv(interim3$obs.time,interim3$status) ~ 1,
                       error="greenwood", conf.type= "plain")
smr3 = summary(interim.surv3, times=pt)

survminer::surv_median(interim.surv3)

p3 = ggsurv(interim.surv3)
p3 = p3+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,40, 2))+theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=11))+
  scale_y_continuous(name ="Proportion surviving")+
  ggtitle("Interim analysis 4 months after the 300th patient's enrollment")+geom_text(aes(x=26, y=0.98, label="PFS at 9 months: 0.77 \n 95% CI: (0.74, 0.81)"), size=3.2)
p3

#####################################################################
#proportion of responses in the interim data
PFS9_interim1 = subset(interim1, obs.time==pt|obs.time>pt)
#Tproportion of responses
alpha_e1 = nrow(PFS9_interim1)/nrow(interim1)
#------------------------------------------------------
#proportion of responses in the interim data
PFS9_interim2 = subset(interim2, obs.time==pt|obs.time>pt)
#Tproportion of responses
alpha_e2 = nrow(PFS9_interim2)/nrow(interim2)
#------------------------------------------------------
#proportion of responses in the interim data
PFS9_interim3 = subset(interim3, obs.time==pt|obs.time>pt)
#Tproportion of responses
alpha_e3 = nrow(PFS9_interim3)/nrow(interim3)
#------------------------------------------------------
#proportion of responses in the full data
PFS9_full = subset(dat, obs.time==pt|obs.time>pt)
#proportion of responses
alpha_e = nrow(PFS9_full)/nrow(dat)

##########################################################################
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

Predictive_probability_no_hist(data=interim1, statistic =predProb_interim_no_hist, R=1000, pt=9,
                               n_max=871, P_T=0.8, theta=0.80, seed=12376)



Predictive_probability_no_hist(data=interim1, statistic =predProb_interim_no_hist, R=1000, pt=9,
                               n_max=871, P_T=0.8, theta=0.75, seed=123)





Predictive_probability_no_hist(data=interim1, statistic =predProb_interim_no_hist, R=1000, pt=9,
                               n_max=871, P_T=0.8, theta=0.70, seed=123)





Predictive_probability_no_hist(data=interim1, statistic =predProb_interim_no_hist, R=1000, pt=9,
                               n_max=871, P_T=0.8, theta=0.65, seed=123)



#------------------------------------------------------------------------
library(ph2bayes)
pt=9
#responders are those whose their PFS is greater than 9 month
PFS91 = subset(interim1, obs.time==pt|obs.time>pt)
y1 = nrow(PFS91)

#number of patient at interim
n1 = nrow(interim1)

#maximum number of patient (interim+future)
nmax= nrow(dat)

alpha_e = 0.5

#The hyperparameter(shape2) of the beta prior: proportion of non-responders in historical data
beta_e = 0.5

#p_s:the efficacy threshold from the standard treatment or target threshold
#theta_t= The prespecified success threshold

predprob(y1, n1, nmax, alpha_e, beta_e, p_s=0.80,theta_t=0.8)


predprob(y1, n1, nmax, alpha_e, beta_e, p_s=0.75,theta_t=0.8)


predprob(y1, n1, nmax, alpha_e, beta_e, p_s=0.70,theta_t=0.8)


predprob(y1, n1, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.8)

#---------------------------------------------------------------------

Predictive_probability_no_hist(data=interim2, statistic =predProb_interim_no_hist, R=1000, pt=9,
                               n_max=871, P_T=0.8, theta=0.80, seed=123)



Predictive_probability_no_hist(data=interim2, statistic =predProb_interim_no_hist, R=1000, pt=9,
                               n_max=871, P_T=0.8, theta=0.75, seed=123)



Predictive_probability_no_hist(data=interim2, statistic =predProb_interim_no_hist, R=1000, pt=9,
                               n_max=871, P_T=0.8, theta=0.70, seed=123)



Predictive_probability_no_hist(data=interim2, statistic =predProb_interim_no_hist, R=1000, pt=9,
                               n_max=871, P_T=0.8, theta=0.65, seed=123)

#------------------------------------------------------------------------
library(ph2bayes)
PFS92 = subset(interim2, obs.time==pt|obs.time>pt)
y2 = nrow(PFS92)

#number of patient at interim
n2 = nrow(interim2)

#maximum number of patient (interim+future)
nmax= nrow(dat)

alpha_e = 0.5

#The hyperparameter(shape2) of the beta prior: proportion of non-responders in historical data
beta_e = 0.5

#p_s:the efficacy threshold from the standard treatment 
#theta_t= The prespecified target probability

predprob(y2, n2, nmax, alpha_e, beta_e, p_s=0.80,theta_t=0.8)

predprob(y2, n2, nmax, alpha_e, beta_e, p_s=0.75,theta_t=0.8)

predprob(y2, n2, nmax, alpha_e, beta_e, p_s=0.70,theta_t=0.8)

predprob(y2, n2, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.8)

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

Predictive_probability_no_hist(data=interim3, statistic =predProb_interim_no_hist, R=1000, pt=9,
                               n_max=871, P_T=0.8, theta=0.80, seed=123)


Predictive_probability_no_hist(data=interim3, statistic =predProb_interim_no_hist, R=1000, pt=9,
                               n_max=871, P_T=0.8, theta=0.75, seed=123)



Predictive_probability_no_hist(data=interim3, statistic =predProb_interim_no_hist, R=1000, pt=9,
                               n_max=871, P_T=0.8, theta=0.70, seed=123)



Predictive_probability_no_hist(data=interim3, statistic =predProb_interim_no_hist, R=1000, pt=9,
                               n_max=871, P_T=0.8, theta=0.65, seed=123)


#------------------------------------------------------------------------
library(ph2bayes)
PFS93 = subset(interim3, obs.time==pt|obs.time>pt)
y3 = nrow(PFS93)

#number of patient at interim
n3 = nrow(interim3)

#maximum number of patient (interim+future)
nmax= nrow(dat)

alpha_e = 0.5

#The hyperparameter(shape2) of the beta prior: proportion of non-responders in historical data
beta_e = 0.5

#p_s:the efficacy threshold from the standard treatment 
#theta_t= The prespecified target probability

predprob(y3, n3, nmax, alpha_e, beta_e, p_s=0.80,theta_t=0.8)

predprob(y3, n3, nmax, alpha_e, beta_e, p_s=0.75,theta_t=0.8)

predprob(y3, n3, nmax, alpha_e, beta_e, p_s=0.70,theta_t=0.8)

predprob(y3, n3, nmax, alpha_e, beta_e, p_s=0.65,theta_t=0.8)

#########################################
#Now Overall survival
dat0= read.csv("~/Desktop/dat2.csv")
colnames(dat0)[1]="USUBJID"
dat1= read.csv("~/Desktop/adsl.csv")
names(dat0)
colnames(dat1)
dat2 = subset(dat1, ARMN==1)
dat3 = dat2[, c(1,2,32,33,35,36,38,39,40,87,138,139,152,153,155,163,165,166,168,169,171,206, 207,208,209)]
dat = merge(dat0, dat3, by="USUBJID")
dat = dat[with(dat, order(TRTSDTC, TRTSTM)),]
names(dat)[32] = "obs.time"
names(dat)[33] = "Censor"
dat$status = ifelse(dat$Censor==1,0, 1)
table(dat$status)
table(dat$status)/sum(table(dat$status))
#table(dat$DSFUREAS)
#censored_dat = subset(dat, status==0)
#table(censored_dat$DSTREAS)
#table(censored_dat$DSTTERM)
#range(dat$obs.time)
surv_dat= survfit(Surv(dat$obs.time,dat$status) ~ 1,
                  error="greenwood", conf.type= "plain")
pt=18
survminer::surv_median(surv_dat)
smr =summary(surv_dat, times=pt)

p0 = ggsurv(surv_dat)
p0 = p0+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,40, 2))+theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=11))+
  scale_y_continuous(name ="Proportion surviving")+
  ggtitle("K-M curve of the treatment arm")+geom_text(x = 28, y=0.95, label="OS at 18 months: 0.82 \n 95% CI: (0.79, 0.84) ", size=3.2)

p0
#------------------------------------------------------
library(lubridate)
dat$study.time = as.Date(dat$TRTSDTC)+ dat$obs.time
# 4 months after 100th patient enrollment
enrol_100 = dat[100,"TRTSDTC"]
cut_off1 = as.Date(enrol_100) %m+% months(4)
interim1 = subset(dat,as.Date(TRTSDTC) <= cut_off1)
dim(interim1)
table(interim1$status)
interim1$status = ifelse(as.Date((interim1$study.time))>as.Date(cut_off1), 0, interim1$status)
table(interim1$status)

interim.surv1= survfit(Surv(interim1$obs.time,interim1$status) ~ 1,
                       error="greenwood", conf.type= "plain")
smr1 = summary(interim.surv1, times=pt)

survminer::surv_median(interim.surv1)

p1 = ggsurv(interim.surv1)
p1 = p1+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,40, 2))+theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=11))+
  scale_y_continuous(name ="Proportion surviving")+
  ggtitle("Interim analysis 4 months after the 100th patient's enrollment")+geom_text(aes(x=26, y=0.98, label="OS at 18 months: 0.82 \n 95% CI: (0.78, 0.86)"), size=3.2)
p1
######################################################################
# 4 months after 200th patient enrollment
enrol_200 = dat[200,"TRTSDTC"]
cut_off2 = as.Date(enrol_200) %m+% months(4)
interim2 = subset(dat,as.Date(TRTSDTC) <= cut_off2)
dim(interim2)
table(interim2$status)
interim2$status = ifelse(as.Date((interim2$study.time))>as.Date(cut_off2), 0, interim2$status)
table(interim2$status)

interim.surv2= survfit(Surv(interim2$obs.time,interim2$status) ~ 1,
                       error="greenwood", conf.type= "plain")
smr2 = summary(interim.surv2, times=pt)

survminer::surv_median(interim.surv2)

p2 = ggsurv(interim.surv2)
p2 = p2+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,40, 2))+theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=11))+
  scale_y_continuous(name ="Proportion surviving")+
  ggtitle("Interim analysis 4 months after the 200th patient's enrollment")+geom_text(aes(x=26, y=0.98, label="OS at 18 months: 0.83 \n 95% CI: (0.79, 0.86)"), size=3.2)
p2

###################################################################
# 4 months after 300th patient enrollment
enrol_300 = dat[300,"TRTSDTC"]
cut_off3 = as.Date(enrol_300) %m+% months(4)
interim3 = subset(dat,as.Date(TRTSDTC) <= cut_off3)
dim(interim3)
table(interim3$status)
interim3$status = ifelse(as.Date((interim3$study.time))>as.Date(cut_off3), 0, interim3$status)
table(interim3$status)

interim.surv3= survfit(Surv(interim3$obs.time,interim3$status) ~ 1,
                       error="greenwood", conf.type= "plain")
smr3 = summary(interim.surv3, times=pt)

survminer::surv_median(interim.surv3)

p3 = ggsurv(interim.surv3)
p3 = p3+theme(text=element_text(size=11)) +theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(0,40, 2))+theme(plot.margin= unit(c(0.7, 0.7, 0.7,0.7), "cm"))+
  theme(plot.title = element_text(hjust = 0.5, size=11))+
  scale_y_continuous(name ="Proportion surviving")+
  ggtitle("Interim analysis 4 months after the 300th patient's enrollment")+geom_text(aes(x=26, y=0.98, label="OS at 18 months: 0.83 \n 95% CI: (0.80, 0.87)"), size=3.2)
p3

#####################################################################
#proportion of responses in the interim data
PFS9_interim1 = subset(interim1, obs.time==pt|obs.time>pt)
#Tproportion of responses
alpha_e1 = nrow(PFS9_interim1)/nrow(interim1)
#------------------------------------------------------
#proportion of responses in the interim data
PFS9_interim2 = subset(interim2, obs.time==pt|obs.time>pt)
#Tproportion of responses
alpha_e2 = nrow(PFS9_interim2)/nrow(interim2)
#------------------------------------------------------
#proportion of responses in the interim data
PFS9_interim3 = subset(interim3, obs.time==pt|obs.time>pt)
#Tproportion of responses
alpha_e3 = nrow(PFS9_interim3)/nrow(interim3)
#------------------------------------------------------
#proportion of responses in the full data
PFS9_full = subset(dat, obs.time==pt|obs.time>pt)
#proportion of responses
alpha_e = nrow(PFS9_full)/nrow(dat)

##########################################################################
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
#P_T=0.8, theta=0.85
Predictive_probability_no_hist(data=interim1, statistic =predProb_interim_no_hist, R=1000, pt=18,
                               n_max=871, P_T=0.8, theta=0.85, seed=123)


#P_T=0.8, theta=0.80
Predictive_probability_no_hist(data=interim1, statistic =predProb_interim_no_hist, R=1000, pt=18,
                               n_max=871, P_T=0.8, theta=0.80, seed=123)




#P_T=0.8 , theta =0.75
Predictive_probability_no_hist(data=interim1, statistic =predProb_interim_no_hist, R=1000, pt=18,
                               n_max=871, P_T=0.8, theta=0.75, seed=12387)




#P_T=0.8, theta =0.70
Predictive_probability_no_hist(data=interim1, statistic =predProb_interim_no_hist, R=1000, pt=18,
                               n_max=871, P_T=0.8, theta=0.70, seed=123)



#------------------------------------------------------------------------
library(ph2bayes)
pt=18
#responders are those whose their PFS is greater than 9 month
PFS91 = subset(interim1, obs.time==pt|obs.time>pt)
y1 = nrow(PFS91)

#number of patient at interim
n1 = nrow(interim1)

#maximum number of patient (interim+future)
nmax= nrow(dat)

alpha_e = 0.5

#The hyperparameter(shape2) of the beta prior: proportion of non-responders in historical data
beta_e = 0.5

#p_s:the efficacy threshold from the standard treatment or target threshold
#theta_t= The prespecified success threshold

predprob(y1, n1, nmax, alpha_e, beta_e, p_s=0.85,theta_t=0.8)
#

predprob(y1, n1, nmax, alpha_e, beta_e, p_s=0.80,theta_t=0.8)
#

predprob(y1, n1, nmax, alpha_e, beta_e, p_s=0.75,theta_t=0.8)
#

predprob(y1, n1, nmax, alpha_e, beta_e, p_s=0.70,theta_t=0.8)
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#P_T=0.8, theta=0.85
Predictive_probability_no_hist(data=interim2, statistic =predProb_interim_no_hist, R=1000, pt=18,
                               n_max=871, P_T=0.8, theta=0.85, seed=123)


Predictive_probability_no_hist(data=interim2, statistic =predProb_interim_no_hist, R=1000, pt=18,
                               n_max=871, P_T=0.8, theta=0.80, seed=123)


Predictive_probability_no_hist(data=interim2, statistic =predProb_interim_no_hist, R=1000, pt=18,
                               n_max=871, P_T=0.8, theta=0.75, seed=123)


Predictive_probability_no_hist(data=interim2, statistic =predProb_interim_no_hist, R=1000, pt=18,
                               n_max=871, P_T=0.8, theta=0.70, seed=123)


#------------------------------------------------------------------------
library(ph2bayes)
PFS92 = subset(interim2, obs.time==pt|obs.time>pt)
y2 = nrow(PFS92)


#number of patient at interim
n2 = nrow(interim2)

#maximum number of patient (interim+future)
nmax= nrow(dat)

alpha_e = 0.5

#The hyperparameter(shape2) of the beta prior: proportion of non-responders in historical data
beta_e = 0.5

#p_s:the efficacy threshold from the standard treatment 
#theta_t= The prespecified target probability

predprob(y2, n2, nmax, alpha_e, beta_e, p_s=0.85,theta_t=0.8)
#
predprob(y2, n2, nmax, alpha_e, beta_e, p_s=0.80,theta_t=0.8)
#
predprob(y2, n2, nmax, alpha_e, beta_e, p_s=0.75,theta_t=0.8)
#
predprob(y2, n2, nmax, alpha_e, beta_e, p_s=0.70,theta_t=0.8)
#
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

Predictive_probability_no_hist(data=interim3, statistic =predProb_interim_no_hist, R=1000, pt=18,
                               n_max=871, P_T=0.8, theta=0.85, seed=123)




Predictive_probability_no_hist(data=interim3, statistic =predProb_interim_no_hist, R=1000, pt=18,
                               n_max=871, P_T=0.8, theta=0.80, seed=123)



Predictive_probability_no_hist(data=interim3, statistic =predProb_interim_no_hist, R=1000, pt=18,
                               n_max=871, P_T=0.8, theta=0.75, seed=12398)


Predictive_probability_no_hist(data=interim3, statistic =predProb_interim_no_hist, R=1000, pt=18,
                               n_max=871, P_T=0.8, theta=0.70, seed=123)




#------------------------------------------------------------------------
library(ph2bayes)
PFS93 = subset(interim3, obs.time==pt|obs.time>pt)
y3 = nrow(PFS93)

#number of patient at interim
n3 = nrow(interim3)

#maximum number of patient (interim+future)
nmax= nrow(dat)

alpha_e = 0.5

#The hyperparameter(shape2) of the beta prior: proportion of non-responders in historical data
beta_e = 0.5

#p_s:the efficacy threshold from the standard treatment 
#theta_t= The prespecified target probability

predprob(y3, n3, nmax, alpha_e, beta_e, p_s=0.85,theta_t=0.8)
#
predprob(y3, n3, nmax, alpha_e, beta_e, p_s=0.80,theta_t=0.8)
#
predprob(y3, n3, nmax, alpha_e, beta_e, p_s=0.75,theta_t=0.8)
#
predprob(y3, n3, nmax, alpha_e, beta_e, p_s=0.70,theta_t=0.8)
#
###################################################################

