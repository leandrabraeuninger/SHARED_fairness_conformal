rm(list=ls())
setwd("~/Turing Projects/Conformal_Prediction/Recidivism/")
read.csv("compas_rates.csv")->data
read.csv("two_yr_recid.csv")->outcome
head(data)
head(outcome)
######
match(paste(tolower(outcome$last),tolower(outcome$first)),paste(tolower(data$LastName),tolower(data$FirstName)))->m
outcome[-which(is.na(m)),]->out_r
data[m[-which(is.na(m))],]->data_r
head(out_r)
head(data_r)
plot(data_r$DecileScore,out_r$v_decile_score)
cor.test(data_r$RawScore,out_r$decile_score)###ok so this is highly correlated
hist(data_r$RawScore[which(data_r$Ethnic_Code_Text=="African-American")])
hist(data_r$RawScore[which(data_r$Ethnic_Code_Text=="Caucasian")])
boxplot(data_r$RawScore[which(data_r$Ethnic_Code_Text=="African-American")],data_r$RawScore[which(data_r$Ethnic_Code_Text=="Caucasian")])
wilcox.test(data_r$RawScore[which(data_r$Ethnic_Code_Text=="African-American")],data_r$RawScore[which(data_r$Ethnic_Code_Text=="Caucasian")])
range(data_r$RawScore)
###############################
colnames(data_r)
colnames(out_r)
head(out_r)
head(data_r)
times<-rep(2,nrow(out_r))
as.numeric(as.Date(out_r$r_offense_date)-as.Date(out_r$screening_date))/356.25->re_of_times
times[which(out_r$two_year_recid==1)]<-re_of_times[which(out_r$two_year_recid==1)]
events<-out_r$two_year_recid
cbind(times,events)
library('survival')
summary(coxph(Surv(times,events)~data_r$Ethnic_Code_Text))###race is not a predictor of re-offence
summary(coxph(Surv(times,events)~data_r$RawScore))###recidivism score is
#############################ok so we need to change the outcome to predicting how soon they will re-offend
###################
##so using this i develop 3 bins
##low risk high risk and medium risk
###################
cbind(times,events,data_r$RawScore)->surv_events
head(surv_events)
hist(surv_events[which(surv_events[,2]==1),1])###ok cool so there is a dist
max(surv_events[,1])##some re-offend in like 3 years
hist(surv_events[,3])
plot(surv_events[which(surv_events[,2]==1),3],surv_events[which(surv_events[,2]==1),1])
####ok so lets just bin 0-6 months, 6 months-1yr, 1yr1.5yr, 1.5yr-2yr, more than 2 yrs
###so this is 5 bins
data.frame(surv_events)->surv_events
####
bins<-rep(0,nrow(surv_events))
table(surv_events$times[which(surv_events$events==0)])###cool
####
bins[which(surv_events$times < 2)]<-1
#bins[which(surv_events$times <= 1.5)]<-2
bins[which(surv_events$times <= 1)]<-2
#bins[which(surv_events$times <= 0.5)]<-4
plot(bins,surv_events$times)###cool
# boxplot(surv_events$V3[which(bins==1)],surv_events$V3[which(bins==2)],surv_events$V3[which(bins==3)],
#         surv_events$V3[which(bins==4)])
boxplot(surv_events$V3[which(bins==0)],surv_events$V3[which(bins==1)],surv_events$V3[which(bins==2)])

hist(surv_events$V3)

summary(lm(bins~surv_events$V3))###ok so it is sig
length(which(surv_events[,2]==1))###so we have 3000+ re-offenders
surv_events$V3+abs(min(surv_events$V3))->ss
hist(ss/max(ss))
ss/max(ss)->sig_score
# boxplot(sig_score[which(bins==1)],sig_score[which(bins==2)],sig_score[which(bins==3)],
#         sig_score[which(bins==4)])
boxplot(sig_score[which(bins==0)],sig_score[which(bins==1)],sig_score[which(bins==2)])
####its just a linear transformation
wilcox.test(sig_score[which(bins==4)],sig_score[which(bins==2)])##ok
wilcox.test(sig_score[which(bins==3)],sig_score[which(bins==2)])##ok
wilcox.test(sig_score[which(bins==4)],sig_score[which(bins==3)])##ok
wilcox.test(sig_score[which(bins==1)],sig_score[which(bins==2)])##ok
wilcox.test(sig_score[which(bins==0)],sig_score[which(bins==1)])##ok
###########
wilcox.test(srat_score[which(bins==0)],srat_score[which(bins==1)])##ok
boxplot(srat_score[which(bins==0)],srat_score[which(bins==1)])##ok
#####
wilcox.test(srat_score[which(bins==2)],srat_score[which(bins==1)])##ok
table(srat_score)
#########################################
####right so the score is weaker for the early vs late offenders
####do this in more detail
# hist(1/(1+exp(-surv_events$V3)))
# ####lets see which of these we can
# #####but they are associated!
# ##################make a subset aa and cauc pull out a set and do group based conformal on the raw score as a prob
# ##########use sigmoid function 1/(1+e^-x)
# 1/(1+exp(-data_r$RawScore))->sig_score
hist(sig_score)
summary(coxph(Surv(times,events)~data_r$Ethnic_Code_Text))###race is not a predictor of re-offence
summary(coxph(Surv(times,events)~sig_score))###recidivism score is
boxplot(sig_score[which(data_r$Ethnic_Code_Text=="African-American")],sig_score[which(data_r$Ethnic_Code_Text=="Caucasian")])
wilcox.test(sig_score[which(data_r$Ethnic_Code_Text=="African-American")],sig_score[which(data_r$Ethnic_Code_Text=="Caucasian")])
#######
library('survminer')
data.frame(times,events,sig_score,data_r$Ethnic_Code_Text)->all
colnames(all)<-c("time","event","score","race")
surv_cutpoint(all,time="time",event = "event",variables ="score")->a1
###################
#make it more refined
all[which(events==1),]->all2
all2[events[which(all2$time>1)],2]<-0
surv_cutpoint(all2,time="time",event = "event",variables ="score")->a2
###################so looks like 3 catagories
unlist(c(a1$cutpoint[1],a2$cutpoint[1]))->cps####this is the prob of non re-offence
cps
srat_score<-rep(2,length(sig_score))
srat_score[which(all$score<cps[1])]<-1
srat_score[which(all$score<cps[2])]<-0
table(srat_score)
plot(srat_score,all$time)
summary(coxph(Surv(times,events)~srat_score))###recidivism score is
cbind(all,srat_score)->all3
data.frame(all3)->all3
colnames(all3)
plot(all3$score,all3$srat_score)

fit <- survfit(Surv(time,event)~srat_score,data=all3)
###################
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  #conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in years",   # customize X axis label
  xlim=c(0,2),
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("low","med","high"),    # change legend labels.
  palette = 
    c("pink","#E7B800", "#2E9FDF") # custom color palettes.
)

#########
table(all3$race)
which(all3$race=="Caucasian")->caus
which(all3$race=="African-American")->aa
all3[c(caus,aa),]->all3
######################
####you have 3 cutoffs but they are not always defined by the training data - so lets think about it and just try with 2
b<-30
train_calib_test<-mat.or.vec(b,14)
for(i in 1:b){
  print(i)
  length(aa)+length(caus)###6005 samples
  #########so now we need three groups - 2000 training points, 2000 calibration test, remaining 2k test
  samp<-sample(nrow(all3),4000,replace = F)
  all3[samp[1:2000],]->train
  all3[samp[2001:4000],]->calib
  all3[-samp,]->test
  #
  # table(train$race)
  # table(test$race)
  # table(calib$race)
  #################
  ###first use training data to determine cut offs
  surv_cutpoint(train,time="time",event = "event",variables ="score")->a1
  # train[which(train$event==1),]->train2
  # #train->train2
  # train2$event[which(train2$time>1)]<-0
  # train2$time[which(train2$time>1)]<-1
  # #cbind(train2$event,train$event)
  # surv_cutpoint(train2,time="time",event = "event",variables ="score")->a2
  # unlist(c(a1$cutpoint[1],a2$cutpoint[1]))->cps
  cps<-as.numeric(a1$cutpoint[1])
  pred_train<-rep(2,nrow(train))
  pred_train[which(as.numeric(train$score)<cps[1])]<-1
  #pred_train[which(train$score<cps[1])]<-0
  #table(pred_train)
  fit <- survfit(Surv(time,event)~pred_train,data=train)
  print(surv_pvalue(fit))
  ###################
  ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = TRUE,         # show confidence intervals for
    # point estimaes of survival curves.
    #conf.int.style = "step",  # customize style of confidence intervals
    xlab = "Time in years",   # customize X axis label
    xlim=c(0,2),
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    risk.table = "abs_pct",  # absolute number and percentage at risk.
    risk.table.y.text.col = T,# colour risk table text annotations.
    risk.table.y.text = FALSE,# show bars instead of names in text annotations
    # in legend of risk table.
    ncensor.plot = TRUE,      # plot the number of censored subjects at time t
    surv.median.line = "hv",  # add the median survival pointer.
    legend.labs =
      c("low","high"),    # change legend labels.
    palette =
      c("#E7B800", "#2E9FDF") # custom color palettes.
  )
  #dev.off()
  #################so we have our 3 groups selected via our cutoff and our 3 outcomes
  #################look at calib data asign outcomes
  calib_outcome<-rep(2,nrow(calib))
  # calib_outcome[which(calib$event==0)]<-0
  # calib_outcome[which(calib$event==1&calib$time<2)]<-1
  # calib_outcome[which(calib$event==1&calib$time<=1)]<-2
  calib_outcome[which(calib$event==1)]<-1
  #table(calib_outcome)
  #################
  #################now define the conformal score
  #cps
  #calib$score
  #hist(abs(cps[2]-calib$score)[which(calib_outcome==0)])
  conformal_score<-rep(NA,nrow(calib))
  ##
  conformal_score[which(calib_outcome==2)]<-abs(cps[1]-calib$score)[which(calib_outcome==2)]
  conformal_score[which(calib_outcome==2 & calib$score>cps[1])]<-0
  ##
  conformal_score[which(calib_outcome==1)]<-abs(cps[1]-calib$score)[which(calib_outcome==1)]
  conformal_score[which(calib_outcome==1 & calib$score<cps[1])]<-0
  # conformal_score[which(calib_outcome==2)]<-abs(cps[2]-calib$score)[which(calib_outcome==2)]
  # conformal_score[which(calib_outcome==2 & calib$score>cps[2])]<-0
  # ####
  # conformal_score[which(calib_outcome==1)]<-apply(cbind(abs(cps[2]-calib$score),abs(cps[1]-calib$score)),1,min)[which(calib_outcome==1)]
  # conformal_score[which(calib_outcome==1 & calib$score<cps[2] & calib$score>cps[1])]<-0
  # ####
  # conformal_score[which(calib_outcome==0)]<-abs(cps[1]-calib$score)[which(calib_outcome==0)]
  # conformal_score[which(calib_outcome==0 & calib$score<cps[1])]<-0
  # #################
  #hist(conformal_score)
  #hist(conformal_score[which(calib_outcome==2)])
  #hist(conformal_score[which(calib_outcome==1)])##different dist so treat seperately
  #################its generally accurate but sometimes not and sometimes its very wrong
  alpha<-0.1
  #n_cauc<-nrow(calib)
  #q_aa0<-q_cauc0<-ceiling((length(which(calib_outcome==0))+1)*(1-alpha))/length(which(calib_outcome==0))
  q_aa1<-q_cauc1<-ceiling((length(which(calib_outcome==1))+1)*(1-alpha))/length(which(calib_outcome==1))
  q_aa2<-q_cauc2<-ceiling((length(which(calib_outcome==2))+1)*(1-alpha))/length(which(calib_outcome==2))
  #s_q_a0<-s_q_c0<-quantile(conformal_score[which(calib_outcome==0)],q_cauc0)
  s_q_a1<-s_q_c1<-quantile(conformal_score[which(calib_outcome==1)],q_cauc1)
  s_q_a2<-s_q_c2<-quantile(conformal_score[which(calib_outcome==2)],q_cauc2)
  # alpha<-0.1
  # n_cauc<-length(which(out_calib$race=="Caucasian"))
  # q_cauc<-ceiling((n_cauc+1)*(1-alpha))/n_cauc
  # s_q_c<-quantile(conformal_score[which(out_calib$race=="Caucasian")],q_cauc)
  #####ok thats our quantile value for caucasions
  ######ok now compute the conformal score globally for all test set in 3 settings
  #conformal_score_test_0<-mapply(function(i){return(if(test$score[i]<cps[1]){0}else{abs(cps[1]-test$score[i])})},1:nrow(test))
  conformal_score_test_1<-mapply(function(i){return(if(test$score[i]<cps[1]){0}else{abs(cps[1]-test$score[i])})},1:nrow(test))
  conformal_score_test_2<-mapply(function(i){return(if(test$score[i]>cps[1]){0}else{abs(cps[1]-test$score[i])})},1:nrow(test))
  #conformal_score_test_1<-mapply(function(i){return(if(test$score[i]<cps[2] & test$score[i]>cps[1]){0}else{min(abs(cps[2]-test$score[i]),abs(cps[1]-test$score[i]))})},1:nrow(test))
  ####################asign classes
  preds<-mat.or.vec(nrow(test),2)
  #preds[which(conformal_score_test_0<s_q_c0),1]<-1
  preds[which(conformal_score_test_1<s_q_c1),1]<-1
  preds[which(conformal_score_test_2<s_q_c2),2]<-1
  # table(apply(preds,1,sum))###nothing is cerain! almost all of them have all 3 and some have 2
  # table(preds[,1])
  # table(preds[,2])
  # table(preds[,3])###ok so every conformal score could be the middle one, because of the alpha choice set to 0.2 for now
  summary(lm(apply(preds,1,sum)~test$race))[4]$coefficients[2,]->reg_cauc_vs_black
  c(mean(apply(preds,1,sum)[which(test$race=="Caucasian")]),
  mean(apply(preds,1,sum)[which(test$race!="Caucasian")]))->un_race###more uncertain for blacks here
  ##################
  test_outcome<-rep(2,nrow(test))
  test_outcome[which(test$event==1)]<-1
  # test_outcome[which(test$event==1&test$time<2)]<-1
  # test_outcome[which(test$event==1&test$time<=1)]<-2
  #table(test_outcome)
  ###################
  #length(which(test_outcome==0 & preds[,1]==1))/length(which(test_outcome==0))
  #length(which(test_outcome==1 & preds[,1]==1))/length(which(test_outcome==1))
  #length(which(test_outcome==2 & preds[,2]==1))/length(which(test_outcome==2))###more wrong for early offenders
  ######################
  coverage<-(#length(which(test_outcome==0 & preds[,1]==1))+
    length(which(test_outcome==1 & preds[,1]==1))+
    length(which(test_outcome==2 & preds[,2]==1)))/length(test_outcome)##correct coverage
  #summary(lm(apply(preds,1,sum)~as.factor(test_outcome)))
  summary(lm(apply(preds,1,sum)~as.factor(test_outcome)))[4]$coefficients[2,]->reg_outcome
  #summary(lm(apply(preds,1,sum)~as.factor(test_outcome)+as.factor(test$race)))###indep of outcome uncerta
  c(#mean(apply(preds,1,sum)[which(test_outcome==0)]),##high uncertainty when dont re-offend
  mean(apply(preds,1,sum)[which(test_outcome==1)]),##moderate uncertainty for medium risk
  mean(apply(preds,1,sum)[which(test_outcome==2)]))->un_out##lower uncertainty for high risk offenders
  ########################################
  #####now we have 4 classes for our test data - certain and predict 1, uncertain and predict 1, certain and predict 2, uncertain and predict 2
  #####lets asign these 4 classes
  apply(preds,1,sum)->unc
  rep(2,length(unc))->pe
  pe[which(test$score<cps)]<-1
  table(data.frame(pe,unc))
  classes<-rep(NA,length(unc))
  classes[which(pe==1 & unc==1)]<-1
  classes[which(pe==1 & unc==2)]<-2
  classes[which(pe==2 & unc==2)]<-3
  classes[which(pe==2 & unc==1)]<-4
  table(classes)
  ###now survival analysis
  fit <- survfit(Surv(time,event)~(classes),data=test)
  print(surv_pvalue(fit))
  ###################
  ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = TRUE,         # show confidence intervals for
    # point estimaes of survival curves.
    #conf.int.style = "step",  # customize style of confidence intervals
    xlab = "Time in years",   # customize X axis label
    xlim=c(0,2),
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    risk.table = "abs_pct",  # absolute number and percentage at risk.
    risk.table.y.text.col = T,# colour risk table text annotations.
    risk.table.y.text = FALSE,# show bars instead of names in text annotations
    # in legend of risk table.
    ncensor.plot = TRUE,      # plot the number of censored subjects at time t
    surv.median.line = "hv",  # add the median survival pointer.
    legend.labs =
      c("low uncert","low cert","high cert","high uncert"
        ),    # change legend labels.
    palette =
      c("red","pink",
        "#E7B800", "#2E9FDF") # custom color palettes.
  )
  ####so output the uncertainty for blacks and whites and by outcome and the association and the coverage
  ####do coverage seperately blacks and whites too
  c(cps,coverage,un_race,un_out,reg_cauc_vs_black,reg_outcome)->out
  out->train_calib_test[i,]}
#################
# colnames(train_calib_test)<-c("lower_cut_off","upper_cut_off","white_uc","black_uc","0_uc","1_uc","2_uc",
#                               "est_r","se_r","t.val_r","p_r","est_o","se_o","t.val_o","p_o")
colnames(train_calib_test)<-c("cut_off","coverage","white_uc","black_uc","1_uc","2_uc",
                              "est_r","se_r","t.val_r","p_r","est_o","se_o","t.val_o","p_o")

head(train_calib_test)
hist(train_calib_test[,1])##cool
hist(train_calib_test[,2])##cool
mean(train_calib_test[,2])###coverage is pretty steady
boxplot(train_calib_test[,3],train_calib_test[,4])##always more uncertain for blacks
wilcox.test(train_calib_test[,3],train_calib_test[,4])##
boxplot(train_calib_test[,5],train_calib_test[,6])##always more uncertain for blacks
wilcox.test(train_calib_test[,5],train_calib_test[,6])##always more uncertain for blacks
hist(train_calib_test[,10])
length(which(train_calib_test[,10]<0.05))/b
hist(train_calib_test[,14])
length(which(train_calib_test[,14]<0.05))/b
#####
#######now look at 3 groups
##############################
###########################
b<-30
train_calib_test<-mat.or.vec(b,30)
for(i in 1:b){
  print(i)
  length(aa)+length(caus)###6005 samples
  #########so now we need three groups - 2000 training points, 2000 calibration test, remaining 2k test
  samp<-sample(nrow(all3),4000,replace = F)
  all3[samp[1:2000],]->train
  all3[samp[2001:4000],]->calib
  all3[-samp,]->test
  #
  # table(train$race)
  # table(test$race)
  # table(calib$race)
  #################
  ###first use training data to determine cut offs
  boxplot(all3$score[which(all3$event==1)],all3$score[which(all3$event==0)])##score is lower if no offence, so high score more risk
  surv_cutpoint(train,time="time",event = "event",variables ="score")->a1
  train[which(train$score>as.numeric(a1$cutpoint)),]->train2
  #train->train2
  surv_cutpoint(train2,time="time",event = "event",variables ="score")->a2
  unlist(c(as.numeric(a1$cutpoint[1]),as.numeric(a2$cutpoint[1])))->cps
  #cps<-as.numeric(a1$cutpoint[1])
  pred_train<-rep(2,nrow(train))
  pred_train[which(as.numeric(train$score)<cps[2])]<-1
  pred_train[which(train$score<cps[1])]<-0
  table(pred_train)
  # table(train$race[which(pred_train==0)])
  # table(train$race[which(pred_train==1)])
  # table(train$race[which(pred_train==2)])
  ####score of 2 is highest risk, score of 0 is lowest
  fit <- survfit(Surv(time,event)~pred_train,data=train)
  print(surv_pvalue(fit))
  ###################
  ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = TRUE,         # show confidence intervals for
    # point estimaes of survival curves.
    #conf.int.style = "step",  # customize style of confidence intervals
    xlab = "Time in years",   # customize X axis label
    xlim=c(0,2),
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    risk.table = "abs_pct",  # absolute number and percentage at risk.
    risk.table.y.text.col = T,# colour risk table text annotations.
    risk.table.y.text = FALSE,# show bars instead of names in text annotations
    # in legend of risk table.
    ncensor.plot = TRUE,      # plot the number of censored subjects at time t
    surv.median.line = "hv",  # add the median survival pointer.
    legend.labs =
      c("low","mid","high"),    # change legend labels.
    palette =
      c("pink","#E7B800", "#2E9FDF") # custom color palettes.
  )
  #dev.off()
  #################so we have our 3 groups selected via our cutoff and our 3 outcomes
  #################look at calib data asign outcomes
  calib_outcome<-rep(2,nrow(calib))
  calib_outcome[which(calib$event==0)]<-0
  calib_outcome[which(calib$event==1&calib$time<2)]<-1
  calib_outcome[which(calib$event==1&calib$time<=1)]<-2
  #calib_outcome[which(calib$event==1)]<-1
  ##outcome of 2 is reoffend in 1st year, outcome of 1 is in 2nd year, 0 is no reoffend
  #table(calib_outcome)
  # table(calib$race[which(calib_outcome==0)])
  # table(calib$race[which(calib_outcome==1)])
  # table(calib$race[which(calib_outcome==2)])
  # #################
  #################now define the conformal score
  #cps
  #calib$score
  #hist(abs(cps[2]-calib$score)[which(calib_outcome==0)])
  conformal_score<-rep(NA,nrow(calib))
  ##conformal score expecting outcome 2=reoffend in 1 year= 0 if > thresh and score-thresh else
  conformal_score[which(calib_outcome==2)]<-abs(cps[2]-calib$score)[which(calib_outcome==2)]
  conformal_score[which(calib_outcome==2 & calib$score>cps[2])]<-0
  ##same for outcome 0 but with earlier thresh
  conformal_score[which(calib_outcome==0)]<-abs(cps[1]-calib$score)[which(calib_outcome==0)]
  conformal_score[which(calib_outcome==0 & calib$score<cps[1])]<-0
  ###for outcome 1 need more to pick the minimum distance from either thresh
  conformal_score[which(calib_outcome==1)]<-apply(cbind(abs(cps[2]-calib$score),abs(cps[1]-calib$score)),1,min)[which(calib_outcome==1)]
  conformal_score[which(calib_outcome==1 & calib$score<cps[2] & calib$score>cps[1])]<-0
  # #################
  # hist(conformal_score)
  # hist(conformal_score[which(calib_outcome==2)])
  # hist(conformal_score[which(calib_outcome==0)])
  # hist(conformal_score[which(calib_outcome==1)])##different dist so treat seperately
  # ################its generally accurate but sometimes not and sometimes its very wrong
  alpha<-0.1
  #n_cauc<-nrow(calib)
  q_aa0<-q_cauc0<-ceiling((length(which(calib_outcome==0))+1)*(1-alpha))/length(which(calib_outcome==0))
  q_aa1<-q_cauc1<-ceiling((length(which(calib_outcome==1))+1)*(1-alpha))/length(which(calib_outcome==1))
  q_aa2<-q_cauc2<-ceiling((length(which(calib_outcome==2))+1)*(1-alpha))/length(which(calib_outcome==2))
  s_q_a0<-s_q_c0<-quantile(conformal_score[which(calib_outcome==0)],q_cauc0)
  s_q_a1<-s_q_c1<-quantile(conformal_score[which(calib_outcome==1)],q_cauc1)
  s_q_a2<-s_q_c2<-quantile(conformal_score[which(calib_outcome==2)],q_cauc2)
  # alpha<-0.1
  # n_cauc<-length(which(out_calib$race=="Caucasian"))
  # q_cauc<-ceiling((n_cauc+1)*(1-alpha))/n_cauc
  # s_q_c<-quantile(conformal_score[which(out_calib$race=="Caucasian")],q_cauc)
  #####ok thats our quantile value for caucasions
  ######ok now compute the conformal score globally for all test set in 3 settings
  conformal_score_test_0<-mapply(function(i){return(if(test$score[i]<cps[1]){0}else{abs(cps[1]-test$score[i])})},1:nrow(test))
  conformal_score_test_1<-mapply(function(i){return(if(test$score[i]<cps[2] & test$score[i]>cps[1]){0}else{min(abs(cps[2]-test$score[i]),abs(cps[1]-test$score[i]))})},1:nrow(test))
  conformal_score_test_2<-mapply(function(i){return(if(test$score[i]>cps[2]){0}else{abs(cps[2]-test$score[i])})},1:nrow(test))
  #
  ####################asign classes
  preds<-mat.or.vec(nrow(test),3)
  preds[which(conformal_score_test_0<s_q_c0),1]<-1
  preds[which(conformal_score_test_1<s_q_c1),2]<-1
  preds[which(conformal_score_test_2<s_q_c2),3]<-1
  table(apply(preds,1,sum))###nothing is cerain! almost all of them have all 3 and some have 2
  table(preds[,1])
  table(preds[,2])
  table(preds[,3])###ok so every conformal score could be the middle one, because of the alpha choice set to 0.2 for now
  summary(lm(apply(preds,1,sum)~test$race))[4]$coefficients[2,]->reg_cauc_vs_black
  c(mean(apply(preds,1,sum)[which(test$race=="Caucasian")]),
    mean(apply(preds,1,sum)[which(test$race!="Caucasian")]))->un_race###more uncertain for blacks here
  ##################
  test_outcome<-rep(2,nrow(test))
  test_outcome[which(test$event==0)]<-0
  test_outcome[which(test$event==1&test$time<2)]<-1
  test_outcome[which(test$event==1&test$time<=1)]<-2
  #table(test_outcome)
  ###################
  length(which(test_outcome==0 & preds[,1]==1))/length(which(test_outcome==0))
  length(which(test_outcome==1 & preds[,2]==1))/length(which(test_outcome==1))
  length(which(test_outcome==2 & preds[,3]==1))/length(which(test_outcome==2))###more wrong for early offenders
  ######################
  coverage<-(length(which(test_outcome==0 & preds[,1]==1))+
    length(which(test_outcome==1 & preds[,2]==1))+
      length(which(test_outcome==2 & preds[,3]==1)))/length(test_outcome)##correct coverage
  ######also want overall accuracy 
  
  #summary(lm(apply(preds,1,sum)~as.factor(test_outcome)))
  summary(lm(apply(preds,1,sum)~as.factor(test_outcome)))[4]$coefficients[2,]->reg_outcome
  #summary(lm(apply(preds,1,sum)~as.factor(test_outcome)+as.factor(test$race)))###indep of outcome uncerta
  c(mean(apply(preds,1,sum)[which(test_outcome==0)]),##high uncertainty when dont re-offend
    mean(apply(preds,1,sum)[which(test_outcome==1)]),##moderate uncertainty for medium risk
    mean(apply(preds,1,sum)[which(test_outcome==2)]))->un_out##lower uncertainty for high risk offenders
  ########################################
  #####now we have 4 classes for our test data - certain and predict 1, uncertain and predict 1, certain and predict 2, uncertain and predict 2
  #####lets asign these 4 classes
  apply(preds,1,sum)->unc
  rep(1,length(unc))->pe
  pe[which(test$score<cps[1])]<-0
  pe[which(test$score>cps[2])]<-2
  ###compare pe with true outcome for accuracy
  length(which(test_outcome==0 & pe==0))/length(which(test_outcome==0))
  length(which(test_outcome==1 & pe==1))/length(which(test_outcome==1))
  length(which(test_outcome==2 & pe==2))/length(which(test_outcome==2))#
  ###much lower
  (length(which(test_outcome==0 & pe==0))+
    length(which(test_outcome==1 & pe==1))+
    length(which(test_outcome==2 & pe==2)))/length(test_outcome)->point_accuracy
  ####accuracy is total number of correct predictions over total number of predictions
  length(which(test_outcome==0 & pe==0))/length(which(pe==0))
  length(which(test_outcome==1 & pe==1))/length(which(pe==1))
  length(which(test_outcome==2 & pe==2))/length(which(pe==2))##basically shite for middle class
  ##focus on top and bottom class
  (length(which(test_outcome==0 & pe==0))+length(which(test_outcome==2 & pe==2)))/length(which(pe!=1))##about 62% acc
  ####
  table(data.frame(pe,unc))
  classes<-pe
  classes[which(pe==0 & unc==1)]<- -1
  classes[which(pe==2 & unc==1)]<- 3
  ###compare accuracy of certain classes and point estimates and via ethnicity
  c(length(which(test_outcome==0 & pe==0))/length(which(pe==0)),
  length(which(test_outcome==0 & classes==-1))/length(which(classes==-1)),
  length(which(test_outcome==0 & classes==0))/length(which(classes==0)))->low_risk_acc
    
  ####
  c(length(which(test_outcome==2 & pe==2))/length(which(pe==2)),
  length(which(test_outcome==2 & classes==3))/length(which(classes==3)),
  length(which(test_outcome==2 & classes==2))/length(which(classes==2)))->high_risk_acc
  
  ####ok so should output this and seperate by race
  c(length(which(test_outcome==0 & pe==0 & test$race=="Caucasian"))/length(which(pe==0 & test$race=="Caucasian")),
  length(which(test_outcome==0 & classes==-1 & test$race=="Caucasian"))/length(which(classes==-1& test$race=="Caucasian")),
  length(which(test_outcome==0 & classes==0 & test$race=="Caucasian"))/length(which(classes==0& test$race=="Caucasian")))->low_risk_white
  ###
  c(length(which(test_outcome==0 & pe==0 & test$race!="Caucasian"))/length(which(pe==0 & test$race!="Caucasian")),
  length(which(test_outcome==0 & classes==-1 & test$race!="Caucasian"))/length(which(classes==-1& test$race!="Caucasian")),
  length(which(test_outcome==0 & classes==0 & test$race!="Caucasian"))/length(which(classes==0& test$race!="Caucasian")))->low_risk_black
  ###
  c(length(which(test_outcome==2 & pe==2 & test$race=="Caucasian"))/length(which(pe==2 & test$race=="Caucasian")),
  length(which(test_outcome==2 & classes==3 & test$race=="Caucasian"))/length(which(classes==3& test$race=="Caucasian")),
  length(which(test_outcome==2 & classes==2 & test$race=="Caucasian"))/length(which(classes==2& test$race=="Caucasian")))->high_risk_white
  ###
  c(length(which(test_outcome==2 & pe==2 & test$race!="Caucasian"))/length(which(pe==2 & test$race!="Caucasian")),
  length(which(test_outcome==2 & classes==3 & test$race!="Caucasian"))/length(which(classes==3& test$race!="Caucasian")),
  length(which(test_outcome==2 & classes==2 & test$race!="Caucasian"))/length(which(classes==2& test$race!="Caucasian")))->high_risk_black
  c(low_risk_acc,low_risk_black,low_risk_white,high_risk_acc,high_risk_black,high_risk_white)->acc_break_down
  table(classes)
  ###now survival analysis
  fit <- survfit(Surv(time,event)~(classes),data=test)
  print(surv_pvalue(fit))
  # table(test$race[which(pe==0)])/length(which(pe==0))
  # table(test$race[which(pe==1)])/length(which(pe==1))
  # table(test$race[which(pe==2)])/length(which(pe==2))
  # table(test$race[which(classes== -1)])/length(which(classes== -1))
  # table(test$race[which(classes== 0)])/length(which(classes== 0))
  # table(test$race[which(classes== 1)])/length(which(classes== 1))
  # table(test$race[which(classes== 2)])/length(which(classes== 2))
  # table(test$race[which(classes== 3)])/length(which(classes== 3))
  # ########
  # (table(test$race[which(classes== 0)])+
  # table(test$race[which(classes== 1)])+
  # table(test$race[which(classes== 2)]))/(length(which(classes== 1))+length(which(classes== 0))+length(which(classes== 2)))
  ########
  # c(table(test$race[which(classes== 3)])/table(test$race[which(pe==2)]),
  # table(test$race[which(classes== -1)])/table(test$race[which(pe==0)]))->per_certain
  ###################
  # pdf("surv_curv_classes.pdf",height = 6,width = 5)
  # ggsurvplot(
  #   fit,                     # survfit object with calculated statistics.
  #   pval = TRUE,             # show p-value of log-rank test.
  #   conf.int = TRUE,         # show confidence intervals for
  #   # point estimaes of survival curves.
  #   #conf.int.style = "step",  # customize style of confidence intervals
  #   xlab = "Time in years",   # customize X axis label
  #   xlim=c(0,2),
  #   ggtheme = theme_light(), # customize plot and risk table with a theme.
  #   risk.table = "abs_pct",  # absolute number and percentage at risk.
  #   risk.table.y.text.col = T,# colour risk table text annotations.
  #   risk.table.y.text = FALSE,# show bars instead of names in text annotations
  #   # in legend of risk table.
  #  # ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  #   surv.median.line = "hv",  # add the median survival pointer.
  #   legend.labs =
  #     c("low cert","low uncert","mid uncert","high uncert","high cert"
  #     ),    # change legend labels.
  #   palette =
  #     c("red","pink","green",
  #       "#E7B800", "#2E9FDF") # custom color palettes.
  # )
  # dev.off()
  # # ######################
  # fit.p <- survfit(Surv(time,event)~pe,data=test)
  # print(surv_pvalue(fit.p))
  # pdf("surv_curv_point_est.pdf",height = 6,width = 5)
  # ggsurvplot(
  #   fit.p,                     # survfit object with calculated statistics.
  #   pval = TRUE,             # show p-value of log-rank test.
  #   conf.int = TRUE,         # show confidence intervals for
  #   # point estimaes of survival curves.
  #   #conf.int.style = "step",  # customize style of confidence intervals
  #   xlab = "Time in years",   # customize X axis label
  #   xlim=c(0,2),
  #   ggtheme = theme_light(), # customize plot and risk table with a theme.
  #   risk.table = "abs_pct",  # absolute number and percentage at risk.
  #   risk.table.y.text.col = T,# colour risk table text annotations.
  #   risk.table.y.text = FALSE,# show bars instead of names in text annotations
  #   # in legend of risk table.
  # #  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  #   surv.median.line = "hv",  # add the median survival pointer.
  #   legend.labs =
  #     c("low","high","mid"
  #     ),    # change legend labels.
  #   palette =
  #     c("pink",
  #       "#E7B800", "#2E9FDF") # custom color palettes.
  # )
  # dev.off()
  ##
  # test[which(test$race=="Caucasian"),]->t_c
  # data.frame(t_c,classes[which(test$race=="Caucasian")])->t_c
  # colnames(t_c)[6]<-"class"
  # test[-which(test$race=="Caucasian"),]->t_a
  # data.frame(t_a,classes[-which(test$race=="Caucasian")])->t_a
  # colnames(t_a)[6]<-"class"
  # fit3 <- survfit(Surv(time,event)~class,data=t_a)
  # print(surv_pvalue(fit3))
  # ###################
  # ggsurvplot(
  #   fit2,                     # survfit object with calculated statistics.
  #   pval = TRUE,             # show p-value of log-rank test.
  #   conf.int = TRUE,         # show confidence intervals for
  #   # point estimaes of survival curves.
  #   #conf.int.style = "step",  # customize style of confidence intervals
  #   xlab = "Time in years",   # customize X axis label
  #   xlim=c(0,2),
  #   ggtheme = theme_light(), # customize plot and risk table with a theme.
  #   risk.table = "abs_pct",  # absolute number and percentage at risk.
  #   risk.table.y.text.col = T,# colour risk table text annotations.
  #   risk.table.y.text = FALSE,# show bars instead of names in text annotations
  #   # in legend of risk table.
  #   ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  #   surv.median.line = "hv",  # add the median survival pointer.
  #   legend.labs =
  #     c("high cert","high uncert","mid uncert","low uncert","low cert"
  #     ),    # change legend labels.
  #   palette =
  #     c("red","pink","green",
  #       "#E7B800", "#2E9FDF") # custom color palettes.
  # )
  # ####
  # ggsurvplot(
  #   fit3,                     # survfit object with calculated statistics.
  #   pval = TRUE,             # show p-value of log-rank test.
  #   conf.int = TRUE,         # show confidence intervals for
  #   # point estimaes of survival curves.
  #   #conf.int.style = "step",  # customize style of confidence intervals
  #   xlab = "Time in years",   # customize X axis label
  #   xlim=c(0,2),
  #   ggtheme = theme_light(), # customize plot and risk table with a theme.
  #   risk.table = "abs_pct",  # absolute number and percentage at risk.
  #   risk.table.y.text.col = T,# colour risk table text annotations.
  #   risk.table.y.text = FALSE,# show bars instead of names in text annotations
  #   # in legend of risk table.
  #   ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  #   surv.median.line = "hv",  # add the median survival pointer.
  #   legend.labs =
  #     c("high cert","high uncert","mid uncert","low uncert","low cert"
  #     ),    # change legend labels.
  #   palette =
  #     c("red","pink","green",
  #       "#E7B800", "#2E9FDF") # custom color palettes.
  # )
  # 
  ####so output the uncertainty for blacks and whites and by outcome and the association and the coverage
  ####do coverage seperately blacks and whites too
  c(cps,s_q_c0,s_q_c1,s_q_c2,coverage,point_accuracy,un_race,un_out,acc_break_down)->out
  out->train_calib_test[i,]}
#################
colnames(train_calib_test)<-c("lower_cut_off","upper_cut_off","q0","q1","q2","coverage","accuracy","white_uc","black_uc","0_uc","1_uc","2_uc",
                              "low_risk_acc","low_risk_acc_c","low_risk_acc_uc","low_risk_black","low_risk_black_c","low_risk_black_uc","low_risk_white","low_risk_white_c","low_risk_white_uc",
                              "high_risk_acc","high_risk_acc_c","high_risk_acc_uc","high_risk_black","high_risk_black_c","high_risk_black_uc","high_risk_white","high_risk_white_c","high_risk_white_uc")
# colnames(train_calib_test)<-c("cut_off","coverage","white_uc","black_uc","1_uc","2_uc",
#                               "est_r","se_r","t.val_r","p_r","est_o","se_o","t.val_o","p_o")
data.frame(train_calib_test)->output
par(mfrow=c(2,2))
hist(output$coverage)
hist(output$lower_cut_off)
hist(output$upper_cut_off)
hist(output$q0)
hist(output$q1)
hist(output$q2)
boxplot(output$low_risk_acc,output$low_risk_acc_c)
boxplot(output$low_risk_black,output$low_risk_black_c)
boxplot(output$low_risk_white,output$low_risk_white_c)

#pdf("Fig1.pdf",height=6, width = 6)
par(mfrow=c(2,2))
hist(output$accuracy,main="Accuracy of COMPAS score risk groups",xlab="Accuracy",col="lightblue");abline(v=1/3,col="red")
hist(output$coverage,main="Coverage of Prediction Sets",xlab="Coverage",col="lightblue");abline(v=0.9,col="red")
hist(output$lower_cut_off,col="lightblue",xlab="c_low",main="c_low distribution")
hist(output$upper_cut_off,col="lightblue",xlab="c_high",main="c_high distribution")
##
hist(output$q0,col="lightblue",xlab="q^low",main="q^low")
hist(output$q1,col="lightblue",xlab="q^med",main="q^med")
hist(output$q2,col="lightblue",xlab="q^high",main="q^high")
#dev.off()
pdf("Fig2.pdf",height=6, width = 6)
par(mfrow=c(2,2))
boxplot(output$low_risk_acc,output$low_risk_acc_c,output$low_risk_acc_uc,
        col=c("lightgreen","pink","lightblue"),names=c("Point-estimate","certain","uncertain"),
        main="Low Risk allocation accuracy",ylab="Accuracy")
wilcox.test(output$low_risk_acc,output$low_risk_acc_c)
wilcox.test(output$low_risk_acc,output$low_risk_acc_uc)
mean(output$low_risk_acc)
mean(output$low_risk_acc_c)
mean(output$low_risk_acc_uc)
###########
boxplot(output$high_risk_acc,output$high_risk_acc_c,output$high_risk_acc_uc,
        col=c("lightgreen","pink","lightblue"),names=c("Point-estimate","certain","uncertain"),
        main="High Risk allocation accuracy",ylab="Accuracy")
wilcox.test(output$high_risk_acc,output$high_risk_acc_c)
wilcox.test(output$high_risk_acc,output$high_risk_acc_uc)
mean(output$high_risk_acc)
mean(output$high_risk_acc_c)
mean(output$high_risk_acc_uc)
dev.off()
####
pdf("Fig3.pdf",height=6, width = 6)
par(mfrow=c(2,2))
#dev.off()
boxplot(output$black_uc,output$white_uc,
        col=c("brown","pink"),names=c("African-American","Caucasian"),
        main="Uncertainty Based on Ethnicity",ylab="Mean uncertainty")
wilcox.test(output$black_uc,output$white_uc)
##########
boxplot(output$low_risk_black,output$low_risk_black_c,output$low_risk_black_uc,
        col=c("lightgreen","pink","lightblue"),names=c("Point-estimate","certain","uncertain"),
        main="Low Risk allocation accuracy: African-American",ylab="Accuracy")
boxplot(output$low_risk_white,output$low_risk_white_c,output$low_risk_white_uc,
        col=c("lightgreen","pink","lightblue"),names=c("Point-estimate","certain","uncertain"),
        main="Low Risk allocation accuracy: Caucasian",ylab="Accuracy")
wilcox.test(output$low_risk_white,output$low_risk_white_c,paired = T)
mean(output$low_risk_white);mean(output$low_risk_white_c)
wilcox.test(output$low_risk_black,output$low_risk_black_c,paired = T)
mean(output$low_risk_black);mean(output$low_risk_black_c)
###
boxplot(output$high_risk_black,output$high_risk_black_c,output$high_risk_black_uc,
        col=c("lightgreen","pink","lightblue"),names=c("Point-estimate","certain","uncertain"),
        main="High Risk allocation accuracy: African-American",ylab="Accuracy")
boxplot(output$high_risk_white,output$high_risk_white_c,output$high_risk_white_uc,
        col=c("lightgreen","pink","lightblue"),names=c("Point-estimate","certain","uncertain"),
        main="High Risk allocation accuracy: Caucasian",ylab="Accuracy")
wilcox.test(output$high_risk_white,output$high_risk_white_c,paired = T)
mean(output$high_risk_white);mean(output$high_risk_white_c)
wilcox.test(output$high_risk_black,output$high_risk_black_c,paired = T)
mean(output$high_risk_black);mean(output$high_risk_black_c)
#####################################
###
dev.off()
# #####
# boxplot(output$high_risk_black,output$high_risk_white,
#         output$high_risk_black_c,output$high_risk_white_c,
#         output$high_risk_black_uc,output$high_risk_white_uc,
#         col=c("green","lightgreen",
#               "red","pink",
#               "blue","lightblue"),
#         names=c("Point-estimate_AA","Point-estimate_C",
#                 "certain_AA","certain_C",
#                 "uncertain_AA","uncertain_C"),
#         main="High Risk allocation accuracy: African-American",ylab="Accuracy")
# ##
wilcox.test(output$high_risk_white,output$high_risk_black)
wilcox.test(output$high_risk_white_c,output$high_risk_black_c)
wilcox.test(output$low_risk_white,output$low_risk_black)
wilcox.test(output$low_risk_white_c,output$low_risk_black_c)
##
t.test(output$low_risk_white/output$low_risk_black,mu=1)
t.test(output$low_risk_white_c/output$low_risk_black_c,mu=1)
###
t.test(output$low_risk_white/output$low_risk_black,
            output$low_risk_white_c/output$low_risk_black_c,paired = T)
##
wilcox.test(output$low_risk_white-output$low_risk_black,
output$low_risk_white_c-output$low_risk_black_c,paired = T)
##
wilcox.test(output$high_risk_white-output$high_risk_black,
            output$high_risk_white_c-output$high_risk_black_c,paired = T)
###
boxplot(output$high_risk_white-output$high_risk_black,
            output$high_risk_white_c-output$high_risk_black_c,paired = T)
###
boxplot(output$low_risk_white-output$low_risk_black,
            output$low_risk_white_c-output$low_risk_black_c,paired = T)

dev.off()
#####
mean(output$lower_cut_off);mean(output$upper_cut_off)
mean(output$coverage)
#boxplot(output$low_risk_black,output$low_risk_black_c,output$low_risk_white,output$low_risk_white_c)
boxplot(output$low_risk_black_c/output$low_risk_black,output$low_risk_white_c/output$low_risk_white)
###
boxplot(output$low_risk_black/output$low_risk_white,output$low_risk_black_c/output$low_risk_white_c)
###
boxplot(output$high_risk_acc,output$high_risk_acc_c)
boxplot(output$high_risk_black,output$high_risk_black_c)
boxplot(output$high_risk_white,output$high_risk_white_c)
##
boxplot(output$high_risk_black/output$high_risk_white,output$high_risk_black_c/output$high_risk_white_c)
#boxplot(output$low_risk_black,output$low_risk_black_c,output$low_risk_white,output$low_risk_white_c)
boxplot(output$high_risk_black_c/output$high_risk_black,output$high_risk_white_c/output$high_risk_white)

wilcox.test(output$white_uc,output$black_uc)
###
boxplot(output$X0_uc,output$X1_uc,output$X2_uc)
####
boxplot(output[,c(12:13)])
wilcox.test(output[,c(12)],output[,13])
##
boxplot(output[,c(14:15)])
wilcox.test(output[,c(14)],output[,15])
boxplot(train_calib_test[,3],train_calib_test[,4])##always more uncertain for blacks
wilcox.test(train_calib_test[,3],train_calib_test[,4])##
boxplot(train_calib_test[,5],train_calib_test[,6])##always more uncertain for blacks
wilcox.test(train_calib_test[,5],train_calib_test[,6])##always more uncertain for blacks
hist(train_calib_test[,10])
length(which(train_calib_test[,10]<0.05))/b
hist(train_calib_test[,14])
length(which(train_calib_test[,14]<0.05))/b
save.image("class_cond_conformal_surv_analysis.rd")
