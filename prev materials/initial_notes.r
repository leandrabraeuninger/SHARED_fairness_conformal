setwd("~/TCGA_conformal/")
rm(list=ls())
list.dirs()->ds
ds[-c(which(ds=="."),which(ds=="./Misc"))]->ds
setwd(ds[1])
read.delim("ba295155-272e-43eb-9d6a-e4c9c392e68b.rna_seq.augmented_star_gene_counts.tsv",skip = 1)->x
head(x)
####extract  the unstranded and the gene id and name
x[,1:3]->annos
x$unstranded->x
i<-2
mapply(function(i){
  setwd("../")
  setwd(ds[i])
  list.files()->l
  if(length(l)!=1){l<-l[-which(l=="annotations.txt")]}
  read.delim(l,skip=1)->x2
  print(c(i,unique(x2[,1]==annos[,1])))
  return(x2$unstranded)
}, 2:length(ds))->y
ds[26]
dim(y)
y[1:10,1:10]
cbind(x,y)->data
rownames(data)<-annos$gene_id
mapply(function(i){return(unlist(strsplit(ds[i],split="./"))[2])},1:length(ds))->colnames(data)
data[-c(1:4),]->data
data[1:10,1:10]
head(annos)
annos[-c(1:4),]->annos
########################################ok thats all the data
dim(data);dim(annos)
setwd("~/TCGA_conformal/Misc/")
read.delim("clinical.tsv")->clin
read.delim("gdc_sample_sheet.2023-05-12.tsv")->mf
head(mf)
match(colnames(data),mf$File.ID)->m
cbind(mf$File.ID[m],mf$Case.ID[m])->ids
head(ids)
colnames(data)[1:6]
match(ids[,2],clin$case_submitter_id)->m2
length(which(is.na(m2)))
clin[m2,]->clin
cbind(clin$case_submitter_id,ids,colnames(data))[1:10,]###all match
colnames(data2)<-ids[,1]
####
length(unique(clin$case_submitter_id))###ok so we use this and only this
match(unique(clin$case_submitter_id),ids[,2])->mm
data[,mm]->data2
match(unique(clin$case_submitter_id),clin$case_submitter_id)->mmm
clin[mmm,]->clin2
cbind(clin2$case_submitter_id,colnames(data2))[1:10,]##
unique(clin2$case_submitter_id==colnames(data2))##good
length(unique(colnames(data2)))
length(unique(colnames(data)))
#
######################################
data2[1:10,1:10]
boxplot(data2[,1:10])
head(clin2)
library('DESeq2')
rownames(clin2)<-clin2$case_submitter_id
dds <- DESeqDataSetFromMatrix(countData = data2, colData = clin2,design = ~1)
dds<-estimateSizeFactors(dds)
x <- counts(dds, normalized=TRUE)
x[1:10,1:10]
boxplot(x[,1:10])
clin<-clin2
data<-x
#save(data,clin,dds,file="TCGA-BRCA_LB.Rd")
load("TCGA-BRCA_LB.Rd")
######################################
colnames(clin)
table(clin$vital_status)###so ~1/6 died, 1 unknown
data[,-which(clin$vital_status=="'--")]->data
clin[-which(clin$vital_status=="'--"),]->clin
head(clin)
table(clin$days_to_diagnosis)##ok Dx is zero days
cbind(as.numeric(clin$days_to_death),clin$days_to_death)##so as.numeric works
cbind(as.numeric(clin$days_to_death),clin$vital_status)###death recorded
unique(as.numeric(clin$days_to_death)[which(clin$vital_status=="Alive")])###all alive have no time til death
which(is.na(as.numeric(clin$days_to_death)[which(clin$vital_status=="Dead")]))###so one of the dead have no time until death, remove them
which(clin$days_to_death=="'--" & clin$vital_status=="Dead")->rm
data[,-rm]->data
clin[-rm,]->clin
table(clin$cause_of_death)###no known cause of death - ass. all cause mortality
cbind(as.numeric(clin$days_to_last_follow_up),clin$days_to_last_follow_up)###works
cbind(as.numeric(clin$days_to_last_follow_up),clin$vital_status)###ok
cbind(as.numeric(clin$days_to_last_follow_up),clin$vital_status,clin$days_to_death)
plot(as.numeric(clin$days_to_death)[which(clin$vital_status=="Dead")],
  as.numeric(clin$days_to_last_follow_up)[which(clin$vital_status=="Dead")],xlab="death",ylab="last fu");abline(a=0,b=1)##good death always after last fu
############################################################survival analysis
library('survival')
event<-rep(NA,ncol(data))
event[which(clin$vital_status=="Dead")]<-1
event[which(clin$vital_status=="Alive")]<-0
table(event)##good
surv_time<-rep(NA,ncol(data))
surv_time[which(clin$vital_status=="Dead")]<-as.numeric(clin$days_to_death)[which(clin$vital_status=="Dead")]
surv_time[which(clin$vital_status=="Alive")]<-as.numeric(clin$days_to_last_follow_up)[which(clin$vital_status=="Alive")]
table(surv_time)###ok so one of the surv times is negative
length(which(is.na(surv_time)))###none are missing
which(surv_time<0)->rm
data[,-rm]->data
clin[-rm,]->clin
surv_time[-rm]->surv_time
event[-rm]->event
############################################################
head(clin)
r<-as.factor(clin$age_at_index)
#table(r[-which(r=="'--")])
table(r)
#summary(coxph(Surv(surv_time[-which(r=="'--")],event[-which(r=="'--")])~r[-which(r=="'--")]))
summary(coxph(Surv(surv_time,event)~as.numeric(r)))###so age is significant
############################################################
head(clin)
r<-as.factor(clin$ethnicity)
#table(r[-which(r=="'--")])
table(r)
#summary(coxph(Surv(surv_time[-which(r=="'--")],event[-which(r=="'--")])~r[-which(r=="'--")]))
summary(coxph(Surv(surv_time,event)~as.factor(r)))###so age is significant
plot(survfit(Surv(surv_time/365.25,event)~r),xlab="Time (years)",
     col=c("lightblue","pink","black"),ylab="survival probability",lwd=2)###the small number of hispanics do better
##############################################################
head(clin)
r<-as.factor(clin$race)
#table(r[-which(r=="'--")])
table(r)
unique(r)->u
u
which(r==u[1]|r==u[2])->f
#summary(coxph(Surv(surv_time[-which(r=="'--")],event[-which(r=="'--")])~r[-which(r=="'--")]))
summary(coxph(Surv(surv_time[f],event[f])~as.character(r[f])))###so age is significant
###############################################################
###stage is interesting
head(clin)
r<-as.factor(clin$ajcc_pathologic_stage)
#table(r[-which(r=="'--")])
table(r)
#######fo ease lets aggregate into stage 1, 2, 3 and 4
stage<-rep(NA,length(r))
# stage[which(r=="Stage I")]<-1
# stage[which(r=="Stage IA")]<-1+1/3
# stage[which(r=="Stage IB")]<-1+2/3
# stage[which(r=="Stage II")]<-2
# stage[which(r=="Stage IIA")]<-2+1/3
# stage[which(r=="Stage IIB")]<-2+2/3
# stage[which(r=="Stage III")]<-3
# stage[which(r=="Stage IIIA")]<-3+1/4
# stage[which(r=="Stage IIIB")]<-3+2/4
# stage[which(r=="Stage IIIC")]<-3+3/4
# stage[which(r=="Stage IV")]<-4
###########
stage[which(r=="Stage I")]<-1
stage[which(r=="Stage IA")]<-1
stage[which(r=="Stage IB")]<-1
stage[which(r=="Stage II")]<-2
stage[which(r=="Stage IIA")]<-2
stage[which(r=="Stage IIB")]<-2
stage[which(r=="Stage III")]<-3
stage[which(r=="Stage IIIA")]<-3
stage[which(r=="Stage IIIB")]<-3
stage[which(r=="Stage IIIC")]<-3
stage[which(r=="Stage IV")]<-4

#summary(coxph(Surv(surv_time[-which(r=="'--")],event[-which(r=="'--")])~r[-which(r=="'--")]))
summary(coxph(Surv(surv_time,event)~stage+as.numeric(clin$age_at_diagnosis)))###so stage sig
summary(coxph(Surv(surv_time,event)~as.factor(stage)))###1&2 similar
plot(survfit(Surv(surv_time/365.25,event)~stage),xlab="Time (years)",
     col=c("lightblue","blue","red","darkred"),ylab="survival probability",lwd=2)###looks correct
#####################
#####################
prolif<-c("MKI67","AURKA","BIRC5","CCNB1","MYBL2")
invasion<-c("MMP11","CTSV")
HER2<-c("GRB7","ERBB2")
oest<-c("ESR1","PGR","BCL2","SCUBE2")
g<-"GSTM1"
b<-"BAG1"
cd<-"CD68"
c(prolif,invasion,HER2,oest,g,b,cd)->og
refs<-c("GAPDH","RPLP0","GUSB","ACTB","TFRC")
########ok so the formula uses RT-qPCR normalised the the refs so lets just normalise by the mean of the refs
0.47*apply(data[match(HER2,annos$gene_name),],2,sum)-
  0.34*apply(data[match(oest,annos$gene_name),],2,sum)+
  1.04*apply(data[match(prolif,annos$gene_name),],2,sum)+
  0.1*apply(data[match(invasion,annos$gene_name),],2,sum)+
  0.05*data[match(cd,annos$gene_name),]-
  0.08*data[match(g,annos$gene_name),]-
  0.07*data[match(b,annos$gene_name),]->onco_un
apply(data[match(refs,annos$gene_name),],2,sum)->norms
hist(onco_un)
hist(norms)
hist(onco_un/norms)
onco_un/norms->onco_n
#########################
summary(coxph(Surv(surv_time,event)~onco_un+norms))###1&2 similar
summary(coxph(Surv(surv_time,event)~t(data[match(og,annos$gene_name),])))###1&2 similar
annos[match(og,annos$gene_name),]
########################
read.csv("mp_genes.csv")->mp
head(mp)
mp[which(mp$correlation>0),]->mp_p
mp[which(mp$correlation<0),]->mp_n
match(mp_p$gene.name,annos$gene_name)->m
summary(coxph(Surv(surv_time,event)~t(data[m[-which(is.na(m))],])))###1&2 similar
######
match(mp_n$gene.name,annos$gene_name)->m2
summary(coxph(Surv(surv_time,event)~t(data[m2[-which(is.na(m2))],])))###1&2 similar
#######################################
apply(data[m2[-which(is.na(m2))],],2,mean)/apply(data[m[-which(is.na(m))],],2,mean)->q
summary(coxph(Surv(surv_time,event)~q))###1&2 similar
summary(coxph(Surv(surv_time,event)~apply(data[m2[-which(is.na(m2))],],2,mean)+apply(data[m[-which(is.na(m))],],2,mean)))###1&2 similar
table(clin$treatment_type)
summary(coxph(Surv(surv_time,event)~apply(data[m2[-which(is.na(m2))],],2,mean)+apply(data[m[-which(is.na(m))],],2,mean)
              +as.numeric(clin$age_at_diagnosis)+stage+as.factor(clin$gender)))###1&2 similar
####
summary(coxph(Surv(surv_time,event)~q
              +as.numeric(clin$age_at_diagnosis)+stage+as.factor(clin$gender)))###1&2 similar
#########
library('survminer')
#########################################
data.frame(surv_time/365.25,event,q,apply(data[m2[-which(is.na(m2))],],2,mean),
           apply(data[m[-which(is.na(m))],],2,mean),stage,as.numeric(clin$age_at_diagnosis),as.factor(clin$gender))->all
head(all)
colnames(all)<-c("time","event","score","pos_score","neg_score","stage","age","sex")
#surv_cutpoint(all,time="time",event = "event",variables =c("score","stage","age"))->a1
surv_cutpoint(all,time="time",event = "event",variables =c("pos_score","neg_score"))->a1
a1
#hist(q)
# all[which(all$score>as.numeric(a1$cutpoint)),]->all2
# surv_cutpoint(all2,time="time",event = "event",variables ="score")->a2
# unlist(c(as.numeric(a1$cutpoint[1]),as.numeric(a2$cutpoint[1])))->cps
# pred<-rep(2,nrow(all))
# pred[which(as.numeric(all$score)<cps[2])]<-1
# pred[which(all$score<cps[1])]<-0
# table(data.frame(pred,event))
as.matrix(a1$cutpoint)[,1]->cps
pred<-rep(1,nrow(all))
pred[which(all$pos_score <=cps[1] & all$neg_score>cps[2])]<-0##lowest risk
# pred[which(all$pos_score>cps[1] & all$neg_score>cps[2])]<-1##highest risk
# pred[which(all$pos_score<=cps[1] & all$neg_score<= cps[2])]<-1##highest risk
pred[which(all$pos_score>cps[1] & all$neg_score<= cps[2])]<-2##highest risk
table(pred)
table(data.frame(pred,event))
summary(coxph(Surv(surv_time,event)~pred
              +as.numeric(clin$age_at_diagnosis)+stage+as.factor(clin$gender)))###1&2 similar
###sig indep of the others
#######
data.frame(all,pred)->all2
colnames(all2)
fit <- survfit(Surv(time,event)~pred,data=all2)
print(surv_pvalue(fit,data=all2))
##################
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for
  # point estimaes of survival curves.
  #conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in days",   # customize X axis label
  xlim=c(0,10),
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
#  legend.labs =
 #   c("low","mid"),#,"high"),    # change legend labels.
  #palette =
   # c("pink","#E7B800"),#, "#2E9FDF"), # custom color palettes.
  data=all2
)
#####################################ok so build the loop
all$time[which(all$time>10)]<-10
all$event[which(all_hist$time>10& all_hist$event==1)]<-0
hist(all$stage)
(all$stage-min(all$stage,na.rm=T))/(max(all$stage,na.rm=T)-1)->stage_n
table(stage_n)
stage_n->all$pos_score##high is bad
hist((all$age-min(all$age,na.rm = T))/(max(all$age,na.rm=T)-min(all$age,na.rm=T)))
(all$age-min(all$age,na.rm = T))/(max(all$age,na.rm=T)-min(all$age,na.rm=T))->age_n
1-age_n->all$neg_score##so high is good
##################################
b<-30
dim(all)
train_calib_test<-mat.or.vec(b,30)
for(i in 1:b){
  # samp<-sample(nrow(all),2*ceiling(nrow(all)/3),replace = F)
  # all[samp[1:ceiling(nrow(all)/3)],]->train
  # all[samp[-c(1:ceiling(nrow(all)/3))],]->calib
  samp<-sample(nrow(all),300,replace = F)
  all[samp[1:200],]->train
  all[samp[-c(1:200)],]->calib
  all[-samp,]->test
  ##
  surv_cutpoint(train,time="time",event = "event",variables =c("pos_score","neg_score"))->a1
  as.matrix(a1$cutpoint)[,1]->cps
  pred_train<-rep(1,nrow(train))
  pred_train[which(train$pos_score <=cps[1] & train$neg_score>cps[2])]<-0##lowest risk
  pred_train[which(train$pos_score>cps[1] & train$neg_score<= cps[2])]<-2##highest risk
  table(pred_train)
  table(data.frame(pred_train,train$event))
  summary(coxph(Surv(train$time,train$event)~pred_train))
                #+train$age+train$stage+train$sex))###1&2 similar
  ###sig indep of the others
  #######
  data.frame(train,pred_train)->train2
  fit <- survfit(Surv(time,event)~pred_train,data=train2)
  print(surv_pvalue(fit,data=train2))
  ##################
  ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = TRUE,         # show confidence intervals for
    # point estimaes of survival curves.
    #conf.int.style = "step",  # customize style of confidence intervals
    xlab = "Time in days",   # customize X axis label
    xlim=c(0,10),
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    risk.table = "abs_pct",  # absolute number and percentage at risk.
    risk.table.y.text.col = T,# colour risk table text annotations.
    risk.table.y.text = FALSE,# show bars instead of names in text annotations
    # in legend of risk table.
    ncensor.plot = TRUE,      # plot the number of censored subjects at time t
    surv.median.line = "hv",  # add the median survival pointer.
    #  legend.labs =
    #   c("low","mid"),#,"high"),    # change legend labels.
    #palette =
    # c("pink","#E7B800"),#, "#2E9FDF"), # custom color palettes.
    data=train2
  )
  #####
  calib_outcome<-rep(2,nrow(calib))
  calib_outcome[which(calib$event==0)]<-0
  calib_outcome[which(calib$event==1&calib$time<5)]<-1
  #calib_outcome[which(calib$event==1&calib$time<=10)]<-2
  table(calib_outcome)
  conformal_score<-rep(NA,nrow(calib))
  ##conformal score expecting outcome 2=reoffend in 1 year= 0 if > thresh and score-thresh else
  #########fix these
  # conformal_score[which(calib_outcome==2)]<-sqrt((cps[1]-calib$pos_score)[which(calib_outcome==2)]^2+(cps[2]-calib$neg_score)[which(calib_outcome==2)]^2)
  # conformal_score[which(calib_outcome==2 & calib$pos_score>cps[1] & calib$neg_score> cps[2])]<-sqrt((cps[2]-calib$neg_score)[which(calib_outcome==2& calib$pos_score>cps[1]& calib$neg_score> cps[2])]^2)
  # conformal_score[which(calib_outcome==2 & calib$neg_score<=cps[2] & calib$pos_score<=cps[1])]<-sqrt((cps[1]-calib$pos_score)[which(calib_outcome==2& calib$neg_score<=cps[2]& calib$pos_score<=cps[1])]^2)
  # conformal_score[which(calib_outcome==2 & calib$pos_score>cps[1] & calib$neg_score<= cps[2])]<-0
  # ##same for outcome 0 but with earlier thresh
  # conformal_score[which(calib_outcome==1)]<-sqrt((cps[1]-calib$pos_score)[which(calib_outcome==1)]^2+(cps[2]-calib$neg_score)[which(calib_outcome==1)]^2)
  # conformal_score[which(calib_outcome==1 & calib$pos_score>cps[1]&calib$neg_score<= cps[2])]<-apply(cbind(sqrt((cps[1]-calib$pos_score)[which(calib_outcome==1& calib$pos_score>cps[1]&calib$neg_score<= cps[2])]^2),
  #                                                                                                 sqrt((cps[2]-calib$neg_score)[which(calib_outcome==1& calib$pos_score>cps[1]&calib$neg_score<= cps[2])]^2)),1,min)
  # conformal_score[which(calib_outcome==1 & calib$pos_score<=cps[1]&calib$neg_score>cps[2])]<-apply(cbind(sqrt((cps[1]-calib$pos_score)[which(calib_outcome==1& calib$pos_score<=cps[1]&calib$neg_score> cps[2])]^2),
  #                                                                                                 sqrt((cps[2]-calib$neg_score)[which(calib_outcome==1& calib$pos_score<=cps[1]&calib$neg_score> cps[2])]^2)),1,min)
  # conformal_score[which(calib_outcome==1 & calib$pos_score<=cps[1] & calib$neg_score<= cps[2])]<-0
  # conformal_score[which(calib_outcome==1 & calib$pos_score>cps[1] & calib$neg_score>cps[2])]<-0
  # ####
  # conformal_score[which(calib_outcome==0)]<-sqrt((cps[1]-calib$pos_score)[which(calib_outcome==0)]^2+(cps[2]-calib$neg_score)[which(calib_outcome==0)]^2)
  # conformal_score[which(calib_outcome==0 & calib$pos_score<=cps[1] & calib$neg_score<cps[2])]<-sqrt((cps[2]-calib$neg_score)[which(calib_outcome==0& calib$pos_score<=cps[1] & calib$neg_score<cps[2])]^2)
  # conformal_score[which(calib_outcome==0 & calib$pos_score>cps[1]& calib$neg_score>cps[2])]<-sqrt((cps[1]-calib$pos_score)[which(calib_outcome==0& calib$pos_score>cps[1]& calib$neg_score>cps[2])]^2)
  # conformal_score[which(calib_outcome==0 & calib$pos_score<=cps[1] & calib$neg_score>cps[2])]<-0
  conformal_score_cal_0<-mapply(function(i){return(if(is.na(calib$pos_score[i]+calib$neg_score[i])){NA}else 
    if(calib$pos_score[i] <= cps[1] & calib$neg_score[i] > cps[2]){0}else 
      if(calib$pos_score[i]>cps[1] & calib$neg_score[i]<=cps[1]){sqrt((cps[1]-calib$pos_score[i])^2)} else 
        if(calib$pos_score[i]<=cps[1] & calib$neg_score[i]<cps[2]){sqrt((cps[2]-calib$neg_score[i])^2)} else{sqrt((cps[1]-calib$pos_score[i])^2+(cps[2]-calib$neg_score[i])^2)}
  )},1:nrow(calib))
  ######
  conformal_score_cal_1<-mapply(function(i){return(if(is.na(calib$pos_score[i]+calib$neg_score[i])){NA}else 
    if(calib$pos_score[i]>cps[1] & calib$neg_score[i]>cps[2]){0}else if(calib$pos_score[i]<=cps[1] & calib$neg_score[i]<=cps[2]){0}else
    if(calib$neg_score[i]>cps[2]& calib$pos_score[i]<=cps[1]){min(sqrt((cps[1]-calib$pos_score[i])^2),
                                                                sqrt((cps[2]-calib$neg_score[i])^2))}else
                                                                  if(calib$pos_score[i]>cps[1] & calib$neg_score[i]<cps[2]){min(sqrt((cps[1]-calib$pos_score[i])^2),
                                                                                                                              sqrt((cps[2]-calib$neg_score[i])^2))}else{sqrt((cps[1]-calib$pos_score[i])^2+(cps[2]-calib$neg_score[i])^2)}
  )},1:nrow(calib))
  ######
  conformal_score_cal_2<-mapply(function(i){return(if(is.na(calib$pos_score[i]+calib$neg_score[i])){NA}else 
    if(calib$pos_score[i]>cps[1] & calib$neg_score[i]<cps[2]){0}else 
    if(calib$neg_score[i]<cps[2]& calib$pos_score[i]<=cps[1]){sqrt((cps[1]-calib$pos_score[i])^2)}else
      if(calib$pos_score[i]>cps[1] & calib$neg_score[i]>cps[2]){sqrt((cps[2]-calib$neg_score[i])^2)}else{sqrt((cps[1]-calib$pos_score[i])^2+(cps[2]-calib$neg_score[i])^2)}
  )},1:nrow(calib))
  conformal_score[which(calib_outcome==2)]<-conformal_score_cal_2[which(calib_outcome==2)]
  conformal_score[which(calib_outcome==1)]<-conformal_score_cal_1[which(calib_outcome==1)]
  conformal_score[which(calib_outcome==0)]<-conformal_score_cal_0[which(calib_outcome==0)]
  hist(conformal_score)
  hist(conformal_score[which(calib_outcome==2)])
  hist(conformal_score[which(calib_outcome==0)])
  hist(conformal_score[which(calib_outcome==1)])##different dist so treat seperately
  #which(is.na(conformal_score))
  # ################its generally accurate but sometimes not and sometimes its very wrong
  alpha<-0.2
  q_aa0<-q_cauc0<-ceiling((length(which(calib_outcome==0))+1)*(1-alpha))/length(which(calib_outcome==0))
  q_aa1<-q_cauc1<-ceiling((length(which(calib_outcome==1))+1)*(1-alpha))/length(which(calib_outcome==1))
  q_aa2<-q_cauc2<-ceiling((length(which(calib_outcome==2))+1)*(1-alpha))/length(which(calib_outcome==2))
  s_q_a0<-s_q_c0<-quantile(conformal_score[which(calib_outcome==0)],q_cauc0,na.rm = T)
  s_q_a1<-s_q_c1<-quantile(conformal_score[which(calib_outcome==1)],q_cauc1,na.rm = T)
  s_q_a2<-s_q_c2<-quantile(conformal_score[which(calib_outcome==2)],q_cauc2,na.rm = T)
  ########################################################################################################
  conformal_score_test_0<-mapply(function(i){return(if(is.na(test$pos_score[i]+test$neg_score[i])){NA}else if(test$pos_score[i] <= cps[1] && test$neg_score[i] > cps[2]){0}else if(test$pos_score[i]>cps[1] && test$neg_score[i]<=cps[1]){sqrt((cps[1]-test$pos_score[i])^2)} else if(test$pos_score[i]<=cps[1] && test$neg_score[i]<cps[2]){sqrt((cps[2]-test$neg_score[i])^2)} else{sqrt((cps[1]-test$pos_score[i])^2+(cps[2]-test$neg_score[i])^2)}
    )},1:nrow(test))
  ######
  conformal_score_test_1<-mapply(function(i){return(if(is.na(test$pos_score[i]+test$neg_score[i])){NA}else if(test$pos_score[i]>cps[1] & test$neg_score[i]>cps[2]){0}else if(test$pos_score[i]<=cps[1] & test$neg_score[i]<=cps[2]){0}else
    if(test$neg_score[i]>cps[2]& test$pos_score[i]<=cps[1]){min(sqrt((cps[1]-test$pos_score[i])^2),
                                                                 sqrt((cps[2]-test$neg_score[i])^2))}else
      if(test$pos_score[i]>cps[1] & test$neg_score[i]<cps[2]){min(sqrt((cps[1]-test$pos_score[i])^2),
                                                                    sqrt((cps[2]-test$neg_score[i])^2))}else{sqrt((cps[1]-test$pos_score[i])^2+(cps[2]-test$neg_score[i])^2)}
  )},1:nrow(test))
  ######
  conformal_score_test_2<-mapply(function(i){return(if(is.na(test$pos_score[i]+test$neg_score[i])){NA}else if(test$pos_score[i]>cps[1] & test$neg_score[i]<cps[2]){0}else 
    if(test$neg_score[i]<cps[2]& test$pos_score[i]<=cps[1]){sqrt((cps[1]-test$pos_score[i])^2)}else
      if(test$pos_score[i]>cps[1] & test$neg_score[i]>cps[2]){sqrt((cps[2]-test$neg_score[i])^2)}else{sqrt((cps[1]-test$pos_score[i])^2+(cps[2]-test$neg_score[i])^2)}
  )},1:nrow(test))
  ####################asign classes
  hist(conformal_score_test_0)
  hist(conformal_score_test_1)
  hist(conformal_score_test_2)
  ######
  preds<-mat.or.vec(nrow(test),3)
  preds[which(conformal_score_test_0<s_q_c0),1]<-1
  preds[which(conformal_score_test_1<s_q_c1),2]<-1
  preds[which(conformal_score_test_2<s_q_c2),3]<-1
  table(apply(preds,1,sum))###nothing is cerain! almost all of them have all 3 and some have 2
  table(preds[,1])
  table(preds[,2])
  table(preds[,3])###ok so every conformal score could be the middle one, because of the alpha choice set to 0.2 for now
  ###############
  test_outcome<-rep(2,nrow(test))
  test_outcome[which(test$event==0)]<-0
  test_outcome[which(test$event==1&test$time<5)]<-1
  table(test_outcome)
  ###################
  length(which(test_outcome==0 & preds[,1]==1))/length(which(test_outcome==0))
  length(which(test_outcome==1 & preds[,2]==1))/length(which(test_outcome==1))
  length(which(test_outcome==2 & preds[,3]==1))/length(which(test_outcome==2))
  ######################
  coverage<-(length(which(test_outcome==0 & preds[,1]==1))+
               length(which(test_outcome==1 & preds[,2]==1))+
               length(which(test_outcome==2 & preds[,3]==1)))/length(test_outcome)##correct coverage
  ######also want overall accuracy
  c(mean(apply(preds,1,sum)[which(test_outcome==0)]),##high uncertainty when dont re-offend
    mean(apply(preds,1,sum)[which(test_outcome==1)]),##moderate uncertainty for medium risk
    mean(apply(preds,1,sum)[which(test_outcome==2)]))->un_out##lower uncertainty for high risk offenders
  ########################################
  #####now we have 4 classes for our test data - certain and predict 1, uncertain and predict 1, certain and predict 2, uncertain and predict 2
  #####lets asign these 4 classes
  apply(preds,1,sum)->unc
  rep(1,length(unc))->pe
  pe[which(test$pos_score>cps[1] & test$neg_score<cps[2])]<-0
  pe[which(test$pos_score<cps[1] & test$neg_score>cps[2])]<-2
  table(pe)
  ###compare pe with true outcome for accuracy
  length(which(test_outcome==0 & pe==0))/length(which(test_outcome==0))
  length(which(test_outcome==1 & pe==1))/length(which(test_outcome==1))
  length(which(test_outcome==2 & pe==2))/length(which(test_outcome==2))#
  ###much lower
  (length(which(test_outcome==0 & pe==0))+
      length(which(test_outcome==1 & pe==1))+
      length(which(test_outcome==2 & pe==2)))/length(test_outcome)->point_accuracy
  ####
  table(data.frame(pe,unc))
  classes<-pe
  classes[which(pe==0 & unc==1)]<- -1
  classes[which(pe==2 & unc==1)]<- 3
  ###compare accuracy of certain classes and point estimates and via ethnicity
  c(length(which(test_outcome==0 & pe==0))/length(which(pe==0)),
    length(which(test_outcome==0 & classes==-1))/length(which(classes==-1)),
    length(which(test_outcome==0 & classes==0))/length(which(classes==0)))->low_risk_acc
  table(classes)
  ####
  fit <- survfit(Surv(time,event)~(classes),data=test)
  print(surv_pvalue(fit))
  ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = TRUE,         # show confidence intervals for
    # point estimaes of survival curves.
    #conf.int.style = "step",  # customize style of confidence intervals
    xlab = "Time in years",   # customize X axis label
    #   xlim=c(0,2),
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    risk.table = "abs_pct",  # absolute number and percentage at risk.
    risk.table.y.text.col = T,# colour risk table text annotations.
    risk.table.y.text = FALSE,# show bars instead of names in text annotations
    # in legend of risk table.
    ncensor.plot = TRUE,      # plot the number of censored subjects at time t
    surv.median.line = "hv",  # add the median survival pointer.
    # legend.labs =
    #   c("low cert","low uncert","mid uncert","high uncert","high cert"
    #   ),    # change legend labels.
    # palette =
    #   c("red","pink","green",
    #     "#E7B800", "#2E9FDF") # custom color palettes.
  )
  fit.p <- survfit(Surv(time,event)~pe,data=test)
  print(surv_pvalue(fit.p))
  ggsurvplot(
    fit.p,                     # survfit object with calculated statistics.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = TRUE,         # show confidence intervals for
    # point estimaes of survival curves.
    #conf.int.style = "step",  # customize style of confidence intervals
    xlab = "Time in years",   # customize X axis label
    #   xlim=c(0,2),
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    risk.table = "abs_pct",  # absolute number and percentage at risk.
    risk.table.y.text.col = T,# colour risk table text annotations.
    risk.table.y.text = FALSE,# show bars instead of names in text annotations
    # in legend of risk table.
    ncensor.plot = TRUE,      # plot the number of censored subjects at time t
    surv.median.line = "hv",  # add the median survival pointer.
  )
  ############something wrong, try with just q
  ####################################
  #######################################
########################################
  ###################################
  #all->all_hist
  all$time[which(all$time>10)]<-10
  all$event[which(all_hist$time>10& all_hist$event==1)]<-0
  hist(all$stage)
  (all$stage-min(all$stage,na.rm=T))/(max(all$stage,na.rm=T)-1)->stage_n
  table(stage_n)
  stage_n->all$pos_score##high is bad
  hist((all$age-min(all$age,na.rm = T))/(max(all$age,na.rm=T)-min(all$age,na.rm=T)))
  (all$age-min(all$age,na.rm = T))/(max(all$age,na.rm=T)-min(all$age,na.rm=T))->age_n
  1-age_n->all$neg_score##so high is good
  ##################################
  samp<-sample(nrow(all),2*ceiling(nrow(all)/3),replace = F)
  all[samp[1:(2*ceiling(nrow(all)/3)-200)],]->train
  all[samp[-c(1:(2*ceiling(nrow(all)/3-200)))],]->calib
  all[-samp,]->test  
  head(all)
  boxplot(all$pos_score[which(all$event==1)],all$pos_score[which(all$event==0)])##score is lower if no offence, so high score more risk
  boxplot(all$neg_score[which(all$event==1)],all$neg_score[which(all$event==0)])##score is lower if no offence, so high score more risk
  surv_cutpoint(train,time="time",event = "event",variables ="score")->a1
  train[which(train$score>as.numeric(a1$cutpoint)),]->train2
  surv_cutpoint(train2,time="time",event = "event",variables ="score")->a2
  unlist(c(as.numeric(a1$cutpoint[1]),as.numeric(a2$cutpoint[1])))->cps
  pred_train<-rep(2,nrow(train))
  pred_train[which(as.numeric(train$score)<cps[2])]<-1
  pred_train[which(train$score<cps[1])]<-0
  table(pred_train)
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
      #xlim=c(0,2),
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
  calib_outcome<-rep(2,nrow(calib))
  calib_outcome[which(calib$event==0)]<-0
  calib_outcome[which(calib$event==1&calib$time<5)]<-1
  table(calib_outcome)
  conformal_score<-rep(NA,nrow(calib))
  conformal_score[which(calib_outcome==2)]<-abs(cps[2]-calib$score)[which(calib_outcome==2)]
  conformal_score[which(calib_outcome==2 & calib$score>cps[2])]<-0
  ##same for outcome 0 but with earlier thresh
  conformal_score[which(calib_outcome==0)]<-abs(cps[1]-calib$score)[which(calib_outcome==0)]
  conformal_score[which(calib_outcome==0 & calib$score<cps[1])]<-0
  ###for outcome 1 need more to pick the minimum distance from either thresh
  conformal_score[which(calib_outcome==1)]<-apply(cbind(abs(cps[2]-calib$score),abs(cps[1]-calib$score)),1,min)[which(calib_outcome==1)]
  conformal_score[which(calib_outcome==1 & calib$score<cps[2] & calib$score>cps[1])]<-0
  # #################
  hist(conformal_score)
  hist(conformal_score[which(calib_outcome==2)])
  hist(conformal_score[which(calib_outcome==0)])
  hist(conformal_score[which(calib_outcome==1)])##different dist so treat seperately
  # ################its generally accurate but sometimes not and sometimes its very wrong
  alpha<-0.1
  #n_cauc<-nrow(calib)
  q_aa0<-q_cauc0<-ceiling((length(which(calib_outcome==0))+1)*(1-alpha))/length(which(calib_outcome==0))
  q_aa1<-q_cauc1<-ceiling((length(which(calib_outcome==1))+1)*(1-alpha))/length(which(calib_outcome==1))
  q_aa2<-q_cauc2<-ceiling((length(which(calib_outcome==2))+1)*(1-alpha))/length(which(calib_outcome==2))
  s_q_a0<-s_q_c0<-quantile(conformal_score[which(calib_outcome==0)],q_cauc0)
  s_q_a1<-s_q_c1<-quantile(conformal_score[which(calib_outcome==1)],q_cauc1)
  s_q_a2<-s_q_c2<-quantile(conformal_score[which(calib_outcome==2)],q_cauc2)
  ######ok now compute the conformal score globally for all test set in 3 settings
  conformal_score_test_0<-mapply(function(i){return(if(test$score[i]<cps[1]){0}else{abs(cps[1]-test$score[i])})},1:nrow(test))
  conformal_score_test_1<-mapply(function(i){return(if(test$score[i]<cps[2] & test$score[i]>cps[1]){0}else{min(abs(cps[2]-test$score[i]),abs(cps[1]-test$score[i]))})},1:nrow(test))
  conformal_score_test_2<-mapply(function(i){return(if(test$score[i]>cps[2]){0}else{abs(cps[2]-test$score[i])})},1:nrow(test))
  ####################asign classes
  preds<-mat.or.vec(nrow(test),3)
  preds[which(conformal_score_test_0<s_q_c0),1]<-1
  preds[which(conformal_score_test_1<s_q_c1),2]<-1
  preds[which(conformal_score_test_2<s_q_c2),3]<-1
  table(apply(preds,1,sum))###nothing is cerain! almost all of them have all 3 and some have 2
  table(preds[,1])
  table(preds[,2])
  table(preds[,3])###ok so every conformal score could be the middle one, because of the alpha choice set to 0.2 for now
  test_outcome<-rep(2,nrow(test))
  test_outcome[which(test$event==0)]<-0
  test_outcome[which(test$event==1&test$time<5)]<-1
  table(test_outcome)
  ##
  ###################
  length(which(test_outcome==0 & preds[,1]==1))/length(which(test_outcome==0))
  length(which(test_outcome==1 & preds[,2]==1))/length(which(test_outcome==1))
  length(which(test_outcome==2 & preds[,3]==1))/length(which(test_outcome==2))###more wrong for early offenders
  ######################
  coverage<-(length(which(test_outcome==0 & preds[,1]==1))+
               length(which(test_outcome==1 & preds[,2]==1))+
               length(which(test_outcome==2 & preds[,3]==1)))/length(test_outcome)##correct coverage
  ######also want overall accuracy
  summary(lm(apply(preds,1,sum)~as.factor(test_outcome)))
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
  table(classes)
  ###now survival analysis
  fit <- survfit(Surv(time,event)~(classes),data=test)
  print(surv_pvalue(fit))
  ggsurvplot(
   fit,                     # survfit object with calculated statistics.
     pval = TRUE,             # show p-value of log-rank test.
     conf.int = TRUE,         # show confidence intervals for
     # point estimaes of survival curves.
     #conf.int.style = "step",  # customize style of confidence intervals
     xlab = "Time in years",   # customize X axis label
  #   xlim=c(0,2),
     ggtheme = theme_light(), # customize plot and risk table with a theme.
     risk.table = "abs_pct",  # absolute number and percentage at risk.
     risk.table.y.text.col = T,# colour risk table text annotations.
     risk.table.y.text = FALSE,# show bars instead of names in text annotations
     # in legend of risk table.
     ncensor.plot = TRUE,      # plot the number of censored subjects at time t
     surv.median.line = "hv",  # add the median survival pointer.
     # legend.labs =
     #   c("low cert","low uncert","mid uncert","high uncert","high cert"
     #   ),    # change legend labels.
     # palette =
     #   c("red","pink","green",
     #     "#E7B800", "#2E9FDF") # custom color palettes.
   )
  fit.p <- survfit(Surv(time,event)~pe,data=test)
  print(surv_pvalue(fit.p))
  ggsurvplot(
     fit.p,                     # survfit object with calculated statistics.
     pval = TRUE,             # show p-value of log-rank test.
     conf.int = TRUE,         # show confidence intervals for
     # point estimaes of survival curves.
     #conf.int.style = "step",  # customize style of confidence intervals
     xlab = "Time in years",   # customize X axis label
  #   xlim=c(0,2),
      ggtheme = theme_light(), # customize plot and risk table with a theme.
   risk.table = "abs_pct",  # absolute number and percentage at risk.
     risk.table.y.text.col = T,# colour risk table text annotations.
     risk.table.y.text = FALSE,# show bars instead of names in text annotations
     # in legend of risk table.
    ncensor.plot = TRUE,      # plot the number of censored subjects at time t
     surv.median.line = "hv",  # add the median survival pointer.
  )
  
  c(cps,s_q_c0,s_q_c1,s_q_c2,coverage,point_accuracy,un_race,un_out,acc_break_down)->out
  
  out->train_calib_test[i,]}

#################

colnames(train_calib_test)<-c("lower_cut_off","upper_cut_off","q0","q1","q2","coverage","accuracy","white_uc","black_uc","0_uc","1_uc","2_uc",
                              
                              "low_risk_acc","low_risk_acc_c","low_risk_acc_uc","low_risk_black","low_risk_black_c","low_risk_black_uc","low_risk_white","low_risk_white_c","low_risk_white_uc",
                              
                              "high_risk_acc","high_risk_acc_c","high_risk_acc_uc","high_risk_black","high_risk_black_c","high_risk_black_uc","high_risk_white","high_risk_white_c","high_risk_white_uc")

# colnames(train_calib_test)<-c("cut_off","coverage","white_uc","black_uc","1_uc","2_uc",
