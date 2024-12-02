
setwd("D:/Documents/Research/TMR/closedLoop/bimanual_closed_loop/data/group/")
library(ez)  
library(ggplot2)  
library(Rmisc)
library(matlab)
library(ungeviz)
library(lme4)
load(file = 'MSLSummary.RData') #  from Behav_11_SRTTAnalyses

library(tidyverse)
library(ggsci)
library(see)
library(cowplot)
.


up1=rgb(252,195,255, maxColorValue = 255)
up2=rgb(172,40,174, maxColorValue = 255)
up3=rgb(83,1,83, maxColorValue = 255)
down1=rgb(197,197,255, maxColorValue = 255)
down2=rgb(92,92,255, maxColorValue = 255)
down3=rgb(1,1,131, maxColorValue = 255)
not1=rgb(176,254,254, maxColorValue = 255)
not2=rgb(0,104,103, maxColorValue = 255)
not3=rgb(0,41,41, maxColorValue = 255)

duration = read.csv(file = 'sleepDuration.csv',header = T,sep = ';');#  from EEG_0_ExtractSleepscore
duration$Sub= as.factor(duration$Sub)
duration$Stage= as.factor(duration$Stage)
tp_fp = read.csv(file = 'trueFalsePositive.csv',header = T,sep = ';');
tp_fp$type = as.factor(tp_fp$type)


offLineGain=matrix(,length(allSub)*length(allSequence),15)
colnames(offLineGain)=c(colnames(MSLSummary)[c(1,4:7)],'Gain_RT','Gain_Acc','N1','N2','N3','REM','WAKE','TP','FP','Tot_stim')
offLineGain=as.data.frame(offLineGain)
preBlocks = c(22:24)
postBlocks = c(25:27)
# preBlocks = c(1:21)
# postBlocks = c(25:45)

counter = 1
for (idx_sub in 1:length(allSub))
{
  for (idx_seq in 1:length(allCond))
  {
    
    meanPreRT   = mean(MSLSummary$Median[MSLSummary$Sub==allSub[idx_sub] &
                                           MSLSummary$Block %in% as.factor(preBlocks) &
                                           MSLSummary$Condition==allCond[idx_seq]],na.rm = T)
    
    meanPostRT  = mean(MSLSummary$Median[MSLSummary$Sub==allSub[idx_sub] &
                                           MSLSummary$Block %in% as.factor(postBlocks) &
                                           MSLSummary$Condition==allCond[idx_seq]],na.rm = T)
    
    meanPreAcc  = mean(MSLSummary$Acc[MSLSummary$Sub==allSub[idx_sub] &
                                         MSLSummary$Block %in% as.factor(preBlocks) &
                                         MSLSummary$Condition==allCond[idx_seq]],na.rm = T)

    meanPostAcc = mean(MSLSummary$Acc[MSLSummary$Sub==allSub[idx_sub] &
                                         MSLSummary$Block %in% as.factor(postBlocks) &
                                         MSLSummary$Condition==allCond[idx_seq]],na.rm = T)
   
    offLineGain$Sub[counter]       = allSub[idx_sub]
    offLineGain$Sequence[counter]  = as.character(MSLSummary$Sequence[MSLSummary$Sub==allSub[idx_sub] &
                                                                        MSLSummary$Block %in% as.factor(preBlocks) &
                                                                        MSLSummary$Condition==allCond[idx_seq]][1])
    
    offLineGain$Condition[counter] = allCond[idx_seq]
    offLineGain$Order[counter] = as.character(MSLSummary$Order[MSLSummary$Sub==allSub[idx_sub] &
                                                                    MSLSummary$Block %in% as.factor(preBlocks) &
                                                                    MSLSummary$Condition==allCond[idx_seq]][1])
    
    if (MSLSummary$Sound[MSLSummary$Sub==allSub[idx_sub] &
                         MSLSummary$Block %in% as.factor(preBlocks) &
                         MSLSummary$Condition==allCond[idx_seq]][1]=="1")
    {offLineGain$Sound[counter] ='Low pitch'}
    
    if (MSLSummary$Sound[MSLSummary$Sub==allSub[idx_sub] &
                         MSLSummary$Block %in% as.factor(preBlocks) &
                         MSLSummary$Condition==allCond[idx_seq]][1]=="2")
    {offLineGain$Sound[counter] ='White noise'}
    
    if (MSLSummary$Sound[MSLSummary$Sub==allSub[idx_sub] &
                         MSLSummary$Block %in% as.factor(preBlocks) &
                         MSLSummary$Condition==allCond[idx_seq]][1]=="3")
    {offLineGain$Sound[counter] ='High pitch'}
    
    
    offLineGain$Sound[counter] = as.character(MSLSummary$Sound[MSLSummary$Sub==allSub[idx_sub] &
                                                                    MSLSummary$Block %in% as.factor(preBlocks) &
                                                                    MSLSummary$Condition==allCond[idx_seq]][1])

    
    offLineGain$Gain_RT[counter]   = ((meanPreRT-meanPostRT)/meanPreRT)*100
    offLineGain$Gain_Acc[counter]  = ((meanPostAcc-meanPreAcc)/meanPreAcc)*100

    if (allSub[idx_sub]=='CL_28')
    {
      offLineGain$N1[counter]   = 'NA'
      offLineGain$N2[counter]   = 'NA'
      offLineGain$N3[counter]   = 'NA'
      offLineGain$REM[counter]  = 'NA'
      offLineGain$WAKE[counter] = 'NA'
    }
    else
    {
      offLineGain$N1[counter]   = duration$DV[duration$Stage==' S1' & duration$Sub ==allSub[idx_sub]]
      offLineGain$N2[counter]   = duration$DV[duration$Stage==' S2' & duration$Sub ==allSub[idx_sub]]
      offLineGain$N3[counter]   = duration$DV[duration$Stage==' S3' & duration$Sub ==allSub[idx_sub]]
      offLineGain$REM[counter]  = duration$DV[duration$Stage==' REM' & duration$Sub ==allSub[idx_sub]]
      offLineGain$WAKE[counter] = duration$DV[duration$Stage==' wake' & duration$Sub ==allSub[idx_sub]]

    }
    counter = counter+1
  }
}


offLineGain$Sub       = as.factor(offLineGain$Sub)
offLineGain$Sequence  = as.factor(offLineGain$Sequence)
offLineGain$Condition = as.factor(offLineGain$Condition)
offLineGain$Order     = as.factor(offLineGain$Order)
offLineGain$Sound     = as.factor(offLineGain$Sound)
offLineGain$TP     = as.numeric(offLineGain$TP)
offLineGain$FP     = as.numeric(offLineGain$FP)
offLineGain$Tot_stim     = as.numeric(offLineGain$Tot_stim)


sorted_Condition  =  c( ' up' , ' down',' not' )
offLineGain$Condition = factor(offLineGain$Condition, levels = sorted_Condition)


#RT
ezANOVA(  offLineGain, 
        dv = .(   Gain_RT),
        wid=.(Sub), 
        within = .(Condition), 
        detailed=T)


upVsDown  = t.test(offLineGain$Gain_RT[offLineGain$Condition==' up' ],offLineGain$Gain_RT[offLineGain$Condition==' down' ] ,paired = T,alternative = 'greater')
upVsNot   = t.test(offLineGain$Gain_RT[offLineGain$Condition==' up' ],offLineGain$Gain_RT[offLineGain$Condition==' not' ],paired = T,alternative = 'greater')
downVsNot = t.test(offLineGain$Gain_RT[offLineGain$Condition==' down' ],offLineGain$Gain_RT[offLineGain$Condition==' not' ],paired = T,alternative = 'less')

cohensD(offLineGain$Gain_RT[offLineGain$Condition==' up'],offLineGain$Gain_RT[offLineGain$Condition==' down'],method = 'paired')
cohensD(offLineGain$Gain_RT[offLineGain$Condition==' up'],offLineGain$Gain_RT[offLineGain$Condition==' not'],method = 'paired')
cohensD(offLineGain$Gain_RT[offLineGain$Condition==' down'],offLineGain$Gain_RT[offLineGain$Condition==' not'],method = 'paired')

p.adjust(c(upVsDown$p.value,upVsNot$p.value,downVsNot$p.value),n= 3, method = 'fdr')

#Acc Fig 2b
ggplot(summarySE(offLineGain, measurevar="Gain_RT", 
                 groupvars=c("Sub","Condition"),na.rm=T), aes( y=Gain_RT, fill = Condition,x = Condition)) +
  geom_violin(position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c(up3,down3,not3))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = "grey")+
  stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Condition), position=position_dodge(1),
               color=c( up1,down1,not1) ) +
  stat_summary(fun=median, geom = "crossbar", width = 0.7,size=0.5,alpha = 0.5, aes(group=Condition), position=position_dodge(1),
               color=c( up1,down1,not1)) +
  geom_point(aes(fill=Condition,group=Sub),size=2,shape=19, position = position_dodge(0.2),
             color = rep(c(up1,down1,not1), times = length(allSub))) +
  # facet_wrap(~ROI_Name)+
  # coord_cartesian(ylim=c(-10,100))+
  theme_classic() 



#Acc Fig S9 b
ezANOVA(offLineGain , 
        dv = .(   Gain_Acc),
        wid=.(Sub), 
        within = .(Condition), 
        detailed=T)


ggplot(summarySE(offLineGain, measurevar="Gain_Acc", 
                 groupvars=c("Sub","Condition"),na.rm=T), aes( y=Gain_Acc, fill = Condition,x = Condition)) +
  geom_violin(position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c(up3,down3,not3))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = "grey")+
  stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Condition), position=position_dodge(1),
               color=c( up1,down1,not1) ) +
  stat_summary(fun=median, geom = "crossbar", width = 0.6,size=0.5,alpha = 0.5, aes(group=Condition), position=position_dodge(1),
               color=c( up1,down1,not1)) +
  geom_point(aes(fill=Condition,group=Sub),size=2,shape=19, position = position_dodge(0.2),
             color = rep(c(up1,down1,not1), times = length(allSub))) +
  # facet_wrap(~ROI_Name)+
  coord_cartesian(ylim=c(-25,50))+
  theme_classic() 

 




