
setwd("D:/Documents/Research/TMR/closedLoop/bimanual_closed_loop/data/group/")
library(ez)  
library(ggplot2)  
library(Rmisc)
library(matlab)

CI_perso<- function (x, ci = 0.95) 
{
  a <- mean(x)
  s <- sd(x)
  n <- length(x)
  error <- qt(ci + (1 - ci)/2, df = n - 1) * s/sqrt(n)
  return(round(c(mean = a, lower = a - error,upper = a + error),1)  )
}



### Impact of the number of cues on offline changes in performance
load(file = 'offLineGain.RData') #  from BEHAV_12_MSL_rofflineGains

TMRIndex=matrix(,length(allSub)*2,3)
colnames(TMRIndex)=c(colnames(MSLSummary)[c(1)],'Condition','index_RT','N1','N2','N3','REM','WAKE')
TMRIndex=as.data.frame(TMRIndex)
counter = 1
for (idx_sub in 1:length(allSub))
{
  
  
  offlineGainUp  = offLineGain$Gain_RT[offLineGain$Sub==allSub[idx_sub] & offLineGain$Condition==' up' ]
  offlineGaiDown = offLineGain$Gain_RT[offLineGain$Sub==allSub[idx_sub] & offLineGain$Condition==' down' ]
  offlineGainNot = offLineGain$Gain_RT[offLineGain$Sub==allSub[idx_sub] & offLineGain$Condition==' not' ]
  
  TMRIndex$Sub[counter]       = allSub[idx_sub]
  TMRIndex$Condition[counter]      = 'up'
  TMRIndex$index_RT[counter]  = (offlineGainUp-offlineGainNot)
  if (allSub[idx_sub]=='CL_28')
  {
    TMRIndex$N1[counter]   = 'NA'
    TMRIndex$N2[counter]   = 'NA'
    TMRIndex$N3[counter]   = 'NA'
    TMRIndex$REM[counter]  = 'NA'
    TMRIndex$WAKE[counter] = 'NA'
  }
  else
  {
    TMRIndex$N1[counter]   = duration$DV[duration$Stage==' S1' & duration$Sub ==allSub[idx_sub]]
    TMRIndex$N2[counter]   = duration$DV[duration$Stage==' S2' & duration$Sub ==allSub[idx_sub]]
    TMRIndex$N3[counter]   = duration$DV[duration$Stage==' S3' & duration$Sub ==allSub[idx_sub]]
    TMRIndex$REM[counter]  = duration$DV[duration$Stage==' REM' & duration$Sub ==allSub[idx_sub]]
    TMRIndex$WAKE[counter] = duration$DV[duration$Stage==' wake' & duration$Sub ==allSub[idx_sub]]
    
  }
  
  counter = counter+1
  
  TMRIndex$Sub[counter]       = allSub[idx_sub]
  TMRIndex$Condition[counter]      = 'down'
  TMRIndex$index_RT[counter]  = (offlineGaiDown-offlineGainNot)
  
  if (allSub[idx_sub]=='CL_28')
  {
    TMRIndex$N1[counter]   = 'NA'
    TMRIndex$N2[counter]   = 'NA'
    TMRIndex$N3[counter]   = 'NA'
    TMRIndex$REM[counter]  = 'NA'
    TMRIndex$WAKE[counter] = 'NA'
  }
  else
  {
    TMRIndex$N1[counter]   = duration$DV[duration$Stage==' S1' & duration$Sub ==allSub[idx_sub]]
    TMRIndex$N2[counter]   = duration$DV[duration$Stage==' S2' & duration$Sub ==allSub[idx_sub]]
    TMRIndex$N3[counter]   = duration$DV[duration$Stage==' S3' & duration$Sub ==allSub[idx_sub]]
    TMRIndex$REM[counter]  = duration$DV[duration$Stage==' REM' & duration$Sub ==allSub[idx_sub]]
    TMRIndex$WAKE[counter] = duration$DV[duration$Stage==' wake' & duration$Sub ==allSub[idx_sub]]
    
  }
  counter = counter+1
  
}


TMRIndex$N1 = as.numeric(TMRIndex$N1)
TMRIndex$N2 = as.numeric(TMRIndex$N2)
TMRIndex$N3 = as.numeric(TMRIndex$N3)
TMRIndex$REM = as.numeric(TMRIndex$REM)
TMRIndex$WAKE = as.numeric(TMRIndex$WAKE)

cor.test(TMRIndex$index_RT[TMRIndex$Condition=='up'],TMRIndex$N1[TMRIndex$Condition=='up'])
cor.test(TMRIndex$index_RT[TMRIndex$Condition=='down'],TMRIndex$N1[TMRIndex$Condition=='down'])
ggplot(TMRIndex , aes(x = index_RT, y = N1, color = Condition )) +
  geom_point(size =2,alpha = 0.5 ) +
  scale_color_manual(values = c(down2,up2,down2,up2)) +
  geom_smooth(method = lm,alpha=0.2,size = 1) +
  theme_minimal() +
  ggtitle("")+
  ylab('Time spent in N1') +
  xlab('TMR Index') +
  theme(axis.line = element_line(color = "grey70"))


cor.test(TMRIndex$index_RT,TMRIndex$N2)
cor.test(TMRIndex$index_RT[TMRIndex$Condition=='up'],TMRIndex$N2[TMRIndex$Condition=='up'])
cor.test(TMRIndex$index_RT[TMRIndex$Condition=='down'],TMRIndex$N2[TMRIndex$Condition=='down'])
ggplot(TMRIndex , aes(x = index_RT, y = N2, color = Condition )) +
  geom_point(size =2,alpha = 0.5 ) +
  scale_color_manual(values = c(down2,up2,down2,up2)) +
  geom_smooth(method = lm,alpha=0.2,size = 1) +
  theme_minimal() +
  ggtitle("")+
  ylab('Time spent in N2') +
  xlab('TMR Index') +
  theme(axis.line = element_line(color = "grey70"))


cor.test(TMRIndex$index_RT,TMRIndex$N3)
cor.test(TMRIndex$index_RT[TMRIndex$Condition=='up'],TMRIndex$N3[TMRIndex$Condition=='up'])
cor.test(TMRIndex$index_RT[TMRIndex$Condition=='down'],TMRIndex$N3[TMRIndex$Condition=='down'])
ggplot(TMRIndex , aes(x = index_RT, y = N3, color = Condition )) +
  geom_point(size =2,alpha = 0.5 ) +
  scale_color_manual(values = c(down2,up2,down2,up2)) +
  geom_smooth(method = lm,alpha=0.2,size = 1) +
  theme_minimal() +
  ggtitle("")+
  ylab('Time spent in N3') +
  xlab('TMR Index') +
  theme(axis.line = element_line(color = "grey70"))


cor.test(TMRIndex$index_RT,TMRIndex$REM)
cor.test(TMRIndex$index_RT[TMRIndex$Condition=='up'],TMRIndex$REM[TMRIndex$Condition=='up'])
cor.test(TMRIndex$index_RT[TMRIndex$Condition=='down'],TMRIndex$REM[TMRIndex$Condition=='down'])
ggplot(TMRIndex , aes(x = index_RT, y = REM, color = Condition )) +
  geom_point(size =2,alpha = 0.5 ) +
  scale_color_manual(values = c(down2,up2,down2,up2)) +
  geom_smooth(method = lm,alpha=0.2,size = 1) +
  theme_minimal() +
  ggtitle("")+
  ylab('Time spent in REM') +
  xlab('TMR Index') +
  theme(axis.line = element_line(color = "grey70"))




cor.test(TMRIndex$index_RT,TMRIndex$WAKE)
cor.test(TMRIndex$index_RT[TMRIndex$Condition=='up'],TMRIndex$WAKE[TMRIndex$Condition=='up'])
cor.test(TMRIndex$index_RT[TMRIndex$Condition=='down'],TMRIndex$WAKE[TMRIndex$Condition=='down'])
ggplot(TMRIndex , aes(x = index_RT, y = WAKE, color = Condition )) +
  geom_point(size =2,alpha = 0.5 ) +
  scale_color_manual(values = c(down2,up2,down2,up2)) +
  geom_smooth(method = lm,alpha=0.2,size = 1) +
  theme_minimal() +
  ggtitle("")+
  ylab('Time spent in WAKE') +
  xlab('TMR Index') +
  theme(axis.line = element_line(color = "grey70"))



# Correlation gain / stim
ggplot(offLineGain[offLineGain$Condition!=" not",] , aes(x = Gain_RT, y = TP, color = Condition )) +
  geom_point(size =2,alpha = 0.5 ) +
  scale_color_manual(values = c(down2,up2,down2,up2)) +
  geom_smooth(method = lm,alpha=0.2,size = 1) +
  theme_minimal() +
  ggtitle("")+
  ylab('Time spent in N1') +
  xlab('TMR Index') +
  theme(axis.line = element_line(color = "grey70"))

cor.test(offLineGain$Gain_RT[offLineGain$Condition==' up'],offLineGain$TP[offLineGain$Condition==' up'])
cor.test(offLineGain$Gain_RT[offLineGain$Condition==' down'],offLineGain$TP[offLineGain$Condition==' down'])



ggplot(offLineGain[offLineGain$Condition!=" not",] , aes(x = Gain_RT, y = FP, color = Condition )) +
  geom_point(size =2,alpha = 0.5 ) +
  scale_color_manual(values = c(down2,up2,down2,up2)) +
  geom_smooth(method = lm,alpha=0.2,size = 1) +
  theme_minimal() +
  ggtitle("")+
  ylab('Time spent in N1') +
  xlab('TMR Index') +
  theme(axis.line = element_line(color = "grey70"))

cor.test(offLineGain$Gain_RT[offLineGain$Condition==' up'],offLineGain$FP[offLineGain$Condition==' up'])
cor.test(offLineGain$Gain_RT[offLineGain$Condition==' down'],offLineGain$FP[offLineGain$Condition==' down'])



ggplot(offLineGain[offLineGain$Condition!=" not",] , aes(x = Gain_RT, y = Tot_stim, color = Condition )) +
  geom_point(size =2,alpha = 0.5 ) +
  scale_color_manual(values = c(down2,up2,down2,up2)) +
  geom_smooth(method = lm,alpha=0.2,size = 1) +
  theme_minimal() +
  ggtitle("")+
  ylab('Time spent in N1') +
  xlab('TMR Index') +
  theme(axis.line = element_line(color = "grey70"))

cor.test(offLineGain$Gain_RT[offLineGain$Condition==' up'],offLineGain$Tot_stim[offLineGain$Condition==' up'])
cor.test(offLineGain$Gain_RT[offLineGain$Condition==' down'],offLineGain$Tot_stim[offLineGain$Condition==' down'])

  
#Equivalent number of discarded trials in the different sessions and conditions

# MSLData from Behav_11_SRTTAnalyses (after outlyier rejection)


DiscardedTrl =matrix(,length(allSub)*(length(allSess)+1)*length(allSequence),5)
colnames(DiscardedTrl)=c('Sub','Condition','Session',"PercDisc",'Nb')
DiscardedTrl=as.data.frame(DiscardedTrl)
counter = 1
for (idx_sub in 1:length(allSub))
{
  for (idx_seq in 1:length(allSequence))
  {
    
    for (idx_sess in 1 : length(allSess))
    {
      
      
      tmp = MSLData[MSLData$Sub==allSub[idx_sub] & 
                      MSLData$Session==allSess[idx_sess] &
                      MSLData$Condition==allCond[idx_seq]  
                    , ]
      DiscardedTrl$Sub[counter]       = allSub[idx_sub]
      DiscardedTrl$Session[counter]   = as.character(tmp$Session[1])
      DiscardedTrl$Condition[counter] = as.character(tmp$Condition[1])
      DiscardedTrl$PercDisc[counter]  = (length(tmp$RT[ tmp$Outlier=='FALSE' ] )/length(tmp$RT))*100
      DiscardedTrl$Nb[counter]  = length(tmp$RT[tmp$Outlier=='FALSE' ]  )
      counter = counter+1
      
      if (idx_sess ==length(allSess))
      {
        tmp = MSLData[MSLData$Sub==allSub[idx_sub] & 
                        MSLData$Session==allSess[idx_sess] &
                        MSLData$Condition==allCond[idx_seq] 
                      & MSLData$Block %in% as.factor(c(25:27)), ]
        DiscardedTrl$Sub[counter]       = allSub[idx_sub]
        DiscardedTrl$Session[counter]   = 'pre3'
        DiscardedTrl$Condition[counter] = as.character(tmp$Condition[1])
        DiscardedTrl$PercDisc[counter]  = (length(tmp$RT[ tmp$Outlier=='FALSE' ] )/length(tmp$RT))*100
        DiscardedTrl$Nb[counter]  = length(tmp$RT[tmp$Outlier=='FALSE']  )
        counter = counter+1
        
      }
    }
  }
}

DiscardedTrl$Session = factor(DiscardedTrl$Session, levels = c('1',  '2',  'pre3','3'))

tmp = tapply(DiscardedTrl$PercDisc,list(DiscardedTrl$Condition,DiscardedTrl$Session),CI_perso)

ezANOVA(DiscardedTrl[!DiscardedTrl$Session %in% c('2','pre3'),], 
        dv = .(PercDisc),
        wid=.(Sub), 
        within = .(Session,Condition), 
        detailed=T)


ezANOVA(DiscardedTrl[DiscardedTrl$Session %in% c('2','pre3'),], 
        dv = .(PercDisc),
        wid=.(Sub), 
        within = .(Session,Condition), 
        detailed=T)


plt = ezPlot(
  data = outlierTrl#[outlierTrl$Session %in% c('2','pre3'),]
  , dv = .(PercDisc)
  , wid = .(Sub)
  , within = .(Session)
  , x = .(Session)
  , x_lab = 'Session'
  , y_lab = 'PercDisc'
)
print(plt)

a = t.test(outlierTrl$PercDisc[outlierTrl$Session=='1'],outlierTrl$PercDisc[outlierTrl$Session=='2'],paired = T)
b = t.test(outlierTrl$PercDisc[outlierTrl$Session=='1'],outlierTrl$PercDisc[outlierTrl$Session=='3'],paired = T)
c = t.test(outlierTrl$PercDisc[outlierTrl$Session=='2'],outlierTrl$PercDisc[outlierTrl$Session=='3'],paired = T)

t.test(outlierTrl$PercDisc[outlierTrl$Session=='2'],outlierTrl$PercDisc[outlierTrl$Session=='pre3'],paired = T)



cohensD(outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' up'],
        outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' down'] )

b = t.test(outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' up'],
           outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' not']  ,paired = T)
cohensD(outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' up'],
        outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' not']  )

c = t.test(outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' down'],
           outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' not']  ,paired = T)
cohensD(outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' down'],
        outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' not']  )





ezANOVA(outlierTrl[outlierTrl$Session %in% c('2','pre3'),], 
        dv = .(PercDisc),
        wid=.(Sub), 
        within = .(Session,Condition), 
        detailed=T)



a = t.test(outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' up'],
           outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' down']  ,paired = T)
cohensD(outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' up'],
        outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' down'] )

b = t.test(outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' up'],
           outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' not']  ,paired = T)
cohensD(outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' up'],
        outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' not']  )

c = t.test(outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' down'],
           outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' not']  ,paired = T)
cohensD(outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' down'],
        outlierTrl$PercDisc[outlierTrl$Session %in% c('2','pre3') & outlierTrl$Condition == ' not']  )

p.adjust(c(a$p.value,b$p.value,c$p.value),method = 'fdr')

plt = ezPlot(
  data = outlierTrl#[outlierTrl$Session %in% c('2','pre3'),]
  , dv = .(PercDisc)
  , wid = .(Sub)
  , within = .(Condition,Session)
  , x = .(Session)
  , split = .(Condition)
  , x_lab = 'Session'
  , y_lab = 'PercDisc'
  , split_lab = 'Condition'
)
print(plt)

  
#Equivalent vigilance during each behavioral session
  ### PVT
PVT = read.csv(file = 'PVT.csv',header = T,sep = ';');

PVT$RT = as.numeric(PVT$RT)*1000
PVT$RT = PVT$RT*1000
PVT$Session= as.factor(PVT$Session)
PVT$Sub= as.factor(PVT$Sub)
CI_perso(PVT$RT[!PVT$Sub %in% c('CL_02','CL_03','CL_04','CL_06','CL_13','CL_17') & PVT$Session==' preNight'])
CI_perso(PVT$RT[!PVT$Sub %in% c('CL_02','CL_03','CL_04','CL_06','CL_13','CL_17') & PVT$Session==' postNight'])

aovPVT = ezANOVA (PVT[!PVT$Sub %in% c('CL_02','CL_03','CL_04','CL_06','CL_13','CL_17'),], dv = .(RT), wid=.(Sub), within = .(Session), detailed=T)


tmpPVT = summarySE(PVT, measurevar="RT", 
                   groupvars=c("Sub","Session"))

tmp =summarySE(PVT,measurevar = 'RT', groupvars=c("Sub","Session"),na.rm=T)

allSub = levels(PVT$Sub)
pvtTsv=matrix(,length(allSub),3)
colnames(pvtTsv)=c('Sub','prenight',"postnight")
pvtTsv=as.data.frame(pvtTsv)
counter = 1
for (idx_sub in 1:length(allSub))
{
  
  
  pvtTsv$Sub[idx_sub]       = allSub[idx_sub]
  pvtTsv$prenight[idx_sub]  = median(PVT$RT[PVT$Sub==allSub[idx_sub] & PVT$Session==' preNight'])
  pvtTsv$postnight[idx_sub] = median(PVT$RT[PVT$Sub==allSub[idx_sub] & PVT$Session==' postNight'])
  
  
}
pvtTsv$Session= as.factor(pvtTsv$Session)
pvtTsv$Sub= as.factor(pvtTsv$Sub)

write.csv(pvtTsv,"PVT.tsv")

  ### SSS
questionnaires = read.csv(file = 'questionnairesSummarize.csv',header = T,sep = '\t');
questionnaires$Qr=as.factor(questionnaires$Qr)

allSub = levels(questionnaires$Qr)
SSS=matrix(,length(allSub)*2 ,3)
colnames(SSS)=c('Sub','Session',"DV")
SSS=as.data.frame(SSS)
counter = 1
for (idx_sub in 1:length(allSub))
{
  for (idx_sess in c(1,2) )
  {
    
    if (idx_sess==1)
    {
      col = 8
    }
    if (idx_sess==2)
    {
      col = 9
    }
    SSS$Sub[counter]     = allSub[idx_sub]
    SSS$Session[counter] = allSess[idx_sess]
    SSS$DV[counter]      = questionnaires[idx_sub,col]
    
    counter = counter +1
  }
}


SSS$Session=as.factor(SSS$Session)
SSS$Sub=as.factor(SSS$Sub)

aovSSS = ezANOVA (SSS, dv = .(DV), wid=.(Sub), within = .(Session), detailed=T)

tapply(SSS$DV,SSS$Session,CI)

### Equivalent baseline performance between the three movement sequences 
#RT

ezANOVA (MSLSummary[MSLSummary$Session =="1",], dv = .(Median), wid = .(Sub),
         within= .(Sequence,Block), detailed=T,type=1)

ezANOVA (MSLSummary[MSLSummary$Session =="2",], dv = .(Median), wid = .(Sub),
         within= .(Sequence,Block), detailed=T)

# Accuracy
ezANOVA (MSLSummary[MSLSummary$Session =="1",], dv = .(Acc), wid = .(Sub),
         within= .(Sequence,Block), detailed=T,type=1)

ezANOVA (MSLSummary[MSLSummary$Session =="2",], dv = .(Acc), wid = .(Sub),
         within= .(Sequence,Block), detailed=T)



### Equivalent baseline performance between the three conditions 
# Median RT
ezANOVA (MSLSummary[MSLSummary$Session =="1",], dv = .(Median), wid = .(Sub),
         within= .(Condition,Block), detailed=T,type=1)

ezANOVA (MSLSummary[MSLSummary$Session =="2" ,], dv = .(Median), wid = .(Sub),
         within= .(Condition,Block), detailed=T)


# Accuracy
ezANOVA (MSLSummary[MSLSummary$Session =="1",], dv = .(Acc), wid = .(Sub),
         within= .(Condition,Block), detailed=T,type=1)

ezANOVA (MSLSummary[MSLSummary$Session =="2",], dv = .(Acc), wid = .(Sub),
         within= .(Condition,Block), detailed=T)



#Equivalent impact of the three different sounds Fig S8
ezANOVA(offLineGain , 
        dv = .(   Gain_RT),
        wid=.(Sub), 
        within = .(Sound), 
        detailed=T)

mod.ez <- ezANOVA(data = datalong,
                  dv = .(rating), 
                  wid = .(id),
                  within = .(factor1, factor2),
                  detailed = TRUE,
                  return_aov = TRUE)

ggplot(summarySE(offLineGain, measurevar="Gain_RT", 
                 groupvars=c("Sub","Sound"),na.rm=T), aes( y=Gain_RT, fill = Sound,x = Sound)) +
  geom_violin(position = position_dodge(width = 0.9))+
  # scale_fill_manual(values = c(up3,down3,not3))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = "grey")+
  stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Sound), position=position_dodge(1) ) +
  stat_summary(fun=median, geom = "crossbar", width = 0.7,size=0.5,alpha = 0.5, aes(group=Sound), position=position_dodge(1)) +
  geom_point(aes(fill=Sound,group=Sub),size=2,shape=19, position = position_dodge(0.2), times = length(allSub)) +
  # facet_wrap(~ROI_Name)+
  # coord_cartesian(ylim=c(-65,50))+
  theme_classic() 


### Motor exectution performances
# randomSummary and #  from Behav_11_SRTTAnalyses  from Behav_11_SRTTAnalyses
learningRate=matrix(,length(allSub)*2,3)
colnames(learningRate)=c('Sub','Task',"Rate")
learningRate=as.data.frame(learningRate)
counter = 1
for (idx_sub in 1:length(allSub))
{
      tmpMSL = MSLSummary[MSLSummary$Sub==allSub[idx_sub] & 
                              MSLSummary$Block %in% as.factor(c(1:4,41:45)),]
      meanMSLpre = mean(tmpMSL$Median[tmpMSL$Session=="1"])
      meanMSLpost = mean(tmpMSL$Median[tmpMSL$Session=="3"])
      
      
      tmpRandom = randomSummary[randomSummary$Sub==allSub[idx_sub] ,]
      meanRandomPre = mean(tmpRandom$Mean[tmpRandom$Condition==" preNight"])
      meanRandomPost = mean(tmpRandom$Mean[tmpRandom$Condition==" postNight"])
      

    learningRate$Sub[counter]   = allSub[idx_sub]
    learningRate$Task[counter]  = 'MSL'
    learningRate$Rate[counter]  =  (meanMSLpre-meanMSLpost)/meanMSLpre*100
    
    learningRate$Sub[counter+1]   = allSub[idx_sub]
    learningRate$Task[counter+1]  = 'Random'
    learningRate$Rate[counter+1]  =  (meanRandomPre-meanRandomPost)/meanRandomPre*100

    
    counter = counter +2
    
  
}


learningRate$Sub       = as.factor(learningRate$Sub)
learningRate$Task  = as.factor(learningRate$Task)
shapiro.test(learningRate$Rate[learningRate$Task=='Random'])
shapiro.test(learningRate$Rate[learningRate$Task=='MSL'])
t.test(learningRate$Rate[learningRate$Task=='Random'],learningRate$Rate[learningRate$Task=='MSL'],paired = T)
CI(learningRate$Rate[learningRate$Task=='Random']);sd(learningRate$Rate[learningRate$Task=='Random'])
CI(learningRate$Rate[learningRate$Task=='MSL']);sd(learningRate$Rate[learningRate$Task=='MSL'])
cor.test(learningRate$Rate[learningRate$Task=='Random'],learningRate$Rate[learningRate$Task=='MSL'],paired = T)

