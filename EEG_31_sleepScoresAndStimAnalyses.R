library('ez')
library('Rmisc')
library('Directional')
library('dplyr')
library('circular')
library('lsr')


setwd("D:/Documents/Research/TMR/closedLoop/bimanual_closed_loop/data/group")

CI_perso<- function (x, ci = 0.95) 
{
  a <- mean(x)
  s <- sd(x)
  n <- length(x)
  error <- qt(ci + (1 - ci)/2, df = n - 1) * s/sqrt(n)
  return(round(c(mean = a, lower = a - error,upper = a + error),1)  )
}

############ SLEEP Duration
duration = read.csv(file = 'sleepDuration.csv',header = T,sep = ';');
duration$Sub = as.factor(duration$Sub)
duration$Stage = as.factor(duration$Stage)
all_sub= levels(duration$Sub)
all_stage = levels(duration$Stage)


WAKE = CI_perso(duration$DV [duration$Stage == ' wake' ])

S1 = CI_perso( duration$DV [duration$Stage == ' S1' ])

S2 = CI_perso(duration$DV [duration$Stage == ' S2' ] )

S3 = CI_perso(duration$DV [duration$Stage == ' S3' ])

REM = CI_perso(duration$DV [duration$Stage == ' REM' ])

NREM = CI_perso(duration$DV [duration$Stage == ' S2' ] + duration$DV [duration$Stage == ' S3' ])  



ICSleepDurationPercentage = CI_perso(duration$Percentage[duration$Stage == ' sleep'])

sleepEff = CI_perso(duration$Percentage[duration$Stage == ' S2' ] + 
     duration$Percentage[ duration$Stage == ' S3' ] + duration$Percentage[ duration$Stage == ' REM' ])

mean(duration$Percentage,na.rm = T)+3*sd(duration$Percentage,na.rm = T)>duration$Percentage
duration$Percentage>mean(duration$Percentage,na.rm = T)-3*sd(duration$Percentage,na.rm = T)

mean(duration$DV [duration$Stage == ' REM' ],na.rm = T)+3*sd(duration$DV [duration$Stage == ' REM' ],na.rm = T)>duration$DV [duration$Stage == ' REM' ]
duration$DV [duration$Stage == ' REM' ]>mean(duration$DV [duration$Stage == ' REM' ],na.rm = T)-3*sd(duration$DV [duration$Stage == ' REM' ],na.rm = T)
mean(duration$DV [duration$Stage == ' S2' ],na.rm = T)+3*sd(duration$DV [duration$Stage == ' S2' ],na.rm = T)>duration$DV [duration$Stage == ' S2' ]
duration$DV [duration$Stage == ' S2' ]>mean(duration$DV [duration$Stage == ' S2' ],na.rm = T)-3*sd(duration$DV [duration$Stage == ' S2' ],na.rm = T)
mean(duration$DV [duration$Stage == ' S1' ],na.rm = T)+3*sd(duration$DV [duration$Stage == ' S1' ],na.rm = T)>duration$DV [duration$Stage == ' S1' ]
duration$DV [duration$Stage == ' S1' ]>mean(duration$DV [duration$Stage == ' S1' ],na.rm = T)-3*sd(duration$DV [duration$Stage == ' S1' ],na.rm = T)
mean(duration$DV [duration$Stage == ' S3' ],na.rm = T)+3*sd(duration$DV [duration$Stage == ' S3' ],na.rm = T)>duration$DV [duration$Stage == ' S3' ]
duration$DV [duration$Stage == ' S3' ]>mean(duration$DV [duration$Stage == ' S3' ],na.rm = T)-3*sd(duration$DV [duration$Stage == ' S3' ],na.rm = T)
mean(duration$DV [duration$Stage == ' wake' ],na.rm = T)+3*sd(duration$DV [duration$Stage == ' wake' ],na.rm = T)>duration$DV [duration$Stage == ' wake' ]
duration$DV [duration$Stage == ' wake' ]>mean(duration$DV [duration$Stage == ' wake' ],na.rm = T)-3*sd(duration$DV [duration$Stage == ' wake' ],na.rm = T)


############ SLEEP Latendy
latency = read.csv(file = 'sleepLatency.csv',header = T,sep = ';');
ICLatencyS1= CI_perso(latency$S1)

mean(latency$S1,na.rm = T)+3*sd(latency$S1,na.rm = T)>latency$S1
latency$S1>mean(latency$S1,na.rm = T)-3*sd(latency$S1,na.rm = T)


############ SLEEP arousal
arousal = read.csv(file = 'sleepArousal.csv',header = T,sep = ';');


mean(arousal$Duration,na.rm = T)+3*sd(arousal$Duration,na.rm = T)>arousal$Duration
arousal$Duration>mean(arousal$Duration,na.rm = T)-3*sd(arousal$Duration,na.rm = T)


############ STIM REPARTITION



detection = read.csv(file = 'logAnalyzer.txt',header = T,sep = '\t')
boxplot(detection$Trough_TP,detection$Peak_TP)

cuesCount = read.csv(file = 'stimEfficiency.csv',header = T,sep = ';');

tapply(cuesCount$Percentage[cuesCount$Stage==' NREM'],cuesCount$type[cuesCount$Stage==' NREM'],CI_perso)

tapply(cuesCount$DV[cuesCount$type==' all'],list(cuesCount$Stage[cuesCount$type==' all']),CI_perso)
tapply(cuesCount$DV[cuesCount$type==' up'],list(cuesCount$Stage[cuesCount$type==' up']),CI_perso)
tapply(cuesCount$DV[cuesCount$type==' down'],list(cuesCount$Stage[cuesCount$type==' down']),CI_perso)

aovCuesCount = ezANOVA (cuesCount[cuesCount$type!=' all' & (cuesCount$Stage!=' total' & cuesCount$Stage!=' NREM'),], dv = .(DV),wid=.(Sub), within = .(type,Stage), detailed=T)



mean(cuesCount$DV[cuesCount$type==' all' & cuesCount$Stage == ' total'],na.rm = T)+3*sd(cuesCount$DV[cuesCount$type==' all' & cuesCount$Stage == ' total'],na.rm = T)>cuesCount$DV[cuesCount$type==' all' & cuesCount$Stage == ' total']
cuesCount$DV[cuesCount$type==' all' & cuesCount$Stage == ' total']>mean(cuesCount$DV[cuesCount$type==' all' & cuesCount$Stage == ' total'],na.rm = T)-3*sd(cuesCount$DV[cuesCount$type==' all' & cuesCount$Stage == ' total'],na.rm = T)

AOVplot = ezPlot(
  data = cuesCount[cuesCount$type!=' all' & (cuesCount$Stage!=' total' & cuesCount$Stage!=' NREM'),],
  , dv = .(DV)
  , wid = .(Sub)
  , within = .(type,Stage)
  , x = .(Stage)
  , split = .(type)
  , x_lab =    'Block'
  , y_lab =    'RT' )
print(AOVplot)

tp_fp = read.csv(file = 'trueFalsePositive.csv',header = T,sep = ';');
tapply(tp_fp$TP ,tp_fp$type,CI_perso)
tapply(tp_fp$percentageTP ,tp_fp$type,CI_perso)

t.test(tp_fp$TP[tp_fp$type==' up'],tp_fp$TP[tp_fp$type==' down'],paired = T)
cor.test(tp_fp$TP[tp_fp$type==' up'],tp_fp$TP[tp_fp$type==' down'],paired = T)

summarySw = read.csv(file = 'summarySW.csv',header = T,sep = ';');
summarySw$SampleTroughEEG = as.numeric(summarySw$SampleTroughEEG)
summarySw$stimEEG = as.numeric(summarySw$stimEEG)
summarySw$stimECHT = as.numeric(summarySw$stimECHT)
summarySw$SampleTroughECHT = as.numeric(summarySw$SampleTroughECHT)
summarySw$SampleTroughECHT2EEG = as.numeric(summarySw$SampleTroughECHT2EEG)

summarySw$StimTroughDelay = (summarySw$stimEEG-summarySw$SampleTroughEEG)*1/500

tmp = summarySE(summarySw[summarySw$isIn==1 & summarySw$isTP==1 
                          & summarySw$SampleTroughECHT2EEG!=0 & summarySw$SampleTroughEEG!=0,], 
                measurevar="StimTroughDelay", 
          groupvars=c("Sub","Cond"))
tmp =tapply(summarySw$isTP[summarySw$isTP==1 | summarySw$isFP==1] ,list(summarySw$Sub[summarySw$isTP==1 | summarySw$isFP==1],summarySw$Cond[summarySw$isTP==1 | summarySw$isFP==1]),length)
tmp2 = tapply(summarySw$isTP ,list(summarySw$Sub,summarySw$Cond),length)
tmp3 = cbind(tmp[,2],tmp[,1],tmp2[,2])


colnames(tmp3) = c('up', 'down', 'not')
percDiscarded = nbTrials/tmp3*100

tapply(tmp$StimTroughDelay*1000, tmp$Cond, CI_perso)

ezPlot (tmp, 
        dv = .(StimTroughDelay), wid=.(Sub), 
        within = .(Cond), x = .(Cond))
          


nbTrialsTroughLocked = read.csv(file = 'nbTrl_anlayzed_onlyTP.txt',header = T,sep = ',');

all_sub= levels(stim_interval_duration$Sub)

df_nbTrialsTroughLocked=matrix(,length(all_sub)*3 ,3)
colnames(df_nbTrialsTroughLocked)=c('Sub','Condition',"nb")
df_nbTrialsTroughLocked=as.data.frame(df_nbTrialsTroughLocked)
counter = 1
for (idx_sub in 1:length(all_sub))
{
  df_nbTrialsTroughLocked$Sub[counter]      = all_sub[idx_sub]
  df_nbTrialsTroughLocked$Condition[counter]     = 'up'
  df_nbTrialsTroughLocked$nb[counter]      = nbTrialsTroughLocked$up[idx_sub]
  df_nbTrialsTroughLocked$Sub[counter+1]      = all_sub[idx_sub]
  df_nbTrialsTroughLocked$Condition[counter+1]     = 'down'
  df_nbTrialsTroughLocked$nb[counter+1]      = nbTrialsTroughLocked$down[idx_sub]
  df_nbTrialsTroughLocked$Sub[counter+2]      = all_sub[idx_sub]
  df_nbTrialsTroughLocked$Condition[counter+2]     = 'not'
  df_nbTrialsTroughLocked$nb[counter+2]      = nbTrialsTroughLocked$not[idx_sub]
 counter  = counter+3
}
sorted_Condition  =  c( 'up' , 'down','not' )
df_nbTrialsTroughLocked$Condition = factor(df_nbTrialsTroughLocked$Condition, levels = sorted_Condition)

ezANOVA (df_nbTrialsTroughLocked, dv = .(nb), wid=.(Sub), within = .(Condition), detailed=T)

t.test(nbTrialsTroughLocked$up,nbTrialsTroughLocked$down,paired=T)
t.test(nbTrialsTroughLocked$up,nbTrialsTroughLocked$not,paired=T)
t.test(nbTrialsTroughLocked$down,nbTrialsTroughLocked$not,paired=T)



nbTrialsCueLocked = read.csv(file = 'nbTrl_Cue-locked_analyzed_onlyTP.txt',header = T,sep = '\t');
t.test(nbTrialsCueLocked$TrueUP,nbTrialsCueLocked$ShamUp,paired=T)
t.test(nbTrialsCueLocked$TrueDown,nbTrialsCueLocked$ShamDown,paired=T)
t.test(nbTrialsCueLocked$TrueUP,nbTrialsCueLocked$TrueDown,paired=T)

df_nbTrialsCueLocked=matrix(,length(all_sub)*4 ,3)
colnames(df_nbTrialsCueLocked)=c('Sub','Condition',"nb")
df_nbTrialsCueLocked=as.data.frame(df_nbTrialsCueLocked)
counter = 1
for (idx_sub in 1:length(all_sub))
{
  df_nbTrialsCueLocked$Condition[counter]     = 'TrueUP'
  df_nbTrialsCueLocked$Sub[counter]      = all_sub[idx_sub]
  df_nbTrialsCueLocked$nb[counter]      = nbTrialsCueLocked$TrueUP[idx_sub]
  
  df_nbTrialsCueLocked$Condition[counter+1]     = 'ShamUp'
  df_nbTrialsCueLocked$Sub[counter+1]      = all_sub[idx_sub]
  df_nbTrialsCueLocked$nb[counter+1]      = nbTrialsCueLocked$ShamUp[idx_sub]
  
  df_nbTrialsCueLocked$Condition[counter+2]     = 'TrueDown'
  df_nbTrialsCueLocked$Sub[counter+2]      = all_sub[idx_sub]
  df_nbTrialsCueLocked$nb[counter+2]      = nbTrialsCueLocked$TrueDown[idx_sub]
  
  
  df_nbTrialsCueLocked$Condition[counter+3]     = 'ShamDown'
  df_nbTrialsCueLocked$Sub[counter+3]      = all_sub[idx_sub]
  df_nbTrialsCueLocked$nb[counter+3]      = nbTrialsCueLocked$ShamDown[idx_sub]
  counter  = counter+4
}

sorted_Condition  =  c( 'TrueUP' , 'ShamUp','TrueDown','ShamDown' )
df_nbTrialsCueLocked$Condition = factor(df_nbTrialsCueLocked$Condition, levels = sorted_Condition)

tapply(df_nbTrialsCueLocked$nb,df_nbTrialsCueLocked$Condition , CI_perso)
ezANOVA (df_nbTrialsCueLocked[df_nbTrialsCueLocked$Condition!='ShamDown',], 
         dv = .(nb), wid=.(Sub), within = .(Condition), detailed=T)



############ SLEEP instrucitons


CI(questionnaires$sleep_duration_estimate_StMary)
CI(questionnaires$sleep_quality_St_Mary)


### Acti + sleep journal
CI(na.omit(questionnaires$N.1))
CI(na.omit(questionnaires$N.2))
CI(na.omit(questionnaires$N.3))


############ Participant characteristics

CI(questionnaires$Edinburgh_handeness)
CI(questionnaires$Daytime_sleepiness)

CI(questionnaires$Beck_depression)
CI(questionnaires$Beck_anxiety)
CI(questionnaires$PQSI)
CI(questionnaires$Chronotype)


### Arousals per condition
allSub = c('sub-01', 'sub-02', 'sub-03', 'sub-04', 'sub-05', 'sub-06', 'sub-07', 'sub-08', 'sub-09', 
           'sub-10', 'sub-11', 'sub-12', 'sub-13', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18', 
           'sub-19', 'sub-20', 'sub-21', 'sub-22', 'sub-23', 'sub-24', 'sub-25', 'sub-26', 'sub-27', 'sub-29', 'sub-30', 'sub-31')
arousalConditions <- read_delim("arousalConditions.csv", 
                                delim = ";", escape_double = FALSE, trim_ws = TRUE)
arousalConditions$Sub = as.factor(arousalConditions$Sub)
arousalConditions$Condition = as.factor(arousalConditions$Condition)

arousalConditions_Count=matrix(,length(allSub)*2,3)
colnames(arousalConditions_Count)=c('Sub','Cond','Count')
arousalConditions_Count=as.data.frame(arousalConditions_Count)

counter = 1
for (idx_sub in 1:length(allSub))
{
  arousalConditions_Count$Cond [counter]  = 'up'
  arousalConditions_Count$Sub [counter]  =allSub[idx_sub]
  arousalConditions_Count$Count [counter]  = length(arousalConditions$Condition[arousalConditions$Sub==allSub[idx_sub] & arousalConditions$Condition=='up'])
  
  arousalConditions_Count$Cond [counter+1]  = 'down'
  arousalConditions_Count$Sub [counter+1]  =allSub[idx_sub]
  arousalConditions_Count$Count [counter+1]  = length(arousalConditions$Condition[arousalConditions$Sub==allSub[idx_sub] & arousalConditions$Condition=='down'])
  counter = counter +2
  
  arousalConditions_Count$Cond [counter+2]  = 'rest'
  arousalConditions_Count$Sub [counter+2]  =allSub[idx_sub]
  arousalConditions_Count$Count [counter+2]  = length(arousalConditions$Condition[arousalConditions$Sub==allSub[idx_sub] & arousalConditions$Condition=='rest'])
}


ezPlot (arousalConditions_Count[!arousalConditions_Count$Cond%in% c('out','rest'),], dv = .(Count), wid=.(Sub), within = .(Cond), x = .(Cond))

DensityUpvsDown = t.test(arousalConditions_Count$Count[arousalConditions_Count$Cond=='up'],arousalConditions_Count$Count[arousalConditions_Count$Cond=='down'], paired=T)
cohensD(arousalConditions_Count$Count[arousalConditions_Count$Cond=='up'],arousalConditions_Count$Count[arousalConditions_Count$Cond=='down'], method = "paired")
CI_perso(arousalConditions_Count$Count[arousalConditions_Count$Cond=='up'])
CI_perso(arousalConditions_Count$Count[arousalConditions_Count$Cond=='down'])

ezANOVA(arousalConditions_Count[!arousalConditions_Count$Cond%in% c('out','rest'),], 
        dv = .(   Count),
        wid=.(Sub), 
        within = .(Cond), 
        detailed=T)

DensityUpvsRest = t.test(arousalConditions_Count$Count[arousalConditions_Count$Cond=='up'],arousalConditions_Count$Count[arousalConditions_Count$Cond=='not'], paired=T)
DensityDownVsRest = t.test(arousalConditions_Count$Count[arousalConditions_Count$Cond=='down'],arousalConditions_Count$Count[arousalConditions_Count$Cond=='not'], paired=T)
p.adjust(c(DensityUpvsDown$p.value,DensityUpvsRest$p.value,DensityDownVsRest$p.value),n=3,method = 'fdr')

#Effectsize computation
cohensD(arousalConditions_Count$Count[arousalConditions_Count$Cond=='up'],arousalConditions_Count$Count[arousalConditions_Count$Cond=='down'], method = "paired")
cohensD(arousalConditions_Count$Count[arousalConditions_Count$Cond=='up'],arousalConditions_Count$Count[arousalConditions_Count$Cond=='rest'], method = "paired")
cohensD(arousalConditions_Count$Count[arousalConditions_Count$Cond=='down'],arousalConditions_Count$Count[arousalConditions_Count$Cond=='rest'], method = "paired")

tapply(arousalConditions_Count$Count, arousalConditions_Count$Cond, CI_perso)


