setwd("D:/Documents/Research/TMR/closedLoop/bimanual_closed_loop/data/group/")
library(ez)  
library(ggplot2)  
library(Rmisc)
library(matlab)

up1=rgb(252,195,255, maxColorValue = 255)
up2=rgb(172,40,174, maxColorValue = 255)
up3=rgb(83,1,83, maxColorValue = 255)
down1=rgb(197,197,255, maxColorValue = 255)
down2=rgb(92,92,255, maxColorValue = 255)
down3=rgb(1,1,131, maxColorValue = 255)
not1=rgb(176,254,254, maxColorValue = 255)
not2=rgb(0,104,103, maxColorValue = 255)
not3=rgb(0,41,41, maxColorValue = 255)

CI_perso<- function (x, ci = 0.95) 
{
  a <- mean(x)
  s <- sd(x)
  n <- length(x)
  error <- qt(ci + (1 - ci)/2, df = n - 1) * s/sqrt(n)
  return(round(c(mean = a, lower = a - error,upper = a + error),1)  )
}

### Random SRTT 


randomSRTT = read.csv(file = 'RandomSRTT.csv',header = T,sep = ';');



randomSRTT$Block   = as.factor(randomSRTT$Block)
randomSRTT$Cue     = as.factor(randomSRTT$Cue)
randomSRTT$Rep     = as.factor(randomSRTT$Rep)
randomSRTT$Sub     = as.factor(randomSRTT$Sub)
randomSRTT$Session = as.factor(randomSRTT$Session)
sorted_Session     = c(" preNight",  " postNight")
randomSRTT$Session = factor(randomSRTT$Session, levels = sorted_Session)

allSess       = levels(randomSRTT$Session)
allBlock      = levels(randomSRTT$Block)
allCue        = levels(randomSRTT$Cue)
allRep        = levels(randomSRTT$Rep)
allSub        = levels(randomSRTT$Sub)


####### OUtlier


for (idx_sub in 1:length(allSub))
{
  for (idx_block in 1 : length(allBlock))
  {
    
    # RT outlyers
    tmpLim = randomSRTT[randomSRTT$Sub==allSub[idx_sub] & randomSRTT$Acc=='1' & 
                          randomSRTT$Block==allBlock[idx_block],]
    
    
    IQR = quantile(tmpLim$RT)[4]-quantile(tmpLim$RT)[2]
    limInf = quantile(tmpLim$RT)[2]-(1.5*IQR)
    limSup = quantile(tmpLim$RT)[4]+(1.5*IQR)
    
    
    randomSRTT$Outlier[randomSRTT$Sub==allSub[idx_sub] & 
                         randomSRTT$Block==allBlock[idx_block]] = randomSRTT$Acc[randomSRTT$Sub==allSub[idx_sub] & 
                                                                                   randomSRTT$Block==allBlock[idx_block]]=='1' & 
      randomSRTT$RT[randomSRTT$Sub==allSub[idx_sub] & 
                      randomSRTT$Block==allBlock[idx_block]]>limInf     & randomSRTT$RT[randomSRTT$Sub==allSub[idx_sub] & 
                                                                                          randomSRTT$Block==allBlock[idx_block]]<limSup 
    
  }
}


####### Compute Accuracy by Block, mean rt per key press PI 

randomSummary=matrix(,length(allSub)*length(allBlock)*length(allSess),9)
colnames(randomSummary)=c("Sub", "Session","Block","Sequence","Condition","Sound","Median","Mean", "Acc")
randomSummary=as.data.frame(randomSummary)
counter = 1

for (idx_sub in 1:length(allSub))
{
  for (idx_sess in 1 : length(allSess))
  {
    
    for (idx_block in 1 : length(allBlock))
    {
      
      # RT outlyers
      tmp = randomSRTT[randomSRTT$Sub==allSub[idx_sub] & 
                         randomSRTT$Session==allSess[idx_sess] & 
                         randomSRTT$Block==allBlock[idx_block] ,]
      IQR = quantile(tmp$RT)[4]-quantile(tmp$RT)[2]
      limInf = quantile(tmp$RT)[2]-(1.5*IQR)
      limSup = quantile(tmp$RT)[4]+(1.5*IQR)
      
      
      randomSummary$Sub[counter]       = allSub[idx_sub]
      if (allSess[idx_sess] ==  " postNight")
      {
        # randomSummary$Block[counter]     = as.numeric(allBlock[idx_block])+41  
        randomSummary$Block[counter]     = as.numeric(allBlock[idx_block]) 
        randomSummary$Session[counter]     = 3
        
      }
      if (allSess[idx_sess] ==  " preNight")
      {
        randomSummary$Block[counter]     = allBlock[idx_block]
        randomSummary$Session[counter]     = 1
      }
      
      # randomSummary$Condition[counter]   = allSess[idx_sess]
      randomSummary$Condition[counter]   = "random"
      randomSummary$Sequence[counter]   = "NA"
      randomSummary$Sound[counter]   = "NA"
      
      
      # % Accuracy per block
      randomSummary$Acc[counter]       = sum(tmp$Acc)/length(tmp$Acc)
      
      #Mean RT per key presses
      randomSummary$Mean[counter]      = mean(tmp$RT[tmp$Acc=='1' & tmp$RT>limInf & tmp$RT<limSup ])
      randomSummary$Median[counter]    = median(tmp$RT[tmp$Acc=='1' & tmp$RT>limInf & tmp$RT<limSup ])
      
      counter = counter+1 
    }
    
  }
}


randomSummary$Sub       = as.factor(randomSummary$Sub)
randomSummary$Condition = as.factor(randomSummary$Condition)
randomSummary$Block     = as.factor(randomSummary$Block)
sorted_Block     = c('1', '2', '3', '4', '42', '43', '44', '45')
randomSummary$Block = factor(randomSummary$Block, levels = sorted_Block)


#### Sequential SRTT (MSL)

MSLData = read.csv(file = 'sequentialSRTT_tmp.csv',header = T,sep = ';')

MSLData$RT=MSLData$RT*1000
MSLData$Acc=MSLData$Cue==MSLData$Rep
MSLData$Acc = as.numeric(MSLData$Acc)
MSLData$Block      = as.factor(MSLData$Block)
MSLData$Session    = as.factor(MSLData$Session)
MSLData$OrdinalPos = as.factor(MSLData$OrdinalPos)
MSLData$Cue        = as.factor(MSLData$Cue)
MSLData$Rep        = as.factor(MSLData$Rep)
MSLData$Sound      = as.factor(MSLData$Sound)
MSLData$Order      = as.factor(MSLData$Order)
MSLData$Repetition = as.factor(MSLData$Repetition)
MSLData$Condition = as.factor(MSLData$Condition)
MSLData$Sub       = as.factor(MSLData$Sub)
MSLData$Sequence  = as.factor(MSLData$Sequence)

sorted_Block  = paste(sort(as.integer(levels(MSLData$Block))))
MSLData$Block = factor(MSLData$Block, levels = sorted_Block)

sorted_Session  = c('1',  '2',  '3')
MSLData$Session = factor(MSLData$Session, levels = sorted_Session)

sorted_Condition  =  c( ' up' , ' down',' not' )
MSLData$Condition = factor(MSLData$Condition, levels = sorted_Condition)

allSub        = levels(MSLData$Sub)
allSess       = levels(MSLData$Session)
allBlock      = levels(MSLData$Block)
allSequence   = levels(MSLData$Sequence)
allCond       = levels(MSLData$Condition)
allOrdinalPos = levels(MSLData$OrdinalPos)
allRepetition = levels(MSLData$Repetition)
allCue        = levels(MSLData$Cue)
allRep        = levels(MSLData$Rep)


####### OUtlier


for (idx_sub in 1:length(allSub))
{
  for (idx_block in 1 : length(allBlock))
  {
      
    # RT outlyers
    tmpLim = MSLData[MSLData$Sub==allSub[idx_sub] & MSLData$Acc=='1' & 
                       MSLData$Block==allBlock[idx_block],]
    

    IQR = quantile(tmpLim$RT)[4]-quantile(tmpLim$RT)[2]
    limInf = quantile(tmpLim$RT)[2]-(1.5*IQR)
    limSup = quantile(tmpLim$RT)[4]+(1.5*IQR)

    
    MSLData$Outlier[MSLData$Sub==allSub[idx_sub] & 
                      MSLData$Block==allBlock[idx_block]] = MSLData$Acc[MSLData$Sub==allSub[idx_sub] & 
                                                                          MSLData$Block==allBlock[idx_block]]=='1' & 
      MSLData$RT[MSLData$Sub==allSub[idx_sub] & 
                   MSLData$Block==allBlock[idx_block]]>limInf     & MSLData$RT[MSLData$Sub==allSub[idx_sub] & 
                   MSLData$Block==allBlock[idx_block]]<limSup 
    
  }
}

#percentage of outlier in the correct trials (acc == 1 and outlier == False)
length(MSLData$RT[MSLData$Acc=='1' & MSLData$Outlier=='FALSE'])/length(MSLData$RT[MSLData$Acc=='1' ])*100
#percentage of incorrect trials in all  trials (acc == 0 )
length(MSLData$RT[MSLData$Acc=='0' ])/length(MSLData$RT)*100

#percentage of discarded trials (acc == 0 and outlier == False) in total
(length(MSLData$RT[ MSLData$Outlier=='FALSE' ] )/length(MSLData$RT))*100




####### Compute Accuracy by Block, mean rt per key press PI 


MSLSummary=matrix(,length(allSub)*length(allBlock)*length(allSequence),10)
colnames(MSLSummary)=c(colnames(MSLData)[c(1:7)],"Mean","Median", "Acc")
MSLSummary=as.data.frame(MSLSummary)
counter = 1
for (idx_sub in 1:length(allSub))
{
  for (idx_block in 1 : length(allBlock))
  {

    for (idx_seq in 1:length(allSequence))
    {
      tmp = MSLData[MSLData$Sub==allSub[idx_sub] & 
                      MSLData$Block==allBlock[idx_block] &
                      MSLData$Sequence==allSequence[idx_seq] 
                       ,]

      MSLSummary$Sub[counter]       = allSub[idx_sub]
      MSLSummary$Block[counter]     = allBlock[idx_block]
      MSLSummary$Sequence[counter]  = allSequence[idx_seq]
      MSLSummary$Session[counter]   = as.character(tmp$Session[1])
      MSLSummary$Condition[counter] = as.character(tmp$Condition[1])
      MSLSummary$Order[counter]     = as.character(tmp$Order[1])
      
      if (tmp$Sound[1]=="1")
      {MSLSummary$Sound[counter] ='Low pitch'}
      
      if (tmp$Sound[1]=="2")
      {MSLSummary$Sound[counter] ='White noise'}
      
      if (tmp$Sound[1]=="3")
      {MSLSummary$Sound[counter] ='High pitch'}
      
      
      # % Accuracy per block
      MSLSummary$Acc[counter]       = sum(tmp$Acc)/length(tmp$Acc)
      MSLData$PerCorr[MSLData$Sub==allSub[idx_sub] & 
                        MSLData$Block==allBlock[idx_block] &
                        MSLData$Sequence==allSequence[idx_seq]] = repmat(sum(tmp$Acc[tmp$Outlier==T& tmp$Acc==1])/length(tmp$Acc),length(tmp$Sequence),1)
      
      #Mean RT per key presses
      MSLSummary$Mean[counter]      = mean(tmp$RT[tmp$Outlier==T & tmp$Acc==1 ])
      MSLSummary$Median[counter]    = median(tmp$RT[tmp$Outlier==T & tmp$Acc==1])

      counter = counter+1
      
    }
  }
}


MSLSummary$Sub       = as.factor(MSLSummary$Sub)
MSLSummary$Session   = as.factor(MSLSummary$Session)
MSLSummary$Block     = as.factor(MSLSummary$Block)
MSLSummary$Sequence  = as.factor(MSLSummary$Sequence)
MSLSummary$Condition = as.factor(MSLSummary$Condition)

sorted_Block  = paste(sort(as.integer(levels(MSLSummary$Block))))
MSLSummary$Block = factor(MSLSummary$Block, levels = sorted_Block)
MSLSummary$Session = factor(MSLSummary$Session, levels = sorted_Session)
MSLSummary$Condition = factor(MSLSummary$Condition, levels = sorted_Condition)



####### Plot raw by block and Condition  Fig 2a
#RT
ggplot(rbind(summarySE(MSLSummary, measurevar="Median", groupvars=c(  "Block", "Condition"),na.rm=T),summarySE(randomSummary, measurevar="Median", groupvars=c(  "Block","Condition"),na.rm=T)),
       aes(x=Block, y=Median, group = Condition )) + 
  geom_point(aes(color=Condition,shape=Condition),size=5)+
  geom_line(aes(color=Condition))+
  scale_shape_manual(values=c(16, 1,18,16,16))+
  scale_size_manual(values=c(2,2,2,2,2))+
  geom_ribbon(aes(ymin = Median-se,
                  ymax = Median+se,fill=factor(Condition)),alpha = 0.2)+
  scale_color_manual(values=c(up2,down2,not2,"Black","Black")) +
  scale_fill_manual(values=c(up1, down2, not2,"Black","Black")) +
  ylab("Reaction time (ms)")+
  xlab("Block of training")+
  # coord_cartesian(xlim=c(22,27),ylim=c(0,450))+
  theme_classic()


#Acc Fig S9 a
ggplot(rbind(summarySE(MSLSummary, measurevar="Acc", groupvars=c(  "Block", "Condition"),na.rm=T),summarySE(randomSummary, measurevar="Acc", groupvars=c(  "Block","Condition"),na.rm=T)),
       aes(x=Block, y=Acc, group = Condition )) + 
  geom_point(aes(color=Condition,shape=Condition),size=5)+
  geom_line(aes(color=Condition))+
  scale_shape_manual(values=c(16, 1,18,16,16))+
  scale_size_manual(values=c(2,2,2,2,2))+
  geom_ribbon(aes(ymin = Acc-se,
                  ymax = Acc+se,fill=factor(Condition)),alpha = 0.2)+
  scale_color_manual(values=c(up2,down2,not2,"Black","Black")) +
  scale_fill_manual(values=c(up1, down2, not2,"Black","Black")) +
  ylab("Reaction time (ms)")+
  xlab("Block of training")+
  theme_classic()


###### Plot raw by block and Sequence  
#RT Fig S7 a
ggplot(summarySE(MSLSummary[MSLSummary$Block %in% as.factor(c(1:24)),], measurevar="Median", 
                 groupvars=c(  "Block","Sequence"),na.rm=T),
       aes(x=Block, y=Median, group = Sequence )       ) + 
  geom_point(aes(color=Sequence,shape=Sequence),size=2)+
  geom_line(aes(color=Sequence))+
  scale_shape_manual(values=c(16, 17,18,16,16))+
  scale_size_manual(values=c(2,2,2,2,2))+
  geom_ribbon(aes(ymin = Median-se,
                  ymax = Median+se,fill=factor(Sequence)),alpha = 0.2)+
  ylab("Reaction time (ms)")+
  coord_cartesian(ylim = c(300,650))+
  xlab("Block of training")+
  theme_classic()

#Acc Fig S7 b
ggplot(summarySE(MSLSummary[MSLSummary$Block %in% as.factor(c(1:24)),], measurevar="Acc", 
                 groupvars=c(  "Block","Sequence"),na.rm=T),
       aes(x=Block, y=Acc, group = Sequence )       ) + 
  geom_point(aes(color=Sequence,shape=Sequence),size=2)+
  geom_line(aes(color=Sequence))+
  scale_shape_manual(values=c(16, 17,18,16,16))+
  scale_size_manual(values=c(2,2,2,2,2))+
  geom_ribbon(aes(ymin = Acc-se,
                  ymax = Acc+se,fill=factor(Sequence)),alpha = 0.2)+
  ylab("Reaction time (ms)")+
  xlab("Block of training")+
  coord_cartesian(ylim = c(0.84,1))+
  theme_classic()

