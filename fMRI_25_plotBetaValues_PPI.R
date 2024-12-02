setwd("D:/Documents/Research/TMR/closedLoop/bimanual_closed_loop/analyses/fmri/Analyses/group/task/rfx_031221/")
library(ggplot2)
library(ggpubr)
library(ez)
library(Rmisc)
library(readr)
library(lme4)
library(REdaS)


up1=rgb(252,195,255, maxColorValue = 255)
up2=rgb(172,40,174, maxColorValue = 255)
up3=rgb(83,1,83, maxColorValue = 255)
down1=rgb(197,197,255, maxColorValue = 255)
down2=rgb(92,92,255, maxColorValue = 255)
down3=rgb(1,1,131, maxColorValue = 255)
not1=rgb(176,254,254, maxColorValue = 255)
not2=rgb(0,104,103, maxColorValue = 255)
not3=rgb(0,41,41, maxColorValue = 255)


### Hippocampus seed
PPI_rightHippo <- read.csv("D:/Documents/Research/TMR/closedLoop/bimanual_closed_loop/analyses/fmri/Analyses/group/task/rfx_031221/PPI_right_hippocampus.csv", sep=";")
sorted_Session  = c('pre','post')
PPI_rightHippo$Session = factor(PPI_rightHippo$Session, levels = sorted_Session)
sorted_Condition  = c('not','down','up')
PPI_rightHippo$Condition = factor(PPI_rightHippo$Condition, levels = sorted_Condition)
PPI_rightHippo$Sub       = as.factor(PPI_rightHippo$Sub)
PPI_rightHippo$ROI_Name  = as.factor(PPI_rightHippo$ROI_Name)
PPI_rightHippo$ROI_Coord = as.factor(PPI_rightHippo$ROI_Coord)
PPI_rightHippo$Condition = as.factor(PPI_rightHippo$Condition)

  # Right PMC within condition (fig 8c)

ggplot(PPI_rightHippo[PPI_rightHippo$ROI_Name=='pmc_right' & PPI_rightHippo$Condition=='up',], aes(x=Session , y=Beta, fill = Session )) + 
  geom_violin(position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c(up1,up3))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color =up2) +
  geom_point(aes(fill=Session,group=Sub),size=2,shape=19, position = position_dodge(0.2),
             color=c(up2)) +
  stat_summary(fun=mean, geom="point", shape=18, size=10,aes(group=Session), position=position_dodge(1),
               color=c(up3,up1)) +
  stat_summary(fun=median, geom = "crossbar", width = 0.5,size=0.5,aes(group=Session), position=position_dodge(0.2),
               color=c(up3,up1)) +
  facet_wrap(~ROI_Name)+
  facet_wrap(~ROI_Name,nrow = 1)+
  # coord_cartesian(ylim = c(-0.5, 0.7))+
  ggtitle('Up')+
  theme_classic() 


allSess       = levels(PPI_rightHippo$Session)
allSub        = levels(PPI_rightHippo$Sub)
allROI        = levels(PPI_rightHippo$ROI_Name)
allROI_Coord  = levels(PPI_rightHippo$ROI_Coord)
allCond       = levels(PPI_rightHippo$Condition)


  # Right M1 Between Conditions (Fig 8b)

diffInteractions=matrix(,length(allCond)*length(allSub)*length(allROI_Coord),5)
colnames(diffInteractions)=c(colnames(PPI_rightHippo)[c(1:4)],"Diff")
diffInteractions=as.data.frame(diffInteractions)
counter = 1
for (idx_sub in 1:length(allSub))
{
   
    for (idx_ROI in 1 : length(allROI_Coord))
    {
      for (idx_Cond in 1 : length(allCond))
      {
        
       tmp = PPI_rightHippo[PPI_rightHippo$Sub == allSub[idx_sub] & PPI_rightHippo$ROI_Coord == allROI_Coord[idx_ROI] & PPI_rightHippo$Condition==allCond[idx_Cond],]
        
        diffInteractions$Sub[counter]       = allSub[idx_sub]
        diffInteractions$Condition[counter] = as.character(tmp$Condition[1])
        diffInteractions$ROI_Name[counter]  = as.character(tmp$ROI_Name[1])
        diffInteractions$ROI_Coord[counter] = as.character(tmp$ROI_Coord[1])
        diffInteractions$Diff[counter] = tmp$Beta[tmp$Session=='post' ]-tmp$Beta[tmp$Session=='pre']
        counter = counter+1
    }  
  }
}



diffInteractions$Sub       = as.factor(diffInteractions$Sub)
diffInteractions$ROI_Name  = as.factor(diffInteractions$ROI_Name)
diffInteractions$ROI_Coord = as.factor(diffInteractions$ROI_Coord)
diffInteractions$Condition = as.factor(diffInteractions$Condition)
diffInteractions$Condition = factor(diffInteractions$Condition, levels = c('up','down','not'))

ggplot(diffInteractions[diffInteractions$ROI_Coord== "52 -24  40"  ,],
       aes(x=Condition , y=Diff,fill=Condition )) + 
  geom_violin(position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c(up3,down3,not3))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = "grey")+
  stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Condition), position=position_dodge(1),
               color=c( up2,down2,not2) ) +
  stat_summary(fun=median, geom = "crossbar", width = 0.7,size=0.5,alpha = 0.5, aes(group=Condition), position=position_dodge(1),
               color=c( up2,down2,not2)) +
  geom_point(aes(fill=Condition,group=Sub),size=2,shape=19, position = position_dodge(0.2),
             color = rep(c(up1,down1,not1), times = length(allSub))) +
  facet_wrap(~ROI_Name)+
  # coord_cartesian(ylim = c(-0.6, 0.7))+
  theme_classic() 



### Putamen seed


PPI_rightPutamen <- read.csv("D:/Documents/Research/TMR/closedLoop/bimanual_closed_loop/analyses/fmri/Analyses/group/task/rfx_031221/PPI_right_putamen.csv", sep=";")
sorted_Session  = c('pre','post')
PPI_rightPutamen$Session = factor(PPI_rightPutamen$Session, levels = sorted_Session)
sorted_Condition  =  c( 'up' , 'down','not' )
PPI_rightPutamen$Condition = factor(PPI_rightPutamen$Condition, levels = sorted_Condition)
PPI_rightPutamen$Sub       = as.factor(PPI_rightPutamen$Sub)
PPI_rightPutamen$ROI_Name  = as.factor(PPI_rightPutamen$ROI_Name)
PPI_rightPutamen$ROI_Coord = as.factor(PPI_rightPutamen$ROI_Coord)
PPI_rightPutamen$Condition = as.factor(PPI_rightPutamen$Condition)

  # right M1 Within condition (Fig 9a)

ggplot(PPI_rightPutamen[PPI_rightPutamen$ROI_Name=='M1_right'  & PPI_rightPutamen$Condition=='down',], aes(x=Session , y=Beta, fill = Session )) + 
  geom_violin(position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c(down1,down3))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color =down2) +
  geom_point(aes(fill=Session,group=Sub),size=2,shape=19, position = position_dodge(0.2),
             color=c(down2)) +
  stat_summary(fun=mean, geom="point", shape=18, size=10,aes(group=Session), position=position_dodge(1),
               color=c(down3,down1)) +
  stat_summary(fun=median, geom = "crossbar", width = 0.5,size=0.5,aes(group=Session), position=position_dodge(0.2),
               color=c(down3,down1)) +
  facet_wrap(~ROI_Name)+
  facet_wrap(~ROI_Name,nrow = 1)+
  coord_cartesian(ylim = c(-0.6, 0.7))+
  ggtitle('Down')+
  theme_classic() +
  guides(fill="none")


  # left M1 Within condition (Fig 9b)


allSess       = levels(PPI_rightPutamen$Session)
allSub        = levels(PPI_rightPutamen$Sub)
allROI        = levels(PPI_rightPutamen$ROI_Name)
allROI_Coord  = levels(PPI_rightPutamen$ROI_Coord)
allCond       = levels(PPI_rightPutamen$Condition)


diffInteractions_puta=matrix(,length(allCond)*length(allSub),5)
colnames(diffInteractions_puta)=c(colnames(Interactions)[c(1:3,5)],"Diff")
diffInteractions_puta=as.data.frame(diffInteractions_puta)
counter = 1
for (idx_sub in 1:length(allSub))
{
  for (idx_Cond in 1 : length(allCond))
  {
    
    for (idx_ROI in 1 : length(allROI_Coord[1]))
    {
      
      tmp = PPI_rightPutamen[PPI_rightPutamen$Sub == allSub[idx_sub] & PPI_rightPutamen$ROI_Coord == allROI_Coord[1] & PPI_rightPutamen$Condition==allCond[idx_Cond],]
      
      diffInteractions_puta$Sub[counter]       = allSub[idx_sub]
      diffInteractions_puta$Condition[counter] = as.character(tmp$Condition[1])
      diffInteractions_puta$ROI_Name[counter]  = as.character(tmp$ROI_Name[1])
      diffInteractions_puta$ROI_Coord[counter] = as.character(tmp$ROI_Coord[1])
      diffInteractions_puta$Diff[counter] = tmp$Beta[tmp$Session=='post' ]-tmp$Beta[tmp$Session=='pre']
      counter = counter+1
    }  
  }
}



diffInteractions_puta$Sub       = as.factor(diffInteractions_puta$Sub)
diffInteractions_puta$ROI_Name  = as.factor(diffInteractions_puta$ROI_Name)
diffInteractions_puta$ROI_Coord = as.factor(diffInteractions_puta$ROI_Coord)
diffInteractions_puta$Condition = as.factor(diffInteractions_puta$Condition)
diffInteractions_puta$Condition = factor(diffInteractions_puta$Condition, levels = sorted_Condition)

ggplot(diffInteractions_puta[diffInteractions_puta$ROI_Coord== "-32 -20  50"  ,],
       aes(x=Condition , y=Diff,fill=Condition )) + 
  geom_violin(position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c(up3,down3,not3))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = "grey")+
  stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Condition), position=position_dodge(1),
               color=c( up1,down1,not1) ) +
  stat_summary(fun=median, geom = "crossbar", width = 0.7,size=1,alpha = 0.5, aes(group=Condition), position=position_dodge(1),
               color=c(up1,down1,not1)) +
  geom_point(aes(fill=Condition,group=Sub),size=2,shape=19, position = position_dodge(0.2),
             color = rep(c(up1,down1,not1), times = length(allSub))) +
  facet_wrap(~ROI_Name)+
  ggtitle("PUTA-MOTOR OVERNIGHT CHANGE UP VS NOT (precentral left -32 -20 50) ")+
  ylab('Overnight change (post - pre)') +
  coord_cartesian(ylim = c(-1, 1))+
  theme_classic() +
  guides(fill="none")

