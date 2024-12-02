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

### Regressions hippocampus seed

Reg_PPI_rightHippo <- read.csv("D:/Documents/Research/TMR/closedLoop/bimanual_closed_loop/analyses/fmri/Analyses/group/task/rfx_031221/regressions_PPI_right_hippocampus.csv", sep=";")
sorted_Session  = c('pre','post')
Reg_PPI_rightHippo$Session = factor(Reg_PPI_rightHippo$Session, levels = sorted_Session)
sorted_Condition  = c('down','up')
Reg_PPI_rightHippo$Condition = factor(Reg_PPI_rightHippo$Condition, levels = sorted_Condition)
Reg_PPI_rightHippo$Sub       = as.factor(Reg_PPI_rightHippo$Sub)
Reg_PPI_rightHippo$ROI_Name  = as.factor(Reg_PPI_rightHippo$ROI_Name)
Reg_PPI_rightHippo$ROI_Coord = as.factor(Reg_PPI_rightHippo$ROI_Coord)
Reg_PPI_rightHippo$Condition = as.factor(Reg_PPI_rightHippo$Condition)



allSess       = levels(Reg_PPI_rightHippo$Session)
allSub        = levels(Reg_PPI_rightHippo$Sub)
allROI        = levels(Reg_PPI_rightHippo$ROI_Name)
allROI_Coord  = levels(Reg_PPI_rightHippo$ROI_Coord)
allCond       = levels(Reg_PPI_rightHippo$Condition) 

Reg__PPI_diffSession_rightHippo=matrix(,length(Reg_PPI_rightHippo[,1])/2,8)
colnames(Reg__PPI_diffSession_rightHippo)=c(colnames(Reg_PPI_rightHippo)[c(1:3,5)],"Diff_Beta",colnames(Reg_PPI_rightHippo)[c(6,9,10)])
Reg__PPI_diffSession_rightHippo=as.data.frame(Reg__PPI_diffSession_rightHippo)
counter = 1
for (idx_sub in 1:length(allSub))
{
  
  for (idx_ROI in 1 : length(allROI_Coord))
  {
    for (idx_Cond in 1 : length(allCond))
    {
      
      tmp = Reg_PPI_rightHippo[Reg_PPI_rightHippo$Sub == allSub[idx_sub] & Reg_PPI_rightHippo$ROI_Coord == allROI_Coord[idx_ROI] & Reg_PPI_rightHippo$Condition == allCond[idx_Cond] ,]
      
      Reg__PPI_diffSession_rightHippo$Sub[counter]       = allSub[idx_sub]
      Reg__PPI_diffSession_rightHippo$Condition[counter] = as.character(tmp$Condition[1])
      Reg__PPI_diffSession_rightHippo$ROI_Name[counter]  = as.character(tmp$ROI_Name[1])
      Reg__PPI_diffSession_rightHippo$ROI_Coord[counter] = as.character(tmp$ROI_Coord[1])
      Reg__PPI_diffSession_rightHippo$Diff_Beta[counter] = na.omit(tmp$Beta[tmp$Session=='post'])-na.omit(tmp$Beta[tmp$Session=='pre'])
      Reg__PPI_diffSession_rightHippo$TMR_idx[counter] = tmp$TMR_idx[1]
      Reg__PPI_diffSession_rightHippo$Sigma_power[counter] = tmp$Sigma_power[1]
      Reg__PPI_diffSession_rightHippo$Peak_amplitude[counter] = tmp$Peak_amplitude[1]
      counter = counter+1
    }  
  }
}



Reg__PPI_diffSession_rightHippo$Sub       = as.factor(Reg__PPI_diffSession_rightHippo$Sub)
Reg__PPI_diffSession_rightHippo$ROI_Name  = as.factor(Reg__PPI_diffSession_rightHippo$ROI_Name)
Reg__PPI_diffSession_rightHippo$ROI_Coord = as.factor(Reg__PPI_diffSession_rightHippo$ROI_Coord)
Reg__PPI_diffSession_rightHippo$Condition = as.factor(Reg__PPI_diffSession_rightHippo$Condition)


# Left Hippocampus connectivity x TMR index  (Fig 8d )

ggplot(Reg__PPI_diffSession_rightHippo[Reg__PPI_diffSession_rightHippo$ROI_Coord=='-32 -18  42' & Reg__PPI_diffSession_rightHippo$Condition=='down',], aes(x = Peak_amplitude, y = Diff_Beta, color = Condition)) +
  geom_point(size =2,alpha = 0.5 ) +
  scale_color_manual(values = c(down2)) +
  geom_smooth(method = lm,alpha=0.2,size = 1) +
  theme_minimal() +
  ggtitle("right hippocampus Within condition overnight change x TMR index in the ")+
  coord_cartesian(ylim = c(-1, 1))+
  # stat_cor(method = "spearman")+
  theme(axis.line = element_line(color = "grey70"))


### Regressions Caudate seed (figures be)
Reg_PPI_rightCaudate <- read.csv("D:/Documents/Research/TMR/closedLoop/bimanual_closed_loop/analyses/fmri/Analyses/group/task/rfx_031221/regressions_PPI_right_caudate.csv", sep=";")
sorted_Session  = c('pre','post')
Reg_PPI_rightCaudate$Session = factor(Reg_PPI_rightCaudate$Session, levels = sorted_Session)
sorted_Condition  =  c( 'up' , 'down', 'not' )
Reg_PPI_rightCaudate$Condition = factor(Reg_PPI_rightCaudate$Condition, levels = sorted_Condition)
Reg_PPI_rightCaudate$Sub       = as.factor(Reg_PPI_rightCaudate$Sub)
Reg_PPI_rightCaudate$ROI_Name  = as.factor(Reg_PPI_rightCaudate$ROI_Name)
Reg_PPI_rightCaudate$ROI_Coord = as.factor(Reg_PPI_rightCaudate$ROI_Coord)
Reg_PPI_rightCaudate$Condition = as.factor(Reg_PPI_rightCaudate$Condition)


allSess       = levels(Reg_PPI_rightCaudate$Session)
allSub        = levels(Reg_PPI_rightCaudate$Sub)
allROI        = levels(Reg_PPI_rightCaudate$ROI_Name)
allROI_Coord  = levels(Reg_PPI_rightCaudate$ROI_Coord)
allCond       = levels(Reg_PPI_rightCaudate$Condition) 

Reg__PPI_diffSession_rightCaudate=matrix(,length(Reg_PPI_rightCaudate[,1])/2,8)
colnames(Reg__PPI_diffSession_rightCaudate)=c(colnames(Reg_PPI_rightCaudate)[c(1:3,6)],"Diff_Beta",colnames(Reg_PPI_rightCaudate)[c(8:10)])
Reg__PPI_diffSession_rightCaudate=as.data.frame(Reg__PPI_diffSession_rightCaudate)
counter = 1
for (idx_sub in 1:length(allSub))
{
  
  for (idx_ROI in 1 : length(allROI_Coord))
  {
    for (idx_Cond in 1 : length(allCond))
    {
      
      tmp = Reg_PPI_rightCaudate[Reg_PPI_rightCaudate$Sub == allSub[idx_sub] & Reg_PPI_rightCaudate$ROI_Coord == allROI_Coord[idx_ROI] & Reg_PPI_rightCaudate$Condition == allCond[idx_Cond] ,]
      
      Reg__PPI_diffSession_rightCaudate$Sub[counter]       = allSub[idx_sub]
      Reg__PPI_diffSession_rightCaudate$Condition[counter] = as.character(tmp$Condition[1])
      Reg__PPI_diffSession_rightCaudate$ROI_Name[counter]  = as.character(tmp$ROI_Name[1])
      Reg__PPI_diffSession_rightCaudate$ROI_Coord[counter] = as.character(tmp$ROI_Coord[1])
      Reg__PPI_diffSession_rightCaudate$Diff_Beta[counter] = na.omit(tmp$Beta[tmp$Session=='post'])-na.omit(tmp$Beta[tmp$Session=='pre'])
      Reg__PPI_diffSession_rightCaudate$TMR_idx[counter] = tmp$TMR_idx[1]
      Reg__PPI_diffSession_rightCaudate$Sigma_power[counter] = tmp$Sigma_power[1]
      Reg__PPI_diffSession_rightCaudate$Peak_amplitude[counter] = tmp$Peak_amplitude[1]
      Reg__PPI_diffSession_rightCaudate$Contrast[counter] = as.character(tmp$Contrast[1])
      counter = counter+1
    }  
  }
}



Reg__PPI_diffSession_rightCaudate$Sub       = as.factor(Reg__PPI_diffSession_rightCaudate$Sub)
Reg__PPI_diffSession_rightCaudate$ROI_Name  = as.factor(Reg__PPI_diffSession_rightCaudate$ROI_Name)
Reg__PPI_diffSession_rightCaudate$ROI_Coord = as.factor(Reg__PPI_diffSession_rightCaudate$ROI_Coord)
Reg__PPI_diffSession_rightCaudate$Condition = as.factor(Reg__PPI_diffSession_rightCaudate$Condition)

# Left Hippocampus connectivity x TMR index  (Fig 9d right panel)

ggplot(Reg__PPI_diffSession_rightCaudate[ Reg__PPI_diffSession_rightCaudate$ROI_Name == 'hippocampus_left' & Reg__PPI_diffSession_rightCaudate$Condition=='down',], 
       aes(x = TMR_idx, y = Diff_Beta, color = Condition)) +
  geom_point(size =2,alpha = 0.5 ) +
  scale_color_manual(values = c(down2)) +
  geom_smooth(method = lm,alpha=0.2,size = 1) +
  theme_minimal() +
  ggtitle("CAUDATE-HIPPO x TMR Down (hippocampus left -16 -40   6)")+
  coord_cartesian(ylim = c(-1.3, 1), xlim = c(-70,30))+
  ylab('Overnight change (post - pre)') +
  # stat_cor(method = "pearson")+
  theme(axis.line = element_line(color = "grey70"))

## Regressions Putamen seed 
Reg_PPI_rightPutamen <- read.csv("D:/Documents/Research/TMR/closedLoop/bimanual_closed_loop/analyses/fmri/Analyses/group/task/rfx_031221/regressions_PPI_right_putamen.csv", sep=";")
sorted_Session  = c('pre','post')
Reg_PPI_rightPutamen$Session = factor(Reg_PPI_rightPutamen$Session, levels = sorted_Session)
sorted_Condition  = c('down','up')
Reg_PPI_rightPutamen$Condition = factor(Reg_PPI_rightPutamen$Condition, levels = sorted_Condition)
Reg_PPI_rightPutamen$Sub       = as.factor(Reg_PPI_rightPutamen$Sub)
Reg_PPI_rightPutamen$ROI_Name  = as.factor(Reg_PPI_rightPutamen$ROI_Name)
Reg_PPI_rightPutamen$ROI_Coord = as.factor(Reg_PPI_rightPutamen$ROI_Coord)
Reg_PPI_rightPutamen$Condition = as.factor(Reg_PPI_rightPutamen$Condition)


allSess       = levels(Reg_PPI_rightPutamen$Session)
allSub        = levels(Reg_PPI_rightPutamen$Sub)
allROI        = levels(Reg_PPI_rightPutamen$ROI_Name)
allROI_Coord  = levels(Reg_PPI_rightPutamen$ROI_Coord)
allCond       = levels(Reg_PPI_rightPutamen$Condition) 

Reg__PPI_diffSession_rightPutamen=matrix(,length(Reg_PPI_rightPutamen[,1])/2,8)
colnames(Reg__PPI_diffSession_rightPutamen)=c(colnames(Reg_PPI_rightPutamen)[c(1:3,5)],"Diff_Beta",colnames(Reg_PPI_rightPutamen)[c(7:9)])
Reg__PPI_diffSession_rightPutamen=as.data.frame(Reg__PPI_diffSession_rightPutamen)
counter = 1
for (idx_sub in 1:length(allSub))
{
  
  for (idx_ROI in 1 : length(allROI_Coord))
  {
    for (idx_Cond in 1 : length(allCond))
    {
      
      tmp = Reg_PPI_rightPutamen[Reg_PPI_rightPutamen$Sub == allSub[idx_sub] & Reg_PPI_rightPutamen$ROI_Coord == allROI_Coord[idx_ROI] & Reg_PPI_rightPutamen$Condition == allCond[idx_Cond] ,]
      
      Reg__PPI_diffSession_rightPutamen$Sub[counter]       = allSub[idx_sub]
      Reg__PPI_diffSession_rightPutamen$Condition[counter] = as.character(tmp$Condition[1])
      Reg__PPI_diffSession_rightPutamen$ROI_Name[counter]  = as.character(tmp$ROI_Name[1])
      Reg__PPI_diffSession_rightPutamen$ROI_Coord[counter] = as.character(tmp$ROI_Coord[1])
      Reg__PPI_diffSession_rightPutamen$Diff_Beta[counter] = na.omit(tmp$Beta[tmp$Session=='post'])-(tmp$Beta[tmp$Session=='pre'])
      Reg__PPI_diffSession_rightPutamen$TMR_idx[counter] = tmp$TMR_idx[1]
      Reg__PPI_diffSession_rightPutamen$Sigma_power[counter] = tmp$Sigma_power[1]
      Reg__PPI_diffSession_rightPutamen$Peak_amplitude[counter] = tmp$Peak_amplitude[1]
      counter = counter+1
    }  
  }
}



Reg__PPI_diffSession_rightPutamen$Sub       = as.factor(Reg__PPI_diffSession_rightPutamen$Sub)
Reg__PPI_diffSession_rightPutamen$ROI_Name  = as.factor(Reg__PPI_diffSession_rightPutamen$ROI_Name)
Reg__PPI_diffSession_rightPutamen$ROI_Coord = as.factor(Reg__PPI_diffSession_rightPutamen$ROI_Coord)
Reg__PPI_diffSession_rightPutamen$Condition = as.factor(Reg__PPI_diffSession_rightPutamen$Condition)

# Left right M1 connectivity x TMR index  (Fig 9d left panel)

ggplot(Reg__PPI_diffSession_rightPutamen[Reg__PPI_diffSession_rightPutamen$ROI_Coord =='26  -8  44' & Reg__PPI_diffSession_rightPutamen$Condition =='down' ,], 
       aes(x = TMR_idx, y = Diff_Beta, color = Condition)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("darkblue","magenta")) +
  geom_smooth(method = lm,alpha=0.2) +
  theme_minimal() +
  ggtitle("Within condition PPI overnight change x TMR index in the connectivity betwwen the right putamen and the right caudate")+
  theme(axis.line = element_line(color = "grey70"))



#Left Hippocampus connectivity x sigma power  (Fig 9c left panel with caudate seed)

ggplot(rbind(Reg__PPI_diffSession_rightCaudate[Reg__PPI_diffSession_rightCaudate$Contrast=='Sigmadown' & Reg__PPI_diffSession_rightCaudate$Condition=='down' & Reg__PPI_diffSession_rightCaudate$ROI_Name=='hippocampus_left',], 
             Reg__PPI_diffSession_rightPutamen[Reg__PPI_diffSession_rightPutamen$Contrast=='Sigmadown' & Reg__PPI_diffSession_rightPutamen$Condition=='down' ,]),
       aes(x = Sigma_power, y = Diff_Beta, color = ROI_Coord)) +
  geom_point(size =2,alpha = 0.5 ) +
  scale_color_manual(values = c(down2,down3)) +
  geom_smooth(method = lm,alpha=0.2,size = 1) +
  theme_minimal() +
  ggtitle("CAUDATE-HIPPO x Sigmma down (hippocampus right 14 -30 -10 ) ")+
  coord_cartesian(ylim = c(-1.7, 1.5))+
  ylab('Overnight change (post - pre)') +
  # stat_cor(method = "spearman")+
  theme(axis.line = element_line(color = "grey70"))


#Left aSPL connectivity x SO peak amplitude  (Fig 9c right panel with caudate seed)
ggplot(rbind(Reg__PPI_diffSession_rightCaudate[Reg__PPI_diffSession_rightCaudate$Contrast=='Peakdown' & Reg__PPI_diffSession_rightCaudate$Condition=='down' & Reg__PPI_diffSession_rightCaudate$ROI_Name=='hippocampus_right_shifted' ,], 
             Reg__PPI_diffSession_rightPutamen[Reg__PPI_diffSession_rightPutamen$Contrast=='Peakdown' & Reg__PPI_diffSession_rightPutamen$Condition=='down' & Reg__PPI_diffSession_rightPutamen$ROI_Name=='hippocampus_right',]),
       aes(x = Peak_amplitude, y = Diff_Beta, color = ROI_Coord)) +
  geom_point(size =2,alpha = 0.5) +
  scale_color_manual(values = c(down2,down3)) +
  geom_smooth(method = lm,alpha=0.2,size = 1) +
  theme_minimal() +
  ggtitle("CAUDATE-HIPPO x Peak down (hippocampus right 14 -32 -8 ) ")+
  coord_cartesian(ylim = c(-1.7, 1.5),xlim = c(0, 17))+
  ylab('Overnight change (post - pre)') +
  # stat_cor(method = "spearman")+
  theme(axis.line = element_line(color = "grey70"))+
  guides(fill="none")

