setwd("D:/Documents/Research/TMR/closedLoop/bimanual_closed_loop/analyses/fmri/Analyses/group/task/rfx_031221/")
library(ggplot2)
library(ggpubr)
library(ez)
library(Rmisc)
library(readr)
library(lme4)
library(REdaS)
library(cocor)

up1=rgb(252,195,255, maxColorValue = 255)
up2=rgb(172,40,174, maxColorValue = 255)
up3=rgb(83,1,83, maxColorValue = 255)
down1=rgb(197,197,255, maxColorValue = 255)
down2=rgb(92,92,255, maxColorValue = 255)
down3=rgb(1,1,131, maxColorValue = 255)
not1=rgb(176,254,254, maxColorValue = 255)
not2=rgb(0,104,103, maxColorValue = 255)
not3=rgb(0,41,41, maxColorValue = 255)


Reg_TMR <- read.csv("D:/Documents/Research/TMR/closedLoop/bimanual_closed_loop/analyses/fmri/Analyses/group/task/rfx_031221/regressions_activation_raw.csv", sep=";")
sorted_Session  = c('pre','post')
Reg_TMR$Session = factor(Reg_TMR$Session, levels = sorted_Session)
sorted_Condition  = c('up','down')
Reg_TMR$Condition = factor(Reg_TMR$Condition, levels = sorted_Condition)

Reg_TMR$Sub       = as.factor(Reg_TMR$Sub)
Reg_TMR$ROI_Name  = as.factor(Reg_TMR$ROI_Name)
Reg_TMR$ROI_Coord = as.factor(Reg_TMR$ROI_Coord)
Reg_TMR$Condition = as.factor(Reg_TMR$Condition)

allSess       = levels(Reg_TMR$Session)
allSub        = levels(Reg_TMR$Sub)
allROI        = levels(Reg_TMR$ROI_Name)
allROI_Coord  = levels(Reg_TMR$ROI_Coord)
allCond       = levels(Reg_TMR$Condition) 

Reg_diffSession=matrix(,length(Reg_TMR[,1])/2,8)
colnames(Reg_diffSession)=c(colnames(Reg_TMR)[c(1:3,5)],"Diff_Beta",colnames(Reg_TMR)[c(7:9)])
Reg_diffSession=as.data.frame(Reg_diffSession)
counter = 1
for (idx_sub in 1:length(allSub))
{
  
  for (idx_ROI in 1 : length(allROI_Coord))
  {
    for (idx_Cond in 1 : length(allCond))
    {
    
      tmp = Reg_TMR[Reg_TMR$Sub == allSub[idx_sub] & Reg_TMR$ROI_Coord == allROI_Coord[idx_ROI] & Reg_TMR$Condition == allCond[idx_Cond] ,]
      
      Reg_diffSession$Sub[counter]       = allSub[idx_sub]
      Reg_diffSession$Condition[counter] = as.character(tmp$Condition[1])
      Reg_diffSession$ROI_Name[counter]  = as.character(tmp$ROI_Name[1])
      Reg_diffSession$ROI_Coord[counter] = as.character(tmp$ROI_Coord[1])
      Reg_diffSession$Diff_Beta[counter] = tmp$Beta[tmp$Session=='post']-tmp$Beta[tmp$Session=='pre']
      Reg_diffSession$TMR_idx[counter] = tmp$TMR_idx[1]
      Reg_diffSession$Sigma_power[counter] = tmp$Sigma_power[1]
      Reg_diffSession$Peak_amplitude[counter] = tmp$Peak_amplitude[1]
      counter = counter+1
    }  
  }
}



Reg_diffSession$Sub       = as.factor(Reg_diffSession$Sub)
Reg_diffSession$ROI_Name  = as.factor(Reg_diffSession$ROI_Name)
Reg_diffSession$ROI_Coord = as.factor(Reg_diffSession$ROI_Coord)
Reg_diffSession$Condition = as.factor(Reg_diffSession$Condition)

# TMR index  
# right Caudate (Fig 6a)
ggplot(Reg_diffSession[(Reg_diffSession$ROI_Name=='caudate_right_up' & Reg_diffSession$Condition=='up') | (Reg_diffSession$ROI_Name=='caudate_right_down' & Reg_diffSession$Condition=='down'),],
       aes(x = TMR_idx, y = Diff_Beta, color = Condition)) +
  geom_point(size =2,alpha = 0.5 , show.legend = FALSE) +
  scale_color_manual(values = c(down2,up2)) +
  geom_smooth(method = lm,alpha=0.2,size = 1, show.legend = FALSE) +
  theme_minimal() +
  ggtitle("right caudate Within condition overnight change x TMR index in the ")+
  coord_cartesian(ylim = c(-0.7, 0.7), xlim = c(-70,30))+
  theme(axis.line = element_line(color = "grey70"))


# right Hippocampus (Fig 7c)
ggplot(Reg_diffSession[(Reg_diffSession$ROI_Name=='hippocampus_right_up' & Reg_diffSession$Condition=='up') | (Reg_diffSession$ROI_Name=='hippocampus_right_down' & Reg_diffSession$Condition=='down'),],
       aes(x = TMR_idx, y = Diff_Beta, color = Condition)) +
  geom_point(size =2,alpha = 0.5 , show.legend = FALSE) +
  scale_color_manual(values = c(down2,up2)) +
  geom_smooth(method = lm,alpha=0.2,size = 1, show.legend = FALSE) +
  theme_minimal() +
  ggtitle("right hippocampus Within condition overnight change x TMR index in the ")+
  coord_cartesian(ylim = c(-0.7, 0.7), xlim = c(-70,30))+
  theme(axis.line = element_line(color = "grey70"))


# Peak Amplitude
# Left M1 (Fig 6b top panel)

ggplot(Reg_diffSession[(Reg_diffSession$ROI_Name=='M1_L_up' & Reg_diffSession$Condition=='up') | (Reg_diffSession$ROI_Name=='M1_L_down' & Reg_diffSession$Condition=='down'),], 
       aes(x = Peak_amplitude , y = Diff_Beta, color = Condition)) +
  geom_point(size =2,alpha = 0.5 , show.legend = FALSE) +
  scale_color_manual(values = c(down2,up2)) +
  geom_smooth(method = lm,alpha=0.2,size = 1, show.legend = FALSE) +
  theme_minimal() +
  ggtitle("left M1 Within condition overnight change x Peak Amplitude  ")+
  coord_cartesian(ylim = c(-1.2, 1.2),xlim = c(0, 17))+
  theme(axis.line = element_line(color = "grey70"))

#right Pallidum (Fig 6b bottom panel)

ggplot(Reg_diffSession[Reg_diffSession$ROI_Name=='Pallidum_R_up' & Reg_diffSession$Condition=='up' ,], 
       aes(x = Peak_amplitude , y = Diff_Beta, color = Condition)) +
  geom_point(size =2,alpha = 0.5 , show.legend = FALSE) +
  scale_color_manual(values = c(up2)) +
  geom_smooth(method = lm,alpha=0.2,size = 1, show.legend = FALSE) +
  theme_minimal() +
  ggtitle("Right Pallidum Within condition overnight change x Peak Amplitude  ")+
  coord_cartesian(ylim = c(-0.5, 0.5),xlim = c(0, 17))+
  theme(axis.line = element_line(color = "grey70"))



