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


# Main effects : Within Conditions

PrevsPost <- read.csv("D:/Documents/Research/TMR/closedLoop/bimanual_closed_loop/analyses/fmri/Analyses/group/task/rfx_031221/PrevsPost_betaValues_ActivationAnalysis.csv", sep=";")
sorted_Session  = c('pre','post')
PrevsPost$Session = factor(PrevsPost$Session, levels = sorted_Session)
sorted_Condition  = c('up','down','not')
PrevsPost$Condition = factor(PrevsPost$Condition, levels = sorted_Condition)
PrevsPost$ROI_Name = as.factor(PrevsPost$ROI_Name)

###Fig 6a Putamen_Right

plotUp = ggplot(PrevsPost[PrevsPost$ROI_Name=='Putamen_Right' & PrevsPost$Condition=='up',], aes(x=Session , y=Beta, fill = Session )) + 
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
  coord_cartesian(ylim = c(-0.5, 0.7))+
  ggtitle('Up')+
  theme_classic() +
  guides(fill="none")

plotDown = ggplot(PrevsPost[PrevsPost$ROI_Name=='Putamen_Right' & PrevsPost$Condition=='down',], aes(x=Session , y=Beta, fill = Session )) + 
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
  coord_cartesian(ylim = c(-0.5, 0.7))+
  ggtitle('Down')+
  theme_classic() +
  guides(fill="none")

plotNot = ggplot(PrevsPost[PrevsPost$ROI_Name=='Putamen_Right' & PrevsPost$Condition=='not',], aes(x=Session , y=Beta, fill = Session )) + 
  geom_violin(position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c(not1,not3))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = c(not2)) +
  geom_point(aes(fill=Session,group=Sub),size=2,shape=19, position = position_dodge(0.2),
             color=c(not2)) +
  stat_summary(fun=mean, geom="point", shape=18, size=10,aes(group=Session), position=position_dodge(1),
               color=c(not3,not1)) +
  stat_summary(fun=median, geom = "crossbar", width = 0.5,size=0.5,aes(group=Session), position=position_dodge(0.2),
               color=c(not3,not1)) +
  facet_wrap(~ROI_Name)+
  facet_wrap(~ROI_Name,nrow = 1)+
  coord_cartesian(ylim = c(-0.5, 0.7))+
  ggtitle('Not')+
  theme_classic() +
  guides(fill="none")

plot = ggarrange(plotUp, plotDown, plotNot , 
          ncol = 3, nrow = 1)
annotate_figure(plot, top = text_grob("Task-related BOLD WITHIN CONDITION (Putamen Right 26  -8  -4) ", 
                                      color = "red", face = "bold", size = 14))


###Fig 6a M1_Left

plotUp =ggplot(PrevsPost[PrevsPost$ROI_Name=='M1_Left' & PrevsPost$Condition=='up',], aes(x=Session , y=Beta, fill = Session )) + 
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
  coord_cartesian(ylim = c(-0.2, 2.6))+
  ggtitle('Up')+
  theme_classic() +
  guides(fill="none")

plotDown = ggplot(PrevsPost[PrevsPost$ROI_Name=='M1_Left' & PrevsPost$Condition=='down',], aes(x=Session , y=Beta, fill = Session )) + 
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
  coord_cartesian(ylim = c(-0.2, 2.6))+
  ggtitle('Down')+
  theme_classic() +
  guides(fill="none")

plotNot = ggplot(PrevsPost[PrevsPost$ROI_Name=='M1_Left' & PrevsPost$Condition=='not',], aes(x=Session , y=Beta, fill = Session )) + 
  geom_violin(position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c(not1,not3))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = c(not2)) +
  geom_point(aes(fill=Session,group=Sub),size=2,shape=19, position = position_dodge(0.2),
             color=c(not2)) +
  stat_summary(fun=mean, geom="point", shape=18, size=10,aes(group=Session), position=position_dodge(1),
               color=c(not3,not1)) +
  stat_summary(fun=median, geom = "crossbar", width = 0.5,size=0.5,aes(group=Session), position=position_dodge(0.2),
               color=c(not3,not1)) +
  facet_wrap(~ROI_Name)+
  facet_wrap(~ROI_Name,nrow = 1)+
  coord_cartesian(ylim = c(-0.2, 2.6))+
  ggtitle('Not')+
  theme_classic() +
  guides(fill="none")

plot = ggarrange(plotUp, plotDown, plotNot , 
          ncol = 3, nrow = 1)
annotate_figure(plot, top = text_grob("Task-related BOLD WITHIN CONDITION (M1 Left -38 -20  52 ) ", 
                                      color = "red", face = "bold", size = 14))


###Fig 6a Caudate_Right

plotUp = ggplot(PrevsPost[PrevsPost$ROI_Name=='Caudate_Right' & PrevsPost$Condition=='up',], aes(x=Session , y=Beta, fill = Session )) + 
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
  coord_cartesian(ylim = c(-1, 1.2))+
  ggtitle('Up')+
  theme_classic() +
  guides(fill="none")

plotDown = ggplot(PrevsPost[PrevsPost$ROI_Name=='Caudate_Right' & PrevsPost$Condition=='down',], aes(x=Session , y=Beta, fill = Session )) + 
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
  coord_cartesian(ylim = c(-1, 1.2))+
  ggtitle('Down')+
  theme_classic() +
  guides(fill="none")

plotNot = ggplot(PrevsPost[PrevsPost$ROI_Name=='Caudate_Right' & PrevsPost$Condition=='not',], aes(x=Session , y=Beta, fill = Session )) + 
  geom_violin(position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c(not1,not3))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = c(not2)) +
  geom_point(aes(fill=Session,group=Sub),size=2,shape=19, position = position_dodge(0.2),
             color=c(not2)) +
  stat_summary(fun=mean, geom="point", shape=18, size=10,aes(group=Session), position=position_dodge(1),
               color=c(not3,not1)) +
  stat_summary(fun=median, geom = "crossbar", width = 0.5,size=0.5,aes(group=Session), position=position_dodge(0.2),
               color=c(not3,not1)) +
  facet_wrap(~ROI_Name)+
  facet_wrap(~ROI_Name,nrow = 1)+
  coord_cartesian(ylim = c(-1, 1.2))+
  ggtitle('Not')+
  theme_classic() +
  guides(fill="none")

plot = ggarrange(plotUp, plotDown, plotNot , 
                 ncol = 3, nrow = 1)
annotate_figure(plot, top = text_grob("Task-related BOLD WITHIN CONDITION (Caudate_Right  10 0 14) ", 
                                      color = "red", face = "bold", size = 14))

                                      
                                      
##Fig 7a Hippocampus right

plotUp =ggplot(PrevsPost[PrevsPost$ROI_Name=='Hippocampus_Right' & PrevsPost$Condition=='up',], aes(x=Session , y=Beta, fill = Session )) + 
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
  coord_cartesian(ylim = c(-0.5, 0.3))+
  ggtitle('Up')+
  theme_classic() +
  guides(fill="none")

plotDown = ggplot(PrevsPost[PrevsPost$ROI_Name=='Hippocampus_Right' & PrevsPost$Condition=='down',], aes(x=Session , y=Beta, fill = Session )) + 
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
  coord_cartesian(ylim = c(-0.5, 0.3))+
  ggtitle('Down')+
  theme_classic() +
  guides(fill="none")

plotNot = ggplot(PrevsPost[PrevsPost$ROI_Name=='Hippocampus_Right' & PrevsPost$Condition=='not',], aes(x=Session , y=Beta, fill = Session )) + 
  geom_violin(position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c(not1,not3))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = c(not2)) +
  geom_point(aes(fill=Session,group=Sub),size=2,shape=19, position = position_dodge(0.2),
             color=c(not2)) +
  stat_summary(fun=mean, geom="point", shape=18, size=10,aes(group=Session), position=position_dodge(1),
               color=c(not3,not1)) +
  stat_summary(fun=median, geom = "crossbar", width = 0.5,size=0.5,aes(group=Session), position=position_dodge(0.2),
               color=c(not3,not1)) +
  facet_wrap(~ROI_Name)+
  facet_wrap(~ROI_Name,nrow = 1)+
  coord_cartesian(ylim = c(-0.5, 0.3))+
  ggtitle('Not')+
  theme_classic() +
  guides(fill="none")

plot = ggarrange(plotUp, plotDown, plotNot + rremove("x.text"), 
          ncol = 3, nrow = 1)
annotate_figure(plot, top = text_grob("Task-related BOLD WITHIN CONDITION (Hippocampus Right 32 -38  -6) ", 
                                      color = "red", face = "bold", size = 14))



#Interaction : between conditions


diffInteractions=matrix(,length(allCond)*length(allSub)*length(allROI),5)
colnames(diffInteractions)=c(colnames(PrevsPost)[c(1:3,5)],"Diff")
diffInteractions=as.data.frame(diffInteractions)
counter = 1
for (idx_sub in 1:length(allSub))
{
  for (idx_Cond in 1 : length(allCond))
  {
   
    for (idx_ROI in 1 : length(allROI_Coord))
    {
      
      tmp = PrevsPost[PrevsPost$Sub == allSub[idx_sub] & PrevsPost$ROI_Coord == allROI_Coord[idx_ROI] & PrevsPost$Condition==allCond[idx_Cond],]
      
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
diffInteractions$Condition = factor(diffInteractions$Condition, levels = sorted_Cond)
diffInteractions$Condition = factor(diffInteractions$Condition, levels = c('up','down','not'))

#Fig 6b  Right caudate cluster 1: Up vs Down
ggplot(diffInteractions[diffInteractions$ROI_Coord== "20  18  12"  & diffInteractions$Condition!='not',],
       aes(x=Condition , y=Diff,fill=Condition )) + 
  geom_violin(position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c(up3,down3))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = "grey")+
  stat_summary(fun=mean, geom="point", shape=18, size=10, aes(group=Condition), position=position_dodge(1),
               color=c( up1,down1) ) +
  stat_summary(fun=median, geom = "crossbar", width = 0.7,size=0.5,alpha = 0.5, aes(group=Condition), position=position_dodge(1),
               color=c( up1,down1)) +
  geom_point(aes(fill=Condition,group=Sub),size=2,shape=19, position = position_dodge(0.2),
             color = rep(c(down1,up1), times = length(allSub))) +
  facet_wrap(~ROI_Name)+
  coord_cartesian(ylim = c(-0.6, 0.7))+
  theme_classic() 

#Fig 6b  Right caudate cluster 2: Up vs Not

ggplot(diffInteractions[diffInteractions$ROI_Name== "Caudate_R_upVsNot" & diffInteractions$Condition!='down',], 
       aes(x=Condition , y=Diff,fill=Condition )) + 
  geom_violin(position = position_dodge(width = 0.9))+
  scale_fill_manual(values=c( up3,not3))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = "grey") +
  geom_point(aes(fill=Condition,group=Sub),size=2,shape=21, position = position_dodge(0.2)) +
  stat_summary(fun=mean, geom="point", shape=18, size=10,alpha = 0.7, aes(group=Condition), position=position_dodge(1),
               color=c(up1,not1) ) +
  stat_summary(fun=median, geom = "crossbar", width = 0.7,size=0.5,alpha = 0.5, aes(group=Condition), position=position_dodge(1),
               color=c(up1,not1) ) +
  geom_point(aes(fill=Condition,group=Sub),size=2,shape=19, position = position_dodge(0.2),
             color = rep(c(not1,up1), times = length(allSub))) +
  coord_cartesian(ylim = c(-0.6, 0.7))+
  facet_wrap(~ROI_Name)+
  theme_classic() 

#Fig 6b  Right caudate cluster 3: Down vs Not

ggplot(diffInteractions[diffInteractions$ROI_Name== "Caudate_R_DownVsNot" & diffInteractions$Condition!='up',], 
       aes(x=Condition , y=Diff,fill=Condition )) + 
  geom_violin(position = position_dodge(width = 0.9))+
  scale_fill_manual(values=c( down3,not3))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = "grey") +
  geom_point(aes(fill=Condition,group=Sub),size=2,shape=21, position = position_dodge(0.2)) +
  stat_summary(fun=mean, geom="point", shape=18, size=10,alpha = 0.7, aes(group=Condition), position=position_dodge(1),
               color=c(down1,not1) ) +
  stat_summary(fun=median, geom = "crossbar", width = 0.7,size=0.5,alpha = 0.5, aes(group=Condition), position=position_dodge(1),
               color=c(down1,not1) ) +
  geom_point(aes(fill=Condition,group=Sub),size=2,shape=19, position = position_dodge(0.2),
             color = rep(c(down1,not1), times = length(allSub))) +
  coord_cartesian(ylim = c(-0.6, 0.7))+
  facet_wrap(~ROI_Name)+
  theme_classic() 


#Fig 7b  Hippocampus right : Up vs Not


ggplot(diffInteractions[diffInteractions$ROI_Name== "Hippocampus_R_upVsNot" & diffInteractions$Condition!='down',], 
       aes(x=Condition , y=Diff,fill=Condition )) + 
  geom_violin(position = position_dodge(width = 0.9))+
  scale_fill_manual(values=c( up3,not3))+
  geom_line(aes(group=Sub), position = position_dodge(0.2),color = "grey") +
  geom_point(aes(fill=Condition,group=Sub),size=2,shape=21, position = position_dodge(0.2)) +
  stat_summary(fun=mean, geom="point", shape=18, size=10,alpha = 0.7, aes(group=Condition), position=position_dodge(1),
               color=c(up1,not1) ) +
  stat_summary(fun=median, geom = "crossbar", width = 0.7,size=0.5,alpha = 0.5, aes(group=Condition), position=position_dodge(1),
               color=c(up1,not1) ) +
  geom_point(aes(fill=Condition,group=Sub),size=2,shape=19, position = position_dodge(0.2),
             color = rep(c(not1,up1), times = length(allSub))) +
  coord_cartesian(ylim = c(-0.45, 0.3))+
  facet_wrap(~ROI_Name)+
  theme_classic() 


