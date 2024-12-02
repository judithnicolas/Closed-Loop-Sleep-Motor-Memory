%% get group data from log file
% Judith Nicolas
% Created 2020 at KU Leuven

[listSub,listSubEEG,listSubBehav] = getFullDatasets;


clc

pathIn = initPath.Exp;

dirOutput= [pathIn 'data\group\'];
dirInput= [pathIn 'data\' ];
nbSessionMSL = 2; %including training and test
keyPressesRandom = 60;
BEHAV_0_getPVT(listSub,dirOutput,dirInput);
BEHAV_0_getRandomSRTT(listSub,dirOutput,dirInput,nbSessionMSL,keyPressesRandom );
BEHAV_0_getSequentialSRTT(listSubBehav,dirInput,dirOutput);


