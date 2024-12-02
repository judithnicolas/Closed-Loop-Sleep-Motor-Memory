%% Initialization
[listSub,listSubEEG,listSubBehav] = getFullDatasets;
[TMRup,TMRdown,offlineGainUp,offlineGainDown,offlineGainNot,offlineGainUp_Training,offlineGainDown_Training,offlineGainNot_Training,TMRup_Training,TMRdown_Training]= getBehaviour();

[powerSigmaUp,powerSigmaDown,peakAmplitudeUp,peakAmplitudeDown]= getEEG();


pathIn=initPath.Exp;

date = '031221';

%% Pre-processing

%Extrtact scans where to  cut
scanExtract(pathIn,sub)

%Extract Onset
extractOnsets(pathIn,sub)

%Copy files from RawOrganized to Analyses
copyfile ([pathIn 'analyses\fmri\RawOrganized\' sub ], [pathIn 'analyses\fmri\Analyses\' sub '\task\' ])

% Run python script
system(['python ' pathIn 'analyses\fmri\scripts\splitter.py ' sub ' ' pathIn ])

%run Preprocessing
fMRI_11_preprocessing_MNI(pathIn,sub)
 
%% Run motion quality check
[rpcheck_trn_min_max,high_motion_volumes ]= fMRI_12_Motion_quality_check(pathIn,listSubBehav);
 
%% model_specification_estimation

pathIn=initPath.Exp;
date = '031221';
 for idx_sub = 1 : length(listSubBehav)
     sub = listSubBehav{idx_sub};
     fprintf('%s\n',sub)
     
     fMRI_13_model_specification_estimation(pathIn,sub,date)
 end
 
 %% Contrasts
 
 for idx_sub = 1 : length(listSubBehav)
     sub = listSubBehav{idx_sub};
     fprintf('%s\n',sub)
     mod = 0; % 0/1 => without/with behavioral modulator
     fMRI_14_extractContrasts(pathIn, sub,date,mod) 
     mod = 1; % 0/1 => without/with behavioral modulator
     fMRI_14_extractContrasts(pathIn, sub,date,mod) 
 end
 

 
  %% Run ffx
 
 mod = 0; % 0/1 => without/with behavioral modulator
 for idx_sub = 1%: length(listSubBehav)
     sub = listSubBehav{idx_sub};
     mkdir([ pathIn '\analyses\fmri\Analyses\JOBS\analysisFFx' date '\'])
     fMRI_21_ffx_MNI(sub,pathIn,date,mod)
 end
 
  %% Run Rfx
mod = 0;
fMRI_0_smoothAllContrastImages(pathIn,listSubBehav(1),date,mod,[])
fMRI_23_rfx_MNI(pathIn,listSubBehav,date,mod)  

 
  %% Run Rfx cov
mod = 0; % 0/1 => without/with behavioral modulator
type = 'behaviour/TMRIndex';
extension = [''];

condition = 'up';
contrast_of_interest = [9]; 
fMRI_24_rfx_MNI_cov(pathIn,listSubBehav,date,mod,TMRup,type,condition,contrast_of_interest,extension)  

condition = 'down';
contrast_of_interest = [11]; 
fMRI_24_rfx_MNI_cov(pathIn,listSubBehav,date,mod,TMRdown,type,condition,contrast_of_interest,extension)  


listSubCorEEGBehav = intersect(listSubEEG,listSubBehav, 'stable');
type = 'eeg/sigmaPower';
condition = 'up';
contrast_of_interest = [9];
fMRI_24_rfx_MNI_cov(pathIn,listSubCorEEGBehav,date,mod,powerSigmaUp(1:length(listSubCorEEGBehav)),type,condition,contrast_of_interest,extension)

condition = 'down';
contrast_of_interest = [11]; 
fMRI_24_rfx_MNI_cov(pathIn,listSubCorEEGBehav,date,mod,powerSigmaDown(1:length(listSubCorEEGBehav)),type,condition,contrast_of_interest,extension) 


type = 'eeg/peakAmplitude';
condition = 'up';
contrast_of_interest = [9]; 
fMRI_24_rfx_MNI_cov(pathIn,listSubCorEEGBehav,date,mod,peakAmplitudeUp(1:length(listSubCorEEGBehav)),type,condition,contrast_of_interest,extension)  

condition = 'down';
contrast_of_interest = [11]; 
fMRI_24_rfx_MNI_cov(pathIn,listSubCorEEGBehav,date,mod,peakAmplitudeDown(1:length(listSubCorEEGBehav)),type,condition,contrast_of_interest,extension) 

 
 %% VOI and PPI Computation
contrastname = 'NOT-UP'; 
voiname = {};
voicoord = {};
voiname{1} = 'right_putamen';   voicoord{1} = [18 12 -2];
voiname{2} = 'right_caudate';   voicoord{2} = [10 14 12];
voiname{3} = 'right_hippocampus';   voicoord{3} = [32 -38 -6];

contrast_nb = [261 1];%{'261_preNightTraining_VS_rest_Against_up_VS_down'} and {'1_postNightTraining_VS_rest_Against_up_VS_down'}         
fMRI_31_VOI_PPI_MNI(pathIn,listSubBehav,date,voiname,voicoord,contrastname,contrast_nb) 

voiname = {};
voicoord = {};
contrastname = 'NOT-UP'; 
contrast_nb = [221 51];%{'221_preNightTraining_VS_rest_Against_not_VS_up'} and {'51_postNightTraining_VS_rest_Against_not_VS_up'}
fMRI_31_VOI_PPI_MNI(pathIn,listSubBehav,date,voiname,voicoord,contrastname,contrast_nb) 


%% Run PPI ffx

contrastname = 'UP-DOWN'; 
voiname = {};
voiname{1} = 'right_putamen';   
voiname{2} = 'right_caudate';  

fMRI_32_FFX_PPI_MNI(pathIn,listSubBehav,date,voiname,contrastname)

contrastname = 'NOT-UP'; 
voiname = {};
voiname{1} = 'right_hippocampus';   
fMRI_32_FFX_PPI_MNI(pathIn,listSubBehav,date,voiname,contrastname)

 %% Run Rfx PPI

 mod = 0;
voiname = {};
voiname{1} = 'right_putamen';   
voiname{2} = 'left_putamen';    
voiname{3} = 'right_caudate';  
voiname{4} = 'left_caudate';    
voiname{5} = 'right_hippocampus';   
voiname{6} = 'left_hippocampus';    


for idx_voi = 1 : length(voiname)
    extension = ['\PPI_analysis\' voiname{idx_voi } ];

    fMRI_0_smoothAllContrastImages(pathIn,listSubBehav,date,mod,extension)
end

for idx_voi = 1 : length(voiname)
    extension = ['\PPI_analysis\' voiname{idx_voi } ];
    
    
    rfx_CLTMR_MNI(pathIn,listSubBehav,date,mod,extension)
    
end


%% Run rfx cov PPI 


mod = 0;
voiname = {};
voiname{1} = 'right_putamen';   
voiname{2} = 'right_caudate';  
voiname{3} = 'right_hippocampus';   

type = 'behaviour/TMRIndex';


condition = 'up';
contrast_of_interest = [13 ];
for idx_voi = 1: length(voiname)
    extension = ['\PPI_analysis\' voiname{idx_voi } ];
    fMRI_32_FFX_PPI_MNI(pathIn,listSubBehav,date,mod,TMRup,type,condition,contrast_of_interest,extension) % cov  = gain per condition to script
end

condition = 'down';
contrast_of_interest = [15];
for idx_voi = 1 : length(voiname)
    extension = ['\PPI_analysis\' voiname{idx_voi } ];
    fMRI_32_FFX_PPI_MNI(pathIn,listSubBehav,date,mod,TMRdown,type,condition,contrast_of_interest,extension) % cov  = gain per condition to script
end


listSubCorEEGBehav = intersect(listSubEEG,listSubBehav, 'stable');

type = 'eeg/sigmaPower';

condition = 'up';
contrast_of_interest = [13 ];
for idx_voi = 1: length(voiname)
    extension = ['\PPI_analysis\' voiname{idx_voi } ];
    fMRI_32_FFX_PPI_MNI(pathIn,listSubBehav,date,mod,powerSigmaUp(1:length(listSubCorEEGBehav)),type,condition,contrast_of_interest,extension) % cov  = gain per condition to script
end

condition = 'down';
contrast_of_interest = [15];
for idx_voi = 1 : length(voiname)
    extension = ['\PPI_analysis\' voiname{idx_voi } ];
    fMRI_32_FFX_PPI_MNI(pathIn,listSubBehav,date,mod,powerSigmaDown(1:length(listSubCorEEGBehav)),type,condition,contrast_of_interest,extension) % cov  = gain per condition to script
end

type = 'eeg/peakAmplitude';

condition = 'up';
contrast_of_interest = [13 ];
for idx_voi = 1: length(voiname)
    extension = ['\PPI_analysis\' voiname{idx_voi } ];
    fMRI_32_FFX_PPI_MNI(pathIn,listSubBehav,date,mod,peakAmplitudeUp(1:length(listSubCorEEGBehav)),type,condition,contrast_of_interest,extension) % cov  = gain per condition to script
end

condition = 'down';
contrast_of_interest = [15];
for idx_voi = 1 : length(voiname)
    extension = ['\PPI_analysis\' voiname{idx_voi } ];
    fMRI_32_FFX_PPI_MNI(pathIn,listSubBehav,date,mod,peakAmplitudeDown(1:length(listSubCorEEGBehav)),type,condition,contrast_of_interest,extension) % cov  = gain per condition to script
end