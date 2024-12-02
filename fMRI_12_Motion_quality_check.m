function [rpcheck_trn_min_max,high_motion_volumes ]= fMRI_12_Motion_quality_check(pathIn,sub)

%% SPM12 Spatial Preprocessing pipeline - Stress Study - Leuven 12/07/2018 - GA

Dir.data  = [pathIn '\analyses\fmri\Analyses\'];

cd(Dir.data)

% data.subjects   = {  sub    };
data.subjects   =   sub    ;

for sub = 1:1:length(data.subjects)   

       data.sessions (sub, 1:3) =     { ...
        'prenight/training', ...
        'prenight/test', ...
        'postnight/training',...
    }';

end


rp_files = [];


for sub_counter = 1 : 1 : length(data.subjects)
    for sess_counter = 1 : 1 : size(data.sessions,2)     
     
        cd([Dir.data data.subjects{sub_counter} '/' data.sessions{sub_counter,sess_counter} '/'])
        files = dir('*.txt');
        rp_files {sub_counter, sess_counter} = [files.folder '/' files.name]; 
        % [rp_files {[files.folder '/' files.name]}]; 
        clear files; 
        
    end
end
clear sub_counter; clear sess_counter;
cd(Dir.data)

rpcheck_trn_min_max = [];
rpcheck_rot_min_max = [];
rpcheck_trn_disp = [];
rpcheck_rot_disp = [];
high_motion_volumes={}
number_high_motion_volumes_trans = [];

for sub_counter = 1:1:length(data.subjects) 
    fig=figure('Units','centimeters','Position',[0 0 25 25]); 
    for sess_counter = 1:1:size(data.sessions,2)  
        [rp] = load(rp_files{sub_counter,sess_counter});
        high_motion_vol_translation = NaN(900,6); 
        nvol = size(rp,1); 
        
        %% identify high-motion volumes (3 mm) for x y z Translation 
        for i = 1:1:3 % voor parameter x,y,z tranlation 
            for ii = 1:nvol % for volume 1:nvol
%                 if rp(ii,i)<-5 || rp(ii,i)>5 
                if rp(ii,i)<-4 || rp(ii,i)>4
                    high_motion_vol_translation(ii,i)=1;
                else
                    high_motion_vol_translation(ii,i)=0;end
            end
        end
        
        for i = 4:1:6 %radians to degrees 
                for ii = 1:length(rp(:,4))
                    rp(ii,i) = radtodeg(rp(ii,i));
                end
        end
        clear i; clear ii; 
        
        subplot (3,1,sess_counter)  
        plot(1:nvol,rp(:,1))
        hold on  
        plot(1:nvol,rp(:,2))
        hold on
        plot(1:nvol,rp(:,3))
        hold on
        ylabel('Translation(mm)');xlim([1 nvol]); xlabel('Timepoints');
        legend('X translation','Y translation','Z translation') 
        title([char(data.subjects(sub_counter)),' ',char(data.sessions(sub_counter,sess_counter))]);
        hold on
        
        %%% save minima and maxima 
        high_motion_volumes.translation{sub_counter,sess_counter} (1,1) = {find(high_motion_vol_translation(:,1)==1)}; 
        high_motion_volumes.translation{sub_counter,sess_counter} (1,2) = {find(high_motion_vol_translation(:,2)==1)}; 
        high_motion_volumes.translation{sub_counter,sess_counter} (1,3) = {find(high_motion_vol_translation(:,3)==1)}; 
        
        number_high_motion_volumes_trans = [number_high_motion_volumes_trans; nansum(high_motion_vol_translation(:,1)) nansum(high_motion_vol_translation(:,2)) nansum(high_motion_vol_translation(:,3))];
       
        rpcheck_trn_min_max = [rpcheck_trn_min_max; sub_counter sess_counter min(min(rp(:,1:3))) max(max(rp(:,1:3))) max(max(abs(rp(:,1:3))))];
        rpcheck_rot_min_max = [rpcheck_rot_min_max; sub_counter sess_counter min(min(rp(:,4:6))) max(max(rp(:,4:6))) max(max(abs(rp(:,4:6))))];
        rpcheck_trn_disp = [rpcheck_trn_disp; sub_counter sess_counter max(max(rp(:,1)))-min(min(rp(:,1))) max(max(rp(:,2)))-min(min(rp(:,2))) max(max(rp(:,3)))-min(min(rp(:,3)))];
        rpcheck_rot_disp = [rpcheck_rot_disp; sub_counter sess_counter max(max(rp(:,4)))-min(min(rp(:,4))) max(max(rp(:,5)))-min(min(rp(:,5))) max(max(rp(:,6)))-min(min(rp(:,6)))];
        clear nvol; clear npar; clear rp; clear high_motion_vol_translation;      
    end
    
    saveas(fig,char(strcat([Dir.data '\Realignment_checks\taskBasedAnalysis\',data.subjects{sub_counter}, '_trn'],'.jpg')));
%     close(fig);
    fig=figure('Units','centimeters','Position',[0 0 25 25]); 
    for sess_counter = 1:1:size(data.sessions,2)  
        [rp] = load(rp_files{sub_counter,sess_counter});
        nvol = size(rp,1); 
        for i = 4:1:6
                for ii = 1:length(rp(:,4))
                    rp(ii,i) = radtodeg(rp(ii,i));
                end
        end
        clear i; clear ii; 
        subplot (3,1,sess_counter)  
        plot(1:nvol,rp(:,4))
        hold on  
        plot(1:nvol,rp(:,5))
        hold on
        plot(1:nvol,rp(:,6))
        hold on
        
        ylabel('Rotation(deg)');xlim([1 nvol]); xlabel('Timepoints');
        legend('X rotation','Y rotation','Z rotation') 
        title([char(data.subjects(sub_counter)),' ',char(data.sessions(sub_counter,sess_counter))]);
        hold on
        clear nvol; clear npar; clear rp;            
    end
    saveas(fig,char(strcat([Dir.data '\Realignment_checks\taskBasedAnalysis\',data.subjects{sub_counter}, '_rot'],'.jpg')));
%     close(fig);
    clear sess_counter; 
end

clear sub_counter; 

save(char(strcat([Dir.data '\Realignment_checks\taskBasedAnalysis\Realignment_check.mat'])),'rpcheck_trn_min_max','rpcheck_trn_disp','rpcheck_rot_min_max','rpcheck_rot_disp','data','number_high_motion_volumes_trans','high_motion_volumes');
% clear all;
fprintf('%s','motion check done\n')

