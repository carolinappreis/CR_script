cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')

clear all
ad=struct; ad.metrics=cell(1,1);ad.metricnames={'subj';'cond';'trials';'mean env';'peak2peak';'variance'};

load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/ACC_zscore_acr_cond.mat'))

%%% normalised to NS condition: load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/ACC_spiral_NSnorm.mat'))

for iii=1:3
    clearvars -except ad iii signal_all map
    for co=1:3
        env=abs(hilbert(signal_all{iii}));
        
        for trial= 1:size(map{iii,co},1)
            env_dum=env(map{iii,co}(trial,1):map{iii,co}(trial,end));
            filt_dum=signal_all{iii}(map{iii,co}(trial,1):map{iii,co}(trial,end));
               
            samplerate2=1000;
            bins=1:samplerate2/2:length(filt_dum)-1;
            for ii=1:length(bins)-1
                pp(1,ii)=peak2peak(filt_dum(bins(ii):bins(ii+1)));
            end
            
%             ad.metrics{iii,co}(trial,:) = [mean(env_dum); nanmean(pp) ; var(filt_dum)];
            
            ad.metrics{iii,co}(trial,:) = [mean(env_dum)];
            clear bins pp  env_dum filt_dum
        end
    end
    close all
end



% % cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA')
% % clearvars -except ad 
% % 
% % filename=strcat('spiral_zsc.mat');
% % % % % % % % save(filename)
% % 
% % load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/spiral_zsc.mat')

%%%% order of conditions here NS/HFS/PSS. Values imported from met variable
%%%% to table in SPSS for stats at group level.
met=cell(3,3);

for i=1:3
    for ii=1:3
        met{i,ii}=ad.metrics{i,ii}; %%%% met{pt,cond}(n_spirals,n_quadrants,metrics)
    end
    subplot(3,1,i)
    bar([mean(met{i,1}) mean(met{i,2}) mean(met{i,3})]);
    change_from_NS(i,:)=([(mean(met{i,1})-mean(met{i,3})) (mean(met{i,1})-mean(met{i,2}))])
 

end

