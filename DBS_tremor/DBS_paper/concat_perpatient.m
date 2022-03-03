

clear all
ad=struct; ad.metrics=cell(1,1); ad.spirals=cell(1,1); ad.sig=cell(1,1); ad.signames={'data1';'data2';'time';'envelope'};ad.metricnames={'subj';'cond';'trials';'nspirals';'nquadrants';'mean env';'peak2peak';'variance'};

close all
cond={'NS';'HF';'C'};
cohort=[1 3 4 6];


for iii=2:4
    clearvars -except ad iii cohort cond signal_all map
    concat=[];
    for co=1:3
        cn=1;
        if (co==1 && iii==2 | iii==3)
            ntrial=[1 2];
        else
            ntrial=1;
        end
        
        for trial=1:length(ntrial)
            load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/P0',num2str(cohort(iii)),'_clean_',num2str(cond{co,1}),num2str(trial),'_SH.mat'));
            
                        
            centre=[535 361];
            signal1(1,:)=signal(1,:)-centre(1);
            signal1(2,:)=signal(2,:)-centre(2);
            signal1(3,:)= sqrt((signal1(1,:).^2)+(signal1(2,:).^2));
            
            [peak]=dindfpeak(signal1(3,:),samplerate2);
            [a,b]=  butter(2, [(peak-2)/(0.5*samplerate2) (peak+2)/(0.5*samplerate2)], 'bandpass'); %15
            filt_cs=(filtfilt(a,b,signal1(3,:)));
            
            if isempty(concat)
                beg=1;
            else
            beg=find(concat==concat(end));
            end
            
            concat=[concat filt_cs];
            
            
           map{iii-1,co}(trial,:)=[beg length(filt_cs)+beg-1];
           
           clear signal1 peak a b filt_cs signal
            
        end
    end
     signal_all{iii-1,1}=zscore(concat);
end
clearvars -except signal_all cond map
cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/')

save('zscore_acr_cond.mat')

