% cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')

clear all; close all
cohort = [1 3 4 6];
cond={'NS';'HFS';'PLS_S'};
folder={'/DATA/';'/DATA_hf/';'/DATA_pls_spiral/'};
norm=1; %%% 0= non zscore 1 = zscore across condtions and 2= normalised to NS condition

for iii=2:4
    concat=[];
    for co=1:size(cond,1)
        load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM',num2str(folder{co,1}),'P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
        
        
        [start, ending, tremor_ds, samplerate]=match(SmrData,iii,co); bg=start; eg=ending;
        clear start ending
        
        run('time_split.m')
        if isempty(t_spi{iii,2,co})
            
            
            if iii==3 && co==3
                che=bg(1);
                clear bg eg
                bg=3522;
                dum=bg+(t_spi{iii,1,co}.*10);
                x=find((dum(:,1))<che);
                if ~isempty(x)
                    dum(1:x(end),:)=[];
                end
                
            else
                dum=bg+(t_spi{iii,1,co}.*10);
            end
            start=dum(:,1); ending=dum(:,2);
            
            
            if iii==4 && co==3
                ending(end)=150600;
            end
            
            if iii==4 && co==2
                start(2)=198300;
            end
            
            if iii==2 && co==3
                clear start ending
                start=[85070 95150 109300 136500 158700 289300 305900 363300 376500];
                ending=[95150 109300 136500 158700 177800 305900 328200 376500 394000];
            end
            %
            
        else
            dum=[bg(1)+(t_spi{iii,1,co}.*10) ; bg(2)+(t_spi{iii,2,co}.*10)];
            start=dum(:,1); ending=dum(:,2);
            
        end
        
        [f1]=spirals_cond(tremor_ds,samplerate,start,ending);
        
        
        y=tremor_ds;
        y(4,:)=sqrt((y(1,:).^2+y(2,:).^2+y(3,:).^2));
        
        if co==1
            norm_NS(iii-1,1)=mean(y(4,:));
            norm_NS(iii-1,2)=std(y(4,:));
        end
        
        handup = [];
        for ix = 1:length(start)
            handup = [handup start(ix):ending(ix)]; %#ok<*AGROW>
        end
        handup = sort(handup,'ascend');
        
        
        [peak]=dindfpeak(y(4,handup),samplerate);
        
        if (peak-2)>=1
            [b,a]=butter(2,[(peak-2)/(0.5*samplerate) (peak+2)/(0.5*samplerate)],'bandpass'); %15
        else
            [b,a]=butter(2,[(1)/(0.5*samplerate) (peak+2)/(0.5*samplerate)],'bandpass'); %15
        end
        %
        %         [b,a]=butter(2,[(1)/(0.5*samplerate) 8/(0.5*samplerate)],'bandpass'); %15
        
        if norm==2
            
            y(4,:)=(y(4,:)-norm_NS(iii-1,1))/norm_NS(iii-1,2);
            signal1=(filtfilt(b,a,y(4,:)));
            
        else
            
            signal1=(filtfilt(b,a,y(4,:)));
            
        end
        
        %         figure(iii-1)
        %         subplot(3,1,1)
        %         plot(signal1)
        
        for trial=1:length(start)
            
            if isempty(concat)
                beg=1;
            else
                beg=find(concat==concat(end));
            end
            
            %             subplot(3,1,2)
            %             plot(signal1(start(trial):ending(trial)))
            
            concat=[concat (signal1(start(trial):ending(trial)))];
            
            map{iii-1,co}(trial,:)=[beg length(signal1(start(trial):ending(trial)))+beg-1];
        end
        clearvars -except map concat signal_all iii cohort norm cond norm_NS folder co
        %         subplot(3,1,3)
        %         plot(concat)
    end
    
    clearvars -except signal_all map norm cond cohort iii concat cond folder
    
    if norm==1
        signal_all{iii-1,1}=zscore(concat);
    else
        signal_all{iii-1,1}=concat;
    end
    
end

clearvars -except signal_all cond map norm

cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA')

if norm==1
    name1='zsCONC';
elseif norm==0
    name1='NONzs';
else
    name1='NSnorm';
end

filename=strcat('ACC_spiral_',num2str(name1),'.mat');
save(filename)

% % save('ACC_zscore_acr_cond.mat')


