cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')

clear all; close all
cohort = [1 3 4 6];
norm=0;

for iii=3
    concat=[];
    for cond=4
        s=struct;
        if cond==1
            load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_NS.mat'))
        elseif cond==2
            load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_RS.mat'))
        elseif cond==3
            load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_spiral/P0',num2str(cohort(iii)),'_PLS_S.mat'))
        else
            load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_hf/P0',num2str(cohort(iii)),'_HFS.mat'))
        end
        
        [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
        data=d.data_raw;  tremor_ds=d.data_ds;
        addon=92; addon_end=35;
        
        if cond==1
            NS_BE_S
            start=hu{iii,:};
            ending=hd{iii,:};
            
        elseif cond==2
            co=2;start=[];ending=[];yy=[];h_up=[];spiral=1;
            [peak_ax, start, ending, yy, h_up, ts, te]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
            start=ts{iii,co};
            ending=te{iii,co};
            
            if iii==3
                ending(end-1)=2329000;
                start(end)=[];
                ending(end)=[];  
            end
            
            %to avoid mov artifacts
            start=start+100;
            ending=ending-100;
        elseif cond==3
            new = find(data(2,:) > 4);
            difp = find((diff(new)) > 100000); % are you trying to threshold at 9.6 seconds?
            ep_1 = [new(difp) new(end)];
            sp_1 = [new(1) new(difp+1)];
            
            dum=find((ep_1-sp_1)<(60000./samplerate)*samplerateold);
            if ~isempty(dum)
                sp_1(1,dum)=NaN;
                ep_1(1,dum)=NaN;
            end
            sp=sp_1(~isnan(sp_1));
            ep=ep_1(~isnan(ep_1)); 
    
            % %                                     time=1:length(data(1,:));
            % %                                     plot(time,data(4,:))
            % %                                     hold on
            % %                                     plot(time,data(3,:))
            % %                                     plot(time(sp),data(4,sp),'r.')
            % %                                     plot(time(ep),data(4,ep),'k.')
            % %
            
            start= floor((sp./samplerateold)*samplerate)+addon; ending = floor((ep./samplerateold)*samplerate)+addon+addon_end;
            
            if iii==4
                ending=[150900];
            end
            
            
            if iii==3
                
                
                ending(1)=109000;
                start(2)=180700;
                ending(2)=256000;
                start(3)=529000;
                ending(3)=595900;
            end
        else
            HFS_BE_S; start=hu{iii,:}; ending=hd{iii,:};
        end
        
        
[f1]=spirals_cond(tremor_ds,samplerate,start,ending);
        
        
        y=d.data_ds;
        y(4,:)=sqrt((y(1,:).^2+y(2,:).^2+y(3,:).^2));
        
        if cond==1 
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
            
            concat=[concat signal1(start(trial):ending(trial))];
            
            map{iii-1,cond}(trial,:)=[beg length(signal1(start(trial):ending(trial)))+beg-1];
        end
        clearvars -except map concat signal_all iii cohort pca1 y norm cond norm_NS
%         subplot(3,1,3)
%         plot(concat)
    end
    
    clearvars -except signal_all map norm cond cohort iii concat cond
    
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


