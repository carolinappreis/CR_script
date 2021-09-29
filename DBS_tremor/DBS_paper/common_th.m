
clear; close
cohort = [1 3 4 6];
match=[[3 3 2];[3 3 2];[1 1 1];[1 1 1]];



for iii=2:4
    if iii==4
        match(4,:)=[3 3 2];
    end
    signal=[]; unf_signal=[];
    for cond=1:4
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
            [peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
            st=start{iii,co};clear start; start=st;
            en=ending{iii,co};clear ending; ending=en;
            
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
            
            if iii==4
           
                ep(2)=ep(1);
                ep(1)=953991;
                sp(2)=1022995;
            end
            %                         time=1:length(data(1,:));
            %                         plot(time,data(4,:))
            %                         hold on
            %                         plot(time(sp),data(4,sp),'r.')
            %                         plot(time(ep),data(4,ep),'k.')
            
            start= floor((sp./samplerateold)*samplerate)+addon;
            ending = floor((ep./samplerateold)*samplerate)+addon+addon_end;
            
        else
            
            HFS_BE_S
            start=hu{iii,:};
            ending=hd{iii,:};
            
        end
        
        y=zscore(d.data_ds(match(cond,iii-1),:));
        samplerate=1000;
        Fpeak=5;
        [afilt, bfilt] = butter(2, [(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
        fdata=filtfilt(afilt,bfilt,y);
        x= abs(hilbert(fdata));
    
        handup = [];
        for i = 1:length(start)
            handup = [handup start(i):ending(i)]; %#ok<*AGROW>
        end
        clear i
        
        ggg=d.data_ds(match(cond,iii-1),:);
        handup = sort(handup,'ascend');
        unf_signal=[ unf_signal ggg(handup)];
        signal= [signal x(handup)];
        clearvars -except match cohort dur_b int_b cond iii Pcurve Ppeak F sig_all signal Pcurve Ppeak unf_signal ufsignal_all
    end
    sig_all{iii-1,1}=signal;
    ufsignal_all{iii-1,1}=unf_signal;
end

for i=1:3
    subplot(1,3,i)
    data=sig_all{i,1};
    pct(i,:)=[prctile(data,25) prctile(data,50)];
    plot(data)
    hold on
    yline(pct(i,1),'r','LineWidth',2)
    yline(pct(i,2),'k','LineWidth',2)
    clear data fdata envdata 
end

clearvars -except pct ufsignal_all unf_signal    %saved as com_pct