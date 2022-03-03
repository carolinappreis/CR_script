clear; close
load('Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper/comb_pct.mat')
cohort = [1 3 4 6];
match=[[NaN 3 NaN 1];[NaN 3  NaN 1];[NaN 1 NaN 1];[NaN 1 NaN 1]];


for iii=1:4
    if (iii~=1 & iii~=3)
        if iii==4
            match(4,:)=[NaN 3 NaN 1];
            signal=[];
            for cond=3
%                 1:4
                s=struct;
                if cond==1
                    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_NS.mat'))
                elseif cond==2
                    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_RS.mat'))
                elseif cond==3
                    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_spiral/P0',num2str(cohort(iii)),'_PLS_P.mat'))
                else
                    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_hf/P0',num2str(cohort(iii)),'_HFS.mat'))
                end
                
                [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
                data=d.data_raw;  tremor_ds=d.data_ds;
                addon=92; addon_end=35;
                
                if cond==1
                    NS_BE_P
                    start=hu{iii,:};
                    ending=hd{iii,:};
                    
                elseif cond==2
                    co=2;start=[];ending=[];yy=[];h_up=[];spiral=0;
                    [peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
                    st=start{iii,co};clear start; start=st;
                    en=ending{iii,co};clear ending; ending=en;
                    
                elseif cond==3
             
                    [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
                    
                    co=3; %%% condifition 3 - phase locked
                    spiral=0;
                    out=struct; start=cell(1,3); ending=cell(1,3); yy=cell(1,3); h_up=cell(1,3); s=struct;
                    
                    [peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
                     st=start{iii,co};clear start; start=st;
                    en=ending{iii,co};clear ending; ending=en;

                else
                    
                    HFS_BE_P
                    start=hu{iii,:};
                    ending=hd{iii,:};
                    
                end
                
                y=zscore(d.data_ds(match(cond,iii),:));
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
                handup = sort(handup,'ascend');
                
                signal= [signal x(handup)];
                
                clearvars -except match cohort dur_b int_b cond iii Pcurve Ppeak F sig_all signal
            end
            sig_all{iii,1}=signal;
        end
    end
end





for i=1:4
    subplot(1,3,i)
    data=sig_all{i,1};
    pct(i,:)=[prctile(data,25) prctile(data,50)];
    plot(data)
    hold on
    yline(pct(i,1),'r','LineWidth',2)
    yline(pct(i,2),'k','LineWidth',2)
    clear data fdata envdata
end

clearvars -except pct