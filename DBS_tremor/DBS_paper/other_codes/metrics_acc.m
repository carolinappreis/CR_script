

clear all; close all
cohort = [1 3 4 6];
pca1=1;


for iii=4
    signal=[]; unf_signal=[];
    for cond=3
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
            
            start= floor((sp./samplerateold)*samplerate)+addon; ending = floor((ep./samplerateold)*samplerate)+addon+addon_end;
            
        else 
            HFS_BE_S; start=hu{iii,:}; ending=hd{iii,:};  
        end
        
        

        y=d.data_ds;
        y(4,:)=sqrt((y(1,:).^2+y(2,:).^2+y(3,:).^2));
        
        concat
        
    end
end

