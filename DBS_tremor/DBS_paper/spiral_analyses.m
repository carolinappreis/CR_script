clear; close
load('Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper/comb_pct.mat')
cohort = [1 3 4 6];
match=[[3 3 2];[3 3 2];[1 1 1];[1 1 1]];


for iii=2
    if iii==4
        match(4,:)=[3 3 2];
    end
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
        
        
        
        handup = [];
        for i = 1:length(start)
            handup = [handup start(i):ending(i)]; %#ok<*AGROW>
        end
        clear i
        handup = sort(handup,'ascend');
        
        
        
        signal= d.data_ds(match(cond,iii-1),handup);
        [Pxx,F] = pwelch(signal, samplerate, [], round(samplerate), samplerate);
        Pxxrange = Pxx(3:10);
         Pcurve(iii-1,cond,:)=Pxx./(sum(Pxx(100:500)));
%          Pcurve(iii-1,cond,:)=Pxx;
%         Ppeak(iii-1,cond,:) = max(Pxxrange);
        Ppeak(iii-1,cond,:)=max(squeeze(Pcurve(iii-1,cond,3:10)));
        
        
        
        %         aa = match(cond,iii-1);
        %         [Pxx,F] = pwelch(tremor_ds(aa,handup), samplerate, [], round(samplerate), samplerate);
        %         frange = F(3:10);
        %         Pxxrange = Pxx(3:10);
        %         Freqpeak= frange(find(Pxxrange == max(Pxxrange)));
        %         Ppeak= max(Pxxrange);
        %         peak_ax =[(Freqpeak(find(Ppeak == max(Ppeak)))) (find(Ppeak == max(Ppeak)))];
        peak_ax=5;
        co=1; s=struct;
        [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);
        
%         [xx]=pha_sta(d,s)
        
                figure(iii-1)
%                 subplot(1,4,cond)
                plot(s.env{iii,1}(match(cond,iii-1),:))
                hold on
                %                 yline(thr)
                if cond~=2
                    for i=1:length(start)
                        xline(start(i),'LineWidth',2)
                        xline(ending(i),'LineWidth',2)
                    end
                end
                box('off')
        
        
%        psi(iii-1,cond,:)=[circ_r([ s.phase{iii,1}(1,handup)-s.phase{iii,1}(2,handup)]') ;circ_r([s.phase{iii,1}(1,handup)- s.phase{iii,1}(3,handup)]'); circ_r([ s.phase{iii,1}(2,handup)- s.phase{iii,1}(3,handup)]')];
%         
       avg_t_sev(iii-1,cond,:)=[[ mean(s.env_acc{iii,1}(match(cond,iii-1),handup))] [std(s.env_acc{iii,1}(match(cond,iii-1),handup))]];
%        avg_t_zsev(iii-1,cond,:)=[[ mean(s.zenv{iii,1}(match(cond,iii-1),handup))] [std(s.zenv{iii,1}(match(cond,iii-1),handup))]];
        
        data_e=s.zenv{iii,1};
        env=data_e(match(cond,iii-1),handup);
% % %         
         thr=pct(iii-1,2);
%          thr=prctile(env,50);
        
        [ab_2,norm_bd,begin2,nr]=anti_burst(env,thr);
% %         
        nr_above(iii-1,cond,:)=nr;
        prct_below(iii-1,cond,:)=(sum(norm_bd)./length(env))*100;
        dur_b(iii-1,cond,:)=[nanmean(norm_bd) (nanmean(norm_bd)+nanstd(norm_bd))./sqrt(length(begin2))];
        int_b(iii-1,cond,:)=[nanmean(ab_2) (nanmean(ab_2)+nanstd(ab_2))./sqrt(length(ab_2))];
        clearvars -except match cohort dur_b int_b cond iii Pcurve Ppeak F pct prct_below avg_t_sev avg_t_zsev psi nr_above
        
    end
end

%%%saved as all_cond_spiral