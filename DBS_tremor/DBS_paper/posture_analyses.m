clear; close
load('Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper/comb_pct.mat')
cohort = [1 3 4 6];
match=[[NaN 3 NaN 1];[NaN 3  NaN 1];[NaN 1 NaN 1];[NaN 1 NaN 1]];


for iii=1:4
    if (iii~=1 & iii~=3)
        if iii==4
            match(4,:)=[NaN 3 NaN 1];
        end
        for cond=1:4
            s=struct;
            if cond==1
                load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_NS.mat'))
            elseif cond==2
                load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_RS.mat'))
            elseif cond==3
                if iii==2
                load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P0',num2str(cohort(iii)),'_PLS_P1.mat'))
                else
                  load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P0',num2str(cohort(iii)),'_PLS_P2.mat')) 
                end
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
                
                % %             for i=1:3
                % %                 figure(i)
                % %                 signal= d.data_ds(i,:);
                % %                 [m,n]=butter(3,[10/(0.5*d.samplerateold)],'low');
                % %                 t3(1,:)=filtfilt(m,n,signal);
                % %
                % %
                % %                 subplot(2,1,1)
                % %                 plot(t3)
                % %                 hold on
                % %                 for i=1:length(start)
                % %                     xline(start(i),'LineWidth',2)
                % %                     xline(ending(i),'LineWidth',2)
                % %                 end
                % %                 subplot(2,1,2)
                % %                 plot(d.ds_dall(i,:))
                % %
                % %             end
                
            elseif cond==2
                
                co=2;start=[];ending=[];yy=[];h_up=[];spiral=0;
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
                
                % for i=1:3
                %     figure(i)
                %     signal= d.data_ds(i,:);
                %     [m,n]=butter(3,[10/(0.5*d.samplerateold)],'low');
                %     t3(1,:)=filtfilt(m,n,signal);
                %
                %     subplot(2,1,1)
                %     plot(t3)
                %     subplot(2,1,2)
                %     plot(d.ds_dall(i,:))
                % end
            end
            
            
            
            handup = [];
            for i = 1:length(start)
                handup = [handup start(i):ending(i)]; %#ok<*AGROW>
            end
            clear i
            handup = sort(handup,'ascend');
            
            
            
            
                    signal= d.data_ds(match(cond,iii),handup);
                    [Pxx,F] = pwelch(signal, samplerate, [], round(samplerate), samplerate);
                    Pxxrange = Pxx(3:10);
                     Pcurve(iii,cond,:)=Pxx./(sum(Pxx(100:500)));
%                      Pcurve(iii,cond,:)=Pxx;
%                     Ppeak(iii,cond,:) = max(Pxxrange);
                    Ppeak(iii,cond,:)=max(squeeze(Pcurve(iii,cond,3:10)));
            
            
            
%                     aa = match(cond,iii);
%                     [Pxx,F] = pwelch(tremor_ds(aa,handup), samplerate, [], round(samplerate), samplerate);
%                     frange = F(3:10);
%                     Pxxrange = Pxx(3:10);
%                     Freqpeak= frange(find(Pxxrange == max(Pxxrange)));
%                     Ppeak= max(Pxxrange);
            %         peak_ax =[(Freqpeak(find(Ppeak == max(Ppeak)))) (find(Ppeak == max(Ppeak)))];
            peak_ax=5;
            co=1; s=struct;
            [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);
            
            
            figure(iii)
            subplot(1,4,cond)
            plot(s.env{iii,1}(match(cond,iii),:))
            hold on
            %                 yline(thr)
            if cond~=2
                for i=1:length(start)
                    xline(start(i),'LineWidth',2)
                    xline(ending(i),'LineWidth',2)
                end
            end
            box('off')
            
            
            %       psi(iii,cond,:)=[circ_r([ s.phase{iii,1}(1,handup)-s.phase{iii,1}(2,handup)]') ;circ_r([s.phase{iii,1}(1,handup)- s.phase{iii,1}(3,handup)]'); circ_r([ s.phase{iii,1}(2,handup)- s.phase{iii,1}(3,handup)]')];
            
                    avg_t_sev(iii,cond,:)=[[ mean(s.env_acc{iii,1}(match(cond,iii),handup))] [std(s.env_acc{iii,1}(match(cond,iii),handup))]];
            %        avg_t_zsev(iii,cond,:)=[[ mean(s.zenv{iii,1}(match(cond,iii),handup))] [std(s.zenv{iii,1}(match(cond,iii),handup))]];
            
            %         data_e=s.zenv{iii,1};
            %         env=data_e(match(cond,iii),handup);
            %
            %          thr=pct(iii,2);
            % %                  thr=prctile(env,50);
            %
            %         [ab_2,norm_bd,begin2]=anti_burst(env,thr);
            % %
            % %
            %           prct_below(iii,cond,:)=(sum(norm_bd)./length(env))*100;
            %         dur_b(iii,cond,:)=[mean(norm_bd) (nanmean(norm_bd)+nanstd(norm_bd))./sqrt(length(begin2))];
            %         int_b(iii,cond,:)=[mean(ab_2) (nanmean(ab_2)+nanstd(ab_2))./sqrt(length(ab_2))];
            clearvars -except match cohort dur_b int_b cond iii Pcurve Ppeak F pct prct_below avg_t_sev avg_t_zsev psi
            
        end
    end
end
end
for i=1:3
    psi_pls(i,:)=[mean(squeeze(psi(i,3,:))) std(squeeze(psi(i,3,:)))]
end


for ii=1:4
    for cc=1:4
        subplot(1,4,ii)
        plot(F(3:10),(squeeze(Pcurve(ii,cc,3:10))))
        hold on
    end
    
    for i=2:4
        change(ii,i-1,:)=round(((Ppeak(ii,i,:)-Ppeak(ii,1,:))./Ppeak(ii,1,:)*100),0);
    end
    
end
figure
bar(change)
box('off')
ylabel('% supression from baseline')
xlabel('patients')


load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
color_b1=[aegean;stone;blushred];

data=int_b;
f1=figure()
for iii=1:3
    subplot(1,3,iii)
    err=squeeze(data(iii,:,2));
    x=squeeze(data(iii,:,1));
    bar([x],'EdgeColor','none','FaceColor',color_b1(1,:),'FaceAlpha',0.5)
    hold on
    %     errorbar([1:4],x,err,'.','Color',color_b1(1,:),'LineWidth',2)
    % ylim([0 500])
    % yticks(0:100:500)
    xticklabels({'NS','RS','PLS','HFS'})
    if data==dur_b
        ylabel('burst duration')
    else
        ylabel('int burst interval')
    end
    box('off')
    set(gca,'FontSize',14)
    clear x err
end
f1.OuterPosition= [1,100,1000,300];
set(f1,'color','w');

