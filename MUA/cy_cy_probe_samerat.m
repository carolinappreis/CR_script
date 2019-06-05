clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\MUA')
load ('SNR_cycle_probe.mat')
same_a=[28;29];
WaveData_DC=ctx(same_a,:);

hp=find((cellfun(@isempty,SUA))==0);
r=[];
for i=1:length(same_a)
r=[r find(hp(:,1)==same_a(i))];
end


for i=1:length(r)
stat_d1{i,1}=stat_d{r(i),1};
end

for i =1:size(stat_d,1)
    %     SSNR.psd(i,:)=pwelch(WaveData_DC(stat_d(i),:),1000,[],1000,1000);
    SSNR.ecog_filt(i,:)=ecogbf_match(i,:);
    SSNR.ecog_env(i,:)=abs(hilbert(ecogbf_match(i,:)));
    env=SSNR.ecog_env(i,:);Ecogfiltered=SSNR.ecog_filt(i,:);
    SSNR.onset_raw{i,1}=bursts(env);
    SSNR.offset_raw{i,1}=bursts_off(env);
    SSNR.onset_phase_al{i,1}=bursts_aligned(env,Ecogfiltered);
    SSNR.offset_phase_al{i,1}=bursts_aligned_off(env,Ecogfiltered);
    SSNR.ecog_phase(i,:)=wrapTo2Pi(angle(hilbert(ecogbf_match(i,:))));
    cy_bursts{i,1}=cycles_10(env,Ecogfiltered);
    clearvars env Ecogfiltered
end


for u=r(1):r(end)
    units_match1=units_match{u,1};
%     close all
    for um=1:size(units_match1,1)
        block = cy_bursts{u,1}{2,1}(any(cy_bursts{u,1}{2,1},2),:);
        for d1=1:size(block,1)
            for d2=1:size(block,2)
                if d2+1<length(block(d1,:))
                    epoch=block(d1,d2):block(d1,d2+1);
                    l=find(units_match1(um,epoch)==1);
                    pha_b{um,1}{d1,d2}=SSNR.ecog_phase(u,epoch(l));
                    pha_b_all{u,um}{d1,d2}=SSNR.ecog_phase(u,epoch(l));
                    idx_spkcycle{um,1}{d1,d2}=epoch(l);
                end
            end
        end
    end
    for cell=1:size(pha_b,1)
        for ii =1:size(pha_b{cell,1},2)
            for i=1:size(pha_b{cell,1},1)
                if ~isempty (pha_b{cell,1}{i,ii})
                    bu{i,1}=pha_b{cell,1}{i,ii}(1);
                end
            end
            bu=bu(~cellfun('isempty',bu));
            vec_lg(cell,ii)=circ_r(cell2mat(bu));
            pref_pha(cell,ii)=circ_mean(cell2mat(bu));
            %         zm(cell,ii) = circ_r(cell2mat(bu)).*(exp(sqrt(-1).*(circ_mean(cell2mat(bu)))));
            clear bu
        end
        cyc_avg(1,:)=mean(vec_lg,1);
        cyc_ang_avg(1,:)=circ_mean(pref_pha);
        
        clearvars -except zm cell pha_b_all pha_b vec_lg  cyc_avg cyc_ang_avg pref_pha units_match cy_bursts SSNR u f units_match1
%     for i=11:15
%         figure(cell+1)
%         if i==11
%             p1=polarplot([0 pref_pha(cell,i)], [0, vec_lg(cell,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%             hold on
%         elseif i==12
%             p2=polarplot([pref_pha(cell,i-1) pref_pha(cell,i)], [vec_lg(cell,i-1), vec_lg(cell,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%         elseif i==13
%             p3=polarplot([pref_pha(cell,i-1) pref_pha(cell,i)], [vec_lg(cell,i-1), vec_lg(cell,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%         elseif i==14
%             p4=polarplot([pref_pha(cell,i-1) pref_pha(cell,i)], [vec_lg(cell,i-1), vec_lg(cell,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%         elseif i==15
%             p5=polarplot([pref_pha(cell,i-1) pref_pha(cell,i)], [vec_lg(cell,i-1), vec_lg(cell,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%             
%         end
%     end
%     rlim([0 0.7])
%     legend([p1 p2 p3 p4 p5],{'1st cycle','2nd cycle','3rd cycle','4th cycle','5th cycle'},'Box','off','Orientation','horizontal','Location','south','FontSize',9)
%     
    end
        %per animal
        for i=11:15
            if i==11
                p1=polarplot([0 cyc_ang_avg(i)], [0, cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
                hold on
            elseif i==12
                p2=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
            elseif i==13
                p3=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
            elseif i==14
                p4=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
            elseif i==15
                p5=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
    
            end
            rlim([0 0.6])
        end
        legend([p1 p2 p3 p4 p5],{'1st cycle','2nd cycle','3rd cycle','4th cycle','5th cycle'},'Box','off','Orientation','horizontal','Location','south','FontSize',9)
%    close all
    clear pha_b
end




% circ_plot(cell2mat(bu),'pretty','bo',true,'linewidth',2,'color','r')

%checking plots
% plot(log(SSNR.psd(:,1:50))')
% env=SSNR.ecog_env(1,:);Ecogfiltered=SSNR.ecog_filt(1,:);
%
% plot(time,env)
% hold on
% plot(time,Ecogfiltered)
% plot(time(SSNR.onset_phase_al{1,1}{1,1}),Ecogfiltered(SSNR.onset_phase_al{1,1}{1,1}),'bo')
% plot(time(SSNR.offset_phase_al{1,1}{1,1}),Ecogfiltered(SSNR.offset_phase_al{1,1}{1,1}),'ko')
% plot(time(SSNR.onset_raw{1,1}{1,1}),env(SSNR.onset_raw{1,1}{1,1}),'b.')
% plot(time(SSNR.offset_raw{1,1}{1,1}),env(SSNR.offset_raw{1,1}{1,1}),'k.')
%
% plot(time,Ecogfiltered)
% hold on
% plot(time(cy_bursts{1,1}{2,1}(2,:)),Ecogfiltered(cy_bursts{1,1}{2,1}(2,:)),'r.')
% plot(time(cell2mat(idx_spkcycle{1,1}(1,:))),Ecogfiltered(cell2mat(idx_spkcycle{1,1}(1,:))),'k.')