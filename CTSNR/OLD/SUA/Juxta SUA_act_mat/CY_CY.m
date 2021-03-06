clear all
 cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
load ('SUA_BZ')
load('BZ_cycle')

for i =1:size(stat_d,1)
    UBZ.psd(i,:)=pwelch(WaveData_DC(stat_d(i),:),1000,[],1000,1000);
    UBZ.ecog_filt(i,:)=ecogbf_match(i,:);
    UBZ.ecog_env(i,:)=abs(hilbert(ecogbf_match(i,:)));
    env=UBZ.ecog_env(i,:);Ecogfiltered=UBZ.ecog_filt(i,:);
    UBZ.onset_raw{i,1}=bursts(env);
    UBZ.offset_raw{i,1}=bursts_off(env);
    UBZ.onset_phase_al{i,1}=bursts_aligned(env,Ecogfiltered);
    UBZ.offset_phase_al{i,1}=bursts_aligned_off(env,Ecogfiltered);
    UBZ.ecog_phase(i,:)=wrapTo2Pi(angle(hilbert(ecogbf_match(i,:))));
    cy_bursts{i,1}=cycles_10(env,Ecogfiltered);
    clearvars env Ecogfiltered
end

for f=1:size(cy_bursts,1)
    block = cy_bursts{f,1}{2,1}(any(cy_bursts{f,1}{2,1},2),:);
    for d1=1:size(block,1)
        for d2=1:size(block,2)
            if d2+1<length(block(d1,:))
                epoch=block(d1,d2):block(d1,d2+1);
                l=find(units_match(f,epoch)==1);
                pha_b{f,1}{d1,d2}=UBZ.ecog_phase(f,epoch(l));
                pha_bl{f,1}(d1,d2)=length(UBZ.ecog_phase(f,epoch(l)));
                idx_spkcycle{f,1}{d1,d2}=epoch(l);
            end
        end
    end
end

% %checking plots
% plot(log(UBZ.psd(:,1:50))')
% env=UBZ.ecog_env(1,:);Ecogfiltered=UBZ.ecog_filt(1,:);
%
% plot(time,env)
% hold on
% plot(time,Ecogfiltered)
% plot(time(UBZ.onset_phase_al{1,1}{1,1}),Ecogfiltered(UBZ.onset_phase_al{1,1}{1,1}),'bo')
% plot(time(UBZ.offset_phase_al{1,1}{1,1}),Ecogfiltered(UBZ.offset_phase_al{1,1}{1,1}),'ko')
% plot(time(UBZ.onset_raw{1,1}{1,1}),env(UBZ.onset_raw{1,1}{1,1}),'b.')
% plot(time(UBZ.offset_raw{1,1}{1,1}),env(UBZ.offset_raw{1,1}{1,1}),'k.')
%
% plot(time,Ecogfiltered)
% hold on
% plot(time(cy_bursts{1,1}{2,1}(2,:)),Ecogfiltered(cy_bursts{1,1}{2,1}(2,:)),'r.')
% plot(time(cell2mat(idx_spkcycle{1,1}(1,:))),Ecogfiltered(cell2mat(idx_spkcycle{1,1}(1,:))),'k.')


for un=1:size(pha_b,1)
    for ii =1:size(pha_b{un,1},2)
        for i=1:size(pha_b{un,1},1)
            if ~isempty (pha_b{un,1}{i,ii})&& numel(find(pha_bl{un,1}(:,ii)>0))>20
                bu{i,1}=pha_b{un,1}{i,ii}(1);
            else
                bu{i,1}=[];
            end
        end
        
        if ~isempty (cell2mat(bu))
            bu1=cell2mat(bu);
            bu1=bu1(randperm(20));
        else
            bu1=NaN;
        end
 
        vec_lg(un,ii)=circ_r(bu1);
        pref_pha(un,ii)=circ_mean(bu1);
        %         zm(cell,ii) = circ_r(cell2mat(bu)).*(exp(sqrt(-1).*(circ_mean(cell2mat(bu)))));
        clear bu
    end    
    clearvars -except zm cell pha_b vec_lg  pref_pha pha_bl
end


    
fig=figure()
cy_plots=6:14;
for i=1:length(cy_plots)
    subplot(1,length(cy_plots),i,polaraxes)
    for un=1:size(vec_lg,1)
        polarplot([0 pref_pha(un,cy_plots(i))], [0, vec_lg(un,cy_plots(i))],'linewidth',2)
        hold on
    end
end
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 40, 10];
fig.Color='w';


% for n=1:size(pref_pha,2)
%     cyc_vavg(1,n)=mean(vec_lg(~isnan(vec_lg(n,:))));
%     cyc_aavg(1,n)=circ_mean(pref_pha(~isnan(pref_pha(n,:))));
% end
% fig=figure()
% cy_plots=6:14;
% for i=1:length(cy_plots)
%     subplot(1,length(cy_plots),i,polaraxes)
%     for un=1:size(vec_lg,1)
%         polarplot([0 cyc_aavg(1,cy_plots(i))], [0, cyc_vavg(1,cy_plots(i))],'linewidth',2)
%         hold on
%     end
% end
% fig.Units = 'centimeters';
% fig.OuterPosition= [10, 10, 40, 10];
% fig.Color='w';
% for n=1:size(pha_b,1)
%     for i=11:15
%         %         subplot(1,size(pha_b,1),n)
%         figure(n+1)
%         if i==11
%             p1=polarplot([0 pref_pha(n,i)], [0, vec_lg(n,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%             hold on
%         elseif i==12
%             p2=polarplot([pref_pha(n,i-1) pref_pha(n,i)], [vec_lg(n,i-1), vec_lg(n,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%         elseif i==13
%             p3=polarplot([pref_pha(n,i-1) pref_pha(n,i)], [vec_lg(n,i-1), vec_lg(n,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%         elseif i==14
%             p4=polarplot([pref_pha(n,i-1) pref_pha(n,i)], [vec_lg(n,i-1), vec_lg(n,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%         elseif i==15
%             p5=polarplot([pref_pha(n,i-1) pref_pha(n,i)], [vec_lg(n,i-1), vec_lg(n,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%             
%         end
%     end
%     rlim([0 0.7])
%     legend([p1 p2 p3 p4 p5],{'1st cycle','2nd cycle','3rd cycle','4th cycle','5th cycle'},'Box','off','Orientation','horizontal','Location','south','FontSize',9)
% end


% % figure(1)
% % for i=11:15
% %
% %     if i==11
% %         p1=polarplot([0 cyc_ang_avg(i)], [0, cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
% %           hold on
% %     elseif i==12
% %         p2=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
% %     elseif i==13
% %         p3=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
% %     elseif i==14
% %         p4=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
% %     elseif i==15
% %         p5=polarplot([cyc_ang_avg(i-1) cyc_ang_avg(i)], [cyc_avg(i-1) cyc_avg(i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
% %
% %     end
% %     rlim([0 0.6])
% % end
% %     legend([p1 p2 p3 p4 p5],{'1st cycle','2nd cycle','3rd cycle','4th cycle','5th cycle'},'Box','off','Orientation','horizontal','Location','south','FontSize',9)
% %
% %

