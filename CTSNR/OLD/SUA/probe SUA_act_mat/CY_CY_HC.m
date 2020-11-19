
clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load ('NEW_SNR_cycle.mat')

stat_d1=stat_d(~cellfun('isempty',stat_d));
units_match1=units_match(~cellfun('isempty',units_match));
ecogbf_match1=ecogbf_match(any(ecogbf_match,2),:);




for i =1:size(units_match1,1)
    %     SBZ.psd(i,:)=pwelch(WaveData_DC(stat_d(i),:),1000,[],1000,1000);
    SBZ.ecog_filt(i,:)=ecogbf_match1(i,:);
    SBZ.ecog_env(i,:)=abs(hilbert(ecogbf_match1(i,:)));
    env=SBZ.ecog_env(i,:);Ecogfiltered=SBZ.ecog_filt(i,:);
    SBZ.onset_raw{i,1}=bursts(env);
    SBZ.offset_raw{i,1}=bursts_off(env);
    SBZ.onset_phase_al{i,1}=bursts_aligned(env,Ecogfiltered);
    SBZ.offset_phase_al{i,1}=bursts_aligned_off(env,Ecogfiltered);
    SBZ.ecog_phase(i,:)=wrapToPi(angle(hilbert(ecogbf_match1(i,:))));
    cy_bursts{i,1}=cycles_10(env,Ecogfiltered);
    cy_bursts{i,1}=cell2mat(cy_bursts{i,1});
    clearvars env Ecogfiltered
end

go=1;
phase_units=[];
for u=1:size(units_match1,1)
    units_match2=units_match1{u,1};
    close all
    block=[];
    %     block = cy_bursts{u,1}{2,1}(any(cy_bursts{u,1}{2,1},2),:);
    block = cy_bursts{u,1}(any(cy_bursts{u,1},2),:);
    for um=1:size(units_match2,1)
            phase_units=[ phase_units  SBZ.ecog_phase(u,find(units_match2(um,:)==1))];
        for d1=1:size(block,1)
            for d2=1:size(block,2)
                if d2+1<length(block(d1,:))
                    epoch=block(d1,d2):block(d1,d2+1);
                    l=find(units_match2(um,epoch)==1);
                    pha_b{um,1}{d1,d2}=SBZ.ecog_phase(u,epoch(l));
                    pha_b_l{um,1}(d1,d2)=length(SBZ.ecog_phase(u,epoch(l)));
                    pha_b_all{u,um}{d1,d2}=SBZ.ecog_phase(u,epoch(l));
                    %                     idx_spkcycle{u,um}{d1,d2}=epoch(l);
                    if isempty (epoch(l))
                        idx_spkcycle{u,um}(d1,d2)=0;
                    else
                        idx_spkcycle{u,um}(d1,d2)=1;
                    end
                end
            end
        end
    end
    
    for ctc=1:size(pha_b,1)
        for ii =1:size(pha_b{ctc,1},2)
            for i=1:size(pha_b{ctc,1},1)
                if ~isempty (pha_b{ctc,1}{i,ii})
                    bu{i,1}=pha_b{ctc,1}{i,ii}(1);
                else
                    bu{i,1}=[];
                    %         zm(ctc,ii) = circ_r(ctc2mat(bu)).*(exp(sqrt(-1).*(circ_mean(ctc2mat(bu)))));
                end
            end
            
            bubu=cell2mat(bu);
            if ~isempty (bubu)&& length(bubu)>=20
                bu1=bubu(randperm(20));
            else
                bu1=NaN;
            end
            vec_lg(go,ii)=circ_r(bu1);
            pref_pha(go,ii)=circ_mean(bu1);
            uni_pha(go,ii)=circ_rtest(bu1);
            clear bu
        end
        clearvars -except uni_pha go zm idx_spkcycle pha_b_all pha_b pha_b_l vec_lg phase_units ctc cyc_avg cyc_ang_avg pref_pha units_match cy_bursts SBZ u f units_match1
        go=go+1;
    end
    clear pha_b
end

%  size(cell2mat(units_match1),1)==size(vec_lg,1)
n=[];
r=1;
for i=1:size(vec_lg,1)
    if ~isnan(sum(vec_lg(i,:)))
        vl(r,:)=vec_lg(i,:);
        pp(r,:)=pref_pha(i,:);
        r=r+1;
    end
end


%%% cycle eleven is burst onset

fig=figure()
cy_plots=7:15;
for i=1:length(cy_plots)
    subplot(1,length(cy_plots),i,polaraxes)
    for un=1:size(vl,1)
        polarplot([0 pp(un,cy_plots(i))], [0, vl(un,cy_plots(i))],'linewidth',2)
        rlim([0 0.5])
        hold on
    end
end

fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 40, 10];
fig.Color='w';



vec_m=nanmean(vl,1);
ang_m=circ_mean(pp);
fig=figure()
cy_plots=6:16; % first 5
for i=1:length(cy_plots)
    if ~isempty(intersect(cy_plots(i), [1:10]))
        cl=[0.5 0.5 0.5];
    elseif cy_plots(i)==11
        cl=[0.5 0 0];
    else
        cl=[0 0 0.5];
    end
    if i==1
        polarplot([0 ang_m(1,cy_plots(i))], [0, vec_m(1,cy_plots(i))],'Color',cl,'linewidth',1,'MarkerIndices',[3],'Marker','d')
    else
        polarplot([ang_m(1,cy_plots(i-1)) ang_m(1,cy_plots(i))], [vec_m(1,cy_plots(i-1)), vec_m(1,cy_plots(i))],'Color',cl,'linewidth',1,'MarkerIndices',[2],'Marker','d')
    end
    hold on
end

color_b=[0.5 0 0];
fig=figure()
% plot(nanmean(vec_lg,1),'-d','LineWidth',1.5,'Color',color_b)
plot(mean(vl),'-d','LineWidth',1.5,'Color',color_b)
hold on
xline(11,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',[0.5 0.5 0.5],'LineWidth',2)
xlim([0 22])
xlabel('Number of {\beta} cycles')
xticks([1:2:21])
xticklabels ({'-10','-8','-6','-4','-2','0','2','4','6','8','10'})
ylabel('Vector length')
box('off')
set(gca,'FontSize',12)
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 10, 10];
fig.Color='w';

% % 
% % color_b=[0.5 0.5 0.5];
% % 
% % vm=circ_mean(phase_units');
% % [vstd,vs]=circ_std(phase_units');
% % fig=figure()
% % polarplot([0 vm-vstd], [0 1])
% % hold on
% % polarplot([0 vm+vstd], [0 1])
% % polarplot([0 vm], [0 1])
% % 
% % 
% % polarhistogram((phase_units'),'FaceColor',[0.5 0.5 0],'FaceAlpha',.8,'EdgeColor','none','BinWidth',0.5)



% fig.Units = 'centimeters';
% fig.OuterPosition= [10, 10, 40, 10];
% fig.Color='w';


% 
% for i=11:15
%     figure(ctc+1)
%     if i==11
%         p1=polarplot([0 pref_pha(ctc,i)], [0, vec_lg(ctc,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%         hold on
%     elseif i==12
%         p2=polarplot([pref_pha(ctc,i-1) pref_pha(ctc,i)], [vec_lg(ctc,i-1), vec_lg(ctc,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%     elseif i==13
%         p3=polarplot([pref_pha(ctc,i-1) pref_pha(ctc,i)], [vec_lg(ctc,i-1), vec_lg(ctc,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%     elseif i==14
%         p4=polarplot([pref_pha(ctc,i-1) pref_pha(ctc,i)], [vec_lg(ctc,i-1), vec_lg(ctc,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%     elseif i==15
%         p5=polarplot([pref_pha(ctc,i-1) pref_pha(ctc,i)], [vec_lg(ctc,i-1), vec_lg(ctc,i)],'linewidth',2,'MarkerIndices',[2],'Marker','d')
%         
%     end
% end
% rlim([0 0.7])
% legend([p1 p2 p3 p4 p5],{'1st cycle','2nd cycle','3rd cycle','4th cycle','5th cycle'},'Box','off','Orientation','horizontal','Location','south','FontSize',9)
% 
% 
