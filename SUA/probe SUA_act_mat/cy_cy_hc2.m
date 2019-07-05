

clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load('SNR_sua_skrate.mat') ; 
subj= SNR.beta_rats;

Ecog=SNR.ctx(subj,:);
for i =1:size(subj,1)
for ii=1:size(SNR.beta_idx{subj(i),1},2)
    hp=SNR.beta_idx{subj(i),1}(1,ii);
    data_all{i,1}(ii,:)=SNR.sua{subj(i),1}(hp,:);
    clear hp
end
end

srn=1000;

for i =1:size(Ecog,1)
    %     SSNR.psd(i,:)=pwelch(WaveData_DC(stat_d(i),:),1000,[],1000,1000);
      [b,a]=butter(2,[15/(0.5*srn) 35/(0.5*srn)],'bandpass');
    Ecogfiltered=filtfilt(b,a,Ecog(i,:));
    SSNR.ecog_filt(i,:)=Ecogfiltered;
    SSNR.ecog_env(i,:)=abs(hilbert(Ecogfiltered));
    env=SSNR.ecog_env(i,:);Ecogfiltered=SSNR.ecog_filt(i,:);
    SSNR.onset_raw{i,1}=bursts(env);
    SSNR.offset_raw{i,1}=bursts_off(env);
    SSNR.onset_phase_al{i,1}=bursts_aligned(env,Ecogfiltered);
    SSNR.offset_phase_al{i,1}=bursts_aligned_off(env,Ecogfiltered);
    SSNR.ecog_phase(i,:)=wrapTo2Pi(angle(hilbert(Ecog(i,:))));
    cy_bursts{i,1}=cycles_10(env,Ecogfiltered); 
% cy_bursts{i,1}=cell2mat(cy_bursts{i,1});
    clearvars env Ecogfiltered
end

go=1;
for u=1:size(data_all,1)
    units_match2=data_all{u,1};
    close all
    block=[];
      block = cy_bursts{u,1}{2,1}(any(cy_bursts{u,1}{2,1},2),:);
%  block = cy_bursts{u,1}(any(cy_bursts{u,1},2),:);
    for um=1:size(units_match2,1)
        for d1=1:size(block,1)
            for d2=1:size(block,2)
                if d2+1<length(block(d1,:))
                    epoch=block(d1,d2):block(d1,d2+1);
                    l=find(units_match2(um,epoch)==1);
                    pha_b{um,1}{d1,d2}=SSNR.ecog_phase(u,epoch(l));
                    pha_b_l{um,1}(d1,d2)=length(SSNR.ecog_phase(u,epoch(l)));
                    pha_b_all{u,um}{d1,d2}=SSNR.ecog_phase(u,epoch(l));
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
                bu1=cell2mat(bu);
                bu1=bu1(randperm(20));
            else
                bu1=NaN;
            end
            vec_lg(go,ii)=circ_r(bu1);
            pref_pha(go,ii)=circ_mean(bu1);
            uni_pha(go,ii)=circ_rtest(bu1)
            clear bu  
        end
        clearvars -except uni_pha go zm idx_spkcycle pha_b_all pha_b pha_b_l vec_lg ctc cyc_avg cyc_ang_avg pref_pha units_match cy_bursts SSNR u f data_all Ecog
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
fig=figure()
cy_plots=7:15;
for i=1:length(cy_plots)
    subplot(1,length(cy_plots),i,polaraxes)
    for un=1:size(vl,1)
        polarplot([0 pp(un,cy_plots(i))], [0, vl(un,cy_plots(i))],'linewidth',2)
        rlim([0 0.9])
        hold on
    end
end

fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 40, 10];
fig.Color='w';

% fig=figure()
% cy_plots=6:14;
% for i=1:length(cy_plots)
%     subplot(1,length(cy_plots),i,polaraxes)
%     for un=1:size(vec_lg,1)
%         polarplot([0 pref_pha(un,cy_plots(i))], [0, vec_lg(un,cy_plots(i))],'linewidth',2)
%         hold on
%     end
% end
% fig.Units = 'centimeters';
% fig.OuterPosition= [10, 10, 40, 10];
% fig.Color='w';
% rlim([0 0.7])
% legend([p1 p2 p3 p4 p5],{'1st cycle','2nd cycle','3rd cycle','4th cycle','5th cycle'},'Box','off','Orientation','horizontal','Location','south','FontSize',9)


% vec_m=nanmean(vec_lg,1);
% ang_m=circ_mean(pref_pha);
% fig=figure()
% cy_plots=6:14;
% for i=1:length(cy_plots)
%     subplot(1,length(cy_plots),i,polaraxes)
%     for un=1:size(vec_lg,1)
%         polarplot([0 ang_m(1,cy_plots(i))], [0, vec_m(1,cy_plots(i))],'linewidth',2)
%         hold on
%     end
% end
% fig.Units = 'centimeters';
% fig.OuterPosition= [10, 10, 40, 10];
% fig.Color='w';