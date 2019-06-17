
clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load ('SNR_cycle_sua.mat')


for i =1:size(SNR.sua,1)
    [b,a]=butter(2,[20/(0.5*srn) 30/(0.5*srn)],'bandpass');
    SSNR.ecog_filt(i,:)=filtfilt(b,a,SNR.ctx(i,:));
    SSNR.ecog_env(i,:)=abs(hilbert(SSNR.ecog_filt(i,:)));
    
    env=SSNR.ecog_env(i,:);Ecogfiltered=SSNR.ecog_filt(i,:);
    
    SSNR.onset_raw{i,1}=bursts(env);
    SSNR.offset_raw{i,1}=bursts_off(env);
    SSNR.onset_phase_al{i,1}=bursts_aligned(env,Ecogfiltered);
    SSNR.offset_phase_al{i,1}=bursts_aligned_off(env,Ecogfiltered);
    SSNR.ecog_phase(i,:)=wrapTo2Pi(angle(hilbert(SNR.ctx(i,:))));
    cy_bursts{i,1}=cycles_10(env,Ecogfiltered);
    clearvars env Ecogfiltered
end


 n=1;
for u=1:size(data_all,1)
    units_match2=data_all{u,1};
    close all
    block=[];
    block = cy_bursts{u,1}{2,1}(any(cy_bursts{u,1}{2,1},2),:);
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
                if ~isempty (pha_b{ctc,1}{i,ii})  && numel(find(pha_b_l{ctc,1}(:,ii)>0))>20
                    bu{i,1}=pha_b{ctc,1}{i,ii}(1);
                else
                    bu{i,1}=[];
                    %         zm(ctc,ii) = circ_r(ctc2mat(bu)).*(exp(sqrt(-1).*(circ_mean(ctc2mat(bu)))));
                end
            end
            if ~isempty (cell2mat(bu))
                bu1=cell2mat(bu);
            else bu1=NaN;
            end
            vec_lg(n,ii)=circ_r(bu1);
            pref_pha(u,ctc,ii)=circ_mean(bu1);
            clear bu  
        end
        n=n+1;
    end
    clear pha_b
end

plot(nanmean(vec_lg,1),'-d')
xlim([0 22])

figure()
err=nanstd(vec_lg);
errorbar(nanmean(vec_lg),err)
xlim([0 22])