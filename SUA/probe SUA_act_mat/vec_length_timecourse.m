
clear all
 cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
%  cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load ('data_SUA_BZ.mat')

BZ.sua=data_region;
BZ.ctx=Ecog_region;

srn=1000;

for i =1:size(BZ.sua,1)
    [b,a]=butter(2,[20/(0.5*srn) 30/(0.5*srn)],'bandpass');
    SBZ.ecog_filt(i,:)=filtfilt(b,a,BZ.ctx(i,:));
    SBZ.ecog_env(i,:)=abs(hilbert(SBZ.ecog_filt(i,:)));
    
    env=SBZ.ecog_env(i,:);Ecogfiltered=SBZ.ecog_filt(i,:);
    
    SBZ.onset_raw{i,1}=bursts(env);
    SBZ.offset_raw{i,1}=bursts_off(env);
    SBZ.onset_phase_al{i,1}=bursts_aligned(env,Ecogfiltered);
    SBZ.offset_phase_al{i,1}=bursts_aligned_off(env,Ecogfiltered);
    SBZ.ecog_phase(i,:)=wrapTo2Pi(angle(hilbert(BZ.ctx(i,:))));
    cy_bursts{i,1}=cycles_10(env,Ecogfiltered);
    clearvars env Ecogfiltered
end


 n=1;
for u=1:size(data_region,1)
    units_match2=data_region{u,1};
    close all
    block=[];
    block = cy_bursts{u,1}{2,1}(any(cy_bursts{u,1}{2,1},2),:);
    for um=1:size(units_match2,1)
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