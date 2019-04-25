clear all
%cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')

load('SNR.mat');
% for ii=1:length(SNR.idrat)
% data{ii,1}=SNR.filt_thal{SNR.idrat(ii),1}
% end
% data=vertcat(data{:});
bins=[50:20:290];


for ik=1:size(SNR.env_ctx,1)
    clearvars dur
    ref1=SNR.onset_raw_all{ik,1};
    ref1_1=SNR.onset_pa_all{ik,1};
    ref2=SNR.offset_raw_all{ik,1};
    if length(ref1) ~= length(ref1_1)
        ref1=ref1(1:length(ref1_1));
        ref2=ref2(1:length(ref1_1));
    end
    %         [dur,dur_idx]=sort(ref2-ref1,'ascend');
    [dur,dur_idx]=sort(ref2-ref1,'ascend');
    %     dur_all(ik,1)=max(dur);
    if max(dur)>=290
        for t=1:length(bins)-1
            for i=1:length(dur)-1
                if (dur(i))>bins(t) && (dur(i))<=bins(t+1)
                    ind_b1{ik,1}(t)=ref1_1(dur_idx(i));
                    ind_d1{ik,1}(t)=dur(i);
                end
            end
        end
    end
end

for i=1:size(ind_b1,1)
    if size(ind_b1{i,1},2)==size(ind_b1{1,1},2)
        b_idx_rat(i,:)=squeeze(ind_b1{i,1});
        b_d_rat(i,:)=squeeze(ind_d1{i,1});
    end
end

ff=1;
for i=1:size(b_idx_rat,1)
    if b_idx_rat(i,1)~=0
        new_idx(1,ff)=i;
        b_pt(ff,:)=b_idx_rat(i,:);
        d_pt(ff,:)=b_d_rat(i,:);
        ff=ff+1;
    end
end
b_pt1=b_pt;
b_pt=b_pt1(:,[1:6 8]);

SNR.idrat=[new_idx];
for ik=1:length(SNR.idrat)
    ref3=b_pt(ik,:);
    for ct=1:size(SNR.phase_thal{SNR.idrat(ik),1},1)
        clearvars -except ik ct SNR epochs_zd1 epochs_zd dur_all ref3 b_pt
        
        non_norm=unwrap(SNR.phase_ctx(SNR.idrat(ik),:))-unwrap(SNR.phase_thal{SNR.idrat(ik),1}(ct,:)); %circdist
        non_norm1=diff(non_norm);
        znon_norm=zscore(non_norm1);
        el=400;
        for ii=1:length(ref3)
            if ref3(ii)>el
                epochs_z(ii,:)=znon_norm(ref3(ii)-el:ref3(ii)+el);
            end
        end
        for ii=1:size(epochs_z,1)
            for ff=1:length(epochs_z(ii,:))
                if epochs_z(ii,ff)>=1.96 | epochs_z(ii,ff)<=-1.96
                    epochs_z1(ii,ff)=1;
                else
                    epochs_z1(ii,ff)=0;
                end
            end
        end
        
%                                 plot(epochs_z1(ii,:),'r.')
%                                 imagesc(epochs_z1)
%                                 xticks([200:200:800])
%                                 xlim([200 800])
%                                 xticklabels ({'-200','0','200','400'})
        epochs_zd(ct,:,:)=epochs_z1;
    end
    epochs_zd1(ik,:,:)=squeeze(mean(epochs_zd,1));
end

% for i=1:length(SNR.idrat)
%  imagesc(squeeze(epochs_zd1(i,:,:)))
%  close all
% end

slip_b=squeeze(mean(epochs_zd1,1));

imagesc(slip_b)
xlabel ('Time(msec)')
ylabel('Bursts # (Sorted by length)')
xticks([200:200:800])
xlim([200 800])
xticklabels ({'-200','0','200','400'})
title('SNR')
%
%
% tempo=1:size(SNR.env_ctx,2);
% plot(tempo,SNR.env_ctx(SNR.idrat(ik),:))
% hold on
% plot(tempo,SNR.filt_ctx(SNR.idrat(ik),:))
% plot(tempo(ref3),SNR.env_ctx(SNR.idrat(ik),ref3),'b.')
% plot(tempo(ref1),SNR.env_ctx(SNR.idrat(ik),ref1),'ro')

