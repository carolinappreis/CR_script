clear all
%cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')

load('BZ.mat');
% for ii=1:length(BZ.idrat)
% data{ii,1}=BZ.filt_thal{BZ.idrat(ii),1}
% end
% data=vertcat(data{:});
bins=[100;200];


for ik=1:size(BZ.env_ctx,1)
    clearvars dur
    ref1=BZ.onset_raw_all{ik,1};
    ref1_1=BZ.onset_pa_all{ik,1};
    ref2=BZ.offset_raw_all{ik,1};
    if length(ref1) ~= length(ref1_1)
        ref1=ref1(1:length(ref1_1));
        ref2=ref2(1:length(ref1_1));
    end
    %         [dur,dur_idx]=sort(ref2-ref1,'ascend');
    dur=ref2-ref1;
    %     dur_all(ik,1)=max(dur);
    ind_b1=[];
    ind_d1=[];
    if max(dur)>=290
        for t=1:length(bins)-1
            for i=1:length(dur)-1
                if (dur(i))>bins(t) && (dur(i))<=bins(t+1)
                    ind_b1=[ind_b1 ref1_1(dur==dur(i))];
                     ind_d1=[ind_d1 dur(i)];
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

BZ.idrat=[new_idx];
for ik=1:length(BZ.idrat)
    ref3=b_pt(ik,:);
    for ct=1:size(BZ.phase_thal{BZ.idrat(ik),1},1)
        clearvars -except ik ct BZ epochs_zd1 epochs_zd dur_all ref3 b_pt
        
        non_norm=unwrap(BZ.phase_ctx(BZ.idrat(ik),:))-unwrap(BZ.phase_thal{BZ.idrat(ik),1}(ct,:)); %circdist
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

% for i=1:length(BZ.idrat)
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
title('BZ')
%
%
% tempo=1:size(BZ.env_ctx,2);
% plot(tempo,BZ.env_ctx(BZ.idrat(ik),:))
% hold on
% plot(tempo,BZ.filt_ctx(BZ.idrat(ik),:))
% plot(tempo(ref3),BZ.env_ctx(BZ.idrat(ik),ref3),'b.')
% plot(tempo(ref1),BZ.env_ctx(BZ.idrat(ik),ref1),'ro')

