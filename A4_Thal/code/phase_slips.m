clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
  cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')

load('SNR_opt.mat');
% for ii=1:length(SNR.idrat)
% data{ii,1}=SNR.filt_thal{SNR.idrat(ii),1}
% end
% data=vertcat(data{:});
f=1;
SNR.idrat=[1:size(SNR.env_ctx,1)];
for ik=1:length(SNR.idrat)
    ref1=SNR.offset_raw{1,(SNR.idrat(ik))}{2,1};
     ref1_1=SNR.onset_phase_al{1,(SNR.idrat(ik))}{2,1};%%% ONSET
%     ref1_1=SNR.offset_phase_al{1,(SNR.idrat(ik))}{2,1};%%% OFFSET
    ref2=SNR.offset_raw{1,SNR.idrat(ik)}{2,1};
    if length(ref1) ~= length(ref1_1)
        ref1=ref1(1:length(ref1_1));
        ref2=ref2(1:length(ref1_1));
        f=f+1;
    end
    [dur,dur_idx]=sort(ref2-ref1,'ascend');
    th(ik,:)=(numel(dur_idx(dur>250)));
    dur_all(ik,:)=dur(1,end-24:end);
    ref3=ref1_1(dur_idx(end-24:end));
    for ct=1:size(SNR.phase_thal{SNR.idrat(ik),1},1)
        clearvars -except ik ct SNR epochs_zd1 epochs_zd dur_all ref3 f th dur_idx
        
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
        
%                         plot(epochs_z1(ii,:),'r.')
%                         imagesc(epochs_z1)
%                         xticks([200:200:800])
%                         xlim([200 800])
%                         xticklabels ({'-200','0','200','400'})
        epochs_zd(ct,:,:)=epochs_z1;
    end
    epochs_zd1(ik,:,:)=squeeze(mean(epochs_zd,1));
end


% for i=1:length(SNR.idrat)
%  imagesc(squeeze(epochs_zd1(i,:,:)))
%  close all
% end

slip_b=squeeze(mean(epochs_zd1,1));

figure()
imagesc(slip_b)
xlabel ('Time(msec)')
ylabel('Bursts # (Sorted by length)')

%%%%% ONSET
xticks([200:200:800])
xlim([200 800])
xticklabels ({'-200','0','200','400'})
title('SNR aligned to burst onset')

%%%% OFFSET
% xticks([0:200:600])
% xlim([0 600])
% xticklabels ({'-400','-200','0','200'})
% title('SNR aligned to burst offset')

figure()
plot(smooth(sum(slip_b)),'LineWidth',1.5)
xlabel ('Time(msec)')
ylabel('Sum zcores across Bursts')
box('off')

%%% ONSET
xticks([200:200:800])
xlim([200 800])
xticklabels ({'-200','0','200','400'})
title('SNR aligned to burst onset')



%%% OFFSET
% xticks([0:200:600])
% xlim([0 600])
% xticklabels ({''-400','-200','0','200'})
% title('SNR aligned to burst offset')

%
% 
% tempo=1:size(SNR.env_ctx,2);
% plot(tempo,SNR.env_ctx(SNR.idrat(ik),:))
% hold on
% plot(tempo,SNR.filt_ctx(SNR.idrat(ik),:))
% plot(tempo(ref3),SNR.env_ctx(SNR.idrat(ik),ref3),'b.')
% plot(tempo(ref1),SNR.env_ctx(SNR.idrat(ik),ref1),'ro')

