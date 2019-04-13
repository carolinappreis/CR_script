clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')

load('SNR.mat');
%  load('SNR.mat');


for ik=1
%     :size(SNR.env_ctx,1)
    for pp=1:2;
    ref=[];
        non_nomr=[];epochs_idx=[];epochs_t=[];
        ref=SNR.onset_raw{1,ik}{pp,1};
        
        el=400;
        for ii=1:length(ref)
            if ref(ii)>el
                epochs_ctx(pp,ii,:)=((SNR.env_ctx(ik,ref(ii)-el:ref(ii)+el)-median(SNR.env_ctx(ik,ref(ii)-el:ref(ii))))./median(SNR.env_ctx(ik,ref(ii)-el:ref(ii)))).*100;
            end
        end
    end
end

time=1:801;
color_b=[0.8 0.8 0.8];
bb=squeeze(epochs_ctx(1,:,:));
y2=mean(bb,1);
y1=mean(bb,1)+std(bb,1)./sqrt(size(bb,1));
y3=mean(bb,1)-std(bb,1)./sqrt(size(bb,1));
p1=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
clearvars -except time epochs_ctx p1 
color_b1=[0.5 0.5 0.5];
bb=squeeze(epochs_ctx(2,:,:));
y2=mean(bb,1);
y1=mean(bb,1)+std(bb,1)./sqrt(size(bb,1));
y3=mean(bb,1)-std(bb,1)./sqrt(size(bb,1));
p2=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b1)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b1],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b1],'FaceAlpha',[0.2],'EdgeColor','none')


legend([p1 p2],{'short burst','long burst'})
ylabel ('ECoG Beta Amplitude (%Baseline)')
xlabel ('Time (msec)')
xticks([0:200:800])
yticks([0:100:200])
xticklabels ({'-400','-200','0','200','400'})
xlim ([200 800])
ylim([-50 200])
box ('off')
title('SNR recordings')