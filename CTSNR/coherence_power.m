clear all
close all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/final_mats')
load('BZ_coh.mat');
raw_coh=BZ_raw_coh;

% clear all
% close all
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/final_mats')
% load('SNr_bua.mat');
% raw_coh=bua;


c=0;
p=0;
samprate=1000;
%   for pr=1:size(raw_coh,1)
  for pr=size(raw_coh,1)-1:size(raw_coh,1)
c=c+1;
for ct=2:size(raw_coh{pr,1},1)
    freq=1:80;
    [Pxx_ind,F_i]=mscohere((raw_coh{pr,1}(1,:)),(raw_coh{pr,1}(ct,:)),samprate,[],samprate,samprate);
    p=p+1;
    coher(p,:)=Pxx_ind(freq)./sum(Pxx_ind);clear Pxx_ind
    [Pxx_ctx,F_ind]=pwelch(raw_coh{pr,1}(1,:),samprate,[],samprate,samprate);
    power_ctx(c,:)=Pxx_ctx(freq)./sum(Pxx_ctx); clear Pxx_ctx
    [Pxx_probe,F_ind]=pwelch((raw_coh{pr,1}(ct,:)),samprate,[],samprate,samprate);
    power_probe(p,:)=Pxx_probe(freq)./sum(Pxx_probe); clear Pxx_probe
    
end
end

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean','grey');
if name=='bz'
    color_b={grey blood (grey+blood)/2};
else
    color_b={grey aegean (grey+aegean)/2};
end

plt={'power_ctx';'power_probe';'coher'};
label={'Normalised Power' ; 'Normalised Power';'Mag-Squared Coherence';};
lims={[0 0.2];[0 0.0050]; [0 0.04]};
fig=figure;
% time=F_i';
time=freq;
for ii=1:size(plt,1)
    subplot(1,size(plt,1),ii)
    data=eval(plt{ii,1});
    y2=mean(data);
%     y1=y2+(std(data)./sqrt(size(data,1)));
%     y3=y2-(std(data)./sqrt(size(data,1)));
    y1=y2+(std(data));
    y3=y2-(std(data));
    %         if ii~=3
    %             y1=log(y1);y2=log(y2);y3=log(y3);
    %         end
    plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',[color_b{1,ii}])
    hold on
    patch([time fliplr(time)], [y1 fliplr(y2)],[color_b{1,ii}],'FaceAlpha',[0.2],'EdgeColor','none')
    patch([time fliplr(time)], [y2 fliplr(y3)],[color_b{1,ii}],'FaceAlpha',[0.2],'EdgeColor','none')
    xlim ([0 80])
    xticks([0:10:80])
    ylim([lims{ii,1}]);
    box('off')
    fig.Units = 'centimeters';
    fig.OuterPosition= [10, 10, 30, 10];
    fig.Color='w';
    set(gca,'FontSize',12)
    ylabel(sprintf(label{ii,1}));
    xlabel('Frequency (Hz)');
end


%%%% COHERENCE 1 RAT FROM SELECT CHANNELS 10% COH CRITERIA
% % % clear all
% % % close all
% % % cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/final_mats')
% % % load('BZ_coh.mat');
% % % load('SNr_coh.mat');
% % % r=0;
% % % for i=size(BZ_raw_coh,1)-1:size(BZ_raw_coh,1)
% % %     r=r+1;
% % %     Ecog{r,1}=BZ_raw_coh{i,1}(1,:);
% % %     SNr{r,1}=SNR_raw_coh{i,1}(2:end,:);
% % %     BZ{r,1}=BZ_raw_coh{i,1}(2:end,:);
% % % end
% % % 
% % % clearvars -except Ecog BZ SNr
% % % region={'Ecog' 'BZ' 'SNr'};
% % % 
% % % 
% % % 
% % % c=0;
% % % p=0;
% % % samprate=1000;
% % %  for pr=1:size(BZ,1)
% % % c=c+1;
% % % for ct=1:size(BZ{pr,1},1)
% % %     for ct2=1:size(SNr{pr,1})
% % %     freq=1:80;
% % %     [Pxx_ind,F_i]=mscohere((BZ{pr,1}(ct,:)),(SNr{pr,1}(ct2,:)),samprate,[],samprate,samprate);
% % %     p=p+1;
% % %     coher(p,:)=Pxx_ind(freq)./sum(Pxx_ind);clear Pxx_ind
% % %     end
% % % end
% % % end
% % % 
% % % load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean','grey');
% % % color_b={(aegean+blood)/2};
% % % plt={'coher'};
% % % label={'Mag-Squared Coherence';};
% % % lims={[0 0.04]};
% % % fig=figure;
% % % time=freq;
% % % ii=1;
% % % data=eval(plt{ii,1});
% % % y2=mean(data);
% % % % y1=y2+(std(data)./sqrt(size(data,1)));
% % % % y3=y2-(std(data)./sqrt(size(data,1)));
% % % y1=y2+std(data);
% % % y3=y2-std(data);
% % % %         if ii~=3
% % % %             y1=log(y1);y2=log(y2);y3=log(y3);
% % % %         end
% % % plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',[color_b{1,ii}])
% % % hold on
% % % patch([time fliplr(time)], [y1 fliplr(y2)],[color_b{1,ii}],'FaceAlpha',[0.2],'EdgeColor','none')
% % % patch([time fliplr(time)], [y2 fliplr(y3)],[color_b{1,ii}],'FaceAlpha',[0.2],'EdgeColor','none')
% % % xlim ([0 80])
% % % xticks([0:10:80])
% % % ylim([lims{ii,1}]);
% % % box('off')
% % % fig.Units = 'centimeters';
% % % fig.OuterPosition= [10, 10, 10, 10];
% % % fig.Color='w';
% % % set(gca,'FontSize',12)
% % % ylabel(sprintf(label{ii,1}));
% % % xlabel('Frequency (Hz)');
