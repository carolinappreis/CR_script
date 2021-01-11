
clear all
close all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/final_mats')
load('BZ_bua.mat');BZ_bua=bua; clear bua
load('SNr_bua.mat');SNR_bua=bua; clear bua
r=0;
for i=size(BZ_bua,1)-1:size(BZ_bua,1)
    r=r+1;
    Ecog{r,1}=BZ_bua{i,1}(1,:);
    SNr{r,1}=SNR_bua{i,1}(2:end,:);
    BZ{r,1}=BZ_bua{i,1}(2:end,:);
end

clearvars -except Ecog BZ SNr
region={'Ecog' 'BZ' 'SNr'};


c=0;
p=0;
samprate=1000;

 for pr=1:size(Ecog,1)
c=c+1;
for ct=1:size(Ecog{pr,1},1)
    for ct2=1:size(BZ{pr,1})
    freq=1:80;
    [Pxx_ind,F_i]=mscohere((Ecog{pr,1}(ct,:)),(BZ{pr,1}(ct2,:)),samprate,[],samprate,samprate);
    p=p+1;
    coher(p,:)=Pxx_ind(freq)./sum(Pxx_ind);clear Pxx_ind
    end
end
end

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean','grey');
color_b={(grey+blood)/2};
plt={'coher'};
label={'Mag-Squared Coherence';};
lims={[0 0.04]};
fig=figure;
time=freq;
ii=1;
data=eval(plt{ii,1});
y2=mean(data);
% y1=y2+(std(data)./sqrt(size(data,1)));
% y3=y2-(std(data)./sqrt(size(data,1)));
y1=y2+std(data);
y3=y2-std(data);
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
fig.OuterPosition= [10, 10, 10, 10];
fig.Color='w';
set(gca,'FontSize',12)
ylabel(sprintf(label{ii,1}));
xlabel('Frequency (Hz)');
