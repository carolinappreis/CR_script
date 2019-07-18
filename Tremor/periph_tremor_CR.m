clear all
iii=[1 2 3 4 5 6 8 10 11];

for numb=1:length(iii);
load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))
start_clean;

%%% re - estimate tremor characteristics
clear handup Pxx F frange Pxxrange Fpeak tremor_or dummy envelope phase frequency 

handup=[];
for i=1:length(start)
    handup=[handup start(i):ending(i)]; %#ok<*AGROW>
end
handup=sort(handup,'ascend');

[Pxx,F]=pwelch(tremor2(handup),samplerate,[],samplerate,samplerate);

frange=F(3:10);
Pxxrange=Pxx(3:10);

Fpeak=frange(find(Pxxrange==max(Pxxrange)));

if (Fpeak-2)>=1
    [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
    [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end
tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;

dummy=hilbert(tremor_or);
envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
phase=angle(dummy);
frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';

tremor=(data(3,:));% %score(:,1)';%
ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremorx(1:size(ts1.data,3))=ts1.data;
tremor=(data(6,:));% %score(:,1)';%
ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremory(1:size(ts1.data,3))=ts1.data;
tremor=(data(7,:));% %score(:,1)';%
ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremorz(1:size(ts1.data,3))=ts1.data;
if (Fpeak-2)>=1
    [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
    [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end
tremorxf=filtfilt(b,a,tremorx);
tremoryf=filtfilt(b,a,tremory);
tremorzf=filtfilt(b,a,tremorz);

for j=1:length(start)
    x=[tremorxf(start(j):ending(j));tremoryf(start(j):ending(j));tremorzf(start(j):ending(j))];
    [pc,score,latent,tsquare] = pca(x');
    xxx(j,1:3)=pc(1:3,1);
    ma(j)=(find(abs(xxx(j,1:3))==max(abs(xxx(j,1:3)))));
end

tremor_or2=NaN(length(start),1);

for i=1:length(start)
    if (~isnan(start(i)))
        tremor_or2(i,1)=(mean(envelope(ending(i)-1000:ending(i)))-mean(envelope(start(i)-1000:start(i))))/mean(envelope(start(i)-1000:start(i)));
    else
        tremor_or2(i,1)=NaN;
    end
end

clear tt tt2

k=1;
tt=NaN(20,12);
yy = xx ;

for s=1:size(tt,2)
    for i =1:1000;
        yy1=xx(randperm(size(xx,2)));
        tt2(1:sum(yy1==s),1)=tremor_or2(find(yy1==s));
        tt3(i,s)=nanmedian(tt2,1);
        clear tt2 
    end
end

for i=1:12
    tt(1:sum(xx==i),i)=tremor_or2(find(xx==i));
end

ttall (numb,:,:)=tt3;
% clearvars -except tt tt3 iii numb
% cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\A_PS')
% save (strcat('P0',num2str(iii(numb)),'_pha_suffle.mat'));
end
A_LS=ttall;
clearvars -except A_LS


% rr(1:size(tt,2))=mean(prctile(tt3,95));
% rr1(1:size(tt,2))=mean(prctile(tt3,25));
% bar(nanmedian(tt))
% hold on
% plot(rr,'k--','LineWidth',1.5)
% plot(rr1,'k--','LineWidth',1.5)
% box('off')
% title 'Significant phasic-stimulation effect' 
% xticklabels({'0' '30' '60' '90' '120' '150' '180' '210' '240' '270' '300' '330'})
% 
% 
% 
% % find((nanmedian(tt))>(prctile(tt3,95))|(nanmedian(tt))<(prctile(tt3,25)))
% 
% figure()
% likhood_amp=sum(tt>prctile(tt3,95)| tt>0)./sum(~isnan(tt));
% likhood_sup=sum(tt<prctile(tt3,25)| tt<0)./sum(~isnan(tt));
% likhood=[likhood_sup ; likhood_amp];
% bar(likhood')
% title ('Likelihood of significant amplification/supression effect')
% xticklabels({'0' '30' '60' '90' '120' '150' '180' '210' '240' '270' '300' '330'})
% legend('supression','amplification')
% legend('boxoff')
% box('off')



% close all
% figure()
% fig=gcf;
% fig.Color=[1 1 1];
% bar(nanmedian(tt))
% hold on
% stem(tt')
% xticklabels({'0','30','60','90','120','150','180','210','240','270','300','330'})
% % ylim([(-max((max(tt)))-0.1).*100 (max(max(tt))+0.1).*100])
% box('off')
% title ('P8')
% % 
% % figure()
% fig=gcf;
% fig.Color=[1 1 1];
% bar(100.*nanmedian(tt2))
% hold on
% stem((100.*tt2)')
% xticklabels({'0','30','60','90','120','150','180','210','240','270','300','330'})
% % ylim([(-max((max(tt2)))-0.1).*100 (max(max(tt2))+0.1).*100])
% box('off')
% title ('P8')
% 
% 
%
