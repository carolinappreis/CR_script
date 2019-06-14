clear all
% cd('/Users/Carolina/Documents/MATLAB/Tremor')
cd('C:\Users\creis\Documents\MATLAB\pt_data_periphstim')
load ('p08_randomstim_cursos.mat')
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
         tremor_sor2(i,1)=(mean(envelope(sur_end(i)-1000:sur_end(i)))-mean(envelope(sur_start(i)-1000:sur_start(i))))/mean(envelope(sur_start(i)-1000:sur_start(i)));

    else
        tremor_or2(i,1)=NaN;
        tremor_sor2(i,1)=NaN;
    end
end


clear tt
k=1;
tt=NaN(20,12);
tt2=NaN(20,12);

for i=1:12
    tt(1:sum(xx==i),i)=tremor_or2(find(xx==i));
    tt2(1:sum(xx==i),i)=tremor_sor2(find(xx==i));
end

close all
figure()
fig=gcf;
fig.Color=[1 1 1];
bar(100.*nanmedian(tt))
hold on
stem((100.*tt)')
xticklabels({'0','30','60','90','120','150','180','210','240','270','300','330'})
% ylim([(-max((max(tt)))-0.1).*100 (max(max(tt))+0.1).*100])
box('off')
title ('P8')

figure()
fig=gcf;
fig.Color=[1 1 1];
bar(100.*nanmedian(tt2))
hold on
stem((100.*tt2)')
xticklabels({'0','30','60','90','120','150','180','210','240','270','300','330'})
% ylim([(-max((max(tt2)))-0.1).*100 (max(max(tt2))+0.1).*100])
box('off')
title ('P8')




st=NaN(1,size((tt),1));
clear A; A=tt; %b1{f,1};
clear B; B=tt2; %s1{f,1}(1:size(A,1),:);
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.01 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)];
    
end


