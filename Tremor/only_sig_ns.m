clear all
close all
%  iii=[1 2 3 4 5 8 10 11 12 13 16 17 18];
% iii=[2 3 4 5 8 10 11 13 16 17];
iii=[3 5 8 17];

gg=[];
main=[1 1 3 1];
for numb= 1:length(iii);
%     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\reference_pls_ns.mat')
clearvars -except iii numb main ns_hu ns_ee ns_ss ns_change
%     load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(iii(numb)),'_baseline.mat'))

data=SmrData.WvData;
samplerateold=SmrData.SR;
ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
samplerate=1000;
tt=0:1/samplerate:(size(ds_data,2)-1)/samplerate;
tre_3=ds_data([3 5 6],:); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK

before_ns

load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\before_PLS.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/before_PLS.mat')
cell_b=[2 3 4 5 8 10 11 13 16 17];
ha=find(iii(numb)==(cell_b));
segmentb=hu{ha,:};
segmente=hd{ha,:};

handup=[];
for i=1:length(segmentb)
handup=[handup segmentb(i):segmente(i)]; %#ok<*AGROW>
end

clear i
handup=sort(handup,'ascend');
for aa=1:3
[Pxx,F]=pwelch(tre_3(aa,handup),samplerate,[],samplerate,samplerate);
frange=F(3:10);
Pxxrange=Pxx(3:10);
Freqpeak(aa,:)=frange(find(Pxxrange==max(Pxxrange)));
Ppeak(aa,:)=max(Pxxrange);
ps_curves(aa,:)=Pxx;
end
peak_ax=[(Freqpeak(find(Ppeak==max(Ppeak)))) (find(Ppeak==max(Ppeak)))];
Fpeak=peak_ax(1);
if (Fpeak-2)>=1
[b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
[b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end
tf_3=filtfilt(b,a,tre_3')*10*9.81/0.5;
tf_3=tf_3';

for i=1:3
ztf_3(i,:)=zscore(tf_3(i,:));
envelope(i,:)=abs(hilbert(tf_3(i,:)));
z_env(i,:)=abs(hilbert(zscore(tf_3(i,:))));
end

load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\pls_prop.mat','pls_b')
for i=1:length(segmentb)
ns_bhu(1,i)=mean(z_env(main(numb),segmentb(i)-1000:segmentb(i)));
ns_start(1,i)=mean(z_env(main(numb),(segmentb(i)+pls_b(numb)-1000):(segmentb(i)+pls_b(numb))));
ns_end(1,i)=mean(z_env(main(numb),segmente(i)-5000:segmente(i)-4000));
n_change(1,i)=(mean(envelope(main(numb),segmente(i)-1000:segmente(i)))-mean(envelope(main(numb),segmentb(i)-1000:segmentb(i))))./mean(envelope(main(numb),segmentb(i)-1000:segmentb(i)));        
end
ns_change(numb,:)=mean(n_change); clear change
ns_hu(numb,:)=mean(ns_bhu,2);
ns_ss(numb,:)=mean(ns_start,2);
ns_ee(numb,:)=mean(ns_end,2);

clearvars -except iii numb main ns_hu ns_ee ns_ss ns_change
clear numb in2
end