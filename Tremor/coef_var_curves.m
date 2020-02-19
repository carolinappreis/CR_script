clear all
close all
%  iii=[1 2 3 4 5 8 10 11 12 13 16 17 18];
iii=[2 3 4 5 8 10 11 13 16 17];

gg=[];
main=[1 1 3 1 3 3 3 3 1 1];
for numb= 1:length(iii);
%     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\reference_pls_ns.mat')
clearvars -except iii numb main
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
segmentb=hu{numb,:};
segmente=hd{numb,:};

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
env(i,:)=abs(hilbert(tf_3(i,:)));
phase(i,:)=angle(hilbert(tf_3(i,:)));
freq(i,:)=(smooth((1000/(2*pi))*diff(unwrap(phase(i,:))),500))';
end

ns_amp=env(main(numb),handup);
ns_freq=freq(main(numb),handup);
min(ns_freq)
amp_bined=[];
dum2=[];
bins=3.5:0.25:9;
    for i=1:length(bins)
        if i+1<length(bins)
         dum=find(ns_freq>bins(i) & ns_freq<=bins(i+1));
         dum2=[dum2 ; i.* ones(length(dum),1)];
         amp_bined=[amp_bined  ns_amp(dum)];
        end
    end
    
    h_amp(numb,1:2,:)=[dum2'; amp_bined];
    figure(1)
    subplot(5,2,numb)
    histogram(squeeze(h_amp(numb,:,:)));
end



clearvars -except iii numb ns_filt_mainax in2
clear numb in2