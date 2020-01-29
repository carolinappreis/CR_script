clear all
close all
%  iii=[1 2 3 4 5 8 10 11 12 13 16 17 18];
iii=[2 3 4 5 8 10 11 13 16 17];
in2=1;
gg=[];
main=[1 1 3 1 3 3 3 3 1 1];
for numb= 1:length(iii);
%     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\reference_pls_ns.mat')
clearvars -except iii numb ns_filt_mainax in2 main
%     load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(iii(numb)),'_baseline.mat'))
if in2==1
in=3;
elseif in2==2 % other axis 1
in=5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
elseif in2==3 % other axis 2
in=6; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
end
data=SmrData.WvData;
samplerateold=SmrData.SR;
tremor=(data(in,:));
addon=92; addon_end=35;
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
%         figure()
%         plot(F(3:50),ps_curves(:,3:50)','LineWidth',2)
%         legend({'x','y','z'})
%         legend('boxoff')
%         box('off')
%
if (Fpeak-2)>=1
[b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
[b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end
tf_3=filtfilt(b,a,tre_3')*10*9.81/0.5;
tf_3=tf_3';
for i=1:3
ztf_3(i,:)=zscore(tf_3(i,:));
end

ns_filt_mainax{numb,1}=ztf_3(main(numb),handup);
end
clearvars -except iii numb ns_filt_mainax in2
clear numb in2