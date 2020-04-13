clear all
% iii=[1 2 3 4 5 8 10 11 12 13 16 17 18];

iii=[2 3 4 5 8 10 11 13 16 17];
main=[1 1 3 1 3 3 3 3 1 1];
for numb=1:length(iii);
    clearvars -except iii PC A1 B1 numb nostim nostimout samplerate  vr main
          load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
%     load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(iii(numb)),'_baseline.mat'))
    
    
    data=SmrData.WvData;
    samplerateold=SmrData.SR;
    ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    tt=0:1/samplerate:(size(ds_data,2)-1)/samplerate;
    tre_3=ds_data([3 5 6],:); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    
    before_ns
    
    segmentb=hu{numb,:};
    segmente=hd{numb,:};
    
    
    handup=[];
    for i=1:length(segmentb)
        handup=[handup segmentb(i):segmente(i)]; %#ok<*AGROW>
    end
    clear i
    handup=sort(handup,'ascend');
    
    
    
    [b,a]=butter(2,[2/(0.5*samplerate) 10/(0.5*samplerate)],'bandpass'); %15
    
    
    tf_3=filtfilt(b,a,tre_3')*10*9.81/0.5;
    tf_3=tf_3';
    % tremor_or=zscore(tremor_or);
    dummy=(hilbert(tf_3))';
    envelope=abs(dummy);
    zenv=(abs(hilbert(zscore(tf_3))))';
    phase=angle(dummy);
    frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
    
    r=[];
    for ii=length(segmentb)
        sig{numb,1}(ii,:)=[ r tre_3(main(numb),segmentb(ii):segmentb(ii)+50000-1)];
    end
   
    
    clearvars -except nostimout iii numb PC A1 B1 iii stim nostim vr main sig
end

 spectrogram(sig,window,noverlap,nfft)
 
 fs=1000;  
 
 Nx = length(n);
nsc = floor(Nx/4.5);
nov = floor(nsc/2);
nff = max(256,2^nextpow2(nsc));

spectrogram(n,hamming(nsc),nov,nff);
 
 spectrogram(n,fs,[],fs,fs) 

 
 [s,ff,tt,p]=spectrogram(n,hanning(nfft),nfft/2,nfft*2,fs); surf(tt,ff,(p),'edgecolor','none'); axis tight;view(0,90)


 spectrogram(n(1,:),[],[],1000)
 
 x=0:1/samplerate:(size(n,2)-1)/samplerate;
n=sig{1,1}(1,:);

cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')

% save('cr_18.mat')
% save 'newnonstim2.mat'

% ANS_group=nostimout; clear nostimout
% AS_group=stim;
% clearvars  -except ANS_group AS_group
% cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
% save 'A_group'
% % NS=nostimout;
% % no_s=nostim;
% % clearvars  -except NS no_s
% % cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
% % % cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
% % % save 'amp_NS.mat'
