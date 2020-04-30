clear all
% iii=[1 2 3 4 5 8 10 11 12 13 16 17 18];

iii=[2 3 4 5 8 10 11 13 16 17];
for numb=1:length(iii);
    clearvars -except iii PC A1 B1 numb nostim nostimout samplerate  vr
        load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
%     load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(iii(numb)),'_baseline.mat'))
    

    data=SmrData.WvData;
    samplerateold=SmrData.SR;
    ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    tt=0:1/samplerate:(size(ds_data,2)-1)/samplerate;
    tre_3=ds_data([3 5 6],:); 
    emg=[ds_data(7,:);ds_data(8,:)];
    
       before_ns
      
    segmentb=hu{numb,:};
    segmente=hd{numb,:};
     handup=[];
    for i=1:length(segmentb)
        handup=[handup segmentb(i):segmente(i)]; %#ok<*AGROW>
    end
    clear i
     handup=sort(handup,'ascend');
     
    for ii=1:size(emg,1) 
     [Pxx,F]=pwelch(emg(ii,handup),samplerate,[],samplerate,samplerate);
     all_pxx(numb,ii,:)=Pxx;
    end
end