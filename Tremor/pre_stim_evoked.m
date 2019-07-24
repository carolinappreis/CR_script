time=0:1/samplerateold:(size(data,2)-1)/samplerateold;
samplerate=1000;

tremor=(data(3,:));% %score(:,1)';%
ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremor1(1:size(ts1.data,3))=ts1.data; clear tremor

tremor=(data(5,:));% %score(:,1)';%
ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremor2(1:size(ts1.data,3))=ts1.data; clear tremor

tremor=(data(6,:));% %score(:,1)';%
ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremor3(1:size(ts1.data,3))=ts1.data; clear tremor

stim1=(data(2,:));
ts=timeseries(stim1,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
stim(1:size(ts1.data,3))=ts1.data;

labels1=data(4,:);
ts=timeseries(labels1,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
labels(1:size(ts1.data,3))=ts1.data;

emg=data(7,:);
ts=timeseries(emg,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
emg1(1:size(ts1.data,3))=ts1.data;

[b,a]=butter(2,[(1)/(0.5*samplerate) (35)/(0.5*samplerate)],'bandpass'); %15
[bb,aa] = butter(3,10/(samplerate/2),'high');

tremor1f=filtfilt(b,a,tremor1);
tremor2f=filtfilt(b,a,tremor2);
tremor3f=filtfilt(b,a,tremor3);
emg1=filtfilt(bb,aa,emg1);


%%% determine stimulation time points
index=[];
for i=2:size(stim,2)-1
    if stim(i-1)<2.5 && stim(i)>2.5
        index=[index i];
    end
end
clear i


dd2=round(labels*100)./100;
for i=1:length(index)
    xx(i)=round(dd2(index(i))./0.1); %#ok<*SAGROW>
end
clear i