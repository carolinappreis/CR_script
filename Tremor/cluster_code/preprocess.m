function [d]=preprocess(SmrData)
d=struct;
d.data_raw=SmrData.WvData;
d.samplerateold=SmrData.SR;
ts=timeseries(d.data_raw,0:(1/d.samplerateold):((size(d.data_raw,2)-1)/d.samplerateold));
ts1=resample(ts,0:0.001:((size(d.data_raw,2)-1)/d.samplerateold),'linear');
ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
d.samplerate=1000;
tt=0:1/d.samplerate:(size(ds_data,2)-1)/d.samplerate;
d.data_ds=ds_data([3 5 6],:);
end