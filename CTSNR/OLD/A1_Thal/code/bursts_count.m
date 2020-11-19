
clear all
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')

load 'animal_lesion_nolesion.mat'

data= A(lesion,:);

data1= {'kjx127AA01@100-200_m.mat'
    'kjx132e01@100-200_m.mat'
    'kjx136b01@0-100_m.mat'
    'kjx140a01@0-100_m.mat'
    'kjx160c01@0-100_m.mat'
    'kjx166a01@0-100_m.mat'
    'kjx167c01@200-300_m.mat'};


bd_short=cell(1,1);
bd_long=cell(1,1);
bd_all=cell(1,1);
med_ctx=cell(1,1);
rr=1;
for i=1:size(data1,1)
    cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/MATLAB/KN')
    load (data1{i})
    
    samprateold=1/IpsiEEG.interval;
    WaveData(1,:)=IpsiEEG.values;
    WaveData=double(WaveData);
    ts=timeseries(WaveData(1,:),0:(1/samprateold):((size(WaveData,2)-1)/samprateold));
    ts1=resample(ts,0:0.001:((size(WaveData,2)-1)/samprateold),'linear');
    WaveData_DC(1,:)=ts1.data;
    cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/code')
    WaveData_DC(1,:)=makemua_CR_1(WaveData(1,:)',0.001,0.003,round(samprateold),1000,4);
    samprate=1000;
    time=0:(1/samprate):(length(WaveData_DC(1,:))-1)*(1/samprate);
    muafilt=[];
    [Pxx_ind,F_ind]=pwelch(WaveData_DC(1,:),samprate,[],samprate,samprate);
    Pxx_ind_beta=Pxx_ind(16:36);
    filtrange=14+find(Pxx_ind_beta==max(Pxx_ind_beta));
    [b,a]=butter(2,[(filtrange-5)/(0.5*samprate) (filtrange+5)/(0.5*samprate)],'bandpass');
    muafilt(1,1:length(time)) = filtfilt(b,a,WaveData_DC(1,:));
    
    run 'bursts_only.m'
    bd_short{rr,:}=duration1{1};
    bd_long{rr,:}=duration1{2};
    bd_all{rr,:}=duration;
    med_ctx{rr,:}=output_ctx;
    rr=rr+1;
end

clearvars -except bd_short bd_long bd_all med_ctx

for j=1:size(bd_all,1)
bd_all1(j,1)=mean(bd_all{j});
bd_short1(j,1)=mean(bd_short{j});
bd_long1(j,1)=mean(bd_long{j});
end


bd_all_s= [mean(bd_all1);std(bd_all1)];
bd_short_s= [mean(bd_short1);std(bd_short1)];
bd_long_s= [mean(bd_long1);std(bd_long1)];


for i=1:size(med_ctx,1)
short_bursts(i,:)=med_ctx{i}(1);
long_bursts(i,:)=med_ctx{i}(2);
end

short_bursts=cell2mat(short_bursts);
long_bursts=cell2mat(long_bursts);

short_bursts_s=prctile(short_bursts,[25 50 75]);
long_bursts_s=prctile(long_bursts,[25 50 75]);

for i=1:size(bd_short,1)
bn_short(i,1)=size(bd_short{i},2);
bn_long(i,1)=size(bd_long{i},2);
bn_all(i,1)=size(bd_all{i},2);
end


bn_short_s=[mean(bn_short);std(bn_short)];
bn_long_s=[mean(bn_long);std(bn_long)];
bn_all_s=[mean(bn_all);std(bn_all)];


bursts_duration=cell(3,1);
bursts_duration={[bd_short_s];[bd_long_s];[bd_all_s]};

bursts_nmbr=cell(3,1);
bursts_nmbr={[bn_short_s];[bn_long_s];[bn_all_s]};

bursts_plot=cell(2,1);
bursts_plot={[short_bursts_s];[long_bursts_s]};

cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')
clearvars -except bursts_plot bursts_nmbr bursts_duration
save 'bursts_lesioned'
