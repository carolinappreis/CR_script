clear all
iii=[1];
numb=1;
Fpeak=4;

% clearvars -except data Fs segmentb segmente a b ns fs iii numb in2 time_n
% load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))
% load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\start_end_RS.mat')
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/start_end_RS.mat')
load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_RS_PS.mat'))
data=SmrData.WvData;
samplerateold=SmrData.SR;
%%------------------------
tremor3=data([3 6 7],:);

[b,a]=butter(2,[2/(0.5*samplerateold) 15/(0.5*samplerateold)],'bandpass'); %15
[p,q]=butter(3,[5/(0.5*samplerateold)],'low');
[m,n]=butter(3,[1.5/(0.5*samplerateold)],'low');


for i=1:size(tremor3,1)
    t1(i,:)=filtfilt(b,a,tremor3(i,:));
    t2(i,:)=filtfilt(p,q,tremor3(i,:));
    t3(i,:)=filtfilt(m,n,tremor3(i,:));
end
data1=vertcat(t1,t2,t3,tremor3);

Fs=1000;
time=0:1/samplerateold:(size(data1,2)-1)/samplerateold;
ts=timeseries(data1,0:(1/samplerateold):((size(data1,2)-1)/samplerateold));
ts1=resample(ts,0:1/Fs:((size(data1,2)-1)/samplerateold),'linear');
data2(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
time2=0:1/Fs:(size(data2,2)-1)/Fs;


filt_t1=data2(1:3,:);
filt_t2=data2(4:6,:);
filt_t3=data2(7:9,:);
tre=data2(10:12,:);

% s=round((start1*Fs)./samplerateold,0);
% e=round((ending1*Fs)./samplerateold,0);
clear RS_tf1 RS_tf2 RS_tf3 RS_tf4
s=round((start2*Fs)./samplerateold,0);
e=round((ending2*Fs)./samplerateold,0);

t=5*Fs;
for tr=1:length(s)
    for th=1:size(tremor3,1)
        
        RS_tf1{tr,1}(th,:)=filt_t1(th,s(tr):(s(tr)+(t/2)-1));
        RS_tf2{tr,1}(th,:)=filt_t2(th,s(tr):(s(tr)+(t/2)-1));
        RS_tf3{tr,1}(th,:)= filt_t3(th,s(tr):(s(tr)+(t/2)-1));
        RS_tf4{tr,1}(th,:)= tre(th,s(tr):(s(tr)+(t/2)-1));
    end
    
end
i=3;
figure()
plot3(RS_tf1{i,1}(1,:),RS_tf1{i,1}(2,:), RS_tf1{i,1}(3,:))
figure()
plot3(RS_tf2{i,1}(1,:),RS_tf2{i,1}(2,:), RS_tf2{i,1}(3,:))
figure()
plot3(RS_tf3{i,1}(1,:),RS_tf3{i,1}(2,:), RS_tf3{i,1}(3,:))
figure()
plot3(time2(1:length(RS_tf2{i,1}(2,:))),RS_tf2{i,1}(2,:), RS_tf2{i,1}(3,:))
figure()
plot3(time2(1:length(RS_tf1{i,1}(2,:))),RS_tf1{i,1}(1,:), RS_tf1{i,1}(3,:))

