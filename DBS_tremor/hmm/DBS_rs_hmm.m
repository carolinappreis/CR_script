clear all
iii=[1];
numb=1;
DBS_Fpeak

% clearvars -except data Fs segmentb segmente a b ns fs iii numb in2 time_n
% load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))
% load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\start_end_RS.mat')
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/start_end_RS.mat')
load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_RS_PS.mat'))
data=SmrData.WvData;
samplerateold=SmrData.SR;
%%------------------------
tremor3=data([3 6 7],:);

[b,a]=butter(2,[(Fpeak-2)/(0.5*samplerateold) (Fpeak+2)/(0.5*samplerateold)],'bandpass');

[m,n]=butter(3,[1.5/(0.5*samplerateold)],'low');

for i=1:size(tremor3,1)
    dc_t1(i,:)=filtfilt(m,n,tremor3(i,:));
    filt_t1(i,:)=filtfilt(b,a,tremor3(i,:));
    env_t1(i,:)=abs(hilbert(filt_t1(i,:)));
end
data1=vertcat(filt_t1,env_t1,dc_t1);

Fs=(Fpeak+2)*2+5;
time=0:1/samplerateold:(size(data1,2)-1)/samplerateold;
ts=timeseries(data1,0:(1/samplerateold):((size(data1,2)-1)/samplerateold));
ts1=resample(ts,0:1/Fs:((size(data1,2)-1)/samplerateold),'linear');
data2(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
time2=0:1/Fs:(size(data2,2)-1)/Fs;



filt_t3=data2(1:3,:);
env_t3=data2(4:6,:);
dc_t3=data2(7:9,:);

% all conditions
% s=round((start*Fs)./samplerateold,0);
% e=round((ending*Fs)./samplerateold,0);

%just posture

s=round((start1*Fs)./samplerateold,0);
e=round((ending1*Fs)./samplerateold,0);


% jj=3;%% axis
% 
% s=s(pca_idx{1,1}(:)==jj);

    
% % just spiral
% s=round((start2*Fs)./samplerateold,0);
% e=round((ending2*Fs)./samplerateold,0);


% t=5*Fs;
% for th=jj
%     RS_raw1=[];
%     RS_t1=[];
%     RS_e1=[];
%     RS_dc1=[];
%     for tr=1:length(s)
%         RS_raw1=[RS_raw1 tremor3(th,s(tr)-Fs:(s(tr)+t-1))];
%         RS_t1=[RS_t1 filt_t3(th,s(tr)-Fs:(s(tr)+t-1))];
%         RS_e1=[RS_e1 env_t3(th,s(tr)-Fs:(s(tr)+t-1))];
%         RS_dc1=[RS_dc1 dc_t3(th,s(tr)-Fs:(s(tr)+t-1))];
%     end
%     
%     
%     RS_r(1,:)=RS_raw1;
%     RS_t(1,:)=RS_t1;
%     RS_e(1,:)=RS_e1;
%     RS_dc(1,:)=RS_dc1;
%     
% end



t=5*Fs;
for th=1:size(tremor3,1)
    RS_raw1=[];
    RS_t1=[];
    RS_e1=[];
    RS_dc1=[];
    seg=[];
    for tr=1:length(s)
        RS_raw1=[RS_raw1 tremor3(th,s(tr)-Fs:(s(tr)+t-1))];
        RS_t1=[RS_t1 filt_t3(th,s(tr)-Fs:(s(tr)+t-1))];
        RS_e1=[RS_e1 env_t3(th,s(tr)-Fs:(s(tr)+t-1))];
        RS_dc1=[RS_dc1 dc_t3(th,s(tr)-Fs:(s(tr)+t-1))];
        seg=[seg s(tr)-Fs:(s(tr)+t-1)];
    end
    
    
    RS_r(th,:)=RS_raw1;
    RS_t(th,:)=RS_t1;
    RS_e(th,:)=RS_e1;
    RS_dc(th,:)=RS_dc1;
    
end

t=t+Fs;

% cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA')
cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/');

clearvars -except RS_t RS_e RS_dc Fs t RS_raw
save 'RS_hmm_posture.mat'

