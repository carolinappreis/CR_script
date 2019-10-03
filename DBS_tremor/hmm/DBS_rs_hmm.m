clear all
iii=[1];
numb=1;
DBS_Fpeak

% clearvars -except data Fs segmentb segmente a b ns fs iii numb in2 time_n
% load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))
% load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\start_end_RS.mat')
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/start_end_RS.mat')
 load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_RS_PS.mat'))

in2=1;

if in2==1
    in=3;
elseif in2==2 % other axis 1
    in=6;
elseif in2==3 % other axis 2
    in=7;
end
data=SmrData.WvData;
samplerateold=SmrData.SR;
tremor=(data(in,:));
addon=92; addon_end=35;

%%------------------------

Fs=(Fpeak+2)*2+5;
tremor3=data([3 6 7],:);

[b,a]=butter(2,[(Fpeak-2)/(0.5*Fs) ((Fpeak+2))/(0.5*Fs)],'bandpass');

[m,n]=butter(2,[1.5/(0.5*Fs) ],'low');

for i=1:size(tremor3,1)
    dc_t3(i,:)=filtfilt(m,n,tremor3(i,:));
    filt_t3(i,:)=filtfilt(b,a,tremor3(i,:));
    env_t3(i,:)=abs(hilbert(filt_t3(i,:)));
end
data1=vertcat(filt_t3,env_t3,dc_t3);

time=0:1/samplerateold:(size(data1,2)-1)/samplerateold;
ts=timeseries(data1,0:(1/samplerateold):((size(data1,2)-1)/samplerateold));
ts1=resample(ts,0:1/Fs:((size(data1,2)-1)/samplerateold),'linear');
data2(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;


filt_t3=data2(1:3,:);
env_t3=data2(4:6,:);
dc_t3=data2(7:9,:);

% all conditions
% s=round((start*Fs)./samplerateold,0);
% e=round((ending*Fs)./samplerateold,0);

%just posture
s=round((start1*Fs)./samplerateold,0);
e=round((ending1*Fs)./samplerateold,0);

%just spiral
% s=round((start2*Fs)./samplerateold,0);
% e=round((ending2*Fs)./samplerateold,0);

%start to length of epochs decided by t;
% t=[Fs 5*Fs];
% for ep=1:size(t,2);
%     for th=1:size(tremor3,1);
%         for i=1:length(t(ep))
%             for tr=1:length(s)
%                 RS_t1{ep,1}(th,tr,1:t(ep))=filt_t3(th,s(tr):(s(tr)+t(ep)-1));
%                 RS_e1{ep,1}(th,tr,1:t(ep))=env_t3(th,s(tr):(s(tr)+t(ep)-1));
%                 RS_dc1{ep,1}(th,tr,1:t(ep))=dc_t3(th,s(tr):(s(tr)+t(ep)-1));
%             end
%         end
%        
%         RS_t{ep,1}=reshape(RS_t1{ep,1},size(RS_t1{ep,1},1),size(RS_t1{ep,1},2)*size(RS_t1{ep,1},3));
%         RS_e{ep,1}=reshape(RS_e1{ep,1},size(RS_e1{ep,1},1),size(RS_e1{ep,1},2)*size(RS_e1{ep,1},3));
%         RS_dc{ep,1}=reshape(RS_dc1{ep,1},size(RS_dc1{ep,1},1),size(RS_dc1{ep,1},2)*size(RS_dc1{ep,1},3));
%         
%     end
% end

%once sec before:one sec after stim
t=5*Fs;
for ep=1:size(t,2);
    for th=1:size(tremor3,1)
            RS_raw1=[];
            RS_t1=[];
            RS_e1=[];
%             RS_dc1=[];
        for i=1:length(t(ep))
            for tr=1:length(s)
                RS_raw1=[RS_raw1 tremor3(th,s(tr)-Fs:(s(tr)+Fs+t(ep)-1))];
                RS_t1=[RS_t1 filt_t3(th,s(tr)-Fs:(s(tr)+Fs+t(ep)-1))];
                RS_e1=[RS_e1 env_t3(th,s(tr)-Fs:(s(tr)+Fs+t(ep)-1))];
%                 RS_dc1=[RS_dc1 dc_t3(th,s(tr)-Fs:(s(tr)+Fs+t(ep)-1))];
            end
        end

        RS_r{ep,1}(th,:)=RS_raw1;
        RS_t{ep,1}(th,:)=RS_t1;
        RS_e{ep,1}(th,:)=RS_e1;
%         RS_dc{ep,1}(th,:)=RS_dc1;
        
    end
end



% cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA')
cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/');

clearvars -except RS_t RS_e RS_dc Fs t RS_raw
save 'RS_hmm_posture.mat'

