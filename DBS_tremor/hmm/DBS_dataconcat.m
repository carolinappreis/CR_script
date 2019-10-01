clear all
iii=[1];
numb=1;
DBS_Fpeak

% clearvars -except data Fs segmentb segmente a b ns fs iii numb in2 time_n
load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))
load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\start_end_RS.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/start_end_RS.mat')
% load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_RS_PS.mat'))

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
Fs=20;

time=0:1/samplerateold:(size(data,2)-1)/samplerateold;
ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:1/Fs:((size(data,2)-1)/samplerateold),'linear');
data2(1:size(data,1),1:size(ts1.data,3))=ts1.data;
tremor3=data2([3 6 7],:);

for u=1:size(tremor3,1)
    z_tremor3(u,:)=zscore(tremor3(u,:));
end

[b,a]=butter(2,[(2)/(0.5*Fs) ((Fs/2)-1)/(0.5*Fs)],'bandpass'); 

[m,n]=butter(2,[1.5/(0.5*Fs) ],'low');
for i=1:size(tremor3,1)
    dc_t3(i,:)=filtfilt(m,n,z_tremor3(i,:));
    filt_t3(i,:)=filtfilt(b,a,z_tremor3(i,:));
    env_t3(i,:)=abs(hilbert(filt_t3(i,:)));
end

% all conditions
% s=round((start*Fs)./samplerateold,0);
% e=round((ending*Fs)./samplerateold,0);

%just posture
% s=round((start1*Fs)./samplerateold,0);
% e=round((ending1*Fs)./samplerateold,0);

%just spiral
s=round((start2*Fs)./samplerateold,0);
e=round((ending2*Fs)./samplerateold,0);

%once sec before or after stim
t=[Fs];
for ep=1:size(t,2);
    for th=1:size(tremor3,1)
            RS_raw1=[];
            RS_t1=[];
            RS_e1=[];
            RS_dc1=[];
        for i=1:length(t(ep))
            for tr=1:length(s)
                RS_raw1=[RS_raw1 tremor3(th,s(tr)-Fs:(s(tr)-1))];
                RS_t1=[RS_t1 filt_t3(th,s(tr)-Fs:(s(tr)-1))];
                RS_e1=[RS_e1 env_t3(th,s(tr)-Fs:(s(tr)-1))];
                RS_dc1=[RS_dc1 dc_t3(th,s(tr)-Fs:(s(tr)-1))];
            end
        end

        RS_r{ep,1}(th,:)=RS_raw1;
        RS_t{ep,1}(th,:)=RS_t1;
        RS_e{ep,1}(th,:)=RS_e1;
        RS_dc{ep,1}(th,:)=RS_dc1;
        
    end
end



cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA')
% cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/');

clearvars -except RS_t RS_e RS_dc Fs t RS_raw
save 'RS_20fs_9ch_spiral_before.mat'

