clear all
iii=[1];
numb=1;
DBS_Fpeak

% clearvars -except data Fs segmentb segmente a b ns fs iii numb in2 time_n
load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))
load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\start_end_RS.mat')
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
% 2*(Fpeak+2);

time=0:1/samplerateold:(size(data,2)-1)/samplerateold;
ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:1/Fs:((size(data,2)-1)/samplerateold),'linear');
data2(1:size(data,1),1:size(ts1.data,3))=ts1.data;
tremor3=data2([3 6 7],:);
time_n=0:1/Fs:(size(data2,2)-1)/Fs;

[b,a]=butter(2,[2/(0.5*Fs) (Fs./2)/(0.5*Fs)],'bandpass'); %15
%         [b,a]=butter(2,[0.8/(0.5*samplerate) ],'low'); %15
%         tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;

[m,n]=butter(2,[3/(0.5*Fs) ],'low'); %15
for i=1:size(tremor3,1)
    dc_t3(i,:)=filtfilt(m,n,zscore(tremor3(i,:)));
    filt_t3(i,:)=filtfilt(b,a,zscore(tremor3(i,:)));
    env_t3(i,:)=abs(hilbert(filt_t3(i,:)));
end

s=round((start*Fs)./samplerateold,0);
e=round((ending*Fs)./samplerateold,0);
t=[Fs 5*Fs];
for ep=1:size(t,2);
    for th=1:size(tremor3,1);
        for i=1:length(t(ep))
            for tr=1:length(s)
                if ((tr)+t(ep)-1)<length(s) && (s(tr)+t(ep)-1)<length(tremor3)
                    RS_t1{ep,1}(th,tr,1:t(ep))=filt_t3(th,s(tr):(s(tr)+t(ep)-1));
                    RS_e1{ep,1}(th,tr,1:t(ep))=env_t3(th,s(tr):(s(tr)+t(ep)-1));
                    RS_dc1{ep,1}(th,tr,1:t(ep))=dc_t3(th,s(tr):(s(tr)+t(ep)-1));
                end
            end
        end
    end
    RS_t{ep,1}=reshape(RS_t1{ep,1},size(RS_t1{ep,1},1),size(RS_t1{ep,1},2)*size(RS_t1{ep,1},3));
    RS_e{ep,1}=reshape(RS_e1{ep,1},size(RS_e1{ep,1},1),size(RS_e1{ep,1},2)*size(RS_e1{ep,1},3));
    RS_dc{ep,1}=reshape(RS_dc1{ep,1},size(RS_dc1{ep,1},1),size(RS_dc1{ep,1},2)*size(RS_dc1{ep,1},3));
end

cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA')
clearvars -except RS_t RS_e RS_dc Fs t
save 'RS_20fs_9ch.mat'

