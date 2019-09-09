
clear all
close all
 load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\01_NS_PS.mat'))
in2=1; % analysing the "main tremor axis"

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

time=0:1/samplerateold:(size(data,2)-1)/samplerateold;

%%% downsample

ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremor2(1:size(ts1.data,3))=ts1.data;
samplerate=1000;
time2=0:1/samplerate:(size(tremor2,2)-1)/samplerate;

Fpeak=4.5;
if (Fpeak-2)>=1
    [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
    [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end
tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
dummy=hilbert(tremor_or);
envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));

clear data


data1=envelope;
length_epoch=60000;
options = struct();
options.K =4; % number of states
Fs=samplerate;
T=repmat(length_epoch,1,round(length(data1)./length_epoch)-1);
data=data1(1:sum(T))';

[hmm, Gamma, Xi, vpath, GammaInit, residuals, fehist] = hmmmar (data,T,options);

Tsubject=length(data);
Y=hmm.state;
options_test = struct();

t = hmmtest(Gamma,T,Tsubject,Y,options_test,hmm);
