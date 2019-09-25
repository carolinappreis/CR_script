clear all
iii=[1];

for numb=1;
    %     :length(iii);
    clearvars -except nostimout iii numb cc cond nostim
    
    cond={'_NS_BL.mat';'_NS_PS.mat';'_HFS_PS.mat'};
    cc=2;
         load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),cond{cc,1}))
%     load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),cond{cc,1}));
    
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
    
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremor2(1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    
    nsfs=segm_nfilt(tremor2,cc);
%     clearvars -except nostimout iii numb cc cond  nostim
end


clear data
% data1=nsfs;
data=nsfs';
length_epoch=30000;
options = struct();
options.filter=[2 10];
options.Fs=1000;
options.onpower=1;
options.downsample=20;
options.K =3; % number of states
T=repmat(length_epoch,1,length(data)./length_epoch);
% data=data1(1:sum(T))';

[hmm, Gamma] = hmmmar(data,T,options);
[viterbipath] = hmmdecode(data,T,hmm,1);
k=hmm.K;

maxFO = getMaxFractionalOccupancy(Gamma,T,options); % useful to diagnose if the HMM  % is capturing dynamics or grand between-subject  % differences (see Wiki)
FO = getFractionalOccupancy (Gamma,T,options); % state fractional occupancies per session
LifeTimes = getStateLifeTimes (Gamma,T,options); % state life times
Intervals = getStateIntervalTimes (Gamma,T,options); % interval times between state visits
SwitchingRate =  getSwitchingRate(Gamma,T,options); % rate of switching between stats

vpath=viterbipath;
K=hmm.K;
getStateOnsets (vpath,T,Hz,K)