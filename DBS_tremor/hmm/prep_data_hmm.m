clear all
iii=[1];

for numb=1;
    %     :length(iii);
    clearvars -except nostimout iii numb cc cond nostim
    DBS_Fpeak
    
    cond={'_NS_BL.mat';'_NS_PS.mat';'_HFS_PS.mat'};
    
    for cc=1:3;
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
        
        if (Fpeak-2)>=1
            [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
        else
            [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
        end
        %         [b,a]=butter(2,[0.8/(0.5*samplerate) ],'low'); %15
        %         tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
        tremor_or=filtfilt(b,a,tremor2);
        tremor_or=zscore(tremor_or);
        envelope=abs(hilbert(tremor_or));
        
        if cc==1
            rt=segm(envelope,cc);
        elseif cc==2
            ns=segm(envelope,cc);
        elseif cc==3
            fs=segm(envelope,cc);
        end
        
    end
    rnsfs=[rt ns fs];
    %     clearvars -except nostimout iii numb cc cond  nostim
end



clear data
data1=rnsfs;
length_epoch=30000;
options = struct();
options.K =6; % number of states
options.order=0;
Fs=samplerate;
T=repmat(length_epoch,1,length(data1)./length_epoch);
data=data1';

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
onset=getStateOnsets (vpath,T,K);