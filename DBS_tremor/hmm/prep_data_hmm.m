clear all
iii=[1];

for numb=1;
    %     :length(iii);
    clearvars -except iii numb cc cond nostim
    DBS_Fpeak
    
    cond={'_NS_BL.mat';'_NS_PS.mat';'_HFS_PS.mat';};
    
    for cc=2:3;
%         load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),cond{cc,1}))
                load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),cond{cc,1}));
        
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
        
        Fs=20;
        time=0:1/samplerateold:(size(data,2)-1)/samplerateold;
        ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
        ts1=resample(ts,0:1/Fs:((size(data,2)-1)/samplerateold),'linear');
        tremor2(1:size(ts1.data,3))=ts1.data;
        time_n=0:1/Fs:(size(tremor2,2)-1)/Fs;
        
        if (Fpeak-2)>=1
            [b,a]=butter(2,[(Fpeak-2)/(0.5*Fs) (Fpeak+2)/(0.5*Fs)],'bandpass'); %15
        else
            [b,a]=butter(2,[(1)/(0.5*Fs) (Fpeak+2)/(0.5*Fs)],'bandpass'); %15
        end
        %         [b,a]=butter(2,[0.8/(0.5*samplerate) ],'low'); %15
        %         tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
        tremor_or=filtfilt(b,a,tremor2);
        tremor_or=zscore(tremor_or);
        envelope=abs(hilbert(tremor_or));
        
        if cc==1
%             rt=segm(tremor_or,cc,samplerate);
        elseif cc==2
            ns=segm(tremor_or,cc,Fs);
        elseif cc==3
            fs=segm(tremor_or,cc,Fs);
        end
        
    end
    DBS_rs_hmm
     %     clearvars -except nostimout iii numb cc cond  nostim
end

clearvars -except Fs rs ns fs

clear data
data=([ns])';


length_epoch=30*Fs;
T=repmat(length_epoch,1,length(data)./length_epoch);

options = struct();
options.K =3; % number of states
options.order=3;
options.useParallel=0;



[hmm, Gamma] = hmmmar(data,T,options);
close all
plot(data)
hold on
plot(Gamma*5)


X=data;
 fit = hmmspectramar(X,T,hmm,Gamma,options)
% fit = hmmspectramt(X,T,Gamma,options)

figure()
for i=1:options.K
plot(fit.state(i).psd)
hold on
end


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