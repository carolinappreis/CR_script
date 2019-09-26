clear all
load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\RS_input_9ch.mat')

i=1;
data=[RS_t{i,1}; RS_e{i,1}; RS_dc{i,1}];
% data=reshape(data1,1,size(data1,1)*size(data1,2));

length_epoch=t(i);
T=repmat(length_epoch,1,((length(data)./length_epoch)*9))';

sum(T)==numel(data)

options = struct();
options.K =3; % number of states
options.order=3;
options.useParallel=1;



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