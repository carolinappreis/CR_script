clear all
% load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\RS_input_9ch.mat')
load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\RS_20fs_9ch.mat')
i=1;
data1=[RS_t{i,1}];
% ; RS_e{i,1}; RS_dc{i,1}];
data=reshape(data1,1,size(data1,1)*size(data1,2));data=data';
T=repmat(t(i),1,((length(data)./t(i))));

time_n=0:1/Fs:(size(data,1)-1)/Fs;

sum(T)==numel(data)


% decide on order
for order=1:10
sample_rate=Fs;
subplot(1,10,order)
        mar = check_mar_spectra( data,size(data,1), order, sample_rate,0);

        [pxx_welch,f_welch] = pwelch( data,hamming(64),[],256,sample_rate, 'psd' );

        [pxx_yule,f_yule] = pyulear( data,order,1024,sample_rate );

        plot(f_welch,pxx_welch,'linewidth',3)
        hold on
        plot(f_yule,pxx_yule,'g--','linewidth',2)
        plot(mar.freq_vect,squeeze(abs(mar.PSD(1,1,:))),'r:','linewidth',2)
        xlabel('Frequency (Hz)')
        ylabel('Power Spectral Density')
end

legend({'pwelch','pyulear','hmmmar'});

options = struct();
options.initrep = 3; % Quin = 3 
options.K =3; % number of states
options.order=5;
options.useParallel=1;
options.Fs=Fs;
% option.win=18;

%Saed's options
% options = [];
% options.initrep = 10; % Quin = 3 
% options.initcyc = 1e3;
% options.K = 3;
% options.standardise = 1; %  SBK = 1; Quin = 0
% % options.leakagecorr=-1;  % New Added by SBK
%  options.verbose = 1;
% options.Fs = Fs;
% options.useMEX = 1;
% options.zeromean = 0;
% options.dropstates = 0;
% options.useParallel = 1; % set to 1 if you have access to the parallel processing toolbox

 
%%% embedded...
% options.order = 0;
% options.embeddedlags = -5:5;
% options.covtype = 'full';



[hmm, Gamma] = hmmmar(data,T,options);
figure()
plot(data)
hold on
plot(Gamma*2)


X=data;
 fit = hmmspectramar(X,T,hmm,Gamma,options)
% fit = hmmspectramt(X,T,Gamma,options)

figure()
for i=1:options.K
plot(fit.state(i).f,fit.state(i).psd)
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