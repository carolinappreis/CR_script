clear all
close all
% load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\RS_input_9ch.mat')
% load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\RS_20fs_9ch.mat')
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/start_end_RS.mat')
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/RS_2xfpeak_9ch.mat')
iii=2;

d= [RS_t{iii,1}; RS_e{iii,1}; RS_dc{iii,1}];
% stackedplot(d');
% xlim([0 1200])

% data1=[RS_t{iii,1}];
data1=d';

% data=reshape(data1,1,size(data1,1)*size(data1,2));data=data';
T=repmat(t(iii),1,((length(data1)./t(iii))));
epochs=zeros(1,size(data,1));
epochs(1:t(iii):end)=1;
% sum(T)==numel(data)

for sm=1:length(xx)
    stim1(sm,1:t(iii))= repmat(xx(sm),1,t(iii));
end
stim=reshape(stim1',1,size(stim1,1)*size(stim1,2));
time_n=0:1/Fs:(size(data,1)-1)/Fs;


% %%% decide on order
% for order=1:10
% sample_rate=Fs;
% subplot(1,10,order)
%         mar = check_mar_spectra( data,size(data,1), order, sample_rate,0);
% 
%         [pxx_welch,f_welch] = pwelch( data,hamming(64),[],256,sample_rate, 'psd' );
% 
%         [pxx_yule,f_yule] = pyulear( data,order,1024,sample_rate );
% 
%         plot(f_welch,pxx_welch,'linewidth',3)
%         hold on
%         plot(f_yule,pxx_yule,'g--','linewidth',2)
%         plot(mar.freq_vect,squeeze(abs(mar.PSD(1,1,:))),'r:','linewidth',2)
%         xlabel('Frequency (Hz)')
%         ylabel('Power Spectral Density')
% end
% legend({'pwelch','pyulear','hmmmar'});
% 
% options = struct();
% options.initrep = 3; % Quin = 3
% options.K =4; % number of states
% options.order=9;
% options.useParallel=0;
% options.Fs=Fs;
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



[hmm, Gamma] = hmmmar(data1,T,options);

figure()
plot(data1)
hold on
plot(Gamma*2)
% plot(epochs,'k','LineWidth',2)
box('off')
legend ({'temor signal','state1','state2','state3','state4'});
legend('boxoff')

X=data;
fit = hmmspectramar(X,T,hmm,Gamma,options)
% fit = hmmspectramt(X,T,Gamma,options)

figure()
for ii=1:options.K
plot(fit.state(ii).f,fit.state(ii).psd)
hold on
end
legend({'state1','state2','state3','state4'})
legend('boxoff')

k=hmm.K;
[viterbipath] = hmmdecode(data,T,hmm,1);

maxFO = getMaxFractionalOccupancy(Gamma,T,options); % useful to diagnose if the HMM  % is capturing dynamics or grand between-subject  % differences (see Wiki)
FO = getFractionalOccupancy (Gamma,T,options); % state fractional occupancies per session
LifeTimes = getStateLifeTimes (Gamma,T,options); % state life times
Intervals = getStateIntervalTimes (Gamma,T,options); % interval times between state visits
SwitchingRate =  getSwitchingRate(Gamma,T,options); % rate of switching between stats

vpath=viterbipath; %vector with the most likely sequence of states
K=hmm.K;
onset=getStateOnsets (vpath,T,K);

ep=1:t(iii):size(data1,2)+1;
vp_ep=[];
for g=1:length(ep)
    if g+1<=length(ep)
%         vp_ep=[ vp_ep  vpath(ep(g):(ep(g+1)-1))]; %% vpath of the 5
%         seconds 
          vp_ep=[ vp_ep  vpath(((ep(g)+4*Fs)-1):(ep(g+1)-1))]; %% vpath of the last sec of stim
    end
end
vp_ep=vp_ep';

figure()
states={'posture';'spiral'};
for ms=1:2; %motor state
    for no=1:size(vp_ep,1)
        for st=1:k
            vp_res(no,st)= sum(vp_ep(no,:)==st)/size(vp_ep,2);
        end
    end
    
    phase=xx(ms:2:length(xx));
    vp=vp_res(ms:2:length(xx),:);
    
    tt=NaN(12,k);
    
    for w=1:12
        tt(w,:)=mean(vp(find(phase==w),:));
    end
    
   subplot(2,1,ms)
   bar(tt)
   box('off')
   title(states{ms,1})
  clear tt
end


% evokedStateProbability