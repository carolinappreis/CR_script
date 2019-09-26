clear all
% load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\RS_input_9ch.mat')
load('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/hmm/RS_input_9ch.mat')
i=2;
data1=[RS_dc{i,1}];
% ; RS_e{i,1}; RS_dc{i,1}];
data=reshape(data1,1,size(data1,1)*size(data1,2));
length_epoch=t(i);
data=data(1:(end-60))';

T=ones(ceil(size(data,1)/length_epoch),1)*length_epoch; 

% T(end)=rem(size(data',1),length_epoch);


% T=repmat(length_epoch,1,((length(data)./length_epoch)*9))';

sum(T)==numel(data)

options = struct();
options.K =3; % number of states
options.order=3;
options.useParallel=0;
options.Fs=12;
option.win=18;



for order=1:10
figure()


sample_rate=12;

        mar = check_mar_spectra( data,size(data,1), order, sample_rate,0);

        [pxx_welch,f_welch] = pwelch( data,hamming(64),[],256,sample_rate, 'psd' );

        [pxx_yule,f_yule] = pyulear( data,order,1024,sample_rate );

        plot(f_welch,pxx_welch,'linewidth',3)

        plot(f_yule,pxx_yule,'g--','linewidth',2)

        plot(mar.freq_vect,squeeze(abs(mar.PSD(1,1,:))),'r:','linewidth',2)

        xlabel('Frequency (Hz)')

        ylabel('Power Spectral Density')

%         xlim([0 50]);

%         title(['subj' num2str(Good_subjects_sub(sub-3))]);

%         ylim([0 0.27]);

end



% hold on;
% 
%        legend({'pwelch','pyulear','hmmmar'});   
%     
% end
% 


[hmm, Gamma] = hmmmar(data,T,options);
close all
plot(data)
hold on
plot(Gamma*5)


X=data;
%  fit = hmmspectramar(X,T,hmm,Gamma,options)
 fit = hmmspectramt(X,T,Gamma,options)

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