
% data_new_cat : data; T_new_cat : T

% hmm : model 

% Gamma : hmm posterior probability 

% options : hmm options 

%

T_cat=T_new_cat;

maxFO = getMaxFractionalOccupancy(Gamma,T,options);  % maximum fractional isotropy of the subjects 

FO = getFractionalOccupancy (Gamma,T,options); % state fractional occupancies per session

Intervals = getStateIntervalTimes (Gamma,T,options); % interval times between state visits

SwitchingRate =  getSwitchingRate(Gamma,T,options); % rate of switching between stats

lifetimes=getStateLifeTimes(Gamma,T,options);   % SBK it was T_new{sub} . , options

 

lifetimes=cellfun(@(x) x/Fs,lifetimes,'UniformOutput', false);

mean_state_lifetimes=cell2mat(cellfun(@mean,lifetimes,'UniformOutput', false));

lifetimes_old=lifetimes;

lifetimes=[];

for kk=1:size(lifetimes_old,2)

    tmp=cell2mat(lifetimes_old(:,kk)');

    lifetimes{kk}=tmp;

end


 % Lifetime

figure;

[h]=distributionPlot(lifetimes,'histOpt',2,'color','c');

h{2}(1).Marker='.';

h{2}(1).MarkerSize=40;

h{2}(1).Color='r';

h{2}(2).Marker='.';

h{2}(2).MarkerSize=40;

h{2}(2).Color='k';

 

plot4paper('State','Dwell Time (s)');
title ('Lifetime');

 

tempX=[output,'1Subj1ch_Lags', num2str(L),'_stateNO',num2str(N_states),'_sub_',num2str(sub),'_Lifetime','_chNO_',num2str(ch),'.png'];

saveas(gca,tempX); 

 

% FO

figure 

h = imagesc(FO);

% title ('FO');

ylabel('Trials');

xlabel('States');

xtick = 1: size (FO,2) ;

ytick = 1:size (FO,1)

    ax = gca;

 ax.XTick = xtick;

%  ax.YTick = ytick;

colormap jet

colorbar

 title ('Fractional Occupancy ');

 

tempX=[output,'1Subj1ch_Lags', num2str(L),'_stateNO',num2str(N_states),'_sub_',num2str(sub),'_FO_','_chNO',num2str(ch),'.png'];

saveas(gca,tempX); 


%  transition probabilies from any state to any other state with considering staying at one state 

 

k=size(Gamma,2);

figure

    [~,pca1] = pca(Gamma','NumComponents',1);

    [~,ord] = sort(pca1); 

    P = hmm.P;

    for j=1:k, P(j,j) = 0; P(j,:) = P(j,:) / sum(P(j,:));  end

    imagesc(P(ord,ord)); colorbar

    axis square

    hold on

    for j=0:13

        plot([0 13] - 0.5,[j j] + 0.5,'k','LineWidth',2)

        plot([j j] + 0.5,[0 13] - 0.5,'k','LineWidth',2)

    end

    hold off

xtick = 1: size(Gamma,2);

    ax = gca;

 ax.XTick = xtick;

 ax.YTick = xtick;

  xlabel('States');

 ylabel('States');

 title ('Trans. Prob. with persistence probabilities');

 

tempX=[output,'1Subj1ch_Lags', num2str(L),'_stateNO',num2str(N_states),'_sub_',num2str(sub),'_TransOrdinSameState','_chNO_',num2str(ch),'.png'];

saveas(gca,tempX); 


% transition probabilies from any state to any other state,

 % without considering the persistence probabilities (i.e. the probability

 % to remain in the same state)

figure

TP = getTransProbs(hmm);

P = TP;

for j=1:k, P(j,j) = 0

%      P(j,:) = P(j,:) / sum(P(j,:));

end

imagesc(P(ord,ord)); colorbar

axis square

hold on

for j=0:13

   plot([0 13] - 0.5,[j j] + 0.5,'k','LineWidth',2)

   plot([j j] + 0.5,[0 13] - 0.5,'k','LineWidth',2)

end

hold off

 title ('Trans. Prob. without persistence probabilities');

 

tempX=[output,'1Subj1ch_Lags', num2str(L),'_stateNO',num2str(N_states),'_sub_',num2str(sub),'_TransWithoutSameState','_chNO_',num2str(ch),'.png'];

saveas(gca,tempX); 


 

% If order==0 (Gaussian distribution), these purely represents functional connectivity; If order>0 (MAR), these refer to the covariance matrix of the residual

figure 

 

for k1 = 1 : size (Gamma,2)

    if isempty(lifetimes{1, k1})  % in case a state is empty; to not consider it fr the FunConn

        continue

    end 

funcon(:,:,k1) = getFuncConn(hmm,k1);

 

subplot(2,3,k1)  % Change

temp = funcon(:,:,k1);

temp = corrcov(temp);  % Compute correlation matrix from covariance matrix

imagesc(temp)

colormap(parula)

colorbar

xtick = 1: length(subjects) ;

    ax = gca;

  title(['States ',num2str(k1)]);

end 

hold off

 

tempX=[output,'1Subj1ch_Lags', num2str(L),'_stateNO',num2str(N_states),'_sub_',num2str(sub),'_CorrCovFunConn','_chNO_',num2str(ch),'.png'];

saveas(gca,tempX); 

 


%   fractional occupancy correlation matrix

figure 

    FOfig = corr(FO);

    P = FOfig;

      for j=1:k

          P(j,j) = 0; 

%        P (j,:) = P(j,:) / sum(P(j,:)); 

      end

    imagesc(P(ord,ord)); colorbar

    axis square

    hold on

    for j=0:13

        plot([0 13] - 0.5,[j j] + 0.5,'k','LineWidth',2)

        plot([j j] + 0.5,[0 13] - 0.5,'k','LineWidth',2)

    end

    hold off

xtick = 1: k ;

    ax = gca;

 ax.XTick = xtick;

 ax.YTick = xtick;

  ylabel('States');

    xlabel('States');

    colormap(jet)

     title ('Corr. ( FO )');

 

    

 tempX=[output,'1Subj1ch_Lags', num2str(L),'_stateNO',num2str(N_states),'_sub_',num2str(sub),'_CorrFO','_chNO_',num2str(ch),'.png'];

saveas(gca,tempX); 

 

 

% Max FO histogram plots

 

figure 

 histogram(maxFO), ylabel('Frequency of Trials');

 title ('Hist. ( maxFO )');

 tempX=[output,'1Subj1ch_Lags', num2str(L),'_stateNO',num2str(N_states),'_sub_',num2str(sub),'_HistMaxFO','_chNO_',num2str(ch),'.png'];

saveas(gca,tempX); 


 data_plot=data_new_cat;

    reg=6;
    options.order = 0;

    options.win = reg * sample_rate; % multitaper window 

 

    T_cat_z = sum(T_cat);
    fit = hmmspectramt(data_plot,T_cat_z,Gamma,options);

% frequency bands of interest 

[wt,wf] = cwt(data_plot(:,1),'amor',sample_rate);

theta_wavelet_freqs = wf > 4 & wf < 7;

alpha_wavelet_freqs = wf > 8 & wf < 12;

beta_wavelet_freqs = wf > 13 & wf < 30;

 

% % tbins = (1*sample_rate:40*sample_rate);  %For selected part of the signal

tbins = 1:length(data_plot);   % For all the signal

time_s = tbins/sample_rate; 

Gamma_R1 = Gamma;

%  

% %%%%%%%%%%%%%%%%%%%%%%%%

% detect which state is alpha, low and high beta

% theta/alpha/Beta oscillations

hilb_theta = mean(abs(wt(theta_wavelet_freqs,:)),1)';

ts=[];

for kk=1:size(Gamma,2)

    tmp=corr(Gamma(:,kk),hilb_theta);

    

    % ols 

    [cope, varcope]=ols(demean(logit(Gamma(:,kk))),demean(hilb_theta),[1]);

    ts(kk)=cope/sqrt(varcope);

end

[max_tstat theta_state]=max(ts);

 

[a b]=max(Gamma');

theta_mask_gamma = (b==theta_state)';

valid_theta_state = max(max(mean_state_lifetimes))<10 && max_tstat>20;

 

 

hilb_alpha = mean(abs(wt(alpha_wavelet_freqs,:)),1)';

ts=[];

for kk=1:size(Gamma,2)

    tmp=corr(Gamma(:,kk),hilb_alpha);

    

    % ols 

    [cope, varcope]=ols(demean(logit(Gamma(:,kk))),demean(hilb_alpha),[1]);

    ts(kk)=cope/sqrt(varcope);

end

[max_tstat alpha_state]=max(ts);

 

[a b]=max(Gamma');

alpha_mask_gamma = (b==alpha_state)';

valid_alpha_state = max(max(mean_state_lifetimes))<10 && max_tstat>20;

 

%  

hilb_beta = mean(abs(wt(beta_wavelet_freqs,:)),1)';

ts=[];

for kk=1:size(Gamma,2)

    tmp=corr(Gamma(:,kk),hilb_beta);

    

    % ols 

    [cope, varcope]=ols(demean(logit(Gamma(:,kk))),demean(hilb_beta),[1]);

    ts(kk)=cope/sqrt(varcope);

end

[max_tstat beta_state]=max(ts);

 

[a b]=max(Gamma');

 valid_beta_state = max(max(mean_state_lifetimes))<10 && max_tstat>20;

% valid_beta_state = max(mean_state_lifetimes)<10 && max_tstat>20;

beta_mask_gamma = (b==beta_state)';

 

%

 

if valid_theta_state

    disp(['theta state is #', num2str(theta_state)]);

else

    disp(['INSUFFICIENT EVIDENCE FOR A THETA STATE']);

    disp(['Best candidate theta state found was #', num2str(theta_state)]);

end

 

if valid_alpha_state

    disp(['alpha state is #', num2str(alpha_state)]);

else

    disp(['INSUFFICIENT EVIDENCE FOR A ALPHA STATE']);

    disp(['Best candidate alpha state found was #', num2str(alpha_state)]);

end

 

 

if valid_beta_state

    disp(['beta state is #', num2str(beta_state)]);

else

    disp(['INSUFFICIENT EVIDENCE FOR A BETA STATE']);

    disp(['Best candidate beta state found was #', num2str(beta_state)]);

end

%%


VE_reg = normalise(data_plot); % data

VE_reg=VE_reg(1:2000);

tbins=tbins(1:2000);

 Gamma=Gamma(1:2000,:);

time_s=time_s(1:2000);

% highbeta_state = lowbeta_state;

cols={'b','g','r','y','p','m','c'};

legs={'1','2','3','4','5','6','7'};

legs{theta_state}='theta';

legs{alpha_state}='alpha';

% legs{highbeta_state}='high beta';

% legs{lowbeta_state}='low beta';

legs{beta_state}='beta';

% legs{lowgamma_state}='low gamma';

% legs{highgamma_state}='high gamma';

legs=legs(1:size(Gamma,2));

 

legs_thresh{1}='theta';

legs_thresh{2}='alpha';

legs_thresh{3}='beta';

 

 %%

figure('units','normalized','outerposition',[0 0 1 1])

 

% subplot_grid(6,2,1,1)

subplot(6,8,[1 2])

hold on;grid on;

% addpath(genpath('DD:\HMM-MAR-master_Feb6th\HMM-MAR-master\spectral\'));

 

mar = check_mar_spectra( VE_reg, size(VE_reg,1), order, sample_rate,0);

[pxx_welch,f_welch] = pwelch( VE_reg,hamming(64),[],256,sample_rate, 'psd' );

 

plot(f_welch,pxx_welch,'linewidth',3)

plot(mar.freq_vect,squeeze(abs(mar.PSD(1,1,:))),'r:','linewidth',2)

        xlabel('Freq., Hz');

% subplot_grid(6,2,1,2)

subplot(6,8,[3:8])

hold on;grid on;

for ii = 1:size(Gamma,2)

%     plot(fit.state(ii).f,abs(fit.state(ii).psd),cols{ii},'linewidth',3);%

%    if its just one channel 

%     plot(fit.state(ii).f,squeeze(abs(fit.state(ii).psd(:,ch,ch))),cols{ii},'linewidth',3); % if ita 

    plot(fit.state(ii).f,squeeze(abs(fit.state(ii).psd)),cols{ii},'linewidth',3); % if ita 

 

end

 legend(legs);

        xlabel('Freq., Hz');

 

%%%%%%%%%%%%%%%%%%%%%%%%

% plot Gamma

 

% subplot(612);

subplot(6,8,9:16)

 

hold on;

for ii = 1:size(Gamma,2)

    plot(time_s, Gamma(tbins,ii),cols{ii},'linewidth',3);

end

xlim([time_s(1), time_s(end)]);

xlabel('time, s');

 

% legend(legs);

 

%%%%%%%%%%%%%%%%%%%%%%%%

% plot data with state tcs

 

% subplot(613)

subplot(6,8,17:24)

 

for kk=1:size(Gamma,2)

 

    [a b]=max(Gamma');

    mask_gamma = (b==kk)';

    burst_vals_state = VE_reg.*mask_gamma;

    burst_vals_state(burst_vals_state == 0) = nan;

 

    hold on; plot(time_s, burst_vals_state(tbins),cols{kk},'LineWidth',2);    

end

% legend(legs);

% set(gca,'FontSize',12); 

xlim([time_s(1), time_s(end)]);

xlabel('time, s');

 

%%%%%%%%%%%%%%%%%%%%%%%%

% plot cwt

 

% subplot(614)

subplot(6,8,25:32)

 

contourf( time_s,wf, sqrt(abs(wt(:,tbins))), 36 ,'linestyle','none' )

grid on;hold on

xlim([time_s(1), time_s(end)]);

% title('Morlet wavelet transform ')

xlabel('time, s');

%%%%%%%%%%%%%%%%%%%%%%%%

% plot state covariance matrices

 

% state specific cov mats

 

clear covmat corrmat;

% if do_embedded

    for kk=1:size(Gamma,2)

        

          if isempty(lifetimes{1, kk})  % in case a state is empty; to not consider it fr the FunConn

        continue

          end 

    

        [covmat{kk},corrmat{kk}]=getFuncConn(hmm,kk);

    end

% end

 

for kk=1:size(Gamma,2)

%     subplot_grid(6,options.K,5,kk);

%         subplot(6,options.K,kk+20);

 

  if isempty(lifetimes{1, kk})  % in case a state is empty; to not consider it fr the FunConn

        continue

    end 

 

subplot(6,8,32+kk)

 

    h=imagesc(covmat{kk},[0 1]);

    colorbar;

            title(['state',num2str(kk)]);

 

end

 

%%%%%%%%%%%%%%%%%%%%%%%%

% compute phase stability from instantaneous hilbert freq

 

freqbands=[4,45];

dat=data_plot;

datbp = bandpass(dat,freqbands,sample_rate);            

env = transpose(abs(hilbert(datbp')));

phase=(angle(hilbert(datbp')));

instfreq = sample_rate/(2*pi)*diff(unwrap(phase));

 

[aa bb]=max(Gamma');

 

clear phase_stability amplitude instantaneous_freqkk;

for kk = 1:size(Gamma,2)

 

    disp(['Computing for state ' num2str(kk)]);

 

    inds = logical(bb==kk)';

 

    instantaneous_freqkk{kk}=instfreq(inds(2:end));

    phase_stability(kk)=1./std(instfreq(inds(2:end)));

    amplitude{kk}=env(inds);    

end

 

%%%%%%%%%%%%%%%%%%%%%%%%

% plot state-wise distributions of instantaneous hilbert freq and amp

 

for kk = 1:size(Gamma,2)

 

%     subplot_grid(6,2,6,1);

%         subplot(616);

subplot(6,8,[41: 44])

 

 

    [nn,xx]=hist(instantaneous_freqkk{kk},50);

    plot(xx,nn/sum(nn),cols{kk},'LineWidth',2);

    hold on;

    xlabel('Hilb Freq.');

%         xlabel('Freq., Hz');

 

    

%     subplot_grid(6,2,6,2);

%             subplot(6,2,12);

subplot(6,8,[45: 48])

 

 

    [nn,xx]=hist(amplitude{kk},50);

    plot(xx,nn/sum(nn),cols{kk},'LineWidth',2);

    xlabel('Hilb Amp.');

%             xlabel('Freq., Hz');  

 

    

    hold on;

end

 Thresh_masks=Thresh_masks';

 

 

%

    

disp(['phase_stability=', num2str(phase_stability)]);

 

 tempX=[output,'1Subj1ch_Lags', num2str(L),'_stateNO',num2str(N_states),'_sub_',num2str(sub),'_Summary','_chNO_',num2str(ch),'.png'];

saveas(gca,tempX); 