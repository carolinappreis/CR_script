%TC PROJECT v.2

%-----------------------BUA-A4 Analysis----------------------------------------
 cd ('C:\Users\creis\Documents\GitHub\CR_script\A4_Thal\code')
 cd('/Users/Carolina/Documents/GitHub/CR_script/A4_Thal/code') 

%phase synchrony index(psi)consistency during ctx/subctx bursts vs surrogates between ctx/subcortical regions
run 'combine2.m'
run 'exploring_psi.m'
run 'plots_psi_time_bursts.m'
% imporatnt functions are bursts.m burts_off.m bursts_aligned.m hayriye_c
% bursts_aligned_off.m

%phase-slips analysis
run 'phase_slips.m' % as described in PNAS - Cagnan et al.2019
run 'during.m' % according to bins of busrt duration

%PSD ctx and subcortical regions
run 'psd_plots'

%-----------------------SUA Juxta_SUA_act  Analysis------------------------------
cd ('C:\Users\creis\Documents\GitHub\CR_script\SUA\Juxta SUA_act_mat')
cd('/Users/Carolina/Documents/GitHub/CR_script/SUA/Juxta SUA_act_mat') 

%phase consistency of subcortical firing of subcortical areas
run 'psi_SUA.m'

% Single unit firing rate and phase consistency of subcortical firing in
% bursts vs surrogates && ISI, firing rate per region 
run 'Filegen_SUA_act_new.m'

% Sum fo spiking across bursts in time.
run 'gaussian2.m'

%Cycle by cycle analysis
run 'CY_CY.m' % runs function 'cycles_10.m'

%-----------------------SUA Probe_SUA_act  Analysis------------------------------
  cd ('C:\Users\creis\Documents\GitHub\CR_script\SUA\Probe SUA_act_mat')
%  cd('/Users/Carolina/Documents/GitHub/CR_script/SUA/Probe SUA_act_mat') 

%Find region, spike sorting , creat files (newSUA_SNR.mat)
run 'spikes_makemua.m' % will run make_CR_spikes.m function

%Firing-rate of units into alpha, beta ,gamma
run 'process_SUA.m'