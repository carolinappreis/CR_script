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
run 'during_severalb.m' % according to bins of busrt duration
run 'idregion_phaseslip.m' % change of instantaneous frequency at ctx and thal independantly
%PSD ctx and subcortical regions
run 'psd_plots'

%----------------------- Juxta_SUA_act  Analysis------------------------------
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

cd('/Users/Carolina/Documents/GitHub/CR_script/SUA/Probe SUA_act_mat') 
cd('C:\Users\creis\Documents\GitHub\CR_script\SUA\probe SUA_act_mat')


 % Select rats at region and that have shown coh with ctx in the bua
 % analysis/ were chosen by KN in his excel file
 run 'filegen_slected_rats.m'
 
 %Sum fo spiking across bursts in time.
 run 'gaussian_probe.m'
 
 % burst triggered avg
 run 'triggered_avg_sua.m' % data will come from 'spikerate_sua.m'
 
 %correct version of thw above
 run 'new_triggered_sua.m'
 
 %correlation spikes snr bz same rat
 run 'correlation_sua.m'
 
%stats on phase alignement for cycle by clycle analysis
run 'cy_cy_sua_stats.m'  % data_SUA_SNR.mat

%cycle by cycle as Andy
run 'CY_CY_sua.m'

%cycle by cycle HC
run 'CY_CY_HC.m'

%-----------------------MUA------------------------------------------------------
cd('C:\Users\creis\Documents\GitHub\CR_script\MUA')
cd('/Users/Carolina/Documents/GitHub/CR_script/MUA')

%Find region, spike sorting , creat files (newSUA_SNR.mat)
run 'spikes_makemua.m' % will run make_CR_spikes.m function

%Firing-rate of units into alpha, beta ,gamma
run 'process_MUA.m'

% triggered averaged from bursts aligned and non-aligned
run 'triggered_avg_mua.m'

%cycle by cucle analysis prbe
run 'cy_cy_probe.m' % comes from 'cy_cy_probe_stats.m'


%------------------------ plots poster

% psd and coherence
 cd ('C:\Users\creis\Documents\GitHub\CR_script\A4_Thal\code')
 cd('/Users/Carolina/Documents/GitHub/CR_script/A4_Thal/code') 
 run('psd_plots.m')
 run('coherence_plots.m')
 
%amplitude change, PSI 
cd('C:\Users\creis\Documents\GitHub\CR_script')
cd('/Users/Carolina/Documents/GitHub/CR_script')
run('plots_TCdata.m')


