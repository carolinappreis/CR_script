%TC PROJECT

%-----------------------BUA-A1 Analysis----------------------------------------
cd ('C:\Users\creis\Documents\GitHub\CR_script\A1_Thal\code')

%Organizes data acording to lesion and non-lesioned rats
run 'sort_lesion_nonlesion.m'

%Organizes lesioned data acording to thalamic regions (all probes)and takes
%the median change in beta activity in both conditions: short and long
%bursts.
run 'filgen_noncoh_new.m'
    % will run 'sort_data.m' to save files to be ploted - ATENTION,match
    % file names.
    %check which region has significant statistical difference with
    %surrogates: run 'testing_clusters_noncoh.m'
run 'sublplot_short.m'%plots
run 'subplot_long.m'%plots

%%Organizes lesioned data acording to thalamic regions (ONLY BETA COHERENT probes)and takes
%the median change in beta activity in both conditions: short and long
%bursts.
run 'A1_coherent_files.m' %saves changes (BZ_complete_new and SNr_complete_new)  in beta activity aligned to short and long bursts 
    % run sort_data.m (uncomment BZ_complete_new and SNr_complete_new) to save
    % files in data format
    
run 'coherent_regions.m'%plots
    
%-----------------------BUA-A2 Analysis----------------------------------------
% cd ('C:\Users\creis\Documents\GitHub\CR_script\A2_Thal\code')
cd('/Users/Carolina/Documents/GitHub/CR_script/A2_Thal/code')

%Organizes data acording to lesion and non-lesioned rats
run 'sort_lesion_nonlesion.m'

%Selects region of interest; generates BUA; filters EEG leaves probe
%recordings unfitered and aligned to bursts (all burts; no duration discrimination)
run 'thal_all.m'
    %will run 'phase_onset_new.m' % alignes unfiltered subcortical signal
    %to phase-aligned burst onset
    
    %will run 'new_thal.m' %organizes and saves data to be plotted: summed beta
    %evoked; filtered/z-scored/enveloped subcortical signal to get phase-consistency
    
%plot beta evoked events and corresponding envelopes against surrogates
run 'A2_2_patches.m'

%obtaines the point in time where beta envelope rises subcortically
run 'envelope_rise_evoke_cr.m'

%stat. difference between starting point and max amplitude of evoked betaof evoked beta between regions 
run %'?'

%-----------------------BUA-A3 Analysis----------------------------------------
cd ('C:\Users\creis\Documents\GitHub\CR_script\A3_Thal\code')
% cd('/Users/Carolina/Documents/GitHub/CR_script/A3_Thal/code') 

% psi across regions in different frequencies
run 'psi_across_frequencies.m'

%phase synchrony index(psi)between ctx and different subcortical regions
%across animals - mat files from thal_all (A2)
run 'psi_BUA_ctxsub_new.m' % run 'psi_BUA_ctxsub_msc.m' - crashing

%phase synchrony index(psi)from between regions and in the same rat
run 'psi_BUA_subsub.m'

%phase synchrony index(psi)from different regions and rats with the
%cortical &/or subcortical signal in bursts vs surrogates
run 'psi_bursts_msc.m'

%phase mean in bursts (onset-median duration) vs before (200msec)
run 'mang_burst_ctxsub.m'

%-----------------------BUA-A4 Analysis----------------------------------------
%  cd ('C:\Users\creis\Documents\GitHub\CR_script\A4_Thal\code')
cd('/Users/Carolina/Documents/GitHub/CR_script/A4_Thal/code') 


%phase synchrony index(psi)consistency during ctx/subctx bursts vs surrogates between ctx/subcortical regions
run 'PSI_betactxburst_5fctxsub.m'

%phase synchrony index(psi)consistency during ctx bursts vs surrogates between subcortical-subcortical regions
run 'subcortical2_same_rat'

%Rats with 10% coherence in beta with cortex in 15-30Hz and sig. increase
%in psi during bursts
run 'selected_rats' % data comes from 'psi_rats_savemat'

%Regions' psi together 
run 'region_PSI'

%BZ and SNr in high and low beta
run 'same_region_dif_freqs.m'

%-----------------------SUA Juxta_SUA_act  Analysis------------------------------
  cd ('C:\Users\creis\Documents\GitHub\CR_script\SUA\Juxta SUA_act_mat')
%  cd('/Users/Carolina/Documents/GitHub/CR_script/SUA/Juxta SUA_act_mat') 

%phase consistency of subcortical firing of subcortical areas
run 'psi_SUA.m'

% Single unit firing rate and phase consistency of subcortical firing in
% bursts vs surrogates && ISI, firing rate per region 
run 'Filegen_SUA_act_new.m'

% Sum fo spiking across bursts in time.
run 'gaussian2.m'
%-----------------------SUA Probe_SUA_act  Analysis------------------------------
  cd ('C:\Users\creis\Documents\GitHub\CR_script\SUA\Probe SUA_act_mat')
%  cd('/Users/Carolina/Documents/GitHub/CR_script/SUA/Probe SUA_act_mat') 

 % Select rats at region and that have shown coh with ctx in the bua
 % analysis/ were chosen by KN in his excel file
 run 'filegen_slected_rats.m'
 
 %Sum fo spiking across bursts in time.
 run 'gaussian_probe.m'



