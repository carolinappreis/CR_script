%TC PROJECT v.2
%--------------------- CTSNR-----------------------------------------------
cd('/Users/Carolina/Documents/GitHub/CR_script/CTSNR')

%%% BUA
run 'new_mastercode.m' 
%amplitude change
        %%% subcortex ref
        run 'sub_to_ctx.m'
        %%% ref cortex
        run 'ctx_to_sub.m'
%psi_change
        %%% subcortex ref
        run 'sub_to_ctx_psi.m'
        %%% ref cortex
        run 'ctx_to_sub_psi.m'
        
%phase slips
run 'pha_slip_ctx_subctx.m'
run 'pha_slip_ctx.m'
run 'pha_slip_sub.m'

%Overlaps
    %%% subcortex ref
    run 'subctx_ov.m'
    %%% ref cortex
    run 'ecog_sub_ovl.m'

%%% SUA/JUXTA
%triggerd-averages
run('PrepData_trg_avg.m') % will call fx trg_avg.m
        
%Modulation of angle and vector length cycle by cycle       
run('PrepData_cy_by_cy.m') % will call fx cy_by_cy.m 
        %uses output from for SUA 'pl_nonunif_sua.m' with data_SUA_SNR.mat
        %selects units non-uniformily locked to
        %one of three sub beta bands - the one with highest vector length.

%-----------------------BUA-A4 Analysis----------------------------------------
  cd('/Users/Carolina/Documents/GitHub/CR_script/CTSNR/OLD/A4_Thal/code')

%phase synchrony index(psi)consistency during ctx/subctx bursts vs surrogates between ctx/subcortical regions
run 'combine2.m'
run 'exploring_psi.m'
run 'plots_psi_time_bursts.m'
% imporatnt functions are bursts.m burts_off.m bursts_aligned.m hayriye_c
% bursts_aligned_off.m

%phase-slips analysis
run 'phase_slip_last.m' % up to date 
% run 'phase_slips.m' % as described in PNAS - Cagnan et al.2019
% run 'during_severalb.m' % according to bins of busrt duration
run 'idregion_phaseslip.m' % change of instantaneous frequency at ctx and thal independantly
%PSD ctx and subcortical regions
run 'psd_plots'

    
%----------------------- Juxta  Analysis------------------------------
cd('/Users/Carolina/Documents/GitHub/CR_script/CTSNR/OLD/SUA/Juxta SUA_act_mat') 

%phase consistency of subcortical firing of subcortical areas
run 'psi_SUA.m'

% Single unit firing rate and phase consistency of subcortical firing in
% bursts vs surrogates && ISI, firing rate per region 
run 'Filegen_SUA_act_new.m'

% Sum fo spiking across bursts in time.
run 'gaussian2.m'

%matching sua triggered avgs
run 'trig_avg_jux_1b' 
run 'trig_avg_jux_4b' 

%Cycle by cycle analysis
run 'CY_CY.m' % runs function 'cycles_10.m'

% cy_cy 4 betas
run 'cy_cy_j_4b.m'

%cy_cy 1 beta
run 'cy_cy_j_1b.m'


%histogram of firing rate and phase locking of units with ECoG
run 'frate_phaselocking.m'


%-----------------------SUA Probe_SUA_act  Analysis------------------------------

cd('/Users/Carolina/Documents/GitHub/CR_script/CTSNR/OLD/SUA/Probe SUA_act_mat') 


 % Select rats at region and that have shown coh with ctx in the bua
 % analysis/ were chosen by KN in his excel file
 run 'filgen_slected_rats.m'
 
 %Sum fo spiking across bursts in time.
 run 'gaussian_probe.m'
 
 % all units in BZ/snr (not filtered by firing rate)
 run 'new_triggered_sua.m' 
 
 %correlation spikes snr bz same rat
 run 'correlation_sua.m'
 
%stats on phase alignement (takes forward units non-uniformily locked to
%one of three sub beta bands - the one with highest vector length.
run 'pl_nonunif_sua.m'  % uses data_SUA_SNR.mat

run 'trig_avg_pl_nonunif_4bands.m'% uses outout ['NEW_BZ_cycle.mat'], from pl_nonunif_sua.m

%%% trig_avg beta locked, non-uni 
run 'trig_avg_pl_nu_1band.m' 

%OLD versions cycle by cycle 
% run 'CY_CY_AS.m' ; run 'CY_CY_HC.m' 

%%%cycle by cycle using ['NEW_BZ_cycle.mat'], from pl_nonunif_sua.m 
run 'cy_cy_4betas.m'

%%%cycle by cycle using data_SUA_SNR.mat and chosing phase-locked units to
%%%1 beta non-uniformly
run 'cy_cy_1beta.m'


%phase locking of units to ECoG and firing rate
run 'phase_consist.m'

%evolurion of vector length cy cy
run 'vec_length_timecourse.m'

%-----------------------MUA------------------------------------------------------
cd('C:\Users\creis\Documents\GitHub\CR_script\MUA')
cd('/Users/Carolina/Documents/GitHub/CR_script/CTSNR/OLD/MUA')

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
 cd('/Users/Carolina/Documents/GitHub/CR_script/CTSNR/OLD/A4_Thal/code') 
 run('psd_plots.m')
 run('coherence_plots.m')
 
%amplitude change, PSI 
cd('C:\Users\creis\Documents\GitHub\CR_script')
cd('/Users/Carolina/Documents/GitHub/CTSNR/OLD/CR_script')
run('plots_TCdata.m')

%%%---------------------- matfiles
% BZ_bua -- raw data each cell is a probe: 1st array ctx, 2:end probe contacts
% BZ_opt -- all metrics (filtered, phase, onset , coherence)
