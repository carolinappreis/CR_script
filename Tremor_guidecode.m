 cd('/Users/Carolina/Documents/GitHub/CR_script/Tremor')
cd('C:\Users\creis\Documents\GitHub\CR_script\Tremor')
cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')

%Phase tracking within an across stim runs (per phase)
run('phase_tracking.m')

%Mat for Amplitude Response curves and Random Stim dist
run ('amp_arc.m') % will run 'start_cleaner.m'
%to get amp_ARC.mat

% Mat for Amplitude Response Curves and no stim surrogate
run ('non_stim_ampsh.m')

%Mat for Phase Response curves and Random Stim dist
run ('freq_frc.m') % will run 'start_cleaner.m'

%Mat for PRC surrogates
run('non_stim_phasesh.m')

% Group stats - realignment
run('group_realign.m')

%Stim-evoked potential on sig.phase effect phase?
%%%% needs updated ARC's if used again
run('stim_evoked.m') % will run ('pre_stim_evoked.m')

%Check contribution of 3 axis to PCA
run('pc_arc.m')
%saves am_ax.mat - arc with runs where pre-defined main axis was the main
%contributor to pcr
% saves s_arc.mat and t_arc.mat - arc ran of th eother 2 axis

run('pc_frc.m')
% saves same as pc_arc but for frc

%Smoothed response curves
run ('smooth_RC.m')
% saves smooth_arc3.mat and smooth_frc3.mat  

%plot arc prc with thresholds and smooth curves with extra axis response
%curves
run ('arc_prc_plots.m')

% GROUP LEVEL response curves for all data, median split a>median &
% a<median
run('group_rc_smooth_before.m')

% correlations between change in amplitude an absolute A over trials (beg
% test :end test)
run('trial_ampchange.m')

%corr amplitude before and after stim
run('corr_amp_before_after_stim.m') % to be changed if used

%median split arc
run('median_split_a.m')
%saves arc_mediansplit.mat

%median split frc
run('median_split_f.m')
%saves frc_mediansplit.mat


%flow of analysis
%determine main axis
%find peak frequency
%find start and end points of stimulation fitiing the criteria
%find changes in amp and freq 
%exclude changes in amp and freq that are outliers (????)
% look at pc coefficents

%-------------------------------------------------------------------------%


%%% original codes from HC
run('peripheraltremor.m') % will run phasedetection.m
run ('nostim.m')
run('phasedetection_general_phase_abby.m')


%%%  OLD VERSIONS

% sig phase-specific stimulation effects 
run ('periph_tremor_CR.m') %will run ('start_clean.m')

% sig effect of stim compared to non stim periods
run ('tremor_code_HC.m')

%plots from the above
run('plots') % change for NS/PS for amplitude effects ot NS_PHS1 and PS_PHS1 for effects on frequency

%phase-specific change on Frequency induced by stim
run ('stim_phasesh.m')

%sig. induced phase-shift vs.baseline phase-shift
run('non_stim_phasesh.m')

%Correlation of Amplitude with change in Amp and Freq
run('trial_corr.m') %% dont use! y axis includes abs values in x axis
run('corr_change_abs.m') %% dont use! y axis includes abs values in x axis




