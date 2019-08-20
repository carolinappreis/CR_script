
% cd('/Users/Carolina/Documents/GitHub/CR_script/Tremor')
 cd('C:\Users\creis\Documents\GitHub\CR_script\Tremor')
cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')


%Phase tracking within an across stim runs (per phase)
run('phase_tracking.m')

%Mat for Amplitude Response curves and Random Stim dist
run ('amp_arc.m') % will run 'start_cleaner.m'

% Mat for Amplitude Response Curves and no stim surrogate
run ('non_stim_ampsh.m')

%Mat for Phase Response curves and Random Stim dist
run ('freq_frc.m') % will run 'start_cleaner.m'

%Mat for PRC surrogates
run('non_stim_phasesh.m')

%Correlation of Amplitude with change in Amp and Freq
run('trial_corr.m')
run('corr_change_abs.m')
% Group stats - realignment
run('group_realign.m')

%Stim-evoked potential on sig.phase effect phase?
run('stim_evoked.m') % will run ('pre_stim_evoked.m')

%Check contribution of 3 axis to PCA
run('pc_arc.m')
run('pc_frc.m')

%Smoothed response curves
run ('smooth_RC.m')

% GROUP LEVEL response curves
run('group_rc_smooth_before.m')

% Compare amplitude begining/end experiment
run('compare_amp2.m')

%corr amplitude before and after stim
run('corr_amp_before_after_stim.m')

%median split arc
run('median_split_a.m')

%median split frc
run('median_split_f.m')



%flow of analysis
%determine main axis
%find peak frequency
%find start and end points of stimulation fitiing the criteria
%find changes in amp and freq 
%exclude changes in amp and freq that are outliers

% look at pc coefficents
% run 



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




