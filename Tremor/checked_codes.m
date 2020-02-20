%%checked_codes


cd('C:\Users\creis\Documents\GitHub\CR_script\Tremor')
cd('/Users/Carolina/Documents/GitHub/CR_script/Tremor')


run('Periph_PhasicStim_smr.m') 

run('peripheraltremor_newcode.m')

run('NS_sigmoid.m')

run('only_sig_clean_pls.m') % runs 'code_before_PLScell.m' 
% old version run('pls_sig.m')


run('NS_matched') %%%% HC to generate .mat baseline do non stim 
run('NS_mediansplit') %%%% HC to generate .mat baseline do non stim 

run('rc_axis_ns_thresh')%%%% HC To generate .mat file with arcs 

run('median_split_a') 

run('all_plots') %%%% to plot the curves of the pateints of interest
% arcs from rc_axis_ns_threshold; thresholds from NS_new; median split from medial_split_a

run('plot_pls_ns') % absolute amplitude bef posture, before stim, end stim for ns and pls -

% data from only_sig_clean_pls & only_sig_ns OLD verison: data from NS_sigmoid and pls_sig
run('plot_pls_ns_demo') % diagram of envelope

run('gauss_sin_fit.m') % fits using fit funtion instead of fitnlm
run('nlm_fits_fx.m')  %to get fits from sine and gaussian models to raw response curves - are ARCs non-uniform?
run('nonlin_lin_AIC.m') %to run gauss sin and k fits and do moel comparison with AIC 

%%%--------------------------GROUP
    
run('group_arc_plot') %% old code with freq names group_rc_smooth_before
run ('group_statistics.m') %% stats same peak stim non stim
run('group_statistics_eachpeak.m') %% stats each main peak 

run('rampup_stim.m') % sigmoid fits to hand up stim trials (10 trials all partients) + ttest fromo stim/ns
run('rampup_NS.m') % sigmoid fits to hand up 
run('rampup_onemin_PLS.m')

run('pwelch_amp_sup_ns.m') %% will run file taken from 'NS_zscore_handuptimedomain.m'
run('pwelch_3plot.m')
%%%phase-locked stim

run('pls_group_sup.m') %%% plotting median tremor severity during sup 
run('pls_group_amp.m') %%% plotting median tremor severity during amp


%%%%----------------------CLEANING DATA
run ('finding_stim_thre.m') % 3 pateints with thriggering issues - tremor characteristics found. 