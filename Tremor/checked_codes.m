%%checked_codes


cd('C:\Users\creis\Documents\GitHub\CR_script\Tremor')
cd('/Users/Carolina/Documents/GitHub/CR_script/Tremor')


run('Periph_PhasicStim_smr.m') 

run('peripheraltremor_newcode.m')

run('NS_sigmoid.m')

run('pls_sig.m')
% runs 'code_before_PLScell.m' 

run('NS_new') %%%% HC to generate .mat baseline do non stim 

run('rc_axis_ns_thresh')%%%% HC To generate .mat file with arcs 

run('median_split_a') 

run('all_plots') %%%% to plot the curves of the pateints of interest
% arcs from rc_axis_ns_threshold; thresholds from NS_new; median split from medial_split_a

run('plot_pls_ns') % absolute amplitude bef posture, before stim, end stim for ns and pls - data from NS_sigmoid and pls_sig
    run('plot_pls_ns_demo') % diagram of envelope



    GROUP
    
run('group_arc_plot') %% old code with freq names group_rc_smooth_before

