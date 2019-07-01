
% cd('/Users/Carolina/Documents/GitHub/CR_script/Tremor')
cd('C:\Users\creis\Documents\GitHub\CR_script\Tremor')
% sig phase-specific stimulation effects 
run ('periph_tremor_CR.m') %will run ('start_clean.m')
% sig effect of stim compared to non stim periods
run ('tremor_code_HC.m')
%plots from the above
run('plots')

%phase-specific change in stim
run ('priph_tremor_phaseshift.m')

