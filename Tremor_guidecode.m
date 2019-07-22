
cd('/Users/Carolina/Documents/GitHub/CR_script/Tremor')
% cd('C:\Users\creis\Documents\GitHub\CR_script\Tremor')

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




