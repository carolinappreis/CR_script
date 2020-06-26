cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data')

cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA')

cd('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA')
%_________________________________________________________________________%
cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/')
cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_posture')
cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_spiral')
cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/old_DBS')


% ARC for posture and spiral conditions 
run 'DBS_arc.m' %will run DBS_cleaner
%saves DBS_amp_ARC.mat
% plots smooth ARC and sin/gaussian fits

% FRC for posture and spiral conditions 
run 'DBS_frc.m' %will run DBS_cleaner
%saves DBS_freq_FRC.mat
% plots smooth FRC and sin/gaussian fits

% ARC for posture and spiral conditions median split of A
run ('DBS_RS_PSmediansplit.m')
%saves arc_mediansplit.mat

% Finding onset points for motor tasks in NS and HFS recordings
run('DBS_3Dposition.m')

% Surrogates of NS HF to deliniate sig. stim - fidley! 
run('DBS_NS_FS_analysis.m')



cd('C:\Users\creis\Documents\GitHub\CR_script\DBS_tremor\hmm')
cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/hmm')