clear all
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A2_Thal/MAT')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A2_Thal\mat')

load 'animal_lesion_nolesion.mat'
'kjx127z01@170-270_m.mat          '

'kjx132f02@180-280_m.mat          '
    
 'kjx136b01@0-100_m.mat            '
 'kjx140d01@0-100_m.mat            '

 'kjx140f01@70-170_m.mat           '
'kjx160c01@0-100_m.mat            '
    
'kjx166b03@0-100_m.mat            '

'kjx167c01@200-300_m.mat          '


T=[7 9 10 15 17 18 22 27];
S=[1:6 8 11:14 16 19:21 23:26 28:29];

lesion_nonrep=lesion(:,S);

clearvars -except lesion_nonrep lesion A nolesion newfile
% save 'animal_lesion_nolesion.mat'