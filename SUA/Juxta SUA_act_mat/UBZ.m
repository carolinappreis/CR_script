clear all
 cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
load ('SUA_BZ')

freq=17.5;
[b,a]=butter(2,[(freq-5)/(0.5*srn) (freq+5)/(0.5*srn)],'bandpass');

for i =1:size(WaveData_DC,1)
    UBZ.psd(i,:)=pwelch(WaveData_DC(i,:),1000,[],1000,1000);
    UBZ.ecog_filt(i,:)=filtfilt(b,a,WaveData_DC(i,:));
    UBZ.ecog_env(i,:)=abs(hilbert(UBZ.ecog_filt(i,:)));
    env=UBZ.ecog_env(i,:);Ecogfiltered=UBZ.ecog_env(i,:);
    UBZ.onset_raw{1,i}=bursts(env);
%     UBZ.onset_phase_al{1,i}=bursts_aligned(env,Ecogfiltered);
    UBZ.ecog_phase(i,:)=angle(hilbert(UBZ.ecog_filt(i,:)));
    clearvars env Ecogfiltered
end

% plot(log(UBZ.psd(:,1:50))')

 