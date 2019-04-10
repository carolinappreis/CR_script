clear all
clear all
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA:act:mat')
load ('BZ_onerat_inburst.mat')

[b,a]=butter(2,[15 36]./(0.5*samprate));
Ecogfiltered=filtfilt(b,a,WaveData_DC(1,:));

ang_data=angle(hilbert(WaveData_DC));
ang_trials=angle(hilbert(output));
ang_diff=ang_data(2,:)-ang_data(8,:);


abs(mean(exp(1i*(ang_data(5,:)-ang_data(6,:)))))
polar(repmat((ang_data(5,:)-ang_data(6,:)),1,2)',repmat([0 1],1, 100000)');
polarhistogram(ang_data(5,:)-ang_data(6,:),12);
histogram(ang_data(5,:)-ang_data(6,:),12);



ispc_trials=squeeze(abs(mean(exp(1i*diff(ang_trials,[],1)),2)));
plot(1:1:size(output,3),ispc_trials)

ispc_time= squeeze(abs(mean(exp(1i*diff(ang_trials,[],1)),3)));
plot(1:1:size(output,2),ispc_time)

for i =1:115
subplot(2,1,1)
polarhistogram(ang_trials(1,i,2000:2200),12);
subplot(2,1,2)
polarhistogram(ang_trials(1,i,1800:2000),12);
end

% 
% for i =1:115
% subplot(2,1,1)
% polarhistogram(ang_trials(1,i,1800:2000)- ang_trials(7,i,1800:2000),12);
% subplot(2,1,2)
% polarhistogram(ang_trials(1,i,2000:2200)- ang_trials(7,i,2000:2200),12);
% end

% 
cd ('/Users/Carolina/Documents/MATLAB/SUA/Juxta SUA:act:mat')
run Filegen_SUA_act

for i=1:7
polarhistogram(angles_inb{i,1},12);
end

% ang_sig1=abs(angle(hilbert(beta_allctx(1,:))));
% ang_sig2=abs(angle(hilbert(beta_allctx(2,:))));
% dif_angles=ang_sig1-ang_sig2;
% abs(mean(exp(1i*(dif_angles))))
% polar(repmat(dif_angles,1,2)',repmat([0 1],1, 100000)');
    