
clear all
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\amp_ARC.mat')
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\amplifying_phases')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/abs_amp.mat')

figure()
subplot(1,2,1)
m=[mean(ampall,2) mean(ttall,2)];
plot(m(:,1),m(:,2),'k+')
y1=lsline
box('off')
ylabel ('Change in amplitude (zscore)')
xlabel('Amplitude (zscore)')

c=corrcoef(m(:,1),m(:,2));
legend(y1,[num2str(c(1,2))],'box','off')

for i=1:9
a1(1,i)=ttall(i,phase_peak(i));
a2(1,i)=ampall(i,phase_peak(i));
end

subplot(1,2,2)
plot(a1,a2,'k+')
y2=lsline
box('off')
ylabel ('Change in amplitude (zscore)')
xlabel('Amplitude (zscore)')
c2=corrcoef(a1',a2')
legend(y2,[num2str(c2(1,2))],'box','off')




clear all
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\freq_FRC.mat')
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\amplifying_phases')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/abs_freq.mat')

figure()
subplot(1,2,1)
m=[mean(freqall,2) mean(ttall,2)];
plot(m(:,1),m(:,2),'k+')
y1=lsline
box('off')
ylabel ('Change in frequency (zscore)')
xlabel('Frequency (zscore)')
c=corrcoef(m(:,1),m(:,2));
legend(y1,[num2str(c(1,2))],'box','off')

for i=1:9
a1(1,i)=ttall(i,phase_peak(i));
a2(1,i)=freqall(i,phase_peak(i));
end

subplot(1,2,2)
plot(a1,a2,'k+')
y2=lsline
box('off')
ylabel ('Change in frequency (zscore)')
xlabel('Frequency (zscore)')
c2=corrcoef(a1',a2')
legend(y2,[num2str(c2(1,2))],'box','off')

