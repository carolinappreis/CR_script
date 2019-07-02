% cd('C:\Users\creis\Documents\GitHub\CR_script\Tremor')
clear all
 load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\NS\p01_stimnosim')
 load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\PS\p01_pha_suffle.mat')
d=nanmedian(tt);

upperthreshold=prctile(nostimout,99.7917); %bonferroni corrected for 12 comparisons
lowerthreshold=prctile(nostimout,0.2083);

% find(d>=upperthreshold | d<=lowerthreshold)

close all
fig=figure(1)
subplot(1,3,1)
bar(d,'FaceColor',[0.5 0.5 0.5],'EdgeColor', [0.5 0.5 0.5]) 
hold on 
rr(1:12)=upperthreshold;
rr2(1:12)=lowerthreshold;
plot(rr,'r--','LineWidth',1)
plot(rr2,'r--','LineWidth',1)
xlim([0.5 12.5])
box('off')
xticklabels({'0' '30' '60' '90' '120' '150' '180' '210' '240' '270' '300' '330'})
ylabel ('Median amplitude change')
xlabel ('Stimulated phase')
ylim ([-1 1])
title 'Significant stimulation effect'
set(gca,'FontSize',12)

rr(1:size(tt,2))=mean(prctile(tt3,99.7917));
rr1(1:size(tt,2))=mean(prctile(tt3,0.2083));

subplot(1,3,2)
bar(nanmedian(tt),'FaceColor',[0.5 0.5 0.5],'EdgeColor', [0.5 0.5 0.5])
hold on
plot(rr,'r--','LineWidth',1)
plot(rr1,'r--','LineWidth',1)
box('off')
title 'Significant phasic-specific stimulation effect' 
xticklabels({'0' '30' '60' '90' '120' '150' '180' '210' '240' '270' '300' '330'})
ylabel ('Median amplitude change')
xlabel ('Stimulated phase')
ylim ([-1 1])
 set(gca,'FontSize',12)

subplot(1,3,3)
likhood_amp=sum(tt>prctile(tt3,97.5) & tt>0)./sum(~isnan(tt));
likhood_sup=sum(tt<prctile(tt3,2.5) & tt<0)./sum(~isnan(tt));
likhood=[likhood_sup ; likhood_amp];
bar(likhood')
title ('Likelihood of significant amplification/supression effect')
xticklabels({'0' '30' '60' '90' '120' '150' '180' '210' '240' '270' '300' '330'})
legend('supression','amplification')
legend('boxoff')
box('off')
xlabel ('Stimulated phase')
ylabel ('Likelihood') 
ylim ([0 1])


fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 50, 15];
set(gca,'FontSize',12)


set(fig,'color','w');


figure()
plot(C)
hold on
plot(AA,C(AA),'r.')
plot(BB,C(BB),'b.')
set(gcf,'color','w');
title ('Baseline')
xlabel('Time (msec)')
ylabel ('Tremor Envelope')
box('off')

cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim')
load ('P011_randstim_cursos.mat')
in2=1; % analysing the "main tremor axis"
if in2==1
    in=3;
elseif in2==2 % other axis 1
    in=6;
elseif in2==3 % other axis 2
    in=7;
end
data=SmrData.WvData;
samplerateold=SmrData.SR;
tremor=(data(in,:));
addon=92; addon_end=35;

phasedetection;

data=SmrData.WvData;
rep=10; % number of trials for random stim - please enter for each patient
clearvars -except Fpeak in2 in rep SmrData data stimout

samplerateold=SmrData.SR;
time=0:1/samplerateold:(size(data,2)-1)/samplerateold;
tremor=data(in,:);
ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremor2(1:size(ts1.data,3))=ts1.data;
samplerate=1000;

if (Fpeak-2)>=1
    [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
    [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end

tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
figure()
plot(tremor_or)
set(gcf,'color','w');
title ('Tremor during random phase stim')
xlabel('Time(msec)')
ylabel ('Filtered signal main axis')
box('off')