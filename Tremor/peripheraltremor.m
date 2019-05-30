% SMR_File_To_Mat;

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

% bar(0:30:330,100.*nanmedian(tt)) % x axis phase - y axis percent change in 
% % tremor severity at the end of 5 seconds with respect to severity right
% % before stimulation began
% box('off')

close all
figure(2)
fig=gcf;
fig.Color=[1 1 1];
bar(100.*nanmedian(tt))
hold on
stem((100.*tt)')
xticklabels({'0','30','60','90','120','150','180','210','240','270','300','330'})
ylim([(-max((max(tt)))-0.1).*100 (max(max(tt))+0.1).*100])
box('off')
title ('P11')


