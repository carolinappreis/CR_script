clear all
all=[2 3 4 5 8 10 11 13 16];
% iiii=[2 5 8];

iiii=all(1:9);

for numb=1:length(iiii);
    clearvars -except iiii numb
    % load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\PLS\p0',num2str(iiii(numb)),'_PLS.mat'))
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_stim/RS/P0',num2str(iiii(numb)),'_RS.mat'))


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
f1=figure;
bar(0:30:330,100.*nanmedian(tt)) 
xlabel('Stimulation phase (degrees)')
ylabel('Change in tremor severity (m/s^2)')
box('off')
set(f1,'color','w');
set(gca,'FontSize',14)
f1.Units = 'centimeters';
set(gca,'XTickLabelRotation',45)
f1.OuterPosition= [10, 10, 12, 12];
cd('/Users/Carolina/OneDrive - Nexus365/arcs_share/peripheral')
saveas(f1,['pheritremor_',num2str(iiii(numb)),'.png'])
close all
end


% hold on
% plot(0:30:330,100.*tt,'.')
% % x axis phase - y axis percent change in 
% % tremor severity at the end of 5 seconds with respect to severity right
% % before stimulation began


