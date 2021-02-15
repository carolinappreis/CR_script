
clear
close all
% cd('/Users/Carolina/Desktop')
cd('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/SH_spirals_xcels')
NS1=xlsread('P03_HF1_SH.xlsx');

plot(NS1(:,1),NS1(:,2))
close
NS1(find(NS1(:,1)==-1),:)=[];
plot(NS1(:,1),NS1(:,2),'Color',[0.8 0.5 0.8],'LineWidth',2)
%  plot(NS1(:,1),NS1(:,2),'Color',[0 0.5 0.5],'LineWidth',2)
% hold on
% plot(NS1(:,1),NS1(:,2),'r.')
box('off')
xticks([])
yticks([])
ax = gca
ax.Visible = 'off'

cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA')
