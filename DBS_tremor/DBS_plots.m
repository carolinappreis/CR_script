clear all
close all
%  cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data')
 cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\')

ax={'DBS_amp_ARC.mat','s_arc.mat','t_arc.mat'};
cond={'NS_PS_result.mat','NS_ax2.mat','NS_ax3.mat'};


k=1; % no stim a1, a2, a3) 
n=2;%posture (1) vs spiral (2)

load(cond{1,k})
load(ax{1,k})


%  load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash')
 load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','aegean');

cl(1,:)=blushred;
cl(2,:)=aegean;



  bar(ttall{n,1},'FaceColor',cl(n,:),'EdgeColor',cl(n,:))
% bar(squeeze(t_arc(n,:,:)),'FaceColor',cl(n,:),'EdgeColor',cl(n,:))
hold on
% yline(prctile(idv_NS(1,n,:),99.7917),'k--','LineWidth',1)
% yline(prctile(idv_NS(1,n,:),0.2083),'k--','LineWidth',1)
yline(prctile(idv_NS(n,1,:),99.7917),'k--','LineWidth',1)
yline(prctile(idv_NS(n,1,:),0.2083),'k--','LineWidth',1)
box ('off')


% clear all
% cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data')
% 
% load ('DBS_amp_ARC.mat')
% 
% n=1;
% 
% load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash')
% cl=blushred;
% 
% corrls=cell2mat(LS);
% 
% LS1{1,1}=corrls(1,(~isnan(corrls(1,:))));
% LS1{2,1}=corrls(2,(~isnan(corrls(2,:))));
% 
% 
% bar(ttall{n,1},'FaceColor',cl,'EdgeColor',cl)
% hold on
% yline(prctile(LS1{n,1}(1,:),99.7917),'k--','LineWidth',1)
% yline(prctile(LS1{n,1}(1,:),0.2083),'k--','LineWidth',1)
% box ('off')
% 
% if  n==1
%     title('Posture phasic stim vs. Random Stim')
%     
% elseif n==2
%         title ('Spiral phasic stim vs. Random Stim')
%         
% 
% end
    