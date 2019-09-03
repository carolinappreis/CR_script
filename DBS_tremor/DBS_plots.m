clear all
% cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data')
cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\')

load ('DBS_amp_ARC.mat')

cond={'NS_PS_result.mat';'HFS_PS_result.mat'};

k=2;
n=2;

load(cond{k,1})
% load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash')
load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');

cl=blushred;



bar(ttall{n,1},'FaceColor',cl,'EdgeColor',cl)
hold on
yline(prctile(NS(1,n,:),99.7917),'k--','LineWidth',1)
yline(prctile(NS(1,n,:),0.2083),'k--','LineWidth',1)
box ('off')

if k==1 && n==1
    title('Posture phasic stim vs. no stim')
    
elseif k==1 && n==2
        title ('Spiral phasic stim vs. no stim')
        
elseif k==2 && n==1
    title('Posture phasic stim vs. Continuous stim')
    
elseif k==2 && n==2
    title('Spiral phasic stim vs. Continuous stim')
end
    


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
    