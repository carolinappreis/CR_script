%%%% Create smooth data and save

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------------AMP
% %%% loading all ARC axis
% clear all
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\A_group')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\f_ax')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Axis2\s_arc')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Axis3\t_arc')
% 
% % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/A_group')
% % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/f_ax')
% % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Axis2/s_arc')
% % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Axis3/t_arc')
% 
% 
% SO=[S S];
% S1=[f_ax f_ax];
% S2=[s_arc s_arc];
% S3=[t_arc t_arc];
% 
% 
% for ii=1:size(SO,1)
%     dum=[];
%     dum1=[];
%     dum2=[];
%     dum3=[];
%     for i=1:size(SO,2)
%         if i+2<size(SO,2)
%        dum=[dum sum(SO(ii,(i:i+2)))./length(SO(ii,(i:i+2)))];
%        dum1=[dum1 sum(S1(ii,(i:i+2)))./length(SO(ii,(i:i+2)))];
%        dum2=[dum2 sum(S2(ii,(i:i+2)))./length(SO(ii,(i:i+2)))];
%        dum3=[dum3 sum(S3(ii,(i:i+2)))./length(SO(ii,(i:i+2)))];
%         end
%     end
%     sm_s(ii,:)=dum; 
%     sm_s1(ii,:)=dum1;
%     sm_s2(ii,:)=dum2;
%     sm_s3(ii,:)=dum3;
% end
% new=[1:5 7:9]; %%without PD patient (i.e., pt number 6)
% 
% smo_s=sm_s(new,1:12);
% smo_s1=sm_s1(new,1:12);
% smo_s2=sm_s2(new,1:12);
% smo_s3=sm_s3(new,1:12);
% 
% clearvars -except smo_s smo_s1 smo_s2 smo_s3
% % cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
% % cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')

% % save('smooth_arc3.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------------FREQ

% % %%% loading all FRC axis
% %  clear all
% % load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\F_group1')
% % load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\f_frc')
% % load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Axis2\s_frc')
% % load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Axis3\t_frc')
% % 
% % % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/F_group1')
% % % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/f_frc')
% % % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Axis2/s_frc')
% % % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Axis3/t_frc')
% % 
% % 
% % SO=[S S];
% % S1=[f_frc f_frc];
% % S2=[s_frc s_frc];
% % S3=[t_frc t_frc];
% % 
% % 
% % for ii=1:size(SO,1)
% %     dum=[];
% %     dum1=[];
% %     dum2=[];
% %     dum3=[];
% %     for i=1:size(SO,2)
% %         if i+2<size(SO,2)
% %        dum=[dum sum(SO(ii,(i:i+2)))./length(SO(ii,(i:i+2)))];
% %        dum1=[dum1 sum(S1(ii,(i:i+2)))./length(SO(ii,(i:i+2)))];
% %        dum2=[dum2 sum(S2(ii,(i:i+2)))./length(SO(ii,(i:i+2)))];
% %        dum3=[dum3 sum(S3(ii,(i:i+2)))./length(SO(ii,(i:i+2)))];
% %         end
% %     end
% %     sm_s(ii,:)=dum; 
% %     sm_s1(ii,:)=dum1;
% %     sm_s2(ii,:)=dum2;
% %     sm_s3(ii,:)=dum3;
% % end
% % % new=[1:5 7:9]; %%without PD patient (i.e., pt number 6)
% % 
% % smo_s=sm_s(:,1:12);
% % smo_s1=sm_s1(:,1:12);
% % smo_s2=sm_s2(:,1:12);
% % smo_s3=sm_s3(:,1:12);
% % 
% % clearvars -except smo_s smo_s1 smo_s2 smo_s3
% % % cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
% % cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
% % save('smooth_frc3.mat')



%%%% use smoothed data

clear all
close all
metric=0;

if metric==0;
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/smooth_arc3')
% load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash')

load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\A_group')
load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\smooth_arc3')   
load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
new=[1:5 7:9]; %%without PD patient (i.e., pt number 6)
S=S(new,:);

cl=blushred;
cl1=squash;

else
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/smooth_arc3')
% load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','aegean','stone');

load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\F_group1')
load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\smooth_frc3')    
load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','aegean','stone');

cl=aegean;
cl1=stone;
end


for i=1:size(smo_s,1)
    f1=figure(1)
    subplot(1,8,i)
    bar(smo_s(i,:),'LineWidth',0.5,'FaceColor',cl,'EdgeColor',cl)
    hold on
     bar(smo_s1(i,:),'LineStyle','--','LineWidth',1,'FaceColor','none','EdgeColor',cl1)
    yline(0,'LineWidth',1)
    box('off')
    
    f2=figure(2)
    subplot(1,8,i)
    bar(smo_s(i,:),'FaceColor',cl,'EdgeColor',cl)
    hold on
    plot(smo_s2(i,:),'LineWidth',1,'Color',cl1)
    plot(smo_s3(i,:),'LineWidth',1,'Color',cl1)
    yline(0,'LineWidth',1)
    box('off')
    
    f3=figure(3)
    subplot(1,8,i)
    bar(S(i,:),'FaceColor',cl,'EdgeColor',cl)
    hold on
    yline(prctile(no_s(i,:),99.7917),'r--','LineWidth',1)
    yline(prctile(no_s(i,:),0.2083),'r--','LineWidth',1)
    
    yline(mean(prctile(LS(i,:),99.7917)),'k--','LineWidth',1)
    yline(mean(prctile(LS(i,:),0.2083)),'k--','LineWidth',1)
    
    box('off')
    
end


f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 60, 6];
set(gca,'FontSize',9)
set(f1,'color','w');

f2.Units = 'centimeters';
f2.OuterPosition= [10, 10, 60, 6];
set(gca,'FontSize',9)
set(f2,'color','w');

f3.Units = 'centimeters';
f3.OuterPosition= [10, 10, 60, 6];
set(gca,'FontSize',9)
set(f3,'color','w');

