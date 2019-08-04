%%%% Create smooth data and save
% clear all
% 
% %%% loading all ARC axis
% 
% % load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\A_group')
% % load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\f_ax')
% % load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Axis2\s_arc')
% % load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Axis3\t_arc')
% 
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/A_group')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/f_ax')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Axis2/s_arc')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Axis3/t_arc')
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
% cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
% save('smooth_arc3.mat')


%%%% use smoothed data
clear all
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/smooth_arc3')


close all

% for i=1:size(smo_s,1)
%     f1=figure(1)
%     subplot(1,9,i)
%     bar(smo_s(i,:),'LineWidth',1,'FaceColor',[0.5 0 0],'EdgeColor',[0.5 0 0])
%     hold on
%     plot(smo_s1(i,:),'LineWidth',1,'Color','k')
% %   bar(f_ax(new(i),:),'LineStyle','--','LineWidth',1,'FaceColor','none','EdgeColor','k')
%     yline(0,'LineWidth',1)
%     box('off')
%     
%     f2=figure(2)
%     subplot(1,9,i)
%     bar(smo_s(i,:),'FaceColor',[0.5 0 0],'EdgeColor',[0.5 0 0])
%     hold on
%     plot(smo_s2(i,:),'LineWidth',1,'Color','[0.5 0.5 0.5]')
%     plot(smo_s3(i,:),'LineWidth',1,'Color','k')
%     yline(0,'LineWidth',1)
%     box('off')
%     
% end


for i=1:size(smo_s,1)
    f1=figure(1)
    subplot(1,9,i)
    bar(smo_s(i,:),'LineWidth',0.5,'FaceColor',[0.5 0 0],'FaceAlpha',0.5,'EdgeColor',[0.5 0 0])
    hold on
     bar(smo_s1(i,:),'LineStyle','--','LineWidth',1,'FaceColor','none','EdgeColor','k')
    yline(0,'LineWidth',1)
    box('off')
    
    f2=figure(2)
    subplot(1,9,i)
    bar(smo_s(i,:),'FaceColor',[0.5 0 0],'FaceAlpha',0.5,'EdgeColor',[0.5 0 0])
    hold on
    plot(smo_s2(i,:),'LineWidth',1,'Color','[0.5 0.5 0.5]')
    plot(smo_s3(i,:),'LineWidth',1,'Color','k')
    yline(0,'LineWidth',1)
    box('off')
    
end



