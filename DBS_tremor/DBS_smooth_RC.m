% % % %%% Create smooth data and save
% % %
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------------AMP
% % % %%% loading all ARC axis
% clear all
% % load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_amp_ARC')
% % load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\am_ax')
% % load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\s_arc')
% % load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\t_arc')
% 
% load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_amp_ARC')
% load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/am_ax')
% load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/s_arc')
% load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/t_arc')
% 
% for pp=1:2
%     
% SO=repmat(ttall{pp,1},1,3);
% S1=repmat(squeeze(am_ax(pp,:,:))',1,3);
% S2=repmat(squeeze(s_arc(pp,:,:))',1,3);
% S3=repmat(squeeze(t_arc(pp,:,:))',1,3);
% 
%     for ii=1:size(SO,1)
%         for i=size(ttall{1,1},2)+1:size(ttall{1,1},2)*2
%             smo_s(pp,ii,i-12)=sum(SO(ii,(i-1:i+1)))./length(SO(ii,(i-1:i+1)));
%             smo_s1(pp,ii,i-12)=sum(S1(ii,(i-1:i+1)))./length(S1(ii,(i-1:i+1)));
%             smo_s2(pp,ii,i-12)=sum(S2(ii,(i-1:i+1)))./length(S2(ii,(i-1:i+1)));
%             smo_s3(pp,ii,i-12)= sum(S3(ii,(i-1:i+1)))./length(S3(ii,(i-1:i+1)));
%         end
%     end
% end
% 
% clearvars -except smo_s smo_s1 smo_s2 smo_s3
%  cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data')
% % cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data')
% 
% save('smooth_arc3.mat')
% %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------------FREQ
% %
% % %%% loading all FRC axis
% clear all
% % load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_freq_FRC')
% % load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\fm_ax')
% % load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\s_frc')
% % load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\t_frc')
% 
% load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_freq_FRC')
% load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/fm_ax')
% load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/s_frc')
% load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/t_frc')
% 
% 
% for pp=1:2
%     
% SO=repmat(ttall{pp,1},1,3);
% S1=repmat(squeeze(fm_ax(pp,:,:))',1,3);
% S2=repmat(squeeze(s_frc(pp,:,:))',1,3);
% S3=repmat(squeeze(t_frc(pp,:,:))',1,3);
% 
%     for ii=1:size(SO,1)
%         for i=size(ttall{1,1},2)+1:size(ttall{1,1},2)*2
%             smo_s(pp,ii,i-12)=sum(SO(ii,(i-1:i+1)))./length(SO(ii,(i-1:i+1)));
%             smo_s1(pp,ii,i-12)=sum(S1(ii,(i-1:i+1)))./length(S1(ii,(i-1:i+1)));
%             smo_s2(pp,ii,i-12)=sum(S2(ii,(i-1:i+1)))./length(S2(ii,(i-1:i+1)));
%             smo_s3(pp,ii,i-12)= sum(S3(ii,(i-1:i+1)))./length(S3(ii,(i-1:i+1)));
%         end
%     end
% end
% 
% clearvars -except smo_s smo_s1 smo_s2 smo_s3
%  cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data')
% % cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data')
% save('smooth_frc3.mat')

% %
% %
% % %%%% use smoothed data

clear all
close all

metric=0;

if metric==0;
%         load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/am_ax')
%         load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_amp_ARC','LS','tt1','ttall')        
%         load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/NS_PS_result','idv_NS')
%         load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/smooth_arc3')
%         load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash');
    %
    
    load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\am_ax')
    load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\NS_PS_result','idv_NS')
    load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_amp_ARC','LS','tt1','ttall')
    load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\smooth_arc3')
    load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
    
    
    S=ttall;
    cl=blushred;
    cl1=squash;
    
else
%         load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/fm_ax')
%         load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_freq_FRC','LS','tt1','ttall')
%         load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/fNS_PS_result','fidv_NS')
%         load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/smooth_frc3')
%         load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','aegean','stone');

    load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\fm_ax')
    load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_freq_FRC','LS','tt1','ttall')
    load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\fNS_PS_result','fidv_NS')
    load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\smooth_frc3')
    load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','aegean','stone');


    S=ttall;
    cl=aegean;
    cl1=stone;
end

for i=1:size(smo_s,1)
    f=figure(1)
    subplot(1,size(smo_s,1),i)
    bar(smo_s(i,:),'LineWidth',0.5,'FaceColor',cl,'EdgeColor',cl)
    hold on
    bar(smo_s1(i,:),'LineStyle','--','LineWidth',1,'FaceColor','none','EdgeColor','k')
    yline(0,'LineWidth',1)
    box('off')

    f2=figure(2)
    subplot(1,size(smo_s,1),i)
    bar(smo_s(i,:),'FaceColor',cl,'EdgeColor',cl)
    hold on
    plot(smo_s2(i,:),'LineWidth',1,'Color','k')
    plot(smo_s3(i,:),'LineWidth',1,'Color','k')
    yline(0,'LineWidth',1)
    box('off')

    f3=figure(3)
    subplot(1,size(S,1),i)
    bar(S{i,1},'FaceColor',cl,'EdgeColor',cl)
    hold on
    yline(prctile(idv_NS(i,1,:),99.7917),'r--','LineWidth',1)
    yline(prctile(idv_NS(i,1,:),0.2083),'r--','LineWidth',1)
    
    yline(prctile(LS{i,1},99.7917),'k--','LineWidth',1)
    yline(prctile(LS{i,1},0.2083),'k--','LineWidth',1)
    
    box('off')
    
end


% f1.Units = 'centimeters';
% f1.OuterPosition= [10, 10, 60, 6];
% set(gca,'FontSize',8)
% set(f1,'color','w');
% 
% f2.Units = 'centimeters';
% f2.OuterPosition= [10, 10, 60, 6];
% set(gca,'FontSize',8)
% set(f2,'color','w');
% %
% f3.Units = 'centimeters';
% f3.OuterPosition= [10, 10, 60, 6];
% set(gca,'FontSize',8)
% set(f3,'color','w');

i=2;
fig=[smo_s(i,1,:);smo_s2(i,1,:);smo_s3(i,1,:)];

%%------- POSTER IMAGES


for i=1:size(fig,1)
f=figure()
bar(0:30:330,fig(i,:),'LineWidth',0.5,'FaceColor',cl,'EdgeColor',cl)
hold on
yline(0,'LineWidth',1)
ylim([-0.5 0.5])
yticks([ -0.5:0.25:0.5])
box('off')
ylabel('Change in tremor severity')
xlabel('Stimulation phase (degrees)')
    
f.Units = 'centimeters';
f.OuterPosition= [10, 10, 10, 10];
set(gca,'XTickLabelRotation',45)
set(gca,'FontSize',12)
set(f,'color','w');
end
