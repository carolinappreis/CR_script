% % %%% Create smooth data and save
% %
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------------AMP
% % % %%% loading all ARC axis
% clear all
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\amp_ARC')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\am_ax')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\s_arc')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\t_arc')
%
% % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/amp_ARC')
% % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/am_ax')
% % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/s_arc')
% % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/t_arc')
%
% SO=repmat(ttall,1,3);
% S1=repmat(am_ax,1,3);
% S2=repmat(s_arc,1,3);
% S3=repmat(t_arc,1,3);
%  for ii=1:size(SO,1)
%     for i=size(ttall,2)+1:size(ttall,2)*2
%         smo_s(ii,i-12)=sum(SO(ii,(i-1:i+1)))./length(SO(ii,(i-1:i+1)));
%         smo_s1(ii,i-12)=sum(S1(ii,(i-1:i+1)))./length(S1(ii,(i-1:i+1)));
%         smo_s2(ii,i-12)=sum(S2(ii,(i-1:i+1)))./length(S2(ii,(i-1:i+1)));
%         smo_s3(ii,i-12)= sum(S3(ii,(i-1:i+1)))./length(S3(ii,(i-1:i+1)));
%     end
% end
%
% clearvars -except smo_s smo_s1 smo_s2 smo_s3
% % cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
% cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
%
% save('smooth_arc3.mat')
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------------FREQ
%
% %%% loading all FRC axis
% clear all
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\freq_FRC')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\fm_ax')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\s_frc')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\t_frc')
%
% % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/freq_FRC')
% % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/fm_ax')
% % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/s_frc')
% % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/t_frc')
%
%
% SO=repmat(ttall,1,3);
% S1=repmat(fm_ax,1,3);
% S2=repmat(s_frc,1,3);
% S3=repmat(t_frc,1,3);
%
%
% for ii=1:size(SO,1)
%     for i=size(ttall,2)+1:size(ttall,2)*2
%         smo_s(ii,i-12)=sum(SO(ii,(i-1:i+1)))./length(SO(ii,(i-1:i+1)));
%         smo_s1(ii,i-12)=sum(S1(ii,(i-1:i+1)))./length(S1(ii,(i-1:i+1)));
%         smo_s2(ii,i-12)=sum(S2(ii,(i-1:i+1)))./length(S2(ii,(i-1:i+1)));
%         smo_s3(ii,i-12)= sum(S3(ii,(i-1:i+1)))./length(S3(ii,(i-1:i+1)));
%     end
% end
%
% clearvars -except smo_s smo_s1 smo_s2 smo_s3
% % cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
% cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
% save('smooth_frc3.mat')
% %
% %
% %
% % %%%% use smoothed data

clear all
close all

metric=1;

if metric==0;
    % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/smooth_arc3')
    % load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash')
    
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\am_ax')
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\amp_NS','no_s')
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\amp_ARC','LS')
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\smooth_arc3')
    load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
    S=am_ax;
    
    cl=blushred;
    cl1=squash;
    
else
    load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/fm_ax')
    load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/freq_FRC','LS')
    load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/freq_NS','no_s')
    load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/smooth_frc3')
    load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','aegean','stone');
%     
%     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\fm_ax')
%     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\freq_FRC','LS')
%     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\freq_NS','no_s')
%     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\smooth_frc3')
%     load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','aegean','stone');
    S=fm_ax;
    
    cl=aegean;
    cl1=stone;
end


for i=1:size(smo_s,1)
%     f1=figure(1)
%     subplot(1,size(smo_s,1),i)
%     bar(smo_s(i,:),'LineWidth',0.5,'FaceColor',cl,'EdgeColor',cl)
%     hold on
%     bar(smo_s1(i,:),'LineStyle','--','LineWidth',1,'FaceColor','none','EdgeColor','k')
%     yline(0,'LineWidth',1)
%     box('off')
%     
%     f2=figure(2)
%     subplot(1,size(smo_s,1),i)
%     bar(smo_s(i,:),'FaceColor',cl,'EdgeColor',cl)
%     hold on
%     plot(smo_s2(i,:),'LineWidth',1,'Color','k')
%     plot(smo_s3(i,:),'LineWidth',1,'Color','k')
%     yline(0,'LineWidth',1)
%     box('off')
%     
%     
    %%%V1 -------------
    %     f3=figure(3)
    %     subplot(1,size(smo_s,1),i)
    %     bar(S(i,:),'FaceColor',cl,'EdgeColor',cl)
    %     hold on
    %     yline(prctile(no_s(i,:),99.7917),'r--','LineWidth',1)
    %     yline(prctile(no_s(i,:),0.2083),'r--','LineWidth',1)
    %
    %     yline(prctile(mean(LS(i,:,:),3),99.7917),'k--','LineWidth',1)
    %     yline(prctile(mean(LS(i,:,:),3),0.2083),'k--','LineWidth',1)
    %
    %     box('off')
    
    %%%V2-----------
        f3=figure(3)
        subplot(1,size(smo_s,1),i)
        bar(S(i,:),'FaceColor',cl,'EdgeColor',cl)
        hold on
        yline(prctile(no_s(i,:),99.7917),'r--','LineWidth',1)
        yline(prctile(no_s(i,:),0.2083),'r--','LineWidth',1)
        sigrn=[];
        for ii=1:12
            if S(i,ii)>prctile(LS(i,:,ii),97.5) | S(i,ii)<prctile(LS(i,:,ii),2.5)
                sigrn=[sigrn 0];
            else
                sigrn=[sigrn NaN];
            end
        end
        plot(sigrn,'kd','MarkerSize',4,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
        box('off')
    %
    %%%V3------------
%     f3=figure(3)
%     subplot(1,size(smo_s,1),i)
%     bar(S(i,:),'FaceColor',cl,'EdgeColor',cl)
%     hold on
%     yline(prctile(no_s(i,:),99.7917),'r--','LineWidth',1)
%     yline(prctile(no_s(i,:),0.2083),'r--','LineWidth',1)
%     box('off')
%     
%     f4=figure(4)
%     subplot(1,size(smo_s,1),i)
%     bar(S(i,:),'FaceColor',cl,'EdgeColor',cl)
%     hold on
%     for ii=1:12
%         lim(ii,:)=[prctile(LS(i,:,ii),2.5) prctile(LS(i,:,ii),97.5)];
%     end
%     stem(lim,'MarkerSize',2,'Color',[0.5 0.5 0.5])
%     box('off')
    
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

f3.Units = 'centimeters';
f3.OuterPosition= [10, 10, 60, 6];
set(gca,'FontSize',8)
set(f3,'color','w');

% f4.Units = 'centimeters';
% f4.OuterPosition= [10, 10, 60, 6];
% set(gca,'FontSize',8)
% set(f4,'color','w');

