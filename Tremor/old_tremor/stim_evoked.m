clear all
 load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\A_group.mat','S','LS','no_s')
% % % load ('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/A_group.mat')
%  load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\F_group.mat','S','LS','no_s')
% % % load ('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/F_group.mat')

amp.s=S; amp.ls=LS; amp.nos= no_s; clear S  LS no_s
ref_ms= amp.s;
ref_text='amp.s';

iii=[1 2 3 4 5 6 8 10 11];
for numb=1:length(iii);    
clearvars -except iii numb ref_ms amp freq ref_text phase_max
% load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))
load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))

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
pre_stim_evoked;

time1=0:1/samplerate:(size(tremor1,2)-1)/samplerate;

% plot(time1,stim)
% hold on
% plot(time1(index),stim(index),'r.')
% plot(time1,tremorx)

epoch=250;
for st=1:length(index)
    if (~isnan(index(st)))
    evk_1(st,(1:2*epoch+1))= tremor1f(index(st)-epoch:index(st)+epoch);
    evk_2(st,(1:2*epoch+1))= tremor2f(index(st)-epoch:index(st)+epoch);
    evk_3(st,(1:2*epoch+1))= tremor3f(index(st)-epoch:index(st)+epoch);
    check(st,(1:2*epoch+1))= stim(index(st)-epoch:index(st)+epoch);
    emg2(st,(1:2*epoch+1))= emg1(index(st)-epoch:index(st)+epoch);
    end
end

clear ttx tty ttz ch


for i=1:12
    ta1(i,:)=mean(evk_1(find(xx==i),:)); 
    ta2(i,:)=mean(evk_2(find(xx==i),:));
    ta3(i,:)=mean(evk_3(find(xx==i),:));
    ch(i,:)=mean(check(find(xx==i),:));
    emg3(i,:)=mean(emg2(find(xx==i),:));
    
    ta1(ta1==0)=NaN;
    ta2(ta2==0)=NaN;
    ta3(ta3==0)=NaN;
    ch(ch==0)=NaN;
    emg3(emg3==0)=NaN;
end

ref=abs(ref_ms(numb,:));

amp.ns_upth(1,numb)=prctile(amp.nos(numb,:),99.7917); %bonferroni corrected for 12 comparisons
amp.ls_upth(1,numb)=prctile(amp.ls(numb,:),99.7917); 

amp.ns_dwth(1,numb)=prctile(amp.nos(numb,:),0.2083); %bonferroni corrected for 12 comparisons
amp.ls_dwth(1,numb)=prctile(amp.ls(numb,:),0.2083);

amp.ns_sig(1,numb)= (max(ref)>=amp.ns_upth(numb));
amp.ls_sig(1,numb)= (max(ref)>=amp.ls_upth(numb));
idx_phase=find(ref==max(ref));
phase_max(1,numb)=idx_phase;

% figure() 
% subplot(3,1,1)
% plot(tremor1f)
% xlim([time(index(2)).*1000-1000 time(index(2)).*1000+20000])
% ylabel('Main axis')
% box('off')
% subplot(3,1,2)
% plot(tremor2f)
% xlim([time(index(2)).*1000-1000 time(index(2)).*1000+20000])
% box('off')
% subplot(3,1,3)
% plot(tremor3f)
% xlim([time(index(2)).*1000-1000 time(index(2)).*1000+20000])
% box('off')
% % 
% fig=figure()
% subplot(5,1,1)
% plot(ta1(idx_phase,:))
% ylabel('Main axis')
% xline(250,'--')
% % ylim([-0.5 0.5])
% set(gca,'FontSize',12)
% xlim([0 500])
% xticks([0:250:500]);
% xticklabels ({'-250','0','250'})
% box('off')
% 
% subplot(5,1,2)
% plot(ta2(idx_phase,:))
% xline(250,'--')
% % ylim([-0.5 0.5])
% set(gca,'FontSize',12)
% xlim([0 500])
% xticks([0:250:500]);
% xticklabels ({'-250','0','250'})
% box('off')
% 
% subplot(5,1,3)
% plot(ta3(idx_phase,:))
% xline(250,'--')
% % ylim([-0.5 0.5])
% set(gca,'FontSize',12)
% xlim([0 500])
% xticks([0:250:500]);
% xticklabels ({'-250','0','250'})
% box('off')
% 
% subplot(5,1,4)
% plot(ch(idx_phase,:))
% xline(250,'--')
% set(gca,'FontSize',12)
% xlim([0 500])
% xticks([0:250:500]);
% xticklabels ({'-250','0','250'})
% box('off')
% 
% subplot(5,1,5)
% plot(emg3(idx_phase,:))
% xline(250,'--')
% set(gca,'FontSize',12)
% xlim([0 500])
% xticks([0:250:500]);
% xticklabels ({'-250','0','250'})
% box('off')
% 
% fig.Units = 'centimeters';
% fig.OuterPosition= [2, 2, 10, 20];
% set(fig,'color','w');
% 
% fig=figure()
% subplot(2,1,1)
% bar(ref_ms(numb,:),'FaceColor',[0.5 0.5 0.5],'EdgeColor', [0.5 0.5 0.5])
% hold on
% yline(amp.ns_upth(numb),'r--','LineWidth',1)
% yline(amp.ns_dwth(numb),'r--','LineWidth',1)
% box('off')
% ylim([-1 1])
% ylabel('non stim thresholds')
% subplot(2,1,2)
% bar(ref_ms(numb,:),'FaceColor',[0.5 0.5 0.5],'EdgeColor', [0.5 0.5 0.5])
% hold on
% yline(amp.ls_upth(numb),'k--','LineWidth',1)
% yline(amp.ls_dwth(numb),'k--','LineWidth',1)
% box('off')
% ylim([-1 1])
% ylabel('random stim thresholds')
% fig.Units = 'centimeters';
% fig.OuterPosition= [2, 2, 10, 20];
% set(fig,'color','w');

close all
end

bumps=[0 1 0 1 1 0 0 0 0];

clearvars -except bumps freq amp ref_ms ref_text phase_max



