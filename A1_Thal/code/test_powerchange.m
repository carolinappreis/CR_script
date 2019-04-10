% clear all
% cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat/all_regions')
% load ('SNr_nc_fig','all_contacts','surr2','surr_s')
% out3_real= all_contacts;
% clear all_contacts
% out3_sur= surr2;
% clear surr2
% 
% ind=2;
% clear A; A(1:size(out3_real,2),1:551)=out3_real(ind,:,300:850);
% clear B; B(1:size(out3_sur,1),1:551)=out3_sur(:,300:850);
% 
% 
% 
% data=[];
% for i=1:size(A,1)
% data(i,:)=smooth((A(i,:)));
% end
% 
% for i=1:size(data,1)
% max_data(i,1)=find(data(i,10:490)==max(data(i,10:490)))
% end
% max_data=max_data+10;
% 
% 
% 
% 
% time=[1:500];
% figure(1)
% for i=1:size(max_data,1)
% plot(data(i,:))
% hold on
% plot(time(max_data(i,1)),data(i,max_data(i,1)),'r*')
% % close all
% end
% 
% boxplot(max_data)
% plot(max_data,'r.')

clear all
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')
load ('SNR_Nullcoh','all_contacts','surr2','surr_s')
out3_real= all_contacts;
clear all_contacts
out3_sur= surr2;
clear surr2

figure(1)
for ind=1:2;
clear A; A(1:size(out3_real,2),1:551)=out3_real(ind,:,300:850);
clear B; B(1:size(out3_sur,1),1:551)=out3_sur(:,300:850);
subplot(3,1,ind)
plot(mean(A))
hold on
plot(mean(B))
end

subplot(3,1,3)
load('SNR_power','power_nc')
plot(mean(power_nc(:,5:80)))
title ('non_coherent contacts')


clear all
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')
load ('SNR_fig','all_contacts','surr2','surr_s')
out3_real= all_contacts;
clear all_contacts
out3_sur= surr2;
clear surr2
figure(2)
for ind=1:2;
clear A; A(1:size(out3_real,2),1:551)=out3_real(ind,:,300:850);
clear B; B(1:size(out3_sur,1),1:551)=out3_sur(:,300:850);
subplot(3,1,ind)
plot(mean(A))
hold on
plot(mean(B))
end

subplot(3,1,3)
load('SNR_power','power_c')
plot(mean(power_c(:,5:80)))
title ('coherent contacts')


clear all
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat/all_regions')
load ('SNR_nc_fig','all_contacts','surr2','surr_s')
out3_real= all_contacts;
clear all_contacts
out3_sur= surr2;
clear surr2
figure(3)
for ind=1:2;
clear A; A(1:size(out3_real,2),1:551)=out3_real(ind,:,300:850);
clear B; B(1:size(out3_sur,1),1:551)=out3_sur(:,300:850);
subplot(3,1,ind)
plot(mean(A))
hold on
plot(mean(B))
end

subplot(3,1,3)
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')
load('SNR_power','power_c','power_nc')
power_all=vertcat(power_c, power_nc)
plot(mean(power_all(:,5:80)))
title ('All contacts')



%%-------------------------------


clear all
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')
load ('SNR_Nullcoh','all_contacts','surr2','surr_s')
out3_real= all_contacts;
clear all_contacts
out3_sur= surr2;
clear surr2

figure(4)
for ind=1:2;
clear A; A(1:size(out3_real,2),1:551)=out3_real(ind,:,300:850);
clear B; B(1:size(out3_sur,1),1:551)=out3_sur(:,300:850);
subplot(3,1,ind)
plot((A)')
hold on
plot((B)')
end

subplot(3,1,3)
load('SNR_power','power_nc')
plot(power_nc(:,5:80)')
ylim ([0 1*10^6])
title ('non_coherent contacts')


clear all
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')
load ('SNR_fig','all_contacts','surr2','surr_s')
out3_real= all_contacts;
clear all_contacts
out3_sur= surr2;
clear surr2
figure(5)
for ind=1:2;
clear A; A(1:size(out3_real,2),1:551)=out3_real(ind,:,300:850);
clear B; B(1:size(out3_sur,1),1:551)=out3_sur(:,300:850);
subplot(3,1,ind)
plot((A)')
hold on
plot((B)')
end

subplot(3,1,3)
load('SNR_power','power_c')
plot(power_c(:,5:80)')
ylim ([0 1*10^6])
title ('coherent contacts')


clear all
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat/all_regions')
load ('SNR_nc_fig','all_contacts','surr2','surr_s')
out3_real= all_contacts;
clear all_contacts
out3_sur= surr2;
clear surr2
figure(6)
for ind=1:2;
clear A; A(1:size(out3_real,2),1:551)=out3_real(ind,:,300:850);
clear B; B(1:size(out3_sur,1),1:551)=out3_sur(:,300:850);
subplot(3,1,ind)
plot((A)')
hold on
plot((B)')
end

subplot(3,1,3)
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')
load('SNR_power','power_c','power_nc')
power_all=vertcat(power_c, power_nc)
plot(power_all(:,5:80)')
ylim ([0 1*10^6])
title ('All contacts')

