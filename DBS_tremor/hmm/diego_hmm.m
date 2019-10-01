clear all
close all
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/tmp_blah.mat')

options0.embeddedlags=0;
  
[~,b] = pca(data1(:,end-2:end),'NumComponents',2); 
figure;scatter(b(:,1),b(:,2))
[hmm,Gamma,~,vpath]=hmmmar(b,T,options0);
figure; area(Gamma(1:5001,:))
figure; plot(vpath(1:5001),'LineWidth',3);
hold on
plot(data1(1:5001,end)+1,'k','LineWidth',3)

