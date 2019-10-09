clear all
close all
% load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\start_end_RS.mat')
% load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\RS_hmm_posture.mat')
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/start_end_RS.mat')
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/RS_hmm_posture.mat')
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','aegean','stone','blushred','ocean','cream','banana');


data1=RS_t;
data=data1';data=data(:,1);
T=length(data1);
% T=repmat(t,1,(round((size(data,1)./t),0)));

options=struct();
options.initrep=3;
options.K=2;
options.useParallel=0;
options.Fs=Fs;
options.embeddedlags=0;
options.covtype = 'diag';
options.zeromean = 0;
options.order=0;

[hmm,Gamma,~,vpath]=hmmmar(data,T,options);
epochs=zeros(1,size(data,1)); epochs(1:t:end)=1;
time=0:1/Fs:(size(data,1)-1)/Fs;



% col(1,:)=banana ;
% col(2,:)=aegean ;

col(1,:)=blushred;
col(2,:)=stone;


f=figure();
for i=1
    %     subplot(1,3,i)
    h=area(Gamma);
    h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
    
    hold on
%      plot(epochs,'r','LineWidth',0.5)
    plot(data(:,i)*10,'k','LineWidth',1)
    xlim([0 size(data,1)])
    box('off')
end
ylabel('Prob of state')
xlabel('Time')
box('off')
f.Units = 'centimeters';
f.OuterPosition= [10, 10, 10, 5];
set(gca,'FontSize',10)
set(f,'color','w');

X=data;
options.win=t;
fit = hmmspectramt(X,T,Gamma,options)



f=figure()
for state = 1:options.K
    for ii=1:size(data,2)
        plot(fit.state(state).f,fit.state(state).psd(:,ii,ii),'LineWidth',2, 'Color',col(state,:))
        hold on
        ylim([0 1])
        box('off')
    end
end
ylabel ('Power spectral Density')
xlabel ('Frequency (Hz)')
f.Units = 'centimeters';
f.OuterPosition= [10, 10, 8, 8];
set(gca,'FontSize',10)
set(f,'color','w');
% legend({'state1','state2'})
% legend('boxoff')



ep=1:t:size(data1,2);
m_ep1=[];
for g=1:length(ep)
    m_ep1=[m_ep1 median(vpath(ep(g)+5*Fs:ep(g)+6*Fs-1))- median(vpath(ep(g):ep(g)+Fs-1))];
end

m_ep1=m_ep1';

load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/pca_ma.mat')


f=figure();
histogram(pca_idx{1,1}(1,:),'LineWidth',0.5,'FaceColor',col(1,:),'EdgeColor','none')
box('off')
xlabel('Axis')
xticks([1:1:3])
ylabel('#trials with max. PCA coef.')
box('off')
f.Units = 'centimeters';
f.OuterPosition= [10, 10, 8, 8];
set(gca,'FontSize',10)
set(f,'color','w');



numel(find((diff(pca_idx{1,1}))~=0))
numel(find((diff(pca_idx{2,1}))~=0))


phase=xx1

for w=1:12
    %      tt1{w,:}=(m_ep1(find(phase==w)));
    tt1{w,:}=(m_ep1(find(phase==w),:));
    %     tt2{w,:}=(vp_ep2(find(phase==w),:));
end

f=figure()
bar(0:30:330,cellfun(@sum,tt1),'LineWidth',0.5,'FaceColor',col(1,:),'EdgeColor',col(1,:))
hold on
yline(0,'LineWidth',1)
ylabel('Prob of change of state')
xlabel('stimulation phase (degrees)')
box('off')
f.Units = 'centimeters';
f.OuterPosition= [10, 10, 8, 8];
set(gca,'XTickLabelRotation',45)
set(gca,'FontSize',10)
set(f,'color','w');


