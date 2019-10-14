clear all
close all
% load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\start_end_RS.mat')
% load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\RS_hmm_spiral.mat')
% load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','aegean','stone','blushred','ocean','cream','banana');

cc=1

if cc==1
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/start_end_RS.mat')
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/RS_hmm_posture.mat')
% load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','aegean','stone','blushred','ocean','cream','banana');
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','m1_1','m1_2','m2_1','m2_2');
else
    
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/start_end_RS.mat')
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/RS_hmm_spiral.mat')
% load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','aegean','stone','blushred','ocean','cream','banana');
end

data1=RS_t;
data=data1';
% T=length(data1);
T=repmat(t,1,(round((size(data,1)./t),0)));
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


% col{1,1}(1,:)=blushred;
% col{1,1}(2,:)=stone;
% col{2,1}(2,:)=banana ;
% col{2,1}(1,:)=aegean ;
col{1,1}(1,:)=m1_1;
col{1,1}(2,:)=m1_2;
col{2,1}(2,:)=m2_1 ;
col{2,1}(1,:)=m2_2 ;




f=figure();
for i=1:size(data1,1)
   subplot(size(data1,1),1,i)
    h=area(Gamma);
    h(1).FaceColor = col{cc,1}(1,:);
    h(2).FaceColor = col{cc,1}(2,:);
    hold on
    plot(epochs,'w','LineWidth',0.5)
    plot(data(:,i)*10,'k','LineWidth',1)
    xlim([0 size(data,1)])
    box('off')
    set(gca,'FontSize',10)
set(f,'color','w');
     legend({'state1','state2'},'EdgeColor','none')
end
% ylabel('Prob of state')
% xlabel('Time')
f.Units = 'centimeters';
f.OuterPosition= [10, 10, 8, 10];


X=data;
options.win=t;
fit = hmmspectramt(X,T,Gamma,options)



f=figure()
for st = 1:options.K
    for ii=1:size(data,2)
        subplot(size(data,2),1,ii)
        plot(fit.state(st).f,squeeze(abs(fit.state(st).psd(:,ii,ii))),'LineWidth',1.5,'Color',col{cc,1}(st,:))
        hold on
        ylabel ('Power spectral Density')
        xlabel ('Frequency (Hz)')
        ylim([0 0.8])
        box('off')
         legend({'sate1','state2'})
        legend('boxoff')
%        title(['State',num2str(st)])
    end
end

f.Units = 'centimeters';
f.OuterPosition= [10, 10, 10, 15];
set(f,'color','w');
% legend({'state1','state2'})
% legend('boxoff')

ep=1:Fs*6:length(data);
for i=1:length(ep)
    if i<length(ep)
st_ep(i,:)=vpath(ep(i)+Fs:ep(i+1)-1);
    end
end

st_ep=[st_ep ;( vpath(ep(i)+Fs:(ep(i)+Fs*6)-1))'];

first=[];scnd=[];none=[];
for ii =1:size(st_ep,1)
    if sum(st_ep(ii,:))==size(st_ep,2)
        first=[first ii];
    elseif sum(st_ep(ii,:))==2*size(st_ep,2)
        scnd=[scnd ii];
    else
        none=[none ii];
    end
end

states{1,1}=first;
states{2,1}=scnd;
states{3,1}=none;

figure()
bar([numel(first) numel(scnd) numel(none)])
ylabel ('number of trials')
xticklabels({'state1', 'state2', 'state1&state2'})
box('off')
% clearvars -except states Gamma data T epochs vpath fit hmm options cc Fs
% cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA')
% save(['Tall_states',num2str(options.K),'_motor',num2str(cc),'_3ch.mat'])













ep=1:t:size(data1,2);
m_ep1=[];
for g=1:length(ep)
    m_ep1=[m_ep1 median(vpath(ep(g)+5*Fs:ep(g)+6*Fs-1))- median(vpath(ep(g):ep(g)+Fs-1))];
end

m_ep1=m_ep1';

load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/pca_ma.mat')
% load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\pca_ma.mat')

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


r=pca_idx{2,1}(1,:);
for i=1:3
prc(i,1)=(numel(find(r==i)))./length(r);
end

f=figure;
pie(prc)
box('off')
f.Units = 'centimeters';
f.OuterPosition= [10, 10, 8, 8];
set(gca,'FontSize',10)
set(f,'color','w');


numel(find((diff(pca_idx{1,1}))~=0))
numel(find((diff(pca_idx{2,1}))~=0))


phase=xx1

for w=1:12
         tt1{w,:}=(m_ep1(find(phase==w)));
    tt1{w,:}=(m_ep1(find(phase==w),:));
        tt2{w,:}=(vp_ep2(find(phase==w),:));
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


clear all

load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/states2_motor2_3ch.mat')
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','m1_1','m1_2','m2_1','m2_2','blushred','stone','aegean');
data1=data';
Fs=17;
t=5*Fs;

% col{1,1}(1,:)=m1_1;
% col{1,1}(2,:)=m1_2;
% col{2,1}(2,:)=m2_1 ;
% col{2,1}(1,:)=m2_2 ;


col{1,1}(1,:)=m1_1;
col{1,1}(2,:)=blushred;
col{2,1}(2,:)=m2_2 ;
col{2,1}(1,:)=m2_1 ;


f=figure();
for i=3
%     :size(data1,1)
  % subplot(size(data1,1),1,i)
    h=area(Gamma);
    h(1).FaceColor = col{cc,1}(1,:);
    h(2).FaceColor = col{cc,1}(2,:);
     h(1).EdgeColor = 'none';
    h(2).EdgeColor = 'none';
    hold on
     plot(epochs,'w','LineWidth',0.5)
    plot(data(:,i)*10,'k','LineWidth',1)
    xlim([0 size(data,1)])
    box('off')
    ylabel('Probability')
xlabel('Time (seconds)')
xtick=([0:1700:10200]);
xticklabel=({'0', '100', '200', '300','400','500','600'});
    set(gca,'FontSize',14)
set(f,'color','w');
%      legend({'state1','state2'},'EdgeColor','none')
end

f.Units = 'centimeters';
f.OuterPosition= [10, 10, 14, 8];


X=data;
options.win=t;
fit = hmmspectramt(X,T,Gamma,options)



f=figure()
for st = 1:options.K
    for ii=1:size(data,2)
%         subplot(size(data,2),1,ii)

        plot(fit.state(st).f,squeeze(abs(fit.state(st).psd(:,ii,ii))),'LineWidth',1.5,'Color',col{cc,1}(st,:))
        hold on
        ylabel ('PSD')
        xlabel ('Frequency (Hz)')
        ylim([0 0.8])
        box('off')
%          legend({'sate1','state2'})
%         legend('boxoff')
%        title(['State',num2str(st)])

    end
end
f.Units = 'centimeters';
f.OuterPosition= [10, 10, 12, 12];
set(f,'color','w');
 set(gca,'FontSize',14)

% legend({'state1','state2'})
% legend('boxoff')

ep=1:Fs*6:length(data);
for i=1:length(ep)
    if i<length(ep)
st_ep(i,:)=vpath(ep(i)+Fs:ep(i+1)-1);
    end
end

st_ep=[st_ep ;( vpath(ep(i)+Fs:(ep(i)+Fs*6)-1))'];

first=[];scnd=[];none=[];
for ii =1:size(st_ep,1)
    if sum(st_ep(ii,:))==size(st_ep,2)
        first=[first ii];
    elseif sum(st_ep(ii,:))==2*size(st_ep,2)
        scnd=[scnd ii];
    else
        none=[none ii];
    end
end

states{1,1}=first;
states{2,1}=scnd;
states{3,1}=none;

f=figure()
bar([numel(first) numel(scnd) numel(none)])
ylabel ('number of trials')
xticklabels({'state1', 'state2', 'state1 & 2'})
box('off')
f.Units = 'centimeters';
f.OuterPosition= [10, 10, 12, 12];
set(f,'color','w');
 set(gca,'FontSize',14)
