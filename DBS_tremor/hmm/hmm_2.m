clear all
close all
load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\start_end_RS.mat')
load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\RS_20fs_9ch_posture.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/start_end_RS.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/RS_20fs_9ch_posture_beforeafter.mat')
iii=2;
%  t(2)=140;

d= [RS_t{iii,1}; RS_e{iii,1}; RS_dc{iii,1}];
stackedplot(d');
% xlim([0 1200])
% dc=RS_dc{iii,1};
% data2=RS_e{iii,1};
% for i=1:size(dc,1)
%   data1(i,:)=data2(i,:)-(mean(dc(i,:)));
% end


data1=RS_e{iii,1};

data=data1(3,:);data=data';
T=repmat(t(iii),1,((length(data)./t(iii))));

options=struct();
options.initrep=3;
options.K=2;
options.useParallel=1;
options.Fs=Fs;
options.embeddedlags=0;
options.covtype = 'diag';
options.zeromean = 0;
options.order=0;

% options=struct();
% options.initrep=3;
% options.K=3;
% options.useParallel=1;
% options.Fs=Fs;
% % options.embeddedlags=0;
% % options.covtype = 'diag';
% options.zeromean = 0;
% options.order=7;

[hmm,Gamma,~,vpath]=hmmmar(data,T,options);
epochs=zeros(1,size(data,1)); epochs(1:t(iii):end)=2;

figure;area(Gamma)
hold on
plot(epochs,'r','LineWidth',1.5)
plot(data','k','LineWidth',2)

X=data;
% fit = hmmspectramar(X,T,hmm,Gamma,options)
fit = hmmspectramt(X,T,Gamma,options)

figure()
for state = 1:options.K
for ii=1:1
    subplot(1,2,ii+((state-1)*1))
    plot(fit.state(state).f,fit.state(state).psd(:,ii,ii))
    hold on
end
end
legend({'state1','state2'})
legend('boxoff')

figure; plot(vpath,'LineWidth',3);ylim([0.8 2.2])


ep=1:t(iii):size(data1,2)+1;
vp_ep1=[];
vp_ep2=[];
for g=1:length(ep)
    %         vp_ep=[ vp_ep  vpath(ep(g):(ep(g+1)-1))]; %% vpath of the 5
    %         seconds
    %           vp_ep=[ vp_ep  vpath(((ep(g)+4*Fs)-1):(ep(g+1)-1))]; %% vpath of the last sec of stim
    if g+1<=length(ep)
        
        vp_ep1=[ vp_ep1  vpath((ep(g)+5*Fs):(ep(g)+6*Fs-1))];
        vp_ep2=[ vp_ep2  vpath(ep(g):(ep(g)+Fs-1))];
    end
end

vp_ep1=vp_ep1';
vp_ep2=vp_ep2';
% 
% vp_ep=vp_ep1;
% 
% figure()
% 
% for no=1:size(vp_ep,1)
%     for st=1:options.K
%         vp_res(no,st)= sum(vp_ep(no,:)==st)/size(vp_ep,2);
%         vp_res1(no,:)=mean(vp_ep(no,:)==st);
%     end
% end

phase=xx1;


%     tt=NaN(12,options.K);


for w=1:12
    tt1{w,:}=(vp_ep1(find(phase==w),:));
    tt2{w,:}=(vp_ep2(find(phase==w),:));
end

figure()
for i=1:12
    subplot(12,1,i)
    bar(tt1{i,1})
end

% bar(tt(:,1))
% box('off')





% evokedStateProbability