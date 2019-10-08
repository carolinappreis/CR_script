clear all
close all
% load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\start_end_RS.mat')
% load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\RS_hmm_posture.mat')
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/start_end_RS.mat')
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/RS_hmm_posture.mat')


% d= [RS_t; RS_e; RS_dc];
% stackedplot(d');
% xlim([0 1200])
% dc=RS_dc{iii,1};
% data2=RS_e{iii,1};
% for i=1:size(dc,1)
%   data1(i,:)=data2(i,:)-(mean(dc(i,:)));
% end


data1=RS_t;
data=data1';data=data(:,1);
T=length(data1);
% T=repmat(t,1,(round((size(data,1)./t),0)));
% if (sum(T))~=size(data,1)
% data=data1(1:sum(T));
% end

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

for i=1
%     subplot(1,3,i)
area(Gamma)
hold on
plot(epochs,'w','LineWidth',1.5)
plot(data(:,i)*10,'k','LineWidth',1)
xlim([0 size(data,1)])
    box('off')
end



figure;
for i=1:3
    subplot(1,3,i)
    plot(data(:,i)*10)
    hold on
    area(Gamma)
    plot(epochs-0.5,'r','LineWidth',0.5)
    xlim([0 Fs*120])
    box('off')
end




hold on
plot((data(:,1)*10)+1.5)


    X=data;
    % fit = hmmspectramar(X,T,hmm,Gamma,options)

    options.win=t;
    fit = hmmspectramt(X,T,Gamma,options)

    figure()
    for state = 1:options.K
        for ii=1:size(data,2)
            subplot(size(data,2),2,ii+((state-1)*size(data,2)))
            plot(fit.state(state).f,fit.state(state).psd(:,ii,ii))
            ylim([0 1])
            box('off')
        end
    end
% legend({'state1','state2'})
% legend('boxoff')



    ep=1:t:size(data1,2);
    m_ep1=[];
    for g=1:length(ep) 
            m_ep1=[m_ep1 median(vpath(ep(g)+5*Fs:ep(g)+6*Fs-1))- median(vpath(ep(g):ep(g)+Fs-1))];
    end

    m_ep1=m_ep1';

    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/pca_ma.mat')

histogram(pca_idx{1,1}(1,:))
box('off')
xlabel('Axis')
xticks([1:1:3])
ylabel('Number of trials with max. PCA coef.')


numel(find((diff(pca_idx{1,1}))~=0))
numel(find((diff(pca_idx{2,1}))~=0))



 

    % phase=xx1(pca_idx{1,1}(:)==3);

    for w=1:12
    %      tt1{w,:}=(m_ep1(find(phase==w)));
         tt1{w,:}=(m_ep1(find(phase==w),:));
    %     tt2{w,:}=(vp_ep2(find(phase==w),:));
    end

    figure;
    bar(0:30:330,cellfun(@sum,tt1))
box('off')
   




    % evokedStateProbability