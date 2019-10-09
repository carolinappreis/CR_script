
   
clear all
iii=[1];
numb=1;
DBS_Fpeak

load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/01_NS_PS.mat')

in2=1;

if in2==1
    in=3;
elseif in2==2 % other axis 1
    in=6;
elseif in2==3 % other axis 2
    in=7;
end
data=SmrData.WvData;
samplerateold=SmrData.SR;
%%------------------------
tremor3=data([3 6 7],:);

[b,a]=butter(2,[(Fpeak-2)/(0.5*samplerateold) (Fpeak+2)/(0.5*samplerateold)],'bandpass');

[m,n]=butter(3,[1.5/(0.5*samplerateold)],'low');

for i=1:size(tremor3,1)
    dc_t1(i,:)=filtfilt(m,n,tremor3(i,:));
    filt_t1(i,:)=filtfilt(b,a,tremor3(i,:));
    env_t1(i,:)=abs(hilbert(filt_t1(i,:)));
end
data1=vertcat(filt_t1,env_t1,dc_t1);

Fs=(Fpeak+2)*2+5;
time=0:1/samplerateold:(size(data1,2)-1)/samplerateold;
ts=timeseries(data1,0:(1/samplerateold):((size(data1,2)-1)/samplerateold));
ts1=resample(ts,0:1/Fs:((size(data1,2)-1)/samplerateold),'linear');
data2(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;


filt_t3=data2(1:3,:);
env_t3=data2(4:6,:);
dc_t3=data2(7:9,:);

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','aegean','stone','blushred','ocean','cream','banana');


data1=filt_t3(1,:);
data=data1';
T=length(data1);
% T=repmat(t,1,(round((size(data,1)./t),0)));

options=struct();
options.initrep=3;
options.K=3;
options.useParallel=0;
options.Fs=Fs;
options.embeddedlags=0;
options.covtype = 'diag';
options.zeromean = 0;
options.order=0;

[hmm,Gamma,~,vpath]=hmmmar(data,T,options);
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
    plot(data*10,'k','LineWidth',1)
    xlim([0 size(data,1)])
    box('off')
end
% ylabel('Prob of state')
% xlabel('Time')
% box('off')
% f.Units = 'centimeters';
% f.OuterPosition= [10, 10, 10, 5];
% set(gca,'FontSize',10)
% set(f,'color','w');

X=data;
options.win=Fs*10;
fit = hmmspectramt(X,T,Gamma,options)



f=figure()
for state = 1:options.K
    for ii=1:size(data,2)
        plot(fit.state(state).f,fit.state(state).psd,'LineWidth',2)
        hold on
%         ylim([0 1])
        box('off')
    end
end
% ylabel ('Power spectral Density')
% xlabel ('Frequency (Hz)')
% f.Units = 'centimeters';
% f.OuterPosition= [10, 10, 8, 8];
% set(gca,'FontSize',10)
% set(f,'color','w');
% legend({'state1','state2'})
% legend('boxoff')
