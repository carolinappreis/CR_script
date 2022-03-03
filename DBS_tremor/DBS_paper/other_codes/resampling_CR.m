
cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')


clear all
close all
cond={'NS';'HF';'C'};
cohort=[ 1 3 4 6];

iii=3;
trial=1;
co=2;

load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),num2str(trial),'_SH.mat'));
data=(NS1)';

plot3(NS1(:,3),NS1(:,1),NS1(:,2))

data(:,find(data(1,:)==-1))=[];
com_sig=sqrt((data(1,:).^2)+(data(2,:).^2));

% [signal,time]=resample(com_sig,NS1(:,3),1);


x = data(3,:) ;
y = [data(1,:);data(2,:);com_sig] ;
[xn, index] = unique(x);
yn=y(:,index);


N=length(yn);
xi = linspace(min(xn),max(xn),N) ;
for i=1:size(yn,1)
yi(i,:) = interp1(xn,yn(i,:),xi,'linear') ;
end

samplerate=1000/(median(diff(xi)));
time=0:1/samplerate:((size(yi,2)-1)/samplerate);


%%% interpolate

samplerate2=100;
res=0:1/samplerate2:((size(yi,2)-1)/samplerate);
[a,b]=  butter(2, [2/(0.5*samplerate2) 8/(0.5*samplerate2)], 'bandpass'); %15

for i=1:size(yi,1)
data1(i,:)=interp1(time,yi(i,:),res);
filt_cs(i,:)=(filtfilt(a,b,data1(i,:)));
env_cs(i,:)=abs(hilbert(filt_cs(i,:)));
end

% figure
% subplot(3,1,1)
% plot(data(3,:),com_sig,'k.')
% title('raw')
% subplot(3,1,2)
% plot(xn,yn(3,:),'r.')
% title('unique')
% subplot(3,1,3)
% plot(xi,yi(3,:),'.b')
% title('unique resamp')
% 
% figure
% plot(time,yi(3,:),'b.')
% title('unique resamp zero')
% hold on
% plot(res,data1(3,:),'g.')
% title('unique resamp zero upsamp')
% 
% figure
% plot3(res,data1(1,:),data1(2,:))
% title('3D spiral resamp upsamp')
% 
% figure
% subplot(2,1,1)
% plot(data(3,:)./1000,data(1,:))
% title('original x')
% subplot(2,1,2)
% plot(res,data1(1,:))
% title('final x')

