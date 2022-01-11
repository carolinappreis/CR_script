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
samplerate=floor(1000/median(diff(data(3,:))));
ti=0:1/samplerate:((size(data,2)-1)/samplerate);
figure
subplot(1,3,1)
plot(ti,data(1,:))
subplot(1,3,2)
plot(ti,data(2,:))
subplot(1,3,3)
plot(data(1,:),data(2,:))
title('before')
%%% filt & var
com_data=sqrt((data(1,:).^2)+(data(2,:).^2));
[a,b]=  butter(2, [2/(0.5*samplerate) 8/(0.5*samplerate)], 'bandpass'); %15
filt_cs=(filtfilt(a,b,com_data));
env_d=abs(hilbert(filt_cs));
outlier=find(env_d>=(5*(mean(env_d))));


plot(ti,env_d)
hold on
yline(5*mean(env_d))
plot(ti,env_d)
hold on
yline(5*mean(env_d))


points=cell(size(cohort,2),2,size(cond,1));
points{3,1,1}=[[1:14] [4927:4953]];
points{3,2,1}=[[1:12]];
points{3,1,2}=[4277:4301];
points{3,1,3}=[[6507:6529]  [11172:11182]];
%
%
%
ti(points{iii,trial,co})=NaN;
tt=ti(~isnan(ti));
res=0:1/samplerate:((size(tt,2)-1)/samplerate);
data(:,points{iii,trial,co})=NaN;
dd1=[];


for ix=1:2
dd1(ix,:)=data(ix,~isnan(data(ix,:)));
end


sig_data=sqrt((dd1(1,:).^2)+(dd1(2,:).^2));
[c,d]=  butter(2, [2/(0.5*samplerate) 8/(0.5*samplerate)], 'bandpass'); %15
filt_sig=(filtfilt(c,d,sig_data));
env=abs(hilbert(filt_sig));
plot(res,filt_sig)
hold on
plot(res,env)