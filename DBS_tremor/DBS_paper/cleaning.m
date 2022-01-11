
cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')

clear all
close all
cond={'NS';'HF';'C'};
cohort=[ 1 3 4 6];

iii=2;
trial=1;
co=2;

load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),num2str(trial),'_SH.mat'));
data=(NS1)';

plot3(NS1(:,3),NS1(:,1),NS1(:,2),'.')

data(:,find(data(1,:)==-1))=[];
data(isnan(data(3,:)))
com_sig=sqrt((data(1,:).^2)+(data(2,:).^2));

x = data(3,:) ;
y = [data(1,:);data(2,:);com_sig] ;
[xn, index] = unique(x);
yn=y(:,index);

yi=[];
N=length(yn);
xi = linspace(min(xn),max(xn),N) ;
for i=1:size(yn,1)
    yi(i,:) = interp1(xn,yn(i,:),xi,'linear') ;
end

samplerate=1000/(median(diff(xi))); %%% samplerate to change from miliseconds to seconds
time=0:1/samplerate:((size(yi,2)-1)/samplerate); %%% adjust x1 to start in 0

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

figure
subplot(2,1,1)
plot(data(3,:),data(1,:))
title('original x')
subplot(2,1,2)
plot(res,data1(1,:))
title('final x')


%%% clean outliers
threshold=(mean(env_cs(3,:))+3*(std(env_cs(3,:))));
indexexceed=find(env_cs(3,:)>threshold);
diffindex=diff(indexexceed);
pnts=find(diffindex>1);
if ~isempty(pnts)
    begin(1,:)=indexexceed(pnts+1);
    ending(1,:)=indexexceed(pnts);
    begin2=[indexexceed(1) begin];
    ending2=[ending indexexceed(end)];
    
    ind_b=[];
    for i=1:(length(begin2))-1
        if (ending2(i)-begin2(i))>=samplerate2
            ind_b=[ind_b i];
        end
    end
    
    begin3=begin2(ind_b);
    ending3=ending2(ind_b);
end

% figure
% plot(res,filt_cs(3,:))
% hold on
% plot(res,env_cs(3,:))
% yline(threshold)
% plot(res(begin3),env_cs(3,begin3),'k.')

plot3(res,data1(1,:),data1(2,:))



%%% adjust matriz with point of necessary here
points=cell(size(cohort,2),2,size(cond,1));

points{2,1,1}=[[1:6886]];
% points{2,1,1}=[[1:4126] [5658:6022] [6452:6886]];
points{2,2,1}=[];
points{2,1,2}=[[1:1978] [5692:length(res)]];
points{2,1,3}=[];

points{3,1,1}=[[1:584] [6234:6343]];
points{3,2,1}=[[1:2339] [8266:8443]];
points{3,1,2}=[[]];
points{3,1,3}=[[]];

points{4,1,1}=[[1:5453] [11680 12590] [19540:length(res)]];
points{4,1,2}=[[1:5494]];
points{4,1,2}=[];



tempo=[];
signal=[];

res(points{iii,trial,co})=NaN;
data1(:,points{iii,trial,co})=NaN;
tt=res(~isnan(res));
tempo=0:1/samplerate2:((size(tt,2)-1)/samplerate2);
for i=1:3
    signal(i,:)=data1(i,~isnan(data1(i,:)));
end

figure
subplot(1,2,1)
plot(tempo,signal(3,:))
subplot(1,2,2)
plot(signal(1,:),signal(2,:),'.')
title('after')

close all

plot3(tempo,signal(1,:),signal(2,:),'.')

t_spi=cell(size(cohort,2),2,size(cond,1));

% t_spi{2,1,1}=[[1 1939];[1960 4304];];
t_spi{2,1,1}=[1 length(tempo)];
t_spi{2,2,1}=[[1 2534];[2534 5228];[5228 length(tempo)]];
t_spi{2,1,2}=[[1 1892];[1894 length(tempo)]];
t_spi{2,1,3}=[[47 2105];[2105 4072];[4072 6236]];

t_spi{3,1,1}=[[1 1687];[1687 2947];[2947 4354];[4354 5601];[5650 length(tempo)]];
t_spi{3,2,1}=[[1 1668];[1668 2979];[2979 4341];[4341 length(tempo)]];
t_spi{3,1,2}=[[1 1130];[1130 2517];[2517 4207];[4207 5605];[5605 length(tempo)]];
t_spi{3,1,3}=[[1 1420];[1420 2647];[2647 3834];[3834 4923];[4923 6132];[6132 7283];[7283 8507];[8507 9827];[9827 10990];[10990 length(tempo)]];

t_spi{4,1,1}=[[1 6260];[7209 length(tempo)]];
t_spi{4,1,2}=[[1 5236];[5236 length(tempo)]];
t_spi{4,1,3}=[[1 5932];[5932 length(tempo)]];






figure
for mi=1:size(t_spi{iii,trial,co}(:,1),1)
     subplot(1,size(t_spi{iii,trial,co}(:,1),1),mi)
%              subplot(2,size(t_spi{iii,trial,co}(:,1),1)/2,mi)
    plot(signal(1,t_spi{iii,trial,co}(mi,1):t_spi{iii,trial,co}(mi,2)),signal(2,t_spi{iii,trial,co}(mi,1):t_spi{iii,trial,co}(mi,2)),'.')
end

cd(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals'))

filename=strcat('P0',num2str(cohort(iii)),'_clean_',num2str(cond{co,1}),num2str(trial),'_SH.mat');

clearvars -except samplerate2 signal tempo co cohort iii trial cond filename

 save(filename)