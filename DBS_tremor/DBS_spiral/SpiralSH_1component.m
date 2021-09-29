% clear
% load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P01_NS1_SH.mat')
% % % bins=find(NS1(:,1)==-1);
% % % for y=1:length(bins)
% % %     if(y+1)<length(bins)
% % %         data1{1,y}=NS1(bins(y)+1:bins(y+1)-1,:);
% % %     end
% % % end
% % % data1{1,length(bins)-1}=NS1(bins(y)+1:end,:);
% % %
% % %
% % %
% % % % find(NS1(:,1)<10) %%% to find turn new whole spiral
% % %
% % % %%% Data by spiral
% % % %%% for nr=1:size(data1,2)
% % % % nr=1;
% % % % data=(data1{1,nr})';
% % %
% % % %%%OR all data
% % %
% % % NS1(bins,:)=[];
% % % data=NS1';
clear
clear
cohort=[ 1 3 4 6];
iii=cohort(3);
trial=1;
close all

% load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P0',num2str(iii),'_HF',num2str(trial),'_SH.mat'))

 load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P0',num2str(iii),'_C',num2str(trial),'_SH.mat'))
data=NS1';
samplerate=floor(1000/median(diff(data(3,:))));

if iii==4 
    Fpeak=3;
else
    Fpeak=4;
end

[a,b]=  butter(2, [(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
for o=1:2
    filt(o,:)=(filtfilt(a,b,data(o,:))).*(10*9.81/0.5);
end


plot(filt(1,:))

x = [filt(1,:); filt(2,:)];
[pc, score, latent, tsquare, explained] = pca(x');
% new_f(1,:)= filt(1,:).*pc(1,1)+ filt(2,:).*pc(2,1);
new_f(1,:)= filt(1,:).*pc(1,1)+ filt(2,:).*pc(2,1);


new_e=abs(hilbert(new_f));

th=prctile(new_e,55);
id_sp1=find(new_e>th);
id_sp2=find(new_e<th);

idx=id_sp1;
    
d_idx=diff(idx);
pnts=find(d_idx>1);
begin=idx(pnts+1);
ending=idx(pnts);

begin2=[idx(1) begin];
ending2=[ending idx(end)];

dur=ending2-begin2;
findi=find(dur>samplerate);
% dur(find(dur>10))
ending3=ending2(findi);
begin3=begin2(findi);



time=1:length(new_e);

figure
% plot(time,new_f)
hold on
plot(time,new_e)
yline(th)
plot(time(begin3),new_e(begin3),'r.')
plot(time(ending3),new_e(ending3),'b.')


figure
plot(data(:,1),data(:,2),'Color',[0.5 0.5 0.5])
for i=1:length(begin3)
hold on
plot(data(begin3(i):ending3(i),1),data(begin3(i):ending3(i),2),'r.','MarkerSize',10)
% plot(data(id_sp2,1),data(id_sp2,2),'b.','MarkerSize',10)
end


time=1:length(new_e);

figure
plot(time,new_f)
hold on
plot(time,new_e)
yline(th)
plot(time(id_sp1),new_e(id_sp1),'r.')
plot(time(id_sp2),new_e(id_sp2),'b.')

data=data';
figure
plot3(time,data(:,1),data(:,2),'Color',[0.5 0.5 0.5])
hold on
plot3(time(id_sp1),data(id_sp1,1),data(id_sp1,2),'r.','MarkerSize',10)
plot3(time(id_sp2),data(id_sp2,1),data(id_sp2,2),'b.','MarkerSize',10)

figure
plot(data(:,1),data(:,2),'Color',[0.5 0.5 0.5])
hold on
plot(data(id_sp1,1),data(id_sp1,2),'r.','MarkerSize',10)
plot(data(id_sp2,1),data(id_sp2,2),'b.','MarkerSize',10)




