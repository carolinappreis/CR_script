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
close
cohort=[ 1 3 4 6];
iii=cohort(3);
trial=1;
close all

%     load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P0',num2str(iii),'_NS',num2str(trial),'_SH.mat'))


%      load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P0',num2str(iii),'_HF',num2str(trial),'_SH.mat'))

   load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P0',num2str(iii),'_C',num2str(trial),'_SH.mat'))
data=NS1';
samplerate=floor(1000/median(diff(data(3,:))));
samplerate2=samplerate*100;



time=1:length(data(3,:));
time1=1:0.01:length(data(3,:));


for o=1:2
data1(o,:)=interp1(data(o,:),time1,'linear','extrap');
end


g1=figure(1)
[Pxx,F] = pwelch(data1(2,:), 2*samplerate2, [], 2*samplerate2, 2*samplerate2);
plot(F,Pxx,'b','LineWidth',1.5)
xlim([0 50])
% ylim([0 35])
mm=2;MM=50;
Fp1=F(find(Pxx==max(Pxx(find(F==mm):find(F==MM)))));
clear Pxx F
hold on

[Pxx,F] = pwelch(data1(1,:), 2*samplerate2, [], 2*samplerate2, 2*samplerate2);
plot(F,Pxx,'r','LineWidth',1.5)
xlim([0 50])
% ylim([0 35])
Fp2=F(find(Pxx==max(Pxx(find(F==mm):find(F==MM)))));
legend('xaxis','yaxis')
box('off')
set(g1,'color','w');
g1.OuterPosition= [1,100,280,380];
legend('boxoff')
xlabel('Frequency (Hz)')
ylabel('Power (u.a)')



Fpeak=(mean([Fp1 Fp2]));
 

% Fpeak=4

% [a,b]=  butter(2, [(Fpeak-2)/(0.5*samplerate2) (Fpeak+2)/(0.5*samplerate2)], 'bandpass'); %15

[a,b]=  butter(2, [3/(0.5*samplerate2) 6/(0.5*samplerate2)], 'bandpass'); %15

for o=1:2
    filt(o,:)=(filtfilt(a,b,data1(o,:))).*(10*9.81/0.5);
end

x = [filt(1,:); filt(2,:)];
[pc, score, latent, tsquare, explained] = pca(x');
new_f(1,:)= filt(1,:).*pc(1,1)+ filt(2,:).*pc(2,1);

% 
subplot(3,1,1)
plot(time1,data1(2,:),'.')
hold on
plot(time,data(2,:),'.')
subplot(3,1,2)
plot(time1,filt(2,:))
subplot(3,1,3)
plot(time1,new_f)

new_e=abs(hilbert(new_f));

th=prctile(new_e,10);
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
findi=find(dur>samplerate2*2);
% dur(find(dur>10))
ending3=ending2(findi);
begin3=begin2(findi);

m1=figure
plot(time1,new_f)
hold on
plot(time1,new_e)
yline(th)
plot(time1(begin3),new_e(begin3),'r.','LineWidth',3)
plot(time1(ending3),new_e(ending3),'b.','LineWidth',3)
xticks([0:samplerate2./10:floor(length(new_f))])
xticklabels([0:10:170])
ylabel('Amplitude(uV)')
xlabel('Time (s)')
box('off')
set(m1,'color','w');
m1.OuterPosition= [1,100,1000,300];
% % 

figure
plot(time1,data1(2,:),'Color',[0.5 0.5 0.5])
for i=1:length(begin3)
    hold on
    plot(time1(begin3(i):ending3(i)),data1(2,begin3(i):ending3(i)),'r.','MarkerSize',10)
    % plot(data(id_sp2,1),data(id_sp2,2),'b.','MarkerSize',10)
end
    
figure
plot(data1(1,:),data1(2,:),'Color',[0.5 0.5 0.5])
for i=1:length(begin3)
hold on
plot(data1(1,begin3(i):ending3(i)),data1(2,begin3(i):ending3(i)),'r.','MarkerSize',5)
% plot(data(id_sp2,1),data(id_sp2,2),'b.','MarkerSize',10)
end

figure
plot(time1(1,:),data1(1,:),'Color',[0.5 0.5 0.5])
for i=1:length(begin3)
hold on
plot(time1(1,begin3(i):ending3(i)),data1(1,begin3(i):ending3(i)),'r.','MarkerSize',5)
% plot(data(id_sp2,1),data(id_sp2,2),'b.','MarkerSize',10)
end
xticks([0:samplerate2./10:floor(length(new_f))])
xticklabels([0:10:170])
ylabel('Amplitude(uV)')
xlabel('Time (s)')
box('off')
set(m1,'color','w');
m1.OuterPosition= [1,100,1000,300];


figure
plot3(time1,data1(1,:),data1(2,:),'Color',[0.5 0.5 0.5])
hold on
for i=1:length(begin3)
hold on
plot3(time1(1,begin3(i):ending3(i)),data1(1,begin3(i):ending3(i)),data1(2,begin3(i):ending3(i)),'r.','MarkerSize',5)
end


f3=figure(15);
[S,F,T,P] = spectrogram(data1(1,:),hanning(512),125,3:50,samplerate2);
imagesc(T,fliplr(F),P);
set(gca,'FontSize',12);
ax = gca;
ax.YDir = 'normal'
xticks(0:10:160)
yticks(5:10:50)
ylabel('Frequency (Hz)')
xlabel('Time (s)')
box('off')
colorbar
set(f3,'color','w');
f3.OuterPosition= [1,100,1000,300];
xlim([0 90])


