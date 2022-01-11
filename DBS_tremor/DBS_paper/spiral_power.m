clear
close
cohort=[ 1 3 4 6];
iii=cohort(3);
trial=1;
close all

%      load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P0',num2str(iii),'_NS',num2str(trial),'_SH.mat'))


%       load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P0',num2str(iii),'_HF',num2str(trial),'_SH.mat'))

load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P0',num2str(iii),'_C',num2str(trial),'_SH.mat'))
data=(NS1)';
samplerate=floor(1000/median(diff(data(3,:))));
time=0:1/samplerate:((size(data,2)-1)/samplerate);
data(find(data(1,:)==-1),:)=[];

% plot(data(1,:),data(2,:))
% close
% data(find(data(1,:)==-1),:)=[];
% plot(data(1,:),data(2,:),'Color',[0.8 0.5 0.8],'LineWidth',2)


%%% interpolate

samplerate2=100;
res=0:1/samplerate2:((size(data,2)-1)/samplerate);


for o=1:2
    data1(o,:)=interp1(time,data(o,:),res);
end


% plot(data1(1,:),data1(2,:),'r.')
% hold on
% plot(data(1,:),data(2,:),'k.');


%%% in radian
%  data1=data;
%  samplerate2=samplerate;
%  res=time;





centre=[535 361];
d1=data1(1,:)-centre(1);
d2=data1(2,:)-centre(2);
[px,py]=cart2pol(d1,d2);
polarplot(px,py)

com_sig=sqrt((data1(1,:).^2)+(data1(2,:).^2));
[a,b]=  butter(2, [2/(0.5*samplerate2) 8/(0.5*samplerate2)], 'bandpass'); %15
filt_cs=(filtfilt(a,b,com_sig));
env=abs(hilbert(filt_cs));


% %%% testing
one=0.361*60*samplerate2;

% plot(data(1,1:one),data(2,1:one),'Color',[0.8 0.5 0.8],'LineWidth',2)
% hold on
% plot(data(1,st:en),data(2,st:en),'k.')

%%% spectrogram
sg = samplerate2;
ov = samplerate2/2;
frq=[];

close
[s,f,t,p] =spectrogram(filt_cs,sg,ov,frq,samplerate2,'yaxis');

idx=(find(f>2));
cut=idx(1); clear idx
idx=(find(f<8));
cut(1,2)=idx(end); clear idx

s=abs(s);
pow=mean(s(cut(1):cut(2),:));



% ar=[];
% ar(1:find(time==(t(1))))=pow(1);
ar2=[];
for ii=1:length(t)-1
    % %     if ~isempty(find(time==(t(ii))):find(time==(t(ii+1))))
    % %         ar2(find(time==(t(ii))):find(time==(t(ii+1))))=pow(ii+1);
    % %     else
    ar2(find(res==(t(ii))):find(res==(t(ii+1))))=pow(ii+1);
    % %     end
end

dum=find(ar2==0);
ar2(1:dum(end))=pow(1);
ar3=repmat(pow(end),1,length(res(find(res==(t(end))):end))-1);
arr=[ar2 ar3];
arra=smoothdata(arr,'movmean',10);


figure
subplot(2,1,1)
plot(res,arra)
subplot(2,1,2)
plot(res,env)

figure
spectrogram(filt_cs,sg,ov,[],samplerate2,'yaxis');
colorbar
colormap jet
figure
polarscatter(px(1:one), py(1:one), [], arra(1:one),'filled')
colorbar
colormap jet
figure
polarscatter(px, py, [], arra,'filled')
colorbar
colormap jet

% scatter(data1(1,1:one), data1(2,1:one), [], arr(1:one),'filled')

