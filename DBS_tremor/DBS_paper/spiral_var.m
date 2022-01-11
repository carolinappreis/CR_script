

clear all
close all
cond={'NS';'HF';'C'};
cohort=[ 1 3 4 6];

iii=3;
trial=1;
co=2;

load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/P0',num2str(cohort(iii)),'_clean_',num2str(cond{co,1}),num2str(trial),'_SH.mat'));

%%% filt & var
[a,b]=  butter(2, [2/(0.5*samplerate2) 8/(0.5*samplerate2)], 'bandpass'); %15
filt_cs=(filtfilt(a,b,signal(3,:)));
env=abs(hilbert(filt_cs));
pp=peak2peak(filt_cs);

%%% spatial
% signal=data;
centre=[535 361];
d1=signal(1,:)-centre(1);
d2=signal(2,:)-centre(2);
[px,py]=cart2pol(d1,d2);

figure(1)
polarplot(px,py)
hold on
figure(2)
plot(tempo,filt_cs)
hold on


pxx=rad2deg(px);
idx{1,1}=find(pxx>0 & pxx< 90);
idx{2,1}=find(pxx<0 & pxx>-90);
idx{3,1}=find(pxx>-180 & pxx<-90);
idx{4,1}=find(pxx>90 & pxx<180);


t_spi=cell(size(cohort,2),2,size(cond,1));

% t_spi{2,1,1}=[[1 1939];[1960 4304];];
t_spi{2,1,1}=[1 length(tempo)];
t_spi{2,2,1}=[[1 2534];[2534 5228];[5228 length(tempo)]];
% t_spi{2,1,2}=[[1 1161];[1189 3040]; [3040 4862]];
t_spi{2,1,2}=[[1 1161];[1189 3040]];
t_spi{2,1,3}=[[47 2105];[2105 4072];[4072 6236]];

t_spi{2,1,2}=[[1 1892];[1894 length(tempo)]];
t_spi{2,1,3}=[[47 2105];[2105 4072];[4072 6236]];


t_spi{3,1,1}=[[1 1687];[1687 2947];[2947 4354];[4354 5601];[5650 length(tempo)]];
t_spi{3,2,1}=[[1 1668];[1668 2979];[2979 4341];[4341 length(tempo)]];
t_spi{3,1,2}=[[1 1130];[1130 2517];[2517 4207];[4207 5605];[5605 length(tempo)]];
t_spi{3,1,3}=[[1 1420];[1420 2647];[2647 3834];[3834 4923];[4923 6132];[6132 7283];[7283 8507];[8507 9827];[9827 10990];[10990 length(tempo)]];

t_spi{4,1,1}=[[1 6260];[7209 length(tempo)]];
t_spi{4,1,2}=[[1 5236];[5236 length(tempo)]];
t_spi{4,1,3}=[[1 5932];[5932 length(tempo)]];

figure(4)
senv=smoothdata(env,'movmean',samplerate2);
polarscatter(px(:), py(:), [], senv(:),'filled')
colorbar
colormap winter
title('envelope')

figure(5)
for mm=1:size(t_spi{iii,trial,co}(:,1),1)
    
    subplot(1,size(t_spi{iii,trial,co}(:,1),1),mm)
    polarscatter(px(t_spi{iii,trial,co}(mm,1):t_spi{iii,trial,co}(mm,2)),py(t_spi{iii,trial,co}(mm,1):t_spi{iii,trial,co}(mm,2)),[],senv(t_spi{iii,trial,co}(mm,1):t_spi{iii,trial,co}(mm,2)),'filled');
    %colorbar
    colormap winter
    caxis([2 8])
end

for i=1:4 
    figure(1)
    polarplot(px(idx{i,1}),py(idx{i,1}),'.')
    hold on
    figure(2)
    plot(tempo(idx{i,1}),filt_cs(idx{i,1}),'.')
    
    mpipi(1,i) = peak2peak(filt_cs(idx{i,1}));
    mvar(1,i) = var(filt_cs(idx{i,1}));
    menv(1,i)= mean(env(idx{i,1}));  
end



%%% all data
f2=figure()
dat_p1=[menv;mvar;mpipi];
for a=1:3
    subplot(1,3,a)
    x=squeeze(dat_p1(a,:));
    bar([x],'EdgeColor','none','FaceColor',[0.5 0.5 1],'FaceAlpha',0.5)
    hold on
    % ylim([0 500])
    % yticks(0:100:500)
    xticklabels({'1Q','2Q','3Q','4Q'})
    if a==1
        ylabel('mean envelope')
    elseif a==2
        ylabel('mean variance')
    else
        ylabel('mean peak to peak (fx)')
    end
    box('off')
    set(gca,'FontSize',14)
    clear x err
end
f2.OuterPosition= [1,100,1000,300];
set(f2,'color','w');


% % clear all
% % close all
% % cond={'NS';'HF';'C'};
% % cohort=[ 1 3 4 6];
% % 
% % iii=3;
% % trial=1;
% % co=2;
% % 
% % load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/P0',num2str(cohort(iii)),'_clean_',num2str(cond{co,1}),num2str(trial),'_SH.mat'));
% % 
% % %%% filt & var
% % com_sig=sqrt((signal(1,:).^2)+(signal(2,:).^2));
% % [a,b]=  butter(2, [2/(0.5*samplerate2) 8/(0.5*samplerate2)], 'bandpass'); %15
% % filt_cs=(filtfilt(a,b,com_sig));
% % shift=abs(min(filt_cs)-1);
% % 
% % env=abs(hilbert(filt_cs));
% % pp=peak2peak(filt_cs);
% % 
% % 
% % %%% spectrogram
% % sg = samplerate2;
% % ov = samplerate2/2;
% % frq=[];
% % 
% % [s,f,t,p]=spectrogram(filt_cs,sg,ov,frq,samplerate2,'yaxis');
% % 
% % idx=(find(f>2));
% % cut=idx(1); clear idx
% % idx=(find(f<8));
% % cut(1,2)=idx(end); clear idx
% % 
% % s=abs(s);
% % pow=mean(s(cut(1):cut(2),:));
% % ar2=[];
% % for ii=1:length(t)-1
% %     ar2(find(tempo==(t(ii))):find(tempo==(t(ii+1))))=pow(ii+1);
% % end
% % 
% % dum=find(ar2==0);
% % ar2(1:dum(end))=pow(1);
% % ar3=repmat(pow(end),1,length(tempo(find(tempo==(t(end))):end))-1);
% % arr=[ar2 ar3];
% % 
% % %%% spatial
% % % signal=data;
% % centre=[535 361];
% % d1=signal(1,:)-centre(1);
% % d2=signal(2,:)-centre(2);
% % [px,py]=cart2pol(d1,d2);
% % 
% % pxx=rad2deg(px);
% % idx{1,1}=find(pxx>0 & pxx< 90);
% % idx{2,1}=find(pxx<0 & pxx>-90);
% % idx{3,1}=find(pxx>-180 & pxx<-90);
% % idx{4,1}=find(pxx>90 & pxx<180);
% % 
% % 
% % t_spi=cell(size(cohort,2),2,size(cond,1));
% % 
% % % t_spi{2,1,1}=[[1 1939];[1960 4304];];
% % t_spi{2,1,1}=[1 length(tempo)];
% % t_spi{2,2,1}=[[1 2534];[2534 5228];[5228 length(tempo)]];
% % % t_spi{2,1,2}=[[1 1161];[1189 3040]; [3040 4862]];
% % t_spi{2,1,2}=[[1 1161];[1189 3040]];
% % t_spi{2,1,3}=[[47 2105];[2105 4072];[4072 6236]];
% % 
% % t_spi{2,1,2}=[[1 1892];[1894 length(tempo)]];
% % t_spi{2,1,3}=[[47 2105];[2105 4072];[4072 6236]];
% % 
% % 
% % t_spi{3,1,1}=[[1 1687];[1687 2947];[2947 4354];[4354 5601];[5650 length(tempo)]];
% % t_spi{3,2,1}=[[1 1668];[1668 2979];[2979 4341];[4341 length(tempo)]];
% % t_spi{3,1,2}=[[1 1130];[1130 2517];[2517 4207];[4207 5605];[5605 length(tempo)]];
% % t_spi{3,1,3}=[[1 1420];[1420 2647];[2647 3834];[3834 4923];[4923 6132];[6132 7283];[7283 8507];[8507 9827];[9827 10990];[10990 length(tempo)]];
% % 
% % t_spi{4,1,1}=[[1 6260];[7209 length(tempo)]];
% % t_spi{4,1,2}=[[1 5236];[5236 length(tempo)]];
% % t_spi{4,1,3}=[[1 5932];[5932 length(tempo)]];
% % 
% % 
% % figure(1)
% % polarplot(px,py)
% % hold on
% % figure(2)
% % plot(tempo,filt_cs)
% % hold on
% % 
% % figure(3)
% % subplot(2,1,1)
% % plot(tempo,arr)
% % subplot(2,1,2)
% % plot(tempo,env)
% % 
% % 
% % figure(4)
% % subplot(1,2,1)
% % sarr=smoothdata(arr,'movmean',100);
% % polarscatter(px(:), py(:), [], sarr(:),'filled')
% % title('power')
% % colorbar
% % colormap jet
% % subplot(1,2,2)
% % senv=smoothdata(env,'movmean',100);
% % polarscatter(px(:), py(:), [], senv(:),'filled')
% % colorbar
% % colormap jet
% % title('envelope')
% % 
% % 
% % figure(5)
% % for mm=1:size(t_spi{iii,trial,co}(:,1),1)
% %     
% %     subplot(1,size(t_spi{iii,trial,co}(:,1),1),mm)
% %     polarscatter(px(t_spi{iii,trial,co}(mm,1):t_spi{iii,trial,co}(mm,2)),py(t_spi{iii,trial,co}(mm,1):t_spi{iii,trial,co}(mm,2)),[],senv(t_spi{iii,trial,co}(mm,1):t_spi{iii,trial,co}(mm,2)),'filled');
% %     %colorbar
% %     colormap winter
% %     caxis([2 8])
% % end
% % 
% % for i=1:4
% %     
% %     ddum=diff(idx{i,1});
% %     nid=idx{i,1}(find(ddum>1));
% %     
% %     figure(1)
% %     polarplot(px(idx{i,1}),py(idx{i,1}),'.')
% %     figure(2)
% %     plot(tempo(idx{i,1}),filt_cs(idx{i,1}),'.')
% %     
% %     
% %     mpipi(1,i) = peak2peak(filt_cs(idx{i,1}));
% %     mvar(1,i) = var(filt_cs(idx{i,1}));
% %     mpower(1,i)=mean(arr(idx{i,1}));
% %     menv(1,i)= mean(env(idx{i,1}));
% %     
% %     for u=1:(length(nid))-1
% %         rr=(nid(u):nid(u+1)-1);
% %         ll{i,1}(1,u)=length(rr);
% %         seg_menv{i,1}(1,u)=mean(env(rr));
% %         seg_var{i,1}(1,u)=var(filt_cs(rr));
% %         seg_pp{i,1}(1,u)=peak2peak(filt_cs(rr));
% %         seg_power{i,1}(1,u)=mean(arr(rr));
% %         %       clear r maxvalM maxidxM
% %     end
% %     clear ddum nid
% % end
% % 
% % 
% % for i=1:4
% %     dat_p(1,i,:)=[mean(seg_menv{i,1}) std(seg_menv{i,1})];
% %     dat_p(2,i,:)=[mean(seg_power{i,1}) std(seg_power{i,1})];
% %     dat_p(3,i,:)= [mean(seg_var{i,1}) std(seg_var{i,1})];
% %     dat_p(4,i,:)=[mean(seg_pp{i,1}) std(seg_pp{i,1})];
% % end
% % 
% % 
% % %%% from segments
% % f1=figure()
% % for a=1:4
% %     subplot(1,4,a)
% %     x=squeeze(dat_p(a,:,1));
% %     err=squeeze(dat_p(a,:,2));
% %     bar([x],'EdgeColor','none','FaceColor',[0.5 0.5 1],'FaceAlpha',0.5)
% %     hold on
% %     errorbar([1:4],x,err,'.','LineWidth',2)
% %     % ylim([0 500])
% %     % yticks(0:100:500)
% %     xticklabels({'1Q','2Q','3Q','4Q'})
% %     if a==1
% %         ylabel('mean envelope')
% %     elseif a==2
% %         ylabel('mean spectral power')
% %     elseif a==3
% %         ylabel('mean variance')
% %     else a==4
% %         ylabel('mean peak to peak (fx)')
% %     end
% %     
% %     
% %     box('off')
% %     set(gca,'FontSize',14)
% %     clear x err
% % end
% % f1.OuterPosition= [1,100,1000,300];
% % set(f1,'color','w');
% % 
% % 
% % %%% all data
% % f2=figure()
% % dat_p1=[menv;mpower;mvar;mpipi];
% % for a=1:4
% %     subplot(1,4,a)
% %     x=squeeze(dat_p1(a,:));
% %     bar([x],'EdgeColor','none','FaceColor',[0.5 0.5 1],'FaceAlpha',0.5)
% %     hold on
% %     % ylim([0 500])
% %     % yticks(0:100:500)
% %     xticklabels({'1Q','2Q','3Q','4Q'})
% %     if a==1
% %         ylabel('mean envelope')
% %     elseif a==2
% %         ylabel('mean spectral power')
% %     elseif a==3
% %         ylabel('mean variance')
% %     else
% %         ylabel('mean peak to peak (fx)')
% %     end
% %     box('off')
% %     set(gca,'FontSize',14)
% %     clear x err
% % end
% % f2.OuterPosition= [1,100,1000,300];
% % set(f2,'color','w');






% % % clear
% % % close
% % % cohort=[ 1 3 4 6];
% % % iii=cohort(3);
% % % trial=1;
% % % close all
% % %
% % % %       load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P0',num2str(iii),'_NS',num2str(trial),'_SH.mat'))
% % %
% % %
% % % %       load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P0',num2str(iii),'_HF',num2str(trial),'_SH.mat'))
% % %
% % % load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P0',num2str(iii),'_C',num2str(trial),'_SH.mat'))
% % % data=(NS1)';
% % % samplerate=floor(1000/median(diff(data(3,:))));
% % % time=0:1/samplerate:((size(data,2)-1)/samplerate);
% % %
% % % plot(data(1,:),data(2,:))
% % % close
% % % data(find(data(1,:)==-1),:)=[];
% % % plot(data(1,:),data(2,:),'Color',[0.8 0.5 0.8],'LineWidth',2)
% % %
% % %
% % %
% % % %%%interpolate
% % % samplerate2=100;
% % % res=0:1/samplerate2:((size(data,2)-1)/samplerate);
% % %
% % %
% % % for o=1:2
% % %     data1(o,:)=interp1(time,data(o,:),res);
% % % end
% % %
% % %
% % % % plot(data1(1,:),data1(2,:),'r.')
% % % % hold on
% % % % plot(data(1,:),data(2,:),'k.');
% % %
% % % %%% filt & var
% % % com_sig=sqrt((data1(1,:).^2)+(data1(2,:).^2));
% % % [a,b]=  butter(2, [2/(0.5*samplerate2) 8/(0.5*samplerate2)], 'bandpass'); %15
% % % filt_cs=(filtfilt(a,b,com_sig));
% % % shift=abs(min(filt_cs)-1);
% % % filt_cs=filt_cs+shift;
% % %
% % % env=abs(hilbert(filt_cs));
% % % pp=peak2peak(filt_cs);
% % %
% % %
% % % %%% spectrogram
% % % sg = samplerate2;
% % % ov = samplerate2/2;
% % % frq=[];
% % %
% % % [s,f,t,p] =spectrogram(filt_cs,sg,ov,frq,samplerate2,'yaxis');
% % %
% % % idx=(find(f>2));
% % % cut=idx(1); clear idx
% % % idx=(find(f<8));
% % % cut(1,2)=idx(end); clear idx
% % %
% % % s=abs(s);
% % % pow=mean(s(cut(1):cut(2),:));
% % % ar2=[];
% % % for ii=1:length(t)-1
% % %     ar2(find(res==(t(ii))):find(res==(t(ii+1))))=pow(ii+1);
% % % end
% % %
% % % dum=find(ar2==0);
% % % ar2(1:dum(end))=pow(1);
% % % ar3=repmat(pow(end),1,length(res(find(res==(t(end))):end))-1);
% % % arr=[ar2 ar3];
% % %
% % % %%% spatial
% % % % data1=data;
% % % centre=[535 361];
% % % d1=data1(1,:)-centre(1);
% % % d2=data1(2,:)-centre(2);
% % % [px,py]=cart2pol(d1,d2);
% % %
% % % pxx=rad2deg(px);
% % % idx{1,1}=find(pxx>0 & pxx< 90);
% % % idx{2,1}=find(pxx<0 & pxx>-90);
% % % idx{3,1}=find(pxx>-180 & pxx<-90);
% % % idx{4,1}=find(pxx>90 & pxx<180);
% % %
% % %
% % % close
% % % mimax=[];
% % % for i=1:200:length(filt_cs)
% % %     array=[];
% % %     array1=[];
% % %     if i+200<length(filt_cs)
% % %         bin=200;
% % %     else
% % %         bin=length(filt_cs)-i;
% % %     end
% % %     array=i:i+bin;
% % %     dum=filt_cs(array);
% % %     [maxval,maxidx]=findpeaks(dum,'MinPeakDistance',25);
% % %     [minval,minidx]=findpeaks(-(dum),'MinPeakDistance',25);
% % %
% % %     amax=array(maxidx);
% % %     amin=array(minidx);
% % %
% % %     if length(amax)~=length(amin)
% % %         if length(amax)<length(amin)
% % %             array1(1,:)=amax;
% % %             array1(2,:)=amin(1:length(amax));
% % %         elseif length(amax)>length(amin)
% % %             array1(1,:)=amax(1:length(amin));
% % %             array1(2,:)=amin;
% % %         end
% % %     else
% % %         array1(1,:)=amax;
% % %         array1(2,:)=amin;
% % %     end
% % %
% % %
% % %     %%% check
% % %     %         plot(res(array),filt_cs(array))
% % %     %         hold on
% % %     %         plot(res(array1(1,:)),filt_cs(array1(1,:)),'r.');
% % %     %         plot(res(array1(2,:)),filt_cs(array1(2,:)),'k.');
% % %     close
% % %     clear amin amax maxval minval maxidx minidx
% % %     mypp(i:i+bin)=mean(abs(array1(1,:)-array1(2,:)));
% % % end
% % %
% % %
% % % figure(1)
% % % polarplot(px,py)
% % % hold on
% % % figure(2)
% % % plot(res,filt_cs)
% % % hold on
% % %
% % % figure(3)
% % % subplot(2,1,1)
% % % plot(res,arr)
% % % subplot(2,1,2)
% % % plot(res,env)
% % %
% % % one=0.325*samplerate2*60;
% % %
% % % figure(4)
% % % subplot(1,3,1)
% % % sarr=smoothdata(arr,'movmean',60);
% % % polarscatter(px(:), py(:), [], sarr(:),'filled')
% % % title('power')
% % % colorbar
% % % colormap jet
% % % subplot(1,3,2)
% % % senv=smoothdata(env,'movmean',60);
% % % polarscatter(px(:), py(:), [], senv(:),'filled')
% % % colorbar
% % % colormap jet
% % % title('envelope')
% % % subplot(1,3,3)
% % % smpp=smoothdata(mypp,'movmean',60);
% % % polarscatter(px(:), py(:), [], smpp(:),'filled')
% % % title('peak2peak (inh)')
% % % colorbar
% % % colormap jet
% % %
% % %
% % %
% % %
% % % for i=1:4
% % %
% % %     ddum=diff(idx{i,1});
% % %     nid=idx{i,1}(find(ddum>1));
% % %
% % %     figure(1)
% % %     polarplot(px(idx{i,1}),py(idx{i,1}),'.')
% % %     figure(2)
% % %     plot(res(idx{i,1}),filt_cs(idx{i,1}),'.')
% % %
% % %
% % %     mpipi(1,i) = peak2peak(filt_cs(idx{i,1}));
% % %     mvar(1,i) = var(filt_cs(idx{i,1}));
% % %     mpower(1,i)=mean(arr(idx{i,1}));
% % %     menv(1,i)= mean(env(idx{i,1}));
% % %     mmpp(1,i)= mean(mypp(idx{i,1}));
% % %
% % %     for u=1:(length(nid))-1
% % %         rr=(nid(u):nid(u+1)-1);
% % %         ll{i,1}(1,u)=length(rr);
% % %         seg_menv{i,1}(1,u)=mean(env(rr));
% % %         seg_var{i,1}(1,u)=var(filt_cs(rr));
% % %         seg_pp{i,1}(1,u)=peak2peak(filt_cs(rr));
% % %         seg_power{i,1}(1,u)=mean(arr(rr));
% % %         seg_mpp{i,1}(1,u)=mean(mypp(rr));
% % %         %       clear r maxvalM maxidxM
% % %     end
% % %     clear ddum nid
% % % end
% % %
% % %
% % % for i=1:4
% % %     dat_p(1,i,:)=[mean(seg_menv{i,1}) std(seg_menv{i,1})];
% % %     dat_p(2,i,:)=[mean(seg_power{i,1}) std(seg_power{i,1})];
% % %     dat_p(3,i,:)= [mean(seg_var{i,1}) std(seg_var{i,1})];
% % %     dat_p(4,i,:)=[mean(seg_pp{i,1}) std(seg_pp{i,1})];
% % %     dat_p(5,i,:)=[mean(seg_mpp{i,1}) std(seg_mpp{i,1})];
% % % end
% % %
% % % f1=figure()
% % % for a=1:5
% % %     subplot(1,5,a)
% % %     x=squeeze(dat_p(a,:,1));
% % %     err=squeeze(dat_p(a,:,2));
% % %     bar([x],'EdgeColor','none','FaceColor',[0.5 0.5 1],'FaceAlpha',0.5)
% % %     hold on
% % %     errorbar([1:4],x,err,'.','LineWidth',2)
% % %     % ylim([0 500])
% % %     % yticks(0:100:500)
% % %     xticklabels({'1Q','2Q','3Q','4Q'})
% % %     if a==1
% % %         ylabel('mean envelope')
% % %     elseif a==2
% % %         ylabel('mean spectral power')
% % %     elseif a==3
% % %         ylabel('mean variance')
% % %     elseif a==4
% % %         ylabel('mean peak to peak (fx)')
% % %     else
% % %         ylabel('mean peak to peak (inh)')
% % %     end
% % %
% % %
% % %     box('off')
% % %     set(gca,'FontSize',14)
% % %     clear x err
% % % end
% % % f1.OuterPosition= [1,100,1000,300];
% % % set(f1,'color','w');
% % %
% % % f2=figure()
% % % dat_p1=[menv;mpower;mvar;mpipi;mmpp];
% % % for a=1:5
% % %     subplot(1,5,a)
% % %     x=squeeze(dat_p1(a,:));
% % %     bar([x],'EdgeColor','none','FaceColor',[0.5 0.5 1],'FaceAlpha',0.5)
% % %     hold on
% % %     % ylim([0 500])
% % %     % yticks(0:100:500)
% % %     xticklabels({'1Q','2Q','3Q','4Q'})
% % %     if a==1
% % %         ylabel('mean envelope')
% % %     elseif a==2
% % %         ylabel('mean spectral power')
% % %     elseif a==3
% % %         ylabel('mean variance')
% % %     elseif a==4
% % %         ylabel('mean peak to peak (fx)')
% % %     else
% % %         ylabel('mean peak to peak (inh)')
% % %     end
% % %     box('off')
% % %     set(gca,'FontSize',14)
% % %     clear x err
% % % end
% % % f2.OuterPosition= [1,100,1000,300];
% % % set(f2,'color','w');
% % %
% % %
% % %
% % %
% % %
% % %
% % %
% % %
% % %
% % %
% % %
% % %
% % %
% % % %%% testing
% % % % one=0.3106*60*66;
% % % % st=3.7*66;
% % % % en=5.4*66;
% % % % %
% % % % plot(data(1,1:one),data(2,1:one),'Color',[0.8 0.5 0.8],'LineWidth',2)
% % % % hold on
% % % % plot(data(1,st:en),data(2,st:en),'k.')
% % %
% % %
