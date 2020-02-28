% 
% clear all
% close all
% %  iii=[1 2 3 4 5 8 10 11 12 13 16 17 18];
% iii=[2 3 4 5 8 10 11 13 16 17];
% 
% gg=[];
% main=[1 1 3 1 3 3 3 3 1 1];
% for numb= 1:length(iii);
%     %     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\reference_pls_ns.mat')
%     clearvars -except iii numb main amp_bins amp_n_bins
%     %     load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
%     load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(iii(numb)),'_baseline.mat'))
%     
%     data=SmrData.WvData;
%     samplerateold=SmrData.SR;
%     ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
%     ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
%     ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
%     samplerate=1000;
%     tt=0:1/samplerate:(size(ds_data,2)-1)/samplerate;
%     tre_3=ds_data([3 5 6],:); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
%     
%     before_ns
%     segmentb=hu{numb,:};
%     segmente=hd{numb,:};
%     
%     handup=[];
%     for i=1:length(segmentb)
%         handup=[handup segmentb(i):segmente(i)]; %#ok<*AGROW>
%     end
%     
%     clear i
%     handup=sort(handup,'ascend');
%     for aa=1:3
%         [Pxx,F]=pwelch(tre_3(aa,handup),samplerate,[],samplerate,samplerate);
%         frange=F(3:10);
%         Pxxrange=Pxx(3:10);
%         Freqpeak(aa,:)=frange(find(Pxxrange==max(Pxxrange)));
%         Ppeak(aa,:)=max(Pxxrange);
%         ps_curves(aa,:)=Pxx;
%     end
%     peak_ax=[(Freqpeak(find(Ppeak==max(Ppeak)))) (find(Ppeak==max(Ppeak)))];
%     Fpeak=peak_ax(1);
%     if (Fpeak-2)>=1
%         [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
%     else
%         [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
%     end
%     tf_3=filtfilt(b,a,tre_3')*10*9.81/0.5;
%     tf_3=tf_3';
%     
%     for i=1:3
%         env(i,:)=abs(hilbert(tf_3(i,:)));
%         phase(i,:)=angle(hilbert(tf_3(i,:)));
%         freq(i,:)=(smooth((1000/(2*pi))*diff(unwrap(phase(i,:))),500))';
%     end
%     
%     ns_amp=env(main(numb),handup);
%     ns_freq=freq(main(numb),handup);
%     %     [min(ns_freq) max(ns_freq)]
%     amp_bined=[];
%     dum2=[];
%     m_amp_b=[];
%     bins=2:0.5:9;
%     for i=1:length(bins)
%         if i+1<length(bins)
%             dum=find(ns_freq>bins(i) & ns_freq<=bins(i+1));
%             dum2=[dum2 ; bins(i).* ones(length(dum),1)];
%             amp_bined=[amp_bined  ((ns_amp(dum)-nanmean(ns_amp))./nanmean(ns_amp))];
%             m_amp_b(i,:)= nanmean(ns_amp(dum));
% %             m_n_amp(i,:)=nanmean((ns_amp(dum)-nanmean(ns_amp))./nanmean(ns_amp));
%             m_n_amp(i,:)=nanmedian(ns_amp(dum)./nanmedian(ns_amp));
%         end
%         
%         %         figure(1)
%         %         subplot(2,5,numb)
%         %         bar(m_amp_b);
%     end
%     
%     amp_bins(numb,:)=m_amp_b;
%     amp_n_bins(numb,:)=m_n_amp;
% end
% 
clearvars -except amp_bins amp_n_bins bins
load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','stone');
cl=stone;

for i=1:size(amp_n_bins,1)
    f1=figure(2)
    subplot(2,5,i)
    y=amp_n_bins(i,:);
    x=bins(1:end-2);
    bar(x,y,'FaceColor',cl,'EdgeColor',cl)
    hold on
    [rsg,rsg_g,rsg_o]=gauss_fit2(x,y)
    ylim([0 1.5])
    xticks([1:2:14])
%     xticklabels({'2','3','4','5','6','7','8'});
    ylabel({'change in ';'tremor severity'})
    xlabel('frequency(Hz)')
%     set(gca,'FontSize',12)
    box('off')
    legend('off')
    cv(i,:)=rsg.c./rsg.b;
    clear y rsg rsg_o rsg_g
end
f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 25, 10];
set(f1,'color','w');

load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\cleaned_rc12_noaddon.mat')
for i=1:size(tt1,1)
    sum_change(i,1)=sum(abs(nanmedian(tt1{i,1})))*100;
    mag_change(i,1)=(max(nanmedian(tt1{i,1}))-min(nanmedian(tt1{i,1})));
    sup_change(i,1)=min(nanmedian(tt1{i,1}))*100;
end

sig=[2 4 5 10];

f1=figure(3)
subplot(1,3,1)
y2=plot(cv,sum_change,'k.','MarkerSize',10);
y3=lsline;
hold on
y2=plot(cv(sig),sum_change(sig),'ko','MarkerSize',5);
set(y3,'LineWidth',1.5,'Color','red')
box('off')
[r,p]=corrcoef(cv,sum_change)
legend(y3,['r= ', num2str(r(1,2))],'box','off')
xlabel('coefficient of variation')
ylabel('total modulatory effect (%)')
set(gca,'FontSize',12)
f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 8, 10];
set(f1,'color','w');

subplot(1,3,2)
y2=plot(cv,mag_change,'k.','MarkerSize',10);
y3=lsline;
hold on
y2=plot(cv(sig),mag_change(sig),'ko','MarkerSize',5);
set(y3,'LineWidth',1.5,'Color','red')
box('off')
[r,p]=corrcoef(cv,mag_change)
legend(y3,['r= ', num2str(r(1,2))],'box','off')
xlabel('coefficient of variation')
ylabel('range of modulation')
set(gca,'FontSize',12)

subplot(1,3,3)
y2=plot(cv,sup_change,'k.','MarkerSize',10);
y3=lsline;
hold on
y2=plot(cv(sig),sup_change(sig),'ko','MarkerSize',5);
set(y3,'LineWidth',1.5,'Color','red')
box('off')
[r,p]=corrcoef(cv,sup_change)
legend(y3,['r= ', num2str(r(1,2))],'box','off')
xlabel('coefficient of variation')
ylabel('supressive modulatory effect (%)')
set(gca,'FontSize',12)

f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 25, 10];
set(f1,'color','w');

% sup=sup_change([2:3 5:10]);
% cv_sup=cv([2:3 5:10]);
% f1=figure(5)
% y2=plot(cv_sup,sup,'k.','MarkerSize',10);
% y3=lsline;
% set(y3,'LineWidth',1.5,'Color','red')
% box('off')
% [r,p]=corrcoef(cv_sup,sup)
% legend(y3,num2str(r(1,2)),'box','off')
% xlabel('coefficient of variation')
% ylabel('supressive effect only (%)')
% set(gca,'FontSize',12)
% f1.Units = 'centimeters';
% f1.OuterPosition= [10, 10, 8, 10];
% set(f1,'color','w');

% clearvars -except iii numb ns_filt_mainax in2
% clear numb in2