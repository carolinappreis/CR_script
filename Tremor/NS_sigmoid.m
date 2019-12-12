clear all
close all
iii=[1 2 3 4 5 8 10 11 12 13 16 17 18];
in2=1;
for numb=1:length(iii);
    clearvars -except iii numb NS NS_i in2 yy psd_curves peaks prm
%        load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
    load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(iii(numb)),'_baseline.mat'))
    
    if in2==1
        in=3;
    elseif in2==2 % other axis 1
        in=5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    elseif in2==3 % other axis 2
        in=6; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    end
    data=SmrData.WvData;
    samplerateold=SmrData.SR;
    tremor=(data(in,:));
    addon=92; addon_end=35;
    
    ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    tt=0:1/samplerate:(size(ds_data,2)-1)/samplerate;
    tre_3=ds_data([3 5 6],:); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    
    hu={round(([7.91 9.731])*(10^5),1);round([1.598 3.55 5.142 6.585 8.159]*(10^5),1);round([ 1.072 2.272 3.485 4.681]*(10^5),1)...
        ;round([1.116 2.315 4.253 5.885 7.026]*(10^5),1);round([5.257 7.746 9.313 1.069*10 1.183*10]*(10^5),1);...
        round([2.159/10 1.4600 2.9230 4.3520 5.6510]*(10^5),1);round([1.47 2.754 4.232 5.503]*(10^5),1);round([2.4/10 1.317 3.012]*(10^5),1);...
        round([ 134450 340484 466600 639582],1);round([139104 277397 401396 551971],1);round(([16972 158835 388948 582772 734164]),1);round([7.175/10 2.61 4.43 ]*(10^5),1);round([2.474 4.46]*(10^5),1)};
    
    hd={round(([8.958 1.097*10])*(10^5),1);round([2.43 4.486 6.051 7.509 9.037]*(10^5),1);round([ 1.774 2.905 4.186 5.348]*(10^5),1)...
        ;round([1.698 3.091 5.132 6.481 7.825]*(10^5),1);round([6.386 8.465 9.991 1.131*10 1.24*10]*(10^5),1);...
        round([8.309/10 2.1070 3.5590 4.9500 6.3050]*(10^5),1);round([ 2.104 3.401 4.898 6.207]*(10^5),1);round([9.98/10 2.385 4.394]*(10^5),1);...
        round([ 230271 431076 548826 747772],1);round([211850 343728 477800 622380],1);round(([76448 254875 504159 650313 814627]),1);round([1.296 3.362 5.071 ]*(10^5),1);round([2.959 5.272]*(10^5),1)};
    
    
    segmentb=hu{numb,:};
    segmente=hd{numb,:};
    
    handup=[];
    for i=1:length(segmentb)
        handup=[handup segmentb(i):segmente(i)]; %#ok<*AGROW>
    end
    clear i
    handup=sort(handup,'ascend');
    
    for aa=1:3
        [Pxx,F]=pwelch(tre_3(aa,handup),samplerate,[],samplerate,samplerate);
        frange=F(3:10);
        Pxxrange=Pxx(3:10);
        Freqpeak(aa,:)=frange(find(Pxxrange==max(Pxxrange)));
        Ppeak(aa,:)=max(Pxxrange);
        ps_curves(aa,:)=Pxx;
    end
    peak_ax=[(Freqpeak(find(Ppeak==max(Ppeak)))) (find(Ppeak==max(Ppeak)))];
    Fpeak=peak_ax(1);
    
%         figure()
%         plot(F(3:50),ps_curves(:,3:50)','LineWidth',2)
%         legend({'x','y','z'})
%         legend('boxoff')
%         box('off')
%     
    if (Fpeak-2)>=1
        [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    else
        [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    end
    tf_3=filtfilt(b,a,tre_3')*10*9.81/0.5;
    % tremor_or=zscore(tremor_or);
    dummy=(hilbert(tf_3))';
    envelope=abs(dummy);
    zenv=(abs(hilbert(zscore(tf_3))))';
    phase=angle(dummy);
    frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
    
    
    [d,e]=butter(2,[0.5/(0.5*samplerate) ],'low'); %15
    C=(filtfilt(d,e,tre_3'));
    figure()
    for jj=1:3
        subplot(3,1,jj)
        plot(zscore(tf_3(:,jj)))
        hold on
        plot(zscore(C(:,jj)))
        for i=1:size(segmentb,2)
            xline(segmentb(i),'r')
            xline(segmente(i),'k')
        end
        box('off')
    end
%     
    
    if min(segmentb)<20000
        segmentbp=segmentb(2:end);
        segmentep=segmente(2:end);
    else
        segmentbp=segmentb;
        segmentep=segmente;
    end
    
    %         ini=15000;
    ini=20000;
    for ax=1:3
        gg=[];
        for i=1:length(segmentbp)
            gg=[gg ;(zenv(ax,(segmentbp(i)-ini):(segmentbp(i)+60000)))];
            amp_start(ax,1,i)=mean(envelope(ax,segmentbp(i)-1000:segmentbp(i)));
            amp_end(ax,1,i)=mean(envelope(ax,segmentep(i)-1000:segmentep(i)));
        end
        
%         f5=figure(5)
%         subplot(3,1,ax)
%         plot(median(gg),'Color',[0.5 0.5 0.5])
%         hold on
%         xline(ini,'g--','LineWidth',2)
%         title(['NS pt',num2str(iii(numb))])
        
%                 yy(numb,:)=smooth(median(gg,1));
%                 y=yy(numb,1:40000);
%                 x=tt(1:length(y));
%                 initial_params=[];
%                 [param]=sigm_fit(x,y,initial_params)        % automatic initial_params
%                 clear x y
    end
    
%         for i=1:length(segmentb)
%             x=[sum(envelope(1,segmentb(i):segmente(i)));sum(envelope(2,segmentb(i):segmente(i)));sum(envelope(3,segmentb(i):segmente(i)))];
%         end
%         f10=figure(10)
%         bar(x)
%         box('off')
%         f10.Units = 'centimeters';
%         f10.OuterPosition= [10, 10, 8, 8];
%         set(f10,'color','w');
% % 
%     
%         prm(numb,:)=param;
%         peaks(1,numb)=Fpeak;
%         psd_curves(numb,1:3,:)=ps_curves;
%     close all

timings=[median(amp_start,3) median(amp_end,3)];
f1=figure;
plot((timings'),'.','MarkerSize',20)
hold on
plot((timings'),'LineWidth',0.5)
xlim([0 4])
xticks([1 2 3])
xticklabels({'bef handUP','bef handDW'})
xtickangle(45)
ylabel('Tremor severity ')
legend('boxoff')
f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 10, 10];
set(gca,'FontSize',14)
set(f1,'color','w');
box('off')
end





% x=tt(1:length(yy));
% y=yy(1,:);
%
% [param]=sigm_fit(x,y,fixed_params)        % automatic initial_params
% [param]=sigm_fit(x,y,[],initial_params)   % use it when the estimation is poor