clear all
close all
%  iii=[1 2 3 4 5 8 10 11 12 13 16 17 18];
iii=[2 3 4 5 8 10 11 13 16 17];
in2=1;
gg=[];

for numb= 1:length(iii);
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\reference_pls_ns.mat')
    clearvars -except iii numb NS NS_i in2 yy psd_curves peaks prm ns_ref time_ns gg
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
    %     load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(iii(numb)),'_baseline.mat'))
    
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
    
    before_ns
    
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
    %     figure()
    %     for jj=1:3
    %         subplot(3,1,jj)
    %         plot(zscore(tf_3(:,jj)))
    %         hold on
    %         plot(zscore(C(:,jj)))
    %         for i=1:size(segmentb,2)
    %             xline(segmentb(i),'r')
    %             xline(segmente(i),'k')
    %         end
    %         box('off')
    %     end
    %     %
    
    if min(segmentb)<20000
        segmentbp=segmentb(2:end);
        segmentep=segmente(2:end);
    else
        segmentbp=segmentb;
        segmentep=segmente;
    end
    
    ini=10000;
    %     ini=20000;
    for ax=1:3
        %         gg=[];
        for i=1:length(segmentbp)
            gg=[gg ;(zenv(ax,(segmentbp(i)-ini):(segmentbp(i)+ini)))];
            %             amp_start(ax,1,i)=mean(envelope(ax,segmentbp(i)-1000:segmentbp(i)));
            %             amp_end(ax,1,i)=mean(envelope(ax,segmentep(i)-1000:segmentep(i)));
            %             amp_start(ax,1,i)=mean(zenv(ax,segmentbp(i)-1000:segmentbp(i)));
            %             amp_end(ax,1,i)=mean(zenv(ax,(segmentbp(i)+ns_ref(numb)-5000):segmentbp(i)+ns_ref(numb)));
            plot(gg)
            hold on
        end
        
        %                 f5=figure(5)
        %                 subplot(3,1,ax)
        %                 plot(median(gg),'Color',[0.5 0.5 0.5])
        %                 hold on
        %                 xline(ini,'g--','LineWidth',2)
        %                 title(['NS pt',num2str(iii(numb))])
        %
        %                         yy(numb,:)=smooth(median(gg,1));
        %                         y=yy(numb,1:40000);
        %                         x=tt(1:length(y));
        %                         initial_params=[];
        %                         [param]=sigm_fit(x,y,initial_params)        % automatic initial_params
        %                         clear x y
        
        
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
    main=[1 1 3 1 3 3 3 3 1 1];
    time_ns(numb,:)=timings(main(numb),:);
    p_all(numb,:)=pr(main(numb),:);
    clear timings
    %     f1=figure;
    %     plot((timings'),'.','MarkerSize',20)
    %     hold on
    %     plot((timings'),'LineWidth',0.5)
    %     xlim([0 4])
    %     xticks([1 2 3])
    %     xticklabels({'bef handUP','bef handDW'})
    %     xtickangle(45)
    %     ylabel('Tremor severity ')
    %     legend('boxoff')
    %     f1.Units = 'centimeters';
    %     f1.OuterPosition= [10, 10, 10, 10];
    %     set(gca,'FontSize',14)
    %     set(f1,'color','w');
    %     box('off')
end




y=smooth(median(gg,1));
x=1:length(yy)
initial_params=[];
[param]=sigm_fit(x,y,initial_params)        % automatic initial_params
clear x y
