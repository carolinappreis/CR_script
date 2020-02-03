clear all
close all
%  iii=[1 2 3 4 5 8 10 11 12 13 16 17 18];
iii=[2 3 4 5 8 10 11 13 16 17];
in2=1;


main=[1 1 3 1 3 3 3 3 1 1];
ax=1:3;
gg=[];
for numb= 1:length(iii);
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\reference_pls_ns.mat')
    clearvars -except iii numb NS NS_i in2 yy psd_curves peaks prm ns_ref time_ns gg ax main id_param
    %     load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
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
%         figure()
%         for jj=1:3
%             subplot(3,1,jj)
%             plot(zscore(tf_3(:,jj)))
%             hold on
%             plot(zscore(C(:,jj)))
%             for i=1:size(segmentb,2)
%                 xline(segmentb(i),'r')
%                 xline(segmente(i),'k')
%             end
%             box('off')
%         end
        %
    
    if min(segmentb)<20000
        segmentbp=segmentb(2:end);
        segmentep=segmente(2:end);
    else
        segmentbp=segmentb;
        segmentep=segmente;
    end
    
    ini=10000;
    numb_p=[];
    for i=1:length(segmentbp)
        gg=[gg smooth(zenv(ax(main(numb)),(segmentbp(i)-ini):(segmentbp(i)+ini-1)))];
        
        numb_p=[numb_p  smooth(zenv(ax(main(numb)),(segmentbp(i)-ini):(segmentbp(i)+ini-1)))];
        %             plot(gg)
        %             hold on
    end
    
    y=median(numb_p');
    x=1:length(y);
    initial_params=[];
    [param]=sigm_fit(x,y,initial_params)        % automatic initial_params  "min", "max", "x50" and "slope"
    clear x y
    %     close all
    
    id_param(numb,:)=param;
    
    
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

m50=median(id_param(:,3))


y=smooth(median(gg'));
x=1:length(y)
initial_params=[];
[param]=sigm_fit(x,y,initial_params)        % automatic initial_params
clear x y
