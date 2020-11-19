clear all
% iii=[1 2 3 4 5 8 10 11 12 13 16 17 18];

iii=[2 3 4 5 8 10 11 13 16 17];
main=[1 1 3 1 3 3 3 3 1 1];
for numb=1:length(iii);
    clearvars -except iii numb ns_msplit main
    %      load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
    load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(iii(numb)),'_baseline.mat'))
    
    data=SmrData.WvData;
    samplerateold=SmrData.SR;
    
    ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    tt=0:1/samplerate:(size(ds_data,2)-1)/samplerate;
    tre_3=ds_data([3 5 6],:); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    
    before_ns
    
    segmentb=hu{numb,:};
    segmente=hd{numb,:};
    
    
    %     [d,e]=butter(2,[0.5/(0.5*samplerate) ],'low'); %15
    %     C=(filtfilt(d,e,tre_3'));
    %     figure()
    %     for jj=1:3
    %         subplot(3,1,jj)
    %         plot(zscore(tre_3(jj,:)))
    %         hold on
    %         plot(zscore(C(:,jj)))
    %         for i=1:size(segmentb,2)
    %             xline(segmentb(i),'r')
    %             xline(segmente(i),'k')
    %         end
    %         box('off')
    %     end
    
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
    cr(numb,:)=Fpeak;
    
    %     figure()
    %     plot(F(3:50),ps_curves(:,3:50)','LineWidth',2)
    %     legend({'z','y','x'})
    %     legend('boxoff')
    %     box('off')
    
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
    
    ax=main(numb);
    time=1:length(envelope(1,:));
    
    for j=1:5e4
        ix=randi(length(segmentb),1);
        segment=randi([segmentb(ix)+1000 segmente(ix)-5000],1);
        begin3=segment;
        end3=floor(begin3+5*samplerate);
        baseline3(1,j)=(mean(envelope(ax,end3-1000:end3))-mean(envelope(ax,begin3-1000:begin3)))./mean(envelope(ax,begin3-1000:begin3)); %#ok<*SAGROW> %
        bl_amp(1,j)=mean(envelope(ax,begin3-1000:end3));
    end
    
    
    
    amp_1=NaN(2,round(size(bl_amp,1)./2));
    ch_a1=NaN(2,round(size(baseline3,1)./2));
    m=1;
    n=1;
    for i=1:length(bl_amp)
        if bl_amp(i)<=nanmedian(bl_amp)
            amp_1(1,n)= bl_amp(i);
            ch_a1(1,n)= baseline3(i);
            n=n+1;
            
        else
            amp_1(2,m)= bl_amp(i);
            ch_a1(2,m)= baseline3(i);
            m=m+1;
            
        end
    end
    
    
    rep=10;
    for k=1:2
        for i=1:1e6/2
            dum=ch_a1(k,randi(5e4./2,1,rep));
            ns_msplit(numb,k,i)=nanmedian(dum);
            clear dum 
        end
    end
    
    
    
    %                 for i=1:12
    %                     dum=baseline3(randi(5e4,1,rep));
    %                     dum2=dum;
    %                     nostimout(numb,ax,i)=nanmedian(dum2);
    %                 end
    
end

cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
% cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
clearvars -except iii ns_msplit


%%%%-------------------




