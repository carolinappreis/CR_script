clear all
% iii=[1 2 3 4 5 8 10 11 12 13 16 17 18];

iii=[2 3 4 5 8 10 11 13 16 17];
for numb=1:length(iii);
    clearvars -except iii  numb nostim nostimout
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
    %     load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(iii(numb)),'_baseline.mat'))
    
    
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
    %
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
    phase=angle(dummy);
    zenv=(abs(hilbert(zscore(tf_3))))';
    for pt=1:3;
        frequency(pt,:)=(smooth((1000/(2*pi))*diff(unwrap(phase(pt,:))),500))';
    end
    
    for ax=1:3
        
        tremor_k= NaN(1,5e4);
        for i=1:5e4
            ix=randi(length(segmentb),1);
            segment=randi([segmentb(ix)+1000 segmente(ix)-5000],1);
            begin3=segment;
            end3=floor(begin3+5*samplerate);
            
            
            if ~isnan(begin3)
                tremor_or2(i,1:(end3-begin3+1))=unwrap(phase(ax,begin3:end3));
                tremor_or22(i,1:(end3-begin3+1))=(phase(ax,begin3)+(0:1:(end3-begin3))*2*pi/(1000./mean(frequency(ax,begin3-1000:begin3))));
                tremor_k(1,i)= (tremor_or2(i,(end3-begin3+1))-tremor_or22(i,(end3-begin3+1)))/(2*pi*0.001*(end3-begin3)); %mean(frequency(end3-1000:end3));%
            else
                tremor_or22(i,1:5001)=NaN;
                tremor_or2(i,1:5001)=NaN;
                tremor_k(1,i)=NaN;
            end
        end
        
        rep=10;
        for i=1:1e6
            dum=tremor_k(randi(5e4,1,rep));
            dum2=dum;
            p(i)=nanmedian(dum2);
        end
        nostim(numb,ax,:)=p;
        clear p
        for i=1:12
            dum=baseline3(randi(5e4,1,rep));
            dum2=dum;
            nostimout(numb,ax,i)=nanmedian(dum2);
        end
        clear dum dum2 baseline3
    end
    
    clearvars -except nostimout iii numb  nostim
end

cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')

