clear all
close all
%  iii=[1 2 3 4 5 8 10 11 12 13 16 17 18];
iii=[2 3 4 5 8 10 11 13 16 17];
main=[1 2 3 ; 1 2 3; 3 2 1; 1 2 3; 3 2 1; 3 2 1; 3 2 1; 3 2 1; 1 2 3; 1 2 3];
for numb= 1:length(iii);
    %          load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\reference_pls_ns.mat')
    clearvars -except iii numb main ns_hu ns_ee ns_ss ns_change ns_befstim
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
    %     load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(iii(numb)),'_baseline.mat'))
    %
    data=SmrData.WvData;
    samplerateold=SmrData.SR;
    ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    tt=0:1/samplerate:(size(ds_data,2)-1)/samplerate;
    tre_3=ds_data([3 5 6],:); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    before_ns
    subj_b=[2 3 4 5 8 10 11 13 16 17];
    ha=find(iii(numb)==(subj_b));
    segmentb=hu{ha,:};
    segmente=hd{ha,:};
    %     [d,e]=butter(2,[0.5/(0.5*samplerate) ],'low'); %15
    %     C=(filtfilt(d,e,tre_3'));
    %     figure()
    %     plot(zscore(tre_3(main(numb),:)))
    %     hold on
    %     plot(zscore(C(:,main(numb))))
    %     for i=1:size(segmentb,2)
    %         xline(segmentb(i),'r')
    %         xline(segmente(i),'k')
    %     end
    %     box('off')
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
    if (Fpeak-2)>=1
        [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    else
        [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    end
    tf_3=filtfilt(b,a,tre_3')*10*9.81/0.5;
    tf_3=tf_3';
    envelope=abs(hilbert(tf_3));
    for ax=1:3
        for i=1:length(segmentb)
            ns_start(ax,i)=mean(envelope(main(numb,ax),handup));
        end
    end
    ns_befstim{numb,1}=ns_start;
    clearvars -except iii numb main ns_befstim
    clear numb in2
end
cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')