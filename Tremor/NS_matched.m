clear all
% iii=[1 2 3 4 5 8 10 11 12 13 16 17 18];

iii=[2 3 4 5 8 10 11 13 16 17];
for numb=1:length(iii);
    clearvars -except iii PC A1 B1 numb nostim nostimout samplerate  vr
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

    
    [d,e]=butter(2,[0.5/(0.5*samplerate) ],'low'); %15
    C=(filtfilt(d,e,tre_3'));
            figure()
            for jj=1:3
                subplot(3,1,jj)
                    plot(zscore(tre_3(jj,:)))
                    hold on
                    plot(zscore(C(:,jj)))
                    for i=1:size(segmentb,2)
                        xline(segmentb(i),'r')
                        xline(segmente(i),'k')
                    end
                    box('off')
            end
    
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
    
    
    for ax=1:3
  
        vr=var(zenv(ax,handup));
        pw=zenv(ax,handup);

       
        for j=1:5e4
            ix=randi(length(segmentb),1);
            segment=randi([segmentb(ix)+1000 segmente(ix)-5000],1);
            begin3=segment;
            end3=floor(begin3+5*samplerate);
            baseline3(1,j)=(mean(envelope(ax,end3-1000:end3))-mean(envelope(ax,begin3-1000:begin3)))./mean(envelope(ax,begin3-1000:begin3)); %#ok<*SAGROW> %
            % baseline4(i,j)=(mean(frequency(end3-1000:end3))); %#ok<*SAGROW>
            % %  time=1:length(envelope(1,:));
            % %             plot(time,envelope(1,:))
            % %             hold on
            % %             xline(segment-1000)
            % %             xline(segment+5000)
            
        end
        
        
        % %         rep=10;
        % %         for i=1:1e6
        % %             dum=baseline3(randi(5e4,1,rep));
        % %             dum2=dum;
        % %             p(i)=nanmedian(dum2);
        % %         end
        % %         nostim(numb,ax,:)=p;
        % %         clear p
        % %         for i=1:12
        % %             dum=baseline3(randi(5e4,1,rep));
        % %             dum2=dum;
        % %             nostimout(numb,ax,i)=nanmedian(dum2);
        % %         end
        % %         clear dum dum2 baseline3
    end
    
    clearvars -except nostimout iii numb PC A1 B1 iii stim nostim vr
end
% cohort=[2 3 4 5 6 7 8 10 11 12];
% vr=vr(cohort,:);
% bar(vr)

cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')

% save('cr_18.mat')
% save 'newnonstim2.mat'

% ANS_group=nostimout; clear nostimout
% AS_group=stim;
% clearvars  -except ANS_group AS_group
% cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
% save 'A_group'
% % NS=nostimout;
% % no_s=nostim;
% % clearvars  -except NS no_s
% % cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
% % % cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
% % % save 'amp_NS.mat'
