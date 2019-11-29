clear all
iii=[1 2 3 4 5 8 10 11 12 13 16 17];
in2=1;
for numb=1:length(iii);
    clearvars -except iii PC A1 B1 numb nostim nostimout samplerate in2
    %      load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
    load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(iii(numb)),'_baseline.mat'))
    
    if in2==1
        in=3;
    elseif in2==2 % other axis 1
        in=5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    elseif in2==3 % other axis 2
        in=6;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
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
    tre_3=ds_data([3 5 6],:); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    
    hu={round(([7.91 9.731])*(10^5),1);round([1.598 3.55 5.142 6.585 8.159]*(10^5),1);round([ 1.072 2.272 3.485 4.681]*(10^5),1)...
        ;round([1.116 2.315 4.253 5.885 7.026]*(10^5),1);round([5.257 7.746 9.313 1.069*10 1.183*10]*(10^5),1);...
        round([2.159/10 1.4600 2.9230 4.3520 5.6510]*(10^5),1);round([1.47 2.754 4.232 5.503]*(10^5),1);round([2.4/10 1.317 3.012]*(10^5),1);...
        round([ 134450 340484 466600 639582],1);round([139104 277397 401396 551971],1);round(([16972 158835 388948 582772 734164]),1);round([7.175/10 2.61 4.43 ]*(10^5),1)};
    
    hd={round(([8.958 1.097*10])*(10^5),1);round([2.43 4.486 6.051 7.509 9.037]*(10^5),1);round([ 1.774 2.905 4.186 5.348]*(10^5),1)...
        ;round([1.698 3.091 5.132 6.481 7.825]*(10^5),1);round([6.386 8.465 9.991 1.131*10 1.24*10]*(10^5),1);...
        round([8.309/10 2.1070 3.5590 4.9500 6.3050]*(10^5),1);round([ 2.104 3.401 4.898 6.207]*(10^5),1);round([9.98/10 2.385 4.394]*(10^5),1);...
        round([ 230271 431076 548826 747772],1);round([211850 343728 477800 622380],1);round(([76448 254875 504159 650313 814627]),1);round([1.296 3.362 5.071 ]*(10^5),1)};
    
    segmentb=hu{numb,:};
    segmente=hd{numb,:};
    
        [d,e]=butter(2,[0.5/(0.5*samplerate) ],'low'); %15
        C=(filtfilt(d,e,tre_3'));
        figure()
        for jj=1:3
            subplot(3,1,jj)
                plot(zscore(tre_3(jj,:)))
                hold on
                plot(zscore(C(jj,:)))
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
            for j=1:5e4
                ix=randi(length(segmentb),1);
                segment=randi([segmentb(ix)+1000 segmente(ix)-5000],1);
                begin3=segment;
                end3=floor(begin3+5*samplerate);
                baseline3(1,j)=(mean(envelope(ax,end3-1000:end3))-mean(envelope(ax,begin3-1000:begin3)))./mean(envelope(ax,begin3-1000:begin3)); %#ok<*SAGROW> %
                % baseline4(i,j)=(mean(frequency(end3-1000:end3))); %#ok<*SAGROW>
   
            end
            
        rep=10;
        for i=1:1e6
            dum=baseline3(randi(5e4,1,rep));
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
    clearvars -except nostimout iii numb PC A1 B1 iii stim nostim in2
end
cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
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
