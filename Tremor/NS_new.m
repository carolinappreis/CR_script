clear all
iii=[1 2 3 4 5 8 10 11 12 13 16 17];
in2=1;
for numb=length(iii);
    Fpeak_tremor
    clearvars -except iii PC A1 B1 numb NS NS_i Fpeak samplerate in2
    %      load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
    load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(iii(numb)),'_baseline.mat'))
    
    
    
    if in2==1
        in=3;
    elseif in2==2 % other axis 1
        in=5;
    elseif in2==3 % other axis 2
        in=6;
    end
    data=SmrData.WvData;
    samplerateold=SmrData.SR;
    tremor=(data(in,:));
    addon=92; addon_end=35;
    
    
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremor2(1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    tt=0:1/samplerate:(size(tremor2,2)-1)/samplerate;
    
    
    if (Fpeak-2)>=1
        [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    else
        [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    end
    tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
    % tremor_or=zscore(tremor_or);
    dummy=hilbert(tremor_or);
    envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
    phase=angle(dummy);
    frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
    
    tremor=(data(3,:));% %score(:,1)';%
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremorx(1:size(ts1.data,3))=ts1.data;
    tremor=(data(5,:));% %score(:,1)';%
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremory(1:size(ts1.data,3))=ts1.data;
    tremor=(data(6,:));% %score(:,1)';%
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremorz(1:size(ts1.data,3))=ts1.data;
    tremorxf=filtfilt(b,a,tremorx);
    tremoryf=filtfilt(b,a,tremory);
    tremorzf=filtfilt(b,a,tremorz);
    envelope=[abs(hilbert(tremorxf));abs(hilbert(tremoryf));abs(hilbert(tremorzf))];
    phase=[angle(hilbert(tremorxf));angle(hilbert(tremoryf));angle(hilbert(tremorzf))];
    t_all=[tremorx;tremory;tremorz];
    t_allf=[tremorxf;tremoryf;tremorzf];
    
    [b,a]=butter(2,[0.5/(0.5*samplerate) ],'low'); %15
    C=(filtfilt(b,a,t_all'))';
    %     close all
    %     for i=1:3
    %         subplot(3,1,i)
    %         plot(zscore(t_allf(i,:)))
    %         hold on
    %         plot(zscore(C(i,:)),'LineWidth',1)
    %     end
    %
    
    hu={round(([6.849 7.91 9.731])*(10^5),1);round([1.598 3.55 5.142 6.585 8.159]*(10^5),1);round([1.171/10 1.072 2.272 3.485 4.681]*(10^5),1)...
        ;round([1.116 2.315 4.253 5.885 7.026]*(10^5),1);round([5.257 7.746 9.313 1.069*10 1.183*10]*(10^5),1);...
        round([2.159/10 1.4600 2.9230 4.3520 5.6510]*(10^5),1);round([6110/10^5 1.47 2.754 4.232 5.503]*(10^5),1);round([2.4/10 1.317 3.012]*(10^5),1);...
        round([10211 134450 340484 466600 639582],1);round([5265 139104 277397 401396 551971],1);round(([16972 158835 388948 582772 734164]),1);round([7.175/10 2.61 4.43 ]*(10^5),1)};
    
    hd={round(([7.244 8.958 1.097*10])*(10^5),1);round([2.43 4.486 6.051 7.509 9.037]*(10^5),1);round([5.121/10 1.774 2.905 4.186 5.348]*(10^5),1)...
        ;round([1.698 3.091 5.132 6.481 7.825]*(10^5),1);round([6.386 8.465 9.991 1.131*10 1.24*10]*(10^5),1);...
        round([8.309/10 2.1070 3.5590 4.9500 6.3050]*(10^5),1);round([7.305/10 2.104 3.401 4.898 6.207]*(10^5),1);round([9.98/10 2.385 4.394]*(10^5),1);...
        round([85248 230271 431076 548826 747772],1);round([84050 211850 343728 477800 622380],1);round(([76448 254875 504159 650313 814627]),1);round([1.296 3.362 5.071 ]*(10^5),1)};
    
    
    
    segmentb=hu{numb,:};
    segmente=hd{numb,:};
    unstable3=[];
    
    %         close all
    %          plot(zscore(tremor_or))
    %          hold on
    %         plot(zscore(C(1,:)))
    %         for i=1:size(segmentb,2)
    %        xline(segmentb(i),'r')
    %        xline(segmente(i),'k')
    %         end
    %
    
    %
    %         zc=zscore(C);
    %         close all
    %         plot(tt,zscore(tremor_or))
    %         hold on
    %        plot(tt,zc)
    %         plot(tt(segmentb),zc(segmentb),'ko')
    %         plot(tt(segmente),zc(segmente),'ro')
    %
    for ax=1:3
        for i=1
            for j=1:5e4
                ix=randi(length(segmentb),1);
                segment=randi([segmentb(ix)+1000 segmente(ix)-5000],1);
                begin3=segment;
                end3=floor(begin3+5*samplerate);
                while ~isempty(intersect(unstable3,begin3:end3))
                    segment=randi([segmentb(ix)+1000 segmente(ix)-5000],1);
                    begin3=segment;
                    end3=floor(begin3+5*samplerate);
                end
                baseline3(i,j)=(mean(envelope(ax,end3-1000:end3))-mean(envelope(ax,begin3-1000:begin3)))./mean(ax,envelope(begin3-1000:begin3)); %#ok<*SAGROW> %
                %                 baseline4(i,j)=(mean(frequency(end3-1000:end3))); %#ok<*SAGROW>
                
            end
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
        clear dum dum2
    end
    clearvars -except nostimout iii numb PC A1 B1 iii stim nostim
end
 cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
 save 'amp_NS.mat'

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
