clear all
iii=[1];

for numb=1;
    %     :length(iii);
    clearvars -except nostimout iii numb cc cond nostim
    
    DBS_Fpeak
    
    cond={'_NS_BL.mat';'_NS_PS.mat';'_HFS_PS.mat'};
    cc=3;
    if cc==2
        segmentb=[8646 88131 174501 234200 315901 393800];
        segmente=[82701 147501 226501 297200 377500 454600];
        start{1,1}=[8646 174501 315901];
        start{2,1}=[88131 234200 393800];
        ending{1,1}=[82701 226501 377500];
        ending{2,1}=[147501 297200 454600];
        
    elseif cc==3
        segmentb=[1621 112001 205601 287501 366701 450401];
        segmente=[76121 179301 263501 350901 422701 517501];
        start{1,1}=segmentb(1:2:end);
        start{2,1}=segmentb(2:2:end);
        ending{1,1}=segmente(1:2:end);
        ending{2,1}=segmente(2:2:end);
    end
    
    
%     load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),cond{cc,1}))
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),cond{cc,1}));
    
    in2=1; % analysing the "main tremor axis"
    
    if in2==1
        in=3;
    elseif in2==2 % other axis 1
        in=6;
    elseif in2==3 % other axis 2
        in=7;
    end
    data=SmrData.WvData;
    samplerateold=SmrData.SR;
    tremor=(data(in,:));
    addon=92; addon_end=35;
    
    time=0:1/samplerateold:(size(data,2)-1)/samplerateold;
    
    %%% downsample
    
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremor2(1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    
    if (Fpeak-2)>=1
        [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    else
        [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    end
    tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
    [b,a]=butter(2,[0.8/(0.5*samplerate) ],'low'); %15
    % tremor_or=zscore(tremor_or);
    dummy=hilbert(tremor_or);
    envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
    phase=angle(dummy);
    frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
    tremor=(data(3,:));% %score(:,1)';%
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremorx(1:size(ts1.data,3))=ts1.data;
    filt_x=filtfilt(b,a,tremorx);
    tremor=(data(6,:));% %score(:,1)';%
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremory(1:size(ts1.data,3))=ts1.data;
    filt_y=filtfilt(b,a,tremory);
    tremor=(data(7,:));% %score(:,1)';%
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremorz(1:size(ts1.data,3))=ts1.data;
    filt_z=filtfilt(b,a,tremorz);
    timeor=0:1/samplerate:(size(tremorx,2)-1)/samplerate;
    
    
    %
    %     subplot(3,1,1)
    %     plot(timeor(1,:), filt_x(1,:));
    %     subplot(3,1,2)
    %     plot(timeor(1,:), filt_y(1,:));
    %     subplot(3,1,3)
    %     plot(timeor(1,:), filt_z(1,:));
    
    
    
    
    
    %     for hh=1:2
    %         figure(6)
    %         subplot(2,1,hh)
    %         plot(timeor,envelope)
    %         hold on
    %         for i=1:3
    %             figure(5)
    %             subplot(2,1,hh)
    %             plot3(timeor(1,start{hh,1}(i):ending{hh,1}(i)),tremorx(1,start{hh,1}(i):ending{hh,1}(i)), tremorz(1,start{hh,1}(i):ending{hh,1}(i)));
    %             hold on
    %             figure(6)
    %             plot(timeor(1,start{hh,1}(i):ending{hh,1}(i)),envelope(1,start{hh,1}(i):ending{hh,1}(i)),'k')
    %             plot(timeor(1,start{hh,1}(i):ending{hh,1}(i)),envelope(1,start{hh,1}(i):ending{hh,1}(i)),'r')
    %         end
    %     end
    %
    unstable3=[];
    rep=10;
    for ss=1:2
        for i=1
            for i=1:5e4
                ix=randi(length(start{ss,1}),1);
                segment=randi([start{ss,1}(ix)+1000 ending{ss,1}(ix)-5000],1);
                begin3=segment;
                end3=floor(begin3+5*samplerate);
                if ~isnan(begin3)
                    tremor_or2(i,1:(end3-begin3+1))=unwrap(phase(begin3:end3));
                    tremor_or22(i,1:(end3-begin3+1))=(phase(begin3)+(0:1:(end3-begin3))*2*pi/(1000./mean(frequency(begin3-1000:begin3))));
                    tremor_k(i,1)= (tremor_or2(i,(end3-begin3+1))-tremor_or22(i,(end3-begin3+1)))/(2*pi*0.001*(end3-begin3)); %mean(frequency(end3-1000:end3));%
                else
                    tremor_or22(i,1:5001)=NaN;
                    tremor_or2(i,1:5001)=NaN;
                    tremor_k(i,1)=NaN;
                end
            end
        end
        
    end
    for rr=1:1e6
        dum=tremor_k(randi(25e3,1,rep));
        dum2=dum;
        NS_i1(rr)=nanmedian(dum2);
    end
    nostim(ss,numb,:)=NS_i1;
    
    for i=1:12
        dum=tremor_k(randi(5e4,1,rep));
        dum2=dum;
        nostimout(ss,numb,i)=nanmedian(dum2);
    end
    
    clearvars -except nostimout iii numb cc cond  nostim
end

fNS=nostimout;
fidv_NS=nostim;
clearvars  -except fNS fidv_NS
% cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data')
cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data')
save 'fHFS_PS_result.mat'

