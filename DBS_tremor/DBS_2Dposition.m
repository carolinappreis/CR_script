clear all
iii=[1];

for numb=1;
    %     :length(iii);
    clearvars -except nostimout iii numb cc cond 
    DBS_Fpeak
    
    cond={'_NS_BL.mat';'_NS_PS.mat';'_HFS_PS.mat'};
    cc=3;
    if cc==2
        segmentb=  [8646 88131 174501 234200 315901 393800];
        segmente= [82701 147501 226501 297200 377500 454600];
        start{1,1}=[8646 174501 315901];
        start{2,1}=[88131 234200 393800];
        ending{1,1}=[82701 226501 377500];
        ending{2,1}=[147501 297200 454600];
        
    elseif cc==3
        segmentb= [1621 112001 205601 287501 366701 450401];
        segmente= [76121 179301 263501 350901 422701 517501];
        start{1,1}=segmentb(1:2:end);
        start{2,1}=segmentb(2:2:end);
        ending{1,1}=segmente(1:2:end);
        ending{2,1}=segmente(2:2:end);
    end
    
    
    load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),cond{cc,1}))
%     load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),cond{cc,1}));
    
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
    % figure()
    % subplot(2,1,1)
    % plot(timeor,C)
    % hold on
    % plot(timeor(1,segmentb), C(1,segmentb),'r.');
    % plot(timeor(1,segmente), C(1,segmente),'ko');
    % subplot(2,1,2)
    % plot(timeor,C)
    % hold on
    % plot(timeor(1,si),C(1,si),'r.');
    % plot(timeor(1,ei),C(1,ei),'ko');
    
    f=2;
    epoch=segmentb(1):segmente(2);
    
    subplot(3,1,1)
    plot(timeor(1,epoch), filt_x(1,epoch));
    subplot(3,1,2)
    plot(timeor(1,epoch), filt_y(1,epoch));
    subplot(3,1,3)
    plot(timeor(1,epoch), filt_z(1,epoch));
    
    figure()
    subplot(3,1,1)
    plot(filt_x(1,epoch), filt_y(1,epoch));
    subplot(3,1,2)
    plot(filt_x(1,epoch), filt_z(1,epoch));
    subplot(3,1,3)
    plot(filt_y(1,epoch), filt_z(1,epoch));
    
    plot3(timeor(1,:),filt_x(1,:), filt_z(1,:));
    hold on
    plot3(timeor(1,segmentb),filt_x(1,segmentb), filt_y(1,segmentb),'rd');
    
    plot3(timeor(1,epoch),tremorx(1,epoch), filt_y(1,epoch));
    
    figure()
    for f=1:6
    epoch=segmentb(f):segmente(f);
    plot3(filt_x(1,epoch),filt_y(1,epoch), filt_z(1,epoch),'color',rand(1,3));
    hold on
    n(f,:)=segmente(f)-segmentb(f);
    end
    xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')
    
    
    figure(1)
    xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')
    for f=1:6
        epoch=segmentb(f):segmente(f);
        h = animatedline;
        numpoints = length(epoch);
        % x = tremorx(1,epoch);
        % y = tremory(1,epoch);
        % z= tremorz(1,epoch);
        x = filt_x(1,epoch);
        y = filt_y(1,epoch);
        z= filt_z(1,epoch);
        a = tic; % start timer
        for k = 1:numpoints
            addpoints(h,x(k),y(k),z(k));
            b = toc(a); % check timer
            if b > (1/1000)
                drawnow % update screen every 1/30 seconds
                a = tic; % reset timer after updating
                k;
            end
        end
         h = animatedline(x,y,z,'color',rand(1,3));
        hold on
    end
    
end