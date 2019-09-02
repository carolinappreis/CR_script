clear all
iii=[1];

for numb=1;
    %     :length(iii);
    clearvars -except iii numb ttall ampall ph_stim LS tt1
    DBS_Fpeak

    cond={'_NS_BL.mat';'_NS_PS.mat';'_HFS_PS.mat'};
    PC=[29;75;75];
    A1={[];[1 3 7 16 21 23];[]};
    B1={[];[2 6 12 20 22 26];[]};
    
    cc=2;
%     load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),cond{cc,1}))
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_HFS_PS.mat'));
    
    
    
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
    % tremor_or=zscore(tremor_or);
    dummy=hilbert(tremor_or);
    envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
    phase=angle(dummy);
    frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';

    [b,a]=butter(2,[0.5/(0.5*samplerate) ],'low'); %15
    C=(filtfilt(b,a,envelope));
    
%     plot(tremor_or)
%     hold on
%     plot(envelope,'LineWidth',1.5)
%     plot(C,'LineWidth',1.5)
    
     A=prctile(C,PC(cc,numb));
    ind_s=[];
    for i=11:length(C)-11
        if C(i-1)<A && C(i+1)>A
            ind_s=[ind_s i]; %#ok<*AGROW>
        end
    end
    for i=1:(length(ind_s)-1)
        if ind_s(i+1)-ind_s(i)==1
            ind_s(i+1)=NaN;
        end
    end
    ind_s2=ind_s(~isnan(ind_s));
    ind_e=[];
    for i=11:length(C)-11
        if C(i-1)>A && C(i+1)<A
            ind_e=[ind_e i];
        end
    end
    for i=1:(length(ind_e)-1)
        if ind_e(i+1)-ind_e(i)==1
            ind_e(i+1)=NaN;
        end
    end
    ind_e2=ind_e(~isnan(ind_e));
    %
    % plot(C)
    % hold on
    % plot(ind_s2,C(ind_s2),'r.')
    % plot(ind_e2,C(ind_e2),'b.')
    
    if isempty (A1{cc,numb})
        AA=ind_s2;
        BB=ind_e2;
    else
        AA=ind_s2(A1{cc,numb});
        BB=ind_e2(B1{cc,numb});
    end

    segmentb=AA;
    segmente=BB;
    unstable3=[];
%     
%     figure()
%     plot(C)
%     hold on
%     plot(AA,C(AA),'r.')
%     plot(BB,C(BB),'b.')
%     
   
tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
    [b,a]=butter(2,[3/(0.5*samplerate) ],'low'); %15
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

    
si=  [8975 79831 151001 225401 301001 374501];
ei= [65001 145001 220501 292001 350001 435501];

figure()
subplot(2,1,1)
plot(timeor,C)
hold on
plot(timeor(1,segmentb), C(1,segmentb),'r.');
plot(timeor(1,segmente), C(1,segmente),'ko');
subplot(2,1,2)
plot(timeor,C)
hold on
plot(timeor(1,si),C(1,si),'r.');
plot(timeor(1,ei),C(1,ei),'ko');

f=2;
epoch=si(f):ei(f);

subplot(3,1,1)
plot(timeor(1,epoch), tremorx(1,epoch));
subplot(3,1,2)
plot(timeor(1,epoch), tremory(1,epoch));
subplot(3,1,3)
plot(timeor(1,epoch), tremorz(1,epoch));

figure()
subplot(3,1,1)
plot(tremorx(1,epoch), tremory(1,epoch));
subplot(3,1,2)
plot(tremorx(1,epoch), tremorz(1,epoch));
subplot(3,1,3)
plot(tremory(1,epoch), tremorz(1,epoch));

plot3(timeor(1,:),tremorx(1,:), tremory(1,:));
hold on
plot3(timeor(1,si),tremorx(1,si), tremory(1,si),'rd');

plot3(timeor(1,epoch),tremorx(1,epoch), tremory(1,epoch));

plot3(tremorx(1,epoch),tremory(1,epoch), tremorz(1,epoch));

    

    
figure(1)    
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


    
end