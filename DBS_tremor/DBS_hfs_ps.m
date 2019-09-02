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
    load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),cond{cc,1}))
    % load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_HFS_PS.mat'));

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

    
    plot(C)
    hold on
    plot(AA,C(AA),'r.')
    plot(BB,C(BB),'b.')
    
    
    if (Fpeak-2)>=1
    [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
    [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end

tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
dummy=hilbert(tremor_or);
envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
phase=angle(dummy);
frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),251))';

amp_e=NaN(length(segmentb),1);
amp_b=NaN(length(segmentb),1);

for i=1:length(segmentb)
    if (~isnan(segmentb(i)))
            amp_b(i,1)=mean(envelope(segmentb(i)-1000:segmentb(i)));
            amp_e(i,1)=mean(envelope(segmente(i)-1000:segmente(i)));
    else
        amp_e(i,1)=NaN;
        amp_b(i,1)=NaN;
    end
end

% figure(1)
% y2=plot(amp_b,amp_e,'k+');
% y3=lsline;
% set(y3,'LineWidth',2,'Color','red')
% box('off')
% c2=corrcoef(amp_b,amp_e)
% legend(y3,[num2str(c2(1,2))],'box','off')
% figure(2)
% bar(amp_b-amp_e)

amp1_p=amp_b(1:2:end);amp1_s=amp_b(2:2:end);
amp2_p=amp_e(1:2:end);amp2_s=amp_e(2:2:end);

figure (3)
subplot(2,1,1)
y2=plot(amp1_p,amp2_p,'k+');
y3=lsline;
set(y3,'LineWidth',2,'Color','red')
box('off')
c2=corrcoef(amp1_p,amp2_p)
legend(y3,[num2str(c2(1,2))],'box','off')

subplot(2,1,2)
y2=plot(amp1_s,amp2_s,'k+');
y3=lsline;
set(y3,'LineWidth',2,'Color','red')
box('off')
c2=corrcoef(amp1_s,amp2_s)
legend(y3,[num2str(c2(1,2))],'box','off')

figure(4)
subplot(2,1,1)
bar(amp1_p-amp2_p)
subplot(2,1,2)
bar(amp1_s-amp2_s)



end