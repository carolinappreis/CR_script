clear all
iii=[1];

for numb=1;
    %     :length(iii);
    clearvars -except iii numb ttall ampall ph_stim LS tt1
    DBS_Fpeak
    
    cond={'_NS_BL.mat';'_NS_PS.mat';'_HFS_PS.mat'};
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
    tremor=(data(3,:));% %score(:,1)';%
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremorx(1:size(ts1.data,3))=ts1.data;
    tremor=(data(6,:));% %score(:,1)';%
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremory(1:size(ts1.data,3))=ts1.data;
    tremor=(data(7,:));% %score(:,1)';%
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremorz(1:size(ts1.data,3))=ts1.data;
    timeor=0:1/samplerate:(size(tremorx,2)-1)/samplerate;
    
    subplot(3,1,1)
    plot(timeor(1,:), tremorx(1,:));
    subplot(3,1,2)
    plot(timeor(1,:), tremory(1,:));
    subplot(3,1,3)
    plot(timeor(1,:), tremorz(1,:));
    
    start{1,1}=[8975 209001 378001];
    start{2,1}=[89281 242600 393001];
    ending{1,1}=[81221 242500 391701];
    ending{2,1}=[161001 301201 454001];
    
    
    
    
    for hh=1:2
           figure(6)
            subplot(2,1,hh)
            plot(timeor,envelope)
        for i=1:3
%             figure(5)
%             subplot(2,1,hh)
%             plot3(timeor(1,start{hh,1}(i):ending{hh,1}(i)),tremorx(1,start{hh,1}(i):ending{hh,1}(i)), tremorz(1,start{hh,1}(i):ending{hh,1}(i)));
%             hold on
         
            hold on
            plot(timeor(1,start{hh,1}(i):ending{hh,1}(i)),envelope(1,start{hh,1}(i):ending{hh,1}(i)),'k')
            plot(timeor(1,start{hh,1}(i):ending{hh,1}(i)),envelope(1,start{hh,1}(i):ending{hh,1}(i)),'r')
        end
    end
    
segmentb=  [8975 89281 242600 209001 378001 393001];
segmente= [81221 161001 242500  301201 391701  454001]; 
    
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






