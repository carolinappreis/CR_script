clear all
all=[2 3 4 5 8 10 11 13 16];
% iiii=[2 5 8];

iiii=all([ 1 2 8]);

for numb=1:length(iiii);
    clearvars -except iiii numb
    % load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\PLS_1\p0',num2str(iiii(numb)),'_PLS2.mat'))
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/PLS_1/p0',num2str(iiii(numb)),'_PLS2.mat'))
    
       
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
%     
% plot(data(2,:))
% hold on
% plot(data(1,:))



%%% determine stimulation time points
index=[];
for i=2:size(data,2)-1
    if data(2,i-1)<2.5 && data(2,i)>2.5
        index=[index i];
    end
end
clear i

indexes4=[index(1) index(find(diff(index)./samplerateold > 0.95)+1)];
indexes3=[index(find(diff(index)./samplerateold > 0.95)) index(end)];

dd2=round(data(4,:)*100)./100;
for i=1:length(indexes4)
    xx(i)=round(dd2(indexes4(i))./0.1); %#ok<*SAGROW>
end
clear i

start=floor((indexes4./samplerateold)*samplerate)+addon;
ending=floor((indexes3./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);

for nn=1:length(start)
    if numel(start(nn):ending(nn))>30000
        start(nn)=start(nn);
        ending(nn)=ending(nn);
    else
        start(nn)=NaN
        ending(nn)=NaN
    end
end

start=start(~isnan(start));
ending=ending(~isnan(ending));

%%% when patient's hand is up
handup=[];
for i=1:length(start)
    handup=[handup start(i):ending(i)]; %#ok<*AGROW>
end
clear i
handup=sort(handup,'ascend');



%%% tremor characteristics
[Pxx,F]=pwelch(tremor2(handup),samplerate,[],samplerate,samplerate);

frange=F(3:10);
Pxxrange=Pxx(3:10);

Fpeak=frange(find(Pxxrange==max(Pxxrange))); %#ok<*FNDSB>

if (Fpeak-2)>=1
    [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
    [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end

tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;

%%% dc_shift
% [d,e]=butter(3,(0.5)/(0.5*samplerateold),'low');
% dc=(filtfilt(d,e,data(6,:))-2)*10;
% plot(data(1,:))
% hold on
% plot(data(2,:))
% plot(dc,'LineWidth',3)
envelope=abs(hilbert(tremor_or));

for rr=1:length(start)
change(1,rr)=(mean(envelope(ending(rr)-1000))-mean(envelope(start(rr)-1000)))./mean(envelope(start(rr)-1000));
amp_b(1,rr)=mean(envelope(start(rr)-1000));
amp_after(1,rr)=mean(envelope(ending(rr)-1000));
trace_1=envelope(start(rr)-1000:ending(rr));
trace(rr,:)=trace_1(1:59115);
clear trace_1
end

[p,h]=ttest(amp_b,amp_after)

% plot(data(2,:))
% hold on
% plot(data(1,:))

f1=figure;
subplot(1,2,1)
time=1:length(tremor_or);
plot(time,tremor_or,'LineWidth',1,'Color',[0.5 0.5 0.5])
xlim([0 300000])
xticks(0:60000:300000)
xticklabels({'0','1','2','3','4','5'})
set(gca,'FontSize',14)
title('filtered tremor')
box('off')
ylabel('Acceleration (m/s^2)')
xlabel('Time (min) ')
subplot(1,2,2)
color_b1=[0.5 0.5 0.5];
y2=median(trace); 
y1=y2+std(trace); 
y3=y2-std(trace);
time=1:length(trace);
p1=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b1)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b1],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b1],'FaceAlpha',[0.2],'EdgeColor','none')
set(gca,'FontSize',14)
box('off')
title('Avg envelope')
xlim([0 61000])
xticks(0:10000:61000)
xticklabels({'0','10','20','30','40','50','60'})
xlabel('Time (sec)')
ylabel('Change in tremor severity (m/s^2)')
f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 24, 12];
set(f1,'color','w');
cd('/Users/Carolina/OneDrive - Nexus365/arcs_share/peripheral')
% savefig(['PLS_plots_1_',num2str(iiii(numb))])
saveas(f1,['PLS_plots_1_',num2str(iiii(numb)),'.png'])

end