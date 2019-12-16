clear all
% SMR_File_To_Mat;
cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\raw_patient_data\P03')
load('P03_PLS.mat')
numb=1;
cond=2;
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

% close all
% plot(time,data(2,:))
% hold on
% plot(time,data(1,:))
% plot(time(1,indexes3),data(1,indexes3),'ko')
% plot(time(1,indexes4),data(1,indexes4),'ro')




start=floor((indexes4./samplerateold)*samplerate)+addon;
ending=floor((indexes3./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);

for nn=1:length(start)
    if numel(start(nn):ending(nn))>60000
        start(nn)=start(nn);
        ending(nn)=ending(nn);
    else
        start(nn)=NaN;
        ending(nn)=NaN;
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

% %% dc_shift
[d,e]=butter(3,(0.5)/(0.5*samplerateold),'low');
dc=(filtfilt(d,e,data(6,:))-2)*10;
% [l1 l2]=findpeaks(dc*-1);
% [h1 h2]=findpeaks(dc);
% [dh1 dh2]=sort(h1,'descend');
%
% for i=1:numel(start)
% h_up(1,i)=h2(dh2(i));
% end
if cond==1;
    start=start(1:3);
    ending=ending(1:3);
    bef_s=floor([56.36 275.8 491.7]*samplerate);
    tap_in(1,:)=floor([112.2 339.7 557.7]*samplerate);
    tap_in(2,:)=floor([143.6 369.4 585.6]*samplerate);
    tap_out(1,:)=floor([116.1 343 560.7]*samplerate);
    tap_out(2,:)=floor([147.4 372 589.1]*samplerate);
    
else
    start=start(4:5);
    ending=ending(4:5);
    bef_s=floor([757.9 1004]*samplerate);
    tap_in(1,:)=NaN;
    tap_in(2,:)=NaN;
    tap_out(1,:)=NaN;
    tap_out(2,:)=NaN;
end

envelope=abs(hilbert(tremor_or));

for rr=1:length(start)
    amp_b(1,rr)=mean(envelope(start(rr)-1000));
    amp_after(1,rr)=mean(envelope(ending(rr)-1000)); 
    trace(rr,:)=envelope(start(rr)-1000:start(rr)+95000-1); %%% clean short spiral and remove 5s cut
    idx_t(rr,:)=[bef_s(rr) ending(rr)-5000]; %%% clean short spiral and remove 5s cut
    change(1,rr)=(mean(envelope(ending(rr)-1000))-mean(envelope(start(rr)-1000)))./mean(envelope(start(rr)-1000));
end

% [p,h]=ttest(amp_b,amp_after)

ts=timeseries(dc,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
dc2(1:size(ts1.data,3))=ts1.data;

close all
f1=figure;
time2=1:length(tremor_or);
plot(time2,tremor_or,'LineWidth',1,'Color',[0.5 0.5 0.5])
hold on
plot(time2,dc2,'LineWidth',1)

    for r=1:size(tap_in,1)
        for i=1:length(bef_s)
            xline(time2(start(i)),'g','LineWidth',2)
            xline(time2(ending(i)),'g','LineWidth',2)
            if ~isnan(tap_in)
            xline(tap_in(r,i),'r--','LineWidth',0.5)
            xline(tap_out(r,i),'k--','LineWidth',0.5)
            end
        end
    end

box('off')
xlim([bef_s(1)-5000 ending(end)+5000])

f2=figure;
subplot(size(trace,1)+1,1,1)
color_b1=[0.5 0.5 0.5];
y2=median(trace);
y1=y2+std(trace);
y3=y2-std(trace);
time3=1:length(trace);
p1=plot(time3, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b1)
patch([time3 fliplr(time3)], [y1 fliplr(y2)],[color_b1],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time3 fliplr(time3)], [y2 fliplr(y3)],[color_b1],'FaceAlpha',[0.2],'EdgeColor','none')
box('off')
xlim([0 length(y2)])

for r=1:size(tap_in,1)
    for i=1:size(trace,1)
        subplot(size(trace,1)+1,1,1+i)
        plot(time2(idx_t(i,1):idx_t(i,2)),envelope(idx_t(i,1):idx_t(i,2)))
        hold on
        xline(time2(start(i)),'g','LineWidth',2)
        if ~isnan(tap_in)
            xline(tap_in(r,i),'r--','LineWidth',0.5)
            xline(tap_out(r,i),'k--','LineWidth',0.5)
        end
        box('off')
        xlim([time2(idx_t(i,1)) time2(idx_t(i,2))])
    end
end
% xlim([0 61000])
% xticks(0:10000:61000)
% xticklabels({'0','10','20','30','40','50','60'})
% xlabel('Time (sec)')
% ylabel('Change in tremor severity (m/s^2)')

f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 24, 12];
set(f1,'color','w');


% 
% start-bef_s
% mean([tap_out(1,:)-tap_in(1,:) tap_out(2,:)-tap_in(2,:)])
% std([tap_out(1,:)-tap_in(1,:) tap_out(2,:)-tap_in(2,:)])




% xlim([0 300000])
% xticks(0:60000:300000)
% xticklabels({'0','1','2','3','4','5'})
% plot(time,data(2,:))
% hold on
% plot(time,data(1,:),'LineWidth',1,'Color','r')
% plot(time,dc,'LineWidth',1)
% ylim([-0.5 0.5])
% set(gca,'FontSize',14)
% title('filtered tremor')
% box('off')
% ylabel('Acceleration (m/s^2)')
% xlabel('Time (min) ')

% cd('/Users/Carolina/OneDrive - Nexus365/arcs_share/peripheral')
% saveas(f1,['PLS_plots',num2str(iiii(numb)),'.png'])