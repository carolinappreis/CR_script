clear all
all=[2 3 4 5 8 10 11 13 16];
% iiii=[2 5 8];
f=1;
for rep=1:2
    if rep==1
        iiii=all([1 3 4 ]);
    else
        iiii=all(1);
    end
    for numb=1:size(iiii,2);
        % load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\PLS\p0',num2str(iiii(numb)),'_PLS.mat'))
        if rep==1
            load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/PLS/p0',num2str(iiii(numb)),'_PLS.mat'))
        else
            load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/PLS_1/p0',num2str(iiii(numb)),'_PLS2.mat'))
        end
        
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
        
        %%% dc_shift
        % [d,e]=butter(3,(0.5)/(0.5*samplerateold),'low');
        % dc=(filtfilt(d,e,data(6,:))-2)*10;
        % plot(data(1,:))
        % hold on
        % plot(data(2,:))
        % plot(dc,'LineWidth',3)
        envelope=zscore(abs(hilbert(tremor_or)));
        
        bin=floor(50000/3);
        
        for rr=1:length(start)
            
            change(rr,:)=[(mean(envelope(start(rr)-1000:start(rr)))) ...
                mean(envelope(start(rr)+bin-500:start(rr)+bin+500)) ...
                mean(envelope(start(rr)+2*bin-500:start(rr)+2*bin+500)) ...
                mean(envelope(ending(rr)-1000:ending(rr)))];
        end
        
        amp(f,:)=mean(change);
        
        clearvars -except iiii all numb amp f rep
        f=f+1;
    end
    
end

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','olive')

f1=figure;
plot(1:4,amp,'.','MarkerSize',7)
hold on
plot(median(amp,1),'lineWidth',3,'Color',olive)
xlim([0 5])
xticks([1 2 3 4])
xticklabels({'before','mid1',' mid2','end'})
ylim([-1 2])
yticks(-1:1:2)
ylabel('Tremor severity (zscore)')
f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 10, 10];
set(gca,'FontSize',14)
set(f1,'color','w');
box('off')
