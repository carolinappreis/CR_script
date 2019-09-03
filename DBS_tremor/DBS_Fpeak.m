
    clearvars -except iii numb ttall ampall ph_stim LS tt1
     load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))
%          load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_RS_PS.mat'));
    
    
    
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
    
    Fpeak=frange(find(Pxxrange==max(Pxxrange)));

clearvars -except iii numb  Fpeak 