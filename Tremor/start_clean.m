
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

start_run=start(1:12:120);
end_run=ending(12:12:120);
e_nonstim=[start_run length(tremor2)];
s_nonstim=[1 end_run];

ns=[];
for i=1:length(start_run)
    ns=[ns s_nonstim(i):e_nonstim(i)];
end

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

dummy=hilbert(tremor_or);
envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));


%%% removal of time segments [ < mean - std for more than 10 seconds ]
unstable=find(envelope<(mean(envelope(handup))-std(envelope(handup))));

if ~isempty(unstable)
    change_e=[unstable(find(diff(unstable)~=1)) unstable(end)];
    change_b=[unstable(1) unstable(find(diff(unstable)~=1)+1)];
    
    change_ind=find((change_e-change_b)>10000);
    
    unstable2=[];
    if ~isempty(change_ind)
        change_e2=change_e(change_ind);
        change_b2=change_b(change_ind);
        for i=1:length(change_ind)
            unstable2=[unstable2 change_b2(i):change_e2(i)];
        end
    end
    
    unstable3=intersect(handup,unstable2);
    
    clear unstable2;
    unstable2=unstable3;
    
    clear i
    
    start2=start;
    ending2=ending;
    
    if ~isempty(unstable2)
        for i=1:length(unstable2)
            start2(max(find(start<unstable2(i))))=NaN; %#ok<*MXFND>
            ending2(max(find(start<unstable2(i))))=NaN;
        end
    end
    clear i
    
    xx2=xx;
    clear start ending xx
    start=start2(find(~isnan(start2)));
    ending=ending2(find(~isnan(ending2)));
    xx=xx2((find(~isnan(start2))));
    clear start2 ending2 xx2
end