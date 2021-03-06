time=0:1/samplerateold:(size(data,2)-1)/samplerateold;

%% downsample

ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremor2(1:size(ts1.data,3))=ts1.data;
samplerate=1000;

%% determine stimulation time points
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

%% when patient's hand is up
handup=[];
for i=1:length(start)
    handup=[handup start(i):ending(i)]; %#ok<*AGROW>
end
clear i
handup=sort(handup,'ascend');

%% tremor characteristics
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


%% removal of time segments [ < mean - std for more than 10 seconds ]
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

%% re - estimate tremor characteristics
clear handup Pxx F frange Pxxrange Fpeak tremor_or dummy envelope phase frequency

handup=[];
for i=1:length(start)
        handup=[handup start(i):ending(i)]; %#ok<*AGROW>
end
handup=sort(handup,'ascend');

[Pxx,F]=pwelch(tremor2(handup),samplerate,[],samplerate,samplerate);

frange=F(3:10);
Pxxrange=Pxx(3:10);

Fpeak=frange(find(Pxxrange==max(Pxxrange)));

if (Fpeak-2)>=1
    [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
    [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end
tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;

dummy=hilbert(tremor_or);
envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
phase=angle(dummy);
frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';

tremor=(data(3,:));% %score(:,1)';%
ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremorx(1:size(ts1.data,3))=ts1.data;
tremor=(data(5,:));% %score(:,1)';%
ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremory(1:size(ts1.data,3))=ts1.data;
tremor=(data(6,:));% %score(:,1)';%
ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremorz(1:size(ts1.data,3))=ts1.data;
if (Fpeak-2)>=1
    [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
    [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end
tremorxf=filtfilt(b,a,tremorx);
tremoryf=filtfilt(b,a,tremory);
tremorzf=filtfilt(b,a,tremorz);

for j=1:length(start)
    x=[tremorxf(start(j):ending(j));tremoryf(start(j):ending(j));tremorzf(start(j):ending(j))];
    [pc,score,latent,tsquare] = pca(x');
    xxx(j,1:3)=pc(1:3,1);
    ma(j)=(find(abs(xxx(j,1:3))==max(abs(xxx(j,1:3)))));
end

tremor_or2=NaN(length(start),1);

for i=1:length(start)
    if (~isnan(start(i))) % && ma==1) % dont forget to delete the ) when removing comment 
         tremor_or2(i,1)=(mean(envelope(ending(i)-1000:ending(i)))-mean(envelope(start(i)-1000:start(i))))/mean(envelope(start(i)-1000:start(i)));
    else
        tremor_or2(i,1)=NaN;
    end
end


clear tt
k=1;
tt=NaN(20,12);

for i=1:12
    tt(1:sum(xx==i),i)=tremor_or2(find(xx==i)); % please check tt does not contain
    % zeros - this can happen if tremor amplitude was too small and
    % stimulation started and stopped within the 5 second window - if so
    % execute
    
    % tt(tt==0)=NaN;
end



