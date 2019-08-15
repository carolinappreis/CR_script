%#ok<*NASGU> 
%#ok<*UDIM>

time=0:1/samplerateold:(size(data,2)-1)/samplerateold;
numberoftrials=15;

%% downsample

ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremor2(1:size(ts1.data,3))=ts1.data;
samplerate=1000;

%% determine stimulation time points ( start and end of 5 second long phase locking between stim and tremor )
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

start=floor((indexes4./samplerateold)*samplerate)+addon; % addon=delay from TTL to stimulation delivery 
ending=floor((indexes3./samplerateold)*samplerate)+addon+addon_end; % addon_end accounts for burst of stimulation delivered
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

variable=hilbert(tremor_or);
envelope=sqrt((real(variable).^2)+(imag(variable).^2));


%% removal of time segments [ < mean - std for more than 10 seconds ]
% I remove segments where tremor dissipates for more than 10 seconds
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
    
    clear start ending
    start=start2;
    ending=ending2;
    clear start2 ending2
end

% patient specific segment removal (due to voluntary movement etc) 
% if remove==1
%     start(find(start<300*1000 ))=NaN;
% elseif remove==2;
%     start(find(start<100*1000 ))=NaN;
%     start(find(start>648*1000 & start<656*1000))=NaN;
% end
start2=start(~isnan(start));
ending2=ending(~isnan(start));
xx2=xx(~isnan(start));
clear start ending xx
start=start2;
ending=ending2;
xx=xx2;

%% re - estimate tremor characteristics
clear handup Pxx F frange Pxxrange Fpeak tremor_or variable envelope phase frequency

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

variable=hilbert(tremor_or);
envelope=sqrt((real(variable).^2)+(imag(variable).^2));
phase=angle(variable);
frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(variable))),251))'; % frequency estimated from hilbert

%% alternative frequency estimation using windowed FFT ( resolution of 0.2 Hz )
%%% there is a problem with both ways of estimating frequency either due to 
%%% using window or smoothing etc ... 

% for i=1:length(tremor_or)-5000
% seg=tremor_or(i:i+5000);
% f_seg=abs(real(fft(seg)));
% F=0:(1000/5000):(1000-0.2);
% f_alt(i)=F(min(find(f_seg==max(f_seg))));
% end
     
%% tremor_or2 - tremor_or22 is the difference in slope (i.e. change in frequency)

% tremor_or2=NaN(20,5001);
% tremor_or22=NaN(20,5001);
% 
% [a,b]=hist(frequency,0:0.05:10);
% 
% 
% for i=1:length(start)
%     if ~isnan(start(i)) 
%         tremor_or2(i,1:(ending(i)-start(i)+1))=unwrap(phase(start(i):ending(i)));
%         tremor_or22(i,1:(ending(i)-start(i)+1))=(phase(start(i))+(0:1:(ending(i)-start(i)))*2*pi/(1000./mean(frequency(start(i)-1000:start(i)))));
%         tremor_k(i,1)= (tremor_or2(i,(ending(i)-start(i)+1))-tremor_or22(i,(ending(i)-start(i)+1)))/(2*pi*0.001*(ending(i)-start(i))); %mean(frequency(ending(i)-1000:ending(i)));%
% %           tremor_k(i,1)= mean(frequency(ending(i)-1000:ending(i)));%
%     else
%         tremor_or22(i,1:5001)=NaN;
%         tremor_or2(i,1:5001)=NaN;
%         tremor_k(i,1)=NaN;
%     end
% end

%% tremor_or2 - tremor_or22 is the difference in frequency (estimated from unwrapped hilbert)


for i=1:length(start)
    if ~isnan(start(i)) 
        tremor_or2(i,1:(ending(i)-start(i)+1))=frequency(start(i):ending(i));
        tremor_or22(i,1:(ending(i)-start(i)+1))=mean(frequency(start(i)-1000:start(i)));
        tremor_k(i,1)=tremor_or2(i,(ending(i)-start(i)+1))-tremor_or22(i,(ending(i)-start(i)+1));
    else
        tremor_or22(i,1:5001)=NaN;
        tremor_or2(i,1:5001)=NaN;
        tremor_k(i,1)=NaN;
    end
end

clear tt
tt=[]; 
k=1;

tt=NaN(20,12);

for i=1:12
    tt(1:sum(xx==i),i)=tremor_k(find(xx==i));
end

%% if you are using it on 2 Hz data
% indexes would be where your pulse was delivered
% this is again very similar to before , i unwrap phase to get frequency
% change from the slow 
% sampling frequency is 1000 Hz , so i start it 200 ms before stim onset
% and look till 200 ms after stimulation and check the deviation from
% "projected" unwrapped phase of a perfect sinusoid vs our signal
% in your case it would also be interesting to see when or if the perfect
% sinuoid deviates from observed data (i.e. when does the frequency change
% start . frequency of the perfect sinusoid is derived from last one second
 indexes=start;
for i=1:length(indexes)
    variable1=((phase(indexes(i)-200)+(0:1:400)*2*pi/(1000./mean(frequency(indexes(i)-1000:indexes(i))))));
    variable2=(unwrap(phase(indexes(i)-200:indexes(i)+200)));
    
    variable=variable1(end);
    if variable>pi
        variable=variable-2*pi;
    elseif variable<-pi
        variable=variable+2*pi;
    end
    d(i)=variable1(end)-variable2(end);
end

d(find(d>pi))=d(find(d>pi))-2*pi;
d(find(d<-pi))=d(find(d<-pi))+2*pi;

ind=-pi:(2*pi)/12:(pi-1/12); % divides stimulation phase to 12 bins, I would start first dividing it into two bins ( 0-180 and 180-360
d2=NaN(12,100);
for i=1:length(ind)
    d2(i,1:length(find(phase(indexes)>ind(i) & phase(indexes)<(ind(i)+2*pi/12))))=(d(find(phase(indexes)>ind(i) & phase(indexes)<(ind(i)+2*pi/12))));
end
figure()
plot(ind+2*pi/20,(d2'),'.')
hold on
plot(ind+2*pi/20,nanmedian(d2')) 