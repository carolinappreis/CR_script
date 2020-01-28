clear all
load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(11),'_RS.mat'))
iii=[2 3 4 5 8 10 11 13 16 17]; 
main=[1 1 3 1 3 3 3 3 1 1];
in2=main(find(iii==11)); % analysing the "main tremor axis"

if in2==1
    in=3;
elseif in2==2 % other axis 1
    in=5;
elseif in2==3 % other axis 2
    in=6;
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

Fpeak=frange(find(Pxxrange==max(Pxxrange))); %#ok<*FNDSB>


clearvars -except Fpeak in2 main iii
load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(11),'_baseline.mat'))

if in2==1
    in=3;
elseif in2==2 % other axis 1
    in=5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
elseif in2==3 % other axis 2
    in=6;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
end
data=SmrData.WvData;
samplerateold=SmrData.SR;
tremor=(data(in,:));
addon=92; addon_end=35;


ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
samplerate=1000;
tt=0:1/samplerate:(size(ds_data,2)-1)/samplerate;
tre_3=ds_data([3 5 6],:);

samplerate=1000;

if (Fpeak-2)>=1
    [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
    [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end

tf_3=filtfilt(b,a,tre_3')*10*9.81/0.5;
% tremor_or=zscore(tremor_or);
dummy=(hilbert(tf_3))';
envelope=abs(dummy);
zenv=(abs(hilbert(zscore(tf_3))))';
phase=angle(dummy);
frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';


[b,a]=butter(2,[0.1/(0.5*samplerate) ],'low'); %15
C=(filtfilt(b,a,envelope(main(find(iii==11)),:)));
A=mean(C);
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
AA=ind_s2(find((ind_e2-ind_s2)>10000));
BB=ind_e2(find((ind_e2-ind_s2)>10000));

close all
plot(C)
hold on
plot(AA,C(AA),'r.')
plot(BB,C(BB),'r.')