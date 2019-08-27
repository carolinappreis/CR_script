
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

Fpeak=frange(find(Pxxrange==max(Pxxrange))); %#ok<*FNDSB>

if (Fpeak-2)>=1
    [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
    [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end

tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;

dummy=hilbert(tremor_or);
envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));

th1=(Fpeak*5)./2;
th2=(Fpeak*5)+round((Fpeak*5)./5);
%-------

new=find(data(2,:)>4);
difp=find((diff(new))>100000);
ep=[new(difp) new(end)];
sp=[new(1) new(difp+1)];

%     plot(time,data(4,:))
%     hold on
%     plot(time(sp),data(4,sp),'r.')
%     plot(time(ep),data(4,ep),'k.')

for ik=1:length(sp) %%find double start and end points in a stimulation run
    
    s=(find(([indexes4>=sp(ik)]+[indexes4<=ep(ik)])==2));
    e=(find(([indexes3>=sp(ik)]+[indexes3<=ep(ik)])==2));
    tks=(find(diff(xx(s))==0))+1;
    tke=(find(diff(xx(e))==0));
    indexes4(s(tks))=NaN;
    indexes3(e(tke))=NaN;
    xx(e(tke))=NaN;
    
end

for it=1:length(indexes4) %% find runs with trigering issues (too few, too many pulses)
        if numel(index(find(index==indexes4(it)):find(index==indexes3(it))))>=th1 && numel(index(find(index==indexes4(it)):find(index==indexes3(it))))<=th2
            indexes4(it)=indexes4(it);
            indexes3(it)=indexes3(it);
            xx(it)=xx(it);
        else
            indexes4(it)=NaN;
            indexes3(it)=NaN;
            xx(it)=NaN;
        end
end

indexes4=indexes4(~isnan(indexes4));
indexes3=indexes3(~isnan(indexes3));
xx=xx(~isnan(xx));

clear start ending
start=floor((indexes4./samplerateold)*samplerate)+addon;
ending=floor((indexes3./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
% 
% plot(time,data(4,:))
% hold on
% plot(time(index),data(4,index),'r.')
% plot(time(indexes4),data(4,indexes4),'ko')
% plot(time(indexes3),data(4,indexes3),'bo')