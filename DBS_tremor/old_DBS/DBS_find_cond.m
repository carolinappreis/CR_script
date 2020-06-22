
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

%-------

new=find(data(2,:)>4);
difp=find((diff(new))>100000);
ep_1=[new(difp) new(end)];
sp_1=[new(1) new(difp+1)];


if numb==1
    dt1_s=[sp_1(2:2:end)];dt1_e=[ep_1(2:2:end)];
    dt2_s=sp_1(3:2:end); dt2_e=ep_1(3:2:end);
    sp=sp_1(1,2:end);
    ep=ep_1(1,2:end);
end
%     %
%             plot(time,data(4,:))
%         hold on
%         plot(time(sp),data(4,sp),'r.')
%         plot(time(ep),data(4,ep),'k.')


for ik=1:length(sp) %%find double start and end points in a stimulation run
    
    s=(find(([indexes4>=sp(ik)]+[indexes4<=ep(ik)])==2));
    e=(find(([indexes3>=sp(ik)]+[indexes3<=ep(ik)])==2));
    tks=(find(diff(xx(s))==0))+1;
    tke=(find(diff(xx(e))==0));
    
    indexes4(s(tks))=NaN;
    indexes3(e(tke))=NaN;
    xx(e(tke))=NaN;
    
end

%%%%% find runs with trigering issues (too few, too many pulses)
%     th1=(Fpeak*5)./2;
%     th2=(Fpeak*5)+round((Fpeak*5)./5);
%         for it=1:length(indexes4)
%             if numel(index(find(index==indexes4(it)):find(index==indexes3(it))))>=th1 && numel(index(find(index==indexes4(it)):find(index==indexes3(it))))<=th2
%                 indexes4(it)=indexes4(it);
%                 indexes3(it)=indexes3(it);
%                 xx(it)=xx(it);
%             else
%                 indexes4(it)=NaN;
%                 indexes3(it)=NaN;
%                 xx(it)=NaN;
%             end
%         end
%%%%%%%%%%%%%%%%
indexes4=indexes4(~isnan(indexes4));
indexes3=indexes3(~isnan(indexes3));
xx=xx(~isnan(xx));


start1=[];
ending1=[];
xx1=[];
for il=1:length(dt1_s)
    start1=[start1 indexes4(find(([indexes4>=dt1_s(il)]+[indexes4<=dt1_e(il)])==2))];
    ending1=[ending1 indexes3(find(([indexes3>=dt1_s(il)]+[indexes3<=dt1_e(il)])==2))];
    xx1=[xx1 xx(find(([indexes4>=dt1_s(il)]+[indexes4<=dt1_e(il)])==2))];
end

start2=[];
ending2=[];
xx2=[];
for il=1:length(dt2_s)
    if numb==1&&il==3 | numb==1&&il==9
        dums=indexes4(find(([indexes4>=dt2_s(il)]+[indexes4<=dt2_e(il)])==2));
        start2=[start2 dums(1,3:end)]; clear dums
        dume=indexes3(find(([indexes3>=dt2_s(il)]+[indexes3<=dt2_e(il)])==2));
        ending2=[ending2 dume(1,3:end)]; clear dume
        dumx=xx(find(([indexes4>=dt2_s(il)]+[indexes4<=dt2_e(il)])==2));
        xx2=[xx2 dumx(1,3:end)];clear dumx
    else
        dums=indexes4(find(([indexes4>=dt2_s(il)]+[indexes4<=dt2_e(il)])==2));
        start2=[start2 dums]; clear dums
        dume=indexes3(find(([indexes3>=dt2_s(il)]+[indexes3<=dt2_e(il)])==2));
        ending2=[ending2 dume]; clear dume
        dumx=xx(find(([indexes4>=dt2_s(il)]+[indexes4<=dt2_e(il)])==2));
        xx2=[xx2 dumx];clear dumx
        
    end
end

figure()
plot(time,data(4,:))
hold on
plot(time(index),data(4,index),'r.')
plot(time(start1),data(4,start1),'ko')
plot(time(ending1),data(4,ending1),'bo')


clear start ending
start{1,1}=floor((start1./samplerateold)*samplerate)+addon;
ending{1,1}=floor((ending1./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
start{2,1}=floor((start2./samplerateold)*samplerate)+addon;
ending{2,1}=floor((ending2./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
clear xx
xx{1,1}=xx1;
xx{2,1}=xx2;