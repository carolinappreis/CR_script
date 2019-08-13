clear all
 iii=[1 2 3 4 5 6 8 10 11 13];

PC=[70 66 47 47 47 50 50 50 50 55];
A1={([1 3 6 8 12 18 23 27 30 32]);[];[];([2 3 5 6 7]);([1 2 4 6 8 9 10 11]);[];[];([1:9 15]);([2 4 7:10 13:15 22 25]);[]};
B1={([2 5 7 11 17 22 26 29 31 34]);[];[];([2 3 5 6 7]);([1 2 4 6 7 9 10 11 12]);[];[];([1:9 15]);([2 5 7 8 9 12 13 14 19 22 25]);[]};
for numb=length(iii);
% load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
% load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))
 load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(iii(numb)),'_baseline.mat'))
 load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))

in2=1; % analysing the "main tremor axis"
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

phasedetection;
% bar(0:30:330,100.*nanmedian(tt)) 
% box('off')% x axis phase - y axis percent change in
% tremor severity at the end of 5 seconds with respect to severity right
% before stimulation began

% stimout=tt;
stim(numb,:)=nanmedian(tt);
% save stim stimout

data=SmrData.WvData;
rep=10; % number of trials for random stim - please enter for each patient
clearvars -except Fpeak in2 in rep SmrData data nostimout iii numb PC A1 B1 iii stim nostim xx

samplerateold=SmrData.SR;
time=0:1/samplerateold:(size(data,2)-1)/samplerateold;
tremor=data(in,:);
ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremor2(1:size(ts1.data,3))=ts1.data;
samplerate=1000;

if (Fpeak-2)>=1
    [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
    [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end

tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
dummy=hilbert(tremor_or);
envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));

[b,a]=butter(2,[0.1/(0.5*samplerate) ],'low'); %15
C=(filtfilt(b,a,envelope));
 A=prctile(C,PC(numb));

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
    
    if isempty (A1{numb,1})
        AA=ind_s2;
        BB=ind_e2;
    else
        AA=ind_s2(A1{numb,1});
        BB=ind_e2(B1{numb,1});
    end
    
    if numb==5
        AA= [1 ind_s2(A1{numb,1})];
        BB=ind_e2(B1{numb,1});
    end
    
    
    segmentb=AA;
    segmente=BB;
  
% plot(C)
% hold on
% plot(AA,C(AA),'r.')
% plot(BB,C(BB),'b.')

%%% removal of time segments [ < mean - std for more than 10 seconds ]
handup=[];

for i=1:length(segmentb)
    handup=[handup segmentb(i):segmente(i)]; %#ok<*AGROW>
end
clear i
handup=sort(handup,'ascend');

unstable=find(envelope<(mean(envelope(handup))-std(envelope(handup))));

if ~isempty(unstable)
    change_e=[unstable(find(diff(unstable)~=1)) unstable(end)]; %#ok<*FNDSB>
    change_b=[unstable(1) unstable(find(diff(unstable)~=1)+1)];
    
    change_ind=find((change_e-change_b)>10000);
    
    unstable2=[];
    if ~isempty(change_ind)
        change_e2=change_e(change_ind);
        change_b2=change_b(change_ind);
        for i=1:length(change_ind)
            unstable2=[unstable2 change_b2(i):change_e2(i)]; %#ok<*AGROW>
        end
    end
    
    [unstable3 ia ib]=intersect(handup,unstable2);
else
    unstable3=[];
end

%%% analysis

if (Fpeak-2)>=1
    [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
    [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end

tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
dummy=hilbert(tremor_or);
envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
phase=angle(dummy);
frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),251))';

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
[b,a]=butter(2,[2/(0.5*samplerate) 7/(0.5*samplerate)],'bandpass'); %15
tremorxf=filtfilt(b,a,tremorx);
tremoryf=filtfilt(b,a,tremory);
tremorzf=filtfilt(b,a,tremorz);

for j=1:length(segmentb)
    x=[tremorxf(segmentb(j):segmente(j));tremoryf(segmentb(j):segmente(j));tremorzf(segmentb(j):segmente(j))];
    [pc,score,latent,tsquare] = pca(x');
    xxx(j,1:3)=pc(1:3,1);
    ma(j)=(find(abs(xxx(j,1:3))==max(abs(xxx(j,1:3)))));
end

for i=1
    for j=1:5e4
        ix=randi(length(segmentb),1);
         while (segmentb(ix)+1000)<(segmente(ix)-5000)
        segment=randi([segmentb(ix)+1000 segmente(ix)-5000],1);  
        begin3=segment;
        end3=floor(begin3+5*samplerate);
         end
        while ~isempty(intersect(unstable3,begin3:end3))
            segment=randi([segmentb(ix)+1000 segmente(ix)-5000],1);
            begin3=segment;
            end3=floor(begin3+5*samplerate);
        end
        baseline3(i,j)=(mean(envelope(end3-1000:end3))-mean(envelope(begin3-1000:begin3)))./mean(envelope(begin3-1000:begin3)); %#ok<*SAGROW> %
        baseline4(i,j)=(mean(frequency(end3-1000:end3))); %#ok<*SAGROW>
        
    end
end

for i=1:1e6
    dum=baseline3(randi(5e4,1,rep));
    dum2=dum;
    p(i)=nanmedian(dum2);
end
nostim(numb,:)=p;

for i=1:12
    dum=baseline3(randi(5e4,1,rep));
    dum2=dum;
    nostimout(numb,i)=nanmedian(dum2);
end
clearvars -except nostimout iii numb PC A1 B1 iii stim nostim 
end
% ANS_group=nostimout; clear nostimout
% AS_group=stim;
% clearvars  -except ANS_group AS_group
% cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
% save 'A_group'
NS=nostimout; 
S=stim;
idv_NS=nostim;
clearvars  -except NS S idv_NS
% cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
save 'A_group13'
