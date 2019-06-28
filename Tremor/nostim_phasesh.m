clear all
 cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline')
% cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline')
load ('P01_baseline.mat')

 cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim')
% cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim')
load ('P01_randstim_cursos.mat')
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
cd('C:\Users\creis\Documents\GitHub\CR_script\Tremor')

run ('priph_tremor_phaseshift.m');

stimout=tt;
% save stim stimout

data=SmrData.WvData;
rep=10; % number of trials for random stim - please enter for each patient
clearvars -except Fpeak in2 in rep SmrData data stimout

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
% A=mean(C);
A=prctile(C,50);
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
% 
plot(C)
hold on
plot(ind_s2,C(ind_s2),'r.')
plot(ind_e2,C(ind_e2),'b.')

AA=ind_s2;
BB=ind_e2;

% pt 1
A=prctile(C,70);
AA=ind_s2([1 3 6 8 12 18 23 27 30 32]);
BB=ind_e2([2 5 7 11 17 22 26 29 31 34]);
AA=ind_s2([1 3 6 8 12 18 23 27 30 ]);
BB=ind_e2([2 5 7 11 17 22 26 29 31 ]);
%pt2
%A=prctile(C,66);
%pt3
%A=prctile(C,47);
%pt4
%A=prctile(C,47);
% AA=ind_s2([2 3 5 6 7]);
% BB=ind_e2([2 3 5 6 7]);
%pt5
%A=prctile(C,47);
% AA= [1 ind_s2([1 2 4 6 8 9 10 11])];
% BB= ind_e2([1 2 4 6 7 9 10 11 12]);
%pt6
%A=prctile(C,50);
%pt8
%A=prctile(C,50);
%pt10
%A=prctile(C,50);
% AA= ind_s2([1:9 15]);
% BB= ind_e2([1:9 15]);
%pt11
% AA= ind_s2([2 4 7:10 13:15 22 25]);
% BB= ind_e2([2 5 7 8 9 12 13 14 19 22 25]);

% AA=ind_s2(find((ind_e2-ind_s2)>10000));
% BB=ind_e2(find((ind_e2-ind_s2)>10000));
close all
plot(C)
hold on
plot(AA,C(AA),'r.')
plot(BB,C(BB),'b.')

% segment_nostim % please check that the output makes sense - sometimes could be too
%sensitive when tremor is not large enough amplitude so you may need to
%enter the start and end points of when the patient put their hand up by
%hand - if so - it is samples or in milliseconds

%%% segments

segmentb=AA;
segmente=BB;

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
tremor=(data(6,:));% %score(:,1)';%
ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremory(1:size(ts1.data,3))=ts1.data;
tremor=(data(7,:));% %score(:,1)';%
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
        segment=randi([segmentb(ix)+1000 segmente(ix)-5000],1);
        begin3=segment;
        end3=floor(begin3+5*samplerate);
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

nostimout=p;
clearvars -except stimout nostimout C AA BB
save 'p011_stimnosim.mat'

% load stim
% load nostim
% d=nanmedian(stimout);
% 
% upperthreshold=prctile(nostimout,99.7917); %bonferroni corrected for 12 comparisons
% lowerthreshold=prctile(nostimout,0.2083);
% 
% find(d>=upperthreshold | d<=lowerthreshold)
% bar(d) 
% hold on 
% rr(1:12)=upperthreshold;
% rr2(1:12)=lowerthreshold;
% plot(rr,'r--')
% plot(rr2,'r--')
% xlim([0.5 12.5])
% box('off')
% xticklabels({'0' '30' '60' '90' '120' '150' '180' '210' '240' '270' '300' '330'})
% 
% 
% % 
% % close all
% % figure(2)
% % fig=gcf;
% % fig.Color=[1 1 1];
% % bar(100.*nanmedian(tt))
% % hold on
% % stem((100.*tt)')
% % xticklabels({'0','30','60','90','120','150','180','210','240','270','300','330'})
% % ylim([(-max((max(tt)))-0.1).*100 (max(max(tt))+0.1).*100])
% % box('off')
% % title ('P1')
