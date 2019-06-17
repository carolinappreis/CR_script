% first convert baseline smr to mat by running SMR_File_To_Mat
% then run peripheraltremor for the random phase stimulation file
% and then run this script without clearing the workspace
clear all
cd('C:\Users\creis\Documents\MATLAB\pt_data_periphstim')
SMR_File_To_Mat;
data=SmrData.WvData;
rep=10; % number of trials for random stim - please enter for each patient
clearvars -except Fpeak in2 in rep SmrData data

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

segment_nostim % please check that the output makes sense - sometimes could be too 
%sensitive when tremor is not large enough amplitude so you may need to
%enter the start and end points of when the patient put their hand up by
%hand - if so - it is samples or in milliseconds 

%% segments

segmentb=AA; 
segmente=BB; 

%% removal of time segments [ < mean - std for more than 10 seconds ]
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
    
end

%% analysis

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
save nostim nostimout