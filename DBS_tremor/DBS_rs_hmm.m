

clearvars -except data Fs segmentb segmente a b ns fs iii numb in2 time_n
load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))

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


ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts2=resample(ts,0:1/Fs:((size(data,2)-1)/samplerateold),'linear');
data2(:,1:size(ts2.data,3))=ts2.data;
time_n=0:1/Fs:(size(data2(in,:),2)-1)/Fs;

new=find(data2(2,:)>4);
difp=find(diff(new)>500);  %(100000)*Fs./1000
sp_1=[new(1) new(difp+1)];
ep_1=[new(difp) new(end)];

plot(time_n,data2(2,:))
hold on
plot(time_n(sp_1),data2(2,sp_1),'ro')

%%------------------------

Fs=20;

ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts2=resample(ts,0:1/Fs:((size(data,2)-1)/samplerateold),'linear');
data2(:,1:size(ts2.data,3))=ts2.data;
time_n=0:1/Fs:(size(data2(in,:),2)-1)/Fs;

tremor_or=filtfilt(b,a,data2(in,:));
tremor_or=zscore(tremor_or);
envelope=abs(hilbert(tremor_or));

t=300; %(15000.*Fs)./1000;
rs=[];
for tr=1:length(sp_1)
%     if (run/2)>t && ((run/2)+t)<length(tremor_or)
        run=round(length(sp_1(tr):ep_1(tr))./2,0);
        rs=[rs (tremor_or(sp_1(tr-((run-t)):((numel(run)/2)+t)-1))];
%     end
end

%%%%% check eppochs!!