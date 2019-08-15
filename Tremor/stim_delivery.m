clear all
iii=[1 2 3 4 5 8 10 11 13];

for numb=1
%     :length(iii);
    clearvars -except iii numb ttall ampall counts counts1
load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))
% load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))

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

time=0:1/samplerateold:(size(data,2)-1)/samplerateold;

%%% downsample

ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremor2(1:size(ts1.data,3))=ts1.data;
samplerate=1000;


new=find(data(2,:)>4);
end_point=[new(find(diff(new)>100000)) new(end)];
start_point=[new(1) new(find(diff(new(1,1:end-1))>100000)+1)];

t=0:1/1000:(numel(start_point(1):end_point(1))-1)/1000;
do=data(2,start_point(1):end_point(1));
d=zeros(1,numel(start_point(1):end_point(1)));
id_s=round(1:length(d)./12:length(d));
id_end=id_s+50000;
for ih=1:12
d(id_s(ih):id_end(ih))=5;
end
plot(do)
hold on
plot(d,'r','LineWidth',1.5)


for i=1
 n=find(do(id_s(i):id_end(i))>4);
end


end_point=[n(find(diff(n)>100)) new(end)];
start_point=[(1) n(find(diff(n(1,1:end-1))>100)+1)];

% plot(time,data(2,:))
% hold on
% plot(time(new),data(2,new),'r.')
% plot (time(end_point),data(2,end_point),'k.')
% plot(time(start_point),data(2,start_point),'y.')





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

end