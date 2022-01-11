load('P06_RS_cursos.mat')

in2=3; % analysing the "main tremor axis"

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

timeold=0:1/samplerateold:(size(data,2)-1)/samplerateold;

%% downsample

ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
tremor2(1:size(ts1.data,3))=ts1.data;
samplerate=1000;
time=0:1/samplerate:(size(tremor2,2)-1)/samplerate;

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

diff_xx=diff(xx);
diff_in=diff(indexes4);

xx_copy=xx;
indexes4_copy=indexes4;
indexes3_copy=indexes3;

for i=1:length(diff_xx)
    if diff_xx(i)==0 && diff_in(i)<5*samplerateold
        xx_copy(i+1)=NaN;
        indexes4_copy(i+1)=NaN;
        indexes3_copy(i)=NaN;
    end
end

indexes4_2=indexes4_copy(~isnan(indexes4_copy));
xx_2=xx_copy(~isnan(xx_copy));
indexes3_2=indexes3_copy(~isnan(indexes3_copy));

start=floor((indexes4_2./samplerateold)*samplerate)+addon;
ending=floor((indexes3_2./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);

indexesblock_b=floor([index(1) index(find(diff(index)./samplerateold > 20)+1)].*samplerate./samplerateold);
indexesblock_e=floor([index(find(diff(index)./samplerateold > 20)) index(end)].*samplerate./samplerateold);

%% when patient's hand is up
handup=[];
for i=1:length(indexesblock_b)
    handup=[handup indexesblock_b(i):indexesblock_e(i)]; %#ok<*AGROW>
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


%% re - estimate tremor characteristics

dummy=hilbert(tremor_or);
envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
phase=angle(dummy);
frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';

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
end
clear j

tremor_or2=NaN(length(start),1);

[b,a]=butter(2,[0.1/(2*samplerate) 2/(2*samplerate)],'bandpass');
dc_level(1,:)=filtfilt(b,a,tremorx);
dc_level(2,:)=filtfilt(b,a,tremory);
dc_level(3,:)=filtfilt(b,a,tremorz);

figure()
time=0:1/samplerate:(size(envelope,2)-1)/samplerate;
for i=1:3
    z_diff_dc(i,:)=zscore(diff(dc_level(i,:)));
    hold on
    plot(time(1:end-1),z_diff_dc(i,:))
    hold on
    plot(time(start),z_diff_dc(i,start),'r.')
    plot(time(ending),z_diff_dc(i,ending),'g.')
end

% clear start2 ending2
% start2=start;
% ending2=ending;
% for i=1:length(start)
%     %     if sum(abs(z_diff_dc(:,ending(i)))<1.96 & abs(z_diff_dc(:,start(i)))<1.96)==3
%     %         start2(i)=start(i);
%     %         ending2(i)=ending(i);
%     %     elseif sum(abs(z_diff_dc(:,ending(i)))>1.96 & abs(z_diff_dc(:,start(i)))<1.96)>0
%     %         start2(i)=start(i);
%     %         ending2(i)=ending(i)-12000;
%     % %     elseif sum(abs(z_diff_dc(:,ending(i)))<1.96 & abs(z_diff_dc(:,start(i)))>1.96)>0
%     % %         start2(i)=start(i)+12000;
%     % %         ending2(i)=ending(i);
%     %     elseif sum(abs(z_diff_dc(:,ending(i)))>1.96 & abs(z_diff_dc(:,start(i)))>1.96)>0
%     start2(i)=start(i);%+12000;
%     ending2(i)=ending(i)-12000;
%     %     end
% end
% 
% figure()
% time=0:1/samplerate:(size(envelope,2)-1)/samplerate;
% for i=1:3
%     z_diff_dc(i,:)=zscore(diff(dc_level(i,:)));
%     hold on
%     plot(time(1:end-1),z_diff_dc(i,:))
%     hold on
%     plot(time(start2),z_diff_dc(i,start2),'r.')
%     plot(time(ending2),z_diff_dc(i,ending2),'g.')
% end
% 
% start3=start2;
% start3(2:4)=start2(1);
% start3(6:7)=start2(5);
% start3(9)=start3(8);
% start3(11:12)=start3(10);
% start3(14:15)=start3(13);
block_interest=[];
angle_block=[];
D=[20000 46000 26000 23000 3000 5000 3000 17000];
for blocks=1:2:16
    block_interest=[block_interest indexesblock_b(blocks):indexesblock_e(blocks)];
%     eval(['AA=(load(' '''rs0' int2str(blocks/2) '_M.txt''' '));']);    
%     AA(find(AA(:,1)==1),1)=NaN;
%     AA(find(AA(:,2)==1),2)=NaN;
%     AA(find(AA(:,3)==1),3)=NaN;
%     AA(find(AA(:,1)==-1),1)=NaN;
%     AA(find(AA(:,2)==-1),2)=NaN;
%     AA(find(AA(:,3)==-1),3)=NaN;
%     BB(:,1)=AA(~isnan(AA(:,1)),1);
%     BB(:,2)=AA(~isnan(AA(:,2)),2);
%     BB(:,3)=AA(~isnan(AA(:,3)),3);
%     theta=atan2((BB(:,2)-mean(BB(:,2))),(BB(:,1)-mean(BB(:,1))));
%     [theta_fix,ty]=resample(theta,BB(:,3),1);
%     dum=dc_level(1,indexesblock_b(blocks):indexesblock_e(blocks));
%     dum2=(dum-mean(dum(1:35000)))./std(dum(1:35000));
%     theta_aligned=(theta_fix(D(blocks/2):(D(blocks/2)+length(dum)-1)))';
%     angle_block=[angle_block theta_aligned];
%     clear AA BB theta theta_fix ty dum dum2 Xa Ya theta_aligned
end

for i=1:length(start)
    if intersect(start(i),block_interest) & envelope(start(i))>median(envelope(block_interest))
%         [c, ia, ib]=intersect(start(i),block_interest);
%         if ib+5000<length(angle_block)
%         s_angle=find(angle_block(ib:ib+5000)>0 & angle_block(ib:ib+5000)<pi);
%         if length(s_angle>1000)
%             new_start=start(i)+s_angle(1);
%             new_end=start(i)+s_angle(end);
%             tremor_or2(i,1)=(mean(envelope(new_end-1000:new_end))-mean(envelope(new_start-1000:new_start)))/mean(envelope(new_start-1000:new_start));
%         else
%             tremor_or2(i,1)=NaN;
%         end
%         else
%             tremor_or2(i,1)=NaN;
%         end
        
       tremor_or2(i,1)=(mean(envelope(ending(i)-1000:ending(i)))-mean(envelope(start(i)-1000:start(i))))/mean(envelope(start(i)-1000:start(i)));
    else
       tremor_or2(i,1)=NaN;
    end
end

clear tt
k=1;
tt=NaN(20,12);

for i=1:12
    tt(1:sum(xx_2==i),i)=tremor_or2(find(xx_2==i));%./baseline2((i+11)/12);
end