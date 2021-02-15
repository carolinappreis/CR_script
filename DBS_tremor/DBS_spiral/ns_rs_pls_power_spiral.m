clear; close
cohort = [1 3 4 6];

for iii =  2:length(cohort)
   
    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_NS.mat'))
    
    
    data_raw=SmrData.WvData;
    samplerateold=SmrData.SR;
    ts=timeseries(data_raw,0:(1/samplerateold):((size(data_raw,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data_raw,2)-1)/samplerateold),'linear');
    ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    tt=0:1/samplerate:(size(ds_data,2)-1)/samplerate;
    data_t=ds_data([3 6 7],:);
    
    NS_BE_S
    out.b_trials=hu{iii,:};
    out.e_trials=hd{iii,:};
    
    
%         figure(1)
%         time=1:length(data_t(1,:));
%         plot(time,data_t(1,:))
%         hold on
    handup = [];
    for i = 1:length(out.b_trials)
        handup = [handup out.b_trials(i):out.e_trials(i)]; %#ok<*AGROW>

    end
    clear i
    handup = sort(handup,'ascend');
    
    for aa = 1:3
        [Pxx,F] = pwelch(data_t(aa,handup), samplerate, [], round(samplerate), samplerate);
        frange = F(3:10);
        Pxxrange = Pxx(3:10);
        Freqpeak(aa,:) = frange(find(Pxxrange == max(Pxxrange)));
        Ppeak(aa,:) = max(Pxxrange);
        ps_curves(aa,:) = Pxx;
    end
    
    peak_ax = [(Freqpeak(find(Ppeak == max(Ppeak)))) (find(Ppeak == max(Ppeak)))];
    
    plotp(1,iii-1)=Ppeak(peak_ax(2));
   pcurve(1,iii-1,:)=ps_curves(peak_ax(2),:);
    clearvars -except iii cohort plotp pcurve
end

clearvars -except cohort plotp pcurve
for iii =2: 4
    %1:length(cohort)
    
    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_RS.mat'))
    
    
    in2 = 1; % analysing the "main tremor axis"
    
    if in2 == 1
        in = 3;
    elseif in2 == 2 % other axis 1
        in = 6;
    elseif in2 == 3 % other axis 2
        in = 7;
    end
    
    %
    data = SmrData.WvData;
    addon=92; addon_end=35;
    
    
    %%
    samplerateold = SmrData.SR;
    tremor = (data(in, :));
    
    %%%%%%%%%%%%%%%%%%%
    
    time = 0:1/samplerateold:(size(data, 2)-1)/samplerateold;
    
    % downsample
    
    ts = timeseries(tremor,0:(1/samplerateold):((size(data,2)-1) / samplerateold));
    ts1 = resample(ts,0:0.001:((size(data,2)-1) / samplerateold), 'linear');
    tremor2(1:size(ts1.data, 3)) = ts1.data;
    samplerate = 1000;
    
    % determine stimulation time points
    index = [];
    for i = 2:size(data, 2)-1
        if data(2,i-1)<2.5 && data(2,i)>2.5
            index = [index i];
        end
    end
    clear i
    
    % Find trigger
    indexes4 = [index(1) index(find(diff(index) ./ samplerateold > 0.95)+1)];
    indexes3 = [index(find(diff(index) ./ samplerateold > 0.95)) index(end)];
    
    dd2 = round(data(4, :)*100) ./ 100;
    for i = 1:length(indexes4)
        xx(i) = round(dd2(indexes4(i)) ./ 0.1); %#ok<*SAGROW>
    end
    clear i
    
    start = floor((indexes4 ./ samplerateold)*samplerate)+addon;
    ending = floor((indexes3 ./ samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
    
    % when patient's hand is up
    handup = [];
    for i = 1:length(start)
        handup = [handup start(i):ending(i)]; %#ok<*AGROW>
    end
    clear i
    handup = sort(handup,'ascend');
    
    
    % tremor characteristics
    [Pxx, F] = pwelch(tremor2(handup), samplerate, [], samplerate, samplerate);
    
    frange = F(3:10);
    Pxxrange = Pxx(3:10);
    
    Fpeak = frange(find(Pxxrange == max(Pxxrange))); %#ok<*FNDSB>
    
    if (Fpeak-2) >= 1
        [b, a] = butter(2, [(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
    else
        [b, a] = butter(2, [(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
    end
    
    tremor_or = filtfilt(b, a, tremor2)*10*9.81/0.5;
    dummy = hilbert(tremor_or);
    % phase=angle(dummy);
    frequency = (smooth((1000/(2*pi))*diff(unwrap(angle(dummy))), 500))';
    
    tremor = (data(3, :));% %score(:,1)';%
    ts = timeseries(tremor, 0:(1/samplerateold):((size(data, 2)-1)/samplerateold));
    ts1 = resample(ts, 0:0.001:((size(data,2)-1)/samplerateold), 'linear');
    tremorx(1:size(ts1.data, 3)) = ts1.data;
    tremor = (data(6, :));% %score(:,1)';%
    ts = timeseries(tremor, 0:(1/samplerateold):((size(data, 2)-1)/samplerateold));
    ts1 = resample(ts, 0:0.001:((size(data, 2)-1)/samplerateold), 'linear');
    tremory(1:size(ts1.data, 3)) = ts1.data;
    tremor = (data(7, :));% %score(:,1)';%
    ts = timeseries(tremor, 0:(1/samplerateold):((size(data, 2)-1)/samplerateold));
    ts1 = resample(ts, 0:0.001:((size(data, 2)-1)/samplerateold), 'linear');
    tremorz(1:size(ts1.data, 3)) = ts1.data;
    tremorxf = filtfilt(b, a, tremorx);
    tremoryf = filtfilt(b, a, tremory);
    tremorzf = filtfilt(b, a, tremorz);
    envelope = [abs(hilbert(tremorxf)); abs(hilbert(tremoryf)); abs(hilbert(tremorzf))];
    phase = [angle(hilbert(tremorxf)); angle(hilbert(tremoryf)); angle(hilbert(tremorzf))];
    
    new = find(data(2, :) > 4);
    difp = find((diff(new)) > 100000); % are you trying to threshold at 9.6 seconds?
    ep_1 = [new(difp) new(end)];
    sp_1 = [new(1) new(difp+1)];
    
    %%% input start all trial
    start_t=1;
    
    sp=sp_1(1,start_t:end);
    ep=ep_1(1,start_t:end);
    
    %%% posture/spiral trials
    dt1_s=[sp_1(start_t:2:end)];dt1_e=[ep_1(start_t:2:end)];
    dt2_s=sp_1(start_t+1:2:end); dt2_e=ep_1(start_t+1:2:end);
    
    %         time=1:length(data(1,:));
    %         plot(time,data(4,:))
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
    
    
    % % %     %%%% find runs with trigering issues (too few, too many pulses)
    % % %     th1=(Fpeak*5)./2;
    % % %     th2=(Fpeak*5)+round((Fpeak*5)./5);
    % % %     for it=1:length(indexes4)
    % % %         if numel(index(find(index==indexes4(it)):find(index==indexes3(it))))>=th1 && numel(index(find(index==indexes4(it)):find(index==indexes3(it))))<=th2
    % % %             indexes4(it)=indexes4(it);
    % % %             indexes3(it)=indexes3(it);
    % % %             xx(it)=xx(it);
    % % %         else
    % % %             indexes4(it)=NaN;
    % % %             indexes3(it)=NaN;
    % % %             xx(it)=NaN;
    % % %         end
    % % %     end
    % % %     %%%%%%%%%%%%%%%
    indexes4=indexes4(~isnan(indexes4));
    indexes3=indexes3(~isnan(indexes3));
    xx=xx(~isnan(xx));
    
    
    start1=[];
    ending1=[];
    xx1=[];
    for il=1:length(dt1_s)
        start1=[start1 indexes4(find(([indexes4>=dt1_s(il)]+[indexes4<=dt1_e(il)])==2))]; % intersect([1 2 3],[3 4 5])
        ending1=[ending1 indexes3(find(([indexes3>=dt1_s(il)]+[indexes3<=dt1_e(il)])==2))];
        xx1=[xx1 xx(find(([indexes4>=dt1_s(il)]+[indexes4<=dt1_e(il)])==2))];
    end
    
    start2=[];
    ending2=[];
    xx2=[];
    for il=1:length(dt2_s)
        dums=indexes4(find(([indexes4>=dt2_s(il)]+[indexes4<=dt2_e(il)])==2));
        start2=[start2 dums]; clear dums
        dume=indexes3(find(([indexes3>=dt2_s(il)]+[indexes3<=dt2_e(il)])==2));
        ending2=[ending2 dume]; clear dume
        dumx=xx(find(([indexes4>=dt2_s(il)]+[indexes4<=dt2_e(il)])==2));
        xx2=[xx2 dumx];clear dumx
    end
    
    
    %         figure()
    %         plot(time,data(4,:))
    %         hold on
    %         plot(time(index),data(4,index),'r.')
    %         plot(time(start1),data(4,start1),'ko')
    %         plot(time(ending1),data(4,ending1),'bo')
    
    
    clear start ending xx
    
    stt=floor((start2./samplerateold)*samplerate)+addon;
    ett=floor((ending2./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
    
    handup=[];
    for i=1:length(stt)
        handup=[handup stt(i):ett(i)];
    end

    clear i
    handup = sort(handup,'ascend');
    
  data_t=[tremorx;tremory;tremorz];
    
    for aa = 1:3
        [Pxx,F] = pwelch(data_t(aa,handup), samplerate, [], round(samplerate), samplerate);
        frange = F(3:10);
        Pxxrange = Pxx(3:10);
        Freqpeak(aa,:) = frange(find(Pxxrange == max(Pxxrange)));
        Ppeak(aa,:) = max(Pxxrange);
        ps_curves(aa,:) = Pxx;
    end
    
    peak_ax = [(Freqpeak(find(Ppeak == max(Ppeak)))) (find(Ppeak == max(Ppeak)))];
    
    plotp(2,iii-1)=Ppeak(peak_ax(2));
    pcurve(2,iii-1,:)=ps_curves(peak_ax(2),:);
    clearvars -except iii cohort plotp pcurve
        
    
end

for iii=2:4
load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_spiral/P0',num2str(cohort(iii)),'_PLS_S.mat'))

[d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
data=d.data_raw;  tremor_ds=d.data_ds;
addon=92; addon_end=35;

new = find(data(2,:) > 4);
difp = find((diff(new)) > 100000); % are you trying to threshold at 9.6 seconds?
ep_1 = [new(difp) new(end)];
sp_1 = [new(1) new(difp+1)];

dum=find((ep_1-sp_1)<(60000./samplerate)*samplerateold);
if ~isempty(dum)
    sp_1(1,dum)=NaN;
    ep_1(1,dum)=NaN;
end
sp=sp_1(~isnan(sp_1));
ep=ep_1(~isnan(ep_1));

if iii==4
    ep(2)=ep(1);
    ep(1)=953991;
    sp(2)=1022995;
end
    

time=1:length(data(1,:));
plot(time,data(4,:))
hold on
plot(time(sp),data(4,sp),'r.')
plot(time(ep),data(4,ep),'k.')

start= floor((sp./samplerateold)*samplerate)+addon;
ending = floor((ep./samplerateold)*samplerate)+addon+addon_end;

   handup = [];
   for i = 1:length(start)
       handup = [handup start(i):ending(i)]; %#ok<*AGROW>
       
   end
   clear i
   handup = sort(handup,'ascend');
   
   for aa = 1:3
       [Pxx,F] = pwelch(tremor_ds(aa,handup), samplerate, [], round(samplerate), samplerate);
       frange = F(3:10);
      Pxxrange = Pxx(3:10);
        Freqpeak(aa,:) = frange(find(Pxxrange == max(Pxxrange)));
        Ppeak(aa,:) = max(Pxxrange);
        ps_curves(aa,:) = Pxx;
    end
    
    peak_ax = [(Freqpeak(find(Ppeak == max(Ppeak)))) (find(Ppeak == max(Ppeak)))];
    
        plotp(3,iii-1)=Ppeak(peak_ax(2));
    pcurve(3,iii-1,:)=ps_curves(peak_ax(2),:);
    clearvars -except iii cohort plotp pcurve 
    
end

for iii=2:4
load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_hf/P0',num2str(cohort(iii)),'_HFS_S.mat'))

    data_raw=SmrData.WvData;
    samplerateold=SmrData.SR;
    ts=timeseries(data_raw,0:(1/samplerateold):((size(data_raw,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data_raw,2)-1)/samplerateold),'linear');
    ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    tt=0:1/samplerate:(size(ds_data,2)-1)/samplerate;
    data_t=ds_data([3 6 7],:);
    
  
    for aa = 1:3
        [Pxx,F] = pwelch(data_t(aa,:), samplerate, [], round(samplerate), samplerate);
        frange = F(3:10);
        Pxxrange = Pxx(3:10);
        Freqpeak(aa,:) = frange(find(Pxxrange == max(Pxxrange)));
        Ppeak(aa,:) = max(Pxxrange);
        ps_curves(aa,:) = Pxx;
    end
    
    peak_ax = [(Freqpeak(find(Ppeak == max(Ppeak)))) (find(Ppeak == max(Ppeak)))];
    
    plotp(4,iii-1)=Ppeak(peak_ax(2));
    pcurve(4,iii-1,:)=ps_curves(peak_ax(2),:);
    clearvars -except iii cohort plotp pcurve F
end

for i=1:3
p1=figure(1)
subplot(1,3,i)
set(p1,'color','w')
plot(F, squeeze(pcurve(1,i,:)),'Color',[0.5 0.5 0.5],'LineWidth',3)
hold on
plot(F, squeeze(pcurve(2,i,:)),'Color',[0 0.5 0.5],'LineWidth',3)
plot(F, squeeze(pcurve(3,i,:)),'Color',[0.8 0.5 0.5],'LineWidth',3)
plot(F, squeeze(pcurve(4,i,:)),'Color',[0.8 0.5 0.8],'LineWidth',3)
xlim([2 10])
xticks([2:2:10])
ylim([0 3.5e-4])
box('off')
xlabel('Frequency (Hz)')
ylabel('Power(\muV^2)')
legend({'NS','RS','PLS','HFS'})
legend('boxoff')
end

for i=2:4
change(i-1,:)=round(((plotp(i,:)-plotp(1,:))./plotp(1,:)*100),0);
end


