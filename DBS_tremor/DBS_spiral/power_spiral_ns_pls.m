clear; close
cohort = [1 3 4 6];
iii=4
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


 for iii =  2:length(cohort)
        close all
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
    

% time=1:length(data(1,:));
% plot(time,data(4,:))
% hold on
% plot(time(sp),data(4,sp),'r.')
% plot(time(ep),data(4,ep),'k.')

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
    
    plotp(2,iii-1)=Ppeak(peak_ax(2));
    pcurve(2,iii-1,:)=ps_curves(peak_ax(2),:);
    clearvars -except iii cohort plotp pcurve F
 end
  
 
i=3
p1=figure(1)
set(p1,'color','w')
plot(F, squeeze(pcurve(1,i,:)),'Color',[0.5 0.5 0.5],'LineWidth',3)
hold on
plot(F, squeeze(pcurve(2,i,:)),'Color',[0 0.5 0.5],'LineWidth',3)
xlim([2.1 10])
xticks([3:2:10])
box('off')
xlabel('Frequency (Hz)')
ylabel('Power(\muV^2)')
legend({'NO STIM','STIM'})
legend('boxoff')

