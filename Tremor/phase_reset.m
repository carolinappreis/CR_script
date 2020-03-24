clear all
iiii=[1 2 3 4 5 8 10 11 12 13 16 17 19 20];
% iiii=[2 5 8]; %% significant one
iiii=[ 2 3 4 5 8 10 11 13 16 17];

for numb=1:length(iiii)
    clearvars -except   iiii numb  tt_raw seg_raw tt_raw2
    %     load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iiii(numb)),'_RS.mat'))
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iiii(numb)),'_RS.mat'))
    
    data=SmrData.WvData;
    samplerateold=SmrData.SR;
    ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    tt=0:1/samplerate:(size(ds_data,2)-1)/samplerate;
    tre_3=ds_data([3 5 6],:);
    
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
    
    start=floor((indexes4./samplerateold)*samplerate);
    ending=floor((indexes3./samplerateold)*samplerate);%floor(5*samplerate);
    
    %%% when patient's hand is up
    handup=[];
    for i=1:length(start)
        handup=[handup start(i):ending(i)]; %#ok<*AGROW>
    end
    clear i
    handup=sort(handup,'ascend');
    
    
    %%% tremor characteristics
    [Pxx,F]=pwelch(tre_3(1,handup),samplerate,[],samplerate,samplerate);
    
    frange=F(3:10);
    Pxxrange=Pxx(3:10);
    
    Fpeak=frange(find(Pxxrange==max(Pxxrange))); %#ok<*FNDSB>
    
    if (Fpeak-2)>=1
        [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    else
        [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    end
    
    tf_3=filtfilt(b,a,tre_3')*10*9.81/0.5;
    % tremor_or=zscore(tremor_or);
    dummy=(hilbert(tf_3))';
    envelope=abs(dummy);
    zenv=(abs(hilbert(zscore(tf_3))))';
    phase=angle(dummy);
    frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
    % figure()
    %     bar(sum(envelope'))
    %     box('off')
    
    close all
    
    
    new=find(data(2,:)>4);
    difp=find((diff(new))>100000); % are you trying to threshold at 9.6 seconds?
    ep_1=[new(difp) new(end)];
    sp_1=[new(1) new(difp+1)];
    
    %%% input start all trial
    start_t=1;
    sp=sp_1(1,start_t:end);
    ep=ep_1(1,start_t:end);
    
    for ik=1:length(sp) %%find double start and end points in a stimulation run
        
        s=(find(([indexes4>=sp(ik)]+[indexes4<=ep(ik)])==2));
        e=(find(([indexes3>=sp(ik)]+[indexes3<=ep(ik)])==2));
        tks=(find(diff(xx(s))==0))+1;
        tke=(find(diff(xx(e))==0));
        
        indexes4(s(tks))=NaN;
        indexes3(e(tke))=NaN;
        xx(e(tke))=NaN;
        
    end
    
    indexes4=indexes4(~isnan(indexes4));
    indexes3=indexes3(~isnan(indexes3));
    xx=xx(~isnan(xx));
    
    %%%% find runs with trigering issues (too few, too many pulses)
    th1=(Fpeak*5)./2;
    th2=(Fpeak*5)+round((Fpeak*5)./2);
    %
    %   th1=(Fpeak*5*5)./2;
    %     th2=(Fpeak*5*5)+round((Fpeak*5*5)./2);
    
    for it=1:length(indexes4)
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
    %%%%%%%%%%%%%%%
    
    start1=[];
    ending1=[];
    xx1=[];
    for il=1:length(sp)
        start1=[start1 indexes4(find(([indexes4>=sp(il)]+[indexes4<=ep(il)])==2))]; % intersect([1 2 3],[3 4 5])
        ending1=[ending1 indexes3(find(([indexes3>=sp(il)]+[indexes3<=ep(il)])==2))];
        xx1=[xx1 xx(find(([indexes4>=sp(il)]+[indexes4<=ep(il)])==2))];
    end
    
    clear start ending
    start{1,1}=floor((start1./samplerateold)*samplerate);
    ending{1,1}=floor((ending1./samplerateold)*samplerate);%floor(5*samplerate);
    clear xx
    xx{1,1}=xx1;
    
    %-----
    hh=1;
    axx=1;
    
    ti=500;
    seg_raw=NaN(length(start{hh,1}),ti*2);
    seg_pha=NaN(length(start{hh,1}),ti*2);
    
    for i=1:length(start{hh,1})
        if (~isnan(start{hh,1}(i)))
            seg_raw(i,1:(ti*2))=ds_data(3,(start{hh,1}(i)-(ti-1)):(start{hh,1}(i)+ti));
            xx{hh,1}(i)= xx{hh,1}(i);
        else
            seg_raw(i,1:(ti*2))=NaN;
            xx{hh,1}(i)= NaN;
        end
    end
    
    %         %% criteria for outliers
    %
    %         idx_outl=find(tremor_or2>(nanmean(tremor_or2)+(2*(nanstd(tremor_or2))))|tremor_or2<(nanmean(tremor_or2)-(2*(nanstd(tremor_or2)))));
    %         tremor_or2(idx_outl,1)=NaN;
    %         tremor_or3(idx_outl,1)=NaN;
    %         xx(1,idx_outl)=NaN;
    
    
    
    yy=xx{hh,1}(:);
    
    
    for i=1:12
        tt_raw2(numb,i,:)=mean(seg_raw(find(yy==i),:));
        
        seg_raw1=NaN(20,2*ti);
        dum=find(yy==i);
        seg_raw1(1:length(dum),:)=seg_raw(dum,:);
        tt_raw{numb,i}=seg_raw1;
    end
    
    
    
    clearvars -except iiii numb  tt_raw seg_raw tt_raw2
    
end



clear all ; close all
load ('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/1sec_trig.mat')
for i=1:10
    figure()
    %      cr=squeeze(tt_raw2(i,ii,(size(tt_raw2,3)/2-200):(size(tt_raw2,3)/2+200-1)));
    
    for ii=1:12
        cr=squeeze(tt_raw2(i,ii,:));
        subplot(12,1,ii)
        plot(cr)
        hold on
        xline(size(tt_raw2,3)/2)
    end
end








clear all
iiii=[ 2 3 4 5 8 10 11 13 16 17];

for numb=1:length(iiii)
    clearvars -except   iiii numb  tt_raw seg_raw tt_raw2
    %     load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iiii(numb)),'_RS.mat'))
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iiii(numb)),'_RS.mat'))
    
    data=SmrData.WvData;
    
    samplerateold=SmrData.SR;
    ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    tt=0:1/samplerate:(size(ds_data,2)-1)/samplerate;
    tre_3=ds_data([3 5 6],:);
    zdata=zscore(data(3,:));
    
    %%% determine stimulation time points
    index=[];
    for i=2:size(data,2)-1
        if data(2,i-1)<2.5 && data(2,i)>2.5
            index=[index i];
        end
    end
    
    dd2=round(data(4,:)*100)./100;
    for i=1:length(index)
        xx(i)=round(dd2(index(i))./0.1); %#ok<*SAGROW>
    end
    clear i
    
    
    
    f=floor(((250/2)./samplerate)*samplerateold);
    seg_raw=NaN(length(index),f*2);
    seg_pha=NaN(length(index),f*2);
    
    for i=1:length(index)
        if (~isnan(index(i)))
            seg_raw(i,1:(f*2))=zdata(1,(index(i)-(f-1)):(index(i)+f));
        else
            seg_raw(i,1:(f*2))=NaN;
        end
    end
    
    
    yy=xx;
    
    
    for i=1:12
        tt_raw2(numb,i,:)=mean(seg_raw(find(yy==i),:));
        
        %         seg_raw1=NaN(20,2*ti);
        %         dum=find(yy==i);
        %         seg_raw1(1:length(dum),:)=seg_raw(dum,:);
        %         tt_raw{numb,i}=seg_raw1;
    end
    
    
end

dn=squeeze(nanmean(tt_raw2,2));
close all
f1=figure(1);
r=1;
for m=1:10
    subplot(1,10,m)
    plot(dn(m,:)-mean(dn(m,:)));
    hold on
    xline(size(tt_raw2,3)/2,'--')
    box('off')
    xlim([0 size(tt_raw2,3)])
     ylim ([-0.25 0.4])
     xticks([size(tt_raw2,3)/2])
     xticklabels(['stim'])
 yticks([])
    r=r+1;
    
    %     for i=1:12
    %         subplot(10,12,r)
    %         plot(squeeze(tt_raw2(m,i,:)))
    %         hold on
    %         xline(size(tt_raw2,3)/2,'--')
    %         box('off')
    %         xlim([0 size(tt_raw2,3)])
    %         r=r+1;
    % %         xticks([])
    % %         yticks([])
    %     end
    
end



f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 30, 10];
set(f1,'color','w');

