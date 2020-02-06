
clear all

% iiii=[2 3 4 5 8 10 11 13 16 17 18 19 20 21 22 23]; %17 is the last pateint; we have 17 with one pulse and 18 with 5 pulses at the same phase; 19:21 are the second visit of pateint 17 stimulation at 3 different phases with 1 pulse ; 22 and 23 are 2nd visit of pt 17 at 2 different phases with 5 pulses; iiii=[2 3 4 5 8 10 11 13 16 17];
 iiii=[2 3 4 5 8 10 11 13 16 17];
% iiii=[2 8 17];
f=1;

all=[];
for numb= 1:length(iiii)
    clearvars -except iiii all f numb id_param_pls param_pls_1min good bad
    %     close all
    load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\PLS_combined\P0',num2str(iiii(numb)),'_PLSc.mat'))
    %              load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/PLS/p0',num2str(iiii(numb)),'_PLS.mat'))
    
    data=SmrData.WvData;
    samplerateold=SmrData.SR;
    addon=92; addon_end=35;
    
    ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    to=0:1/samplerateold:(size(data,2)-1)/samplerateold;
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
    
    start=floor((indexes4./samplerateold)*samplerate)+addon;
    ending=floor((indexes3./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
    
    handup=[];
    for i=1:length(start)
        handup=[handup start(i):ending(i)]; %#ok<*AGROW>
    end
    clear i
    handup=sort(handup,'ascend');
    
    for aa=1:3
        [Pxx,F]=pwelch(tre_3(aa,handup),samplerate,[],samplerate,samplerate);
        frange=F(3:10);
        Pxxrange=Pxx(3:10);
        Freqpeak(aa,:)=frange(find(Pxxrange==max(Pxxrange)));
        Ppeak(aa,:)=max(Pxxrange);
        ps_curves(aa,:)=Pxx;
    end
    peak_ax=[(Freqpeak(find(Ppeak==max(Ppeak)))) (find(Ppeak==max(Ppeak)))];
    Fpeak=peak_ax(1);
    
    %     figure(2)
    %     plot(F(3:50),ps_curves(:,3:50)','LineWidth',2)
    %     legend({'CED2','CED4','CED5'})
    %     legend('boxoff')
    %     box('off')
    %
    if (Fpeak-2)>=1
        [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    else
        [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    end
    for i=1:3
        tf_3(i,:)=(filtfilt(b,a,tre_3(i,:)))*10*9.81/0.5;
        envelope(i,:)=abs(hilbert(tf_3(i,:)));
        z_env(i,:)=abs(hilbert(zscore(tf_3(i,:))));
    end
    
    
    new=find(data(2,:)>4);
    difp=find((diff(new))>100000); % are you trying to threshold at 9.6 seconds?
    ep_1=[new(difp) new(end)];
    sp_1=[new(1) new(difp+1)];
    
    %%% input start all trial
    start_t=1;
    sp=sp_1(1,start_t:end);
    ep=ep_1(1,start_t:end);
    
    
    %     plot(to,data(1,:))
    %     hold on
    %     plot(to,data(2,:))
    %     plot(to(indexes4),data(1,indexes4),'k.','MarkerSize',10)
    %     plot(to(indexes3),data(1,indexes3),'g.','MarkerSize',10)
    %
    for ik=1:length(sp) %%find double start and end points in a stimulation run
        
        s=(find(([indexes4>=sp(ik)]+[indexes4<=ep(ik)])==2));
        e=(find(([indexes3>=sp(ik)]+[indexes3<=ep(ik)])==2));
        tks=(find(diff(xx(s))==0))+1;
        tke=(find(diff(xx(e))==0));
        
        indexes4(s(tks))=NaN;
        indexes3(e(tke))=NaN;
        
        
    end
    
    indexes4=indexes4(~isnan(indexes4));
    indexes3=indexes3(~isnan(indexes3));
    
    
    %%% find runs with trigering issues (too few, too many pulses)
    th1=(Fpeak*90)./2;
    th2=(Fpeak*90)+round((Fpeak*90)./2);
    
    for it=1:length(indexes4)
        if numel(index(find(index==indexes4(it)):find(index==indexes3(it))))>=th1 && numel(index(find(index==indexes4(it)):find(index==indexes3(it))))<=th2
            indexes4(it)=indexes4(it);
            indexes3(it)=indexes3(it);
        else
            indexes4(it)=NaN;
            indexes3(it)=NaN;
        end
    end
    
    indexes4=indexes4(~isnan(indexes4));
    indexes3=indexes3(~isnan(indexes3));
    
    clear start ending
    start=floor((indexes4./samplerateold)*samplerate);
    ending=floor((indexes3./samplerateold)*samplerate);%floor(5*samplerate);
    
    if numb==3
        dum=([1:4 6:length(start)]);
        start=start(dum);
        ending=ending(dum);
        
    elseif numb==7
        dum=([1:3 5:length(start)]);
        start=start(dum);
        ending=ending(dum);
    end
    
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\before_PLS.mat')
    % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/before_PLS.mat')
    segmentb=round((dc_s{numb,:})*samplerate,1);
    
    %     [d,e]=butter(2,[0.5/(0.5*samplerate) ],'low'); %15
%         C=(filtfilt(d,e,tre_3(1,:)));
%         figure
%         plot(zscore(C))
%         hold on
%         plot(zscore(tf_3(1,:)))
%         hold on
%         for i=1:size(start,2)
%             xline(segmentb(i),'k')
%             xline(start(i),'r')
%         end
    %     box('off')
    %
    %         figure(2)
    %         plot(data(1,:))
    %         hold on
    %         plot(data(2,:))
    %     %     close all
    
    
    
    
    axx=1;
    seg=5000;
    for i=1:length(start)
%         figure
%         plot(tf_3(1,segmentb(i):segmentb(i)+seg2-1))
%         hold on 
%         plot(envelope(1,segmentb(i):segmentb(i)+seg2-1))
%         xline(start(i),'LineWidth',2)
       
        if   segmentb(i)>seg && (segmentb(i)+seg-1)<=size(z_env,2)
            y=smooth(z_env(axx,segmentb(i)-seg:start(i)));

        else
            y=smooth(z_env(axx,segmentb(i):start(i)));

        end
%          all=[ all  smooth(z_env(axx,segmentb(i)-seg:segmentb(i)+seg-1))];
        
        x=1:length(y);
        initial_params=[];
        [param1]=sigm_fit(x,y,initial_params)        % automatic initial_params  "min", "max", "x50" and "slope"
        clear x y
        
        
        y=smooth(z_env(axx,segmentb(i):ending(i)));
        x=1:length(y);
        initial_params=[];
        [param2]=sigm_fit(x,y,initial_params)        % automatic initial_params  "min", "max", "x50" and "slope"
        clear x y
        hold on
        xline(start(i)-segmentb(i))
 
        if param1(3)>0 && param1(3)<start(1) && param2(3)>0 && param2(3)<ending(1) && param2(3)>(start(i)-segmentb(i))
        good{numb,1}(:,i)=i;
        bad{numb,1}(:,i)=NaN;
        else
        bad{numb,1}(:,i)=i;
        good{numb,1}(:,i)=NaN;
        end
    
        close all
    end
    id_param_pls(numb,:)=param1;
    param_pls_1min(numb,:)=param2;
    
    clearvars -except iiii numb id_param_pls all param_pls_1min good bad
    
end
% m50=median(id_param_pls(:,3));