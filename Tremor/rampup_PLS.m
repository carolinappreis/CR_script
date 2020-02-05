clear all
iiii=[2 3 4 5 8 10 11 13 16 17 18 19 20 21 22 23]; %17 is the last pateint; we have 17 with one pulse and 18 with 5 pulses at the same phase; 19:21 are the second visit of pateint 17 stimulation at 3 different phases with 1 pulse ; 22 and 23 are 2nd visit of pt 17 at 2 different phases with 5 pulses; iiii=[2 3 4 5 8 10 11 13 16 17];
%   iiii=[2 3 13 17];
f=1;
in2=1;
all=[];
for numb= 1:length(iiii)
    close all
    load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\PLS_combined\P0',num2str(iiii(numb)),'_PLSc.mat'))
    %              load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/PLS/p0',num2str(iiii(numb)),'_PLS.mat'))
    
    if in2==1
        in=3;
    elseif in2==2 % other axis 1
        in=5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    elseif in2==3 % other axis 2
        in=6; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    end
    data=SmrData.WvData;
    samplerateold=SmrData.SR;
    tremor=(data(in,:));
    addon=92; addon_end=35;
    
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
    
    st=floor((indexes4./samplerateold)*samplerate)+addon;
    en=floor((indexes3./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
    
    for nn=1:length(st)
        if nn>1
            if (numel(st(nn):en(nn))>30000) && (numel(en(nn-1):st(nn))>25000)
                start(nn)=st(nn);
                ending(nn)=en(nn);
            else
                start(nn)=NaN;
                ending(nn)=NaN;
            end
        else
            start(nn)=st(nn);
            ending(nn)=en(nn);
        end
    end
    
    start=start(~isnan(start));
    ending=ending(~isnan(ending));
    
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
    tf_3=filtfilt(b,a,tre_3')*10*9.81/0.5;
    % tremor_or=zscore(tremor_or);
    dummy=(hilbert(tf_3))';
    envelope=abs(dummy);
    z_env=(abs(hilbert(zscore(tf_3))))';
    phase=angle(dummy);
    frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
    
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\before_PLS.mat')
    %                   load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/before_PLS.mat')
    
    
    segmentb=round((dc_s{numb,:})*samplerate,1);
    
    [d,e]=butter(2,[0.5/(0.5*samplerate) ],'low'); %15
    C=(filtfilt(d,e,tre_3'));
    f3=figure(3)
    plot(zscore(tf_3(:,1)))
    hold on
    plot(zscore(C(:,1)))
    for i=1:size(segmentb,2)
        xline(segmentb(i),'r')
        %                 xline(segmente(i),'k')
        xline(start(i),'k')
    end
    box('off')
    close all
    
    % %      numb_p=[];
    % %     axx=1;
    % %     seg=10000;
    % %     for i=1:length(start)
    % %         if   start(i)>seg && (start(i)+seg-1)<=size(z_env,2)
    % %             %                     figure(1)
    % %             all=[ all  smooth(z_env(axx,start(i)-seg:start(i)+seg-1))];
    % %             %                      plot(z_env(axx,start(i)-seg:start(i)+seg-1))
    % %             %                      hold on
    % %             numb_p=[numb_p  smooth(z_env(axx,start(i)-seg:start(i)+seg-1))];
    % %         end
    % %     end
    % %
    % %     if size(numb_p,2)>1
    % %     y=median(numb_p');
    % %     else
    % %     y=numb_p;
    % %     end
    % %     x=1:length(y);
    % %     initial_params=[];
    % %     [param]=sigm_fit(x,y,initial_params)        % automatic initial_params  "min", "max", "x50" and "slope"
    % %     clear x y
    % %     %     close all
    % %
    % %     id_param_pls(numb,:)=param;
    % %
    clearvars -except iiii numb id_param_pls all in2
    
end
% m50=median(id_param_pls(:,3));