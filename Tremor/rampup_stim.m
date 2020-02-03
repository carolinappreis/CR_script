% 
cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
load('param_ns_s.mat')
[p,h]=ttest(id_param_NS,id_param_stim);
h(3)
median(id_param_NS(:,3))
median(id_param_stim(:,3))


clear all
iiii=[1 2 3 4 5 8 10 11 12 13 16 17 19 20];
% iiii=[2 5 8]; %% significant one
iiii=[ 2 3 4 5 8 10 11 13 16 17];
all=[];
for numb=1:length(iiii)
    clearvars -except iiii numb tt1 LS cr all id_param_stim
    load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iiii(numb)),'_RS.mat'))
    %     load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iiii(numb)),'_RS.mat'))
    
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
    
    %     f1=figure(1)
    %     stackedplot(data')
    %     title(['pt_',num2str(iiii(numb))])
    %
    %%% downsample
    
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremor2(1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    
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
    
    %%% when patient's hand is up
    handup=[];
    for i=1:length(start)
        handup=[handup start(i):ending(i)]; %#ok<*AGROW>
    end
    clear i
    handup=sort(handup,'ascend');
    
    
    %%% tremor characteristics
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
    dummy=hilbert(tremor_or);
    % phase=angle(dummy);
    frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
    
    tremor=(data(3,:));% %score(:,1)';%
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremorx(1:size(ts1.data,3))=ts1.data;
    tremor=(data(5,:));% %score(:,1)';%
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremory(1:size(ts1.data,3))=ts1.data;
    tremor=(data(6,:));% %score(:,1)';%
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremorz(1:size(ts1.data,3))=ts1.data;
    tremorxf=filtfilt(b,a,tremorx);
    tremoryf=filtfilt(b,a,tremory);
    tremorzf=filtfilt(b,a,tremorz);
    envelope=[abs(hilbert(tremorxf));abs(hilbert(tremoryf));abs(hilbert(tremorzf))];
    phase=[angle(hilbert(tremorxf));angle(hilbert(tremoryf));angle(hilbert(tremorzf))];
    z_env=[abs(hilbert(zscore(tremorxf)));abs(hilbert(zscore(tremoryf)));abs(hilbert(zscore(tremorzf)))];
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
    
    
    st=floor((sp./samplerateold)*samplerate);
    en=floor((ep./samplerateold)*samplerate);
    
    
    numb_p=[];
    axx=1;
    seg=10000;
    for i=1:length(st)
        if   st(i)>seg && (st(i)+seg-1)<=size(z_env,2)
            %                     figure(1)
            all=[ all  smooth(z_env(axx,st(i)-seg:st(i)+seg-1))];
            %                      plot(z_env(axx,st(i)-seg:st(i)+seg-1))
            %                      hold on
            numb_p=[numb_p  smooth(z_env(axx,st(i)-seg:st(i)+seg-1))];
        end
    end
    
    
    y=median(numb_p');
    x=1:length(y);
   
    if numb==4
        initial_params=[ NaN NaN length(x)/2 NaN];
    elseif numb==6
        initial_params=[ y(1) y(length(x)) length(x)/2 NaN];
    else
         initial_params=[];
    end
    [param]=sigm_fit(x,y,initial_params)        % automatic initial_params  "min", "max", "x50" and "slope"
    clear x y
    %     close all
    
    id_param_stim(numb,:)=param;
    
    clearvars -except PSI_ax pca_ax tt1  iiii numb ampall ph_stim LS cr all id_param_stim
    
end
 m50=median(id_param_stim(:,3));
% 
% %     clearvars -except iiii tt1 LS
% %     cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
% 
% 
% 
% y=median(all');
% x=1:length(y);
% initial_params=[];
% [param]=sigm_fit(x,y,initial_params)        % automatic initial_params  "min", "max", "x50" and "slope"
% clear x y
% 

