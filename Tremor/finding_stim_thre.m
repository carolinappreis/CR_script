clear all
iiii=[1 2 3 4 5 8 10 11 12 13 16 17 19 20];
% iiii=[2 5 8]; %% significant one
iiii=[ 2 3 4 5 8 10 11 13 16 17];

for numb=1:length(iiii)
    
    clearvars -except iiii numb m_all m_notrig m_lowtrig
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
    
    
    
    new=find(data(2,:)>4);
    difp=find((diff(new))>100000); % are you trying to threshold at 9.6 seconds?
    ep_1=[new(difp) new(end)];
    sp_1=[new(1) new(difp+1)];
    
    %%% input start all trial
    start_t=1;
    sp=sp_1(1,start_t:end);
    ep=ep_1(1,start_t:end);
    
    ss=floor((indexes4./samplerateold)*samplerate);
    ee=floor((indexes3./samplerateold)*samplerate);
    
    pp_all=[];
    for jk=1:length(ss)
        pp_all=[pp_all z_env(1,ss(jk):ee(jk))];
    end
    
    dd1=indexes3;
    dd2=indexes4;
    
    low_tre=[];
    for ik=1:length(sp) %%find double start and end points in a stimulation run
        
        s=(find(([indexes4>=sp(ik)]+[indexes4<=ep(ik)])==2));
        e=(find(([indexes3>=sp(ik)]+[indexes3<=ep(ik)])==2));
        tks=(find(diff(xx(s))==0))+1;
        tke=(find(diff(xx(e))==0));
        
        indexes4(s(tks))=NaN;
        indexes3(e(tke))=NaN;
        xx(e(tke))=NaN;
        
    end
    
    d1=dd1(find(isnan(indexes3)));
    d2=dd2(find(isnan(indexes4)));
    d11=floor((d1./samplerateold)*samplerate);
    d22=floor((d2./samplerateold)*samplerate);
    
    low_tre=[];
    for i =1:length(d11)
        low_tre=[low_tre z_env(1,d11:d22)];
    end
    
    %%%% find runs with trigering issues (too few, too many pulses)
    th1=(Fpeak*5)./2;
    th2=(Fpeak*5)+round((Fpeak*5)./2);
    
    ind_out=[];
    for it=1:length(indexes4)
        
        strial=index(find(index==indexes4(it)):find(index==indexes3(it)));
        nontrig=find(diff(strial)>(((Fpeak*5)/2)*samplerate)); %%% lack of triggering for more than half of stim time
        
        if numel(index(find(index==indexes4(it)):find(index==indexes3(it))))>=th1 && numel(index(find(index==indexes4(it)):find(index==indexes3(it))))<=th2
            indexes4(it)=indexes4(it);
            indexes3(it)=indexes3(it);
            xx(it)=xx(it);
        else
            indexes4(it)=NaN;
            indexes3(it)=NaN;
            xx(it)=NaN;
            ind_out=[ind_out it];
        end
        
        
    end
    %%%%%%%%%%%%%%%
    
    
    %%% all rejected 5 sec
    pp=[];
    for jk=1:length(ind_out)
        pp=[pp z_env(1,ss(ind_out(jk)):ee(ind_out(jk)))];
    end
    
%     figure(1)
%     subplot(5,2,numb)
%     histogram(pp_all,'FaceColor','w')
%     hold on
%     histogram(pp,'FaceColor','r')
    m_notrig(numb,:)=median(low_tre);
    m_lowtrig(numb,:)=nanmedian(pp);
    m_all(numb,:)=nanmedian(pp_all);
end


save('median_z_env.mat')



