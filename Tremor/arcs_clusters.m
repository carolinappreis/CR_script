clear all
iiii=[1 2 3 4 5 8 10 11 12 13 16 17 19 20];
% iiii=[2 5 8]; %% significant one
iiii=[ 2 3 4 5 8 10 11 13 16 17];
load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/clusters_BA'))

p=C_NS;
for numb=1:length(iiii)
    u{numb,1}=unique(p{numb,1});
    dum=hist(p{numb,1});
    h{numb,1}=nonzeros(dum)';
end
    
 for numb=1:length(iiii)   
    
    clearvars -except iiii numb tt1 LS seg_env seg_filt
%     load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iiii(numb)),'_RS.mat'))
        load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iiii(numb)),'_RS.mat'))
    
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
    for pt=1:3;
    frequency=(smooth((1000/(2*pi))*diff(unwrap(phase(pt,:))),500))';
    end
    z_env=[abs(hilbert(zscore(tremorxf)));abs(hilbert(zscore(tremoryf)));abs(hilbert(zscore(tremorzf)))];
    z_filt=[zscore(tremorxf);zscore(tremoryf);zscore(tremorzf)];
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
    for j=1:length(start{hh,1})
        if (~isnan(start{hh,1}(j)))
            x=[sum(envelope(1,start{hh,1}(j):ending{hh,1}(j)));sum(envelope(2,start{hh,1}(j):ending{hh,1}(j)));sum(envelope(3,start{hh,1}(j):ending{hh,1}(j)))];
            y=[tremorxf(start{hh,1}(j):ending{hh,1}(j));tremoryf(start{hh,1}(j):ending{hh,1}(j));tremorzf(start{hh,1}(j):ending{hh,1}(j))];
            [pc,score,latent,tsquare] = pca(y');
            yyy(j,1:3)=pc(1:3,1);
            PSI(j,1:2)=[abs(sum(exp(sqrt(-1).*(phase(1,start{hh,1}(j):ending{hh,1}(j))-phase(2,start{hh,1}(j):ending{hh,1}(j)))))./length(phase(1,start{hh,1}(j):ending{hh,1}(j))));abs(sum(exp(sqrt(-1).*(phase(1,start{hh,1}(j):ending{hh,1}(j))-phase(3,start{hh,1}(j):ending{hh,1}(j)))))./length(phase(1,start{hh,1}(j):ending{hh,1}(j))))];
            ma(j)=find(x==max(x));
            ma2(j)=(find(abs(yyy(j,1:3))==max(abs(yyy(j,1:3)))));
            clear x y
            %             x=[tremorxf(start{hh,1}(j):ending{hh,1}(j));tremoryf(start{hh,1}(j):ending{hh,1}(j));tremorzf(start{hh,1}(j):ending{hh,1}(j))];
            %             [pc,score,latent,tsquare] = pca(x');
            %             xxx(j,1:3)=pc(1:3,1);
            %             ma(j)=(find(abs(xxx(j,1:3))==max(abs(xxx(j,1:3)))));
        end
    end
    
    tremor_or2=NaN(length(start{hh,1}),1);
    tremor_pc=NaN(length(start{hh,1}),1);
    
    for axx=1:3
        for i=1:length(start{hh,1})
            if (~isnan(start{hh,1}(i)))
                tremor_or2(axx,i,1)=(mean(envelope(axx,ending{hh,1}(i)-1000:ending{hh,1}(i)))-mean(envelope(axx,start{hh,1}(i)-1000:start{hh,1}(i))))/mean(envelope(axx,start{hh,1}(i)-1000:start{hh,1}(i)));
                tremor_pc(1,i)=(mean(envelope(ma(i),ending{hh,1}(i)-1000:ending{hh,1}(i)))-mean(envelope(ma(i),start{hh,1}(i)-1000:start{hh,1}(i))))/mean(envelope(ma(i),start{hh,1}(i)-1000:start{hh,1}(i)));
                xx{hh,1}(i)= xx{hh,1}(i);
                env_pha(i,1:5000)=z_env(1,start{hh,1}(i):start{hh,1}(i)+5000-1);
                filt_pha(i,1:5000)=z_filt(1,start{hh,1}(i):start{hh,1}(i)+5000-1);
                
            else
                tremor_or2(axx,i,1)=NaN;
                tremor_pc(i,1)=NaN;
                xx{hh,1}(i)= NaN;
                env_pha(i,1:5000)=NaN;
                filt_pha(i,1:5000)=NaN;
            end
        end
        
        %         %% criteria for outliers
        %
        %         idx_outl=find(tremor_or2>(nanmean(tremor_or2)+(2*(nanstd(tremor_or2))))|tremor_or2<(nanmean(tremor_or2)-(2*(nanstd(tremor_or2)))));
        %         tremor_or2(idx_outl,1)=NaN;
        %         tremor_or3(idx_outl,1)=NaN;
        %         xx(1,idx_outl)=NaN;
        
        tt=NaN(20,12);
        yy=xx{hh,1}(:);
        tt_pc=NaN(20,12);
        
        for i=1:12
            tt(1:sum(yy==i),i)=tremor_or2(axx,find(yy==i));
            stim_pha(1,i)=numel(tremor_or2(axx,find(yy==i)));
            tt_pc(1:sum(yy==i),i)=tremor_pc(1,find(yy==i));
            tt(tt==0)=NaN;
            tt_pc(tt_pc==0)=NaN;
            seg_env{numb,i}=env_pha(find(yy==i),:);
            seg_filt{numb,i}=filt_pha(find(yy==i),:);
        end
        tt1{numb,axx}=tt;
        
        for rr=1:100000
            LS(numb,axx,rr)=nanmedian(tt(randi(length(start1),1,10)));
        end
        
        
        %     close all
       
    end
    stim_pha 
    clearvars -except PSI_ax pca_ax tt1  iiii numb ampall ph_stim LS seg_env seg_filt
    
end
clearvars -except iiii tt1 seg_env seg_filt
cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')




