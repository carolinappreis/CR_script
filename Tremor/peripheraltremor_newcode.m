clear all
% iiii=[1 2 3 4 5 8 10 11 12 13 16];
iiii=[2 5 8];

for numb=1:length(iiii);
    clearvars -except iiii numb
    % load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iiii(numb)),'_RS.mat'))
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iiii(numb)),'_RS.mat'))
    
    in2=1; % analysing the "main tremor axis"
    
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
    
    time=0:1/samplerateold:(size(data,2)-1)/samplerateold;
    
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
    
    
    
    %%
    
    new=find(data(2,:)>4);
    difp=find((diff(new))>100000); % are you trying to threshold at 9.6 seconds?
    ep_1=[new(difp) new(end)];
    sp_1=[new(1) new(difp+1)];
    
    %%% input start all trial
    start_t=1;
    sp=sp_1(1,start_t:end);
    ep=ep_1(1,start_t:end);
    
    
    
    
    %plot(time,data(4,:))
    %hold on
    %plot(time(sp),data(4,sp),'r.')
    %plot(time(ep),data(4,ep),'k.')
    
    
    for ik=1:length(sp) %%find double start and end points in a stimulation run
        
        s=(find(([indexes4>=sp(ik)]+[indexes4<=ep(ik)])==2));
        e=(find(([indexes3>=sp(ik)]+[indexes3<=ep(ik)])==2));
        tks=(find(diff(xx(s))==0))+1;
        tke=(find(diff(xx(e))==0));
        
        indexes4(s(tks))=NaN;
        indexes3(e(tke))=NaN;
        xx(e(tke))=NaN;
        
    end
    
    %%%% find runs with trigering issues (too few, too many pulses)
    th1=(Fpeak*5)./2;
    th2=(Fpeak*5)+round((Fpeak*5)./5);
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
    %%%%%%%%%%%%%%%
    indexes4=indexes4(~isnan(indexes4));
    indexes3=indexes3(~isnan(indexes3));
    xx=xx(~isnan(xx));
    
    
    start1=[];
    ending1=[];
    xx1=[];
    for il=1:length(sp)
        start1=[start1 indexes4(find(([indexes4>=sp(il)]+[indexes4<=ep(il)])==2))]; % intersect([1 2 3],[3 4 5])
        ending1=[ending1 indexes3(find(([indexes3>=sp(il)]+[indexes3<=ep(il)])==2))];
        xx1=[xx1 xx(find(([indexes4>=sp(il)]+[indexes4<=ep(il)])==2))];
    end
    
    
    %     figure()
    %     plot(time,data(4,:))
    %     hold on
    %     plot(time(index),data(4,index),'r.')
    %     plot(time(start1),data(4,start1),'ko')
    %     plot(time(ending1),data(4,ending1),'bo')
    
    
    clear start ending
    start{1,1}=floor((start1./samplerateold)*samplerate)+addon;
    ending{1,1}=floor((ending1./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
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
            else
                tremor_or2(axx,i,1)=NaN;
                tremor_pc(i,1)=NaN;
                xx{hh,1}(i)= NaN;
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
            tt_pc(1:sum(yy==i),i)=tremor_pc(1,find(yy==i));
            tt(tt==0)=NaN;
            tt_pc(tt_pc==0)=NaN;
        end
        tt1{hh,axx}=tt;
    end
    
    pca_ax{hh,1}=ma;
    PSI_ax{hh,1}=mean(PSI,1);
    clear yy ma;
    
    
    
    
    for i=1
        a=hist(pca_ax{i,1},3);
        axmax=find(a==max(a));
        f2=figure()
        subplot(1,5,1)
        bar(0:30:330,100.*nanmedian(tt1{i,axmax}))
        hold on
        plot(0:30:330,100.*tt1{i,axmax},'.')
        box('off')
        subplot(1,5,2)
        dum=repmat(nanmedian(tt1{i,axmax}),1,3);
        for iii=size(tt1{1,1},2)+1:size(tt1{1,1},2)*2
            dum1(1,iii-12)=sum(dum(1,(iii-1:iii+1)))./length(dum(1,(iii-1:iii+1)));
        end
        bar(0:30:330,100.*dum1)
        box('off')
        
        subplot(1,5,3)
        cr=repmat(nanmedian(tt_pc),1,3);
        for iii=size(tt1{1,1},2)+1:size(tt1{1,1},2)*2
            cr1(1,iii-12)=sum(cr(1,(iii-1:iii+1)))./length(cr(1,(iii-1:iii+1)));
        end
        bar(0:30:330,100.*cr1)
        box('off')
        
        subplot(1,5,4)
        bar(1:3,a)
        names = {'CED2'; 'CED5';'CED6'};
        set(gca,'xtick',[1:3],'xticklabel',names)
        box('off')
        subplot(1,5,5)
        bar([1 2],[PSI_ax{i}(1,1);PSI_ax{i}(1,2)]);
        names = {'PSI 2 5'; 'PSI 2 6'};
        set(gca,'xtick',[1:2],'xticklabel',names)
        f2.Units = 'centimeters';
        f2.OuterPosition= [10, 10, 60, 10];
        set(f2,'color','w');
        
    end
    
    
    %     clearvars -except PSI_ax pca_ax tt1 tt  iiii numb ampall ph_stim LS
    
end


