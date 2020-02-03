clear all
iiii=[1 2 3 4 5 8 10 11 12 13 16 17 19 20];
% iiii=[2 5 8]; %% significant one
iiii=[ 2 3 4 5 8 10 11 13 16 17];

for numb=1:length(iiii)
    
    clearvars -except iiii numb m_all m_notrig m_lowtrig th t_below var_th t_var
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
    
    initial_idx4=round(data(4,:)*100)./100;
    for i=1:length(indexes4)
        xx(i)=round(initial_idx4(indexes4(i))./0.1); %#ok<*SAGROW>
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
    env=[abs(hilbert(tremorxf));abs(hilbert(tremoryf));abs(hilbert(tremorzf))];
    phase=[angle(hilbert(tremorxf));angle(hilbert(tremoryf));angle(hilbert(tremorzf))];
    z_env=[abs(hilbert(zscore(tremorxf)));abs(hilbert(zscore(tremoryf)));abs(hilbert(zscore(tremorzf)))];
    % figure()
    %     bar(sum(env'))
    %     box('off')
    
    all_idxs=floor((index./samplerateold)*samplerate);
    trg=zeros(1,length(env));
    trg(all_idxs)=1;
    ti=0:1/samplerate:(length(env)-1)/samplerate;
    
%     figure(1)
%     plot(ti,env(1,:))
%     hold on
%     plot(ti,trg)
%     plot(ti(start),trg(1,start),'r*','MarkerSize',6)
%     plot(ti(ending),trg(1,ending),'k*','MarkerSize',6)
    
    
    new=find(data(2,:)>4);
    difp=find((diff(new))>100000); % are you trying to threshold at 9.6 seconds?
    % 104166 do change to see if we can move to 10 sec instead of 9.6
    ep_1=[new(difp) new(end)];
    sp_1=[new(1) new(difp+1)];
    
    %%% input start all trial
    start_t=1;
    sp=sp_1(1,start_t:end);
    ep=ep_1(1,start_t:end);
    
    pp_all=[];
    for jk=1:length(start)
        %         pp_all=[pp_all env(1,start(jk):ending(jk))];
        pp_all=[pp_all env(1,start(jk):ending(jk))];
    end
    
    
    
    for ik=1:length(sp) %%find double start and end points in a stimulation run
        
        s=(find(([indexes4>=sp(ik)]+[indexes4<=ep(ik)])==2));
        e=(find(([indexes3>=sp(ik)]+[indexes3<=ep(ik)])==2));
        tks=(find(diff(xx(s))==0))+1;
        tke=(find(diff(xx(e))==0));
        
        indexes4(s(tks))=NaN;
        indexes3(e(tke))=NaN;
        xx(e(tke))=NaN;
        
    end
    
    d1=find(isnan(indexes3));
    d2=find(isnan(indexes4));
    d11=ending(d1);
    d22=start(d2);
    
    pt_s=[];pt_e=[];
    low_tre=[];
    for i =1:length(d11)
        low_tre=[low_tre env(1,d11(i):d22(i))];
        %         low_tre=[low_tre env(1,d11(i):d22(i))];
        pt_s=[pt_s d11(i)];
        pt_e=[pt_e d22(i)];
    end
    
    
%         figure(1)
%         plot(ti,trg)
%         hold on
%         plot(ti,env(1,:))
%         plot(ti(start),trg(1,start),'r.','MarkerSize',6)
%         plot(ti(ending),trg(1,ending),'k.','MarkerSize',6)
%         plot(ti(d22),trg(d22),'bd','MarkerSize',5)
%         plot(ti(d11),trg(d1),'gd','MarkerSize',5)
    
    if isempty(pt_e)
        t_below(numb,:)=NaN;
        th(numb,:)=NaN;
        low_tre=NaN;
    else
        x=env(1,:);
        nn=(env(1,round(pt_s+1000/Fpeak)));
        mm=(env(1,round(pt_e+1000/Fpeak)));
        q=[nn mm];
        %         x=env(1,:);
        %         q=[env(pt_s) env(pt_e)];
        prtiles = invprctile(x,q);
        th(numb,:)=max(prtiles);
        t_below(numb,:)=mean(d11-d22);
        
    end
    
    clear d1 d2 d11 d22
    %    check
    %        plot(trg)
    %        hold on
    %        plot(tremorxf,'LineWidth',0.5)
    %        plot(env(1,:),'LineWidth',2)
    %        yline(prctile(env(1,:),th))
    %        xline(pt_s)
    %        xline(pt_e)
    %        xline(pt_s+1000/Fpeak,'r','LineWidth',1.5)
    %        qvalues = prctile(x,prtiles) % check if same values
    
    indexes4=indexes4(~isnan(indexes4));
    indexes3=indexes3(~isnan(indexes3));
    xx=xx(~isnan(xx));
    
    %%%% find runs with trigering issues (too few, too many pulses)
    th1=(Fpeak*5)./2;
    th2=(Fpeak*5)+round((Fpeak*5)./2);
    
    
    for it=1:length(indexes4)
        if numel(index(find(index==indexes4(it)):find(index==indexes3(it))))>=th1 && numel(index(find(index==indexes4(it)):find(index==indexes3(it))))<=th2
            indexes4(it)=indexes4(it);
            indexes3(it)=indexes3(it);
            xx(it)=xx(it);
        else
            indexes4(it)=NaN;
            indexes3(it)=NaN;
            xx(it)=NaN;
            %         nontrig=find(diff(strial)>(((Fpeak*5)/2)*samplerate)); %%% lack of triggering for more than half of stim time
        end
    end
    %%%%%%%%%%%%%%%
    d1=find(isnan(indexes3));
    d2=find(isnan(indexes4));
    d11=ending(d1);
    d22=start(d2);
    
    
    
    if ~isempty(d2)
        
        pp=[];
        for jk=1:length(d11)
            pp=[pp env(1,d22(jk):d11(jk))];
            varb(1,jk)=var(tremorxf(1,d22(jk):d11(jk)));
            t_var(numb,1)=mean(d11-d22);
        end
    else
        pp=NaN;
        varb=NaN;
        t_var(numb,1)=NaN;
    end
    
    
%         figure(1)
%         hold on
%         plot(ti([d11 d22]),trg([d11 d22]),'g.','MarkerSize',5)
%         plot(ti,tremorxf)
    
    %--------------------------------------------------
    
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
    
    clear start ending
    start=floor((start1./samplerateold)*samplerate);
    ending=floor((ending1./samplerateold)*samplerate);%floor(5*samplerate);
    clear xx
    xx{1,1}=xx1;
    
    
%     plot(ti(start),trg(1,start),'g.','MarkerSize',10)
%     plot(ti(ending),trg(1,ending),'g.','MarkerSize',10)
    
    
    m_notrig(numb,:)=nanmedian(low_tre);
    m_lowtrig(numb,:)=nanmedian(pp);
    m_all(numb,:)=nanmedian(pp_all);
    var_th(numb,:)=nanmedian(varb);
    
    figure(1)
   subplot(5,2,numb)
    histogram(pp_all,'FaceColor','w')
    hold on
    if ~isnan(pp)
        histogram(pp,'FaceColor','r')
    end
    if  ~isnan(low_tre)
        histogram(low_tre,'FaceColor','g')
    end

% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\newnonstim2.mat')
% cohort=[2 3 4 5 8 10 11 13 16 17];
% dum=intersect(iiii,cohort);
% 
% pt=[];
% for i=1:length(dum)
%     pt=[pt find(iiii==dum(i))];
% end
% 
%   main=[1 1 3 1 3 3 3 3 1 1];
%    figure(numb)
%    histogram(nostim(pt(numb),main(numb),:),'FaceColor','w')
%    hold on
%    histogram(pp_all,'FaceColor','r')
   
end


clearvars -except iiii  m_all m_notrig m_lowtrig th t_below var_th t_var
 cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
% cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
%  save('cleaned_tremor_prop.mat')



