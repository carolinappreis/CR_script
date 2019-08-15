clear all

iii=[1 2 3 4 5 8 10 11 13];

for numb= 4
%     1:length(iii);
    clearvars -except iii numb idx_trg
    % load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))
    
    
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
    
    %%% downsample
    
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremor2(1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    
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
    
    Fpeak=frange(find(Pxxrange==max(Pxxrange)));
    th1=(Fpeak*5)./2;
    th2=(Fpeak*5)+round((Fpeak*5)./5);
    %-------
    
    new=find(data(2,:)>4);
    difp=find((diff(new))>100000);
    ep=[new(difp) new(end)];
    sp=[new(1) new(difp+1)];
    
    %     plot(time,data(2,:))
    %     hold on
    %     plot(time(sp),data(2,sp),'r.')
    %     plot(time(ep),data(2,ep),'k.')
    
    for ik=1:length(sp)
        
        s=(find(([indexes4>=sp(ik)]+[indexes4<=ep(ik)])==2));
        e=(find(([indexes3>=sp(ik)]+[indexes3<=ep(ik)])==2));
        tks=(find(diff(xx(s))==0))+1;
        tke=(find(diff(xx(e))==0));
        indexes4(s(tks))=NaN;
        indexes3(e(tke))=NaN;
        
        if length(sp(ik):ep(ik))>700000
            do=data(2,sp(ik):ep(ik));
            d=zeros(1,length(sp(ik):ep(ik)));
            td=0:1/1000:((length(sp(ik):ep(ik)))-1)/1000;
            
            bins= [1:(length(sp(ik):ep(ik))./12):length(sp(ik):ep(ik))];
            for l=1:length(bins)
                
                d(bins(l):bins(l)+50000)=5;
                
                z=do(bins(l):bins(l)+50000);
                n_trg(1,l)=numel(find(diff(z)<-3));
                
                 if numel(find(diff(z)<-3))>=th1 & numel(find(diff(z)<-3))<th2
%                                     if numel(find(diff(z)<-3))>=th1 

                    good_t(ik,l)=1;
                else
                    good_t(ik,l)=NaN;
                end
                clear z
            end
            %         plot(td,do)
            %         hold on
            %         plot(td,d,'r','LineWidth',2)
        else
            good_t(ik,:)=NaN;
        end
        clearvars do d dt
    end
    
    idx_trg=reshape(good_t',1,size(good_t,1)*size(good_t,2));
    
    
   

end


for i=1:length(idx_trg)
    
    if isnan (idx_trg(i))
        indexes3(i)=NaN;
        indexes4(i)=NaN;
    end
end
     
    indexes4=indexes4(~isnan(indexes4));
    indexes3=indexes3(~isnan(indexes3));
    figure(1)
    plot(time,data(4,:))
    hold on
    plot(time(index),data(4,index),'r.')
    plot(time(indexes4),data(4,indexes4),'ko')
    plot(time(indexes3),data(4,indexes3),'bo')
%     xlim([0 100])
clearvars  start end
start=floor((indexes4./samplerateold)*samplerate)+addon;
ending=floor((indexes3./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);