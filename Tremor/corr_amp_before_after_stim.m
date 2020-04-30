clear all
close all

iii=[2 3 4 5 8 10 11 13 16 17];
for numb=1:length(iii);
    clearvars -except iii numb a_bf a_af tt_al
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))
    %     load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))
    
    in2=1; % analysing the "main tremor axis"
    if in2==1
        in=3;
    elseif in2==2 % other axis 1
        in=5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    elseif in2==3 % other axis 2
        in=6;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    end
    data=SmrData.WvData;
    samplerateold=SmrData.SR;
    tremor=(data(in,:));
    
    
    ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    tt=0:1/samplerate:(size(ds_data,2)-1)/samplerate;
    tre_3=ds_data([3 5 6],:); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    
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
    
    %     figure()
    %     plot(F(3:50),ps_curves(:,3:50)','LineWidth',2)
    %     legend({'z','y','x'})
    %     legend('boxoff')
    %     box('off')
    
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
    
    % amplitude
    tremor_or2=NaN(length(start),1);
    tremor_or3=NaN(length(start),1);
    for ax=1:3
        for i=1:length(start)
            if (~isnan(start(i)))
                tremor_or1(i,1)=mean(envelope(ax,start(i)-1000:start(i)));
                tremor_or2(i,1)=mean(envelope(ax,ending(i)-1000:ending(i)));
                tremor_or3(i,1)=(mean(envelope(ax,ending(i)-1000:ending(i)))-mean(envelope(ax,start(i)-1000:start(i))))/mean(envelope(ax,start(i)-1000:start(i)));
                
                
            else
                tremor_or1(i,1)=NaN;
                tremor_or2(i,1)=NaN;
                tremor_or3(i,1)=NaN;
            end
        end
        
        clear tt tt2
        a_b=NaN(20,12);
        a_a=NaN(20,12);
        tt_c=NaN(20,12);
        
        for i=1:12
            a_b(1:sum(xx==i),i)=tremor_or1(find(xx==i));
            a_a(1:sum(xx==i),i)=tremor_or2(find(xx==i));
            tt_c(1:sum(xx==i),i)=tremor_or2(find(xx==i));
        end
        
        a_bf{numb,ax}=a_b;
        a_af{numb,ax}=a_a;
        tt_all{numb,ax}=tt_c;
        
    end
    %     tt=abs(tt_c);
    %
    %     ttall(numb,:)=nanmedian(tt,2);
    %     ampall1(numb,:)=nanmedian(a_b,2);
    %
    %
    %    %%%% to be changed if used
    %
    %     n=[amp(:,max_ef(numb))  tt(:,max_ef(numb))];
    %
    %     figure(1)
    %     subplot(1,2,1)
    %     plot(n(:,1),n(:,2),'k+');
    %     y1=lsline;
    %     hold on
    %     ylabel('Change in amplitude ')
    %     box('off')
    %
    %
    %     xlabel('Amplitude')
    %     box('off')
    %
end



clear all
cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
load('ns_beforestim.mat')
load('rs_amp_bstim.mat')
load('pls_amp_bstim.mat')



f=[];
for i=1:10
    figure(1)
    subplot(2,5,i)
    
    c=[ns_befstim{i,1}(:,1) ; a_bf{i,1}(:,1) ;pls_amp_t{i,1}(:,1)];
    c=c(~isnan(c));
    d=(1:length(c))';
    y2=plot(d,c,'k.');
    y3=lsline;
    set(y3,'LineWidth',2,'Color','red')
    xlabel('Trials across conditions')
    ylabel('Tremor severity (m/s^2)')
    box('off')
    [c2 p2]=corrcoef(d,c)
     legend(y3,{num2str(round(c2,3))},'box','off','FontSize',12)
    if p2<0.05
        f=[f i];        
    end
end

% 
% figure(1)
% subplot(1,2,2)
% y2=plot(a1,a2,'k+');
% y3=lsline;
% set(y3,'LineWidth',2,'Color','red')
% box('off')
% c2=corrcoef(a1',a2')
% legend(y3,[num2str(c2(1,2))],'box','off')
% 
% for i=1:size(ttfall,1)
%     f1(1,i)=freqall(i,max_ef(i));
%     f2(1,i)=ttfall(i,max_ef(i));
% end
