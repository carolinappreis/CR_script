clear all
close all

iii=[2 3 4 5 8 10 11 13 16 17];
for numb=1:length(iii);
    clearvars -except iii numb arc1 arc2 
         load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))
%     load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))
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
    
    %     %%%% find runs with trigering issues (too few, too many pulses)
    th1=(Fpeak*5)./2;
    th2=(Fpeak*5)+round((Fpeak*5)./2);
    
    %         th1=(Fpeak*5*5)./2;
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
    start=floor((start1./samplerateold)*samplerate);
    ending=floor((ending1./samplerateold)*samplerate);%floor(5*samplerate);
    clear xx
    xx=xx1;
    
    % amplitude
    tremor_or2=NaN(length(start),1);
    tremor_or3=NaN(length(start),1);
    for i=1:length(start)
        if (~isnan(start(i)))
            tremor_or3(i,1)=mean(envelope(1,start(i)-1000:start(i)));
            tremor_or2(i,1)=(mean(envelope(1,ending(i)-1000:ending(i)))-mean(envelope(1,start(i)-1000:start(i))))/mean(envelope(1,start(i)-1000:start(i)));
            env_pha(i,1:5000)=zenv(1,start(i):start(i)+5000-1);
        else
            tremor_or2(i,1)=NaN;
            tremor_or3(i,1)=NaN;
            env_pha(i,1:5000)=NaN;
        end
    end
    
    %%% criteria for outliers
    
    %     idx_outl=find(tremor_or2>(mean(tremor_or2+2*(std(tremor_or2))))|tremor_or2<(mean(tremor_or2-2*(std(tremor_or2)))));
    %     tremor_or2(idx_outl,1)=NaN;
    %     tremor_or3(idx_outl,1)=NaN;
    %     xx(1,idx_outl)=NaN;
    %
    
    amp_1=NaN(2,round(size(tremor_or3,1)./2));
    ch_a1=NaN(2,round(size(tremor_or3,1)./2));
    pha_idx=NaN(2,round(size(tremor_or3,1)./2));
    m=1;
    n=1;
    for i=1:length(start)
        if tremor_or3(i)<=nanmedian(tremor_or3)
            amp_1(1,n)= tremor_or3(i);
            ch_a1(1,n)= tremor_or2(i);
            pha_idx(1,n)=xx(i);
            n=n+1;
            
        else
            amp_1(2,m)= tremor_or3(i);
            ch_a1(2,m)= tremor_or2(i);
            pha_idx(2,m)=xx(i);
            m=m+1;
            
        end
    end
    
    clear tt1 tt2 amp1 amp2
    tt1=NaN(20,12);
    amp1=NaN(20,12);
    tt2=NaN(20,12);
    amp2=NaN(20,12);
    
    for i=1:12
        tt1(1:sum(pha_idx(1,:)==i),i)=ch_a1(1,find(pha_idx(1,:)==i));
        amp1(1:sum(pha_idx(1,:)==i),i)=amp_1(1,find(pha_idx(1,:)==i));
        tt2(1:sum(pha_idx(2,:)==i),i)=ch_a1(2,find(pha_idx(2,:)==i));
        amp2(1:sum(pha_idx(2,:)==i),i)=amp_1(2,find(pha_idx(2,:)==i));
    end
    
    effect1=nanmedian(tt1);
    effect2=nanmedian(tt2);
    
    ef_1=repmat(effect1,1,3);
    ef_2=repmat(effect2,1,3);
    
    %     for i=size(effect1,2)+1:size(effect1,2)*2
    %         arc1(numb,i-12)=nansum(ef_1(1,(i-1:i+1)))./length(ef_1(1,(i-1:i+1)));
    %         arc2(numb,i-12)=nansum(ef_2(1,(i-1:i+1)))./length(ef_2(1,(i-1:i+1)));
    %     end
    
    for i=size(effect1,2)+1:size(effect1,2)*2
        arc1(numb,i-12)=sum(ef_1(1,(i-1:i+1)))./length(ef_1(1,(i-1:i+1)));
        arc2(numb,i-12)=sum(ef_2(1,(i-1:i+1)))./length(ef_2(1,(i-1:i+1)));
    end
end

% cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
% cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
% clearvars -except arc1 arc2 iii
% save('arc_mediansplit_0220.mat')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/arc_mediansplit_0120.mat')
% load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash')

cl=blushred;
cl1=squash;

for i=1:size(arc1,1)
    f1=figure(1)
    subplot(2,5,i)
    bar(arc2(i,:),'FaceColor',cl,'EdgeColor',cl)
    box('off')
    hold on
    bar(arc1(i,:),'LineStyle','--','LineWidth',1,'FaceColor','none','EdgeColor','k')
    box('off')
end

f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 60, 6];
set(gca,'FontSize',9)
set(f1,'color','w');




load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','grey')
%  load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','aegean','grey');


cl2=[0.7 0.7 0.7];


for hh=1:size(arc1,1)
    f1=figure(1)
    subplot(2,5,hh)
    y=[arc1(hh,:);arc2(hh,:)]';
    b = bar(0:30:330,y);
    yline(0,'LineWidth',1)
    ylim([-0.75 0.75])
    yticks([ -1:0.25:1])
    box('off')
    ylabel('Change in tremor severity')
    xlabel('Stimulation phase (degrees)')
    
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 12, 12];
    set(gca,'XTickLabelRotation',45)
    set(gca,'FontSize',14)
    set(f1,'color','w');
    
    
end



f1=figure()
y=[arc1{hh,1};arc2{hh,1}]';
b = bar(0:30:330,y);
yline(0,'LineWidth',1)
ylim([-0.75 0.75])
yticks([ -1:0.25:1])
box('off')
ylabel('Change in tremor severity')
xlabel('Stimulation phase (degrees)')
    
f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 12, 12];
set(gca,'XTickLabelRotation',45)
set(gca,'FontSize',14)
set(f1,'color','w');