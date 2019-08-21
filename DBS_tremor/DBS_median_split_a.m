clear all
iii=[1];

for numb=1;
%     :length(iii);
    clearvars -except iii numb ttall ampall ph_stim LS
    load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))


DBS_cleaner;

%%% re - estimate tremor characteristics
clear handup Pxx F frange Pxxrange Fpeak tremor_or dummy envelope phase frequency 

handup=[];
for i=1:length(start)
    handup=[handup start(i):ending(i)]; %#ok<*AGROW>
end
handup=sort(handup,'ascend');

[Pxx,F]=pwelch(tremor2(handup),samplerate,[],samplerate,samplerate);

frange=F(3:10);
Pxxrange=Pxx(3:10);

Fpeak=frange(find(Pxxrange==max(Pxxrange)));

if (Fpeak-2)>=1
    [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
    [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end
tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
% tremor_or=zscore(tremor_or);
dummy=hilbert(tremor_or);
envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
phase=angle(dummy);
frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
    
    % amplitude
    tremor_or2=NaN(length(start),1);
    for i=1:length(start)
        if (~isnan(start(i)))
            tremor_or3(i,1)=mean(envelope(start(i)-1000:start(i)));
            tremor_or2(i,1)=(mean(envelope(ending(i)-1000:ending(i)))-mean(envelope(start(i)-1000:start(i))))/mean(envelope(start(i)-1000:start(i)));
            
        else
            tremor_or2(i,1)=NaN;
            tremor_or3(i,1)=NaN;
        end
    end
    
    %%% criteria for outliers

idx_outl=find(tremor_or2>(nanmean(tremor_or2)+(2*(nanstd(tremor_or2))))|tremor_or2<(nanmean(tremor_or2)-(2*(nanstd(tremor_or2)))));
    tremor_or2(idx_outl,1)=NaN;
    tremor_or3(idx_outl,1)=NaN;
    xx(1,idx_outl)=NaN;
    
    
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

    for i=size(effect1,2)+1:size(effect1,2)*2
        arc1(numb,i-12)=nansum(ef_1(1,(i-1:i+1)))./length(ef_1(1,(i-1:i+1)));
        arc2(numb,i-12)=nansum(ef_2(1,(i-1:i+1)))./length(ef_2(1,(i-1:i+1)));
    end
end

cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data')
clearvars -except arc1 arc2
save('arc_mediansplit.mat')

load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
cl=blushred;
cl1=squash;

for i=1:size(arc1,1)
    f1=figure(1)
    subplot(1,size(arc1,1),i)
    bar(arc2(i,:),'FaceColor',cl,'EdgeColor',cl)
    box('off')
    hold on
    bar(arc1(i,:),'LineStyle','--','LineWidth',1,'FaceColor','none','EdgeColor','k')
    box('off')
end

% f1.Units = 'centimeters';
% f1.OuterPosition= [10, 10, 60, 6];
% set(gca,'FontSize',9)
% set(f1,'color','w');