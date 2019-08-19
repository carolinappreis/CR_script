clear all
close all
iii=[1 2 3 4 5 8 10 11 12 13];

for numb=1:length(iii);
    clearvars -except iii numb frc1 frc2
load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))
% load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))

start_cleaner;

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
    
    
    tremor_for2=NaN(20,5001);
    tremor_for22=NaN(20,5001);
    
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
    
    [a,b]=hist(frequency,0:0.05:10);
    
    
    for i=1:length(start)
        if ~isnan(start(i))
            tremor_for2(i,1:(ending(i)-start(i)+1))=unwrap(phase(start(i):ending(i)));
            tremor_for22(i,1:(ending(i)-start(i)+1))=(phase(start(i))+(0:1:(ending(i)-start(i)))*2*pi/(1000./mean(frequency(start(i)-1000:start(i)))));
            tremor_k(i,1)= (tremor_for2(i,(ending(i)-start(i)+1))-tremor_for22(i,(ending(i)-start(i)+1)))/(2*pi*0.001*(ending(i)-start(i))); %mean(frequency(ending(i)-1000:ending(i)));%
            %           tremor_k(i,1)= mean(frequency(ending(i)-1000:ending(i)));%
        else
            tremor_for22(i,1:5001)=NaN;
            tremor_for2(i,1:5001)=NaN;
            tremor_k(i,1)=NaN;
        end
    end
    
    amp_1=NaN(2,round(size(tremor_or3,1)./2));
    ch_f1=NaN(2,round(size(tremor_k,1)./2));
    pha_idx=NaN(2,round(size(tremor_k,1)./2));
    m=1;
    n=1;
    for i=1:length(start)
        if tremor_or3(i)<=nanmedian(tremor_or3)
            amp_1(1,n)= tremor_or3(i);
            ch_f1(1,n)= tremor_k(i);
            pha_idx(1,n)=xx(i);
            n=n+1;
            
        else
            amp_1(2,m)= tremor_or3(i);
            ch_f1(2,m)= tremor_k(i);
            pha_idx(2,m)=xx(i);
            m=m+1;
            
        end
    end
 
    clear tt1 tt2 freq1 freq2
    tt1=NaN(20,12);
    freq1=NaN(20,12);
    tt2=NaN(20,12);
    freq2=NaN(20,12);

    for i=1:12
        tt1(1:sum(pha_idx(1,:)==i),i)=ch_f1(1,find(pha_idx(1,:)==i));
        freq1(1:sum(pha_idx(1,:)==i),i)=amp_1(1,find(pha_idx(1,:)==i));
        tt2(1:sum(pha_idx(2,:)==i),i)=ch_f1(2,find(pha_idx(2,:)==i));
        freq2(1:sum(pha_idx(2,:)==i),i)=amp_1(2,find(pha_idx(2,:)==i));
    end

    effect1=nanmedian(tt1);
    effect2=nanmedian(tt2);

    ef_1=repmat(effect1,1,3);
    ef_2=repmat(effect2,1,3);

    for i=size(effect1,2)+1:size(effect1,2)*2
        frc1(numb,i-12)=nansum(ef_1(1,(i-1:i+1)))./length(ef_1(1,(i-1:i+1)));
        frc2(numb,i-12)=nansum(ef_2(1,(i-1:i+1)))./length(ef_2(1,(i-1:i+1)));
    end
end

% clearvars -except frc1 frc2
% save('frc_mediansplit.mat')

close all

load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','aegean','stone');
cl=aegean;
cl1=stone;

for i=1:size(frc1,1)
    f1=figure(1)
    subplot(1,size(frc1,1),i)
    bar(frc2(i,:),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5])
    hold on
    bar(frc1(i,:),'LineStyle','--','LineWidth',1,'FaceColor','none','EdgeColor',cl)
    
end




f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 60, 6];
set(gca,'FontSize',9)
set(f1,'color','w');