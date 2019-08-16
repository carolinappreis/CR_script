clear all
load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\bline_amp.mat')
iii=[1 2 3 4 5 8 10 11 12 13];

for numb=1:length(iii);
    clearvars -except iii numb bamp_ns
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

tremor_or2=NaN(length(start),1);
tremor_or3=NaN(length(start),1);

for i=1:length(start)
    if (~isnan(start(i)))
        tremor_or3(i,1)=mean(envelope(start(i)-1000:start(i)));
        tremor_or2(i,1)=(mean(envelope(ending(i)-1000:ending(i)))-mean(envelope(start(i)-1000:start(i))))/mean(envelope(start(i)-1000:start(i)));

    else
        tremor_or2(i,1)=NaN;
        tremor_or3(i,1)=NaN;
    end
end


% figure(1)
% subplot(1,length(iii),numb)
% plot(tremor_or3,'k+')
% y1=lsline
% set(y1,'LineWidth',2,'Color','red')
% c2=corrcoef([1:length(tremor_or3)]',tremor_or3)
% legend(y1,[num2str(c2(1,2))],'box','off')

dummy=bamp_ns(numb,:); dummy=dummy(find(~isnan(dummy)));
tremor_nss=[dummy tremor_or3'];
figure(1)
subplot(1,length(iii),numb)
plot(tremor_nss,'k+')
y1=lsline
set(y1,'LineWidth',2,'Color','red')
c2=corrcoef((1:length(tremor_nss))',(tremor_nss)');
legend(y1,[num2str(c2(1,2))],'box','off')

end



