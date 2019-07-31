% clear all
iii=[1 2 3 4 5 8 10 11];

for numb=1:length(iii);
    clearvars -except iii numb ttall ampall
load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))
% load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))

start_clean;

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

tt=NaN(20,12);
amp=NaN(20,12);
yy = xx ;

for i=1:12
    tt(1:sum(xx==i),i)=tremor_or2(find(xx==i));
    amp(1:sum(xx==i),i)=tremor_or3(find(xx==i));
end

ttall (numb,:)=nanmedian(tt);
ampall (numb,:)=nanmedian(amp);
% label_sh(numb,:,:)= label_sh1;
clearvars -except ttall iii numb ampall 

end
clearvars -except ttall ampall 
cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim')
% cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim')
save ('abs_amp.mat')
plot(mean(ampall,2),mean(ttall,2),'k+')
lsline
corrcoef(m(:,1),m(:,2))


for i=1:9
a1(1,i)=ttall(i,phase_peak(i));
a2(1,i)=ampall(i,phase_peak(i));
end
plot(a1,a2,'k+')
lsline
corrcoef(a1',a2')