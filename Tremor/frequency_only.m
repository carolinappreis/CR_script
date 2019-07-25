clear all
iii=[1 2 3 4 5 6 8 10 11];

for numb=1:length(iii);
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

dummy=hilbert(tremor_or);
envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
env_z=zscore(envelope);
phase=angle(dummy);
frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
% freq_z=zscore(abs(frequency));

for i=1:length(start)
    if ~isnan(start(i)) 
        tremor_or2(i,1:(ending(i)-start(i)+1))=freq_z(start(i):ending(i));
        tremor_or22(i,1:(ending(i)-start(i)+1))=mean(freq_z(start(i)-1000:start(i)));
        tremor_k(i,1)=tremor_or2(i,(ending(i)-start(i)+1))-tremor_or22(i,(ending(i)-start(i)+1));
    else
        tremor_or22(i,1:5001)=NaN;
        tremor_or2(i,1:5001)=NaN;
        tremor_k(i,1)=NaN;
    end
end

clear tt
tt=[]; 
k=1;

tt=NaN(20,12);

for i=1:12
    tt(1:sum(xx==i),i)=tremor_k(find(xx==i));
    freq(1:sum(xx==i),i)=tremor_or22(find(xx==i));
end

clear tt2

ttall (numb,:)=nanmedian(tt);
freqall (numb,:)=nanmedian(freq);

end
clearvars -except ttall freqall 
cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim')
% cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim')
save ('abs_freq.mat')


plot(mean(freqall,2),mean(ttall,2),'k+')
lsline
m=[mean(freqall,2) mean(ttall,2)];
corrcoef(m(:,1),m(:,2))


for i=1:9
a1(1,i)=ttall(i,phase_peak(i));
a2(1,i)=ampall(i,phase_peak(i));
end
plot(a1,a2,'k+')
lsline
corrcoef(a1',a2')