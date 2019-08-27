clear all
iii=[1];

for numb=1;
%     :length(iii);
    clearvars -except iii numb ttall ampall ph_stim LS tt1
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

tremor_or2=NaN(20,5001);
tremor_or22=NaN(20,5001); 
% [a,b]=hist(frequency,0:0.05:10);

for i=1:length(start)
    if ~isnan(start(i)) 
        tremor_or2(i,1:(ending(i)-start(i)+1))=unwrap(phase(start(i):ending(i)));
        tremor_or22(i,1:(ending(i)-start(i)+1))=(phase(start(i))+(0:1:(ending(i)-start(i)))*2*pi/(1000./mean(frequency(start(i)-1000:start(i)))));
        tremor_k(i,1)= (tremor_or2(i,(ending(i)-start(i)+1))-tremor_or22(i,(ending(i)-start(i)+1)))/(2*pi*0.001*(ending(i)-start(i))); %mean(frequency(ending(i)-1000:ending(i)));%
    else
        tremor_or22(i,1:5001)=NaN;
        tremor_or2(i,1:5001)=NaN;
        tremor_k(i,1)=NaN;
    end
end

% idx_outl=find(tremor_k>(nanmean(tremor_k)+(2*(nanstd(tremor_k))))|tremor_k<(nanmean(tremor_k)-(2*(nanstd(tremor_k)))));
% tremor_k(idx_outl,1)=NaN;
% xx(1,idx_outl)=NaN;

clear tt
clear freq

tt=NaN(25,12);
% freq=NaN(25,12);

for i=1:12
    tt(1:sum(xx==i),i)=tremor_k(find(xx==i));
%     freq(1:sum(xx==i),i)=tremor_or22(find(xx==i));
end

tt1{numb,1}=tt;
ttall (numb,:)=nanmedian(tt);
% freqall (numb,:)=nanmedian(freq);

% for s=1:size(tt,2)
%      for i =1:100000;
%         yy1=xx(randperm(size(xx,2)));
%         tt2(1:sum(yy1==s),1)=tremor_k(find(yy1==s));
%         tt3(i,s)=nanmedian(tt2,1);
%         clear tt2 
%     end
% end
% LS (numb,:,:)=tt3;

for rr=1:100000
LS(numb,rr)=nanmedian(tt(randi(length(start),1,10)));
end

clearvars -except ttall iii numb freqall ph_stim LS tt1

end
clearvars -except ttall  ph_stim LS tt1
cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data')
save('DBS_freq_FRC.mat')

