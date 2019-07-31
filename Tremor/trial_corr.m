 clear all
% iii=[1 2 3 4 5 6 8 10 11];
iii=[1 2 3 4 5 8 10 11];

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
z_tremor=zscore(tremor_or);
dummy=hilbert(z_tremor);
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

clear tt tt2
tt1=NaN(20,12);
amp=NaN(20,12);

for i=1:12
    tt1(1:sum(xx==i),i)=tremor_or2(find(xx==i));
    amp(1:sum(xx==i),i)=tremor_or3(find(xx==i));
end
tt=abs(tt1);

ttall(numb,:)=nanmedian(tt);
ampall(numb,:)=nanmedian(amp);
max_ef(1,numb)=find(ttall(numb,:)==max(ttall (numb,:)));
n=[amp(:,max_ef(numb))  tt(:,max_ef(numb))];

figure(1)
plot(n(:,1),n(:,2),'k+');
y1=lsline;
hold on
ylabel ('Change in amplitude (zscore)')
xlabel('Amplitude (zscore)')

%%% frequency 

% for i=1:length(start)
%     if ~isnan(start(i)) 
%         tremor_for2(i,1:(ending(i)-start(i)+1))=frequency(start(i):ending(i));
%         tremor_for22(i,1:(ending(i)-start(i)+1))=mean(frequency(start(i)-1000:start(i)));
%         tremor_k(i,1)=tremor_for2(i,(ending(i)-start(i)+1))-tremor_for22(i,(ending(i)-start(i)+1));
%     else
%         tremor_for22(i,1:5001)=NaN;
%         tremor_for2(i,1:5001)=NaN;
%         tremor_k(i,1)=NaN;
%     end
% end

%%% frequency 

tremor_for2=NaN(20,5001);
tremor_for22=NaN(20,5001);

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


ttf=[]; 
nf=[];
k=1;

ttf1=NaN(20,12);
freq=NaN(20,12);

for i=1:12
    ttf1(1:sum(xx==i),i)=tremor_k(find(xx==i));
    freq(1:sum(xx==i),i)=tremor_for22(find(xx==i));
end
ttf=abs(ttf1);

ttfall (numb,:)=nanmedian(ttf);
freqall (numb,:)=nanmedian(freq);
max_fef(1,numb)=find(ttfall(numb,:)==max(ttfall (numb,:)));

nf=[freq(:,max_fef(numb))  ttf(:,max_fef(numb))];

% figure(2)
% plot(nf(:,1),nf(:,2),'k+');
% lsline
% hold on
% ylabel ('Change in frequency (zscore)')
% xlabel('Frequency (zscore)')

figure(3)
plot(n(:,1),nf(:,2),'k+');
lsline
hold on
ylabel ('Change in frequency (zscore)')
xlabel('Amplitude (zscore)')
% 
% figure(4)
% plot(nf(:,1),n(:,2),'k+');
% lsline
% hold on
% ylabel ('Change in amplitude (zscore)')
% xlabel('Frequency (zscore)')

end

for i=1:size(ttall,1)
a1(1,i)=ttall(i,max_ef(i));
a2(1,i)=ampall(i,max_ef(i));
end

figure(5)
y2=plot(a1,a2,'k+');
y3=lsline;
set(y3,'LineWidth',2,'Color','red')
box('off')
ylabel ('Change in amplitude (zscore)')
xlabel('Amplitude (zscore)')
c2=corrcoef(a1',a2')
legend(y3,[num2str(c2(1,2))],'box','off')

for i=1:size(ttfall,1)
f1(1,i)=ttfall(i,max_fef(i));
f2(1,i)=freqall(i,max_fef(i));
end

% figure(6)
% y2=plot(f1,f2,'k+');
% y3=lsline;
% set(y3,'LineWidth',2,'Color','red')
% box('off')
% ylabel ('Change in frequency (zscore)')
% xlabel('Frequency (zscore)')
% c2=corrcoef(f1',f2')
% legend(y3,[num2str(c2(1,2))],'box','off')

figure(7)
y2=plot(a1,f2,'k+');
y3=lsline;
set(y3,'LineWidth',2,'Color','red')
box('off')
ylabel ('Change in frequency (zscore)')
xlabel('Amplitude (zscore)')
c2=corrcoef(a1',f2')
legend(y3,[num2str(c2(1,2))],'box','off')

% figure(8)
% y2=plot(f1,a2,'k+');
% y3=lsline;
% set(y3,'LineWidth',2,'Color','red')
% box('off')
% ylabel ('Change in amplitude (zscore)')
% xlabel('Frequency (zscore)')
% c2=corrcoef(f1',a2')
% legend(y3,[num2str(c2(1,2))],'box','off')
