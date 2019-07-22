clear all
iii=[1 2 3 4 5 6 8 10 11];
for numb=1:length(iii);
clearvars -except iii numb P_S ttall 
load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))
% load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))
start_clean;

%%%% re - estimate tremor characteristics
clear handup Pxx F frange Pxxrange Fpeak tremor_or variable envelope phase frequency

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

variable=hilbert(tremor_or);
envelope=sqrt((real(variable).^2)+(imag(variable).^2));
phase=angle(variable);
frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(variable))),251))'; % frequency estimated from hilbert

%%%% alternative frequency estimation using windowed FFT ( resolution of 0.2 Hz )
%%% there is a problem with both ways of estimating frequency either due to 
%%% using window or smoothing etc ... 

% for i=1:length(tremor_or)-5000
% seg=tremor_or(i:i+5000);
% f_seg=abs(real(fft(seg)));
% F=0:(1000/5000):(1000-0.2);
% f_alt(i)=F(min(find(f_seg==max(f_seg))));
% end
     
%%%% tremor_or2 - tremor_or22 is the difference in slope (i.e. change in frequency)

% tremor_or2=NaN(20,5001);
% tremor_or22=NaN(20,5001);
% 
% [a,b]=hist(frequency,0:0.05:10);
% 
% 
% for i=1:length(start)
%     if ~isnan(start(i)) 
%         tremor_or2(i,1:(ending(i)-start(i)+1))=unwrap(phase(start(i):ending(i)));
%         tremor_or22(i,1:(ending(i)-start(i)+1))=(phase(start(i))+(0:1:(ending(i)-start(i)))*2*pi/(1000./mean(frequency(start(i)-1000:start(i)))));
%         tremor_k(i,1)= (tremor_or2(i,(ending(i)-start(i)+1))-tremor_or22(i,(ending(i)-start(i)+1)))/(2*pi*0.001*(ending(i)-start(i))); %mean(frequency(ending(i)-1000:ending(i)));%
% %           tremor_k(i,1)= mean(frequency(ending(i)-1000:ending(i)));%
%     else
%         tremor_or22(i,1:5001)=NaN;
%         tremor_or2(i,1:5001)=NaN;
%         tremor_k(i,1)=NaN;
%     end
% end

%%%% tremor_or2 - tremor_or22 is the difference in frequency (estimated from unwrapped hilbert)


for i=1:length(start)
    if ~isnan(start(i)) 
        tremor_or2(i,1:(ending(i)-start(i)+1))=frequency(start(i):ending(i));
        tremor_or22(i,1:(ending(i)-start(i)+1))=mean(frequency(start(i)-1000:start(i)));
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
end

clear tt2

k=1;
tt3=NaN(20,12);
yy = xx ;

for s=1:size(tt,2)
    for i =1:1000;
        yy1=xx(randperm(size(xx,2)) );
        tt2(1:sum(yy1==s),1)=tremor_k(find(yy1==s));
        tt3(i,s)=nanmedian(tt2,1);
        clear tt2 
    end
end
P_S (numb,:)=nanmedian(tt);
ttall (numb,:,:)=tt3;
clearvars -except P_S tt3 ttall iii numb tt
cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\P_PS')
save (strcat('P0',num2str(iii(numb)),'_pha_suffle_phashift.mat'));
end
P_LS=ttall;
clearvars -except P_LS P_S
cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
save 'F_group.mat'



% rr(1:size(tt,2))=mean(prctile(tt3,95));
% rr1(1:size(tt,2))=mean(prctile(tt3,25));
% bar(nanmedian(tt))
% hold on
% plot(rr,'k--','LineWidth',1.5)
% plot(rr1,'k--','LineWidth',1.5)
% box('off')
% title 'Significant phasic-stimulation effect' 
% xticklabels({'0' '30' '60' '90' '120' '150' '180' '210' '240' '270' '300' '330'})
% 
% 
% 
% % find((nanmedian(tt))>(prctile(tt3,95))|(nanmedian(tt))<(prctile(tt3,25)))
% 
% figure()
% likhood_amp=sum(tt>prctile(tt3,95)| tt>0)./sum(~isnan(tt));
% likhood_sup=sum(tt<prctile(tt3,25)| tt<0)./sum(~isnan(tt));
% likhood=[likhood_sup ; likhood_amp];
% bar(likhood')
% title ('Likelihood of significant amplification/supression effect')
% xticklabels({'0' '30' '60' '90' '120' '150' '180' '210' '240' '270' '300' '330'})
% legend('supression','amplification')
% legend('boxoff')
% box('off')



% % % close all
% % % figure()
% % % fig=gcf;
% % % fig.Color=[1 1 1];
% % % bar(nanmedian(tt),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5])
% % % hold on
% % % plot(tt','.')
% % % xticklabels({'0','30','60','90','120','150','180','210','240','270','300','330'})
% % % % ylim([(-max((max(tt)))-0.1).*100 (max(max(tt))+0.1).*100])
% % % box('off')
% % % title ('P11')
% % % ylabel( 'Frequency change (Hz)')
% % % xlabel ('Stimulated phase')
% % % ylim ([-max(max(tt)) max(max(tt))])


% % figure()
% fig=gcf;
% fig.Color=[1 1 1];
% bar(100.*nanmedian(tt2))
% hold on
% stem((100.*tt2)')
% xticklabels({'0','30','60','90','120','150','180','210','240','270','300','330'})
% % ylim([(-max((max(tt2)))-0.1).*100 (max(max(tt2))+0.1).*100])
% box('off')
% title ('P11')
% 
% 
%