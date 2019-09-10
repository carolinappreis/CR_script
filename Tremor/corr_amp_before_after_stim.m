clear all
close all
iii=[1 2 3 4 5 8 10 11 12 13];
for numb=1:length(iii);
    clearvars -except iii numb ttall ampall  max_ef  ttfall freqall
%     load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))
    load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))
    
    in2=1
    start_cleaner;

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
    phase=angle(dummy);
    frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
    
    % amplitude
    tremor_or2=NaN(length(start),1);
    tremor_or3=NaN(length(start),1);
    for i=1:length(start)
        if (~isnan(start(i)))
            tremor_or1(i,1)=mean(envelope(start(i)-1000:start(i)));
            tremor_or2(i,1)=mean(envelope(ending(i)-1000:ending(i)));
            tremor_or3(i,1)=(mean(envelope(ending(i)-1000:ending(i)))-mean(envelope(start(i)-1000:start(i))))/mean(envelope(start(i)-1000:start(i)));


        else
            tremor_or1(i,1)=NaN;
            tremor_or2(i,1)=NaN;
            tremor_or3(i,1)=NaN;
        end
    end
    
    clear tt tt2
    tt1=NaN(20,12);
    tt2=NaN(20,12);
    tt3=NaN(20,12);
    
    for i=1:12
        tt1(1:sum(xx==i),i)=tremor_or1(find(xx==i));
        tt2(1:sum(xx==i),i)=tremor_or2(find(xx==i));
        tt3(1:sum(xx==i),i)=tremor_or2(find(xx==i));
    end
    tt=abs(tt3);
    
    ttall(numb,:)=nanmedian(tt);
    ampall1(numb,:)=nanmedian(tt1);
    ampall2(numb,:)==nanmedian(tt2);
    max_ef(1,numb)=find(ttall(numb,:)==max(ttall (numb,:)));
    
   %%%% to be changed if used
    
    n=[amp(:,max_ef(numb))  tt(:,max_ef(numb))];
    
    figure(1)
    subplot(1,2,1)
    plot(n(:,1),n(:,2),'k+');
    y1=lsline;
    hold on
    ylabel('Change in amplitude ')
    box('off')
    
    
    xlabel('Amplitude')
    box('off')
    
end

for i=1:size(ttall,1)
    a1(1,i)=ampall(i,max_ef(i));
    a2(1,i)=ttall(i,max_ef(i));
end

figure(1)
subplot(1,2,2)
y2=plot(a1,a2,'k+');
y3=lsline;
set(y3,'LineWidth',2,'Color','red')
box('off')
c2=corrcoef(a1',a2')
legend(y3,[num2str(c2(1,2))],'box','off')

for i=1:size(ttfall,1)
    f1(1,i)=freqall(i,max_ef(i));
    f2(1,i)=ttfall(i,max_ef(i));
end
