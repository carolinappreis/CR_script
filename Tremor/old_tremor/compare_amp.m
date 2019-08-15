% clear all
iii=[1 2 3 4 5 8 10 11];

for numb=1:length(iii);
    clearvars -except iii numb amp_c tremor_c
    % load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))
    
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
    
    
    
    start1=NaN(2,round(numel(start)./2));
    
    start1(1,:)=start(1:round(numel(start)./2));
    start1(2,1:length((round(numel(start)./2)+1):numel(start)))=start((round(numel(start)./2)+1):numel(start));
    start=start1;
    
    ending1=NaN(2,round(numel(ending)./2));
    
    ending1(1,:)=ending(1:round(numel(ending)./2));
    ending1(2,1:length((round(numel(ending)./2)+1):numel(ending)))=ending((round(numel(ending)./2)+1):numel(ending));
    clear ending;
    ending=ending1;
    
    
    for ss=1:size(start,1)
        for i=1:length(start(ss,:))
            if (~isnan(start(ss,i)))
                amp(ss,i,1)=mean(envelope(start(ss,i)-1000:start(ss,i)));
                tremor_or2(ss,i,1)=(mean(envelope(ending(ss,i)-1000:ending(ss,i)))-mean(envelope(start(ss,i)-1000:start(ss,i))))/mean(envelope(start(ss,i)-1000:start(ss,i)));
                
            else
                amp(ss,i,1)=NaN;
                tremor_or2(ss,i,1)=NaN;
            end
        end
        
    end
    
    amp_c(numb,:)=(nanmean(amp,2))';
    tremor_c(numb,:)=(nanmean(tremor_or2,2))';
end

ttest(amp_c(:,2),amp_c(:,1))