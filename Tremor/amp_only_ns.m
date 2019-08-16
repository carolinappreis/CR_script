clear all
iii=[1 2 3 4 5 8 10 11 12 13];
PC=[70 66 47 47 47 50 50 50 50 55];
A1={([1 3 6 8 12 18 23 27 30 32]);[];[];([2 3 5 6 7]);([1 2 4 6 8 9 10 11]);[];([1:9 15]);([2 4 7:10 13:15 22 25]);[1 4 7 12 25 31 45 47];[]};
B1={([2 5 7 11 17 22 26 29 31 34]);[];[];([2 3 5 6 7]);([1 2 4 6 7 9 10 11 12]);[];([1:9 15]);([2 5 7 8 9 12 13 14 19 22 25]);[2 5 10 14 26 37 46 50];[]};

for numb=1:length(iii);
    
    %      load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
    %      load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))
    load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(iii(numb)),'_baseline.mat'))
    load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))
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
    
    [b,a]=butter(2,[0.1/(0.5*samplerate) ],'low'); %15
    C=(filtfilt(b,a,envelope));
    A=prctile(C,PC(numb));
    ind_s=[];
    for i=11:length(C)-11
        if C(i-1)<A && C(i+1)>A
            ind_s=[ind_s i]; %#ok<*AGROW>
        end
    end
    for i=1:(length(ind_s)-1)
        if ind_s(i+1)-ind_s(i)==1
            ind_s(i+1)=NaN;
        end
    end
    ind_s2=ind_s(~isnan(ind_s));
    ind_e=[];
    for i=11:length(C)-11
        if C(i-1)>A && C(i+1)<A
            ind_e=[ind_e i];
        end
    end
    for i=1:(length(ind_e)-1)
        if ind_e(i+1)-ind_e(i)==1
            ind_e(i+1)=NaN;
        end
    end
    ind_e2=ind_e(~isnan(ind_e));
    
    if isempty (A1{numb,1})
        AA=ind_s2;
        BB=ind_e2;
    else
        AA=ind_s2(A1{numb,1});
        BB=ind_e2(B1{numb,1});
    end
    
    if numb==5
        AA= [1 ind_s2(A1{numb,1})];
        BB=ind_e2(B1{numb,1});
    end
    
    
    segmentb=AA;
    segmente=BB;
    
    %     plot(C)
    %     hold on
    %     plot(AA,C(AA),'r.')
    %     plot(BB,C(BB),'b.')
    
    %%% analysis
    
    for i=1:size(segmentb,2)
        amp_ns(numb,i)=median(envelope(segmentb(i):segmente(i)));
    end
    
    clearvars -except nostimout iii numb PC A1 B1 iii amp_ns
end
amp_ns(amp_ns==0)= NaN;
bline_amp=amp_ns;
clearvars -except bline_amp
cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
save ('bline_amp.mat')
