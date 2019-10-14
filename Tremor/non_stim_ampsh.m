clear all
iii=[1 2 3 4 5 8 10 11 12 13 16];
in2=1;
for numb=11
    %     1:length(iii);
    Fpeak_tremor
    clearvars -except iii PC A1 B1 numb NS NS_i Fpeak samplerate in2
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
    PC=[70 66 47 47 47 50 50 50 50 55 45];
    A1={([1 3 6 8 12 18 23 27 30 32]);[];[];([2 3 5 6 7]);([1 2 4 6 8 9 10 11]);[];([1:9 15]);([2 4 7:10 13:15 22 25]);[1 4 7 12 25 31 45 47];[];[1 3 4 7 8]};
    B1={([2 5 7 11 17 22 26 29 31 34]);[];[];([2 3 5 6 7]);([1 2 4 6 7 9 10 11 12]);[];([1:9 15]);([2 5 7 8 9 12 13 14 19 22 25]);[2 5 10 14 26 37 46 50];[];[1 3 4 7 10]};
    
    if in2==1
        in=3;
    elseif in2==2 % other axis 1
        in=5;
    elseif in2==3 % other axis 2
        in=6;
    end
    data=SmrData.WvData;
    samplerateold=SmrData.SR;
    tremor=(data(in,:));
    addon=92; addon_end=35;
    
    
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremor2(1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    
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
    %  A=mean(C);
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
%     
%     plot(C)
%     hold on
%     plot(ind_s2,C(ind_s2),'r.')
%     plot(ind_e2,C(ind_e2),'b.')
    
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
    unstable3=[];
    
%     
%     plot(C)
%     hold on
%     plot(AA,C(AA),'r.')
%     plot(BB,C(BB),'b.')
%     
    %%% analysis
    
    if (Fpeak-2)>=1
        [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    else
        [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
    end
    
    tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
    dummy=hilbert(tremor_or);
    envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
    phase=angle(dummy);
    frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),251))';
    
    
    for i=1
        for j=1:5e4
            ix=randi(length(segmentb),1);
            segment=randi([segmentb(ix)+1000 segmente(ix)-5000],1);
            begin3=segment;
            end3=floor(begin3+5*samplerate);
            while ~isempty(intersect(unstable3,begin3:end3))
                segment=randi([segmentb(ix)+1000 segmente(ix)-5000],1);
                begin3=segment;
                end3=floor(begin3+5*samplerate);
            end
            baseline3(i,j)=(mean(envelope(end3-1000:end3))-mean(envelope(begin3-1000:begin3)))./mean(envelope(begin3-1000:begin3)); %#ok<*SAGROW> %
            baseline4(i,j)=(mean(frequency(end3-1000:end3))); %#ok<*SAGROW>
            
        end
    end
    
   rep=10; 
    for i=1:1e6
        dum=baseline3(randi(5e4,1,rep));
        dum2=dum;
        p(i)=nanmedian(dum2);
    end
    nostim(numb,:)=p;
    
    for i=1:12
        dum=baseline3(randi(5e4,1,rep));
        dum2=dum;
        nostimout(numb,i)=nanmedian(dum2);
    end
    clearvars -except nostimout iii numb PC A1 B1 iii stim nostim
end
% ANS_group=nostimout; clear nostimout
% AS_group=stim;
% clearvars  -except ANS_group AS_group
% cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
% save 'A_group'
NS=nostimout;
no_s=nostim;
clearvars  -except NS no_s
cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
% cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
save 'amp_NS.mat'
