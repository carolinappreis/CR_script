clear all
 iii=[1 2 3 4 5 8 10 11 12 13];

PC=[70 66 47 47 47 50 50 50 50 55];
A1={([1 3 6 8 12 18 23 27 30 32]);[];[];([2 3 5 6 7]);([1 2 4 6 8 9 10 11]);[];([1:9 15]);([2 4 7:10 13:15 22 25]);[];[]};
B1={([2 5 7 11 17 22 26 29 31 34]);[];[];([2 3 5 6 7]);([1 2 4 6 7 9 10 11 12]);[];([1:9 15]);([2 5 7 8 9 12 13 14 19 22 25]);[];[]};
for numb=length(iii)-1;
    clearvars -except iii PC A1 B1 numb NS NS_i
    
         load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(iii(numb)),'_baseline.mat'))
         load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))
%     load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(iii(numb)),'_baseline.mat'))
%     load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))
    
    
    in2=1; % analysing the "main tremor axis"
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
         cd('C:\Users\creis\Documents\GitHub\CR_script\Tremor\old_tremor')
%     cd('/Users/Carolina/Documents/GitHub/CR_script/Tremor')

    run ('phasedetection.m');
    data=SmrData.WvData;
    rep=10; % number of trials for random stim - please enter for each patient
    clearvars -except Fpeak in2 in rep SmrData data nostim iii numb PC A1 B1 NS NS_i xx
    
    samplerateold=SmrData.SR;
    time=0:1/samplerateold:(size(data,2)-1)/samplerateold;
    tremor=data(in,:);
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
    dummy=hilbert(tremor_or);
    envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
    
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
    % plot(C)
    % hold on
    % plot(ind_s2,C(ind_s2),'r.')
    % plot(ind_e2,C(ind_e2),'b.')
    
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
    %%% analysis
    plot(C)
hold on
plot(AA,C(AA),'r.')
plot(BB,C(BB),'b.')
    
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
    
    
    for i=1:5e4
        ix=randi(length(segmentb),1);
        if segmentb(ix)+1000 < segmente(ix)-5000
            segment=randi([segmentb(ix)+1000 segmente(ix)-5000],1);
        else
            segment=NaN;
        end
        begin3=segment;
        end3=floor(begin3+5*samplerate);
        while ~isempty(intersect(unstable3,begin3:end3))
            segment=randi([segmentb(ix)+1000 segmente(ix)-5000],1);
            begin3=segment;
            end3=floor(begin3+5*samplerate);
        end
        
        if ~isnan(begin3)
            tremor_or2(i,1:(end3-begin3+1))=unwrap(phase(begin3:end3));
            tremor_or22(i,1:(end3-begin3+1))=(phase(begin3)+(0:1:(end3-begin3))*2*pi/(1000./mean(frequency(begin3-1000:begin3))));
            tremor_k(i,1)= (tremor_or2(i,(end3-begin3+1))-tremor_or22(i,(end3-begin3+1)))/(2*pi*0.001*(end3-begin3)); %mean(frequency(end3-1000:end3));%
        else
            tremor_or22(i,1:5001)=NaN;
            tremor_or2(i,1:5001)=NaN;
            tremor_k(i,1)=NaN;
        end
    end
    
    %%%% old version to calculate frequency change
    %         if ~isnan(begin3)
    %             tremor_or2(i,1:(end3-begin3+1))=frequency(begin3:end3);
    %             tremor_or22(i,1:(end3-begin3+1))=mean(frequency(begin3-1000:begin3));
    %             tremor_k(i,1)=tremor_or2(i,(end3-begin3+1))-tremor_or22(i,(end3-begin3+1));
    %         else
    %             tremor_or22(i,1:5001)=NaN;
    %             tremor_or2(i,1:5001)=NaN;
    %             tremor_k(i,1)=NaN;
    %         end
    %
    %     end
    
    for i=1:1e6
        dum=tremor_k(randi(25e3,1,rep));
        dum2=dum;
        NS_i1(i)=nanmedian(dum2);
    end
    NS_i(numb,:)=NS_i1;
    
    for i=1:12
        dum=tremor_k(randi(5e4,1,rep));
        dum2=dum;
        NS(numb,i)=nanmedian(dum2);
    end
    
    clearvars -except NS C AA BB numb iii PC A1 B1 nostim stim NS NS_i
    
end
clearvars  -except NS NS_i
 cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
% cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
save 'F_group12'





% load stim
% load nostim
% d=nanmedian(stimout);
%
% upperthreshold=prctile(nostimout,99.7917); %bonferroni corrected for 12 comparisons
% lowerthreshold=prctile(nostimout,0.2083);
%
% find(d>=upperthreshold | d<=lowerthreshold)
% bar(d)
% hold on
% rr(1:12)=upperthreshold;
% rr2(1:12)=lowerthreshold;
% plot(rr,'r--')
% plot(rr2,'r--')
% xlim([0.5 12.5])
% box('off')
% xticklabels({'0' '30' '60' '90' '120' '150' '180' '210' '240' '270' '300' '330'})
%
%
% %
% % close all
% % figure(2)
% % fig=gcf;
% % fig.Color=[1 1 1];
% % bar(100.*nanmedian(tt))
% % hold on
% % stem((100.*tt)')
% % xticklabels({'0','30','60','90','120','150','180','210','240','270','300','330'})
% % ylim([(-max((max(tt)))-0.1).*100 (max(max(tt))+0.1).*100])
% % box('off')
% % title ('P1')
