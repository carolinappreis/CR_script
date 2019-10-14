clear all
iii=[1];

for numb=1;
    %     :length(iii);
    clearvars -except nostimout iii numb cc cond nostim
    DBS_Fpeak
    
    cond={'_HFS_PS.mat'};
    cc=1;
    
    segmentb=[1621 112001 205601 287501 366701 450401];
    segmente=[76121 179301 263501 350901 422701 517501];
    start{1,1}=segmentb(1:2:end);
    start{2,1}=segmentb(2:2:end);
    ending{1,1}=segmente(1:2:end);
    ending{2,1}=segmente(2:2:end);
    
    
    load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),cond{cc,1}))
    %    load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),cond{cc,1}));
    
    in2=1; % analysing the "main tremor axis"
    
    if in2==1
        in=3;
    elseif in2==2 % other axis 1
        in=6;
    elseif in2==3 % other axis 2
        in=7;
    end
    data=SmrData.WvData;
    samplerateold=SmrData.SR;
    tremor=(data(in,:));
    addon=92; addon_end=35;
    
    time=0:1/samplerateold:(size(data,2)-1)/samplerateold;
    
    %%% downsample
    
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
    [b,a]=butter(2,[0.8/(0.5*samplerate) ],'low'); %15
    % tremor_or=zscore(tremor_or);
    dummy=hilbert(tremor_or);
    envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
    phase=angle(dummy);
    frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
    tremor=(data(3,:));% %score(:,1)';%
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremorx(1:size(ts1.data,3))=ts1.data;
    filt_x=filtfilt(b,a,tremorx);
    tremor=(data(6,:));% %score(:,1)';%
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremory(1:size(ts1.data,3))=ts1.data;
    filt_y=filtfilt(b,a,tremory);
    tremor=(data(7,:));% %score(:,1)';%
    ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    tremorz(1:size(ts1.data,3))=ts1.data;
    filt_z=filtfilt(b,a,tremorz);
    timeor=0:1/samplerate:(size(tremorx,2)-1)/samplerate;
    
    for ss=1:2
        for ix=1:length(start{ss,1});
            baseline3(1,ix)=mean(envelope(start{ss,1}(ix):ending{ss,1}(ix)));
        end
        HFS_median(ss,:)=median(baseline3);
        clear baseline3
    end
    
    
%     for i=1:1e6
%         dum=baseline3(randi(5e4,1,rep));
%         dum2=dum;
%         p(i)=nanmedian(dum2);
%     end
%     nostim(ss,numb,:)=p;
%     
%     for i=1:12
%         dum=baseline3(randi(5e4,1,rep));
%         dum2=dum;
%         nostimout(numb,ss,i)=nanmedian(dum2);
%     end
clearvars -except nostimout iii numb cc cond  nostim HFS_median
end

clearvars -except HFS_median
cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data')
% save ('HFS_median')



clear all
%  new=[1:8 10]; %%without PD patient (i.e., pt number 6)
 load ('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_amp_ARC.mat');load ('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\NS_PS_result.mat')
% load ('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/DBS_amp_ARC.mat'); load ('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/DBS_amp_NS.mat')

a.ns=NS; a.s=cell2mat(ttall); clearvars -except a


ref=a.s; %%% max amplitude change vs. max frequecy change
iii=1; %%%%% amp (=0) vs. supressive effect

idmi=[];
idma=[];
for i =1:size(ref,1)
    if (numel(find(ref(i,:)<0)))~=0 % all -except all amplifying subjects
        idmi=[idmi i];
    end
    if (numel(find(ref(i,:)>0)))~=0 % all -except all supressive subjects
        idma=[idma i];
    end
end

if iii==0;
    sub=ref(idma,:);  % iii=0 amplifying effect;
else
    sub=ref(idmi,:);  % iii~=0 supressive effect;
end

for i=1:size(sub,1);
    
    if iii==0
        phase_peak(1,i)=find(sub(i,:)==max(sub(i,:)));
    else
        phase_peak(1,i)=find(sub(i,:)==min(sub(i,:)));
    end
    
    if phase_peak(i)==1;
        a_s_al(i,:)=a.s(i,:);

        
    else
        a_s_al(i,:)=[a.s(i,phase_peak(i):end) a.s(i,1:phase_peak(i)-1)];
        check=[phase_peak(i):size(a_s_al,2) 1:phase_peak(i)-1];
    end
    
    if (phase_peak(i)+5)<= size(a_s_al,2)
        
        a_s2(i,:)= [a.s(i,phase_peak(i)+5:end) a.s(i,1:phase_peak(i)+5-1)];    
    else
        
        dum=(phase_peak(i)+5)-size(a_s_al,2);
        
        a_s2(i,:)= [a.s(i,dum+5:end) a.s(i,1:dum+5-1)];
        clear dum
        
    end
    
end


if kstest(a_s_al(:,1)-a_ns_al(:,1))==1
    
    for i=1:12
        [p1(1,i),h1(1,i)]=signrank(a_s_al(:,i),a_ns_al(:,i));
        a.s_ns=[h1(1) p1(1)];
        
        [p2(1,i),h2(1,i)]=signrank(a_s_al(:,i),a_s2(:,i));
        a.s_180=[h2(1) p2(1)];
        
        [p3(1,i),h3(1,i)]=signrank(a_ns_al(:,i),a_ns2(:,i));
        a.ns_180=[h3(1) p3(1)];
        
        
        
    end
    test_a='wilcoxon';
    
    
else
    
    [p1,h1]=ttest(a_s_al,a_ns_al);
    a.s_ns=[p1(1) h1(1)];
    
    [p2,h2]=ttest(a_s_al,a_s2);
    a.s_180=[p2(1) h2(1)];
    
    [p3,h3]=ttest(a_ns_al,a_ns2);
    a.ns_180=[p3(1) h3(1)];
    test_a='ttest';
end

