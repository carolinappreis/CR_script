clear all, close all
cohort = [ 2 3 4 5 8 10 11 13 16 17];

load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_trials.mat')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/clusters_CR.mat')


main=[1 1 3 1 3 3 3 3 1 1];
ns_mat=[[1 2 3]; [1 2 3]; [3 2 1]; [1 2 3];[3 2 1]; [3 2 1]; [3 2 1]; [3 2 1]; [1 2 3]; [1 2 3]];

for iii = 3
    %     1:length(cohort)
    %%% Baseline
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(cohort(iii)),'_baseline.mat'))
    %    load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(cohort(iii)),'_baseline.mat'))
    rng('default') % set random seed for consistency
    gen=rng;
    
    in2 = 1; % analysing the "main tremor axis"
    
    if in2 == 1
        in = 3;
    elseif in2 == 2 % other axis 1
        in = 5;
    elseif in2 == 3 % other axis 2
        in = 6;
    end
    
    data = SmrData.WvData;
    samplerateold = SmrData.SR;
    tremor = (data(in, :));
    
    ts = timeseries(data, 0:(1 / samplerateold):((size(data, 2)-1) / samplerateold));
    ts1 = resample(ts, 0:0.001:((size(data, 2)-1)/samplerateold), 'linear');
    ds_data(1:size(ts1.data, 1), 1:size(ts1.data, 3)) = ts1.data;
    samplerate = 1000;
    tre_3 = ds_data([3 5 6],:); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    
    before_ns
    
    segmentb=round(hu{iii,:});
    segmente=round(hd{iii,:});
    
    % when patient's hand is up
    handup = [];
    for i = 1:length(segmentb)
        handup = [handup segmentb(i):segmente(i)]; %#ok<*AGROW>
    end
    clear i
    handup = sort(handup,'ascend');
    
    for aa = 1:3
        [Pxx,F] = pwelch(tre_3(aa,handup), samplerate, [], samplerate, samplerate);
        frange = F(3:10);
        Pxxrange = Pxx(3:10);
        Freqpeak(aa,:) = frange(find(Pxxrange == max(Pxxrange)));
        Ppeak(aa,:) = max(Pxxrange);
        ps_curves(aa,:) = Pxx;
    end
    
    peak_ax = [(Freqpeak(find(Ppeak == max(Ppeak)))) (find(Ppeak == max(Ppeak)))];
    Fpeak = peak_ax(1);
    
    if (Fpeak-2) >= 1
        [b,a] = butter(2, [(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
    else
        [b,a] = butter(2, [(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
    end
    
    tremor = (data(3, :));% %score(:,1)';%
    ts = timeseries(tremor, 0:(1/samplerateold):((size(data, 2)-1)/samplerateold));
    ts1 = resample(ts, 0:0.001:((size(data, 2)-1)/samplerateold), 'linear');
    tremorx(1:size(ts1.data, 3)) = ts1.data;
    tremor = (data(5, :));% %score(:,1)';%
    ts = timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1 = resample(ts,0:0.001:((size(data,2)-1)/samplerateold), 'linear');
    tremory(1:size(ts1.data,3)) = ts1.data;
    tremor = (data(6,:));% %score(:,1)';%
    ts = timeseries(tremor, 0:(1/samplerateold):((size(data, 2)-1)/samplerateold));
    ts1 = resample(ts, 0:0.001:((size(data, 2)-1)/samplerateold), 'linear');
    tremorz(1:size(ts1.data, 3)) = ts1.data;
    tremorxf = filtfilt(b, a, tremorx);
    tremoryf = filtfilt(b, a, tremory);
    tremorzf = filtfilt(b, a, tremorz);
    
    envelope = [abs(hilbert(tremorxf)); abs(hilbert(tremoryf)); abs(hilbert(tremorzf))];
    baseline = [tremorxf; tremoryf; tremorzf];
    
    ztremor=[zscore(tremorx);zscore(tremory);zscore(tremorz)];
    %     bline=[filtfilt(b, a, ztremor(1,:)); filtfilt(b, a, ztremor(2,:));filtfilt(b, a, ztremor(3,:))];
    
    % Pre-allocating for speed
    for_cluster = NaN(5e4,5001);
    idx_cseg= NaN(5e4,5001);
    
    for j = 1:5e4
        ix=randi(length(segmentb),1);
        segment=randi([round(segmentb(ix)+1000) round(segmente(ix)-5000)],1);
        begin3=segment;
        end3=floor(begin3+5*samplerate);
        idx_cseg(j,:)=begin3:end3;
        for_cluster(j,:) = ztremor(main(iii), begin3:end3);
    end    
    
    hp=(ztremor(main(iii),handup))';
     
    for gg=1:2
        hh=find(C_NS{iii,1}==gg);
        
        id1=idx_cseg(hh,:)';
        val1=for_cluster(hh,:)';
        
        indi=id1(:);
        val=val1(:);
        
        [indi_c, ia, ic] = unique(indi);
        
        clu{gg,1}=indi_c;
        cia{gg,1}=ia;
        
        clear hh  id1 val1 indi cal indi_c ia ic
    end   
end

p=[clu{1,1} ; clu{2,1}];
g1=clu{1,1}';
g2=clu{2,1}';


num=[];
for i=1:length(g1)
    if ~isempty(find(g2==g1(i)))
       num=[num i];
    end
end

bins=1:5001:length(indi);

for i=1:length(bins)
    if bins(i+1)<length(indi)
    in(i,1:5001)=indi(bins(i):bins(i+1)-1);
    end
end
    
plot(indi,val)
hold on
plot(indi(cia{2,1}(num)),val(cia{1,1}(num)),'.')

mr=sort(cia{1,1}(num),'ascend');
dum=find(diff(cia{1,1}(num))>1);

[ii,jj]=ismember(mr,id1(1,:)); 

C=id1(jj(ii),:)



dum=g1(num);
time=1:size(ztremor,2);
figure()
plot(time,ztremor(1,:))
hold on
plot(time(dum),ztremor(1,dum))

subplot(3,1,1)
plot(time,ztremor(1,:))
hold on
plot(time(g1),ztremor(1,g1))
subplot(3,1,2)
plot(time,ztremor(1,:))
hold on
plot(time(g2),ztremor(1,g2))
subplot(3,1,3)
plot(time,ztremor(1,:))
hold on
plot(time(g1),ztremor(1,g1))
plot(time(id1(:,guardo)),val1(:,guardo),'y');


    
    


