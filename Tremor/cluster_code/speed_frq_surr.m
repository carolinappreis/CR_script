%%% 5000 surrogate frequency change was calculated to be inseted into mod_c
%%% and get 1000000 distribution acording to the 5000 that were in the
%%% cluster chosen


clear all; close all
cohort = [ 2 3 4 5 8 10 11 13 16 17];

% nostim= NaN(length(cohort),3,1e6);
% nostim_f=NaN(length(cohort),3,1e6);
freq_bl= NaN(10,3,5e4);


%%% main axis is adjusted in random stim according to psd calculated from
%%% the non stim condition, and accelerometer axis are swaped on the
%%% amplifyier accordingly - main and ns_mat code for this change
main=[1 1 3 1 3 3 3 3 1 1];
ns_mat=[[1 2 3]; [1 2 3]; [3 2 1]; [1 2 3];[3 2 1]; [3 2 1]; [3 2 1]; [3 2 1]; [1 2 3]; [1 2 3]];


for iii = 1:size(cohort,2)
    clearvars -except iii cohort main method tt1_all ns_mat freq_bl amp_bbl bs_begin bs_end pc_trials change_bl x_all pc1_exp explained_rs
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
    
    for aa = 1:3
        [Pxx,F] = pwelch(tre_3(aa,:), samplerate, [], samplerate, samplerate);
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
    tremorxf = filtfilt(b,a,tremorx);
    tremoryf = filtfilt(b,a,tremory);
    tremorzf = filtfilt(b,a,tremorz);
    
    ztremor=[zscore(tremorx);zscore(tremory);zscore(tremorz)];
    envelope = [abs(hilbert(tremorxf)); abs(hilbert(tremoryf)); abs(hilbert(tremorzf))];
    baseline = [tremorxf; tremoryf; tremorzf];
    phase = [angle(hilbert(tremorxf)); angle(hilbert(tremoryf)); angle(hilbert(tremorzf))];
    for i=1:3
        freqi(i,:)=(smooth((1000/(2*pi))*diff(unwrap(phase(i,:))),500))';
    end
    
    load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/aux_out.mat','bs_end','bs_begin')
    
    mmm=NaN(3,5e4);
    parfor j = 1:length(bs_end)
        begin3=bs_begin(iii,j);
        end3=bs_end(iii,j);
        %             ix=randi(length(segmentb),1);
        %             segment=randi([round(segmentb(ix)+1000) round(segmente(ix)-5000)],1);
        %             begin3=segment;
        %             end3=floor(begin3+5*samplerate);
        for ax = 1:3
            [mm]=simp(ax,iii,end3,begin3,freqi,phase,ns_mat);
            mmm(ax,j)=mm;
        end
    end
    freq_bl(iii,:,:)=mmm;
    
    
    clearvars -except  cohort iii nostim tt1_all main ns_mat amp_bbl bs_begin bs_end  change_bl x_all pc1_exp freq_bl
end

clearvars -except tt1_all amp_bbl bs_begin bs_end change_bl x_all pc1_exp freq_bl
