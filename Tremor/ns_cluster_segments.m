clear all, close all
cohort = [ 2 3 4 5 8 10 11 13 16 17];


main=[1 1 3 1 3 3 3 3 1 1];
method = 'ward'; %% Cluster method - change it here
% Options: ward, average, complete, single, weighted.

for iii = 1:length(cohort)
    
    % load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(cohort(iii)),'_RS.mat'))
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(cohort(iii)),'_RS.mat'))
    
    
    in2 = 1; % analysing the "main tremor axis"
    
    if in2 == 1
        in = 3;
    elseif in2 == 2 % other axis 1
        in = 5;
    elseif in2 == 3 % other axis 2
        in = 6;
    end
    
    %
    data = SmrData.WvData;
    
    %%
    samplerateold = SmrData.SR;
    tremor = (data(in, :));
    
    %%%%%%%%%%%%%%%%%%%
    
    time = 0:1/samplerateold:(size(data, 2)-1)/samplerateold;
    
    % downsample
    
    ts = timeseries(tremor,0:(1/samplerateold):((size(data,2)-1) / samplerateold));
    ts1 = resample(ts,0:0.001:((size(data,2)-1) / samplerateold), 'linear');
    tremor2(1:size(ts1.data, 3)) = ts1.data;
    samplerate = 1000;
    
    % determine stimulation time points
    index = [];
    for i = 2:size(data, 2)-1
        if data(2,i-1)<2.5 && data(2,i)>2.5
            index = [index i];
        end
    end
    clear i
    
    % Find trigger
    indexes4 = [index(1) index(find(diff(index) ./ samplerateold > 0.95)+1)];
    indexes3 = [index(find(diff(index) ./ samplerateold > 0.95)) index(end)];
    
    dd2 = round(data(4, :)*100) ./ 100;
    for i = 1:length(indexes4)
        xx(i) = round(dd2(indexes4(i)) ./ 0.1); %#ok<*SAGROW>
    end
    clear i
    
    start = floor((indexes4 ./ samplerateold)*samplerate);%+addon;
    ending = floor((indexes3 ./ samplerateold)*samplerate);%+addon+addon_end;%floor(5*samplerate);
    
    % when patient's hand is up
    handup = [];
    for i = 1:length(start)
        handup = [handup start(i):ending(i)]; %#ok<*AGROW>
    end
    clear i
    handup = sort(handup,'ascend');
    
    
    % tremor characteristics
    [Pxx, F] = pwelch(tremor2(handup), samplerate, [], samplerate, samplerate);
    
    frange = F(3:10);
    Pxxrange = Pxx(3:10);
    
    Fpeak = frange(find(Pxxrange == max(Pxxrange))); %#ok<*FNDSB>
    
    if (Fpeak-2) >= 1
        [b, a] = butter(2, [(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
    else
        [b, a] = butter(2, [(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
    end
    
    tremor_or = filtfilt(b, a, tremor2)*10*9.81/0.5;
    dummy = hilbert(tremor_or);
    % phase=angle(dummy);
    frequency = (smooth((1000/(2*pi))*diff(unwrap(angle(dummy))), 500))';
    
    tremor = (data(3, :));% %score(:,1)';%
    ts = timeseries(tremor, 0:(1/samplerateold):((size(data, 2)-1)/samplerateold));
    ts1 = resample(ts, 0:0.001:((size(data,2)-1)/samplerateold), 'linear');
    tremorx(1:size(ts1.data, 3)) = ts1.data;
    tremor = (data(5, :));% %score(:,1)';%
    ts = timeseries(tremor, 0:(1/samplerateold):((size(data, 2)-1)/samplerateold));
    ts1 = resample(ts, 0:0.001:((size(data, 2)-1)/samplerateold), 'linear');
    tremory(1:size(ts1.data, 3)) = ts1.data;
    tremor = (data(6, :));% %score(:,1)';%
    ts = timeseries(tremor, 0:(1/samplerateold):((size(data, 2)-1)/samplerateold));
    ts1 = resample(ts, 0:0.001:((size(data, 2)-1)/samplerateold), 'linear');
    tremorz(1:size(ts1.data, 3)) = ts1.data;
    tremorxf = filtfilt(b, a, tremorx);
    tremoryf = filtfilt(b, a, tremory);
    tremorzf = filtfilt(b, a, tremorz);
    envelope = [abs(hilbert(tremorxf)); abs(hilbert(tremoryf)); abs(hilbert(tremorzf))];
    phase = [angle(hilbert(tremorxf)); angle(hilbert(tremoryf)); angle(hilbert(tremorzf))];
    
    new = find(data(2, :) > 4);
    difp = find((diff(new)) > 100000); % are you trying to threshold at 9.6 seconds?
    ep_1 = [new(difp) new(end)];
    sp_1 = [new(1) new(difp+1)];
    
    %%% input start all trial
    start_t = 1;
    sp = sp_1(1, start_t:end);
    ep = ep_1(1, start_t:end);
    
    for ik = 1:length(sp) %%find double start and end points in a stimulation run
        s = (find(([indexes4 >= sp(ik)] + [indexes4 <= ep(ik)]) == 2));
        e = (find(([indexes3 >= sp(ik)] + [indexes3 <= ep(ik)]) == 2));
        tks = (find(diff(xx(s)) == 0)) + 1;
        tke = (find(diff(xx(e)) == 0));
        
        indexes4(s(tks)) = NaN;
        indexes3(e(tke)) = NaN;
        xx(e(tke)) = NaN;
    end
    
    %%%% find runs with trigering issues (too few, too many pulses)
    th1 = (Fpeak * 5) ./ 2;
    th2 = (Fpeak * 5) + round((Fpeak * 5) ./ 2);
    for it = 1:length(indexes4)
        if numel(index(find(index == indexes4(it)):find(index == indexes3(it)))) >= th1 && numel(index(find(index == indexes4(it)):find(index == indexes3(it)))) <= th2
            indexes4(it) = indexes4(it);
            indexes3(it) = indexes3(it);
            xx(it) = xx(it);
        else
            indexes4(it) = NaN;
            indexes3(it) = NaN;
            xx(it) = NaN;
        end
    end
    
    %%%%%%%%%%%%%%%
    indexes4 = indexes4(~isnan(indexes4));
    indexes3 = indexes3(~isnan(indexes3));
    xx = xx(~isnan(xx));
    
    start1 = [];
    ending1 = [];
    xx1 = [];
    for il = 1:length(sp)
        start1 = [start1 indexes4(find(([indexes4 >= sp(il)] + [indexes4 <= ep(il)]) == 2))]; % intersect([1 2 3],[3 4 5])
        ending1 = [ending1 indexes3(find(([indexes3 >= sp(il)] + [indexes3 <= ep(il)]) == 2))];
        xx1 = [xx1 xx(find(([indexes4 >= sp(il)] + [indexes4 <= ep(il)]) == 2))];
    end
    
    clear start ending
    start{1, 1} = floor((start1 ./ samplerateold) * samplerate);%+addon;
    ending{1, 1} = floor((ending1 ./ samplerateold) * samplerate);%+addon+addon_end;%floor(5*samplerate);
    clear xx
    xx{1, 1} = xx1;
    
    %% PCA - Random Stim
    
    % close all
    hh = numel(start);
    for j = 1:length(start{hh, 1})
        x = [tremorxf(start{hh,1}(j):ending{hh,1}(j)); tremoryf(start{hh,1}(j):ending{hh,1}(j)); tremorzf(start{hh,1}(j):ending{hh,1}(j))];
        [pc, score, latent, tsquare] = pca(x');
        pc_trials(j, 1:3) = pc(1:3, 1);
    end
    
    %% Baseline
    clearvars -except iii cohort pxx_z_zpad fig fig2 pc_trials C_all C_RS C_NS PC A1 B1 nostimout_c1 nostim_c1 nostimout_c2 nostim_c2 main method
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Baseline/P0',num2str(cohort(iii)),'_baseline.mat'))
    %    load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Baseline\P0',num2str(cohort(iii)),'_baseline.mat'))
    rng('default') % set random seed for consistency
    
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
    
    [d,e] = butter(2, [0.5/(0.5*samplerate) ],'low'); %15
    C = (filtfilt(d, e, tre_3'));
    
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
    
    % tf_3 = filtfilt(b, a, tre_3')*10*9.81/0.5;
    % dummy = (hilbert(tf_3))';
    % envelope = abs(dummy);
    
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
%     baseline = [tremorxf; tremoryf; tremorzf];
    ztremor=[zscore(tremorx);zscore(tremory);zscore(tremorz)];
    baseline=[filtfilt(b, a, ztremor(1,:)); filtfilt(b, a, ztremor(2,:));filtfilt(b, a, ztremor(3,:))];
    
    before_ns
    
    segmentb=hu{iii,:};
    segmente=hd{iii,:};
    
    % Pre-allocating for speed
    baseline3 = NaN(3,5e4);
    for_cluster = NaN(3,5e4,5001);
    env_cluster = NaN(3,5e4);
    
    for j = 1:5e4
        ix=randi(length(segmentb),1);
        segment=randi([segmentb(ix)+1000 segmente(ix)-5000],1);
        begin3=segment;
        end3=floor(begin3+5*samplerate);
        for ax = 1:3
            %         begin_idx(ax,j) = begin3;
            %         end_idx(ax,j) = end3;
            %             baseline3(ax,j) = (mean(envelope(ax,end3-1000:end3))-mean(envelope(ax, begin3-1000:begin3)))./mean(envelope(ax, begin3-1000:begin3));
            for_cluster(ax,j,:) = baseline(ax, begin3:end3);
            env_cluster(ax,j) = mean(envelope(ax, begin3-1000:begin3));
        end
    end
    
    %
    for j = 1:5e4 % in pca, rows are observations and columns are variables
        for_pca = squeeze(for_cluster(:,j,:)); % should be 3 vs length(segments)
        [pc, score, latent, tsquare] = pca(for_pca');
        pc_trials_ns(j, 1:3) = pc(1:3, 1);
    end
    
    %
    Z = linkage([pc_trials_ns; pc_trials], method);
    c = cluster(Z, 'Maxclust', 2);
    %     c_ns = c(1:5e4); % cluster indices for random stim
    
    idx_b_c1 = []; idx_b_c2 = [];
    for i = 1:5e4 % taking the baseline cluster indices
        if c(i) == 1
            idx_b_c1 = [idx_b_c1 i];
        else % c(i) == 2
            idx_b_c2 = [idx_b_c2 i];
        end
    end
    
    chosen=[2;2;1;2;2;1;2;2;2;2]; %%%% from adapted_clusters_data; chosen cluster
    qq=eval(['idx_b_c' num2str(chosen(iii))]);
    dum=squeeze(for_cluster(main(iii),qq,:));
    
%     [Pxx,F] = pwelch(dum',[],[],[],samplerate);
%     c=mean(Pxx,2);
%     plot(F,c)
%     xlim([0 50])
    
% %     dum2=reshape(dum',1,(size(dum,1)*size(dum,2)));
%          qqc=mean(squeeze(for_cluster(main(iii),qq,:)),1);

% %     [Pxx,F] = pwelch(qqc,samplerate*2,[],[],samplerate);
%         [Pxx,F] = pwelch(qqc,[],[],[],samplerate);
%         plot(F,Pxx)
%         hold on


    
%     baseline_comp(iii,:)=c;
    dumdrns=padarray(dum,[0 2000],0,'both');
    [Pxx_ns,F]=pwelch(dumdrns',[],[],[],samplerate);
    Pxx_ns=(mean(Pxx_ns,2))';

    pxx_z_zpad(iii,:)=Pxx_ns;
    
    clearvars -except pxx_z_zpad  cohort iii C_all C_RS C_NS PC A1 B1 iii nostimout_c1 nostim_c1 nostimout_c2 nostim_c2 randstim_change method main
    
end