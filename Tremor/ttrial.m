clear all, close all
cohort = [ 2 3 4 5 8 10 11 13 16 17];

% nostim= NaN(length(cohort),3,1e6);
% nostim_f=NaN(length(cohort),3,1e6);
tt1_all=cell(10,3);
x_all=cell(10,1);
bs_begin=NaN(10,5e4);
bs_end=NaN(10,5e4);
amp_bbl = NaN(10,3,5e4);
change_bl = NaN(10,3,5e4);


main=[1 1 3 1 3 3 3 3 1 1];
ns_mat=[[1 2 3]; [1 2 3]; [3 2 1]; [1 2 3];[3 2 1]; [3 2 1]; [3 2 1]; [3 2 1]; [1 2 3]; [1 2 3]];


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
    
    hh=1;
    
    hh = numel(start);
    for j = 1:length(start{hh, 1})
        x = [tremorxf(start{hh,1}(j):ending{hh,1}(j)); tremoryf(start{hh,1}(j):ending{hh,1}(j)); tremorzf(start{hh,1}(j):ending{hh,1}(j))];
        [pc, score, latent, tsquare, explained] = pca(x');
        pc_trials(j, 1:3) = pc(1:3, 1);
        explained_rs(j, 1:3) = explained;
    end
    
    
    % % %     tremor_or2=NaN(length(start{hh,1}),1);
    % % %
    % % %     for axx=1:3
    % % %         for i=1:length(start{hh,1})
    % % %             if (~isnan(start{hh,1}(i)))
    % % %                 tremor_or2(axx,i,1)=(mean(envelope(axx,ending{hh,1}(i)-1000:ending{hh,1}(i)))-mean(envelope(axx,start{hh,1}(i)-1000:start{hh,1}(i))))/mean(envelope(axx,start{hh,1}(i)-1000:start{hh,1}(i)));
    % % %                 xx{hh,1}(i)= xx{hh,1}(i);
    % % %
    % % %             else
    % % %                 tremor_or2(axx,i,1)=NaN;
    % % %                 xx{hh,1}(i)= NaN;
    % % %             end
    % % %         end
    % % %
    % % %         %         %% criteria for outliers
    % % %         %
    % % %         %         idx_outl=find(tremor_or2>(nanmean(tremor_or2)+(2*(nanstd(tremor_or2))))|tremor_or2<(nanmean(tremor_or2)-(2*(nanstd(tremor_or2)))));
    % % %         %         tremor_or2(idx_outl,1)=NaN;
    % % %         %         tremor_or3(idx_outl,1)=NaN;
    % % %         %         xx(1,idx_outl)=NaN;
    % % %
    % % %         tt=NaN(20,12);
    % % %         yy=xx{hh,1}(:);
    % % %         tt_pc=NaN(20,12);
    % % %
    % % %         for i=1:12
    % % %             tt(1:sum(yy==i),i)=tremor_or2(axx,find(yy==i));
    % % %
    % % %
    % % %             tt(tt==0)=NaN;
    % % %
    % % %         end
    % % %         tt1_all{iii,axx}=tt;
    % % %
    % % %     end
    % % %
    %% Baseline
    clearvars -except iii cohort main method tt1_all ns_mat amp_bbl bs_begin bs_end pc_trials change_bl x_all
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
    
    before_ns
    
    segmentb=hu{iii,:};
    segmente=hd{iii,:};
    
    
    % Pre-allocating for speed
    for_cluster = NaN(3,5e4,5001);
    
    
    for j = 1:5e4
        ix=randi(length(segmentb),1);
        segment=randi([round(segmentb(ix)+1000) round(segmente(ix)-5000)],1);
        begin3=segment;
        end3=floor(begin3+5*samplerate);
        bs_begin(iii,j)=begin3;
        bs_end(iii,j)=end3 ;
        
        for ax = 1:3
            
            %             tremor_f2(j,1:(end3-begin3+1))=unwrap(phase(ns_mat(iii,ax),begin3:end3));
            %             tremor_f22(j,1:(end3-begin3+1))=(phase(ns_mat(iii,ax),begin3)+(0:1:(end3-begin3))*2*pi/(1000./mean(freqi(ns_mat(iii,ax),begin3-1000:begin3))));
            %             tremor_k(iii,ax,j)= (tremor_f2(j,(end3-begin3+1))-tremor_f22(j,(end3-begin3+1)))/(2*pi*0.001*(end3-begin3)); %mean(frequency(end3-1000:end3));%
            %             clear tremor_f2 tremor_f22
            %%% ----------for frquency seg_bl=tremor_k;
            change_bl(iii,ax,j)=(mean(envelope(ns_mat(iii,ax),end3-1000:end3))-mean(envelope(ns_mat(iii,ax), begin3-1000:begin3)))./mean(envelope(ns_mat(iii,ax), begin3-1000:begin3));
            amp_bbl(iii,ax,j)=mean(envelope(ns_mat(iii,ax), begin3-1000:begin3));
            for_cluster(ax,j,:) = baseline(ns_mat(iii,ax), begin3:end3); % 50000 segments of 5 sec of filtered data
        end
    end
    
    for j = 1:5e4 % in pca, rows are observations and columns are variables
        for_pca = squeeze(for_cluster(:,j,:)); % should be 3 vs length(segments)
        [pc, score, latent, tsquare, explained] = pca(for_pca');
        pc_trials_ns(j, 1:3) = pc(1:3, 1);
        explained_ns(j, 1:3) = explained;
    end
    
    %
    x_all{iii,1}=[pc_trials_ns; pc_trials];
    
    clearvars -except  cohort iii nostim tt1_all main ns_mat amp_bbl bs_begin bs_end  change_bl x_all
end





%% now main script above has baseline axis already matched
%output called old_5000 , was first saved without  matching non stim
%
% % clear all
% % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/old_5000.mat','seg_bl')
% % rng('default') % set random seed for consistency
% % ns_mat=[[1 2 3]; [1 2 3]; [3 2 1]; [1 2 3];[3 2 1]; [3 2 1]; [3 2 1]; [3 2 1]; [1 2 3]; [1 2 3]];
% %
% % for iii=1:10
% %     for ax = 1:3
% %         rep = 10;
% %         baseline3_temp = seg_bl(iii,ns_mat(iii,ax),:);
% %         seg_bl(iii,ax,:)=baseline3_temp;
% %         dum = baseline3_temp(randi(length(baseline3_temp), 1e6, rep));
% %         dum2 = dum;
% %         p = nanmedian(dum2,2);
% %
% %         nostim1(iii,ax,:) = p;
% %
% %         clear dum dum2 p baseline3_temp
% %     end
% % end
% %

%%%% median split

% % %
% % % clear all
% % % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/NS_all.mat','ns_ba','seg_bl');
% % % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_trials.mat');
% % %
% % %
% % % ns_a=NaN(10,2,1e6);
% % % for iii=1:10
% % %     seg=[]; ba=[];
% % %     if isnan (clust_win(iii))
% % %         seg=squeeze(seg_bl(iii,1,:));
% % %         ba=ns_ba(iii,:);
% % %     else
% % %         seg=squeeze(seg_bl(iii,1,clust_trials{iii,clust_win(iii)}));
% % %         ba=ns_ba(iii,clust_trials{iii,clust_win(iii)});
% % %     end
% % %
% % %     a_1=[];
% % %     a_2=[];
% % %     for i=1:length(seg)
% % %         if  ba(i)<median(ba)
% % %             a_1=[a_1 seg(i)];
% % %         else
% % %             a_2=[a_2 seg(i)];
% % %         end
% % %     end
% % %
% % %     rep = 10;
% % %     dum1 = a_1(randi(length(a_1), 1e6, rep));
% % %     dum2 = a_2(randi(length(a_2), 1e6, rep));
% % %     ns_a(iii,1,:) = nanmedian(dum1,2);
% % %     ns_a(iii,2,:) = nanmedian(dum2,2);
% % %     clear dum1 dum2
% % % end