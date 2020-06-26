clear all, close all
rng('default') % set random seed for consistency

% cohort = [1 2 3 4 5 6 7 8 10 11 12 14]; 
% randstim_change = [1 1 2 1 2 1 2 2 1 2 1 2];
cohort = [1 4 5 6 7 8 10 11 12 14];
randstim_change = [1 1 2 1 2 2 1 2 1 2];
main = [1 1 3 1 3 3 1 3 1 3];

method = 'ward'; %% Cluster method - change it here
% Options: ward, average, complete, single, weighted.

pc_idx = 1; %% 1st, 2nd or 3rd principal component - change it here
% Options: 1, 2, 3.

evaluate = 0; %% If evaluate = 1, do silhoutte
test_frequencies = 0;

shuffle_data = 1;
frequency_amplitude = 1; group_level = 0; group_level_full = 1;
power_amplitude = 1;

test_baseline = 0;
variance_explained = 0; 

tuning_curve = 0;

if group_level_full == 1
    cohort = [1 4 5 6 7 8 10 11 12 14 7 14];
    randstim_change = [1 1 2 1 2 2 1 2 1 2 1 2];
end

if evaluate == 1 && pc_idx == 2
    cohort = [7 10 14];
    randstim_change = [2 1 2];
end

if tuning_curve == 1 && pc_idx == 2
    cohort = [7 10 14];
    randstim_change = [2 1 2];
end

for iii = 1:length(cohort)
% cd('C:\Users\beatrizsda\Desktop\Peripheral Stim\Analysis\Cursors\B')
% SMR_File_To_Mat;
load(strcat('/Users/beatrizarruda/OneDrive - Nexus365/Peripheral Stim/Analysis/Mat/RS/P0',num2str(cohort(iii)),'_random_stim_cursors.mat'))
% load(strcat('C:\Users\beatrizsda\OneDrive - Nexus365\Peripheral Stim\Analysis\Mat\RS\P0',num2str(cohort(iii)),'_random_stim_cursors.mat'))

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
% addon=92; addon_end=35;

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

% tremor_or = filtfilt(b, a, tremor2)*10*9.81/0.5;
% dummy = hilbert(tremor_or);
% % phase=angle(dummy);
% frequency = (smooth((1000/(2*pi))*diff(unwrap(angle(dummy))), 500))';

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
% phase = [angle(hilbert(tremorxf)); angle(hilbert(tremoryf)); angle(hilbert(tremorzf))];

% Vector
vec = [tremorxf; tremoryf; tremorzf];
for i = 1:size(vec, 1)
    min_val = abs(min(vec(i, :)));
    vec(i, :) = vec(i, :) + min_val + 0.01;
end

for j = 1:size(vec, 2)
    vector(j, 1) = sqrt(vec(1, j)^2 + vec(2, j)^2 + vec(3, j)^2);
end
% vectorf = filtfilt(b, a, vector);
[ahp, bhp] = butter(2, 2/(0.5*samplerate), 'high');
vectorf = filtfilt(ahp, bhp, vector)';

envelope = [abs(hilbert(tremorxf)); abs(hilbert(tremoryf));...
    abs(hilbert(tremorzf)); abs(hilbert(vectorf))];

% Frequency amplitude envelope
freq_amp = [(smooth((1000/(2*pi))*diff(unwrap(angle(hilbert(tremorxf*10*9.81/0.5)))), 500))';...
    (smooth((1000/(2*pi))*diff(unwrap(angle(hilbert(tremoryf*10*9.81/0.5)))), 500))';...
    (smooth((1000/(2*pi))*diff(unwrap(angle(hilbert(tremorzf*10*9.81/0.5)))), 500))';
    (smooth((1000/(2*pi))*diff(unwrap(angle(hilbert(vectorf*10*9.81/0.5)))), 500))'];

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
th1 = (Fpeak * 5 * 5) ./ 2;
th2 = (Fpeak * 5 * 5) + round((Fpeak * 5 * 5) ./ 2);
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

%% PCA attempt


if group_level_full == 1
    if iii == 7 || iii == 11 || iii == 12  % Patients 10, 7 (2nd time), 14 (2nd time)
        pc_idx = 2;
    else
        pc_idx = 1;
    end
end

% close all
hh = numel(start);
for j = 1:length(start{hh, 1})
    if randstim_change(iii) == 1
        x = [tremorxf(start{hh,1}(j):ending{hh,1}(j)); tremoryf(start{hh,1}(j):ending{hh,1}(j)); tremorzf(start{hh,1}(j):ending{hh,1}(j))];
    else % randstim_change(iii) == 2
        x = [tremorzf(start{hh,1}(j):ending{hh,1}(j)); tremoryf(start{hh,1}(j):ending{hh,1}(j)); tremorxf(start{hh,1}(j):ending{hh,1}(j))];
    end
    [pc, score, latent, tsquare, explained] = pca(x');
    pc_trials(j, 1:3) = pc(1:3, pc_idx);
    explained_rs(j, 1:3) = explained;
end

if variance_explained == 1
    mean_explained_rs(iii, 1) = mean(explained_rs(:, 1));
    mean_explained_rs(iii, 2) = mean(explained_rs(:, 2));
    mean_explained_rs(iii, 3) = mean(explained_rs(:, 3));
    
    median_explained_rs(iii, 1) = median(explained_rs(:, 1));
    median_explained_rs(iii, 2) = median(explained_rs(:, 2));
    median_explained_rs(iii, 3) = median(explained_rs(:, 3));
    
    std_explained_rs(iii, 1) = std(explained_rs(:, 1));
    std_explained_rs(iii, 2) = std(explained_rs(:, 2));
    std_explained_rs(iii, 3) = std(explained_rs(:, 3));
end
%% Frequency amplitude - Random Stim

if frequency_amplitude == 1
    
hh = numel(start);
for j = 1:length(start{hh, 1})
    if randstim_change(iii) == 1
        amp = [freq_amp(1, start{hh,1}(j):ending{hh,1}(j));...
            freq_amp(2, start{hh,1}(j):ending{hh,1}(j));...
            freq_amp(3, start{hh,1}(j):ending{hh,1}(j));...
            freq_amp(4, start{hh,1}(j):ending{hh,1}(j))];
    else % randstim_change(iii) == 2
        amp = [freq_amp(3, start{hh,1}(j):ending{hh,1}(j));...
            freq_amp(2, start{hh,1}(j):ending{hh,1}(j));...
            freq_amp(1, start{hh,1}(j):ending{hh,1}(j));...
            freq_amp(4, start{hh,1}(j):ending{hh,1}(j))];
    end
    mean_freq_amp_rs(1, j) = round(mean(amp(1,:)), 1);
    mean_freq_amp_rs(2, j) = round(mean(amp(2,:)), 1);
    mean_freq_amp_rs(3, j) = round(mean(amp(3,:)), 1);
    mean_freq_amp_rs(4, j) = round(mean(amp(4,:)), 1);
end
       
end

%% Frequencies - Random Stim

if test_frequencies == 1

hh = numel(start);
for j = 1:length(start{hh, 1})
    if randstim_change(iii) == 1
        x = tremorxf(start{hh,1}(j):ending{hh,1}(j));
        y = tremoryf(start{hh,1}(j):ending{hh,1}(j)); 
        z = tremorzf(start{hh,1}(j):ending{hh,1}(j));
    else % randstim_change(iii) == 2
        z = tremorxf(start{hh,1}(j):ending{hh,1}(j));
        y = tremoryf(start{hh,1}(j):ending{hh,1}(j)); 
        x = tremorzf(start{hh,1}(j):ending{hh,1}(j));
    end
  
    [Pxx_x, F_x] = pwelch(x, samplerate, [], samplerate, samplerate);
    frange_x = F_x(3:10);
    Pxxrange_x = Pxx_x(3:10);
    F_x_rs(j, :) = frange_x(find(Pxxrange_x == max(Pxxrange_x)));
    P_x_rs(j,:) = max(Pxxrange_x);
        
    [Pxx_y, F_y] = pwelch(y, samplerate, [], samplerate, samplerate);
	frange_y = F_y(3:10);
	Pxxrange_y = Pxx_y(3:10);
	F_y_rs(j, :) = frange_y(find(Pxxrange_y == max(Pxxrange_y)));
	P_y_rs(j,:) = max(Pxxrange_y);
        
	[Pxx_z, F_z] = pwelch(z, samplerate, [], samplerate, samplerate);
	frange_z = F_z(3:10);
	Pxxrange_z = Pxx_z(3:10);
	F_z_rs(j, :) = frange_z(find(Pxxrange_z == max(Pxxrange_z)));
	P_z_rs(j,:) = max(Pxxrange_z);
end

end

%% Power Amplitude

if power_amplitude == 1
    
    envelopex = envelope(1, :); 
    envelopey = envelope(2, :); 
    envelopez = envelope(3, :); 
    envelopev = envelope(4, :); 
    
hh = numel(start);
for j = 1:length(start{hh, 1})
    
    if randstim_change(iii) == 1
        x = envelopex(start{hh,1}(j):ending{hh,1}(j));
        y = envelopey(start{hh,1}(j):ending{hh,1}(j)); 
        z = envelopez(start{hh,1}(j):ending{hh,1}(j));
        v = envelopev(start{hh,1}(j):ending{hh,1}(j));
    else % randstim_change(iii) == 2
        z = envelopex(start{hh,1}(j):ending{hh,1}(j));
        y = envelopey(start{hh,1}(j):ending{hh,1}(j)); 
        x = envelopez(start{hh,1}(j):ending{hh,1}(j));
        v = envelopev(start{hh,1}(j):ending{hh,1}(j));
    end
  
    envelope_median_rs(j, 1) = median(x);
    envelope_median_rs(j, 2) = median(y);
    envelope_median_rs(j, 3) = median(z);
    envelope_median_rs(j, 4) = median(v);
end
       
end
%% Baseline
clearvars -except test_frequencies iii cohort fig fig2 pc_trials C_RS randstim_change method pc_idx evaluate f1 optk per_neg critval mean_coeff  F_x_rs F_y_rs F_z_rs P_x_rs P_y_rs P_z_rs frequency_amplitude test_baseline Fpeak mean_freq_amp_rs KS_p KS_h T_p T_h test_baseline power_amplitude envelope_median_rs variance_explained median_explained_ns mean_explained_ns std_explained_ns median_explained_rs mean_explained_rs std_explained_rs counts_ns counts_half1 counts_half2 median_freq SR_freq group_level RSum_p RSum_h SR_amp median_amp group_level_full KS_p_freq KS_h_freq T_p_freq T_h_freq RSum_p_freq RSum_h_freq KS_p_power KS_h_power T_p_power T_h_power RSum_p_power RSum_h_power tuning_curve amp_n_bins1 amp_n_bins2 bins main shuffle_data

load(strcat('/Users/beatrizarruda/OneDrive - Nexus365/Peripheral Stim/Analysis/Mat/B/P0',num2str(cohort(iii)),'_baseline_cursors.mat'))
% load(strcat('C:\Users\beatrizsda\OneDrive - Nexus365\Peripheral Stim\Analysis\Mat\B\P0',num2str(cohort(iii)),'_baseline_cursors.mat'))

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
tre_3 = ds_data([3 5 6], :); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    
[d, e] = butter(2, [0.5/(0.5*samplerate) ],'low'); %15
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
Fpeak_ns = peak_ax(1);
    
if (Fpeak_ns-2) >= 1
    [b,a] = butter(2, [(Fpeak_ns-2)/(0.5*samplerate) (Fpeak_ns+2)/(0.5*samplerate)], 'bandpass'); %15
else
    [b,a] = butter(2, [(1)/(0.5*samplerate) (Fpeak_ns+2)/(0.5*samplerate)], 'bandpass'); %15
end

% tf_3 = filtfilt(b, a, tre_3')*10*9.81/0.5;
% dummy = (hilbert(tf_3))';
% envelope = abs(dummy); 

if tuning_curve == 1
    tf_3 = filtfilt(b, a, tre_3') * 10 * 9.81 / 0.5;
    tf_3 = tf_3';
     
    for i = 1:3
        env(i, :) = abs(hilbert(tf_3(i, :)));
        phase(i, :) = angle(hilbert(tf_3(i, :)));
        freq(i, :) = (smooth((1000 / (2 * pi)) * diff(unwrap(phase(i, :))), 500))';
    end
    
    ns_amp = env(main(iii), :);
    ns_freq = freq(main(iii), :);
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
ts1 = resample(ts, 0:0.001:((size(data,   2)-1)/samplerateold), 'linear');
tremorz(1:size(ts1.data, 3)) = ts1.data;
tremorxf = filtfilt(b, a, tremorx);
tremoryf = filtfilt(b, a, tremory);
tremorzf = filtfilt(b, a, tremorz);

% Vector
vec = [tremorxf; tremoryf; tremorzf];
for i = 1:size(vec, 1)
    min_val = abs(min(vec(i, :)));
    vec(i, :) = vec(i, :) + min_val + 0.01;
end

for j = 1:size(vec, 2)
    vector(j, 1) = sqrt(vec(1, j)^2 + vec(2, j)^2 + vec(3, j)^2);
end
% vectorf = filtfilt(b, a, vector);
[ahp, bhp] = butter(2, 2/(0.5*samplerate), 'high');
vectorf = filtfilt(ahp, bhp, vector)';


envelope = [abs(hilbert(tremorxf)); abs(hilbert(tremoryf));...
    abs(hilbert(tremorzf)); abs(hilbert(vectorf))];

% Frequency amplitude envelope
freq_amp = [(smooth((1000/(2*pi))*diff(unwrap(angle(hilbert(tremorxf*10*9.81/0.5)))), 500))';...
    (smooth((1000/(2*pi))*diff(unwrap(angle(hilbert(tremoryf*10*9.81/0.5)))), 500))';...
    (smooth((1000/(2*pi))*diff(unwrap(angle(hilbert(tremorzf*10*9.81/0.5)))), 500))';
    (smooth((1000/(2*pi))*diff(unwrap(angle(hilbert(vectorf*10*9.81/0.5)))), 500))'];

%% PCA - baseline
% Split in 5-sec segments for calculating PCA over time

L = 5000;
tremorxf = split_vector_length(tremorxf, L);
tremoryf = split_vector_length(tremoryf, L);
tremorzf = split_vector_length(tremorzf, L);

if tuning_curve == 1
    ns_amp_split = split_vector_length(ns_amp, L);
    ns_freq_split = split_vector_length(ns_freq, L);
end


if group_level_full == 1
    if iii == 7 || iii == 11 || iii == 12  % Patients 10, 7 (2nd time), 14 (2nd time)
        pc_idx = 2;
    else
        pc_idx = 1;
    end
end

for j = 1:size(tremorxf,1)
    x = [tremorxf(j,:); tremoryf(j,:); tremorzf(j,:)];
    [pc, score, latent, tsquare, explained] = pca(x');
    pc_trials_ns(j, 1:3) = pc(1:3, pc_idx);
    explained_ns(j, 1:3) = explained;
end

if variance_explained == 1
    mean_explained_ns(iii, 1) = mean(explained_ns(:, 1));
    mean_explained_ns(iii, 2) = mean(explained_ns(:, 2));
    mean_explained_ns(iii, 3) = mean(explained_ns(:, 3));
    
    median_explained_ns(iii, 1) = median(explained_ns(:, 1));
    median_explained_ns(iii, 2) = median(explained_ns(:, 2));
    median_explained_ns(iii, 3) = median(explained_ns(:, 3));
    
    std_explained_ns(iii, 1) = std(explained_ns(:, 1));
    std_explained_ns(iii, 2) = std(explained_ns(:, 2));
    std_explained_ns(iii, 3) = std(explained_ns(:, 3));
end
%% Frequency amplitude - Baseline

if frequency_amplitude == 1 || test_baseline == 1
    
    amp_x = split_vector_length(freq_amp(1, :), L);
    amp_y = split_vector_length(freq_amp(2, :), L);
    amp_z = split_vector_length(freq_amp(3, :), L);
    amp_v = split_vector_length(freq_amp(4, :), L);
    
    for j = 1:size(tremorxf,1)
        mean_freq_amp_ns(1, j) = round(mean(amp_x(j, :)), 1);
        mean_freq_amp_ns(2, j) = round(mean(amp_y(j, :)), 1);
        mean_freq_amp_ns(3, j) = round(mean(amp_z(j, :)), 1);
        mean_freq_amp_ns(4, j) = round(mean(amp_v(j, :)), 1);
    end
end

if test_baseline == 1
    metade = floor(size(mean_freq_amp_ns, 2)/2);
    
    half1 = ones(size(mean_freq_amp_ns(1, 1:metade)));
    half2 = 2*ones(size(mean_freq_amp_ns(1, metade+1:end)));
    
    baseline_halves = [half1, half2];
end

%% Frequencies - Baseline

for j = 1:size(tremorxf,1)
    [Pxx_x, F_x] = pwelch(tremorxf(j,:), samplerate, [], samplerate, samplerate);
    frange_x = F_x(3:10);
	Pxxrange_x = Pxx_x(3:10);
    F_x_ns(j, :) = frange_x(find(Pxxrange_x == max(Pxxrange_x)));
    P_x_ns(j,:) = max(Pxxrange_x);
    
    [Pxx_y, F_y] = pwelch(tremoryf(j,:), samplerate, [], samplerate, samplerate);
    frange_y = F_y(3:10);
	Pxxrange_y = Pxx_y(3:10);
    F_y_ns(j, :) = frange_y(find(Pxxrange_y == max(Pxxrange_y)));
    P_y_ns(j,:) = max(Pxxrange_y);
    
    [Pxx_z, F_z] = pwelch(tremorzf(j,:), samplerate, [], samplerate, samplerate);
    frange_z = F_z(3:10);
	Pxxrange_z = Pxx_z(3:10);
    F_z_ns(j, :) = frange_z(find(Pxxrange_z == max(Pxxrange_z)));
    P_z_ns(j,:) = max(Pxxrange_z);
end


%% Power Amplitude - Baseline

envelopex = split_vector_length(envelope(1, :), L);
envelopey = split_vector_length(envelope(2, :), L);
envelopez = split_vector_length(envelope(3, :), L);
envelopev = split_vector_length(envelope(4, :), L);

for j = 1:size(envelopex, 1)
    x = envelopex(j, :);
    y = envelopey(j, :);
    z = envelopez(j, :);
    v = envelopev(j, :);
    envelope_median_ns(j, 1) = median(x);
    envelope_median_ns(j, 2) = median(y);
    envelope_median_ns(j, 3) = median(z);
    envelope_median_ns(j, 4) = median(v);
end

%% Plot PCA

fig = figure(1)
set(gcf, 'color', 'w', 'Position', [80 100 2000 1000])

if iii <= 6
    subplot(6, 2, iii*2) % Odd plots - random stim
    plot(pc_trials(:,1), 'r.')
    hold on
    plot(pc_trials(:,2), 'b.')
    plot(pc_trials(:,3), 'k.')
    xlim([0, size(pc_trials,1)])
    if iii == 1
        title('Random Stim - Patient 1')
    else
        title(sprintf('Patient %d', cohort(iii)))
    end 
    
    if iii == 1 || iii == 7
        legend('CED2', 'CED5', 'CED6')
    end
    
    subplot(6, 2, iii*2-1) % Even plots - baseline
    plot(pc_trials_ns(:,1), 'r.')
    hold on
    plot(pc_trials_ns(:,2), 'b.')
    plot(pc_trials_ns(:,3), 'k.')
    xlim([0, size(pc_trials_ns, 1)])
    
    if iii == 1
        title('Baseline - Patient 1')
    else
        title(sprintf('Patient %d', cohort(iii)))
    end 
    
else % iii > 6
    fig2 = figure(2)
    set(gcf, 'color', 'w', 'Position', [80 100 2000 1000])

    subplot(6, 2, (iii-6)*2) % Odd plots - random stim
    plot(pc_trials(:,1), 'r.')
    hold on
    plot(pc_trials(:,2), 'b.')
    plot(pc_trials(:,3), 'k.')
    xlim([0, size(pc_trials,1)])
    if iii == 7
        title('Random Stim - Patient 7')
    else
        title(sprintf('Patient %d', cohort(iii)))
    end 
    
    if iii == 1 || iii == 7
        legend('CED2', 'CED5', 'CED6')
    end
    
    subplot(6, 2, (iii-6)*2-1) % Even plots - baseline
    plot(pc_trials_ns(:,1), 'r.')
    hold on
    plot(pc_trials_ns(:,2), 'b.')
    plot(pc_trials_ns(:,3), 'k.')
    xlim([0, size(pc_trials_ns,1)])
    
    if iii == 7
        title('Baseline - Patient 7')
    else
        title(sprintf('Patient %d', cohort(iii)))
    end 
    
end

%% Linkage

rng('default') % set random seed for consistency

figure(iii + 2)
set(gcf, 'color', 'w', 'Position', [300,300,1200,300])

Z = linkage([pc_trials_ns; pc_trials], method);
c = cluster(Z, 'Maxclust', 2);


% HERE FIX CLUSTER NOTATION TO MATCH RESAMPLE

if pc_idx == 1
    if cohort(iii) == 1 || cohort(iii) == 8 || cohort(iii) == 12
        c_temp = NaN(size(c));
        c_temp(find(c == 1)) = 2;
        c_temp(find(c == 2)) = 1;
        c = c_temp;
    end
end

if pc_idx == 2
    if cohort(iii) == 1 || cohort(iii) == 6 || cohort(iii) == 12 || cohort(iii) == 14
        c_temp = NaN(size(c));
        c_temp(find(c == 1)) = 2;
        c_temp(find(c == 2)) = 1;
        c = c_temp;
    end
end



c_ns = c(1:size(pc_trials_ns, 1));
c_rs = c(size(pc_trials_ns, 1)+1 : end);


if tuning_curve == 1
    ns_amp_c1 = reshape(ns_amp_split(find(c_ns == 1), :)', [], 1)';
    ns_amp_c2 = reshape(ns_amp_split(find(c_ns == 2), :)', [], 1)';
    
    ns_freq_c1 = reshape(ns_freq_split(find(c_ns == 1), :)', [], 1)';
    ns_freq_c2 = reshape(ns_freq_split(find(c_ns == 2), :)', [], 1)';
    
    bins = 2:0.5:9;
    
    for i = 1:length(bins)
        if i + 1 < length(bins)
            dum1 = find(ns_freq_c1 > bins(i) & ns_freq_c1 <= bins(i+1));

            if ~isempty(dum1) 
                m_n_amp_c1(1, i) = nanmedian(ns_amp_c1(dum1));
            else
                m_n_amp_c1(1, i) = NaN;
            end
        end
     end
    
     max_c1 = max(m_n_amp_c1(1,:));
     m_n_amp_c1 = m_n_amp_c1 ./ max_c1;
     amp_n_bins1(iii, :) = m_n_amp_c1;
     
     
     for i = 1:length(bins)
        if i + 1 < length(bins)
            dum2 = find(ns_freq_c2 > bins(i) & ns_freq_c2 <= bins(i+1));

            if ~isempty(dum2) 
                m_n_amp_c2(1, i) = nanmedian(ns_amp_c2(dum2));
            else
                m_n_amp_c2(1, i) = NaN;
            end
        end
     end
    
     max_c2 = max(m_n_amp_c2(1,:));
     m_n_amp_c2 = m_n_amp_c2 ./ max_c2;
     amp_n_bins2(iii, :) = m_n_amp_c2;
end

ns1 = length(find(c_ns == 1));
ns2 = length(find(c_ns == 2));

counts_ns(iii, 1) = ns1;
if ns1 ~= 0
    counts_ns(iii, 2) = ns1/length(c_ns);
else
    counts_ns(iii, 2) = 0;
end

counts_ns(iii, 3) = ns2;
if ns2 ~= 0
    counts_ns(iii, 4) = ns2/length(c_ns);
else
    counts_ns(iii, 4) = 0;
end

% Z_rs = linkage(pc_trials, 'ward');
% c_rs = cluster(Z_rs, 'Maxclust', 2);
% 
% Z_ns = linkage(pc_trials_ns, 'ward');
% c_ns = cluster(Z_ns, 'Maxclust', 2);


subplot(1,3,1)
scatter3([pc_trials_ns(:, 1); pc_trials(:, 1)], [pc_trials_ns(:, 2); pc_trials(:, 2)], [pc_trials_ns(:, 3); pc_trials(:, 3)], 10, c)
title(sprintf('Patient %d - Combined', cohort(iii)))

subplot(1,3,2)
scatter3(pc_trials_ns(:, 1), pc_trials_ns(:, 2), pc_trials_ns(:, 3), 10, c_ns)
title('Baseline')

subplot(1,3,3)
scatter3(pc_trials(:, 1), pc_trials(:, 2), pc_trials(:, 3), 10, c_rs)
title('Random Stim')

C_RS{iii,1} = c_rs;


%% Test Clusters

if evaluate == 1
    rng('default')
%     [f1, optk, per_neg, critval] = evaluate_clusters([pc_trials_ns; pc_trials], c);
    data = [pc_trials_ns; pc_trials];
    clust = c;
    f1 = figure(length(cohort) + 3)
    set(gcf, 'color','w','Position', [80 100 2000 1000])

    subplot(2, 5, iii)

    [p, h] = silhouette(data, clust);

    for i = 1:2
    	clus{i, 1} = find(clust == i);
        silh(1, i) = (numel(find((p(clus{i, 1})) < 0))) ./ length(clus{i, 1}(:));
        ave_val(1,i) = mean(p(clus{i,1}));
    end

    E1 = evalclusters(data, 'linkage', 'Silhouette', 'klist', 1:2) 
    % E=evalclusters(X,c,'Silhouette')
    optk(iii,1) = E1.OptimalK;
    critval(iii,:) = E1.CriterionValues;
    per_neg(iii,:) = silh;
    mean_coeff(iii,:) = ave_val;
    
    title(sprintf('Patient %d', cohort(iii)))
    xlim([-1 1])
    
    clearvars -except fig fig2 cohort iii C_RS randstim_change method pc_idx evaluate f1 optk per_neg critval mean_coeff test_frequencies frequency_amplitude test_baseline power_amplitude variance_explained median_explained_ns mean_explained_ns std_explained_ns median_explained_rs mean_explained_rs std_explained_rs counts_ns counts_half1 counts_half2 median_freq SR_freq group_level SR_amp median_amp group_level_full KS_p_freq KS_h_freq T_p_freq T_h_freq RSum_p_freq RSum_h_freq KS_p_power KS_h_power T_p_power T_h_power RSum_p_power RSum_h_power tuning_curve amp_n_bins1 amp_n_bins2 bins main shuffle_data


end

%% Frequencies
if test_frequencies == 1
figure(2*length(cohort) + iii), set(gcf, 'color', 'w', 'Position', [80 100 1200 500])

subplot(2,6,1)
histogram(F_x_ns(find(c_ns == 1)),5)
title('X1 - Baseline')
xlabel('Frequency [Hz]')
ylabel('Counts')


subplot(2,6,2)
histogram(F_y_ns(find(c_ns == 1)),5)
title('Y1 - Baseline')
xlabel('Frequency [Hz]')
ylabel('Counts')


subplot(2,6,3)
histogram(F_z_ns(find(c_ns == 1)),5)
title('Z1 - Baseline')
xlabel('Frequency [Hz]')
ylabel('Counts')


subplot(2,6,4)
histogram(F_x_rs(find(c_rs == 1)),5)
title('X1 - Random Stim')
xlabel('Frequency [Hz]')
ylabel('Counts')


subplot(2,6,5)
histogram(F_y_rs(find(c_rs == 1)),5)
title('Y1 - Random Stim')
xlabel('Frequency [Hz]')
ylabel('Counts')


subplot(2,6,6)
histogram(F_z_rs(find(c_rs == 1)),5)
title('Z1 - Random Stim')
xlabel('Frequency [Hz]')
ylabel('Counts')


subplot(2,6,7)
histogram(F_x_ns(find(c_ns == 2)),5)
title('X2 - Baseline')
xlabel('Frequency [Hz]')
ylabel('Counts')


subplot(2,6,8)
histogram(F_y_ns(find(c_ns == 2)),5)
title('Y2 - Baseline')
xlabel('Frequency [Hz]')
ylabel('Counts')


subplot(2,6,9)
histogram(F_z_ns(find(c_ns == 2)),5)
title('Z2 - Baseline')
xlabel('Frequency [Hz]')
ylabel('Counts')


subplot(2,6,10)
histogram(F_x_rs(find(c_rs == 2)),5)
title('X2 - Random Stim')
xlabel('Frequency [Hz]')
ylabel('Counts')


subplot(2,6,11)
histogram(F_y_rs(find(c_rs == 2)),5)
title('Y2 - Random Stim')
xlabel('Frequency [Hz]')
ylabel('Counts')


subplot(2,6,12)
histogram(F_z_rs(find(c_rs == 2)),5)
title('Z2 - Random Stim')
xlabel('Frequency [Hz]')
ylabel('Counts')
end

%% Frequency Amplitudes

if frequency_amplitude == 1

figure(3*length(cohort) + iii), set(gcf, 'color', 'w', 'Position', [80 100 1200 500])

edges_rs = (Fpeak-2:0.2:Fpeak+2); % Create 20 bins.
edges_ns = (Fpeak_ns-2:0.2:Fpeak_ns+2); % Create 20 bins.

xlim_rs = [Fpeak-2, Fpeak+2];
xlim_ns = [Fpeak_ns-2, Fpeak_ns+2];

subplot(2,6,1)

% Plot the histogram.
histogram(mean_freq_amp_ns(1, find(c_ns == 1)), 'BinEdges', edges_ns)
title('X1 - Baseline')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_ns)

subplot(2,6,2)
histogram(mean_freq_amp_ns(2, find(c_ns == 1)), 'BinEdges', edges_ns)
title('Y1 - Baseline')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_ns)

subplot(2,6,3)
histogram(mean_freq_amp_ns(3, find(c_ns == 1)), 'BinEdges', edges_ns)
title('Z1 - Baseline')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_ns)

subplot(2,6,4)
histogram(mean_freq_amp_rs(1, find(c_rs == 1)), 'BinEdges', edges_rs)
title('X1 - Random Stim')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)

subplot(2,6,5)
histogram(mean_freq_amp_rs(2, find(c_rs == 1)), 'BinEdges', edges_rs)
title('Y1 - Random Stim')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)

subplot(2,6,6)
histogram(mean_freq_amp_rs(3, find(c_rs == 1)), 'BinEdges', edges_rs)
title('Z1 - Random Stim')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)

subplot(2,6,7)
histogram(mean_freq_amp_ns(1, find(c_ns == 2)), 'BinEdges', edges_ns)
title('X2 - Baseline')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_ns)

subplot(2,6,8)
histogram(mean_freq_amp_ns(2, find(c_ns == 2)), 'BinEdges', edges_ns)
title('Y2 - Baseline')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_ns)

subplot(2,6,9)
histogram(mean_freq_amp_ns(3, find(c_ns == 2)), 'BinEdges', edges_ns)
title('Z2 - Baseline')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_ns)

subplot(2,6,10)
histogram(mean_freq_amp_rs(1, find(c_rs == 2)), 'BinEdges', edges_rs)
title('X2 - Random Stim')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)

subplot(2,6,11)
histogram(mean_freq_amp_rs(2, find(c_rs == 2)), 'BinEdges', edges_rs)
title('Y2 - Random Stim')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)

subplot(2,6,12)
histogram(mean_freq_amp_rs(3, find(c_rs == 2)), 'BinEdges', edges_rs)
title('Z2 - Random Stim')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)


%% Frequency Amplitudes - VECTOR


figure(4*length(cohort) + iii), set(gcf, 'color', 'w', 'Position', [80 100 500 450])

edges_rs = (Fpeak-2:0.2:Fpeak+2);
edges_ns = (Fpeak_ns-2:0.2:Fpeak_ns+2);

xlim_rs = [Fpeak-2, Fpeak+2];
xlim_ns = [Fpeak_ns-2, Fpeak_ns+2];

subplot(2, 2, 1)
% Plot the histogram.
histogram(mean_freq_amp_ns(4, find(c_ns == 1)), 'BinEdges', edges_ns)
title(sprintf('Patient %d - V1 Baseline', cohort(iii)))
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_ns)

subplot(2, 2, 2)
% Plot the histogram.
histogram(mean_freq_amp_rs(4, find(c_rs == 1)), 'BinEdges', edges_rs)
title(sprintf('Patient %d - V1 Random Stim', cohort(iii)))
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)

subplot(2, 2, 3)
% Plot the histogram.
histogram(mean_freq_amp_ns(4, find(c_ns == 2)), 'BinEdges', edges_ns)
title('V2 Baseline')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_ns)

subplot(2, 2, 4)
% Plot the histogram.
histogram(mean_freq_amp_rs(4, find(c_rs == 2)), 'BinEdges', edges_rs)
title('V2 Random Stim')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)


%% Frequency amplitude - merged baseline and random stim

mean_freq_amp = [mean_freq_amp_ns, mean_freq_amp_rs];
freq_c1 = mean_freq_amp(:, find(c == 1));
freq_c2 = mean_freq_amp(:, find(c == 2));


figure(5*length(cohort) + iii), set(gcf, 'color', 'w', 'Position', [80 100 1000 500])

xlim_rs = [Fpeak-2, Fpeak+2];
xlim_ns = [Fpeak_ns-2, Fpeak_ns+2];

subplot(2,4,1)
% Plot the histogram.
histogram(mean_freq_amp(1, find(c == 1)), 'BinEdges', edges_rs)
title('X1')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)

subplot(2,4,2)
histogram(mean_freq_amp(2, find(c == 1)), 'BinEdges', edges_rs)
title('Y1')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)

subplot(2,4,3)
histogram(mean_freq_amp(3, find(c == 1)), 'BinEdges', edges_rs)
title('Z1')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)

subplot(2,4,4)
histogram(mean_freq_amp(4, find(c == 1)), 'BinEdges', edges_rs)
title('V1')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)

subplot(2,4,5)
histogram(mean_freq_amp(1, find(c == 2)), 'BinEdges', edges_rs)
title('X2')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)

subplot(2,4,6)
histogram(mean_freq_amp(2, find(c == 2)), 'BinEdges', edges_rs)
title('Y2')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)

subplot(2,4,7)
histogram(mean_freq_amp(3, find(c == 2)), 'BinEdges', edges_rs)
title('Z2')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)

subplot(2,4,8)
histogram(mean_freq_amp(4, find(c == 2)), 'BinEdges', edges_rs)
title('V2')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_rs)

%% Statistics - Two-sample Kolmogorov-Smirnov test; two-sample t-test


zeros_vector = zeros(size(freq_c1, 2));


[ks_h_x, ks_p_x] = kstest2(freq_c1(1, :), freq_c2(1, :));
[ks_h_y, ks_p_y] = kstest2(freq_c1(2, :), freq_c2(2, :));
[ks_h_z, ks_p_z] = kstest2(freq_c1(3, :), freq_c2(3, :));
[ks_h_v, ks_p_v] = kstest2(freq_c1(4, :), freq_c2(4, :));

[t_h_x, t_p_x] = ttest2(freq_c1(1, :), freq_c2(1, :));
[t_h_y, t_p_y] = ttest2(freq_c1(2, :), freq_c2(2, :));
[t_h_z, t_p_z] = ttest2(freq_c1(3, :), freq_c2(3, :));
[t_h_v, t_p_v] = ttest2(freq_c1(4, :), freq_c2(4, :));

[rsum_p_x, rsum_h_x] = ranksum(freq_c1(1, :), freq_c2(1, :));
[rsum_p_y, rsum_h_y] = ranksum(freq_c1(2, :), freq_c2(2, :));
[rsum_p_z, rsum_h_z] = ranksum(freq_c1(3, :), freq_c2(3, :));
[rsum_p_v, rsum_h_v] = ranksum(freq_c1(4, :), freq_c2(4, :));

KS_h_freq(iii, 1) = ks_h_x;
KS_h_freq(iii, 2) = ks_h_y;
KS_h_freq(iii, 3) = ks_h_z;
KS_h_freq(iii, 4) = ks_h_v;

KS_p_freq(iii, 1) = ks_p_x;
KS_p_freq(iii, 2) = ks_p_y;
KS_p_freq(iii, 3) = ks_p_z;
KS_p_freq(iii, 4) = ks_p_v;

T_h_freq(iii, 1) = t_h_x;
T_h_freq(iii, 2) = t_h_y;
T_h_freq(iii, 3) = t_h_z;
T_h_freq(iii, 4) = t_h_v;

T_p_freq(iii, 1) = t_p_x;
T_p_freq(iii, 2) = t_p_y;
T_p_freq(iii, 3) = t_p_z;
T_p_freq(iii, 4) = t_p_v;

RSum_p_freq(iii, 1) = rsum_p_x;
RSum_p_freq(iii, 2) = rsum_p_y;
RSum_p_freq(iii, 3) = rsum_p_z;
RSum_p_freq(iii, 4) = rsum_p_v;

RSum_h_freq(iii, 1) = rsum_h_x;
RSum_h_freq(iii, 2) = rsum_h_y;
RSum_h_freq(iii, 3) = rsum_h_z;
RSum_h_freq(iii, 4) = rsum_h_v; 
end


%% Test baseline

if test_baseline == 1
    half1_c1 = []; half2_c1 = []; half1_c2 = []; half2_c2 = [];
    
    for i = 1:length(c_ns)
        
        if c_ns(i) == 1 && baseline_halves(i) == 1
            half1_c1 = [half1_c1, i];
        end
        
        if c_ns(i) == 2 && baseline_halves(i) == 1
            half1_c2 = [half1_c2, i];
        end
        
        if c_ns(i) == 1 && baseline_halves(i) == 2
            half2_c1 = [half2_c1, i];
        end
        
        if c_ns(i) == 2 && baseline_halves(i) == 2
            half2_c2 = [half2_c2, i];
        end
        
    end
    
counts_h1c1 = length(half1_c1);
counts_h1c2 = length(half1_c2);
counts_h2c1 = length(half2_c1);
counts_h2c2 = length(half2_c2);

counts_half1(iii, 1) = counts_h1c1;
if counts_h1c1 ~= 0
    counts_half1(iii, 2) = counts_h1c1/length(find(baseline_halves == 1));
else
    counts_half1(iii, 2) = 0;
end

counts_half1(iii, 3) = counts_h1c2;
if counts_h1c2 ~= 0
    counts_half1(iii, 4) = counts_h1c2/length(find(baseline_halves == 1));
else
    counts_half1(iii, 4) = 0;
end

counts_half2(iii, 1) = counts_h2c1;
if counts_h2c1 ~= 0
    counts_half2(iii, 2) = counts_h2c1/length(find(baseline_halves == 2));
else
    counts_half2(iii, 2) = 0;
end

counts_half2(iii, 3) = counts_h2c2;
if counts_h2c2 ~= 0
    counts_half2(iii, 4) = counts_h2c2/length(find(baseline_halves == 2));
else
    counts_half2(iii, 4) = 0;
end

figure(5*length(cohort) + iii), set(gcf, 'color', 'w', 'Position', [80 100 500 500])

edges_ns = (Fpeak_ns-2:0.2:Fpeak_ns+2); % Create bin edges
xlim_ns = [Fpeak_ns-2, Fpeak_ns+2];

subplot(2,2,1)
% Plot the histogram.
histogram(mean_freq_amp_ns(4, half1_c1), 'BinEdges', edges_ns)
title('Half 1, Cluster 1')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_ns)

subplot(2,2,2)
% Plot the histogram.
histogram(mean_freq_amp_ns(4, half1_c2), 'BinEdges', edges_ns)
title('Half 1, Cluster 2')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_ns)

subplot(2,2,3)
% Plot the histogram.
histogram(mean_freq_amp_ns(4, half2_c1), 'BinEdges', edges_ns)
title('Half 2, Cluster 1')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_ns)

subplot(2,2,4)
% Plot the histogram.
histogram(mean_freq_amp_ns(4, half2_c2), 'BinEdges', edges_ns)
title('Half 2, Cluster 2')
xlabel('Frequency [Hz]')
ylabel('Counts')
xlim(xlim_ns)


%% Baseline test - Statistics - Two-sample Kolmogorov-Smirnov test; two-sample t-test

if ~isempty(mean_freq_amp_ns(4, half1_c1)) && ~isempty(mean_freq_amp_ns(4, half2_c1))
    [ks_h_h1c1_h2c1, ks_p_h1c1_h2c1] = kstest2(mean_freq_amp_ns(4, half1_c1), mean_freq_amp_ns(4, half2_c1));
    [t_h_h1c1_h2c1, t_p_h1c1_h2c1] = ttest2(mean_freq_amp_ns(4, half1_c1), mean_freq_amp_ns(4, half2_c1));
    KS_h(iii, 1) = ks_h_h1c1_h2c1;
    KS_p(iii, 1) = ks_p_h1c1_h2c1;
    T_h(iii, 1) = t_h_h1c1_h2c1;
    T_p(iii, 1) = t_p_h1c1_h2c1;
end

if ~isempty(mean_freq_amp_ns(4, half1_c2)) && ~isempty(mean_freq_amp_ns(4, half2_c2))
    [ks_h_h1c2_h2c2, ks_p_h1c2_h2c2] = kstest2(mean_freq_amp_ns(4, half1_c2), mean_freq_amp_ns(4, half2_c2));
    [t_h_h1c2_h2c2, t_p_h1c2_h2c2] = ttest2(mean_freq_amp_ns(4, half1_c2), mean_freq_amp_ns(4, half2_c2));
    KS_h(iii, 2) = ks_h_h1c2_h2c2;
    KS_p(iii, 2) = ks_p_h1c2_h2c2;
    T_h(iii, 2) = t_h_h1c2_h2c2;
    T_p(iii, 2) = t_p_h1c2_h2c2;
end

if ~isempty(mean_freq_amp_ns(4, half1_c1)) && ~isempty(mean_freq_amp_ns(4, half1_c2))
    [ks_h_h1c1_h1c2, ks_p_h1c1_h1c2] = kstest2(mean_freq_amp_ns(4, half1_c1), mean_freq_amp_ns(4, half1_c2));
    [t_h_h1c1_h1c2, t_p_h1c1_h1c2] = ttest2(mean_freq_amp_ns(4, half1_c1), mean_freq_amp_ns(4, half1_c2));
    KS_h(iii, 3) = ks_h_h1c1_h1c2;
    KS_p(iii, 3) = ks_p_h1c1_h1c2;
    T_h(iii, 3) = t_h_h1c1_h1c2;
    T_p(iii, 3) = t_p_h1c1_h1c2;
end

if ~isempty(mean_freq_amp_ns(4, half2_c1)) && ~isempty(mean_freq_amp_ns(4, half2_c2))
    [ks_h_h2c1_h2c2, ks_p_h2c1_h2c2] = kstest2(mean_freq_amp_ns(4, half2_c1), mean_freq_amp_ns(4, half2_c2));
    [t_h_h2c1_h2c2, t_p_h2c1_h2c2] = ttest2(mean_freq_amp_ns(4, half2_c1), mean_freq_amp_ns(4, half2_c2));
    KS_h(iii, 4) = ks_h_h2c1_h2c2;
    KS_p(iii, 4) = ks_p_h2c1_h2c2;
    T_h(iii, 4) = t_h_h2c1_h2c2;
    T_p(iii, 4) = t_p_h2c1_h2c2;
end

end

%% Power amplitude
if power_amplitude == 1
    
power_merged = [envelope_median_ns; envelope_median_rs]';
power_c1 = power_merged(:, find(c == 1));
power_c2 = power_merged(:, find(c == 2));


figure(7*length(cohort) + iii), set(gcf, 'color', 'w', 'Position', [80 100 1000 500])

% xlim_rs = [Fpeak-2, Fpeak+2];
% xlim_ns = [Fpeak_ns-2, Fpeak_ns+2];
nbins = 10;

subplot(2,4,1)
% Plot the histogram.
histogram(power_merged(1, find(c == 1)), nbins)
% histogram(mean_freq_amp(2, find(c == 1)), 'BinEdges', edges_rs)
title('X1')
xlabel('Power')
ylabel('Counts')
% xlim(xlim_rs)

subplot(2,4,2)
histogram(power_merged(2, find(c == 1)), nbins)
title('Y1')
xlabel('Power')
ylabel('Counts')
% xlim(xlim_rs)

subplot(2,4,3)
histogram(power_merged(3, find(c == 1)), nbins)
title('Z1')
xlabel('Power')
ylabel('Counts')
% xlim(xlim_rs)

subplot(2,4,4)
histogram(power_merged(4, find(c == 1)), nbins)
title('V1')
xlabel('Power')
ylabel('Counts')
% xlim(xlim_rs)

subplot(2,4,5)
histogram(power_merged(1, find(c == 2)), nbins)
title('X2')
xlabel('Power')
ylabel('Counts')
% xlim(xlim_rs)

subplot(2,4,6)
histogram(power_merged(2, find(c == 2)), nbins)
title('Y2')
xlabel('Power')
ylabel('Counts')
% xlim(xlim_rs)

subplot(2,4,7)
histogram(power_merged(3, find(c == 2)), nbins)
title('Z2')
xlabel('Power')
ylabel('Counts')
% xlim(xlim_rs)

subplot(2,4,8)
histogram(power_merged(4, find(c == 2)), nbins)
title('V2')
xlabel('Power')
ylabel('Counts')
% xlim(xlim_rs)

%% Statistics - Two-sample Kolmogorov-Smirnov test; two-sample t-test; Wilcoxon rank sum test

zeros_vector = zeros(size(power_c1, 2));

[ks_h_x, ks_p_x] = kstest2(power_c1(1, :), power_c2(1, :));
[ks_h_y, ks_p_y] = kstest2(power_c1(2, :), power_c2(2, :));
[ks_h_z, ks_p_z] = kstest2(power_c1(3, :), power_c2(3, :));
[ks_h_v, ks_p_v] = kstest2(power_c1(4, :), power_c2(4, :));

[t_h_x, t_p_x] = ttest2(power_c1(1, :), power_c2(1, :));
[t_h_y, t_p_y] = ttest2(power_c1(2, :), power_c2(2, :));
[t_h_z, t_p_z] = ttest2(power_c1(3, :), power_c2(3, :));
[t_h_v, t_p_v] = ttest2(power_c1(4, :), power_c2(4, :));

% 
% [rsum_p_x, rsum_h_x] = ranksum(abs(power_c1(1, :) - power_c2(1, :)), zeros_vector);
% [rsum_p_y, rsum_h_y] = ranksum(abs(power_c1(2, :) - power_c2(2, :)), zeros_vector);
% [rsum_p_z, rsum_h_z] = ranksum(abs(power_c1(3, :) - power_c2(3, :)), zeros_vector);
% [rsum_p_v, rsum_h_v] = ranksum(abs(power_c1(4, :) - power_c2(4, :)), zeros_vector);

KS_h_power(iii, 1) = ks_h_x;
KS_h_power(iii, 2) = ks_h_y;
KS_h_power(iii, 3) = ks_h_z;
KS_h_power(iii, 4) = ks_h_v;

KS_p_power(iii, 1) = ks_p_x;
KS_p_power(iii, 2) = ks_p_y;
KS_p_power(iii, 3) = ks_p_z;
KS_p_power(iii, 4) = ks_p_v;

T_h_power(iii, 1) = t_h_x;
T_h_power(iii, 2) = t_h_y;
T_h_power(iii, 3) = t_h_z;
T_h_power(iii, 4) = t_h_v;

T_p_power(iii, 1) = t_p_x;
T_p_power(iii, 2) = t_p_y;
T_p_power(iii, 3) = t_p_z;
T_p_power(iii, 4) = t_p_v;

RSum_p_power(iii, 1) = rsum_p_x;
RSum_p_power(iii, 2) = rsum_p_y;
RSum_p_power(iii, 3) = rsum_p_z;
RSum_p_power(iii, 4) = rsum_p_v;

RSum_h_power(iii, 1) = rsum_h_x;
RSum_h_power(iii, 2) = rsum_h_y;
RSum_h_power(iii, 3) = rsum_h_z;
RSum_h_power(iii, 4) = rsum_h_v; 
end

%% Reshuffle data

if shuffle_data == 1
    
    nsurrogates = 50000;
    c1_length = length(find(c == 1));
    c2_length = length(find(c == 2));
    envelope_median = [envelope_median_ns', envelope_median_rs'];
    % mean_freq_amp = mean_freq_amp;
    
    if iii == 1 % pre-allocate for speed
        SR_freq = NaN(nsurrogates, 4);
        SR_amp = NaN(nsurrogates, 4);
        median_freq = NaN(nsurrogates, length(cohort), 2, 4);
        median_amp = NaN(nsurrogates, length(cohort), 2, 4);
    end
    
    
    for j = 1:nsurrogates
        idx_all = 1:length(c);
        idx_c1 = randperm(length(c), c1_length); %randperm(n,k) returns a row vector containing k unique integers selected randomly from 1 to n.
        idx_all(idx_c1) = [];
        idx_c2 = idx_all;
        envelope_median_c1 = envelope_median(:, idx_c1);
        envelope_median_c2 = envelope_median(:, idx_c2);
        mean_freq_amp_c1 = mean_freq_amp(:, idx_c1);
        mean_freq_amp_c2 = mean_freq_amp(:, idx_c2);
        
        if group_level == 1 || group_level_full == 1
            median_freq(j, iii, 1, 1) = median(mean_freq_amp_c1(1, :));
            median_freq(j, iii, 2, 1) = median(mean_freq_amp_c2(1, :));
            median_freq(j, iii, 1, 2) = median(mean_freq_amp_c1(2, :));
            median_freq(j, iii, 2, 2) = median(mean_freq_amp_c2(2, :));
            median_freq(j, iii, 1, 3) = median(mean_freq_amp_c1(3, :));
            median_freq(j, iii, 2, 3) = median(mean_freq_amp_c2(3, :));
            median_freq(j, iii, 1, 4) = median(mean_freq_amp_c1(4, :));
            median_freq(j, iii, 2, 4) = median(mean_freq_amp_c2(4, :));
            
            median_amp(j, iii, 1, 1) = median(envelope_median_c1(1, :));
            median_amp(j, iii, 2, 1) = median(envelope_median_c2(1, :));
            median_amp(j, iii, 1, 2) = median(envelope_median_c1(2, :));
            median_amp(j, iii, 2, 2) = median(envelope_median_c2(2, :));
            median_amp(j, iii, 1, 3) = median(envelope_median_c1(3, :));
            median_amp(j, iii, 2, 3) = median(envelope_median_c2(3, :));
            median_amp(j, iii, 1, 4) = median(envelope_median_c1(4, :));
            median_amp(j, iii, 2, 4) = median(envelope_median_c2(4, :));
        end
        
        
        
        if group_level == 1
            if iii == length(cohort)
                idx_delete = [1 2 3 6 7];
                median_freq(:, idx_delete, :, :) = [];
                median_amp(:, idx_delete, :, :) = [];
                
                SR_freq(j, 1) = signrank(abs(median_freq(j, :, 1, 1) - median_freq(j, :, 2, 1))); % X
                SR_freq(j, 2) = signrank(abs(median_freq(j, :, 1, 2) - median_freq(j, :, 2, 2))); % Y
                SR_freq(j, 3) = signrank(abs(median_freq(j, :, 1, 3) - median_freq(j, :, 2, 3))); % Z
                SR_freq(j, 4) = signrank(abs(median_freq(j, :, 1, 4) - median_freq(j, :, 2, 4))); % Vector
                
                SR_amp(j, 1) = signrank(abs(median_amp(j, :, 1, 1) - median_amp(j, :, 2, 1))); % X
                SR_amp(j, 2) = signrank(abs(median_amp(j, :, 1, 2) - median_amp(j, :, 2, 2))); % Y
                SR_amp(j, 3) = signrank(abs(median_amp(j, :, 1, 3) - median_amp(j, :, 2, 3))); % Z
                SR_amp(j, 4) = signrank(abs(median_amp(j, :, 1, 4) - median_amp(j, :, 2, 4))); % Vector
            end
        end

        if group_level_full == 1
            if iii == length(cohort)
                idx_delete = [1 2 3 6];
                median_freq(:, idx_delete, :, :) = [];
                median_amp(:, idx_delete, :, :) = [];
                
                SR_freq(j, 1) = signrank(abs(median_freq(j, :, 1, 1) - median_freq(j, :, 2, 1))); % X
                SR_freq(j, 2) = signrank(abs(median_freq(j, :, 1, 2) - median_freq(j, :, 2, 2))); % Y
                SR_freq(j, 3) = signrank(abs(median_freq(j, :, 1, 3) - median_freq(j, :, 2, 3))); % Z
                SR_freq(j, 4) = signrank(abs(median_freq(j, :, 1, 4) - median_freq(j, :, 2, 4))); % Vector
                
                SR_amp(j, 1) = signrank(abs(median_amp(j, :, 1, 1) - median_amp(j, :, 2, 1))); % X
                SR_amp(j, 2) = signrank(abs(median_amp(j, :, 1, 2) - median_amp(j, :, 2, 2))); % Y
                SR_amp(j, 3) = signrank(abs(median_amp(j, :, 1, 3) - median_amp(j, :, 2, 3))); % Z
                SR_amp(j, 4) = signrank(abs(median_amp(j, :, 1, 4) - median_amp(j, :, 2, 4))); % Vector
            end
        end 
        
    end
     
end

%% Group level

if shuffle_data == 0
if group_level == 1 || group_level_full == 1
    envelope_median = [envelope_median_ns', envelope_median_rs'];
    median_freq(iii, 1, 1) = median(mean_freq_amp(1, find(c == 1)));
    median_freq(iii, 2, 1) = median(mean_freq_amp(1, find(c == 2)));
    median_freq(iii, 1, 2) = median(mean_freq_amp(2, find(c == 1)));
    median_freq(iii, 2, 2) = median(mean_freq_amp(2, find(c == 2)));
    median_freq(iii, 1, 3) = median(mean_freq_amp(3, find(c == 1)));
    median_freq(iii, 2, 3) = median(mean_freq_amp(3, find(c == 2)));
    median_freq(iii, 1, 4) = median(mean_freq_amp(4, find(c == 1)));
    median_freq(iii, 2, 4) = median(mean_freq_amp(4, find(c == 2)));
    
    median_amp(iii, 1, 1) = median(envelope_median(1, find(c == 1)));
    median_amp(iii, 2, 1) = median(envelope_median(1, find(c == 2)));
    median_amp(iii, 1, 2) = median(envelope_median(2, find(c == 1)));
    median_amp(iii, 2, 2) = median(envelope_median(2, find(c == 2)));
    median_amp(iii, 1, 3) = median(envelope_median(3, find(c == 1)));
    median_amp(iii, 2, 3) = median(envelope_median(3, find(c == 2)));
    median_amp(iii, 1, 4) = median(envelope_median(4, find(c == 1)));
    median_amp(iii, 2, 4) = median(envelope_median(4, find(c == 2)));
end


% cohort = [1 4 5 6 7 8 10 11 12 14];
% final_cohort = [6 7 11 12 14];
% idx_final = [4 5 8 9 10];
if group_level == 1
if iii == length(cohort)
    idx_delete = [1 2 3 6 7];
    median_freq(idx_delete, :, :) = [];
    median_amp(idx_delete, :, :) = [];
    
    [SR_freq(1, 1), SR_freq(1, 2)] = signrank(abs(median_freq(:, 1, 1) - median_freq(:, 2, 1))); % X
    [SR_freq(2, 1), SR_freq(2, 2)] = signrank(abs(median_freq(:, 1, 2) - median_freq(:, 2, 2))); % Y
    [SR_freq(3, 1), SR_freq(3, 2)] = signrank(abs(median_freq(:, 1, 3) - median_freq(:, 2, 3))); % Z
    [SR_freq(4, 1), SR_freq(4, 2)] = signrank(abs(median_freq(:, 1, 4) - median_freq(:, 2, 4))); % Vector
    
    [SR_amp(1, 1), SR_amp(1, 2)] = signrank(abs(median_amp(:, 1, 1) - median_amp(:, 2, 1))); % X
    [SR_amp(2, 1), SR_amp(2, 2)] = signrank(abs(median_amp(:, 1, 2) - median_amp(:, 2, 2))); % Y
    [SR_amp(3, 1), SR_amp(3, 2)] = signrank(abs(median_amp(:, 1, 3) - median_amp(:, 2, 3))); % Z
    [SR_amp(4, 1), SR_amp(4, 2)] = signrank(abs(median_amp(:, 1, 4) - median_amp(:, 2, 4))); % Vector
end
end

if group_level_full == 1
if iii == length(cohort)
    idx_delete = [1 2 3 6];
    median_freq(idx_delete, :, :) = [];
    median_amp(idx_delete, :, :) = [];
    
    [SR_freq(1, 1), SR_freq(1, 2)] = signrank(abs(median_freq(:, 1, 1) - median_freq(:, 2, 1))); % X
    [SR_freq(2, 1), SR_freq(2, 2)] = signrank(abs(median_freq(:, 1, 2) - median_freq(:, 2, 2))); % Y
    [SR_freq(3, 1), SR_freq(3, 2)] = signrank(abs(median_freq(:, 1, 3) - median_freq(:, 2, 3))); % Z
    [SR_freq(4, 1), SR_freq(4, 2)] = signrank(abs(median_freq(:, 1, 4) - median_freq(:, 2, 4))); % Vector
    
    [SR_amp(1, 1), SR_amp(1, 2)] = signrank(abs(median_amp(:, 1, 1) - median_amp(:, 2, 1))); % X
    [SR_amp(2, 1), SR_amp(2, 2)] = signrank(abs(median_amp(:, 1, 2) - median_amp(:, 2, 2))); % Y
    [SR_amp(3, 1), SR_amp(3, 2)] = signrank(abs(median_amp(:, 1, 3) - median_amp(:, 2, 3))); % Z
    [SR_amp(4, 1), SR_amp(4, 2)] = signrank(abs(median_amp(:, 1, 4) - median_amp(:, 2, 4))); % Vector
end
end 
end
%%
clearvars -except fig fig2 cohort iii C_RS randstim_change method pc_idx evaluate f1 optk per_neg critval mean_coeff test_frequencies frequency_amplitude test_baseline power_amplitude variance_explained median_explained_ns mean_explained_ns std_explained_ns median_explained_rs mean_explained_rs std_explained_rs counts_ns counts_half1 counts_half2 median_freq SR_freq group_level SR_amp median_amp group_level_full KS_p_freq KS_h_freq T_p_freq T_h_freq RSum_p_freq RSum_h_freq KS_p_power KS_h_power T_p_power T_h_power RSum_p_power RSum_h_power tuning_curve amp_n_bins1 amp_n_bins2 bins main shuffle_data

end

%%
if tuning_curve == 1
    clearvars -except amp_n_bins1 amp_n_bins2 bins cohort
    cl = 'r';
    for i = 1:size(amp_n_bins1, 1)
            
        f1 = figure(100)
        subplot(2, 5, i)
        y = amp_n_bins1(i, :);
        x = bins(1 : end-2);
        bar(x, y, 'FaceColor', cl, 'EdgeColor', cl)
        hold on
        if length(find(isnan(amp_n_bins1(i, :)) == 1)) < size(amp_n_bins1(i, :), 2) - 2
            [rsg, rsg_g, rsg_o] = gauss_fit2(x, y)
            ylim([0 1.5])
            xticks([1 : 2 : 14])
    
            ylabel({'Absolute amplitude (\muV^2)'})
            xlabel('Frequency (Hz)')

            cv1(i, :) = rsg.c ./ rsg.b;
        end
        
        box('off')
        legend('off')
        title(sprintf('Patient %d', cohort(i)))
        clear x y rsg rsg_o rsg_g
    end

   f1.Units = 'centimeters';
    f1.OuterPosition = [10, 10, 32, 15];
    set(f1, 'color', 'w')
    
    for i = 1:size(amp_n_bins2, 1)
        f2 = figure(101)
        subplot(2, 5, i)
        y = amp_n_bins2(i, :);
        x = bins(1 : end-2);
        bar(x, y, 'FaceColor', cl, 'EdgeColor', cl)
        hold on
        if length(find(isnan(amp_n_bins2(i, :)) == 1)) < size(amp_n_bins2(i, :), 2)
            [rsg, rsg_g, rsg_o] = gauss_fit2(x, y)
            ylim([0 1.5])
            xticks([1 : 2 : 14])
        
            ylabel({'Absolute amplitude (\muV^2)'})
            xlabel('Frequency (Hz)')
           
            cv2(i, :) = rsg.c ./ rsg.b;
        end
        
        box('off')
        legend('off')
        title(sprintf('Patient %d', cohort(i)))
        clear x y rsg rsg_o rsg_g       
    end
       
    f2.Units = 'centimeters';
    f2.OuterPosition = [10, 10, 32, 15];
    set(f2, 'color', 'w')
end


%%



function v_new = split_vector(v, n)
L = floor(length(v) / n); % Length of each block
idx_st = NaN(100, 1);
idx_st(1, 1) = 1;
v_new = NaN(100, L);
for i = 2:n
   idx_st(i) = idx_st(i-1) + L;
   v_new(i-1,:) = v(1, idx_st(i-1) : (idx_st(i)-1));
end
v_new(n,:) = v(1, idx_st(n) : idx_st(n) + L - 1);
end

function v_new = split_vector_length(v, L)
n = floor(length(v) / L); % Number of segments
idx_st = NaN(n, 1);
idx_st(1, 1) = 1;

idx_end = NaN(n, 1);
idx_end(1, 1) = L;

v_new = NaN(n, L);
for i = 2:n
   idx_st(i) = idx_st(i-1) + L;
   idx_end(i) = idx_end(i-1) + L;
end

for i = 1:n
   v_new(i, :) = v(1, idx_st(i, 1) : idx_end(i, 1));
end
end


