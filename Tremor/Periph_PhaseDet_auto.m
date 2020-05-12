clear all, close all
cohort = [ 2 3 4 5 8 10 11 13 16 17];

method = 'ward-pca2';

load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/clusters_BA'),'C_RS')
% % load(strcat('C:\Users\beatrizsda\OneDrive - Nexus365\Peripheral Stim\Matlab Outputs\',method,'\C_RS'))
% load(strcat('/Users/beatrizarruda/OneDrive - Nexus365/Peripheral Stim/Matlab Outputs/',method,'/C_RS'))
TT1_C1 = cell(length(cohort),3);
TT1_C2 = cell(length(cohort),3);
TT1_A1 = cell(length(cohort),3);
TT1_A2 = cell(length(cohort),3);

if strcmp(method, 'ward-power') || strcmp(method, 'ward-pca2') || strcmp(method, 'ward-pca3')
    method = 'ward';
end

Nc1 = NaN(length(cohort),1);
Nc2 = NaN(length(cohort),1);
seg_zenv1=cell(length(cohort),12);
seg_env1=cell(length(cohort),12);
seg_zfilt1=cell(length(cohort),12);
seg_filt1=cell(length(cohort),12);
seg_zenv2=cell(length(cohort),12);
seg_env2=cell(length(cohort),12);
seg_zfilt2=cell(length(cohort),12);
seg_filt2=cell(length(cohort),12);


for iii = 1:length(cohort)
    
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(cohort(iii)),'_RS.mat'))
    
    c_rs = C_RS{iii,1};
    
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
    filtt=[tremorxf;tremoryf; tremorzf];
    envelope = [abs(hilbert(tremorxf)); abs(hilbert(tremoryf)); abs(hilbert(tremorzf))];
    phase = [angle(hilbert(tremorxf)); angle(hilbert(tremoryf)); angle(hilbert(tremorzf))];
    tremor_zs = [zscore(tremorx);zscore(tremory); zscore(tremorz)];
    filt_zs = [filtfilt(b, a,tremor_zs(1,:)); filtfilt(b, a,tremor_zs(2,:)); filtfilt(b, a,tremor_zs(3,:));];
    env_zs = [abs(hilbert(filt_zs(1,:))); abs(hilbert(filt_zs(2,:))); abs(hilbert(filt_zs(3,:)))];
    
    new = find(data(2, :) > 4);
    difp = find((diff(new)) > 100000); % are you trying to threshold at 9.6 seconds?
    ep_1 = [new(difp) new(end)];
    sp_1 = [new(1) new(difp+1)];
    
    %%% input start all trial
    start_t = 1;
    sp = sp_1(1, start_t:end);
    ep = ep_1(1, start_t:end);
    
    %     plot(time,data(4,:))
    %     hold on
    %     plot(time(sp),data(4,sp),'r.')
    %     plot(time(ep),data(4,ep),'k.')
    
    
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
    
    
    % figure()
    % plot(time, data(4, :))
    % hold on
    % plot(time(index), data(4, index), 'r.')
    % plot(time(start1), data(4, start1), 'ko')
    % plot(time(ending1), data(4, ending1), 'bo')
    
    clear start ending
    start{1, 1} = floor((start1 ./ samplerateold) * samplerate);%+addon;
    ending{1, 1} = floor((ending1 ./ samplerateold) * samplerate);%+addon+addon_end;%floor(5*samplerate);
    clear xx
    xx{1, 1} = xx1;
    
    %%
    
    hh = numel(start);
    
    
    c_1 = find(c_rs == 1);
    c_2 = find(c_rs == 2);
    
    
    
    %     f2 = figure()
    
    
    
    
    %% %%%%%% Cluster 1 %%%%%%
    if ~isempty(c_1)
        
        c = c_1;
        Nc1(iii,1) = length(c);
        start_c = start{hh, 1}(c);
        ending_c = ending{hh, 1}(c);
        xx_c = xx{1, 1}(c);
        
        
        %%
        
        %-----
        for j = 1:length(start_c)
            if (~isnan(start_c(j)))
                x = [sum(envelope(1, start_c(j):ending_c(j))); sum(envelope(2, start_c(j):ending_c(j))); sum(envelope(3, start_c(j):ending_c(j)))];
                y = [tremorxf(start_c(j):ending_c(j)); tremoryf(start_c(j):ending_c(j)); tremorzf(start_c(j):ending_c(j))];
                [pc, score, latent, tsquare] = pca(y');
                yyy(j, 1:3) = pc(1:3, 1);
                
                
                PSI_c(j, 1:2) = [abs(sum(exp(sqrt(-1) .* (phase(1, start_c(j) : ending_c(j))-phase(2, start_c(j):ending_c(j))))) ./ length(phase(1, start_c(j):ending_c(j)))); abs(sum(exp(sqrt(-1) .* (phase(1, start_c(j):ending_c(j))-phase(3, start_c(j):ending_c(j))))) ./ length(phase(1, start_c(j):ending_c(j))))];
                ma_c(j) = find(x == max(x));
                ma2_c(j) = (find(abs(yyy(j, 1:3)) == max(abs(yyy(j, 1:3)))));
                
                clear x y
            end
        end
        
        tremor_or2_c = NaN(length(start_c), 1);
        tremor_pc_c = NaN(length(start_c), 1);
        bs_env = NaN(length(start_c), 1);
        
        
        for axx = 1:3
            for i = 1:length(start_c)
                if (~isnan(start_c(i)))
                    tremor_or2_c(axx, i, 1) = (mean(envelope(axx, ending_c(i)-1000:ending_c(i))) - mean(envelope(axx, start_c(i) - 1000:start_c(i)))) / mean(envelope(axx, start_c(i) - 1000:start_c(i)));
                    tremor_pc_c(1, i) = (mean(envelope(ma_c(i), ending_c(i)-1000:ending_c(i))) - mean(envelope(ma_c(i), start_c(i) - 1000:start_c(i)))) / mean(envelope(ma_c(i), start_c(i) - 1000:start_c(i)));
                    z_env(axx,i,1:5000)=env_zs(axx,start{hh,1}(i):start{hh,1}(i)+5000-1);
                    z_filt(axx,i,1:5000)=filt_zs(axx,start{hh,1}(i):start{hh,1}(i)+5000-1);
                    n_env(axx,i,1:5000)=envelope(axx,start{hh,1}(i):start{hh,1}(i)+5000-1);
                    n_filt(axx,i,1:5000)=filtt(axx,start{hh,1}(i):start{hh,1}(i)+5000-1);
                    
                    bs_env(axx,i,1)=mean(envelope(axx,start{hh,1}(i)-1000:start{hh,1}(i)));
                    
                    
                    xx_c(i) = xx_c(i);
                else
                    tremor_or2_c(axx, i, 1) = NaN;
                    tremor_pc_c(i, 1) = NaN;
                    z_env(axx,i,1:5000)=NaN;
                    z_filt(axx,i,1:5000)=NaN;
                    n_env(axx,i,1:5000)=NaN;
                    n_filt(axx,i,1:5000)=NaN;
                    
                    bs_env(axx,i,1)=NaN;
                    
                    xx_c(i) = NaN;
                end
            end
            
            tt_c = NaN(20, 12);
            yy_c = xx_c(:);
            tt_pc_c = NaN(20, 12);
            tt_abstim = NaN(20, 12);
            
            
            for i = 1:12
                tt_c(1:sum(yy_c == i), i) = tremor_or2_c(axx, find(yy_c == i));
                tt_abstim(1:sum(yy_c == i), i) = bs_env(axx, find(yy_c == i));
                tt_pc_c(1:sum(yy_c == i), i) = tremor_pc_c(1, find(yy_c == i));
                if length(c)>=10
                    seg_zenv1{iii,i}(axx,:,:)=squeeze(z_env(axx,find(yy_c==i),:));
                    seg_zfilt1{iii,i}(axx,:,:)=squeeze(z_filt(axx,find(yy_c==i),:));
                    seg_env1{iii,i}(axx,:,:)=squeeze(n_env(axx,find(yy_c==i),:));
                    seg_filt1{iii,i}(axx,:,:)=squeeze(n_filt(axx,find(yy_c==i),:));
                end
                
                tt_c(tt_c == 0) = NaN;
                tt_abstim(tt_abstim==0)=NaN;
                tt_pc_c(tt_pc_c == 0) = NaN;
            end
            
            tt1_c{hh, axx} = tt_c;
            tt1_abs{hh,axx}=tt_abstim;
        end
        
        
        
        pca_ax_c{hh, 1} = ma_c;
        PSI_ax_c{hh, 1} = mean(PSI_c, 1);
        
        
        clear yy ma;
        
        %%
        %         % Plots
        %         for i = 1
        %             a_c = hist(pca_ax_c{i, 1}, 1:3);
        %             axmax_c = find(a_c == max(a_c));
        %             subplot(2, 4, 1)
        %             temp = tt1_c{i, axmax_c};
        %             bar(0:30:330, 100 .* nanmedian(temp))
        %             hold on
        %             plot(0:30:330, 100 .* temp,'.')
        %             set(gca, 'XTickLabelRotation', 45)
        %             box('off')
        %             ylabel('Cluster 1')
        %             subplot(2, 4, 2)
        %             bar(0:30:330, 100 .* nanmedian(temp))
        %             set(gca,'XTickLabelRotation', 45)
        %             box('off')
        %             subplot(2, 4, 3)
        %             bar(1:3, a_c)
        %             names = {'CED2'; 'CED5';'CED6'};
        %             set(gca, 'xtick', [1:3], 'xticklabel', names)
        %             box('off')
        %             subplot(2, 4, 4)
        %             bar([1 2], [PSI_ax_c{i}(1, 1); PSI_ax_c{i}(1, 2)]);
        %             names = {'PSI 2 5'; 'PSI 2 6'};
        %             box('off')
        %
        %             set(gca, 'xtick', [1:2], 'xticklabel', names)
        %             f2.Units = 'centimeters';
        %             f2.OuterPosition = [4, 5, 45, 20];
        %             set(f2, 'color', 'w');
        %         end
        
        if length(c) >= 10
            TT1_C1{iii, 1} = tt1_c{1, 1};
            TT1_C1{iii, 2} = tt1_c{1, 2};
            TT1_C1{iii, 3} = tt1_c{1, 3};
            TT1_A1{iii, 1} = tt1_abs{1, 1};
            TT1_A1{iii, 2} = tt1_abs{1, 2};
            TT1_A1{iii, 3} = tt1_abs{1, 3};
        end
        
    end
    
    clearvars -except cohort iii start ending envelope phase tremorxf tremoryf tremorzf xx hh C_RS c_2 f2 TT1_A1 TT1_A2 TT1_C1 TT1_C2 Nc1 Nc2 seg_env1 seg_filt1 seg_env2 seg_filt2 seg_zenv1 seg_zfilt1 seg_zenv2 seg_zfilt2 env_zs filt_zs filtt
    
    %% %%%%%% Cluster 2 %%%%%%
    if ~isempty(c_2)
        
        c = c_2;
        Nc2(iii, 1) = length(c);
        
        start_c = start{hh, 1}(c);
        ending_c = ending{hh, 1}(c);
        xx_c = xx{1, 1}(c);
        
        
        %%
        
        %-----
        for j = 1:length(start_c)
            if (~isnan(start_c(j)))
                x = [sum(envelope(1, start_c(j):ending_c(j))); sum(envelope(2, start_c(j):ending_c(j))); sum(envelope(3, start_c(j):ending_c(j)))];
                y = [tremorxf(start_c(j):ending_c(j)); tremoryf(start_c(j):ending_c(j)); tremorzf(start_c(j):ending_c(j))];
                [pc, score, latent, tsquare] = pca(y');
                yyy(j, 1:3) = pc(1:3, 1);
                
                
                PSI_c(j, 1:2) = [abs(sum(exp(sqrt(-1) .* (phase(1, start_c(j) : ending_c(j))-phase(2, start_c(j):ending_c(j))))) ./ length(phase(1, start_c(j):ending_c(j)))); abs(sum(exp(sqrt(-1) .* (phase(1, start_c(j):ending_c(j))-phase(3, start_c(j):ending_c(j))))) ./ length(phase(1, start_c(j):ending_c(j))))];
                ma_c(j) = find(x == max(x));
                ma2_c(j) = (find(abs(yyy(j, 1:3)) == max(abs(yyy(j, 1:3)))));
                
                clear x y
            end
        end
        
        tremor_or2_c = NaN(length(start_c), 1);
        tremor_pc_c = NaN(length(start_c), 1);
        bs_env = NaN(length(start_c), 1);
        
        
        for axx = 1:3
            for i = 1:length(start_c)
                if (~isnan(start_c(i)))
                    tremor_or2_c(axx, i, 1) = (mean(envelope(axx, ending_c(i)-1000:ending_c(i))) - mean(envelope(axx, start_c(i) - 1000:start_c(i)))) / mean(envelope(axx, start_c(i) - 1000:start_c(i)));
                    tremor_pc_c(1, i) = (mean(envelope(ma_c(i), ending_c(i)-1000:ending_c(i))) - mean(envelope(ma_c(i), start_c(i) - 1000:start_c(i)))) / mean(envelope(ma_c(i), start_c(i) - 1000:start_c(i)));
                    z_env(axx,i,1:5000)=env_zs(axx,start{hh,1}(i):start{hh,1}(i)+5000-1);
                    z_filt(axx,i,1:5000)=filt_zs(axx,start{hh,1}(i):start{hh,1}(i)+5000-1);
                    n_env(axx,i,1:5000)=envelope(axx,start{hh,1}(i):start{hh,1}(i)+5000-1);
                    n_filt(axx,i,1:5000)=filtt(axx,start{hh,1}(i):start{hh,1}(i)+5000-1);
                    
                    bs_env(axx,i,1)=mean(envelope(axx,start{hh,1}(i)-1000:start{hh,1}(i)));
                    
                    
                    xx_c(i) = xx_c(i);
                else
                    tremor_or2_c(axx, i, 1) = NaN;
                    tremor_pc_c(i, 1) = NaN;
                    z_env(axx,i,1:5000)=NaN;
                    z_filt(axx,i,1:5000)=NaN;
                    n_env(axx,i,1:5000)=NaN;
                    n_filt(axx,i,1:5000)=NaN;
                    
                    bs_env(axx,i,1)=NaN;
                    
                    xx_c(i) = NaN;
                end
            end
            
            tt_c = NaN(20, 12);
            yy_c = xx_c(:);
            tt_pc_c = NaN(20, 12);
            tt_abstim = NaN(20, 12);
            
            
            for i = 1:12
                tt_c(1:sum(yy_c == i), i) = tremor_or2_c(axx, find(yy_c == i));
                tt_abstim(1:sum(yy_c == i), i) = bs_env(axx, find(yy_c == i));
                tt_pc_c(1:sum(yy_c == i), i) = tremor_pc_c(1, find(yy_c == i));
                if length(c)>=10
                    seg_zenv2{iii,i}(axx,:,:)=squeeze(z_env(axx,find(yy_c==i),:));
                    seg_zfilt2{iii,i}(axx,:,:)=squeeze(z_filt(axx,find(yy_c==i),:));
                    seg_env2{iii,i}(axx,:,:)=squeeze(n_env(axx,find(yy_c==i),:));
                    seg_filt2{iii,i}(axx,:,:)=squeeze(n_filt(axx,find(yy_c==i),:));
                end
                
                tt_c(tt_c == 0) = NaN;
                tt_abstim(tt_abstim==0)=NaN;
                tt_pc_c(tt_pc_c == 0) = NaN;
            end
            
            tt2_c{hh, axx} = tt_c;
            tt2_abs{hh,axx}=tt_abstim;
        end
        
        
        
        pca_ax_c{hh, 1} = ma_c;
        PSI_ax_c{hh, 1} = mean(PSI_c, 1);
        
        
        clear yy ma;
        
        %%
        % Plots
        %         for i = 1
        %             a_c = hist(pca_ax_c{i, 1}, 1:3);
        %             axmax_c = find(a_c == max(a_c));
        %             subplot(2, 4, 5)
        %             temp = tt1_c{i, axmax_c};
        %             bar(0:30:330, 100 .* nanmedian(temp))
        %             hold on
        %             plot(0:30:330, 100 .* temp,'.')
        %             set(gca, 'XTickLabelRotation', 45)
        %             box('off')
        %             ylabel('Cluster 2')
        %             subplot(2, 4, 6)
        %             bar(0:30:330, 100 .* nanmedian(temp))
        %             set(gca,'XTickLabelRotation', 45)
        %             box('off')
        %             subplot(2, 4, 7)
        %             bar(1:3, a_c)
        %             names = {'CED2'; 'CED5';'CED6'};
        %             set(gca, 'xtick', [1:3], 'xticklabel', names)
        %             box('off')
        %             subplot(2, 4, 8)
        %             bar([1 2], [PSI_ax_c{i}(1, 1); PSI_ax_c{i}(1, 2)]);
        %             names = {'PSI 2 5'; 'PSI 2 6'};
        %             box('off')
        %
        %             set(gca, 'xtick', [1:2], 'xticklabel', names)
        %             f2.Units = 'centimeters';
        %             f2.OuterPosition = [4, 5, 45, 20];
        %             set(f2, 'color', 'w');
        %         end
        %
        if length(c) >= 10
            TT1_C2{iii, 1} = tt2_c{1, 1};
            TT1_C2{iii, 2} = tt2_c{1, 2};
            TT1_C2{iii, 3} = tt2_c{1, 3};
            TT1_A2{iii, 1} = tt2_abs{1, 1};
            TT1_A2{iii, 2} = tt2_abs{1, 2};
            TT1_A2{iii, 3} = tt2_abs{1, 3};
        end
        
    end
    %     filename=['arc_clusters',num2str(iii),'.fig'];
    %     saveas(gcf,filename)
    %     filename=['arc_clusters',num2str(iii),'.svg'];
    %     saveas(gcf,filename)
    
    clearvars -except cohort iii C_RS TT1_C1 TT1_C2  TT1_A1 TT1_A2 Nc1 Nc2 seg_env1 seg_filt1 seg_env2 seg_filt2 seg_zenv1 seg_zfilt1 seg_zenv2 seg_zfilt2
    close all
end
