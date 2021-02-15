clear; close all
cohort = [ 1 3 4 6];
spiral=1;

% nostim= NaN(length(cohort),3,1e6);
% nostim_f=NaN(length(cohort),3,1e6);
tt1_all=cell(length(cohort),3);
x_all=cell(length(cohort),1);
bs_begin=NaN(length(cohort),5e4);
bs_end=NaN(length(cohort),5e4);
amp_bbl = NaN(length(cohort),3,5e4);
change_bl = NaN(length(cohort),3,5e4);
pc1_exp =cell(length(cohort),1);
m_ax=NaN(1,length(cohort));
ns_mat=NaN(length(cohort),3);
rs_mat=NaN(length(cohort),3);
st_sp=cell(length(cohort),1); 
et_sp=cell(length(cohort),1);


 for iii = 1:length(cohort)
    
    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_RS.mat'))
    
    
    in2 = 1; % analysing the "main tremor axis"
    
    if in2 == 1
        in = 3;
    elseif in2 == 2 % other axis 1
        in = 6;
    elseif in2 == 3 % other axis 2
        in = 7;
    end
    
    %
    data = SmrData.WvData;
    addon=92; addon_end=35;
    
    
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
    
    start = floor((indexes4 ./ samplerateold)*samplerate)+addon;
    ending = floor((indexes3 ./ samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
    
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
    tremor = (data(6, :));% %score(:,1)';%
    ts = timeseries(tremor, 0:(1/samplerateold):((size(data, 2)-1)/samplerateold));
    ts1 = resample(ts, 0:0.001:((size(data, 2)-1)/samplerateold), 'linear');
    tremory(1:size(ts1.data, 3)) = ts1.data;
    tremor = (data(7, :));% %score(:,1)';%
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
    start_t=1;
    
    sp=sp_1(1,start_t:end);
    ep=ep_1(1,start_t:end);
    
    %%% posture/spiral trials
    dt1_s=[sp_1(start_t:2:end)];dt1_e=[ep_1(start_t:2:end)];
    dt2_s=sp_1(start_t+1:2:end); dt2_e=ep_1(start_t+1:2:end);
    
    %         time=1:length(data(1,:));
    %         plot(time,data(4,:))
    %         hold on
    %         plot(time(sp),data(4,sp),'r.')
    %         plot(time(ep),data(4,ep),'k.')
    
    
    for ik=1:length(sp) %%find double start and end points in a stimulation run
        
        s=(find(([indexes4>=sp(ik)]+[indexes4<=ep(ik)])==2));
        e=(find(([indexes3>=sp(ik)]+[indexes3<=ep(ik)])==2));
        tks=(find(diff(xx(s))==0))+1;
        tke=(find(diff(xx(e))==0));
        
        indexes4(s(tks))=NaN;
        indexes3(e(tke))=NaN;
        xx(e(tke))=NaN;
        
    end
    
    
    % % %     %%%% find runs with trigering issues (too few, too many pulses)
    % % %     th1=(Fpeak*5)./2;
    % % %     th2=(Fpeak*5)+round((Fpeak*5)./5);
    % % %     for it=1:length(indexes4)
    % % %         if numel(index(find(index==indexes4(it)):find(index==indexes3(it))))>=th1 && numel(index(find(index==indexes4(it)):find(index==indexes3(it))))<=th2
    % % %             indexes4(it)=indexes4(it);
    % % %             indexes3(it)=indexes3(it);
    % % %             xx(it)=xx(it);
    % % %         else
    % % %             indexes4(it)=NaN;
    % % %             indexes3(it)=NaN;
    % % %             xx(it)=NaN;
    % % %         end
    % % %     end
    % % %     %%%%%%%%%%%%%%%
    indexes4=indexes4(~isnan(indexes4));
    indexes3=indexes3(~isnan(indexes3));
    xx=xx(~isnan(xx));
    
    
    start1=[];
    ending1=[];
    xx1=[];
    for il=1:length(dt1_s)
        start1=[start1 indexes4(find(([indexes4>=dt1_s(il)]+[indexes4<=dt1_e(il)])==2))]; % intersect([1 2 3],[3 4 5])
        ending1=[ending1 indexes3(find(([indexes3>=dt1_s(il)]+[indexes3<=dt1_e(il)])==2))];
        xx1=[xx1 xx(find(([indexes4>=dt1_s(il)]+[indexes4<=dt1_e(il)])==2))];
    end
    
    start2=[];
    ending2=[];
    xx2=[];
    for il=1:length(dt2_s)
        dums=indexes4(find(([indexes4>=dt2_s(il)]+[indexes4<=dt2_e(il)])==2));
        start2=[start2 dums]; clear dums
        dume=indexes3(find(([indexes3>=dt2_s(il)]+[indexes3<=dt2_e(il)])==2));
        ending2=[ending2 dume]; clear dume
        dumx=xx(find(([indexes4>=dt2_s(il)]+[indexes4<=dt2_e(il)])==2));
        xx2=[xx2 dumx];clear dumx
    end
    
    
%         figure()
%         plot(time,data(4,:))
%         hold on
%         plot(time(index),data(4,index),'r.')
%         plot(time(start1),data(4,start1),'ko')
%         plot(time(ending1),data(4,ending1),'bo')
    
    
    clear start ending xx
    pstart{1,1}=floor((start1./samplerateold)*samplerate)+addon;
    pending{1,1}=floor((ending1./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
    pstart{2,1}=floor((start2./samplerateold)*samplerate)+addon;
    pending{2,1}=floor((ending2./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
    
    %%%% saving for cluster analysis of spirals
    st_sp{iii,1}=pstart{2,1};
    et_sp{iii,1}=pending{2,1};
    
    
    xx{1,1}=xx1;
    xx{2,1}=xx2;
    
    %%% choosing start{1,1}/ending{1,1}/xx{1,1} to get posture only _ check!
    if spiral ==0
    start{1,1}= pstart{1,1};
    ending{1,1} = pending{1,1};
    yy{1,1}= xx{1,1};
    else
    start{1,1}= pstart{2,1};
    ending{1,1} = pending{2,1};
    yy{1,1}= xx{2,1};
    end
    
     
    for j = 1:length(start{1, 1})
        if (~isnan(start{1,1}(j)))
        yd=[sum(envelope(1,start{1,1}(j):ending{1,1}(j)));sum(envelope(2,start{1,1}(j):ending{1,1}(j)));sum(envelope(3,start{1,1}(j):ending{1,1}(j)))];
        ma_c(j) = find(yd == max(yd)); clear yd
        x = [tremorxf(start{1,1}(j):ending{1,1}(j)); tremoryf(start{1,1}(j):ending{1,1}(j)); tremorzf(start{1,1}(j):ending{1,1}(j))];
        [pc, score, latent, tsquare, explained] = pca(x');
        pc_trials(j, 1:3) = pc(1:3, 1);
        explained_rs(j, 1:3) = explained;
        end
    end
    
    a=hist(ma_c,1:3);
    m_ax(1,iii)=find(a==max(a));
   if m_ax(1,iii)==1
       rs_mat(iii,:)=[1 2 3];
   elseif m_ax(1,iii)==2
       rs_mat(iii,:)=[2 3 1];
   else
       rs_mat(iii,:)=[3 2 1];
   end
    
    tremor_or2=NaN(3,length(start{1,1}));
    
    yy=cell2mat(yy);
    for axx=1:3
        for i=1:length(start{1,1})
            if (~isnan(start{1,1}(i)))
                tremor_or2(axx,i)=(mean(envelope(rs_mat(iii,axx),ending{1,1}(i)-1000:ending{1,1}(i)))-mean(envelope(rs_mat(iii,axx),start{1,1}(i)-1000:start{1,1}(i))))/mean(envelope(rs_mat(iii,axx),start{1,1}(i)-1000:start{1,1}(i)));
                yy(i)= yy(i);
                
            else
                tremor_or2(axx,i,1)=NaN;
                yy(i)= NaN;
            end
        end
        
        %         %% criteria for outliers
        %
        %         idx_outl=find(tremor_or2>(nanmean(tremor_or2)+(2*(nanstd(tremor_or2))))|tremor_or2<(nanmean(tremor_or2)-(2*(nanstd(tremor_or2)))));
        %         tremor_or2(idx_outl,1)=NaN;
        %         tremor_or3(idx_outl,1)=NaN;
        %         xx(1,idx_outl)=NaN;
        
        tt=NaN(20,12);
        tt_pc=NaN(20,12);
        
        
        for i=1:12
            tt(1:sum(yy==i),i)=tremor_or2(axx,find(yy==i));
            
            
            tt(tt==0)=NaN;
            
        end
        tt1_all{iii,axx}=tt;
        
    end
    
 %% Baseline
    clearvars -except iii cohort main method tt1_all ns_mat amp_bbl bs_begin bs_end pc_trials change_bl x_all pc1_exp explained_rs m_ax ns_mat rs_mat spiral
    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_NS.mat'))
    
    
    rng('default') % set random seed for consistency
    gen=rng;
    
    in2 = 1; % analysing the "main tremor axis"
    
    if in2 == 1
        in = 3;
    elseif in2 == 2 % other axis 1
        in = 6;
    elseif in2 == 3 % other axis 2
        in = 7;
    end
    
    data = SmrData.WvData;
    samplerateold = SmrData.SR;
    tremor = (data(in, :));
    
    ts = timeseries(data, 0:(1 / samplerateold):((size(data, 2)-1) / samplerateold));
    ts1 = resample(ts, 0:0.001:((size(data, 2)-1)/samplerateold), 'linear');
    ds_data(1:size(ts1.data, 1), 1:size(ts1.data, 3)) = ts1.data;
    samplerate = 1000;
    tre_3 = ds_data([3 6 7],:); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK
    
    
    for aa = 1:3
        [Pxx,F] = pwelch(tre_3(aa,:), samplerate, [], samplerate, samplerate);
        frange = F(3:10);
        Pxxrange = Pxx(3:10);
        Freqpeak(aa,:) = frange(find(Pxxrange == max(Pxxrange)));
        Ppeak(aa,:) = max(Pxxrange);
        ps_curves(aa,:) = Pxx;
    end
    
    peak_ax = [(Freqpeak(find(Ppeak == max(Ppeak)))) (find(Ppeak == max(Ppeak)))];
    
    
   f1=peak_ax(2);
   
   if (m_ax(1,iii)==1 && f1==1 || m_ax(1,iii)==3 && f1==3)
       ns_mat(iii,:)=[1 2 3];
   elseif (m_ax(1,iii)==1 && f1==3 || m_ax(1,iii)==3 && f1==1)
       ns_mat(iii,:)=[3 2 1];
   elseif (m_ax(1,iii)==1 && f1==2 || m_ax(1,iii)==2 && f1==1)
       ns_mat(iii,:)=[2 3 1];
   end
    
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
    
    if spiral==0
        NS_BE_P
    else
        NS_BE_S
    end
    
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
    pc1_exp{iii,1}(1:3,:)=[(explained_ns(:,1:3))' (explained_rs(:,1:3))'];
   clearvars -except  cohort iii nostim tt1_all main ns_mat amp_bbl bs_begin bs_end  change_bl x_all pc1_exp m_ax rs_mat st_sp et_sp spiral
end

clearvars -except tt1_all  amp_bbl bs_begin bs_end change_bl x_all pc1_exp m_ax ns_mat rs_mat spiral
