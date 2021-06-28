    

clear; close all
cohort = [ 1 3 4 6];
spiral=1;

freq_bl = NaN(length(cohort),3,5e4);
ns_mat=NaN(length(cohort),3);


 for iii = 1:length(cohort)
    clearvars -except iii cohort ns_mat spiral freq_bl
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
   
   
   if spiral==0
       load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_out.mat','bs_end','bs_begin','amp_bbl','change_bl','m_ax')
   else
       load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_spiral.mat','bs_end','bs_begin','amp_bbl','change_bl','m_ax')
   end
   
   
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
    

    mmm=NaN(3,5e4);
    parfor j = 1:length(bs_end)
        begin3=bs_begin(iii,j);
        end3=bs_end(iii,j);
        for ax = 1:3
            [mm]=simp(ax,iii,end3,begin3,freqi,phase,ns_mat);
            mmm(ax,j)=mm;
        end
    end
    freq_bl(iii,:,:)=mmm;
       
 end