clear; close
cohort = [ 1 3 4 6];
cond={'NS';'RS'};
pc_t=cell(size(cohort,2),size(cond,1));
start=cell(10,3); ending=cell(10,3); yy=cell(10,3); out=struct; clust=struct;

for iii =  1
    clearvars -except  cohort cond iii pc_t start ending yy out clust
    co=1
    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
    
  
    data_raw=SmrData.WvData;
    samplerateold=SmrData.SR;
    ts=timeseries(d.data_raw,0:(1/d.samplerateold):((size(d.data_raw,2)-1)/d.samplerateold));
    ts1=resample(ts,0:0.001:((size(d.data_raw,2)-1)/d.samplerateold),'linear');
    ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
    d.samplerate=1000;
    tt=0:1/d.samplerate:(size(ds_data,2)-1)/d.samplerate;
    data_t=ds_data([3 6 7],:);
    
    NS_BE_S
    out.b_trials=hu{iii,:};
    out.e_trials=hd{iii,:};
    
    
    %     figure(1)
    %     time=1:length(data_t(1,:));
    %     plot(time,data_t(1,:))
    %     hold on
    handup = [];
    for i = 1:length(out.b_trials)
        handup = [handup out.b_trials(i):out.e_trials(i)]; %#ok<*AGROW>
        gi{i,1}= [out.b_trials(i):out.e_trials(i)];
        %         plot(time(gi{i,1}),data_t(1,gi{i,1}))
        %         hold on
    end
    clear i
    handup = sort(handup,'ascend');
    out.h_up=gi;
    
    for aa = 1:3
        [Pxx,F] = pwelch(data_t(aa,handup), samplerate, [], round(samplerate), samplerate);
        frange = F(3:10);
        Pxxrange = Pxx(3:10);
        Freqpeak(aa,:) = frange(find(Pxxrange == max(Pxxrange)));
        Ppeak(aa,:) = max(Pxxrange);
        ps_curves(aa,:) = Pxx;
    end
    
    peak_ax = [(Freqpeak(find(Ppeak == max(Ppeak)))) (find(Ppeak == max(Ppeak)))];
    Fpeak = round(peak_ax(1));
    
    if (Fpeak-2) >= 1
        [afilt, bfilt] = butter(2, [(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
    else
        [afilt, bfilt] = butter(2, [(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)], 'bandpass'); %15
    end
    
    [dd, cc] = butter(2,[5/(0.5*samplerate) ],'low'); %15
    
    [b,a]=butter(2,[0.8/(0.5*samplerate) ],'low'); %15
    
    
    for i=1:3
        out.filt{iii,co}(i,:)=filtfilt(afilt,bfilt,data_t(i,:));
        out.l_filt{iii,co}(i,:)=filtfilt(b,a,data_t(i,:));
        out.w_filt{iii,co}(i,:)=filtfilt(dd,cc,data_t(i,:));
        out.env{iii,co}=abs(hilbert(out.filt{iii,co}(i,:)));
    end
    
    out.data_t{iii,co}=data_t;
    x_signal=out.filt{iii,co};
    
    data=out.w_filt{iii,1};
    data1=out.l_filt{iii,1};
    time=1:length(data(1,:));
    figure(20)
    plot(time,out.env{iii,1}(1,:),'Color',[0.5 0.5 0.5])
    hold on
    figure(10)
    %         plot3(time(1,out.h_up{iii,1}{lh,1}),data(2,out.h_up{iii,1}{lh,1}),data(3,out.h_up{iii,1}{lh,1}),'Color',[0.5 0.5 0.5])
    plot3(time,data(2,:),data(3,:),'Color',[0.5 0.5 0.5])
    hold on
    xlabel('time')
    ylabel('y axis')
    zlabel('z axis')
    tt={'NS1';'NS2';'NS3'};
    

    L = 1000;
    
    for o=1:3
        [v_new,out]= split_vector_length(x_signal(o,out.h_up), L,out,iii,co,o);
        ns_sseg{1,o}=v_new; clear v_new
    end
    
    
    for j = 1:size(ns_sseg{lh,o},1)
        x = [ns_sseg{1,1}(j,:); ns_sseg{1,2}(j,:); ns_sseg{1,3}(j,:)];
        [pc, score, latent, tsquare, explained] = pca(x');
        maxx=find(pc(:, 1)==max(pc(:,1)));
% new_seg{
    end
    figure(4)
    subplot(size(out.h_up{1,1},1),1,lh)
    plot(pc_trials_ns(:,1), 'r.')
    hold on
    plot(pc_trials_ns(:,2), 'b.')
    plot(pc_trials_ns(:,3), 'k.')
    xlim([0,size(pc_trials_ns, 1)])
    title(sprintf(tt{lh,1}))
    %
    
    pc_t{1,lh}=pc_trials_ns;
    
end