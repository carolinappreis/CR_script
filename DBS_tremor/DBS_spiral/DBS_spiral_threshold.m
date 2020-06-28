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
    ts=timeseries(data_raw,0:(1/samplerateold):((size(data_raw,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data_raw,2)-1)/samplerateold),'linear');
    ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    tt=0:1/samplerate:(size(ds_data,2)-1)/samplerate;
    data_t=ds_data([3 6 7],:);
    
    NS_BE_S
    out.b_trials=hu{iii,:};
    out.e_trials=hd{iii,:};
    
    
%         figure(1)
%         time=1:length(data_t(1,:));
%         plot(time,data_t(1,:))
%         hold on
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
    
    [dd, cc] = butter(2,[1/(0.5*samplerate) ],'low'); %15
    
    [b,a]=butter(2,[0.8/(0.5*samplerate) ],'low'); %15
    
    
    for i=1:3
        out.filt{iii,co}(i,:)=filtfilt(afilt,bfilt,data_t(i,:));
        out.l_filt{iii,co}(i,:)=filtfilt(b,a,data_t(i,:));
        out.w_filt{iii,co}(i,:)=filtfilt(dd,cc,data_t(i,:));
        out.env{iii,co}(i,:)=abs(hilbert(out.filt{iii,co}(i,:)));
    end
    
    x_signal=out.filt{iii,co};

    L = 1000;
    
    for o=1:3
        [v_new,out,e_new]= split_vector_length(x_signal(o,out.h_up{1,1}),out.env{iii,co}(o,out.h_up{1,1}),L,out,iii,co);
        ns_sseg{1,o}=v_new; clear v_new
        e_sseg{1,o}=e_new;
    end
    
    new_sig=[];
    new_e=[];
    for j = 1:size(ns_sseg{1,o},1)
        x = [ns_sseg{1,1}(j,:); ns_sseg{1,2}(j,:); ns_sseg{1,3}(j,:)];
        [pc, score, latent, tsquare, explained] = pca(x');
        comb_pc(j,1:3)=pc(1:3,1);
        maxx=find(pc(:, 1)==max(pc(:,1)));
        new_sig=[new_sig ns_sseg{1,maxx}(j,:)];
        new_e=[new_e e_sseg{1,maxx}(j,:)];
    end
    
    figure()
    plot(comb_pc(:,1), 'r.')
    hold on
    plot(comb_pc(:,2), 'b.')
    plot(comb_pc(:,3), 'k.')
    xlim([0,size(comb_pc, 1)])
    
    
    th=prctile(new_e,50);
    id_sp1=find(new_e>th);
    id_sp2=find(new_e<th);
    
    tt= 1:length(new_e);
    
    figure()
    plot(tt,new_sig,'Color', [0.5 0.5 0.5])
    hold on
    plot(tt,new_e,'Color', [0.5 0.5 0.5])
    plot(tt(id_sp1),new_e(id_sp1),'r.')
    plot(tt(id_sp2),new_e(id_sp2),'b.')


    figure()
    chosen_d=out.l_filt{iii,co};
    time=1:length(chosen_d(1,:));
    plot3(time(1,out.h_up{1,1}),chosen_d(2,out.h_up{1,1}),chosen_d(3,out.h_up{1,1}),'Color',[0.5 0.5 0.5])
    hold on
    plot3(time(1,out.h_up{1,1}(id_sp1)),chosen_d(2,out.h_up{1,1}(id_sp1)),chosen_d(3,out.h_up{1,1}(id_sp1)),'r.')
    plot3(time(1,out.h_up{1,1}(id_sp2)),chosen_d(2,out.h_up{1,1}(id_sp2)),chosen_d(3,out.h_up{1,1}(id_sp2)),'b.')
    
     
end