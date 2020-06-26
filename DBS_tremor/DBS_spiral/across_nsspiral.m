clear; close
cohort = [ 1 3 4 6];
cond={'NS';'RS'};
pc_t=cell(size(cohort,2),size(cond,1));
start=cell(10,3); ending=cell(10,3); yy=cell(10,3); out=struct; clust=struct;

for iii =  1
    clearvars -except  cohort cond iii pc_t start ending yy out clust
    co=1
    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
    
    data=SmrData.WvData;
    samplerate=SmrData.SR;
    data_t=data([3 6 7],:);
    samplerateold=1000;
    
    NS_BE_S
    out.b_trials{iii,co}=floor((hu{iii,:}./samplerateold)*samplerate);
    out.e_trials{iii,co}=floor((hd{iii,:}./samplerateold)*samplerate);
    
    
    %     figure(1)
    %     time=1:length(data_t(1,:));
    %     plot(time,data_t(1,:))
    %     hold on
    handup = [];
    for i = 1:length(out.b_trials{iii,co})
        handup = [handup out.b_trials{iii,co}(i):out.e_trials{iii,co}(i)]; %#ok<*AGROW>
    end
    clear i
    handup = sort(handup,'ascend');
    out.h_up{iii,co}=handup;
    
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
    
    [dd, cc] = butter(2, [0.5/(0.5*samplerate) (6)/(0.5*samplerate)], 'bandpass'); %15
    
    [b,a]=butter(2,[0.8/(0.5*samplerate) ],'low'); %15
    
    
    for i=1:3
        out.filt{iii,co}(i,:)=filtfilt(afilt,bfilt,data_t(i,:));
        out.l_filt{iii,co}(i,:)=filtfilt(b,a,data_t(i,:));
        out.w_filt{iii,co}(i,:)=filtfilt(dd,cc,data_t(i,:));
    end
    
    out.data_t{iii,co}=data_t;
    x_signal=out.filt{iii,co};
    
    for lh=1:size(out.h_up{1,1},1)
        L = floor((1000./samplerateold)*samplerate);
        
        for o=1:3
            [v_new,out]= split_vector_length(x_signal(o,out.h_up{iii,co}), L,out,iii,co,o);
            ns_sseg{lh,o}=v_new; clear v_new
        end
        
        
        for j = 1:size(ns_sseg{lh,o},1)
            x = [ns_sseg{lh,1}(j,:); ns_sseg{lh,2}(j,:); ns_sseg{lh,3}(j,:)];
            [pc, score, latent, tsquare, explained] = pca(x');
            pc_trials_ns(j, 1:3) = pc(1:3, 1);
            explained_ns(j, 1:3) = explained;
        end
        figure(4)
        plot(pc_trials_ns(:,1), 'r.')
        hold on
        plot(pc_trials_ns(:,2), 'b.')
        plot(pc_trials_ns(:,3), 'k.')
        xlim([0, size(pc_trials_ns, 1)])
        title('NS')
        
        pc_t{1,lh}=pc_trials_ns;
        clear  pc_trials_ns v_new
        
    end
    
    
    x_all{1,1}=[pc_t{1,1}];
    runs{1,1}=1:length(x_all{1,1});
    k=3;
    Z = linkage(x_all{iii,1},'ward');
    c = cluster(Z, 'Maxclust', k);
    
    
    figure()
    scatter3([pc_t{iii,1}(:, 1)], [pc_t{iii,1}(:, 2)], [pc_t{iii,1}(:, 3)], 10, c)
    xlabel('x axis');ylabel('y axis');zlabel('z axis')
    
    
    figure()
    [l h]=silhouette(x_all{iii,1},c);
    title('ALL')
    
    for i=1:1
        clust.C{iii,i}=c(runs{1,i});
        clust.idx{iii,i}=[{find(clust.C{iii,i}==1)}; {find(clust.C{iii,i}==2)}; {find(clust.C{iii,i}==3)}];
        clust.count{iii,1}(1:3,i)=[numel(find(clust.C{iii,i}==1)); numel(find(clust.C{iii,i}==2)); numel(find(clust.C{iii,i}==3))];
        clust.percent{iii,1}(1:3,i)=[numel(find(clust.C{iii,i}==1))./length(runs{1,i})*100; numel(find(clust.C{iii,i}==2))./length(runs{1,i})*100; numel(find(clust.C{iii,i}==3))./length(runs{1,i})*100];
    end
    
    p=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];
    tt={'NS1';'NS2';'NS3'};
    
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
   
    
    for cl=1:3
        id=[clust.idx{iii,co}{cl,1}];
        if ~isempty(id)
            for ii=1:length(id)
                dummy=out.h_up{iii,1}{lh,1}(out.start{iii,1}(id(ii),1):out.ending{iii,1}(id(ii),1));
                figure(10)
                plot3(time(dummy),data1(2,(dummy)),data1(3,(dummy)),'.','Color',p(cl,:))
                %                 plot3(data(1,(dummy)),data(2,(dummy)),data(3,(dummy)),'.','Color',p(cl,:))
                figure(20)
                plot(time(dummy),out.env{iii,1}(1,(dummy)),'.','Color',p(cl,:))
                
                clear dummy
            end
            clear id
        end
    end
    view(-30,15)
end