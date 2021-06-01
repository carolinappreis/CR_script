
%%% choose between 1) patient 3 (DT) 3) patient 6 (ET) stim at 120 deg 4)
%%% pateint 6 (ET) at 240
% 
clear; close
  load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P03_PLS_P')
  load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_out.mat','rs_mat')
rs_max=rs_mat(2,1); clear rs_mat; cr=[];

% clear; close
% load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P06_PLS_P1')
%   load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_out.mat','rs_mat')
% rs_max=rs_mat(4,1);cr=[];


% clear; close
%  load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P06_PLS_P2')
%   load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_out.mat','rs_mat')%%% to get the main axis and points of tp_s(tapping start) and tp_e(tapping ending)
% rs_max=rs_mat(4,1); cr=1;
%% PLS06 P2

%%% -------

[d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;

co=3;
iii=1;
clust=struct; out=struct; start=cell(1,3); ending=cell(1,3); yy=cell(1,3); h_up=cell(1,3); s=struct;

[peak_ax, start, ending, yy, h_up]= dbs_startend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up);

[s]=dbs_zfiltenv(d,peak_ax,co,iii,s,samplerate);

dum=s.filt{1,3};
dum2=s.env_acc{1,3};

r=0;
for t=1:size(tp_s,1)
    for j = 1:length(tp_s{t,1})
        r=r+1;
        if (~isnan(tp_s{1,1}(j)))
            for ii=1:3
                x(ii,:) = dum(ii,(tp_s{t,1}(j)-3000):tp_s{t,1}(j));
                y(ii,:) = mean(dum2(ii,(tp_s{t,1}(j)-3000):tp_s{t,1}(j)));
            end
            [pc, score, latent, tsquare, explained] = pca(x');
            plot_pc(:,r)=pc(1:3,1);
            acc(:,r)=y;
            [k,kk]=sort(pc(1:3,1),'descend');
            pc_before(t,j) = kk(1);
            dif_before(t,j) = k(1)-k(2);
        end
        clearvars -except t j pc_before dum tp_s tp_e start dif_before ending cr plot_pc r acc dum2
    end
end



plot(plot_pc','.')
box('off')
xlim([0 size(plot_pc,2)+1])


if cr==1
ending{1,3}(1,2)=258848;
end

dr=ending{1,3};
for j=1:size(dr,2)
    if (~isnan(dr(j)))
        for ii=1:3
            x(ii,:) = dum(ii,(dr(j)-3000):dr(j));
        end
        [pc, score, latent, tsquare, explained] = pca(x');
        [k,kk]=sort(pc(1:3,1),'descend');
        pc_last(1,j) = kk(1);
        dif_last(1,j)=k(1)-k(2);
    end
    clearvars -except t j pc_before dum tp_s tp_e pc_last cr dif_last dif_before start ending dr
end


for t=1:size(tp_s,1)
    for j = 1:length(tp_s{t,1})
        if (~isnan(tp_s{1,1}(j)))
            for ii=1:3
                x(ii,:) = dum(ii,(tp_e{t,1}(j)):tp_e{t,1}(j)+3000);
            end
            [pc, score, latent, tsquare, explained] = pca(x');
            [k,kk]=sort(pc(1:3,1),'descend');
            pc_after(t,j) = kk(1);
            dif_after(t,j) = k(1)-k(2);
        end
        clearvars -except t j pc_after pc_before dum tp_s tp_e pc_last dif_before dif_last dif_after start ending
    end
end



dr=start{1,3};
for j=1:size(dr,2)
    if (~isnan(dr(j)))
        for ii=1:3
            x(ii,:) = dum(ii,(dr(j)-3000):dr(j));
        end
        [pc, score, latent, tsquare, explained] = pca(x');
        [k,kk]=sort(pc(1:3,1),'descend');
        pc_first(1,j) = kk(1);
        dif_first(1,j)=k(1)-k(2);
    end
        clearvars -except pc_first dif_first cr t j pc_after pc_before dum tp_s tp_e pc_last dif_before dif_last dif_after start ending dr
end

map=[pc_first' pc_before(:,1) pc_after(:,1) pc_before(:,2) pc_after(:,2) pc_last'];

map1=[ dif_first' dif_before(:,1) dif_after(:,1) dif_before(:,2) dif_after(:,2) dif_last'];


figure(1)
for i=1:3
    subplot(1,3,i)
    bar(map(i,:))
    xticklabels({'bef stim on';'bef tap';'after tap';'bef tap';'after tap';'bef stim off'})
    xtickangle(45)
    ylabel('tremor axis')
    box('off')
    ylim([0 3])
    yticks([1:3])
end

figure(2)
for i=1:3
    subplot(1,3,i)
    bar(map1(i,:))
    xticklabels({'bef stim on';'bef tap';'after tap';'bef tap';'after tap';'bef stim off'})
    xtickangle(45)
    ylabel('difference of pca contrib')
    box('off')
end
