

% --- to find animals

clear all
close all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load ('data_all')
load 'BZ_ctx_probe'
idx=[];
for i=1:size(WaveData_DCall,1)
if size(WaveData_DCall{i,1},1)~=1
idx=[idx i];
end
end

% cd('C:\Users\creis\Documents\GitHub\CRcode\codes_thal\A4_Thal')\
cd('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/A4_Thal')
load ('bz_allrats.mat','animals')
% load ('snr_lowbeta.mat','animals')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
load 'animal_lesion_nolesion.mat'
A(lesion(idx(animals)),4:9)
cd('/Users/Carolina/Documents/MATLAB/KN')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\MATLAB\KN')
% 
% %-----



clear all
close all
cd ('C:\Users\creis\Documents\GitHub\CRcode\codes_thal\A4_Thal')
% cd('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/A4_Thal')
% %       load 'bz_allrats.mat'; color_b=[0.5 0 0.5];id_rat=([1 7 9 11])'; rats=size(id_rat,1);
%      load 'snr_allrats.mat'; color_b=[0 0 0.5];id_rat=([1 3 11 12])'; rats=size(id_rat,1);
color_s=[0 0 0];
color_ctxb=[0.5 0.5 0.5];

for i =1:length(id_rat)
clust_b{id_rat(i),:};
clust_b_m(i,:)=mean(ans);
clust_b_sd(i,:)=std(ans)./sqrt(size(ans,1));
clust_s{id_rat(i),:};
clust_s_m(i,:)=mean(ans);
clust_s_sd(i,:)=std(ans)./sqrt(size(ans,1));
power_coh{id_rat(i),:};
power_coh_m(i,:)=mean(ans);
power_coh_sd(i,:)=std(ans)./sqrt(size(ans,1));
power_ctx1(i,:)=power_ctx{id_rat(i),1}(1,:);
end

ctx_b=ctx_b1(id_rat,:);
ctx_b_m=mean(ctx_b,1);
ctx_b_sd=std(ctx_b)./sqrt(size(ctx_b,1));

power_ctx_m=mean(power_ctx1);
power_ctx_sd=std(power_ctx1)./sqrt(size(power_ctx1,1));

cd('C:\Users\creis\Documents\GitHub\CRcode\codes_thal\A4_Thal\code')
for r=1:size(id_rat,1);
    if size(clust_b{id_rat(r),:},1)>1
        st=NaN(1,401);
        clear A; A=clust_b{id_rat(r),:}; %b1{f,1};
        clear B; B=clust_s{id_rat(r),:}; %s1{f,1}(1:size(A,1),:);
        
        hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
        beg=find(st(1,:)<0.01 & st2(1,:)~=0);
        if ~isempty(beg)
            sig_rise(r,:)=[beg(1) beg(end)];
        end
    end
end
% ---------- rats with coherent beta ctx-sub and sign dif with surrogates
% -individual

 time= [1:2*200+1];
for i =1:(rats) 
    subplot(rats+1,1,1)
    y7=ctx_b_m; y6=y7+ctx_b_sd; y8=y7-ctx_b_sd;
    plot(time, y7, 'DisplayName','BZ aligned to short bursts')
    set(plot(time,y7),'LineWidth',1.5,'Color',color_ctxb);
    patch([time fliplr(time)], [y6 fliplr(y7)],[color_ctxb],'FaceAlpha',[0.2],'EdgeColor','none')
    patch([time fliplr(time)], [y7 fliplr(y8)],[color_ctxb],'FaceAlpha',[0.2],'EdgeColor','none')
    xlim ([0 400])
    xticks([0:100:400])
    xticklabels ({'-200','-100','0','100','200'})
    hold on
    
    subplot(rats+1,1,i+1)
    y2=clust_b_m(i,:); y1=y2+clust_b_sd(i,:); y3=y2-clust_b_sd(i,:);
    y5=clust_s_m(i,:); y4=y5+clust_s_sd(i,:); y6=y5-clust_s_sd(i,:);
    plot(time, y2, 'DisplayName','BZ aligned to short bursts')
    set(plot(time,y2),'LineWidth',1.5,'Color',color_b);
    patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
    patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
    hold on
    plot(time, y5, 'DisplayName','BZ aligned to short bursts')
    set(plot(time,y5),'LineWidth',1.5,'Color',color_s);
    patch([time fliplr(time)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
    patch([time fliplr(time)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
    xlim ([0 400])
    ylim ([0 0.6])
    xticks([0:100:400])
    xticklabels ({'-200','-100','0','100','200'})
    
    patch([sig_rise(i,1) sig_rise(i,2) sig_rise(i,2) sig_rise(i,1)],[-0.1 -0.1 0.6 0.6],color_b,'FaceAlpha',[0.1],'EdgeColor','none')
    
end

figure()
    time1=10:50;
for i =1:rats
    subplot(2,1,1)
    y1=power_ctx1(i,time1);
    plot(time1,log(y1), 'LineWidth',1.5,'Color',color_ctxb)
    hold on
    xlabel ('Frequency (Hz)')
    ylabel('Log(power) M1')
    
    subplot(2,1,2)
    y2=power_coh_m(i,time1);
    plot(time1, log(y2), 'LineWidth',1.5,'Color',color_b)
    hold on
    xlabel ('Frequency (Hz)')
    ylabel('Log(power) BZ')
    
end


% ----- mean


% ylabel('Phase Synchrony Index')
% title('BZ-CTX coupling 15-30Hz') 
% title('CTX {\beta} burst')
