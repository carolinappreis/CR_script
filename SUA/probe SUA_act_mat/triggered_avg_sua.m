clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
load('SNR_sua_skrate.mat') ; 
 

% subj= SNR.beta_rats(ismember(SNR.beta_rats,BZ.beta_rats));

subj= SNR.beta_rats;
color_b=[0 0 0.5]; %snr
% color_b=[0.5 0 0.5]; %bz



Ecog=SNR.ctx(subj,:);
for i =1:size(subj,1)
    for ii=1:size(SNR.beta_idx{subj(i),1},2)
        hp=SNR.beta_idx{subj(i),1}(1,ii);
        data_all{i,1}(ii,:)=SNR.sua{subj(i),1}(hp,:);
        clear hp
    end
end


srn=1000;

for j=1:size(data_all,1)
    
    data=data_all{j,1};
    ecog=Ecog(j,:);
    clearvars -except j ii SNR Ecog data ecog data_all srn rec_pa1 rec_npa1 rec_pa2 rec_npa2 region_pl region_spl region_npl region_nspl subj color_b
    
    for ii=1:size(data,1)
        
        [b,a]=butter(2,[15/(0.5*srn) 35/(0.5*srn)],'bandpass');
        Ecogfiltered=filtfilt(b,a,ecog);
        env=abs(hilbert(Ecogfiltered));
        onset1=bursts(env);
        onset1=horzcat(onset1{:});
        onset=bursts_aligned(env,Ecogfiltered);
        onset=horzcat(onset{:});
        
        data_g=smoothdata(data(ii,:),'gaussian',25);
        
        for jj=1:size(onset,2)
            if onset(jj)>200 && onset(jj)+200<length(data_g)
                output_pa(ii,jj,:)=data_g(onset(jj)-200:onset(jj)+200);
                output_npa(ii,jj,:)= data_g(onset1(jj)-200:onset1(jj)+200);
                output_count(ii,jj,:)= data(onset1(jj)-200:onset1(jj)+200);
            end
        end
        rec_pa1{j,1}=reshape(sum(output_pa,2),(size(sum(output_pa,2),1)),(size(sum(output_pa,2),3)));
        rec_npa1{j,1}=reshape(sum(output_npa,2),(size(sum(output_npa,2),1)),(size(sum(output_npa,2),3)));
    end
    rec_pa2(j,:)=mean(rec_pa1{j,1},1);
    rec_npa2(j,:)=mean(rec_npa1{j,1},1);
    clear data
end
rec_pa=cell2mat(rec_pa1);
rec_npa=cell2mat(rec_npa1);

region_pl=zscore(mean(rec_pa,1));
region_spl=zscore(std(rec_pa)./sqrt(size(rec_pa,1)));
region_npl=zscore(mean(rec_npa,1));
region_snpl=zscore(std(rec_npa)./sqrt(size(rec_npa,1)));

time2=[1:401];
color_s=[0.5 0.5 0.5];



subplot(1,2,1)
y2=region_pl; y1=y2+region_spl; y3=y2-region_spl;
y5=region_npl; y4=y5+region_snpl; y6=y5-region_snpl;
p1=plot(time2, y2, 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
xlim ([0 400])
ylim ([-5 5])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
legend([p1],{'phase-aligned'},'FontSize',12)
box ('off')
%         xlabel ('Time (msec)')
%         ylabel('Firing-rate(z-score)')


subplot(1,2,2)
p2=plot(time2, y5, 'LineWidth',1.5,'Color',color_s)
patch([time2 fliplr(time2)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
xlim ([0 400])
ylim ([-5 5])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
legend([p2],{'non-phase aligned'},'FontSize',12)
box ('off')



% clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
% load('rec_pa2_snr.mat')
% s=rec_pa2;
% clear rec_pa2
% load('rec_pa2_bz.mat')
% b=rec_pa2;
% clear rec_pa2
% for i=1:22
%     subplot(3,10,i)
%     plot(b(i,:))
%     hold on
%     plot(s(i,:))
% box('off')
% end


for nn=3
for i =1:size(s{nn,1},1)
subplot(size(s{nn,1},1),1,i)
plot(s{1,1}(i,:))
xlim ([0 10000])
end
figure(2)
for i =1:size(b{nn,1},1)
subplot(size(b{nn,1},1),1,i)
plot(b{1,1}(i,:))
xlim ([0 10000])
end
end