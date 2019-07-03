clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')

for rg = 1:2;
    if rg==1;
        load('data_SUA_BZ.mat');
    elseif rg==2;
        load('data_SUA_SNR.mat');
    end
%     freq=[10;22;32;43;77]; 
    freq={[1 7];[8 12];[15 35];[36 80];[85 100]};
    srn=1000;
    dataregion=cell2mat(data_region);
    for t=1:size(freq,1)
        for  j=1:size(data_region,1)
            for r=1:size(data_region{j,1},1)
                clear data; data=data_region{j,1}(r,:);
                sk_ctc=[];
                for i =1:srn:(length(data)-srn);
                    sk_ctc=[sk_ctc numel(find(data(1,i:i+srn)==1))];
                end
                sk_ctc1(r,:)=mean(sk_ctc);
                data_ones=find(data==1);
                [b,a]=butter(2,[(freq{t,1}(1))/(0.5*srn) (freq{t,1}(2))/(0.5*srn)],'bandpass');
                Ecogfiltered=filtfilt(b,a,Ecog_region(j,:));
                ang=angle(hilbert(Ecogfiltered(data_ones)))';
                vec_an(r,:)=circ_r(ang);
                clear ang;
            end
            vec(t,j,:)=mean(vec_an);
            spikerate{rg,1}(j,:)=mean(sk_ctc1);
            clear vec_an  sk_ctc1
        end
    end
    reg(rg,:)=mean(vec,2);
end

fig=figure()
histogram(spikerate{1,1},5,'FaceColor',[0.5 0 0],'EdgeColor',[0.4 0 0],'FaceAlpha',0.8)
ylabel('Unit count');
ylim([0 5])
xlabel('Frequency (Hz)');
title('Firing rate BZ')
box('off')
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 6, 6];
fig.Color='w';
set(gca,'FontSize',12)

fig=figure()
histogram(spikerate{2,1},5,'FaceColor',[0 0 0.5],'EdgeColor',[0 0 0.4],'FaceAlpha',0.8)
ylabel('Unit count');
ylim([0 5])
xlabel('Frequency (Hz)');
title('Firing rate SNR')
box('off')
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 6, 6];
fig.Color='w';
set(gca,'FontSize',12)

fig=figure()
bar(reg')
title('Phase locking to ECoG')
legend('BZ','SNR','box','off')
ylabel('Vector length');
xlabel('Frequecy (Hz)');
xticklabels({'5-15','16-26','27-37','38-48','49-100'})
box('off')
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 12,10];
fig.Color='w';
set(gca,'FontSize',12)


% cd('/Users/Carolina/Desktop/TC_data')
% savefig('sua_phaseconsist')
% saveas(gca,'sua_phaseconsist.png')

% for t=1:4
%     subplot (1,4,t)
%     polarhistogram(ang1{t,3},12)
% end
% title ('mean phase VM')

% figure()
% for i=1:size(data_all,2)
%     subplot(size(data_all,2),1,i)
%     plot(time,data_all{3,i})
%     ylim ([-2 2])
% %     xlim ([0 10])
% end


% plot (time,Ecogfiltered)
% xlim ([0 0.1])
% figure()
% plot(time,data_all{1,1})
% xlim ([0 0.1])