clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')

for rg = 1:2;
    if rg==1;
        load('data_SUA_BZ.mat');
    elseif rg==2;
        load('data_SUA_SNR.mat');
    end
    freq=[10;22;32;43;77]; srn=1000;
    dataregion=cell2mat(data_region);
    for t=1:size(freq,1)
        for  j=1:size(data_region,1)
            for r=1:size(data_region{j,1},1)
                clear data; data=data_region{j,1}(r,:);
                data_ones=find(data==1);
                [b,a]=butter(2,[(freq(t)-5)/(0.5*srn) (freq(t)+5)/(0.5*srn)],'bandpass');
                if t==length(freq)
                    [b,a]=butter(2,[49/(0.5*srn) 100/(0.5*srn)],'bandpass');
                end
                Ecogfiltered=filtfilt(b,a,Ecog_region(j,:));
                ang=angle(hilbert(Ecogfiltered(data_ones)))';
                vec_an(r,:)=circ_r(ang);
                clear ang;
            end
            vec(t,j,:)=mean(vec_an);
            clear vec_an;
        end
    end
    reg(rg,:)=mean(vec,2);
end

bar(abs(euler1))
title ('Phase consistency of unit firing in cortical signal')
legend('VA','VL','VM')
xticklabels({'5-15Hz','16-26Hz','27-37Hz','38-48Hz','49-100Hz'})


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