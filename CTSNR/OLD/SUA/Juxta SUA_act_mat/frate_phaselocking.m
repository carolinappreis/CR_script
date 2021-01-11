% % % 
% % % 
% % % clear all
% % % % cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
% % % cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
% % % 
% % % for rg = 1
% % %     if rg==1
% % %         load('SUA_BZ.mat');
% % %     elseif rg==2
% % %         load('SUA_CZ.mat');
% % %     end
% % %     
% % %     data_region=data_all;
% % %     Ecog_region=WaveData_DC;
% % %     %     freq=[10;22;32;43;77];
% % %     freq={[1 7];[8 12];[15 35];[36 80];[81 100]};
% % %     
% % %     srn=1000;
% % %     for t=1:size(freq,1)
% % %         for r=1:size(data_region,1)
% % %             data=data_region(r,:);
% % %             sk_ctc=[];
% % %             for i =1:srn:(length(data)-srn);
% % %                 sk_ctc=[sk_ctc numel(find(data(1,i:i+srn)==1))];
% % %             end
% % %             data_ones=find(data==1);
% % %             [b,a]=butter(2,[(freq{t,1}(1))/(0.5*srn) (freq{t,1}(2))/(0.5*srn)],'bandpass');
% % %             Ecogfiltered=filtfilt(b,a,Ecog_region(r,:));
% % %             ang=angle(hilbert(Ecogfiltered(data_ones)))';
% % %             vec_an(r,:)=circ_r(ang);
% % %             clear ang;
% % %             spikerate{rg,1}(r,:)=mean(sk_ctc);
% % %         end
% % %         vec(t,:)=mean(vec_an);
% % %         clear vec_an data
% % %     end
% % %     reg(rg,:)=mean(vec,2);
% % % end
% % % 
% % % if rg~=1
% % %     
% % %     fig=figure()
% % %     subplot(1,3,1)
% % %     title('BZ')
% % %     histogram(spikerate{1,1},5)
% % %     ylabel('Unit count');
% % %     xlabel('Frequency (Hz)');
% % %     title('Firing rate BZ')
% % %     box('off')
% % %     
% % %     subplot(1,3,2)
% % %     title('CZ')
% % %     histogram(spikerate{2,1},5)
% % %     ylabel('Unit count');
% % %     xlabel('Frequency (Hz)');
% % %     title('Firing rate CZ')
% % %     box('off')
% % %     
% % %     subplot(1,3,3)
% % %     bar(reg')
% % %     title ('Phase locking to ECoG')
% % %     legend('BZ','CZ')
% % %     xticklabels({'1-7','16-26','27-37','38-48','49-100'})
% % %     ylabel('Vector length');
% % %     xlabel('Frequecy (Hz)');
% % %     fig.Color='w';
% % %     box('off')
% % %     fig.Units = 'centimeters';
% % %     fig.OuterPosition= [10, 10, 25, 8];
% % %     fig.Color='w';
% % %     
% % % else
% % %     
% % %     fig=figure()
% % %     subplot(1,2,1)
% % %     title('BZ')
% % %     histogram(spikerate{1,1},5)
% % %     ylabel('Unit count');
% % %     xlabel('Frequency (Hz)');
% % %     title('Firing rate BZ')
% % %     box('off')
% % %     
% % %     
% % %     subplot(1,2,2)
% % %     bar(reg')
% % %     title ('Phase locking of BZ units to ECoG')
% % %     xticklabels({'1-7','16-26','27-37','38-48','49-100'})
% % %     ylabel('Vector length');
% % %     xlabel('Frequecy (Hz)');
% % %     
% % %     fig.Color='w';
% % %     box('off')
% % %     fig.Units = 'centimeters';
% % %     fig.OuterPosition= [10, 10, 26, 10];
% % %     fig.Color='w';
% % %     
% % %     
% % %     
% % % end



clear all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load('data_SUA_SNR.mat')


rg=1;

%     freq=[10;22;32;43;77];
freq={[1 7];[8 12];[15 35];[36 80];[81 100]};
srn=1000;
for t=1:size(freq,1)
    p=0;
    for r=1:size(data_region,1)
        for ii=1:size(data_region{r,1},1)
            p=p+1;
            data=data_region{r,1}(ii,:);
            sk_ctc=[];
            for i =1:srn:(length(data)-srn);
                sk_ctc=[sk_ctc numel(find(data(1,i:i+srn)==1))];
            end
            data_ones=find(data==1);
            [b,a]=butter(2,[(freq{t,1}(1))/(0.5*srn) (freq{t,1}(2))/(0.5*srn)],'bandpass');
            Ecogfiltered=filtfilt(b,a,Ecog_region(r,:));
            ang=angle(hilbert(Ecogfiltered(data_ones)))';
            vec_an(p,:)=circ_r(ang);
            clear ang;
            spikerate{rg,1}(p,:)=mean(sk_ctc);
        end
    end
    vec(t,:)=mean(vec_an);
    clear vec_an data
end
reg(rg,:)=mean(vec,2);


fig=figure()
subplot(1,2,1)
title('SNr')
histogram(spikerate{1,1},5)
ylabel('Unit count');
xlabel('Frequency (Hz)');
title('Firing rate SNr')
box('off')


subplot(1,2,2)
bar(reg')
title ('Phase locking of SNr units to ECoG')
xticklabels({'1-7','16-26','27-37','38-48','49-100'})
ylabel('Vector length');
xlabel('Frequecy (Hz)');

fig.Color='w';
box('off')
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 26, 10];
fig.Color='w';





