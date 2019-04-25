clear all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
load('output_count_bz.mat')
load('rec_pa_bz.mat')
%
% data=region_pl{2,1};
time=1:401;
% [maxvalM,maxidxM] = findpeaks(data);
% plot(time,data)
% hold on
% plot(time(maxidxM(1,11:14)),data(maxidxM(1,11:14)),'r.')
% idx_int=maxidxM(11:14);
% spikes=output_count{9,1};
%
% join=[];
% for i=1:size(spikes,1)
%    if sum(spikes(i,idx_int(1)-50:idx_int(1)+50))>=1
%        join=[join i];
%    end
% end
%
% new_spikes=spikes(join,:);
%
% for i=1:size(new_spikes,1)
%     steps=-8:round(14./size(new_spikes,1),1):8;
%     r=find(new_spikes(i,:)==1);
%     new_spikes1(i,r)=steps(i);
% end
%
%
%
% plot(time(join),spikes(join,:),'r.')



for f=4;
%     1:size(output_count,1)
    plot(time,rec_pa(f,:),'k','LineWidth',1.5)
    xlim ([0 400])
    xticks([0:100:400])
    xticklabels ({'-200','-100','0','100','200'})
    box ('off')
    hold on
    for i=1:size(output_count{f,1},1)
        steps=-8:round(13./size(output_count{f,1},1),1):8;
        m=[];
        m=output_count{f,1}(i,:);
        idx=find(m==1);
        m(idx)=steps(i)+2;
        plot(time(idx),m(idx),'r.')
        hold on
    end
    xlim ([0 400])
    ylim ([-6 6])
    close all
end



xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
box ('off')
