clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load('data_SUA_SNR.mat') ; 
SNR.ctx=Ecog_all;
SNR.sua=data_all;
clear data_all

srn=1000;
time=0:0.001:((size(SNR.ctx(1,:),2)-1)./srn);

n=[17.5:7.5:32.5];
freq=cell(1,1);
for i=1:length(n)
    freq{1,i}=[n(i)-5 n(i)+5];
end

idx=[];
for i=1:size(SNR.sua,1)
    if ~isempty(SNR.sua{i,1})
        idx=[idx i];
    end
end

for t=1:size(freq,2)
    for i=1:length(idx)
        for  j=1:size(SNR.sua{idx(i),1},1)
            data=SNR.sua{idx(i),1}(j,:);
            data_all{i,1}(j,:)=data;
            data_ones{i,j}=find(data==1);
            data_one=find(data==1);
            [b,a]=butter(2,[freq{1,t}(1)/(0.5*srn) freq{1,t}(end)/(0.5*srn)],'bandpass');
            Ecogfiltered=filtfilt(b,a,SNR.ctx(idx(i),:));
            ecog_units(i,t,:)=Ecogfiltered;
            hp=wrapToPi(angle(hilbert(Ecogfiltered)));
            ang{i,j}(:,t)=hp(data_one); clear hp;
            vect_length{i,j}(:,t)=circ_r(ang{i,j}(:,t));
            ray_test{i,j}(:,t)=circ_rtest(ang{i,j}(:,t));
            if circ_rtest(ang{i,j}(:,t))<0.05
                r{i,j}(:,t)=1;
            else
                r{i,j}(:,t)=0;
            end
        end
    end
    clear b a
end



for i=1:size(ray_test,1)
    stay=[];
    for j=1:size(ray_test,2)
        if sum(r{i,j})>0
            stay=[stay j];
        end
    end
    stat_d{i,1}=stay;
end



for i =1:size(stat_d,1)
    for ii= 1:size(stat_d{i,1},2)
        if find(cell2mat(vect_length(stat_d{i,1}(ii)))==max(cell2mat(vect_length(stat_d{i,1}(ii)))))==1
            ecogbf_match(i,:)=ecog_units(i,1,:);
        elseif find(cell2mat(vect_length(stat_d{i,1}(ii)))==max(cell2mat(vect_length(stat_d{i,1}(ii)))))==2
            ecogbf_match(i,:)=ecog_units(i,2,:);
        elseif find(cell2mat(vect_length(stat_d{i,1}(ii)))==max(cell2mat(vect_length(stat_d{i,1}(ii)))))==3
            ecogbf_match(i,:)=ecog_units(i,3,:);
        end
    end
    units_match{i,1}=data_all{i,1}(stat_d{i,1}(1:end),:);
end



clearvars -except units_match ecogbf_match stat_d data_all data_ones srn time SNR
% save 'NEW_SNR_cycle'


