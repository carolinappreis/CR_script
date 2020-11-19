%saving subcortical data (median beta change) to plot on files:
%'subplot_long';'subplot_short';'coherent regions'

% clear all
% cd ('\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A1_Thal\mat')
% load 'SNR_complete_new'
% load 'BZ_complete_new'

med_thal1=cell2mat([med_thal_cell(:,1)]);
med_thal2=cell2mat([med_thal_cell(:,2)]);
med_thal1_s=prctile(med_thal1,[25 50 75]);
med_thal2_s=prctile(med_thal2,[25 50 75]);
med_thal_s=cell(2,1);
med_thal_s{1}=med_thal1_s; %short bursts
med_thal_s{2}=med_thal2_s; %long bursts

med_probe=[];
for j=1:2;
    for k=1:size(med_thal_cell,1);
        if ~isempty (med_thal_cell{k,j}) || size(med_thal_cell{k,j},1)~=1  
            med_probe(j,k,1:1001)= median(med_thal_cell{k,j});
        end
    end
end

short_thal=cell(1,1);
long_thal=cell(1,1);
m=1;
for i =1:size(med_thal_cell,1)
    for l=1:size(med_thal_cell{i,1},1)
        short_thal{m,1}(l,:)=med_thal_cell{i,1}(l,:);
    end
    short_thal=short_thal(~cellfun('isempty',short_thal));
    
    for l=1:size(med_thal_cell{i,2},1)
        long_thal{m,1}(l,:)=med_thal_cell{i,2}(l,:);
    end
    long_thal=long_thal(~cellfun('isempty',long_thal));
    
    m=m+1;
end
surr1=surr(~cellfun('isempty',surr));
surr2=cell2mat(surr1);

med_short_thal=[];
med_long_thal=[];
surr=[];

for k=1:size(short_thal,1);
    if size(short_thal{k},1)~=1
        med_short_thal(k,1:1001)= median(short_thal{k});
        med_long_thal(k,1:1001)= median(long_thal{k});
        surr(k,1:1001)= median(surr1{k});
        
    else
        med_short_thal(k,1:1001)= short_thal{k};
        med_long_thal(k,1:1001)= long_thal{k};
        surr(k,1:1001)= surr1{k};
    end
    
end

med_contacts(1,1:(size(med_short_thal,1)),1:1001)=med_short_thal;
med_contacts(2,1:(size(med_short_thal,1)),1:1001)=med_long_thal;

all_contacts(1,1:(size(med_thal1,1)),1:1001)=med_thal1;
all_contacts(2,1:(size(med_thal1,1)),1:1001)=med_thal2;

surr_s=prctile(surr2,[25 50 75]);

plot(med_thal_s{2}(2,:))
hold on
plot(surr_s(2,:))
plot(med_thal_s{1}(2,:))

% to save non coherent saveas BZ_nc_fig change to cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat/all_regions')

clearvars -except med_contacts med_ctx med_thal_cell med_thal_s surr surr_s onset all_contacts surr2
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')
% save 'PC_nc_fig'




