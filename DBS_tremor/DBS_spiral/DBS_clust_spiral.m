clear; close
cohort = [ 1 3 4 6];
cond={'NS';'RS'};
pc_t=cell(size(cohort,2),size(cond,1));
start=cell(10,3); ending=cell(10,3); yy=cell(10,3); out=struct; clust=struct;

for iii =  1
%     :length(cohort)
    clearvars -except  cohort cond iii pc_t start ending yy out clust
    for co=1:size(cond,1)
        
        load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
        
        [out]=s_preprocess(SmrData, out, iii, co);
  
        [pc_t,out]=stitch_xall(out,iii,co,pc_t);
        
    end
    
    [clust]= cluster_s(pc_t,iii);
    
    time_dom(out,clust,iii)

    
    close all
end

