clear; close
cohort = [ 2 3 4 5 8 10 11 13 16 17];
cond={'NS';'RS';'PLS'};
clust=struct; out=struct; start=cell(10,3); ending=cell(10,3); yy=cell(10,3); h_up=cell(10,3); s=struct;
rng('default')
gen=(rng);
%
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat');
for iii = 1:length(cohort)
    clearvars -except  cohort cond iii clust s start ending yy out gen h_up avg
    
    for co= 1
%         size(cond,1)
        
        load(strcat('/Users/Carolina/OneDrive - Nexus365/PERI-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
        
        [d]=preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
        
%         [h]=twitch(iii,d,samplerateold,samplerate,co); clear h;
        
                [afilt, bfilt, start, ending, yy, h_up]= startend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up);
        
                [s]=zfiltenv(d,bfilt,afilt,co,iii,s); clear afilt bfilt
              

%                          avg(iii,:)=[mean(s.env_acc{iii,1}(1,h_up{iii,co}),2) (std(s.env_acc{iii,1}(1,h_up{iii,co})'))'];
        
        %         if ~isnan(clust.win(iii,1))
        %             [out]=cluster_intime(clust,s,iii,co,out);
        %         end        

        %  [start,ending,out,yy]=mod_nc(start,ending,co,samplerate,iii,s,yy,out,gen); %%%tremor amplitude change without clustering
    end
    
%                [clust,out]=clustering2(out,iii,clust,start,ending,yy); %% clustering analysis
% %     
%     for co=2
% %         size(cond,1)
%         [out]=mod_c(clust,out,co,iii,s,h_up);   %%%% tremor amplitude change with clustering
%     end
%    
%     
end
clearvars -except out clust

load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat','out');


%%% 1) plots & stats stim ARCs and mediansplit stim ARCs with clustering 
[r,p]=plots_c(out);


%%% 2) plots & stats tremor power during non-stim / superssive stim effects and amplifying stim effects 
[out]=psd_rs_indseg(out); 
% avg pwelch of 5 sec segments amp and sup - 
        %%% pwelch from 5 sec seg of 50 000 non-stim surrogate comes from
        %%% ppx_nseg.m

 %%% 3) plots & stats group stim effects vs non stim effects      
[stats]=group_arcs(out,stats);

 %%% 3) plots & stats group stim effects vs non stim effects after median
 %%% split
run('median_split_group.m');



%%%%%%%%%%%%%%%% NOT USED IN PAPER / EXTRA ANALYSIS

[out]=polyfit(out); %non uniformity of arc's polynomial fits with r2 and f-stats
[out]=polyfit_nc(out);
% % [nc]=plots_nc(out); %plots ARCs without clustering ---- needs fx mod_nc


[out]=cvar(out);
stats=struct;



