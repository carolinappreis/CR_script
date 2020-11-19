clear all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
% cd ('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
load('probeselect_SUA_SNR.mat')


srn=1000;
for ii=1:size(data_all,1)
    clearvars -except A lesion ii j r region_pl region_npl region_spl region_snpl output_pa output_npa spikerate_all ISI_all data_all Ecog_all srn time
    for  j=1:size(data_all{ii,:},1)
        
        data=data_all{ii,1}(j,:);
        data_ones=find(data==1);
        [b,a]=butter(2,[15/(0.5*srn) 35/(0.5*srn)],'bandpass');
        Ecogfiltered=filtfilt(b,a,Ecog_all(ii,:));
        env=abs(hilbert(Ecogfiltered));
        spkrate_1=[];
        for i =1:srn:(length(data)-srn);
            spkrate_1=[spkrate_1 numel(find(data(i:i+srn)==1))];
            ISI_all{ii,1}(1,j)=mean(diff(time(data==1)));
            %                 psi_all{ii,1}(j,:)=abs(mean(find(ang(data_ones)));
        end
        
        spikerate_all{ii,1}(1,j)=mean(nonzeros(spkrate_1));
        if mean(nonzeros(spkrate_1))-std(nonzeros(spkrate_1))>=13 && mean(nonzeros(spkrate_1))+std(nonzeros(spkrate_1))<35
            %             cd('C:\Users\creis\Documents\GitHub\CR_script\A4_Thal\code')
            cd('/Users/Carolina/Documents/GitHub/CR_script/A4_Thal/code')
            onset1=bursts(env);
            onset=bursts_aligned(env,Ecogfiltered);
            
            
            data_g=smoothdata(data,'gaussian',25);
            
            
            for nb=1:size(onset1,1)
                for jj=1:size(onset{nb,1},2)
                    if onset{nb,1}(jj)>200 && onset{nb,1}(jj)+200<length(data_g)
                        output_pa{j,:}(jj,:)=data_g(onset{nb,1}(jj)-200:onset{nb,1}(jj)+200);
                        output_npa{j,:}(jj,:)= data_g(onset1{nb,1}(jj)-200:onset1{nb,1}(jj)+200);
                    end
                end
                rec_pa{nb,1}(j,:)=sum(output_pa{j,:},1); rec_pa{nb,1}=rec_pa{nb,1}(any(rec_pa{nb,1},2),:);
                rec_npa{nb,1}(j,:)=sum(output_npa{j,:},1); rec_npa{nb,1}=rec_npa{nb,1}(any(rec_npa{nb,1},2),:);
                
                region_pl{ii,nb}=zscore(mean(rec_pa{nb,1},1));
                region_spl{ii,nb}=zscore(std(rec_pa{nb,1})./sqrt(size(rec_pa{nb,1},1)));
                region_npl{ii,nb}=zscore(mean(rec_npa{nb,1},1));
                region_snpl{ii,nb}=zscore(std(rec_npa{nb,1})./sqrt(size(rec_npa{nb,1},1)));
            end
        end
    end
end

small_b=mean(vertcat(region_pl{:,1}),1);
long_b=mean(vertcat(region_pl{:,2}),1);
% plot(small_b)
hold on
plot(long_b)

% rats_id=str2num(A(lesion,4:6));
% unNum     = unique(rats_id(:));
% [n,bin]   = histc(rats_id(:),unNum);



%
% x 13 14 23