clear all;close all
iii=[1 2 3 4 5 6 8 10 11];
for nnn=2;
    clearvars -except iii nnn effect_phase effect_general effect_ind
    for numb=1:length(iii);
        if nnn==1;
            %  load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\PS_PSH\P0',num2str(iii(numb)),'_pha_suffle_phashift.mat'))
            load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/PS_PSH/P0',num2str(iii(numb)),'_pha_suffle_phashift.mat'))
        elseif nnn==2;
            %  load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\PS\P0',num2str(iii(numb)),'_pha_suffle.mat'))
            load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/PS/P0',num2str(iii(numb)),'_pha_suffle.mat'))
        end
        likhood_amp(numb,:)=sum(tt>prctile(tt3,97.5) & tt>0)./sum(~isnan(tt));
        likhood_sup(numb,:)=sum(tt<prctile(tt3,2.5) & tt<0)./sum(~isnan(tt));
%         amp_effect=[]; sup_effect=[];
%         amp_effect=tt(tt>prctile(tt3,97.5) & tt>0);
%         sup_effect=tt(tt<prctile(tt3,2.5) & tt<0);
%         stats(numb,:)=ttest(amp_effect,sup_effect);
        like_effect=sum(tt>prctile(tt3,97.5) & tt>0 | tt<prctile(tt3,2.5) & tt<0) ./ sum(~isnan(tt));
        
        likhood=[likhood_sup ; likhood_amp];
%         bar(likhood')
        alleffect_ind(1,numb)=mean(like_effect)*100;
        ampeffect_ind(1,numb)=mean(likhood_amp(numb,:))*100;
        supeffect_ind(1,numb)=mean(likhood_sup(numb,:))*100;
    end
    figure(nnn)
    bar([mean(alleffect_ind) mean(ampeffect_ind) mean(supeffect_ind) ])  
    effect_phase(nnn,:)=ttest(likhood_amp,likhood_sup)
    effect_general(nnn,1)=ttest(reshape(likhood_amp,1,size(likhood_amp,1)*size(likhood_amp,2)),reshape(likhood_sup,1,size(likhood_sup,1)*size(likhood_sup,2)))
end


    xticklabels({'likelihood of effect','likelihood of amp','likelihood of sup'})
    ylabel('Size effect')