clear all;close all
id=[1 2 3 4 5 6 8 10 11];
for nnn=1:4;
    clearvars -except id numb nnn effect_phase effect_general effect_ind g_effect stats
    for numb=1:length(id);
        
        %%% AMPLITUDE CHANGE
        %%%% stim vs. nno stim
        if nnn==1;
            load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\A_NS\P0',num2str(id(numb)),'_stimnosim.mat'))
            %load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/A_NS/P0',num2str(iii(numb)),'_stimnostim.mat'))
            tt=stimout;
            tt3=nostimout;
            
            %%%% stim vs. phase specific
        elseif nnn==2;
            load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\A_PS\P0',num2str(id(numb)),'_pha_suffle.mat'))
            
            
            %%% FREQUENCY CHANGE
            %%%% stim vs. nno stim
        elseif nnn==3;
            load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\P_PS\P0',num2str(id(numb)),'_pha_suffle_phashift.mat'))
            clear tt3
            load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\P_NS\P0',num2str(id(numb)),'_stimnosim_phashift.mat'))
            tt3=nostimout;
            %load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/P_NS/P0',num2str(iii(numb)),'_stimnosim_phashift.mat'))
            
            %%%% stim vs. phase specific
        elseif nnn==4;
            load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\P_PS\P0',num2str(id(numb)),'_pha_suffle_phashift.mat'))
            %load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/P_PS/P0',num2str(iii(numb)),'_pha_suffle_phashift.mat'))
        end
        
        
        
        likhood_amp(numb,:)=sum(tt>prctile(tt3,97.5) & tt>0)./sum(~isnan(tt));
        likhood_sup(numb,:)=sum(tt<prctile(tt3,2.5) & tt<0)./sum(~isnan(tt));
        %         amp_effect=[]; sup_effect=[];
        amp_effect(numb,:)=mean(tt(tt>prctile(tt3,97.5) & tt>0));
        sup_effect(numb,:)=mean(abs(tt(tt<prctile(tt3,2.5) & tt<0)));
        %         stats(numb,:)=ttest(amp_effect,sup_effect);
        like_effect(numb,:)=sum(tt>prctile(tt3,97.5) & tt>0 | tt<prctile(tt3,2.5) & tt<0) ./ sum(~isnan(tt));
        
        likhood=[likhood_sup ; likhood_amp];
        %         bar(likhood')
        alleffect_ind(1,numb)=mean(like_effect(numb,:))*100;
        ampeffect_ind(1,numb)=mean(likhood_amp(numb,:))*100;
        supeffect_ind(1,numb)=mean(likhood_sup(numb,:))*100;
    end
    
    g_effect(nnn,:)=alleffect_ind;
    [p]=ttest(amp_effect,sup_effect);
    stats(nnn,:)=p;
    %   effect_phase(nnn,:)=ttest(likhood_amp,likhood_sup)
    %   effect_general(nnn,1)=ttest(reshape(likhood_amp,1,size(likhood_amp,1)*size(likhood_amp,2)),reshape(likhood_sup,1,size(likhood_sup,1)*size(likhood_sup,2)))
    
%     figure(nnn)
%     boxplot((like_effect.*100)')
%     xlabel('# patients')
%     box('off')
%     if nnn==1 | nnn==3
%         ylabel('Sig. stim effect (%)')
%     else
%         ylabel('Sig. phase-specific stim effect (%)')
%     end
%     
end

figure()
boxplot((g_effect(1:2,:))')
ylabel('Prob of sig. stim effect (%)')
xticklabels({'Stim vs. No Stim','Phase-specific stim vs. Random Stim'})
box('off')
title('Significant amplitude modulation') 


figure()
boxplot((g_effect(3:4,:))')
ylabel('Prob of sig. stim effect (%)')
xticklabels({'Stim vs. No Stim','Phase-specific stim vs. Random Stim'})
box('off')
title('Significant frequency modulation')

