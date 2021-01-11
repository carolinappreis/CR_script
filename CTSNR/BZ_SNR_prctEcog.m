
% clear all
% close all
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/final_mats')
% load('BZ_SNR_cohfilts.mat');
%
% for i=1:size(BZ,1)
%     Ecog{i,1}=BZ{i,1}(1,:);
%     SNr{i,1}=SNR{i,1}(2:end,:);
%     BZ{i,1}=BZ{i,1}(2:end,:);
% end

% % clear all
% % close all
% % cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/final_mats')
% % load('BZ_SNR_filtcoh_nocrit.mat');
% % BZ=new_data{1,1}(:,2:end,:);
% % SNr=new_data{1,2}(:,2:end,:);
% % Ecog=new_data{1,2}(:,1,:);
% % Ecog1=new_data{1,1}(:,1,:);

clear all
close all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/final_mats')
load('BZ_SNR.mat');

for i=1:size(rat_labels,1)
    Ecog{i,1}=BZ_bua{i,1}(1,:);
    SNr{i,1}=SNr_bua{i,1}(2:end,:);
    BZ{i,1}=BZ_bua{i,1}(2:end,:);
end

% % clear all
% % close all
% % cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/final_mats')
% % load('BZ_coh.mat');
% % load('SNr_coh.mat');
% %
% % r=0;
% % for i=size(BZ_raw_coh,1)-1:size(BZ_raw_coh,1)
% %     r=r+1;
% %     Ecog{r,1}=BZ_raw_coh{i,1}(1,:);
% %     SNr{r,1}=SNR_raw_coh{i,1}(2:end,:);
% %     BZ{r,1}=BZ_raw_coh{i,1}(2:end,:);
% % end

clearvars -except Ecog BZ SNr
region={'Ecog' 'BZ' 'SNr'};


samprate=1000;
for i=1:size(Ecog,1)
    [Pxx_ind,F_i]=pwelch(Ecog(i,:),samprate,[],samprate,samprate);
    frange=find(F_i==15):find(F_i==35);
    Pxx_ind_beta=Pxx_ind(frange);
    filtrange(i,1)=14+find(Pxx_ind_beta==max(Pxx_ind_beta));
end

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean','black','grey');
col={[grey ; black];[squash ; blood];[sky ; aegean]};
data={}; nr_b={}; dur_b={}; dur_nb={};



for co1=1:size(region,2)
    clear dat1
    dat1=eval(region{1,co1});
    for iii=1:size(Ecog,1)
        clearvars -except region Ecog BZ SNr iii co1 data nr_b dur_b dur_nb filtrange dat1
        dat=dat1{iii,1};
        [data,nr_b,dur_b,dur_nb]=disc_bursts(dat,data,co1,iii,nr_b,dur_b,dur_nb,filtrange);
    end
end



ovl_4comb={};  ovl_4comb_shuff={};
comb={[1 2];[1 3];[2 3];[1 2 3]};

for combi=1:size(comb,1)
    m=0;n=0;
    for iii=1:size(data,2)
        for co1=1:size(data{comb{combi,1}(1),iii},1)
            clear ref1;
            ref1=data{comb{combi,1}(1),iii}(co1,:);
            for co2=1:size(data{comb{combi,1}(2),iii},1)
                m=m+1;
                clear ref2 ref_b
                ref2=data{comb{combi,1}(2),iii}(co2,:);
                
                if combi~=4
                    dum=[ref1;ref2];
                    ref_b=find(dum(1,:)==1);
                    [number_ovl,dur_ovl]=overlap(dum,combi);
                    ovl_4comb{combi,1}(m,:)=((([(length(ref_b)-dur_ovl) dur_ovl])./length(ref_b)).*100);
                    [dur_chance_ovl]=chance_ovl(dum,ref_b,combi);
                    ovl_4comb_shuff{combi,1}(m,:)=dur_chance_ovl;
                    clear ref3 clear dum dur_ovl ref_b
                else
                    for co3=1:size(data{comb{combi,1}(3),iii},1)
                        n=n+1;
                        clear ref3
                        ref3=data{comb{combi,1}(3),iii}(co3,:);
                        dum=[ref1;ref2;ref3];
                        ref_b=find(dum(1,:)==1);
                        [number_ovl,dur_ovl]=overlap(dum,combi);
                        ovl_4comb{combi,1}(n,:)=((([(length(ref_b)-dur_ovl) dur_ovl])./length(ref_b)).*100);
                        [dur_chance_ovl]=chance_ovl(dum,ref_b,combi);
                        ovl_4comb_shuff{combi,1}(n,:)=dur_chance_ovl;
                        clear ref3 clear dum dur_ovl ref_b
                    end
                end
            end
        end
    end
end


clearvars -except ovl_4comb  ovl_4comb_shuff comb region

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean','grey');
color_b=[blood;aegean;grey];

m=0;
for i=1:size(comb,1)
    if i~=3
        fig=figure;
        m=m+1;
        data=ovl_4comb{i,1}(:,2);
        surr=ovl_4comb_shuff{i,1}(:,2);
        [p,h]= ttest(data,surr);
        stats(i,:)=h;  clear h
        bar([1:2],[mean(data) mean(surr)],'EdgeColor','none','FaceColor',color_b(m,:),'FaceAlpha',0.5)
        hold on
        plot(1,data,'.','Color',[0.5 0.5 0.5])
        plot(2,surr,'.','Color',[0.5 0.5 0.5])
        if stats(i,:)<0.05
            plot(1.5,50,'k*','LineWidth',1)
        end
        ylim([0 50])
        xticklabels({'OVL data','OVL surr'})
        ylabel('Total time EcogB (%)')
        box('off')
        set(gca,'FontSize',12)
        fig.Units='centimeters';
        fig.OuterPosition= [10, 10, 12, 10];
        set(fig,'color','w');
        if i~=size(comb,1)
            title([eval(sprintf('region{1,(comb{i,1}(2))}'))])
            
        else
            title([eval(sprintf('region{1,(comb{i,1}(2))}')) eval(sprintf('region{1,(comb{i,1}(3))}'))])
        end
        clear data surr p h
    end
end

re=1;
    color_b={{blood,aegean};{grey, aegean};{grey,blood}};
    labels=[{'BZ', 'SNr'} ;{'ECoG','SNr'}; {'ECoG', 'BZ'}];
    Z= [mean(ovl_4comb{1,1}(:,2)) mean(ovl_4comb{2,1}(:,2)) mean(ovl_4comb{4,1}(:,2))];
    fig=figure;
    col=color_b{re,1};
    venn(Z,'FaceColor',col,'FaceAlpha',{0.4,0.4},'EdgeColor','none')
    legend([{labels{re,1}} ; {labels{re,2}}]);
    legend('boxoff')
    
    ovlap=[mean(ovl_4comb{1,1}(:,2))-mean(ovl_4comb{4,1}(:,2))  mean(ovl_4comb{2,1}(:,2))-mean(ovl_4comb{4,1}(:,2)) mean(ovl_4comb{4,1}(:,2))];
    no_ovlap=100-sum(ovlap);
    r=round([no_ovlap ovlap],1);
    
    str = [num2str(r(1)),'%'];
    t = text(-3.5,2,str);
    s = t.FontSize;
    t.FontSize = 14;
    
    str = [num2str(r(2)),'%'];
    t = text(-0.5,0,str);
    s = t.FontSize;
    t.FontSize = 14;
    
    str = [num2str(r(3)),'%'];
    t = text(3.5,0,str);
    s = t.FontSize;
    t.FontSize = 14;
    
    str = [num2str(r(4)),'%'];
    t = text(1.5,0,str);
    s = t.FontSize;
    t.FontSize = 14;
    
    set(gca,'visible','off')
    fig.Units = 'centimeters';
    fig.OuterPosition= [10, 10, 16, 12];


