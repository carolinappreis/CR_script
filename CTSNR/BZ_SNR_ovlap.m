

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

% % clear all
% % close all
% % cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/final_mats')
% % load('BZ_SNR_cohfilts.mat');
% %
% % for i=1:size(BZ,1)
% %     Ecog{i,1}=BZ{i,1}(1,:);
% %     SNr{i,1}=SNR{i,1}(2:end,:);
% %     BZ{i,1}=BZ{i,1}(2:end,:);
% % end

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



dur_ovl={}; all={}; prct_all={}; all_shuff={}; dur_shuff={};
comb={[1 2];[1 3];[2 3];[1 2 3]};
for iii=1:size(data,2)
    for combi=1:size(comb,1)
        m=0; n=0;
        for co1=1:size(data{comb{combi,1}(1),iii},1)
            clear ref1;
            ref1=data{comb{combi,1}(1),iii}(co1,:);
            for co2=1:size(data{comb{combi,1}(2),iii},1)
                m=m+1;
                clear ref2;
                ref2=data{comb{combi,1}(2),iii}(co2,:);
                
                if combi~=4
                    dum=[ref1;ref2];
                    [number_ovl,dur_ovl]=overlap(dum,combi);
                    ref_b=find(dum(1,:)==1);
                    [dur_chance_ovl]=chance_ovl(dum,ref_b,combi);
                    dur_shuff{combi,iii}(m,:)=dur_chance_ovl;
                    dur_ovlap{combi,iii}(1,m)=dur_ovl;
                    
                else
                    for co3=1:size(data{comb{combi,1}(3),iii},1)
                        clear ref3
                        ref3=data{comb{combi,1}(3),iii}(co3,:);
                        n=n+1;
                        dum=[ref1;ref2;ref3];
                        [number_ovl,dur_ovl]=overlap(dum,combi);
                        all{1,iii}(1,n)=dur_ovl;
                        
                        [dur_chance_ovl]=chance_ovl(dum,ref_b,combi);
                        all_shuff{1,iii}(n,:)=dur_chance_ovl;
                    end
                end
                clear dur_ovl
            end
        end
    end
end

comb={[1 2];[1 3];[2 3];[1 2 3]};
for i=1:3
    dur_b1{i,1}=horzcat(dur_b{i,:});
    dur_ovlap1{i,1}=horzcat(dur_ovlap{i,:});
    dur_shuff1{i,1}=horzcat(dur_shuff{i,:});
end
all1{1,1}=horzcat(all{1,:});
all_shuff1{1,1}=horzcat(all_shuff{1,:});


pr=1;
for ii=1:3
    p1(ii)=nanmean(dur_b1{ii,pr});
    int1(ii)=nanmean(dur_ovlap1{ii,pr});
end

total_b=[p1(1)+p1(2)-int1(1)-nanmean(all1{1,pr})+p1(3)-int1(2)-int1(3)-nanmean(all1{1,pr})]/1000;

int2=int1-nanmean(all1{1,pr});
r(1)=p1(1)-(int2(1)+int2(2));
r(2)=p1(2)-(int2(1)+int2(3));
r(3)=p1(3)-(int2(2)+int2(3));

int1(1,4)=nanmean(all1{1,pr});

bursts=[r int1]./1000;
% total_n_b=nanmean(none)./1000;
total_n_b=100-sum(bursts);

Z=[p1 int1]./1000;

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean','grey');
color_b={grey,blood,aegean};
fig=figure;
venn(Z,'FaceColor',color_b,'FaceAlpha',{0.4,0.4,0.4},'EdgeColor','none')
legend({'Ecog'; 'BZ'; 'SNr'});
legend('boxoff')

prcts=round([total_n_b bursts],1);

str = [num2str(prcts(1)),'%'];
t = text(-3,6,str);
s = t.FontSize;
t.FontSize = 14;

str = [num2str(prcts(2)),'%'];
t = text(-0.5,0,str);
s = t.FontSize;
t.FontSize = 14;

str = [num2str(prcts(3)),'%'];
t = text(3.5,0,str);
s = t.FontSize;
t.FontSize = 14;

str = [num2str(prcts(4)),'%'];
t = text(1.5,4,str);
s = t.FontSize;
t.FontSize = 14;

str = [num2str(prcts(5)),'%'];
t = text(2,0,str);
s = t.FontSize;
t.FontSize = 12;

str = [num2str(prcts(6)),'%'];
t = text(0.5,2,str);
s = t.FontSize;
t.FontSize = 12;

str = [num2str(prcts(7)),'%'];
t = text(2.5,2,str);
s = t.FontSize;
t.FontSize = 12;

str = [num2str(prcts(8)),'%'];
t = text(2,1.5,str);
s = t.FontSize;
t.FontSize = 12;

set(gca,'visible','off')
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 16, 18];




figure
bar(Z)