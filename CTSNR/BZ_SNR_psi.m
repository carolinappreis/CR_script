
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

combo=[1 2 3 ; 2 1 3; 3 1 2];

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean','black','grey');
col={[grey ; black];[squash ; blood];[sky ; aegean]};

for iii=1:size(combo,1)
for co1=combo(iii,1)
    clearvars -except region Ecog BZ SNr col co1 iii combo
    dat1=eval(region{1,co1});
    for co2=combo(iii,2)
        clearvars -except region Ecog BZ SNr col co1 co2 dat1 iii combo 
        dat2=eval(region{1,co2});
        for co3=combo(iii,3)
            clearvars -except region Ecog BZ SNr col co1 co2 co3 dat1 dat2 iii combo
            dat3=eval(region{1,co3});
            if co2~=co3
                [data,baseline]=change_psi(dat1,dat2,dat3,co1,co2,co3);
                el=500;
                time=1:el+1;
                for ii=1:2
                    fig=figure;
                    for on=1
%                         on=1:2
                        clear A; A=squeeze(data{1,ii}(on,:,:));
                        if size(A,1)>1
                            val=mean(A);
                        else
                            val=A;
                        end
%                         y2=zscore(val);
%                         y1=zscore(val+std(A)./sqrt(size(A,1)));
%                         y3=zscore(val-std(A)./sqrt(size(A,1)));
                        y2=(val);
                        y1=(val+std(A)./sqrt(size(A,1)));
                        y3=(val-std(A)./sqrt(size(A,1)));
                        %cl= eval(sprintf('col{co%d,1}',ii));
                        cl=(col{co2,1}(on,:)+col{co3,1}(on,:))./2;
                        plot(time, y2,'LineStyle','-','LineWidth',1.5,'Color',cl)
                        hold on
                        patch([time fliplr(time)], [y1 fliplr(y2)],cl,'FaceAlpha',0.2,'EdgeColor','none')
                        patch([time fliplr(time)], [y2 fliplr(y3)],cl,'FaceAlpha',0.2,'EdgeColor','none')
                        
                        site=([[2+0.3 2+0.3 2+0.4 2+0.4];[2+0.6 2+0.6 2+0.7 2+0.7]]);
                        st=NaN(1,length(time));
                        clear B; B=data{1,ii+2};
                        hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
                        beg=[];
                        beg=find(st(1,:)<0.05 & st2(1,:)~=0);
                        sig_rise_all=[];
                        if ~isempty(beg)
                            sig_rise_all=[beg(1) beg(end)];
                        end
                        clear st st2
                        if ~isempty (beg)
                            patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],site(on,:),col{co1,1}(on,:),'EdgeColor','none')
                            clear beg
                        end
                    end
                    hold on
                    clear cl
                    
                    %         xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)
                    xline(200,'--','Color',[0.5 0.5 0.5],'LineWidth',2)
                    box('off')
                    xlim ([0 500])
                    ylim([-3 3])
                    xticks([0:100:500])
                    xticklabels ({'-200','-100','0','100','200','300'})
                    fig.Units = 'centimeters';
                    fig.OuterPosition= [10, 10, 8, 10];
                    fig.Color='w';
                    set(gca,'FontSize',12)
                     if ii==1
                         ylabel('PSI within bursts (zscore)');
                     else
                         ylabel('PSI across bursts (zscore)');
                     end
                    xlabel('Time (ms)');
                    title([eval(sprintf('region{1,co%d}',2)) eval(sprintf('region{1,co%d}',3))])
                    
                    str = ['psi:',num2str(baseline(ii))];
                    t = text(-1.5,-2.5,str);
                    s = t.FontSize;
                    t.FontSize = 12;
                end
            end
        end
    end
end
end
