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


load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean','black','grey');
col={[grey ; black];[squash ; blood];[sky ; aegean]};

for co1=1:size(region,2)
    clearvars -except region Ecog BZ SNr col co1
    dat1=eval(region{1,co1});
    for co2=1:size(region,2)
        clearvars -except region Ecog BZ SNr col co1 co2 dat1
        dat2=eval(region{1,co2});
        [data]=change_amp(dat1,dat2,co1,co2);
        
        if ~isempty(data)
            el=500;
            time=1:el+1;
            for ii=1:size(data,2)-1
                fig=figure;
                for on=1:2
                    clear A; A=squeeze(data{1,ii}(on,:,:));
                    y2=100.*(mean(A));
                    y1=100.*(mean(A)+std(A)./sqrt(size(A,1)));
                    y3=100.*(mean(A)-std(A)./sqrt(size(A,1)));
                    cl= eval(sprintf('col{co%d,1}',ii));
                    plot(time, y2,'LineStyle','-','LineWidth',1.5,'Color',cl(on,:))
                    hold on
                    patch([time fliplr(time)], [y1 fliplr(y2)],cl(on,:),'FaceAlpha',0.2,'EdgeColor','none')
                    patch([time fliplr(time)], [y2 fliplr(y3)],cl(on,:),'FaceAlpha',0.2,'EdgeColor','none')
                    if co1~=co2 && ii==2
                        site=([[15+1 15+1 15+2 15+2];[15+4 15+4 15+5 15+5]]);
                        st=NaN(1,length(time));
                        clear B; B=data{1,3};
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
                end
                                    %         xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)
                    xline(200,'--','Color',[0.5 0.5 0.5],'LineWidth',2)
                    box('off')
                    xlim ([0 500])
                    if max(mean(A).*100)>80
                        ylim ([-150 150])
                    else
                        ylim ([-20 20])
                    end
                    xticks([0:100:500])
                    xticklabels ({'-200','-100','0','100','200','300'})
                    fig.Units = 'centimeters';
                    fig.OuterPosition= [10, 10, 8, 10];
                    fig.Color='w';
                    set(gca,'FontSize',12)
                    ylabel('change in beta amplitude (%)');
                    xlabel('Time (ms)');
                    title(eval(sprintf('region{1,co%d}',ii)))
            end
        end
    end
end

