
% cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')

% clear all
close all
load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/spiral_zsc_4bins.mat'));

ad.metrics{1,1,1}=[ad.metrics{1,1,1}; ad.metrics{1,1,2}];
ad.metrics{2,1,1}=[ad.metrics{2,1,1}; ad.metrics{2,1,2}];


met=cell(3,3);
env_all=cell(3,3);
pol_dat=cell(3,3);

for i=1:3
    for ii=1:3
        met{i,ii}=ad.metrics{i,ii,1}; %%%% met{pt,cond}(n_spirals,n_quadrants,metrics)
        env_all{i,ii}=ad.sig{i,ii,1}(end,:);
        pol_dat{i,ii}=ad.sig{i,ii,1}(1:2,:);
    end
end
env_all{1,1}=[ad.sig{1,1,1}(end,:) ad.sig{1,1,2}(end,:)];
env_all{2,1}=[ad.sig{2,1,1}(end,:) ad.sig{2,1,2}(end,:)];

pol_dat{1,1}=[ad.sig{1,1,1}(1:2,:) ad.sig{1,1,2}(1:2,:)];
pol_dat{2,1}=[ad.sig{2,1,1}(1:2,:) ad.sig{2,1,2}(1:2,:)];


me=1;
for ss= 1:size(met,1)
    for cc=1:size(met,2)
        met{ss,cc}(:,:,me)
    end
end

%%% Figure paper 4.C. ss= patients cc=condition

for ss=3 
    %     1:size(ad.metrics,1)
    for cc=1:size(ad.metrics,2)
        % %         for me=1
        % %             [P,ANOVATAB,STATS] = anova1(met{ss,cc}(:,:,me));
        % %             close
        % %             [c,m,h,nms] = multcompare(STATS);
        % %             close
        % %             results{ss,cc}=find(c(:,end)<0.05);
        % %
        % %             [val,hh]=ttest([met{ss,cc}(:,1,me) ; met{ss,cc}(:,3,me)], [met{ss,cc}(:,2,me) ; met{ss,cc}(:,4,me)]);
        % %             test13{ss,cc}=val; clear hh val
        % %         end
        
        limits(1,:)=[min(env_all{ss,cc}) max(env_all{ss,cc})-(0.5*max(env_all{ss,cc}))];
        limits(2,:)=[min(horzcat(env_all{ss,:})) (max(horzcat(env_all{ss,:})))/mean(horzcat(env_all{ss,:}))];
        %         limits(2,:)=[min(horzcat(env_all{ss,:})) (max(horzcat(env_all{ss,:})))];
        
        
        % %                         f1=figure(ss)
        % %                         senv=smoothdata(env_all{ss,cc},'movmean',samplerate2/2);
        % %                         subplot(1,3,cc)
        % %                         n=polarscatter(pol_dat{ss,cc}(1,:),pol_dat{ss,cc}(2,:),8,senv,'filled')
        % %                         colorbar
        % %                         colormap winter
        % %                         rticks([0 400])
        % %                         rticklabels([])
        % %                         caxis([limits(2,1) limits(2,2)])
        % %                         f1.OuterPosition= [100,1000,1100,600];
        % %                         set(f1,'color','w');
        % %
        % %
        % % me=1;
        % %
        % %
        % % if ~isempty(results{1,me}{ss,cc})
        % %     dum=results{1,me}{ss,cc};
        % %     name=strcat('Q',num2str(c(dum,1)),'-Q',num2str(c(dum,2)))
        % %     title(name,'FontSize',8)
        % % else
        % %     title('none')
        % % end
        
        
        CT=cbrewer('seq','YlOrRd',9);
        CT(1,:)=[];
       
        for i=1:size(ad.spirals{ss,cc},1)
            pdat=cell2mat(ad.spirals{ss,cc}(i));
            senv=smoothdata(pdat(3,:),'movmean',samplerate2/2);
            f1=figure();
            
            %                 if size(ad.spirals{ss,cc},1)>5
            %                 subplot(2,5,i)
            %                 else
            %                 subplot(1,(size(ad.spirals{ss,cc},1)),i)
            %                 end
            ll= 350;
            polarscatter(pdat(1,:),pdat(2,:),12,senv,'filled')
            hold on
            y=linspace(-ll,ll,100);
            x=repmat(degtorad(90),1,100);
            polarplot(x,y,'k','LineWidth',1.5)
            hold on
            y=linspace(-ll,ll,100);
            x=repmat(degtorad(180),1,100);
            polarplot(x,y,'k','LineWidth',1.5)
            
            colorbar 
            colormap (CT)
            rlim([0 ll])
            rticks([0 ll])
            rticklabels([])
            caxis([limits(1,1) limits(1,2)])
            set(gca,'FontSize',14)
            f1.OuterPosition= [1,100,300,300];
            set(f1,'color','w');
            
        end
        %         f1.OuterPosition= [1,100,300,300];
        %             set(f1,'color','w');
        %                     close all
        clear limits pdata senv f1
        
    end
end

% % for me=1
% %     for ss= 1:size(met,1)
% %         for cc=1:size(met,2)
% %             datt(cc,:)=nanmean(met{ss,cc}(:,:,me));
% %             dst(cc,:)=nanstd(met{ss,cc}(:,:,me));
% %             
% %             %%%realign to max amplitude
% %             %             norm=mean(met{ss,cc}(:,:,1));
% %             %             amax=find(norm==max(norm));
% %             %             if amax==1
% %             %                 datt1=datt;
% %             %             else
% %             %                 datt1=[datt(1,amax:end) datt(1,1:amax-1)];
% %             %             end
% %             
% %             datt1(cc,:)=datt(cc,:)./max(datt(cc,:));
% %             
% %             f2=figure(ss)
% %             subplot(size(met,2),1,cc)
% %             bar([datt(cc,:)],'EdgeColor','none','FaceColor',[0.8 0.5 1],'FaceAlpha',0.5)
% %             hold on
% %             errorbar([1:numb_bins],datt(cc,:),dst(cc,:),'.','Color',[0.8 0.5 1],'LineWidth',2)
% %             % ylim([0 500])
% %             % yticks(0:100:500)
% %             if numb_bins==4
% %                 xticklabels({'0','90','180','270'})
% %             elseif numb_bins==6
% %                 xticklabels({'0','60','120','180','240','300'})
% %             else xticklabels({'0','45','90','135','180','225','270','315'})
% %             end
% %             %
% %             box('off')
% %             xlabel('degrees')
% %             name=ad.metricnames{5+me, 1};
% %             ylabel(name)
% %             f2.OuterPosition= [1,100,200,500];
% %             set(f2,'color','w');
% %             
% %         end
% %         f1=figure(size(met,1)+1)
% %         bar(datt','EdgeColor','none')
% %         box('off')
% %         legend({'NS','HFS','PLS'})
% %         legend('boxoff')
% %         f1.OuterPosition= [1,100,500,200];
% %         set(f1,'color','w');
% %         
% %         %
% %         
% %         close all
% %     end
% %     
% %     
% % end

