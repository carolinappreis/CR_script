clear all
ad=struct; ad.metrics=cell(1,1); ad.spirals=cell(1,1); ad.sig=cell(1,1); ad.signames={'data1';'data2';'time';'envelope'};ad.metricnames={'subj';'cond';'trials';'nspirals';'nquadrants';'mean env';'peak2peak';'variance';'freq'};

close all
cond={'NS';'HF';'C'};
cohort=[1 3 4 6];
pca1=0; %%% use the sqrt or the pca as the signal to extract the envelope from
norm=1; %%% norm=1:sqrt(x2+y2) is zscored across condition ;


for iii=2:4
    clearvars -except ad iii cohort cond pca1 norm gg avg_freq
    for co=1:3
        cn=1;
        if (co==1 && iii==2 | iii==3)
            ntrial=[1 2];
        else
            ntrial=1;
        end
        alli=0;
        for trial=1:length(ntrial)
            load(strcat('/Users/Carolina/Desktop/clean_SH_spirals/P0',num2str(cohort(iii)),'_clean_',num2str(cond{co,1}),num2str(trial),'_SH.mat'));
            load(strcat('/Users/Carolina/Desktop/clean_SH_spirals/zscore_acr_cond.mat'))

            centre=[535 361];
            signal1(1,:)=signal(1,:)-centre(1);
            signal1(2,:)=signal(2,:)-centre(2);

            if norm==0
                signal1(3,:)= sqrt((signal1(1,:).^2)+(signal1(2,:).^2));
            else
                signal1(3,:)=signal_all{iii-1}(map{iii-1,co}(trial,1):map{iii-1,co}(trial,end));
            end

            si=signal1(3,:);
            %%% find peak tremor
            [peak,Pxx,F]=dindfpeak(si,samplerate2);
            if (peak-2)>=1
                [a,b]=butter(2,[(peak-2)/(0.5*samplerate2) (peak+2)/(0.5*samplerate2)],'bandpass'); %15
            else
                [a,b]=butter(2,[(1)/(0.5*samplerate2) (peak+2)/(0.5*samplerate2)],'bandpass'); %15
            end

            %             [a,b]=  butter(2, [2/(0.5*samplerate2) 10/(0.5*samplerate2)], 'bandpass'); %15


            if pca1==0
                filt_cs=(filtfilt(a,b,signal1(3,:)));

            else

                for i=1:2
                    to_pc(i,:)=(filtfilt(a,b,signal1(i,:)));
                end

                dx = [to_pc(1,:); to_pc(2,:)];
                [pc, score, latent, tsquare, explained] = pca(dx');
                filt_cs= to_pc(1,:).*pc(1,1)+ to_pc(2,:).*pc(2,1);
                clear a b to_pc pc comb

            end

            %%% filt & var
            env=abs(hilbert(filt_cs));
            phase=angle(hilbert(filt_cs));
            freq=(smooth((100/(2*pi))*diff(unwrap(phase)),50))';


            [px,py]=cart2pol(signal1(1,:),signal1(2,:));

            run('time_split.m')

            %             figure(2)
            %             polarplot(px,py)
            %             hold on
            %             figure(3)
            %             plot(tempo,filt_cs)
            %             hold on


            for ns=1:size(t_spi{iii,trial,co}(:,1),1)
                alli=alli+1;

                % % check power/space
                % %                 [mama,maid]=sort(env,'descend');
                % %                 figure(1)
                % %                 if size(t_spi{iii,trial,co}(:,1),1)>5
                % %                 subplot(2,(size(t_spi{iii,trial,co}(:,1),1))/2,ns)
                % %                 else
                % %                 subplot(1,(size(t_spi{iii,trial,co}(:,1),1)),ns)
                % %                 end
                % %                 plot3(tempo,signal1(1,:),signal1(2,:))
                % %                 hold on
                % %                 plot3(tempo(maid(1:length(maid)/2)),signal1(1,(maid(1:length(maid)/2))),signal1(2,(maid(1:length(maid)/2))),'r.')
                % %                 clear mama maid
                % %                 [mama,maid]=sort(env,'ascend');
                % %                 plot3(tempo(maid(1:length(maid)/2)),signal1(1,(maid(1:length(maid)/2))),signal1(2,(maid(1:length(maid)/2))),'k.')
                % %                 xlim([tempo(t_spi{iii,trial,co}(ns,1)) tempo(t_spi{iii,trial,co}(ns,2))])
                % %
                % %
                px1=px(t_spi{iii,trial,co}(ns,1):t_spi{iii,trial,co}(ns,2));
                py1=py(t_spi{iii,trial,co}(ns,1):t_spi{iii,trial,co}(ns,2));
                filt_cs1=filt_cs(t_spi{iii,trial,co}(ns,1):t_spi{iii,trial,co}(ns,2));
                env1=env(t_spi{iii,trial,co}(ns,1):t_spi{iii,trial,co}(ns,2));
                if length(freq)<(t_spi{iii,trial,co}(ns,2))

                    frq1=freq(t_spi{iii,trial,co}(ns,1):t_spi{iii,trial,co}(ns,2)-1);
                else
                    frq1=freq(t_spi{iii,trial,co}(ns,1):t_spi{iii,trial,co}(ns,2));

                end
                tempo1=tempo(t_spi{iii,trial,co}(ns,1):t_spi{iii,trial,co}(ns,2));
                pxx=rad2deg(px1);

                numb_bins=4;
                [idx,cl]=resolt_bins(numb_bins,pxx);


                for q=1:size(idx,1)


                    %                                         figure(2)
                    %                                         polarplot(px1(idx{q,1}),py1(idx{q,1}),'.','Color',cl(q,:))
                    %                                         hold on
                    %                                         figure(3)
                    %                                         plot(tempo1(idx{q,1}),filt_cs1(idx{q,1}),'.','Color',cl(q,:))
                    %                                         hold on

                    dum=filt_cs1(idx{q,1});
                    bins=1:samplerate2/2:length(dum)-1;
                    for ii=1:length(bins)-1
                        pp(1,ii)=peak2peak(dum(bins(ii):bins(ii+1)));
                    end

                    ad.metrics{iii-1,co,trial}(ns,q,:) = [mean(frq1(idx{q,1}(1:end-1)))];

                    gg{iii-1,co}(alli,q)=[mean(frq1(idx{q,1}(1:end-1)))];


                    clear dum bins pp
                end

                
                if length(px1)~=length(frq1)
                frr1=[frq1 NaN];
                else
                    frr1=frq1;
                end

                ad.spirals{iii-1,co}{cn,1}=[px1;py1;frr1];
                cn=cn+1;
                clear px1 py1 pxx idx frq1 frr1
            end

            close all
            frr=[freq NaN];
            ad.sig{iii-1,co,trial}=[px;py;tempo;frr]; clear frr

          clear signal1 s1

        end
    end
end

cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/')
clearvars -except ad numb_bins pca1 norm gg
samplerate2=100;
if pca1==1
    name='pca';
else
    name='';
end

if norm==1
    name1='zsc';
else
    name1='';
end

%% HEAT MAP


peaki=[4 4 3.5];

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

for ss=1:size(ad.metrics,1)
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

        %         limits(1,:)=[min(env_all{ss,cc}) max(env_all{ss,cc})-(0.5*max(env_all{ss,cc}))];
        %         limits(2,:)=[min(horzcat(env_all{ss,:})) (max(horzcat(env_all{ss,:})))/mean(horzcat(env_all{ss,:}))];
        %
        %          limits(2,:)=[min(abs((peaki(ss)-horzcat(env_all{ss,:})))) max(abs(peaki(ss)-horzcat(env_all{ss,:})))];
        limits(1,:)=[ 0 2];
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
            senv=smoothdata(abs(peaki(ss)-pdat(3,:)),'movmean',samplerate2/2);
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





%% ------------------

% % %%%%%% PER QUADRANT
%%% clearvars -except gg
% % %%% gg (patient, condition)(spiral,quadrants)
% %
% %
% %
% % % filename=strcat('frqspiral_',num2str(name1),'_',num2str(name),num2str(numb_bins),'bins.mat');
% % % % % save(filename)
% % r=0;
% % for iii=1:3
% %     for co=1:3
% %
% %         mgg{iii,co}= mean(gg{iii,co});
% %         peaki=[4 4 3.5];
% %         ns1=1;
% %         n1=abs(mgg{iii,co}-peaki(iii));
% %         if max(n1)<4
% %             for ii=1:length(n1)
% %                 if n1(ii)<0.5
% %                     p1(ns1,ii)=1;
% %                 elseif n1(ii)>=0.5 && n1(ii)<1
% %                     p1(ns1,ii)=2;
% %                 elseif n1(ii)>=1 && n1(ii)<2
% %                     p1(ns1,ii)=3;
% %                 else
% %                     p1(ns1,ii)=4;
% %                 end
% %             end
% %         else
% %             error('outlier')
% %         end
% %
% %         code_mgg{iii,co}=p1;
% %
% %         for ns=1:size(gg{iii,co},1)
% %             y= gg{iii,co}(ns,:);
% %             n=abs(y-peaki(iii));
% %             if max(n)<4
% %                 for ii=1:length(n)
% %                     if n(ii)<0.5
% %                         p(ns,ii)=1;
% %                     elseif n(ii)>=0.5 && n(ii)<1
% %                         p(ns,ii)=2;
% %                     elseif n(ii)>=1 && n(ii)<2
% %                         p(ns,ii)=3;
% %                     else
% %                         p(ns,ii)=4;
% %                     end
% %                 end
% %             else
% %                 error('outlier')
% %             end
% %
% %         end
% %
% %         for ip=1:size(p,2)
% %             total{iii,co}(ip,:)=[numel(find(p(:,ip)==1))./size(p,1) numel(find(p(:,ip)==2))./size(p,1) numel(find(p(:,ip)==3))./size(p,1) numel(find(p(:,ip)==4))./size(p,1)];
% %
% %             f1=figure(co)
% %             subplot(1,4,ip)
% %             pie(total{iii,co}(ip,:))
% %             box('off')
% %
% %         end
% %
% %         f1.Units = 'centimeters';
% %         f1.OuterPosition= [10, 10, 35, 8];
% %         set(f1,'color','w')
% %         legend('{<0.5HZ}','{0.5-1Hz}','{1-2Hz}','{>2Hz}')
% %         legend('boxoff')
% %
% %         clear p y p1 ns1 n
% %     end
% % end
% %
