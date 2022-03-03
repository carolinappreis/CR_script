clear all
ad=struct; ad.metrics=cell(1,1); ad.spirals=cell(1,1); ad.sig=cell(1,1); ad.signames={'data1';'data2';'combsig';'time';'envelope'};ad.metricnames={'subj';'cond';'trials';'nspirals';'nquadrants';'mean';'peak2peak';'variance'};

% close all
cond={'NS';'HF';'C'};
cohort=[ 1 3 4 6];


for iii=2:4
    clearvars -except ad iii cohort cond
    for co=1:3
        cn=1;
        if (co==1 && iii==2 | iii==3)
            ntrial=[1 2];
        else
            ntrial=1;
        end
        
        for trial=1:length(ntrial)
            load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/P0',num2str(cohort(iii)),'_clean_',num2str(cond{co,1}),num2str(trial),'_SH.mat'));
            
            
            %            centre=[ mean(signal(1,:)) mean(signal(2,:))];
            centre=[535 361];
            signal1(1,:)=signal(1,:)-centre(1);
            signal1(2,:)=signal(2,:)-centre(2);
            signal1(3,:)=(sqrt((signal1(1,:).^2)+(signal1(2,:).^2)));
            
            
            %%% find peak tremor
            [peak]=dindfpeak(signal1,samplerate2);
            
            %%% filt & var
            [a,b]=  butter(2, [(peak-2)/(0.5*samplerate2) (peak+2)/(0.5*samplerate2)], 'bandpass'); %15
            filt_cs=(filtfilt(a,b,zscore(signal1(3,:))));
            env=abs(hilbert(filt_cs));
            
            %             plot(signal1(1,:))
            %             hold on
            %             plot(signal1(2,:))
            %             plot(signal1(3,:))
            run('time_split.m')
            
            
            [px,py]=cart2pol(signal1(1,:),signal1(2,:));
            for ns=1:size(t_spi{iii,trial,co}(:,1),1)
                
                %  check power/space
                %                 [mama,maid]=sort(env,'descend');
                %                 figure(1)
                %                 subplot(2,(size(t_spi{iii,trial,co}(:,1),1))/2,ns)
                %                 plot3(tempo,signal1(1,:),signal1(2,:))
                %                 hold on
                %                 plot3(tempo(maid(1:length(maid)/2)),signal1(1,(maid(1:length(maid)/2))),signal1(2,(maid(1:length(maid)/2))),'r.')
                %                 clear mama maid
                %                 [mama,maid]=sort(env,'ascend');
                %                 plot3(tempo(maid(1:length(maid)/2)),signal1(1,(maid(1:length(maid)/2))),signal1(2,(maid(1:length(maid)/2))),'k.')
                %                 xlim([tempo(t_spi{iii,trial,co}(ns,1)) tempo(t_spi{iii,trial,co}(ns,2))])
                %
                
                px1=px(t_spi{iii,trial,co}(ns,1):t_spi{iii,trial,co}(ns,2));
                py1=py(t_spi{iii,trial,co}(ns,1):t_spi{iii,trial,co}(ns,2));
                filt_cs1=filt_cs(t_spi{iii,trial,co}(ns,1):t_spi{iii,trial,co}(ns,2));
                env1=env(t_spi{iii,trial,co}(ns,1):t_spi{iii,trial,co}(ns,2));
                pxx=rad2deg(px1);
                
                idx{1,1}=find(pxx>0 & pxx< 90);
                idx{2,1}=find(pxx<0 & pxx>-90);
                idx{3,1}=find(pxx>-180 & pxx<-90);
                idx{4,1}=find(pxx>90 & pxx<180);
                
                for q=1:4
                    
                    dum=filt_cs1(idx{q,1});
                    bins=1:samplerate2/2:length(dum)-1;
                    for ii=1:length(bins)-1
                        pp(1,ii)=peak2peak(dum(bins(ii):bins(ii+1)));
                    end
                    
                    ad.metrics{iii-1,co,trial}(ns,q,:) = [mean(env1(idx{q,1})); mean(pp) ; var(filt_cs1(idx{q,1}))];
                    clear dum bins pp
                end
                
                ad.spirals{iii-1,co}{cn,1}=[px1;py1;env1];
                cn=cn+1;
                clear px1 py1 pxx idx env1 
            end
            ad.sig{iii-1,co,trial}=[px;py;tempo;env];
            clear signal1
        end
    end
end

cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/')
clearvars -except ad
samplerate2=100;
% % save ('toplot_spirals1')