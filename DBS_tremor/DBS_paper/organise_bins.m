clear all
ad=struct; ad.metrics=cell(1,1); ad.spirals=cell(1,1); ad.sig=cell(1,1); ad.signames={'data1';'data2';'time';'envelope'};ad.metricnames={'subj';'cond';'trials';'nspirals';'nquadrants';'mean env';'peak2peak';'variance';'freq'};

close all
cond={'NS';'HF';'C'};
cohort=[1 3 4 6];
pca1=0; %%% use the sqrt or the pca as the signal to extract the envelope from
norm=1; %%% norm=1:sqrt(x2+y2) is zscored across condition ;


for iii=2:4
    clearvars -except ad iii cohort cond pca1 norm
    for co=1:3
        cn=1;
        if (co==1 && iii==2 | iii==3)
            ntrial=[1 2];
        else
            ntrial=1;
        end
        
        for trial=1:length(ntrial)
            load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/P0',num2str(cohort(iii)),'_clean_',num2str(cond{co,1}),num2str(trial),'_SH.mat'));
            load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/zscore_acr_cond.mat'))
            
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
                    
                    ad.metrics{iii-1,co,trial}(ns,q,:) = [mean(env1(idx{q,1})); nanmean(pp) ; var(filt_cs1(idx{q,1}));mean(frq1(idx{q,1}(1:end-1))) ];
                    clear dum bins pp
                end
                
                ad.spirals{iii-1,co}{cn,1}=[px1;py1;env1];
                cn=cn+1;
                clear px1 py1 pxx idx env1
            end
            close all
            ad.sig{iii-1,co,trial}=[px;py;tempo;env];
            clear signal1 s1
        end
    end
end

cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/')
clearvars -except ad numb_bins pca1 norm
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

filename=strcat('frqspiral_',num2str(name1),'_',num2str(name),num2str(numb_bins),'bins.mat');
% % %  save(filename)