

clear all; close all
cohort = [1 3 4 6];
pca1=1;


for iii=4
    signal=[]; unf_signal=[];
    for cond=2
        s=struct;
        if cond==1
            load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_NS.mat'))
        elseif cond==2
            load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_RS.mat'))
        elseif cond==3
            load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_spiral/P0',num2str(cohort(iii)),'_PLS_S.mat'))
        else
            load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_hf/P0',num2str(cohort(iii)),'_HFS.mat'))
        end
        
        [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
        data=d.data_raw;  tremor_ds=d.data_ds;
        addon=92; addon_end=35;
        
        if cond==1
            NS_BE_S
            start=hu{iii,:};
            ending=hd{iii,:};
            
        elseif cond==2
            co=2;start=[];ending=[];yy=[];h_up=[];spiral=1;
            [peak_ax, start, ending, yy, h_up, ts, te]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
            start=ts{iii,co};
            ending=te{iii,co};
            
        elseif cond==3
            new = find(data(2,:) > 4);
            difp = find((diff(new)) > 100000); % are you trying to threshold at 9.6 seconds?
            ep_1 = [new(difp) new(end)];
            sp_1 = [new(1) new(difp+1)];
            
            dum=find((ep_1-sp_1)<(60000./samplerate)*samplerateold);
            if ~isempty(dum)
                sp_1(1,dum)=NaN;
                ep_1(1,dum)=NaN;
            end
            sp=sp_1(~isnan(sp_1));
            ep=ep_1(~isnan(ep_1));
            
            if iii==4
                
                ep(2)=ep(1);
                ep(1)=953991;
                sp(2)=1022995;
            end
            %                         time=1:length(data(1,:));
            %                         plot(time,data(4,:))
            %                         hold on
            %                         plot(time(sp),data(4,sp),'r.')
            %                         plot(time(ep),data(4,ep),'k.')
            
            start= floor((sp./samplerateold)*samplerate)+addon;
            ending = floor((ep./samplerateold)*samplerate)+addon+addon_end;
            
        else
            
            HFS_BE_S
            start=hu{iii,:};
            ending=hd{iii,:};
            
        end
        
        for i=1:3
        y(i,:)=zscore(d.data_ds(i,:));
        end
        
        y(4,:)=sqrt((y(1,:).^2+y(2,:).^2+y(3,:).^2));
        

        handup = [];
        for ix = 1:length(start)
            handup = [handup start(ix):ending(ix)]; %#ok<*AGROW>
        end
        handup = sort(handup,'ascend');
        
        samplerate=1000;
        [a,b]=  butter(3,[1/(0.5*samplerate)], 'high'); %15
        
        for i=1:4
            xyzc(i,:)=(filtfilt(a,b,y(i,handup)));
%             figure(1)
%             plot(xyzc(i,:))
%             hold on
%             [peak,Pxx,F,cut]=dindfpeak(xyzc(i,:),samplerate);

        end
        
        [peak,Pxx,F,cut]=dindfpeak(xyzc(4,:),samplerate);

        comb_signal{iii-1,cond,:}=xyzc(4,:);

        Pcurve(iii-1,cond,:)=Pxx./(sum(Pxx(100:500)));
        Ppeak(iii-1,cond,:)=max(squeeze(Pcurve(iii-1,cond,cut(1):cut(2))));
        fpeak(iii-1,cond)=peak;

         
        if pca1==0
            next_sig=(filtfilt(a,b,xyzc(4,:)));
        else
            [a,b]=  butter(2, [2/(0.5*samplerate) 8/(0.5*samplerate)], 'bandpass'); %15

            for ip=1:3
                to_pc(ip,:)=(filtfilt(a,b,xyzc(ip,:)));
            end
            
            dx = [to_pc(1,:); to_pc(2,:);to_pc(3,:)];
            [pc, score, latent, tsquare, explained] = pca(dx');
            
%             plot(pc(1,1)*score(:,1))
%             hold on
%             plot(to_pc(1,:),'.')
            
            next_sig= to_pc(1,:).*pc(1,1)+ to_pc(2,:).*pc(2,1)+ to_pc(3,:).*pc(3,1);
            pca_signal{iii-1,cond,:}=next_sig;

            
%        close all
       [peak,Pxx,F,cut]=dindfpeak(next_sig,samplerate);
       fpeak(iii-1,cond)=peak;
       Pcurve(iii-1,cond,:)=Pxx./(sum(Pxx(100:500)));
       Ppeak(iii-1,cond,:)=max(squeeze(Pcurve(iii-1,cond,cut(1):cut(2))));
            
        end
        

       
       [a,b]=  butter(2, [(peak-2)/(0.5*samplerate) (peak+2)/(0.5*samplerate)], 'bandpass'); %15

      sig_fh=(filtfilt(a,b,next_sig));
       
 
       env_h=abs(hilbert(sig_fh));
       
        %%% clean outliers
        idp=find(sig_fh>(mean(sig_fh)+5*(std(sig_fh))));
        idm=find(sig_fh<(mean(sig_fh)-5*(std(sig_fh))));
        indexexceed=[idp idm];
        indexexceed = sort(indexexceed,'ascend');

        
        diffindex=diff(indexexceed);
         pnts=find(diffindex>1);
%         pnts=[];
        if ~isempty(pnts)
            bb(1,:)=indexexceed(pnts+1);
            ee(1,:)=indexexceed(pnts);
            begin2=[indexexceed(1) bb];
            ending2=[ee indexexceed(end)];
            
            ind_b=[];
            for i=1:(length(begin2))-1
                if (ending2(i)-begin2(i))>1
                    ind_b=[ind_b i];
                end
            end
            
            begin3=begin2(ind_b);
            ending3=ending2(ind_b);
            
            sig_hdum=sig_fh;
            for ix = 1:length(begin3)
                sig_hdum(begin3(ix):ending3(ix))=NaN;
            end
            sig_hc=sig_hdum(~isnan(sig_hdum));
            env_hc=env_h((~isnan(sig_hdum)));
        else
            sig_hc=sig_fh;
            env_hc=env_h;
        end

        
% %         figure
% %         subplot(2,1,1)
% %         plot(sig_fh)
% %         hold on
% %         plot(env_h)
% %         yline((mean(sig_fh)+5*(std(sig_fh))))
% %         yline((mean(sig_fh)-5*(std(sig_fh))))
% %         subplot(2,1,2)
% %         plot(sig_hc)
% %         hold on
% %         plot(env_hc)
% %         hold on
% %         yline((mean(sig_fh)+5*(std(sig_fh))))
% %         yline((mean(sig_fh)-5*(std(sig_fh))))
% %         
        
%         close all

         avg_t_sev(iii-1,cond,:)=[mean(env_h) std(env_h)];
       
        unf_signal=[unf_signal sig_hc];
        signal= [signal env_h];
        clearvars -except cohort cond iii fpeak Pcurve Ppeak pca1 signal unf_signal ufsignal_all avg_t_sev comb_sig pca_sig pca1 env_all
    end
    env_all{iii-1,1}=signal;
    ufsignal_all{iii-1,1}=unf_signal;
end

for i=1:3
    subplot(1,3,i)
    data=env_all{i,1};
    pct(i,:)=[prctile(data,25) prctile(data,50)];
    plot(data)
    hold on
    yline(pct(i,1),'r','LineWidth',2)
    yline(pct(i,2),'k','LineWidth',2)
    clear data fdata envdata
end

clearvars -except pct ufsignal_all Pcurve Ppeak avg_t_sev comb_sig pca_sig  sig_all env_all pca1 fpeak

if pca1==1
    name='pca';
else
    name='sqrt';
end

filename=strcat('ACCspiral_',num2str(name),'.mat');
save(filename)






for i=1:3
    bar_tremor(i,:)=avg_t_sev(i,:,1)./(pct(i,2).*(10*9.81/0.5));
    
    signal=ufsignal_all{i,1}; samplerate =1000;
    [Pxx,F] = pwelch(signal, samplerate, [], round(samplerate*2), samplerate);
    cur(i,:)=Pxx./(sum(Pxx(100:500)));
    peak_concat(i,:)=max(squeeze(cur(i,7:21)));
    
end

figure
bar(bar_tremor)

figure
for ii=1:3
    for cc=1:4
        subplot(1,3,ii)
        plot(F,(squeeze(Pcurve(ii,cc,:))./ peak_concat(ii)))
        hold on
    end
    point(ii,:)=Ppeak(ii,:)./peak_concat(ii);
    for i=2:4
        change(ii,i-1,:)=round(((point(ii,i,:)-point(ii,1,:))./point(ii,1,:)*100),0);
    end
    box('off')
    xlim([2 10])
end

figure
bar(-change)
box('off')
ylabel('% supression from baseline')
xlabel('patients')
