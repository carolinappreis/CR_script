function  [out]=mod_c(clust,out,co,iii,s,cohort,h_up,spiral)
% rng(gen)
m_ax=1;
z_sig=s.z{iii,co};
envelope=s.env{iii,co};
freqi=s.freq{iii,co};
phase=s.phase{iii,co};
% zenv=s.zenv{iii,co}.*(10*9.81/0.5);
zenv=s.zenv{iii,co};
env_acc=s.env_acc{iii,co};


if co==1
    
        if spiral==0
            load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_out.mat','amp_bbl','change_bl','freq_bl')
        else
            load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_spiral.mat','amp_bbl','change_bl','freq_bl')
        end

    
    
    
    
    %      t=1:length(env_acc(m_ax,:));
    %         plot(t,env_acc(m_ax,:))
    
%     for m=1:length(out.sns{iii,1})
%         out.all_amp{iii,co}(1,m)=mean(envelope(m_ax,out.sns{iii,1}(1,m)-1000:out.sns{iii,1}(1,m)));
% %                 nsa{iii,1}(1,m)=mean(env_acc(m_ax,out.sns{iii,co}(1,m):out.ens{iii,1}(1,m)));
%         %         hold on
%         %        plot(t(out.sns{iii,co}(1,m):out.ens{iii,1}(1,m)),env_acc(m_ax,out.sns{iii,co}(1,m):out.ens{iii,1}(1,m)));
%     end
    
%        out.nsa=mean(nsa{iii,1});
%     
%        out.std_nsa=std(nsa{iii,1});
    
    
    
%     tremor properties hist- frequency and amplitude
% % % % %         if ~isnan(clust.win(iii,1))
% % % % %             idx_ns{iii,1}=h_up{iii,co};
% % % % %             out.ns{iii,1}=z_sig(m_ax,h_up{iii,co});
% % % % %             out.env_ns{iii,1}=envelope(m_ax,h_up{iii,co});
% % % % %     
% % % % %         else
% % % % %             st= out.start_c{iii,co};
% % % % %             ed= out.ending_c{iii,co};
% % % % %             for j=1:length(st)
% % % % %                 idx_cseg(:,j)=st(j):ed(j);
% % % % %                 val_cluster(:,j)=z_sig(m_ax,st(j):ed(j));
% % % % %                 val_cluster1(:,j)=envelope(m_ax,st(j):ed(j));
% % % % %             end
% % % % %             indi=idx_cseg(:);
% % % % %             val=val_cluster(:);
% % % % %             val1=val_cluster1(:);
% % % % %     
% % % % %             [indi_c, ia] = unique(indi);
% % % % %             idx_ns{iii,1}=indi_c;
% % % % %             out.ns{iii,1}=val(ia,1);
% % % % %             out.env_ns{iii,1}=val1(ia,1);
% % % % %         end
% % % % %     
    
    
% % % % %         ns_amp=envelope(m_ax,idx_ns{iii,co});
% % % % %         ns_freq=freqi(m_ax,idx_ns{iii,co});
% % % % %         [min(ns_freq) max(ns_freq)];
% % % % %         bins=2:0.5:9;
% % % % %         bins=2:0.25:9;
% % % % %         for i=1:length(bins)
% % % % %             if i+1<length(bins)
% % % % %                 dum=find(ns_freq>bins(i) & ns_freq<=bins(i+1));
% % % % %                 if ~isempty (dum)
% % % % %                     m_n_amp(1,i)=nanmedian(ns_amp(dum));
% % % % %                 else
% % % % %                     m_n_amp(1,i)=NaN;
% % % % %                 end
% % % % %                 clear dum
% % % % %             end
% % % % %         end
% % % % %     
% % % % %         out.amp_n_bins(iii,:)=m_n_amp./max(m_n_amp) %% /nanmedian(ns_amp);
% % % % %         out.bins=bins;
    
    
% % % % % %     % % % % %     %%% baseline 50000 - taken from aux_master.m
% % % % % %     %%%amp surr
    rep = out.trials(iii,1);
%     rep=10;
    for ax=1:3
        bb(ax,:)=change_bl(iii,ax,(clust.idx{iii,co}));
        ba(ax,:)=amp_bbl(iii,ax,(clust.idx{iii,co}));
        dum = bb(randi(length(bb), 1e6, rep));
        out.change_c{iii,co}(ax,:)=nanmedian(dum,2); clear dum
        for i = 1:12
            dum = bb(randi(length(bb), 1e6, rep));
            out.ns_arc{iii,ax}(:,i)=nanmedian(dum,2); clear dum
        end
        
        clear dum dum2 baseline3
    end
% % % % % %     
% % % % % %     %%%freq surr
    rep = out.trials(iii,1);
    %rep=10;
    for ax=1:3
        bb(ax,:)=freq_bl(iii,ax,(clust.idx{iii,co}));
        ba(ax,:)=amp_bbl(iii,ax,(clust.idx{iii,co}));
        dum = bb(randi(length(bb), 1e6, rep));
        out.fchange{iii,co}(ax,:)=nanmedian(dum,2); clear dum
        for i = 1:12
            dum = bb(randi(length(bb), 1e6, rep));
            out.ns_prc{iii,ax}(:,i)=nanmedian(dum,2); clear dum
        end
        
        clear dum dum2 baseline3
    end
    
    
    
    %%%baseline median split
    
    ns_a=NaN(2,1e6);
    
    seg=bb(m_ax,:);
    ba1=ba(m_ax,:);
    
    a_1=seg(1,find(ba1<=median(ba1)));
    a_2=seg(1,find(ba1>median(ba1)));
    
    
    %     rep=5;
    dum1 = a_1(randi(length(a_1), 1e6, out.trials(iii,2)));
    dum2 = a_2(randi(length(a_2), 1e6, out.trials(iii,3)));
    ns_a(1,:) = (nanmedian(dum1,2))';
    ns_a(2,:) = (nanmedian(dum2,2))';
    clear dum1 dum2
    
    
    
        out.ns_ms{iii,1}=ns_a;
    
    % % % % %     %%%% median power ns --- because it used to be compared with PLS no
    % % % % %     %%%% need for clustering
    % % % % %     % % %     st=out.sns{iii,co}; et=out.ens{iii,co};
    % % % % %     % % %     sig=[2 4 5 7 10];
    % % % % %     % % %     if ismember(iii,sig)
    % % % % %     % % %         z_seg1=NaN(length(st),3);
    % % % % %     % % %
    % % % % %     % % %         for i=1:length(st)
    % % % % %     % % %             if (~isnan(st(i))) && (st(i)-5000)>0
    % % % % %     % % %                 ghost=[st(i)+5000:et(i)-5000];
    % % % % %     % % %                 dum=st(i)+5000+round(length(ghost)/2);
    % % % % %     % % %                 z_seg1(i,:)=[mean(zenv(m_ax,st(i):st(i)+5000)) mean(zenv(m_ax,dum-5000:dum)) mean(zenv(m_ax,et(i)-5000:et(i)))];
    % % % % %     % % %                 clear ghost dum dum1
    % % % % %     % % %             end
    % % % % %     % % %         end
    % % % % %     % % %         out.ns3(iii,:)=median(z_seg1);
    % % % % %     % % %         out.ns_hu{iii,1}=z_sig(m_ax,h_up{iii,co});
    % % % % %     % % %     else
    % % % % %     % % %         out.ns3(iii,:)=NaN(1,3);
    % % % % %     % % %         out.ns_hu{iii,1}=[];
    % % % % %     % % %     end
    
else
    
    if co==2
        
        st=out.start_c{iii,co};
        et=out.ending_c{iii,co};
        yyt=out.yy{iii,co};
        
        for ax=1:3
            tremor_or2=NaN(1,length(st));
            tremor_or3=NaN(1,length(st));
            tremor_or4=NaN(1,length(st));
            tremor_f2=NaN(20,5001);
            tremor_f22=NaN(20,5001);
            z_seg1=NaN(length(st),5000);
            e_env1=NaN(length(st),6000);
            env_5=NaN(1,length(st));
            
            
            for i=1:length(st)
                if (~isnan(st(i)))
                    tremor_or2(1,i)=(mean(envelope(ax,et(i)-1000:et(i)))-mean(envelope(ax,st(i)-1000:st(i))))/mean(envelope(ax,st(i)-1000:st(i)));
                    tremor_or3(1,i)=mean(envelope(ax,st(i)-1000:st(i)));
                    tremor_or4(1,i)=mean(env_acc(ax,et(i)-1000:et(i)));
                    
                    env_5(1,i)=mean(((envelope(ax,st(i):st(i)+5000-1))./(max(envelope(ax,st(i):st(i)+5000-1)))));
                    
                    tremor_f2(i,1:(et(i)-st(i)+1))=unwrap(phase(ax,st(i):et(i)));
                    tremor_f22(i,1:(et(i)-st(i)+1))=(phase(ax,st(i))+(0:1:(et(i)-st(i)))*2*pi/(1000./mean(freqi(ax,st(i)-1000:st(i)))));
                    tremor_k(i,1)= (tremor_f2(i,(et(i)-st(i)+1))-tremor_f22(i,(et(i)-st(i)+1)))/(2*pi*0.001*(et(i)-st(i))); %mean(frequency(et(i)-1000:et(i)));%
                    
                    z_seg1(i,1:5000)=z_sig(m_ax,st(i):st(i)+5000-1);
%                     e_env1(i,1:5000)=envelope(m_ax,st(i):st(i)+5000-1);
                    e_env1(i,1:6000)=zenv(m_ax,st(i):st(i)+6000-1);
                    
                    
                    
                else
                    env_5(1,i)=NaN;
                    tremor_or2(1,i)=NaN;
                    tremor_or3(1,i)=NaN;
                    tremor_or4(1,i)=NaN;
                    
                    tremor_f22(i,1:5001)=NaN;
                    tremor_f2(i,1:5001)=NaN;
                    tremor_k(i,1)=NaN;
                    
                    
                    z_seg1(i,1:5000)=NaN;
                    e_env1(i,1:6000)=NaN;
                    yyt(i)= NaN;
                end
            end
            
            out.all_amp{iii,co}=tremor_or3;
            out.env5_time{iii,1}(ax,:)=env_5;
            
            amp_1=NaN(2,round(size(tremor_or3,2)./2));
            ch_a1=NaN(2,round(size(tremor_or3,2)./2));
            pha_idx=NaN(2,round(size(tremor_or3,2)./2));
            m=1;
            n=1;
            for i=1:length(st)
                if tremor_or3(i)<=nanmedian(tremor_or3)
                    amp_1(1,n)= tremor_or3(i);
                    ch_a1(1,n)= tremor_or2(i);
                    pha_idx(1,n)=yyt(i);
                    n=n+1;
                    
                else
                    amp_1(2,m)= tremor_or3(i);
                    ch_a1(2,m)= tremor_or2(i);
                    pha_idx(2,m)=yyt(i);
                    m=m+1;
                    
                end
            end
            
            env5_p=NaN(20,12);
            
            tt=NaN(20,12);
            ee=NaN(20,12);
            
            tt1=NaN(20,12);
            amp1=NaN(10,12);
            tt2=NaN(20,12);
            amp2=NaN(10,12);
            ttf=NaN(20,12);
            
            for i=1:12
                
                env5_p(1:sum(yyt==i),i)=env_5(1,find(yyt==i));
                tt(1:sum(yyt==i),i)=tremor_or2(1,find(yyt==i));
                ee(1:sum(yyt==i),i)=tremor_or4(1,find(yyt==i));
                tt(tt==0)=NaN;
                ttf(1:sum(yyt==i),i)=tremor_k(find(yyt==i));
                
                out.z_seg{iii,i}= z_seg1(find(yyt==i),:);
                out.zenv_seg{iii,i}= e_env1(find(yyt==i),:);
                
                
                tt1(1:sum(pha_idx(1,:)==i),i)=ch_a1(1,find(pha_idx(1,:)==i));
                amp1(1:sum(pha_idx(1,:)==i),i)=amp_1(1,find(pha_idx(1,:)==i));
                tt2(1:sum(pha_idx(2,:)==i),i)=ch_a1(2,find(pha_idx(2,:)==i));
                amp2(1:sum(pha_idx(2,:)==i),i)=amp_1(2,find(pha_idx(2,:)==i));
            end
            
            out.env5_p{iii,co}{ax,1}=env5_p;
            
            out.arc1{iii,co}{ax,1}=tt1;
            out.arc2{iii,co}{ax,1}=tt2;
            out.change_c{iii,co}{ax,1}=tt;
            out.end_env{iii,co}{ax,1}=ee;
            out.fchange{iii,co}{ax,1}=ttf;
            out.amp1{iii,co}{ax,1}=amp1;
            out.amp2{iii,co}{ax,1}=amp2;
            
            clear tremor_or2 tremor_or3 tt tt1 tt2 ttf tremor_f2 tremor_f22 env5_p env_5
        end
        
    else
        
        sig=[2 4 5 7 10];
        if ismember(iii,sig)
            load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/hu_hd_pls.mat','st')
            st=st{find(sig==iii),1};
            et=out.ending_c{iii,co};
            z_seg1=NaN(length(st),3);
            
            for i=1:length(st)
                if (~isnan(st(i))) && (st(i)-5000)>0
                    ghost=[st(i)+5000:et(i)-5000];
                    dum=st(i)+5000+round(length(ghost)/2);
                    z_seg1(i,:)=[mean(zenv(m_ax,st(i):st(i)+5000)) mean(zenv(m_ax,dum-5000:dum)) mean(zenv(m_ax,et(i)-5000:et(i)))];
                    clear ghost dum dum1
                end
            end
            out.pls3(iii,:)=median(z_seg1);
            out.pls_hu{iii,1}=z_sig(m_ax,h_up{iii,co});
        else
            out.pls3(iii,:)=NaN(1,3);
            out.pls_hu{iii,1}=[];
        end
    end
end

