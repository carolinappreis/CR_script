function  [out,clust,freq_bl,amp_bl]=mod2(out,co,iii,s,freq_bl,amp_bl,start,ending,yy,clust,spiral,h_up)

if spiral==0
    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_out.mat','bs_end','bs_begin')
else
    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_spiral.mat','bs_end','bs_begin')
end

envelope=s.env{iii,co};
freqi=s.freq{iii,co};
phase=s.phase{iii,co};
zenv=s.zenv{iii,co};
env_acc=s.env_acc{iii,co};
z_sig=s.z{iii,co};


if co==1
    
    out.start_c{iii,co}=bs_begin(iii,clust.idx{iii,co});
    out.ending_c{iii,co}=bs_end(iii,clust.idx{iii,co});
    out.sns{iii,co}=start{iii,co};
    out.ens{iii,co}=ending{iii,co};
    
%     
%     %%%%% amp & freq surr
%     load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/trials_for_newcode.mat')
%     %%% trials correspond to avg number of stim runs delivered in the
%     %%% phase_search condition (1st col all, 2n col low med split and 3rd
%     %%% col high med split
%     rep = trials{spiral+1,1}(iii,1);
%     for ax=1:3
%         a1=amp_bl(iii,ax,(clust.idx{iii,co}));
%         b1=freq_bl(iii,ax,(clust.idx{iii,co}));
%         duma1 = a1(randi(length(a1), 1e6, rep));
%         dumb1 = b1(randi(length(b1), 1e6, rep));
%         out.mod_amp{iii,co}(ax,:)=nanmedian(duma1,2); clear duma1 a1
%         out.mod_freq{iii,co}(ax,:)=nanmedian(dumb1,2); clear dumb1 b1
%     end
%     
    
    %%%%% tremor properties hist- frequency and amplitude
    
    [match_ax]=link_ax(spiral);
    
    
    if ~isnan(clust.win(iii,1))
        idx_ns{iii,1}=h_up{iii,co};
        out.ns{iii,1}=z_sig(match_ax(co,iii,1),h_up{iii,co});
        out.env_ns{iii,1}=envelope(match_ax(co,iii,1),h_up{iii,co});
        out.zenv_ns{iii,1}=zenv(match_ax(co,iii,1),h_up{iii,co});

        
    else
        st= out.start_c{iii,co};
        ed= out.ending_c{iii,co};
        for j=1:length(st)
            idx_cseg(:,j)=st(j):ed(j);
            val_cluster(:,j)=z_sig(match_ax(co,iii,1),st(j):ed(j));
            val_cluster1(:,j)=envelope(match_ax(co,iii,1),st(j):ed(j));
            val_cluster2(:,j)=zenv(match_ax(co,iii,1),st(j):ed(j));

        end
        indi=idx_cseg(:);
        val=val_cluster(:);
        val1=val_cluster1(:);
        val2=val_cluster2(:);
        
        [indi_c, ia] = unique(indi);
        idx_ns{iii,1}=indi_c;
        out.ns{iii,1}=val(ia,1);
        out.env_ns{iii,1}=val1(ia,1);
        out.zenv_ns{iii,1}=val2(ia,1);
    end
    
    ns_amp=envelope(match_ax(co,iii,1),idx_ns{iii,co});
    ns_freq=freqi(match_ax(co,iii,1),idx_ns{iii,co});
    [min(ns_freq) max(ns_freq)];
    bins=2:0.5:9;
    bins=2:0.25:9;
    for i=1:length(bins)
        if i+1<length(bins)
            dum=find(ns_freq>bins(i) & ns_freq<=bins(i+1));
            if ~isempty (dum)
                m_n_amp(1,i)=nanmedian(ns_amp(dum));
            else
                m_n_amp(1,i)=NaN;
            end
            clear dum
        end
    end
    
    out.amp_n_bins(iii,:)=m_n_amp./max(m_n_amp) %% /nanmedian(ns_amp);
    out.bins=bins;
    
else
    
    out.start_c{iii,co}=start{iii,co}(clust.idx{iii,co});
    out.ending_c{iii,co}=ending{iii,co}(clust.idx{iii,co});
    out.yy{iii,co}=yy{iii,co}(clust.idx{iii,co});
    
    st=out.start_c{iii,co};
    et=out.ending_c{iii,co};
    yyt=out.yy{iii,co};
    
    
    phase_difs=[wrapToPi(phase(1,:)-phase(2,:));wrapToPi(phase(1,:)-phase(3,:));wrapToPi(phase(2,:)-phase(3,:))];
    
    for i=1:length(st)
        if (~isnan(st(i)))
            cont=[sum(envelope(1,st(i):et(i))) sum(envelope(2,st(i):et(i))) sum(envelope(3,st(i):et(i)))];
            ma(i)=find(cont==max(cont));
        end
    end
    
    out.axe_sum(iii,:)=hist(ma,1:3);
    
        for ax=1:3
    
            tremor_or2=NaN(1,length(st));
            tremor_or3=NaN(1,length(st));
            tremor_f2=NaN(20,5001);
            tremor_f22=NaN(20,5001);
            env_5=NaN(1,length(st));
            e_env1=NaN(length(st),6000);
    
            for i=1:length(st)
                if (~isnan(st(i)))
                    out.psi{iii,1}(ax,i)=circ_r((phase_difs(ax,st(i):et(i)))');
                    tremor_or2(1,i)=(mean(envelope(ax,et(i)-1000:et(i)))-mean(envelope(ax,st(i)-1000:st(i))))/mean(envelope(ax,st(i)-1000:st(i)));
                    tremor_or3(1,i)=mean(envelope(ax,st(i)-1000:st(i)));
    
                    env_5(1,i)=mean(((envelope(ax,st(i):st(i)+5000-1))./(max(envelope(ax,st(i):st(i)+5000-1)))));

    
                    tremor_f2(i,1:(et(i)-st(i)+1))=unwrap(phase(ax,st(i):et(i)));
                    tremor_f22(i,1:(et(i)-st(i)+1))=(phase(ax,st(i))+(0:1:(et(i)-st(i)))*2*pi/(1000./mean(freqi(ax,st(i)-1000:st(i)))));
                    tremor_k(i,1)= (tremor_f2(i,(et(i)-st(i)+1))-tremor_f22(i,(et(i)-st(i)+1)))/(2*pi*0.001*(et(i)-st(i))); %mean(frequency(et(i)-1000:et(i)));%
    
                    e_env1(i,:)=env_acc(ax,st(i)-1000:st(i)+5000-1);
    
                else
                    env_5(1,i)=NaN;
                    e_env1(i,:)=NaN(1:5000);
                    tremor_or2(1,i)=NaN;
                    tremor_or3(1,i)=NaN;
    
                    tremor_f22(i,1:5001)=NaN;
                    tremor_f2(i,1:5001)=NaN;
                    tremor_k(i,1)=NaN;
                end
            end
    
            out.all_amp{iii,co}=tremor_or3;
            out.env5_time{iii,1}(ax,:)=env_5;
    
    
            env5_p=NaN(20,12);
            tt=NaN(20,12);
            ttf=NaN(20,12);
            time_zenv=NaN(20,12,6000);
    
            for i=1:12
                env5_p(1:sum(yyt==i),i)=env_5(1,find(yyt==i)); env5_p(env5_p==0)=NaN;
                tt(1:sum(yyt==i),i)=tremor_or2(1,find(yyt==i)); tt(tt==0)=NaN;
                ttf(1:sum(yyt==i),i)=tremor_k(find(yyt==i)); ttf(ttf==0)=NaN;
                time_zenv(1:sum(yyt==i),i,:)=e_env1((find(yyt==i)),:);
            end
            out.seg_zenv{iii,1}{ax,1}=time_zenv;
            out.norm_env{iii,1}{ax,1}=env5_p;
            out.mod_amp{iii,co}{ax,1}=tt;
            out.mod_freq{iii,co}{ax,1}=ttf;
    
            clear tremor_or2 tremor_or3 tt tt1 tt2 ttf tremor_f2 tremor_f22 env5_p env_5
        end
end
end

