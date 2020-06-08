function [out]=mod_c(clust,out,co,iii,s,h_up)
% rng(gen)
m_ax=1;
z_sig=s.z{iii,co};
envelope=s.env{iii,co};
freqi=s.freq{iii,co};
phase=s.phase{iii,co};
% zenv=s.zenv{iii,co}.*(10*9.81/0.5);
zenv=s.zenv{iii,co};

if co==1
    %%% tremor properties hist- frequency and amplitude
    
    ns_amp=envelope(m_ax,h_up{iii,co});
    ns_freq=freqi(m_ax,h_up{iii,co});
    %     [min(ns_freq) max(ns_freq)];
    dum2=[];
    bins=2:0.5:9;
    for i=1:length(bins)
        if i+1<length(bins)
            dum=find(ns_freq>bins(i) & ns_freq<=bins(i+1));
            dum2=[dum2 ; bins(i).* ones(length(dum),1)];
            if ~isempty (dum)
                m_n_amp(1,i)=nanmedian(ns_amp(dum)./nanmedian(ns_amp));
            else
                m_n_amp(1,i)=NaN;
            end
        end
    end
    
    out.amp_n_bins(iii,:)=m_n_amp;
    out.bins=bins;
    
    %%%% median power ns
    st=out.sns{iii,co}; et=out.ens{iii,co};
    sig=[2 4 5 7 10];
    if ismember(iii,sig)
        z_seg1=NaN(length(st),3);
        
        for i=1:length(st)
            if (~isnan(st(i))) && (st(i)-5000)>0
                ghost=[st(i)+5000:et(i)-5000];
                dum=st(i)+5000+round(length(ghost)/2);
                z_seg1(i,:)=[mean(zenv(m_ax,st(i):st(i)+5000)) mean(zenv(m_ax,dum-5000:dum)) mean(zenv(m_ax,et(i)-5000:et(i)))];
                clear ghost dum dum1
            end
        end
        out.ns3(iii,:)=median(z_seg1);
        out.ns_hu{iii,1}=z_sig(m_ax,h_up{iii,co});
    else
        out.ns3(iii,:)=NaN(1,3);
        out.ns_hu{iii,1}=[];
    end
    
    %%% baseline 50000 - taken from code ttrials
    rep = 10;
    for ax=1:3
        load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/NS_all','seg_bl')
        baseline3=seg_bl(iii,ax,(clust.idx{iii,co}));
        dum = baseline3(randi(length(baseline3), 1e6, rep));
        out.change_c{iii,co}(ax,:)=nanmedian(dum,2); clear dum
        clear dum dum2 baseline3
    end
    
    %             rep = 10;
    %         baseline3_temp = baseline3(ns_mat(iii,ax),:);
    %         dum = baseline3_temp(randi(length(baseline3_temp), 1e6, rep));
    %         dum2 = dum;
    %         p = nanmedian(dum2,2);
    %
    %         nostim(iii,ax,:) = p;
    %
    %         clear dum dum2 p baseline3_temp
    
    
    
    % %     rng('default') % set random seed for consistency
    % %     envelope=s.env{iii,co};
    % %
    % %
    % %     st=start{iii,co}(clust.idx{iii,co});
    % %
    % %
    % %     baseline3 = NaN(3,size(st,2));
    % %     for j = 1:size(st,2)
    % %         begin3=st(j);
    % %         end3=floor(begin3+5*samplerate);
    % %         for ax = 1:3
    % %             baseline3(ax,j) = (mean(envelope(ax,end3-1000:end3))-mean(envelope(ax, begin3-1000:begin3)))./mean(envelope(ax, begin3-1000:begin3));
    % %         end
    % %     end
    % %     rep=10;
    % %     for ax=1:size(s.raw,1)
    % %         dum = baseline3(randi(length(baseline3), 1e6, rep));
    % %         out.change{iii,co}(ax,:)=nanmedian(dum,2); clear dum
    % %     end
    
else
    
    if co==2
        
        st=out.start_c{iii,co};
        et=out.ending_c{iii,co};
        yyt=out.yy{iii,co};
        
        for ax=1:3
            tremor_or2=NaN(1,length(st));
            tremor_or3=NaN(1,length(st));
            tremor_f2=NaN(20,5001);
            tremor_f22=NaN(20,5001);
            z_seg1=NaN(length(st),5000);
            
            for i=1:length(st)
                if (~isnan(st(i)))
                    tremor_or2(1,i)=(mean(envelope(ax,et(i)-1000:et(i)))-mean(envelope(ax,st(i)-1000:st(i))))/mean(envelope(ax,st(i)-1000:st(i)));
                    tremor_or3(1,i)=mean(envelope(ax,st(i)-1000:st(i)));
                    
                    tremor_f2(i,1:(et(i)-st(i)+1))=unwrap(phase(ax,st(i):et(i)));
                    tremor_f22(i,1:(et(i)-st(i)+1))=(phase(ax,st(i))+(0:1:(et(i)-st(i)))*2*pi/(1000./mean(freqi(ax,st(i)-1000:st(i)))));
                    tremor_k(i,1)= (tremor_f2(i,(et(i)-st(i)+1))-tremor_f22(i,(et(i)-st(i)+1)))/(2*pi*0.001*(et(i)-st(i))); %mean(frequency(et(i)-1000:et(i)));%
                    
                    z_seg1(i,1:5000)=z_sig(m_ax,st(i):st(i)+5000-1);
                else
                    tremor_or2(1,i)=NaN;
                    tremor_or3(1,i)=NaN;
                    
                    tremor_f22(i,1:5001)=NaN;
                    tremor_f2(i,1:5001)=NaN;
                    tremor_k(i,1)=NaN;
                    
                    z_seg1(i,1:5000)=NaN;
                    yyt(i)= NaN;
                end
            end
            
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
            
            tt=NaN(20,12);
            
            tt1=NaN(20,12);
            amp1=NaN(10,12);
            tt2=NaN(20,12);
            amp2=NaN(10,12);
            ttf=NaN(20,12);
            
            for i=1:12
                tt(1:sum(yyt==i),i)=tremor_or2(1,find(yyt==i));
                tt(tt==0)=NaN;
                ttf(1:sum(yyt==i),i)=tremor_k(find(yyt==i));
                
                out.z_seg{iii,i}= z_seg1(find(yyt==i),:);
                
                tt1(1:sum(pha_idx(1,:)==i),i)=ch_a1(1,find(pha_idx(1,:)==i));
                amp1(1:sum(pha_idx(1,:)==i),i)=amp_1(1,find(pha_idx(1,:)==i));
                tt2(1:sum(pha_idx(2,:)==i),i)=ch_a1(2,find(pha_idx(2,:)==i));
                amp2(1:sum(pha_idx(2,:)==i),i)=amp_1(2,find(pha_idx(2,:)==i));
            end
            
            out.arc1{iii,co}{ax,1}=tt1;
            out.arc2{iii,co}{ax,1}=tt2;
            out.change_c{iii,co}{ax,1}=tt;
            out.fchange{iii,co}{ax,1}=ttf;
            
            clear tremor_or2 tremor_or3 tt tt1 tt2 ttf tremor_f2 tremor_f22
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


