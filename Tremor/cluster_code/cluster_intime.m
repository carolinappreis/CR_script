function [p1]=cluster_intime(clust,s,iii,co,out)
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/aux_out.mat','bs_end','bs_begin')


if co==1 | co==2
%     for j=1:size(out.start_c{iii,co},2)
%         y = [sum(s.env{iii,co}(1, out.start_c{iii,co}:out.ending_c{iii,co})); sum(s.env{iii,co}(2, out.start_c{iii,co}:out.ending_c{iii,co})); sum(s.env{iii,co}(3, out.start_c{iii,co}:out.ending_c{iii,co}))];
%         ma_c(j) = find(y == max(y));
%     end
%     figure(co)
%     subplot(2,5,iii)
%     hist(ma_c, 1:3)
%     box('off')
%     
    
    if co==1
        z_sig=s.z{iii,co};
        m_ax=1;
        
        comp=1:50000;
        comp(1,find(ismember(comp,clust.idx{iii,co}')))=NaN;
        st_alt=bs_begin(iii,~isnan(comp));
        et_alt=bs_end(iii,~isnan(comp));
        st_clu=bs_begin(iii,clust.idx{iii,co});
        et_clu=bs_end(iii,clust.idx{iii,co});
        
        for ii=1:2
            
            if ii==1
                st=st_clu;
                ed=et_clu;
            else
                st=st_alt;
                ed=et_alt;
            end
            
            for j=1:length(st)
                idx_cseg(:,j)=st(j):ed(j);
                val_cluster(:,j)=z_sig(m_ax,st(j):ed(j));
            end
            indi=idx_cseg(:);
            val=val_cluster(:);
            
            [indi_c, ia] = unique(indi);
            clu{ii,1}=indi_c;
            cia{ii,1}=ia;
            
            clear idx_cseg val_cluster indi_c ia
        end
        
        p=[clu{1,1} ; clu{2,1}];
        g1=clu{1,1}';
        g2=clu{2,1}';
        
        num=[];
        for i=1:length(g1)
            if ~isempty(find(g2==g1(i)))
                num=[num i];
            end
        end
        
        
        dum=g1(num);
        time=1:size(z_sig,2);
        
        
        figure(3+iii)
        subplot(3,1,1)
        plot(time,z_sig(1,:))
        hold on
        plot(time(dum),z_sig(1,dum))
        title ('overlap')
        subplot(3,1,2)
        plot(time,z_sig(1,:))
        hold on
        plot(time(g1),z_sig(1,g1))
        title ('cluster chosen')
        subplot(3,1,3)
        plot(time,z_sig(1,:))
        hold on
        plot(time(g2),z_sig(1,g2))
        title ('cluster other')
    end
end

p1=0;
end