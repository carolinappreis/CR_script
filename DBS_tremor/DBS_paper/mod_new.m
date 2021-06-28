function  [out]=mod_new(clust,out,co,iii,s,spiral,start,ending)

envelope=s.env{iii,co};
freqi=s.freq{iii,co};
phase=s.phase{iii,co};
tfilt=s.filt{iii,co};


if co==1
    
    if spiral==0
        load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_out.mat','bs_begin','bs_end')
    else
        load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_spiral.mat','bs_begin','bs_end')
    end
    
    
    fff=NaN(3,5e4); aaa=NaN(3,5e4);clt = NaN(3,5e4,5001);
    
    parfor j = 1:length(bs_end)
        begin3=bs_begin(iii,j);
        end3=bs_end(iii,j);
        for ax = 1:3
            [ffr,aam]=simp_new(ax,end3,begin3,freqi,phase,envelope)
            fff(ax,j)=ffr;
            aaa(ax,j)=aam;
            clt(ax,j,:)=tfilt(ax,begin3:end3);
        end
    end
    
    freq_bl(iii,:,:)=fff;
    amp_bl(iii,:,:)=aaa;
    
    
    for j = 1:5e4 % in pca, rows are observations and columns are variables
        for_pca = squeeze(clt(:,j,:)); % should be 3 vs length(segments)
        [pc, score, latent, tsquare, explained] = pca(for_pca');
        out.x_all{iii,1}(co,j,1:3) = pc(1:3, 1);
        out.pc1_exp{iii,1}(co,j,1:3) = explained;
    end
    
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/trials_for_newcode.mat')
    %%% trials correspond to avg number of stim runs delivered in the
    %%% phase_search condition (1st col all, 2n col low med split and 3rd
    %%% col high med split
    
    
    %%%%%%%%% amp & freq surr
    rep = trials{spiral+1,1}(iii,1);
    for ax=1:3
        a1=amp_bl(iii,ax,(clust.idx{iii,co}));
        b1=freq_bl(iii,ax,(clust.idx{iii,co}));
        duma1 = a1(randi(length(a1), 1e6, rep));
        dumb1 = b1(randi(length(b1), 1e6, rep));
        out.amp_surr{iii,co}(ax,:)=nanmedian(duma1,2); clear duma1 a1
        out.freq_surr{iii,co}(ax,:)=nanmedian(dumb1,2); clear dumb1 b1
    end
    
else
    x=[];
    for j = 1:length(start{iii,co})
        if (~isnan(start{iii,co}(j)))
            for i=1:3
                x(i,:)= tfilt(i,(start{iii,co}(j):ending{iii,co}(j)));
            end
            [pc, score, latent, tsquare, explained] = pca(x');
            out.x_all{iii,1}(co,j,1:3) = pc(1:3, 1);
            out.pc1_exp{iii,1}(co,j,1:3) = explained;
        end
    end
    
    
end


end

