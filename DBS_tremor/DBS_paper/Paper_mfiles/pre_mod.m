function  [freq_bl,amp_bl,out]=pre_mod(co,iii,s,spiral,start,ending,out,freq_bl,amp_bl,samplerate)

envelope=s.env{iii,co};
freqi=s.freq{iii,co};
phase=s.phase{iii,co};
tfilt=s.filt{iii,co};

if co==1
 
    
%%%% To keep the same "random" indexes from Beatriz's code  
   
    if spiral==0
        load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_out.mat','bs_begin','bs_end')
    else
        load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_aux_spiral.mat','bs_begin','bs_end')
    end

%%%% To generate random indexes from begining    
% %     if spiral==0
% %         NS_BE_P
% %     else
% %         NS_BE_S
% %     end  
% %     segmentb=hu{iii,:};
% %     segmente=hd{iii,:};
% %      rng('default') % set random seed for consistency
% %      gen=rng;
% %    
% %     for j = 1:5e4
% %         ix=randi(length(segmentb),1);
% %         segment=randi([round(segmentb(ix)+1000) round(segmente(ix)-5000)],1);
% %         begin3=segment;
% %         end3=floor(begin3+5*samplerate);
% %         bs_begin(iii,j)=begin3;
% %         bs_end(iii,j)=end3 ;
% %     end
         
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
%     
    freq_bl(iii,:,:)=fff;
    amp_bl(iii,:,:)=aaa;
%     load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/a1')
%     load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/f1')
    
    for j = 1:5e4 % in pca, rows are observations and columns are variables
        for_pca = squeeze(clt(:,j,:)); % should be 3 vs length(segments)
        [pc, score, latent, tsquare, explained] = pca(for_pca');
        out.x_all{iii,co}(j,1:3) = pc(1:3, 1);
        out.pc1_exp{iii,co}(j,1:3) = explained;
    end
    
else
    
    x=[];
    for j = 1:length(start{iii,co})
        if (~isnan(start{iii,co}(j)))
            for i=1:3
                x(i,:)= tfilt(i,(start{iii,co}(j):ending{iii,co}(j)));
            end
            [pc, score, latent, tsquare, explained] = pca(x'); clear x
            out.x_all{iii,co}(j,1:3) = pc(1:3, 1);
            out.pc1_exp{iii,co}(j,1:3) = explained;
        end
    end
end

end

