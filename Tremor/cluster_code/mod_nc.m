function [start,ending,out,yy]=mod_nc(start,ending,co,samplerate,iii,s,yy,out,gen)
% [pc_trials,start,ending,out,yy]=mod_nc(start,ending,co,samplerate,iii,s,yy,out)



if co==1
    for ax=1:3
        load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/new_cluster_CR.mat','nostim')
        out.change_nc{iii,co}(ax,:)=nostim(iii,ax,:);
    end
    

    
    %%%%% ------ below should match above
    
    %     rng(gen);
    %     envelope=s.env{iii,co};
    %     filt_tremor=s.filt{iii,co};
    %     baseline3 = NaN(3,5e4);
    %     for_cluster = NaN(3,5e4,5001);
    %     idx_b=zeros(1,5e4);
    %     idx_e=zeros(1,5e4);
    %
    %     for j = 1:5e4
    %         ix=randi(length(start{iii,co}),1);
    %         segment=randi([round(start{iii,co}(ix)+1000) round(ending{iii,co}(ix)-5000)],1);
    %         begin3=segment;
    %         end3=floor(begin3+5*samplerate);
    %         idx_b(1,j)= begin3;
    %         idx_e(1,j)= end3;
    %         for ax = 1:3
    %             baseline3(ax,j) = (mean(envelope(ax,end3-1000:end3))-mean(envelope(ax, begin3-1000:begin3)))./mean(envelope(ax, begin3-1000:begin3));
    %             for_cluster(ax,j,:) = filt_tremor(ax, begin3:end3);
    %         end
    %     end
    %
    %
    %
    %     for j = 1:5e4 % in pca, rows are observations and columns are variables
    %         seg_pc = squeeze(for_cluster(:,j,:)); % should be 3 vs length(segments)
    %         [pc, score, latent, tsquare] = pca(seg_pc');
    %         pc_trials(j, 1:3) = pc(1:3, 1);
    %     end
    %
    %     start{iii,1}=[]; start{iii,1}= idx_b;
    %     ending{iii,1}=[];ending{iii,1}= idx_e;
    %
    %     rep=10;
    %     for ax=1:3
    %         dum = baseline3(randi(length(baseline3), 1e6, rep));
    %         dum2 = dum;
    %         out.change_nc{iii,co}(ax,:)=nanmedian(dum2,2); clear dum dum2
    %     end
    
    
else %_______________________________________________________________________
    
    envelope=s.env{iii,co};
    filt_tremor=s.filt{iii,co};
    
    for j = 1:length(start{iii,co})
        seg_pc = [filt_tremor(1,start{iii,co}(j):ending{iii,co}(j)); filt_tremor(2,start{iii,co}(j):ending{iii,co}(j)) ;filt_tremor(3,start{iii,co}(j):ending{iii,co}(j))];
        [pc, score, latent, tsquare] = pca(seg_pc');
        pc_trials(j, 1:3) = pc(1:3, 1);
    end
    

tremor_or2=NaN(length(start{iii,co}),1);

for axx=1:3
    for i=1:length(start{iii,co})
        if (~isnan(start{iii,co}(i)))
            tremor_or2(axx,i,1)=(mean(envelope(axx,ending{iii,co}(i)-1000:ending{iii,co}(i)))-mean(envelope(axx,start{iii,co}(i)-1000:start{iii,co}(i))))/mean(envelope(axx,start{iii,co}(i)-1000:start{iii,co}(i)));
        else
            tremor_or2(axx,i,1)=NaN;
            tremor_pc(i,1)=NaN;
            yy{iii,co}(i)= NaN;
        end
    end
    
    tt=NaN(20,12);
    
    for i=1:12
        tt(1:sum(yy{iii,co}==i),i)=tremor_or2(axx,find(yy{iii,co}==i));
        tt(tt==0)=NaN;
    end
    out.change_nc{iii,co}{axx,1}=tt;
    
end

end

end