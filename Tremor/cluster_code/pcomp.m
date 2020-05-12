function [pc_trials,runs,start,ending]=pcomp(filt_tremor,start,ending,co,samplerate,iii)


if co==1
    %        baseline3 = NaN(3,5e4);
    for_cluster = NaN(3,5e4,5001);
    idx_b=zeros(1,5e4);
    idx_e=zeros(1,5e4);
    for j = 1:5e4
        ix=randi(length(start{iii,co}),1);
        segment=randi([start{iii,co}(ix)+1000 ending{iii,co}(ix)-5000],1);
        begin3=segment;
        end3=floor(begin3+5*samplerate);
        idx_b(1,j)= begin3;
        idx_e(1,j)= end3;
        for ax = 1:3
            %             baseline3(ax,j) = (mean(envelope(ax,end3-1000:end3))-mean(envelope(ax, begin3-1000:begin3)))./mean(envelope(ax, begin3-1000:begin3));
            for_cluster(ax,j,:) = filt_tremor(ax, begin3:end3);
        end
    end
    
    %
    for j = 1:5e4 % in pca, rows are observations and columns are variables
        seg_pc = squeeze(for_cluster(:,j,:)); % should be 3 vs length(segments)
        [pc, score, latent, tsquare] = pca(seg_pc');
        pc_trials(j, 1:3) = pc(1:3, 1);
    end
    
    start{iii,1}=[]; start{iii,1}= idx_b;
    ending{iii,1}=[];ending{iii,1}= idx_e;
    
else %_______________________________________________________________________
    
    for j = 1:length(start{iii,co})
        seg_pc = [filt_tremor(1,start{iii,co}(j):ending{iii,co}(j)); filt_tremor(2,start{iii,co}(j):ending{iii,co}(j)) ;filt_tremor(3,start{iii,co}(j):ending{iii,co}(j))];
        [pc, score, latent, tsquare] = pca(seg_pc');
        pc_trials(j, 1:3) = pc(1:3, 1);
    end
    
end

runs(1,co)=length(start{iii,co});

end