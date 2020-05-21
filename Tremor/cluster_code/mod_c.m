function [out,s]=mod_c(clust,s,samplerate,out,yy,start,ending,co,iii)

if co==1
    rng('default') % set random seed for consistency
    envelope=s.env{iii,co};
    
    
    st=start{iii,co}(clust.idx{iii,co});
    
    
    baseline3 = NaN(3,size(st,2));
    for j = 1:size(st,2)
        begin3=st(j);
        end3=floor(begin3+5*samplerate);
        for ax = 1:3
            baseline3(ax,j) = (mean(envelope(ax,end3-1000:end3))-mean(envelope(ax, begin3-1000:begin3)))./mean(envelope(ax, begin3-1000:begin3));
        end
    end
    rep=10;
    for ax=1:size(s.raw,1)
        dum = baseline3(randi(length(baseline3), 1e6, rep));
        out.change{iii,co}(ax,:)=nanmedian(dum,2); clear dum
    end
    
else
    
    envelope=s.env{iii,co};
    st=start{iii,co}(clust.idx{iii,co});
    et=ending{iii,co}(clust.idx{iii,co});
    yyt=yy{iii,co}(clust.idx{iii,co});
    
    
    tremor_or2=NaN(3,length(st));
    
    for axx=1:3
        for i=1:length(st)
            if (~isnan(st(i)))
                tremor_or2(axx,i)=(mean(envelope(axx,et(i)-1000:et(i)))-mean(envelope(axx,st(i)-1000:st(i))))/mean(envelope(axx,st(i)-1000:st(i)));
            else
                tremor_or2(axx,i)=NaN;
                yyt(i)= NaN;
            end
        end
        
        tt=NaN(20,12);
        
        for i=1:12
            tt(1:sum(yyt==i),i)=tremor_or2(axx,find(yyt==i));
            tt(tt==0)=NaN;
        end
        out.change_nc{iii,co}{axx,1}=tt;
        
    end
    
    
end
end
