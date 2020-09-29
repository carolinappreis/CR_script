
for iii=1:10
    
    bb=out.change_c{iii,1}(1,:);
    
    parfor i = 1:1000000
        dum = bb(randi(length(bb),1,12));  
        q=find(dum==max(dum));
        if q==1
         s_al_a(i,:) =dum  ;
        else
         s_al_a(i,:)=[dum(q:end) dum(1:q-1)];   
        end
    end
    
    
   pp(iii,:)=prctile(s_al_a,95);
   clear bb s_al_a
end


load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat','out');


for ii=1:10
        dist=out.change_c{ii,1}(1,:);
    for g=1:12
        maxs=high(ii,g);
        prc(ii,g)=numel(find(dist<maxs))/numel(dist)*100;
        clear max
    end
    clear max dist
end
