function [tocomp tocomp_r]=pref_slip(epochs_z1,b_ps,a_ps,after_on,r,tocomp,a_r,b_r,tocomp_r)
ini=350;
slip=epochs_z1(:,ini:400);
bp=b_ps(:,ini:400);
ap=a_ps(:,ini:400);
br=b_r(:,ini:400);
ar=a_r(:,ini:400);
mmatch=NaN(size(slip,1),2);
rmatch=NaN(size(slip,1),2);
for i=1:size(slip,1)
    dum=find(slip(i,:)==1);
    if ~isempty(dum)
    mmatch(i,1)=circ_mean((bp(i,dum))');
    mmatch(i,2)=circ_mean((ap(i,dum))');
    rmatch(i,1)=nanmean(br(i,dum));
    rmatch(i,2)=nanmean(ar(i,dum));
    end
    clear dum  
end

for ii=1:2
comb(ii,1)=circ_mean(mmatch(find(~isnan(mmatch(:,ii))),ii));
end
comb(3,1)=circ_mean(after_on(find(~isnan(mmatch(:,ii))),1));

tocomp(r,:)=[comb(3)-comb(1) comb(3)-comb(2)];

tocomp_r(r,:)=nanmean(rmatch);

end