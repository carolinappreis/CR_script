function [number_ovl,dur_ovl]=overlap(dum)    
ovlap=(dum==2); 
      d_ovl=diff(ovlap);
      d1=find(d_ovl~=0);
      if d_ovl(d1(1))
          beging=find(d_ovl==1);
          ending=find(d_ovl==-1);
      else
          beging=find(d_ovl==-1);
          ending=find(d_ovl==1);
      end

        dur1=ending-beging;
        number_ovl=numel(find(dur1>50));
        dur_ovl=sum(dur1(find(dur1>50)));
        
end
        
        %       subplot(2,1,1)
%       plot(probe(1,:))
%       hold on
%       subplot(2,1,2)
%       plot(probe(ct,:))
%       hold on
%       for p=1:length(beging)
%           subplot(2,1,1)
%           plot(probe(1,:))
%           hold on
%           xline((beging(p)),'k')
%           xline((ending(p)),'k')
%           ylim([0 2])
%           
%           subplot(2,1,2)
%           plot(probe(ct,:))
%           hold on
%           xline((beging(p)),'k')
%           xline((ending(p)),'k')
%           ylim([0 2])
%       end
%  