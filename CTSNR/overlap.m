function [number_ovl,dur_ovl]=overlap(dum,combi)
dumi=sum(dum);
if combi~=4
    ovlap=(dumi==2);
else
    ovlap=(dumi==3);
end

d_ovl=diff(ovlap);
d1=find(d_ovl~=0);

if ~isempty(d1)
    beging=find(d_ovl==1)+1;
    ending=find(d_ovl==-1);
    
    if d_ovl(d1(1))==-1
        ending=ending(2:end);
    end
    
    if d_ovl(d1(end))==1
        beging=beging(1:end-1);
    end
    
    %     time=1:length(dumi);
    %     plot(time,dumi)
    %     hold on
    %     plot(time(beging),dumi(beging)','r.')
    %     plot(time(ending),dumi(ending)','b.')
    idx_num=[];
    temp=zeros(1,size(dumi,2));
    for i=1:length(beging)
        if ending(i)-beging(i)>50
            temp(1,beging(i):ending(i))=1;
            idx_num(1,i)=i;
        end
    end
    
    dur_ovl=sum(temp);
    number_ovl=numel(idx_num);
else
    dur_ovl=NaN;
    number_ovl=NaN;
end
%%% old version
% % %       d_ovl=diff(ovlap);
% % %       d1=find(d_ovl~=0);
% % %       if d_ovl(d1(1))
% % %           beging=find(d_ovl==1);
% % %           ending=find(d_ovl==-1);
% % %       else
% % %           beging=find(d_ovl==-1);
% % %           ending=find(d_ovl==1);
% % %       end
% % %
% % %       if numel(ending)>numel(beging)
% % %           ending=ending(2:end);
% % %       elseif numel(ending)<numel(beging)
% % %           beging=beging(1:length(ending));
% % %       end
% % %
% % %         dur1=ending-beging;
% % %         number_ovl=numel(find(dur1>50));
% % %         dur_ovl=sum(dur1(find(dur1>50)));
% % %         number_ovl=numel(dur1);
% % %         dur_ovl=sum(dur1);
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