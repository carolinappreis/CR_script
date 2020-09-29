function [low,high,m_arc_a,m_arc_s]=ns_th(out,ref_ns,iii,low,high,m_arc_a,m_arc_s)
[p,q]=(max(ref_ns'));
[m,n]=(min(ref_ns'));

parfor ii=1:size(ref_ns,1)
    % main peak at 0 deg
    if q(ii)==1;
        s_al_a(ii,:)=ref_ns(ii,:);
        check(1,ii)=NaN;
    else
        s_al_a(ii,:)=[ref_ns(ii,q(ii):end) ref_ns(ii,1:q(ii)-1)];
        %         check(ii,:)=[q(ii):size(s_al,2) 1:q(ii)-1];
        %         check(1,ii)=sum(diff([q(ii):size(s_al_a,2) 1:q(ii)-1]));
    end
    
    % main peak at 0 dega
    if n(ii)==1;
        s_al_s(ii,:)=ref_ns(ii,:);
        check(1,ii)=NaN;
    else
        s_al_s(ii,:)=[ref_ns(ii,n(ii):end) ref_ns(ii,1:n(ii)-1)];
        %         check(ii,:)=[q(ii):size(s_al,2) 1:q(ii)-1];
        %         check(1,ii)=sum(diff([n(ii):size(s_al_s,2) 1:n(ii)-1]));
    end
   
end


m_arc_a(iii,:)=median(s_al_a);
m_arc_s(iii,:)=median(s_al_s);
high(iii,:)=prctile(s_al_a,95);
low(iii,:)=prctile(s_al_s,5);

end

% % function [low,high]=ns_th(ref_ns,iii,low,high)
% % [p,q]=(max(ref_ns'));
% % for ii=1:size(ref_ns,1)
% %     % main peak at 0 deg
% %     if q(ii)==1;
% %         s_al(ii,:)=ref_ns(ii,:);
% %         check(1,ii)=NaN;
% %     else
% %         s_al(ii,:)=[ref_ns(ii,q(ii):end) ref_ns(ii,1:q(ii)-1)];
% %         %         check(ii,:)=[q(ii):size(s_al,2) 1:q(ii)-1];
% %         check(1,ii)=sum(diff([q(ii):size(s_al,2) 1:q(ii)-1]));
% %     end
% % end
% % high(iii,:)=prctile(s_al,95); clear s_al check q p ii
% %
% % [p,q]=(min(ref_ns'));
% % for ii=1:size(ref_ns,1)
% %     % main peak at 0 dega
% %     if q(ii)==1;
% %         s_al(ii,:)=ref_ns(ii,:);
% %         check(1,ii)=NaN;
% %     else
% %         s_al(ii,:)=[ref_ns(ii,q(ii):end) ref_ns(ii,1:q(ii)-1)];
% %         %         check(ii,:)=[q(ii):size(s_al,2) 1:q(ii)-1];
% %         check(1,ii)=sum(diff([q(ii):size(s_al,2) 1:q(ii)-1]));
% %     end
% % end
% % low(iii,:)=prctile(s_al,5);
% %
% % end