
function[seg_env,seg_pp,seg_var]=segmentation(idx,q,ns,env1,tempo1,filt_cs1,seg_env,seg_pp,seg_var,cl)



dseg=find(diff(idx{q,1})>1);

seg_env(ns,q,1)=mean(env1(idx{q,1}(1):idx{q,1}(dseg(1))));
seg_pp(ns,q,1)=peak2peak(filt_cs1(idx{q,1}(1):idx{q,1}(dseg(1))));
seg_var(ns,q,1)=var(filt_cs1(idx{q,1}(1):idx{q,1}(dseg(1))));

% figure(4)
% hold on
% plot(tempo1(idx{q,1}(1):idx{q,1}(dseg(1))),filt_cs1(idx{q,1}(1):idx{q,1}(dseg(1))),'.','Color',cl(q,:))


for ii=1:length(dseg)-1
    seg_env(ns,q,ii+1)=mean(env1(idx{q,1}(dseg(ii)+1):idx{q,1}(dseg(ii+1))));
    seg_pp(ns,q,ii+1)=peak2peak(filt_cs1(idx{q,1}(dseg(ii)+1):idx{q,1}(dseg(ii+1))));
    seg_var(ns,q,ii+1)=var(filt_cs1(idx{q,1}(dseg(ii)+1):idx{q,1}(dseg(ii+1))));
    
%     figure(4)
%     plot(tempo1(idx{q,1}(dseg(ii)+1):idx{q,1}(dseg(ii+1))),filt_cs1(idx{q,1}(dseg(ii)+1):idx{q,1}(dseg(ii+1))),'.','Color',cl(q,:))
%     
end
seg_env(ns,q,length(dseg)+1)=mean(env1(idx{q,1}(dseg(end)+1):idx{q,1}(end)));
seg_pp(ns,q,length(dseg)+1)=peak2peak(filt_cs1(idx{q,1}(dseg(end)+1):idx{q,1}(end)));
seg_var(ns,q,length(dseg)+1)=var(filt_cs1(idx{q,1}(dseg(end)+1):idx{q,1}(end)));


% figure(4)
% plot(tempo1(idx{q,1}(dseg(end)+1):idx{q,1}(end)),filt_cs1(idx{q,1}(dseg(end)+1):idx{q,1}(end)),'.','Color',cl(q,:))

end