function [ecog_shuf]=suffle_data(ctx)
dur=200;
segm=(reshape(ctx,dur,length(ctx)/dur))';
idx=randperm(length(ctx)/dur);
new=[];
for ii=1:length(idx)
   new=[new segm(idx(ii),:)];
end
ecog_shuf=new;

%%%% to check if duration of segments (dur) is enough to mantain the
%%%% spectral contect of the signal 
% % srn=1000;
% % [Pxx_ind,F_ind]=pwelch(ctx,srn,[],srn,srn);
% % plot(F_ind(1:50),Pxx_ind(1:50))
% % hold on
% % clear Pxx_ind F_ind
% % [Pxx_ind,F_ind]=pwelch(ecog_shuf,srn,[],srn,srn);
% % plot(F_ind(1:50),Pxx_ind(1:50),'r')
% % ylim([0 20000])
end