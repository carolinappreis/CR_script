load stim
load nostim

%% this is meant to normalise the data but for some reason didnot work on the first
% patient data so I commented it out

% [transdat,lambda] = boxcox(1+nostimout);
% 
% d=1+nanmedian(stimout);
% trans_stim=((d.^lambda)-1)./lambda;
% 
% upperthreshold=prctile(transdat,99.7917);
% lowerthreshold=prctile(transdat,0.2083);
% 
% find(trans_stim>=upperthreshold | trans_stim<=lowerthreshold)

%% this is without normalising the data 

d=nanmedian(stimout);

upperthreshold=prctile(nostimout,99.7917); %bonferroni corrected for 12 comparisons
lowerthreshold=prctile(nostimout,0.2083);

find(d>=upperthreshold | d<=lowerthreshold)