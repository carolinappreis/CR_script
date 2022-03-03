function [idx,cl]=resolt_bins(numb_bins,pxx)

if numb_bins==4
    idx{1,1}=find(pxx>0 & pxx< 90);
    idx{2,1}=find(pxx<0 & pxx>-90);
    idx{3,1}=find(pxx>-180 & pxx<-90);
    idx{4,1}=find(pxx>90 & pxx<180);
    
    cl=[[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560]];
    
    
elseif numb_bins==6
    idx{1,1}=find(pxx>0 & pxx< 60);
    idx{2,1}=find(pxx<0 & pxx>-60);
    idx{3,1}=find(pxx<-60 & pxx>-120);
    idx{4,1}=find(pxx>-180 & pxx<-120);
    idx{5,1}=find(pxx>120 & pxx<180);
    idx{6,1}=find(pxx>60 & pxx< 120);
    
    
    cl=[[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];...
        [0, 0.6070, 0.7710];[0.4000, 0.3650, 0.0980]];
    
elseif numb_bins==8
    idx{1,1}=find(pxx>45 & pxx< 90);
    idx{2,1}=find(pxx>0 & pxx< 45);
    idx{3,1}=find(pxx<0 & pxx>-45);
    idx{4,1}=find(pxx<-45 & pxx>-90);
    idx{5,1}=find(pxx>-135 & pxx<-90);
    idx{6,1}=find(pxx>-180 & pxx<-135);
    idx{7,1}=find(pxx<180 & pxx>135);
    idx{8,1}=find(pxx<135 & pxx>90);
    
    cl=[[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];...
        [0, 0.6070, 0.7710];[0.4000, 0.3650, 0.0980];[0.490, 0.6970, 0.1250];[0.4940, 0.1840, 0.8060]];
end
end