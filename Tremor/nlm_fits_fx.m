
close all
clear all
load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\cleaned_rc12.mat')

cohort=[2 3 4 5 8 10 11 13 16 17];
dum=intersect(iiii,cohort);

pt=[];
for i=1:length(dum)
    pt=[pt find(iiii==dum(i))];
end

main=[1 1 3 1 3 3 3 3 1 1];
%   main=[1 1 1 1 1 1 1 1 1 1];
for pp=1:size(pt,2)
    for kk=1:size(tt1,2)
        
        curve=nanmedian(tt1{pt(pp),kk},1);
        SO=repmat(curve,1,3);
        for i=size(tt1{pt(pp),kk},2)+1:size(tt1{pt(pp),kk},2)*2
            smooth_c(1,i-12)=sum(SO(1,(i-1:i+1)))./length(SO(1,(i-1:i+1)));
        end
        smoo_all(pp,kk,:)=smooth_c;  clear smooth_c
    end

    smoo_main(pp,:)=squeeze(smoo_all(pp,1,:));
    s_main(pp,:)=nanmedian(tt1{pt(pp),1},1);
end


close all
for i=1:10
x = 1:12;
y = s_main(i,:);

beta0 = [(max(y)-min(y)) length(x)/2 2];
mg = @(F,x)(F(1)*exp(-((x-F(2))/F(3)).^2));
mdg = fitnlm(x,y,mg,beta0)


beta1 = [(max(y))];
ms = @(F,x)(F(1)*sin(0.5*x));
mds = fitnlm(x,y,ms,beta1)


subplot(2,5,i)
bar(x,y)
hold on
plot(x, mdg.Fitted,'k','LineWidth',2)
plot(x,mds.Fitted,'r','LineWidth',2)
box('off')
end
%%%%--------------------






