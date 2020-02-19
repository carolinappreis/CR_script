
close all
clear all
load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\cleaned_rc12_noaddon.mat')
load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
cl=blushred;

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

cl=[0.5 0.5 0.5];
close all
for i=1:10
    
x = 1:12;
y = s_main(i,:);


beta0 = [(max(y)-min(y)) length(x)/2 1];
mg = @(F,x)(F(1)*exp(-((x-F(2))/F(3)).^2));
mdg = fitnlm(x,y,mg,beta0)


beta1 = [(max(y)) length(x)/2];
ms = @(F,x)(F(1)*sin(0.5*x+F(2)));
mds = fitlm(x,y,ms,beta1)

f1=figure(1);
subplot(2,5,i)
bar(0:30:330,y,'FaceColor',cl,'EdgeColor',cl)
hold on
plot(0:30:330,mdg.Fitted,'k','LineWidth',2)
plot(0:30:330,mds.Fitted,'r','LineWidth',2)
ylabel('Change in tremor severity')
xlabel('Stimulation phase (degrees)')
set(gca,'XTickLabelRotation',45)
set(gca,'FontSize',12)
box('off')
end
%%%%--------------------
f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 55, 15];
set(f1,'color','w');





