clear all
cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\raw_data')
SMR_File_To_Mat;

in2=1; % analysing the "main tremor axis" 

if in2==1
    in=3;
elseif in2==2 % other axis 1
    in=6;
elseif in2==3 % other axis 2
    in=7;
end
data=SmrData.WvData;
samplerateold=SmrData.SR;
tremor=(data(in,:));
addon=92; addon_end=35;

DBS_phasedetection;


for i=1:2
a=hist(pca_ax{i,1},3);
axmax=find(a==max(a));
f2=figure(i)
subplot(1,3,1)
bar(0:30:330,100.*nanmedian(tt1{i,axmax})) 
hold on
plot(0:30:330,100.*tt1{i,axmax},'.')
box('off')
subplot(1,3,2)
bar(0:30:330,100.*nanmedian(tt1{i,axmax}))  
box('off')
subplot(1,3,3)
histogram(pca_ax{i,1})
box('off')
f2.Units = 'centimeters';
f2.OuterPosition= [10, 10, 40, 10];
set(f2,'color','w');
if i==1
title('posture')
else
title('spiral')
end
end


% x axis phase - y axis percent change in 
% tremor severity at the end of 5 seconds with respect to severity right
% before stimulation began

