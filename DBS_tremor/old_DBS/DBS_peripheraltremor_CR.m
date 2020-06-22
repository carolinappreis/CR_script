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

f1=figure(1)
subplot(1,3,1)
bar(0:30:330,100.*nanmedian(tt1{1,1})) 
hold on
plot(0:30:330,100.*tt1{1,1},'.')
box('off')
title('trials posture')
subplot(1,3,2)
bar(0:30:330,100.*nanmedian(tt1{1,1}))  
box('off')
title('arc posture')
subplot(1,3,3)
histogram(pca_ax{1,1}) 
box('off')
title('pca posture')
f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 40, 10];
set(f1,'color','w');

f2=figure(2)
subplot(1,3,1)
bar(0:30:330,100.*nanmedian(tt1{1,1})) 
hold on
plot(0:30:330,100.*tt1{2,1},'.')
box('off')
title('trials spiral')
subplot(1,3,2)
bar(0:30:330,100.*nanmedian(tt1{2,1}))  
box('off')
title('arc spiral')
subplot(1,3,3)
histogram(pca_ax{2,1})
box('off')
f2.Units = 'centimeters';
f2.OuterPosition= [10, 10, 40, 10];
set(f2,'color','w');
title('pca spiral')



% x axis phase - y axis percent change in 
% tremor severity at the end of 5 seconds with respect to severity right
% before stimulation began

