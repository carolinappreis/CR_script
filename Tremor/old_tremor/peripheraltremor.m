clear all
cd('C:\Users\creis\Documents\MATLAB\pt_data_periphstim')
load ('P013_randomstim_cursos.mat')
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

phasedetection;

bar(0:30:330,100.*nanmedian(tt)) % x axis phase - y axis percent change in 
% tremor severity at the end of 5 seconds with respect to severity right
% before stimulation began

stimout=tt;
save stim stimout