clear all
cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA')
load ('01_RS_PS.mat')

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

bar(0:30:330,100.*nanmedian(tt)) 
hold on
plot(0:30:330,100.*tt,'.')
% x axis phase - y axis percent change in 
% tremor severity at the end of 5 seconds with respect to severity right
% before stimulation began


