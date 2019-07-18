clear all
iii=[1 2 3 4 5 6 8 10 11];
PC=[70 66 47 47 47 50 50 50 50];
A1={([1 3 6 8 12 18 23 27 30 32]);[];[];([2 3 5 6 7]);([1 2 4 6 8 9 10 11]);[];[];([1:9 15]);([2 4 7:10 13:15 22 25])};
B1={([2 5 7 11 17 22 26 29 31 34]);[];[];([2 3 5 6 7]);([1 2 4 6 7 9 10 11 12]);[];[];([1:9 15]);([2 5 7 8 9 12 13 14 19 22 25])};
for numb=5
    %1:length(iii);
% load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))
load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))

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


% [b,a]=butter(2,[1/(0.5*samplerate) 15/(0.5*samplerate)],'bandpass'); %15
% tremorx1=filtfilt(b,a,tremorx);
% tremory1=filtfilt(b,a,tremory);
% tremorz1=filtfilt(b,a,tremorz);
stim=data(2,:);

epoch=250;
for st=1:length(start)
    if (~isnan(start(st)))
    evk_x(st,(1:2*epoch+1))= tremorxf(start(st)-epoch:start(st)+epoch);
    evk_y(st,(1:2*epoch+1))= tremoryf(start(st)-epoch:start(st)+epoch);
    evk_z(st,(1:2*epoch+1))= tremorzf(start(st)-epoch:start(st)+epoch);
%     check(st,(1:2*epoch+1))= stim(start(st)-epoch:start(st)+epoch);
    end
end

clear ttx tty ttz ch


for i=1:12
    ttx(i,:)=mean(evk_x(find(xx==i),:)); 
    tty(i,:)=mean(evk_y(find(xx==i),:));
    ttz(i,:)=mean(evk_z(find(xx==i),:));
%     ch(i,:)=mean(check(find(xx==i),:));
    
    ttx(ttx==0)=NaN;
    tty(tty==0)=NaN;
    ttz(ttz==0)=NaN;
%     ch(ch==0)=NaN;
end



for i=1:12;
fig=figure(1)
subplot(3,1,1)
plot(ttx(i,:))
hold on
xline(250,'--')
% ylim([-0.5 0.5])
set(gca,'FontSize',12)
box('off')

subplot(3,1,2)
plot(tty(i,:))
xline(250,'--')
% ylim([-0.5 0.5])
set(gca,'FontSize',12)
box('off')

subplot(3,1,3)
plot(ttz(i,:))
xline(250,'--')
% ylim([-0.5 0.5])
set(gca,'FontSize',12)
box('off')

% subplot(4,1,4)
% plot(ch(i,:))
% xline(250,'--')
% set(gca,'FontSize',12)
% box('off')

fig.Units = 'centimeters';
fig.OuterPosition= [5, 5, 8, 15];
set(fig,'color','w');

close all
end

end