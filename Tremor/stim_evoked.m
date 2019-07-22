clear all
iii=[1 2 3 4 5 6 8 10 11];
for numb=5
%     1:length(iii);    
clearvars -except iii numb
load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))
% load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))

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
pre_stim_evoked;

time1=0:1/samplerate:(size(tremorx,2)-1)/samplerate;

% plot(time1,stim)
% hold on
% plot(time1(index),stim(index),'r.')
% plot(time1,tremorx)

epoch=250;
for st=1:length(index)
    if (~isnan(index(st)))
    evk_x(st,(1:2*epoch+1))= tremorxf(index(st)-epoch:index(st)+epoch);
    evk_y(st,(1:2*epoch+1))= tremoryf(index(st)-epoch:index(st)+epoch);
    evk_z(st,(1:2*epoch+1))= tremorzf(index(st)-epoch:index(st)+epoch);
    check(st,(1:2*epoch+1))= stim(index(st)-epoch:index(st)+epoch);
    emg2(st,(1:2*epoch+1))= emg1(index(st)-epoch:index(st)+epoch);
    end
end

clear ttx tty ttz ch


for i=1:12
    ttx(i,:)=mean(evk_x(find(xx==i),:)); 
    tty(i,:)=mean(evk_y(find(xx==i),:));
    ttz(i,:)=mean(evk_z(find(xx==i),:));
    ch(i,:)=mean(check(find(xx==i),:));
    emg3(i,:)=mean(emg2(find(xx==i),:));
    
    ttx(ttx==0)=NaN;
    tty(tty==0)=NaN;
    ttz(ttz==0)=NaN;
    ch(ch==0)=NaN;
    emg3(emg3==0)=NaN;
end



for i=1:12;
fig=figure(1)
subplot(5,1,1)
plot(ttx(i,:))
xline(250,'--')
% ylim([-0.5 0.5])
set(gca,'FontSize',12)
box('off')

subplot(5,1,2)
plot(tty(i,:))
xline(250,'--')
% ylim([-0.5 0.5])
set(gca,'FontSize',12)
box('off')

subplot(5,1,3)
plot(ttz(i,:))
xline(250,'--')
% ylim([-0.5 0.5])
set(gca,'FontSize',12)
box('off')

subplot(5,1,4)
plot(ch(i,:))
xline(250,'--')
set(gca,'FontSize',12)
box('off')

subplot(5,1,5)
plot(emg3(i,:))
xline(250,'--')
set(gca,'FontSize',12)
box('off')

fig.Units = 'centimeters';
fig.OuterPosition= [5, 5, 20, 27];
set(fig,'color','w');

% close all
end

end