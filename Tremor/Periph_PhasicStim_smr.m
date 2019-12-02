clear all
cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
SMR_File_To_Mat;


in2=1; % analysing the "main tremor axis"

if in2==1
    in=3;
elseif in2==2 % other axis 1
    in=5;
elseif in2==3 % other axis 2
    in=6;
end
data=SmrData.WvData;
samplerateold=SmrData.SR;
tremor=(data(in,:));
addon=92; addon_end=35;

Periph_PhaseDet;


for i=1
    a=hist(pca_ax{i,1},3);
    axmax=find(a==max(a));
    f2=figure(i)
    subplot(1,4,1)
    bar(0:30:330,100.*nanmedian(tt1{i,axmax}))
    hold on
    plot(0:30:330,100.*tt1{i,axmax},'.')
    set(gca,'XTickLabelRotation',45)
    box('off')
    subplot(1,4,2)
    bar(0:30:330,100.*nanmedian(tt1{i,axmax}))
    set(gca,'XTickLabelRotation',45)
    box('off')
    subplot(1,4,3)
    bar(1:3,a)
    names = {'CED2'; 'CED5';'CED6'};
    set(gca,'xtick',[1:3],'xticklabel',names)
    box('off')
    subplot(1,4,4)
    bar([1 2],[PSI_ax{i}(1,1);PSI_ax{i}(1,2)]);
    names = {'PSI 2 5'; 'PSI 2 6'};
    box('off')
    set(gca,'xtick',[1:2],'xticklabel',names)
    f2.Units = 'centimeters';
    f2.OuterPosition= [10, 10, 60, 10];
    set(f2,'color','w');
    if i==1
        title('condition1')
    else
        title('condition2')
    end
end


% x axis phase - y axis percent change in
% tremor severity at the end of 5 seconds with respect to severity right
% before stimulation began

