%run pls_sig numb=2 (pt 3); breakpoint 175

gg=[zenv(1,(start(1)-30000):(start(1)+80000))];
time=1:length(gg);
plot(time,gg(1,:),'Color',[0 0.5 0.5],'LineWidth',1)
hold on
for i=1:3
    f1=figure(1)
    c1=[time(segmentb(1)-1000) time(start(1)-5000) time(ending(1)-5000)]
    c2=[time(segmentb(1)) time(start(1)) time(ending(1))]
    x = [c1(i) c2(i) c2(i) c1(i)];
    y = [0 0 3 3];
    patch(x,y,[ 1 0.8 0.5],'FaceAlpha',[0.2],'EdgeColor','none')
    xline(c2(i),'--','Color',[0.8 0 0 ],'LineWidth',1.5)
end

box('off')
xlabel('time (sec)')
xlim([0 100000])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
ylabel('tremor severity (m/s^2)')
set(gca,'FontSize',12)
f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 30, 8];
%run NS_sigmoid numb=2 (pt 3); breakpoint 124



gg=zenv(1,(segmentbp(1)-3500):(103700+25150+60000));
time=1:length(gg);
plot(time,gg(1,:),'Color',[0 0.5 0.5],'LineWidth',1)
hold on
for i=1:2
    c1=[time(segmentbp(1)-1000) (103700+25150-5000) ]
    c2=[(segmentbp(1)) (25150) ]
    x = [c1(i) c2(i) c2(i) c1(i)];
    y = [0 0 2.5 2.5];
    patch(x,y,[0.5 0.5 0.5],'FaceAlpha',[0.2],'EdgeColor','none')
end
box('off')
