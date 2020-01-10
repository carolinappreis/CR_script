%run pls_sig numb=2 (pt 3); breakpoint 175

gg=[zenv(1,(start(1)-30000):(start(1)+60000))];
time=1:length(gg);
plot(time,gg(1,:),'Color',[0 0.5 0.5],'LineWidth',1)
hold on
for i=1:3
    c1=[time(segmentb(1)-1000) time(start(1)-5000) time(ending(1)-5000)]
    c2=[time(segmentb(1)) time(start(1)) time(ending(1))]
    x = [c1(i) c2(i) c2(i) c1(i)];
    y = [0 0 2.5 2.5];
    patch(x,y,[0.5 0.5 0.5],'FaceAlpha',[0.2],'EdgeColor','none')
end
box('off')

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
