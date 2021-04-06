

function [fitobj,goodness,output] = WBFit(x, y)
[xData, yData] = prepareCurveData( x, y );
fo = fitoptions('Method','NonlinearLeastSquares', 'StartPoint',[0 0])
%     'Lower',[0  mean(x)-2 ],...
%      'Upper',[max(y) mean(x)+1 ],...
% 'StartPoint',[0 0]);
% ft = fittype('a*b*x^(b-1)*exp(-a*x^b)','options',fo);
ft = fittype('x./(b^2))*exp((x^2)./2*b^2))','options',fo);

[fitobj,goodness,output] = fit( xData, yData, ft);

h2=plot(fitobj);
set(h2,'LineWidth',1.5,'Color','r')
