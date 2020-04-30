function [fitobj,goodness,output] = poli_fit(x,y)

y=y';x=x';
xData=x(find(~isnan(y))); yData=y(~isnan(y));

ft = fittype( 'a*x^3 + b*x^2 + c*x + d');
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Display = 'Off';
opts.Upper = [Inf Inf Inf];
opts.Lower = [0 0 0];
opts.StartPoint = [max(y) nanmean(x) nanmean(x)+3];


[fitobj,goodness,output] = fitlm( xData, yData, ft );
h2=plot(fitobj);
set(h2,'LineWidth',1.5,'Color','k')




end

