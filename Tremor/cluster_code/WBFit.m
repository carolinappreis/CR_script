

function [fitobj,goodness,output] = WBFit(x, y)
  ft = fittype( 'a*b*x^(b-1)*exp(-a*x^b)' );
% %   ft = fittype( 'c*a*b*x^(b-1)*exp(-a*x^b)');  
% % %  b=k 
% % %  a=beta
% %  
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Upper = [Inf Inf];
opts.Lower = [0 0];
opts.StartPoint = [1 0];

% % % opts.Upper = [Inf Inf Inf];
% % % opts.Lower = [0 0 0];
% % % opts.StartPoint = [1 0 5];

 xData=x(find(~isnan(y))); yData=y(~isnan(y));
% Fit model to data.
[fitobj,goodness,output] = fit( xData, yData, ft, opts );
