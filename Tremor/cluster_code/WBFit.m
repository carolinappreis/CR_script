
function [fitobj,goodness,output] = WBFit(x, y, iii,spiral)
  
%  [xData, yData] = prepareCurveData( [], y );
 xData=x(find(~isnan(y))); yData=y(~isnan(y));
 
 

  
  %%% use this if 2 terms for all ft=fittype('(b/a)*((x/a)^(b-1))*exp(-((x/a)^b))');
   %    ft = fittype( 'a*b*x^(b-1)*exp(-a*x^b)' );
% %   ft = fittype( 'c*a*b*x^(b-1)*exp(-a*x^b)');  
% % %  b=k 
% % %  a=beta
% opts.Display = 'Off';
% opts.Upper = [Inf Inf];
% opts.Lower = [1 1];
% opts.StartPoint = [1 2];
% %  


if spiral==0 && iii==3
ft=fittype('c+(b/a)*((x/a)^(b-1))*exp(-((x/a)^b))');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Upper = [Inf Inf 1];
opts.Lower = [0 0 0];
opts.StartPoint = [1 0.5 0.2];
    else
ft=fittype('(b/a)*((x/a)^(b-1))*exp(-((x/a)^b))');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Upper = [Inf Inf];
opts.Lower = [1 1];
opts.StartPoint = [1 2];
end

% Fit model to data.

[fitobj,goodness,output] = fit( xData, yData, ft, opts );
