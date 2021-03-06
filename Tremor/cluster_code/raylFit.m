function [fitobj,goodness,output] = raylFit(x, y)
[xData, yData] = prepareCurveData( x, y );
%CREATEFIT(Y)
%  Create a fit.
%
%  Data for '1_sin' fit:
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 05-Sep-2019 17:45:42


%% Fit: '1_sin'.
% [xData, yData] = prepareCurveData( [], y );
% Set up fittype and options.
%  ft = fittype( 'a*exp(-((x).^2)/(2*b.^2))', 'independent', 'x', 'dependent', 'y' );
 ft = fittype( 'a*exp(-((x-b).^2)./(2*(c.^2)))');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [max(y) mean(x) 0.5];
opts.Upper = [max(y) mean(x)+2 0.75];
opts.Lower = [0  mean(x)-2 0.1];

 y=y';x=x';
 xData=x(find(~isnan(y))); yData=y(~isnan(y));
 
% Fit model to data.
[fitobj,goodness,output] = fit( xData, yData, ft, opts );
% Plot fit with data.
% figure( 'Name', '1_sin' );
% h = plot( fitresult, xData, yData );
 h2=plot(fitobj);
 set(h2,'LineWidth',1.5,'Color','r')

% legend( h, 'y', '1_sin', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
% ylabel( 'y', 'Interpreter', 'none' );



