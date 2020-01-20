function [gof,fitresult1] = sin_fit(y)
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

%  Auto-generated by MATLAB on 05-Sep-2019 17:24:46


%% Fit: '1_sin'.
[xData, yData] = prepareCurveData( [], y );

% Set up fittype and options.
ft = fittype( 'a*sin(b*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.Lower = [-Inf 0.5 -Inf];
% opts.StartPoint = [0.5 0.5543 0.658594407192643];
% opts.Upper = [Inf 0.5 Inf];
opts.Lower = [-Inf 0.5 -0.2];
opts.StartPoint = [2 0.5543 0];
opts.Upper = [Inf 0.5 0.2];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', '1_sin' );
% h = plot( fitresult, xData, yData );
h1=plot(fitresult);
set(h1,'LineWidth',1.5,'Color','k')

% legend( h, 'y', '1_sin', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
% ylabel( 'y', 'Interpreter', 'none' );



