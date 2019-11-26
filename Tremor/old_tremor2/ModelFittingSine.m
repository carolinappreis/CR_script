function ModResults = ModelFittingSine(RawData,ParamSet,stPlot)

% RawData is a vector: 1*N
% Modelling: Data = Amp*sin(2*pi*t + Phi)  (with a constant frequency of 2Hz/2 cycles)
% ModParam.Amp = Amp; ModParam.Phi = Phi;
ParamInit = ParamSet.Init;
LowParam = ParamSet.Min;
UpParam = ParamSet.Max;

Opts = optimset('DerivativeCheck','off','Display','off','TolX',1e-8,'TolFun',1e-8, ...
    'Diagnostics','off','MaxIter',10000);
ModParam = fmincon(@(Param)myfun(Param, RawData),ParamInit,[],[],[],[],LowParam,UpParam,[],Opts);
%ModParam = fmincon(@(Param)myfun(Param, ErrData, RotAng),ParamInit,ConstrA,ConstrB,[],[],LowParam,UpParam,[]);

nn = length(RawData);
tt = (1:nn)./nn;
FitData = ModParam(1).*sin(2*pi*ModParam(3)*tt + ModParam(2));

ModResults.FitData = FitData;
ModResults.ModParam = ModParam;
ModResults.VAF = 1-var(FitData-RawData)/var(RawData);

if stPlot == 1
    figure;
    plot(RawData,'b');
    hold on;
    plot(FitData,'r');
    legend('RawData','FitData');
end

end

function SqrDiff = myfun(Param, RawData)

% Param(1) = Amp; Param(2) = Phi;

nn = length(RawData);
tt = (0:nn-1)./nn;

XX = Param(1).*sin(2*pi*Param(3)*tt + Param(2));
 
Err = RawData - XX;

SqrDiff = sum(Err.^2); 

end


