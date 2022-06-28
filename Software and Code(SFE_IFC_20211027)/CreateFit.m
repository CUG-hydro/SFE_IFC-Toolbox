function [fit1,fit2,fit3,gof1,gof2,gof3] = CreateFit(COR_MRCdata)
%CreateFit 针对校正后的MRC数据进行曲线拟合
%   此处显示详细说明
[fit1,fit2,fit3,gof1,gof2,gof3]=deal([]);
xData = COR_MRCdata(:,1);
yData = COR_MRCdata(:,2);
opts1 = fitoptions('Method','NonlinearLeastSquares');
  opts1.StartPoint = 0;%max(yData)];
  ft1 = fittype( @(a,  x )yData(1)*exp(-a*x), 'independent', 'x', ...
          'dependent', 'y','Options',opts1); 
 % Fit model to data.
[fit1, gof1] = fit(xData, yData, ft1);

 opts2 = fitoptions('Method','NonlinearLeastSquares');
 opts2.StartPoint = [0.01 1.5];% max(yData)];
 opts2.Lower = [0 1.001 ];
 opts2.Upper = [1000 1.999];

 
  ft2 = fittype( @(a, b, x )(((yData(1))^(1-b))-(1-b)*a*x).^(1/(1-b)), 'independent', {'x'}, ...
              'dependent', {'y'}, 'coefficients',{'a','b'},'Options',opts2);
    
 [fit2, gof2] = fit(xData, yData, ft2);

opts3 = fitoptions('Method','NonlinearLeastSquares','MaxFunEvals',1000,'MaxIter',1000);    
opts3.Lower = [0.0001 0.0001];% 0];
opts3.StartPoint = [1  0.005];% max(yData)];
opts3.Upper = [10000 0.999];% 5000]; 

 ft3 = fittype(@(a,b,x) yData(1)*(1+(((1-b)*yData(1)^(1-b)).*x./(a*b))).^(1/(b-1)), 'independent', {'x'}, ...
           'dependent', 'y', 'coefficients',{'a','b'},'Options',opts3);
[fit3, gof3] = fit(xData, yData, ft3);

end

