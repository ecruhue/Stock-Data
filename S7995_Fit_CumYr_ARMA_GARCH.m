bplist00�_WebMainResource�	
^WebResourceURL_WebResourceFrameName_WebResourceData_WebResourceMIMEType_WebResourceTextEncodingName_�https://collab.its.virginia.edu/access/content/attachment/baf13f30-073e-4377-85b6-08f9a53082f3/_anon_/a47346c5-4112-4ee2-ad0a-b4dbb4f22aae/S7995_Fit_CumYr_ARMA_GARCH.mPO	�<html><head></head><body><pre style="word-wrap: break-word; white-space: pre-wrap;">%



TDaysYr=253;
NumYr=16;
NumPredict=126;
% 6 month prediction
BigL=NumYr*TDaysYr+NumPredict;

% SPchosen is one of the numbers you were given, for a particular Type.
SPchosen=10;

% Let Tp be the Type: CD=1; Finance=2; Hcare=3; Indstrls=4
% You need to specify Tp
Tp=2;
if (Tp==1)
   BigData=Price_ACl_CD(1:BigL,SPchosen);
elseif (Tp==2)
   BigData=Price_ACl_Finance(1:BigL,SPchosen);
elseif (Tp==3)
    BigData=Price_ACl_Hcare(1:BigL,SPchosen);
elseif (Tp==4)
    BigData=Price_ACl_Indstrls(1:BigL,SPchosen);
end

Mdl=arima('ARlags',1:2,'D',1,'MAlags',1:2,'Variance',garch('GARCHlags',1,'ARCHlags',1));
%Mdl=arima('ARlags',1:2,'D',1,'MAlags',1:2,'Variance',garch('ARCHlags',1:5));

PredictMat=zeros(NumYr,NumPredict);
PredictCVar=zeros(NumYr,NumPredict);

for r=1:NumYr 
    IV=1:r*TDaysYr;
    IVL=length(IV);
    X=log(BigData(IV));
    IV1=1:r*TDaysYr+NumPredict;
    IVL1=length(IV1);
    Y=log(BigData(IV1));
    
    
    % plot acf, pacf for diff(X) and for squared.
    
    figure(1);
    S7995_acfpacfGARCH(diff(X),30);
    
    % Perform est under model Mdl
    [EstMdl,EstParamCov,logL,info]=estimate(Mdl,X);

    % look at the estimates, P-values to get a sense of the quality of the fit

    % Get the standardized residuals, fit and the cond variance from fit:
    [Res,V,LogL]=infer(EstMdl,X);
    SRes=Res./sqrt(V);
    
    figure(2);
    S7995_acfpacfGARCH(SRes,30);
    % acf, pacf, GARCH 
     
    figure(3);
    subplot(4,1,1);
    plot(X);
    title(['YrNumber: ',num2str(r),',','Data']);
    subplot(4,1,2);
    plot(diff(X));
    title(['YrNumber: ',num2str(r),',','(diff log) Data']);
    subplot(4,1,3);
    plot(V);
    title('Cond Variance');
    subplot(4,1,4);
    plot(SRes);
    title('Standardized Residuals');
    xlabel('time');

   
    [E0,V0]=infer(EstMdl,X);
    [X1,XMSE1,V1]=forecast(EstMdl,NumPredict,'Y0',X(1:IVL),'E0',E0,'V0',V0);

    PredictMat(r,:)=X1;
    PredictCVar(r,:)=V1;
    
    
    figure(4);
    subplot(2,1,1);
    plot(1:IVL1,Y,'k-');
    hold on
    plot(IVL+1:IVL+NumPredict,X1,'r-','Linewidth',2);
    plot(IVL+1:IVL+NumPredict,X1+1.96*sqrt(XMSE1),'k--','Linewidth',2);
    plot(IVL+1:IVL+NumPredict,X1-1.96*sqrt(XMSE1),'k--','Linewidth',2);
    hold off

    subplot(2,1,2);
    plot(V0,'k-');
    hold on
    plot(IVL+1:IVL+NumPredict,V1,'r-','LineWidth',2);
    hold off

    

    pause;
       
end








</pre></body></html>]text/x-matlabUUTF-8    ( 7 N ` v �>?
�                           