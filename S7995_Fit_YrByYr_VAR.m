bplist00�_WebMainResource�	
^WebResourceURL_WebResourceFrameName_WebResourceData_WebResourceMIMEType_WebResourceTextEncodingName_�https://collab.its.virginia.edu/access/content/attachment/baf13f30-073e-4377-85b6-08f9a53082f3/_anon_/b3e64289-fcb5-48d6-a1a8-fc69f3bc8783/S7995_Fit_YrByYr_VAR.mPO
t<html><head></head><body><pre style="word-wrap: break-word; white-space: pre-wrap;">%



TDaysYr=253;
NumYr=16;
%NumPredict=126;
% 6 month prediction
NumPredict=20;

BigL=1+NumYr*TDaysYr+NumPredict;

% SPchosen is one of the numbers you were given, for a particular Type.
SPchosenV=[ 4 10 12];
L_S=length(SPchosenV);

% Let Tp be the Type: CD=1; Finance=2; Hcare=3; Indstrls=4
% You need to specify Tp
Tp=2;
if (Tp==1)
   BigData=diff(log(Price_ACl_CD(1:BigL,SPchosenV)));
elseif (Tp==2)
   BigData=diff(log(Price_ACl_Finance(1:BigL,SPchosenV)));
elseif (Tp==3)
    BigData=diff(log(Price_ACl_Hcare(1:BigL,SPchosenV)));
elseif (Tp==4)
    BigData=diff(log(Price_ACl_Indstrls(1:BigL,SPchosenV)));
end

% Here, in VAR, we need to directly take the diff of the log,
    % unlike in arima where it can be built in

p=3;
Mdl=varm(L_S,p); 

PredictMat=zeros(NumYr,L_S,NumPredict);


for r=1:NumYr 
    IV=(r-1)*TDaysYr+1:r*TDaysYr;
    IVL=length(IV);
    X=BigData(IV,:);
    nX=length(X);
    
    IV1=(r-1)*TDaysYr+1:r*TDaysYr+NumPredict;
    IVL1=length(IV1);
    Y=BigData(IV1,:);
    nY=length(Y);
    
    [EstMdl,EstParamCov,logL,info]=estimate(Mdl,X);
    
    %summarize(EstMdl)
    
    %NOTE: If you uncomment the line above, all the parameters will be written
    %out - they are worth seeing, at least occasionally
    
    
    
    
    Res=infer(EstMdl,X);
    
    
  


    [Forecast,ForecastMSE] = forecast(EstMdl,NumPredict,X);
    extractMSE = @(x)diag(x)';
    MSE = cellfun(extractMSE,ForecastMSE,'UniformOutput',false);
    SE = sqrt(cell2mat(MSE));

    ForecastFI = zeros(NumPredict,Mdl.NumSeries,2);
    ForecastFI(:,:,1) = Forecast - 2*SE;
    ForecastFI(:,:,2) = Forecast + 2*SE;

    subplot(3,1,1);
    h1 = plot(1:nY,Y(:,1),'-k');
    hold on;
    h2 = plot(nX+1:nX+NumPredict,Forecast(:,1),'r-','LineWidth',2);
    h3 = plot(nX+1:nX+NumPredict,ForecastFI(:,1,1),'k--','LineWidth',2);
    plot(nX+1:nX+NumPredict,ForecastFI(:,1,2),'k--','LineWidth',2);
    hold off
    
    subplot(3,1,2);
    h4 = plot(1:nY,Y(:,2),'-k');
    hold on;
    h5 = plot(nX+1:nX+NumPredict,Forecast(:,2),'r-','LineWidth',2);
    h6 = plot(nX+1:nX+NumPredict,ForecastFI(:,2,1),'k--','LineWidth',2);
    plot(nX+1:nX+NumPredict,ForecastFI(:,2,2),'k--','LineWidth',2);
    hold off 
    
    subplot(3,1,3);
    h7 = plot(1:nY,Y(:,3),'-k');
    hold on;
    h8 = plot(nX+1:nX+NumPredict,Forecast(:,3),'r-','LineWidth',2);
    h9 = plot(nX+1:nX+NumPredict,ForecastFI(:,3,1),'k--','LineWidth',2);
    plot(nX+1:nX+NumPredict,ForecastFI(:,3,2),'k--','LineWidth',2);    
    hold off  
    
    PredictMat(r,:,:)=Forecast';

    pause;
       
end








</pre></body></html>]text/x-matlabUUTF-8    ( 7 N ` v �89��                           �