bplist00�_WebMainResource�	
^WebResourceURL_WebResourceFrameName_WebResourceData_WebResourceMIMEType_WebResourceTextEncodingName_�https://collab.its.virginia.edu/access/content/attachment/baf13f30-073e-4377-85b6-08f9a53082f3/_anon_/27329c59-a848-4bea-aa29-71c26d36812b/S7995_LSTM_MTS_CumYr_041119.mPO/<html><head></head><body><pre style="word-wrap: break-word; white-space: pre-wrap;">%


TDaysYr=253;
NumYr=16;
NumPredict=126;
% 6 month prediction

BigL=NumYr*TDaysYr+NumPredict;


% 400 firms in 2000-2017 SP500 data-SPchosenV identifies a subset of interest
% 87 firms in 2000-2017 SP100 data-SPchosenV identifies a subset of interest

% Let Tp be the Type: CD=1; Finance=2; Hcare=3; Indstrls=4
% You need to specify Tp
Tp=2;

% For S7995-the 87 firms were reduced to 4 Types-14,14,14,13 in number
% SPchosen be the numbers you were given, for a particular Type.
SPchosenV=[4 10 12];
L_S=length(SPchosenV);

if (Tp==1)
   BigData=Price_ACl_CD(1:BigL,SPchosenV)';
elseif (Tp==2)
   BigData=Price_ACl_Finance(1:BigL,SPchosenV)';
elseif (Tp==3)
    BigData=Price_ACl_Hcare(1:BigL,SPchosenV)';
elseif (Tp==4)
    BigData=Price_ACl_Indstrls(1:BigL,SPchosenV)';
end


% NOTE: For multiple firms, each firm is in a row (not column as before)



%the following 2 are to be chosen - some choices better than others
NumInPerResp=3;
NumHiddenUnits = 200;


% figure(1);
% plot(BigData');

NumResponses=L_S;
NumInput=NumInPerResp*NumResponses;
NumFeatures = NumInput;



PredictMat=zeros(NumYr,L_S,NumPredict);
PredictRMSE=zeros(NumYr,L_S,1);


% mypool=parpool;

% parfor r=1:NumYr 
for r=1:NumYr 
    BData=BigData;
    IV=1:r*TDaysYr+NumPredict;
    IVL=length(IV);
    data=BData(:,IV);

     
    NumTimeStepsTrain=IVL-NumPredict;
   

    dataTrain = data(:,1:NumTimeStepsTrain+1);
    dataTest = data(:,NumTimeStepsTrain+1:end);
    Nd=length(dataTrain(1,:));
    Nt=length(dataTest(1,:));

    dataTrainStandardized=zeros(L_S,Nd);
    dataTestStandardized=zeros(L_S,Nt);

    muV=zeros(L_S,1);
    sigV=zeros(L_S,1);

    for jj=1:L_S
       muV(jj)=mean(dataTrain(jj,:));
       sigV(jj)=std(dataTrain(jj,:));
       dataTrainStandardized(jj,:) = (dataTrain(jj,:) - muV(jj))/sigV(jj);
       dataTestStandardized(jj,:) = (dataTest(jj,:) - muV(jj))/sigV(jj);
    end


    XTrainN=zeros(L_S,NumInPerResp,Nd-NumInPerResp);
    XTestN=zeros(L_S,NumInPerResp,Nt-NumInPerResp);
    
    for gg=1:L_S
      for k=1:NumInPerResp
        XTrainN(gg,k,:)=dataTrainStandardized(gg,NumInPerResp-k+1:Nd-k);
        XTestN(gg,k,:)=dataTestStandardized(gg,NumInPerResp-k+1:Nt-k);
      end
    end
    
    XTrain=reshape(XTrainN,L_S*NumInPerResp,Nd-NumInPerResp);
    XTest=reshape(XTestN,L_S*NumInPerResp,Nt-NumInPerResp);
    
    YTrain = dataTrainStandardized(:,NumInPerResp+1:Nd);
%     YTest = dataTestStandardized(:,NumInPerResp+1:Nt);
    
    



layers = [ ...
sequenceInputLayer(NumFeatures)
lstmLayer(NumHiddenUnits)
fullyConnectedLayer(NumResponses)
regressionLayer];

%The following two, I set (different from Matlab Doc)
MEpoch=125;
LRDropP=60;

options = trainingOptions('adam', ...
'MaxEpochs',MEpoch, ...
'GradientThreshold',1, ...
'InitialLearnRate',0.005, ...
'LearnRateSchedule','piecewise', ...
'LearnRateDropPeriod',LRDropP, ...
'LearnRateDropFactor',0.2, ...
'Verbose',0);

net = trainNetwork(XTrain,YTrain,layers,options);





net = predictAndUpdateState(net,XTrain);
[net,YPredN] = predictAndUpdateState(net,reshape(YTrain(:,end:-1:end-NumInPerResp+1),NumInput,1));

YPredNN=[YPredN ; reshape(YTrain(:,end:-1:end-NumInPerResp+2),NumInput-NumResponses,1)];

numTimeStepsTest =NumPredict;
for i = 2:numTimeStepsTest
    [net,YPredN] = predictAndUpdateState(net,YPredNN(:,i-1),'ExecutionEnvironment','cpu');
    YPredNN(:,i)=[YPredN ; YPredNN(1:NumInput-NumResponses,i-1)];
end


YTest=dataTest(:,NumInPerResp+1:Nt);

PMat=zeros(L_S,NumPredict);
PRMSE=zeros(L_S,1);
for jj=1:L_S
    PMat(jj,:) = sigV(jj)*YPredNN(jj,:)+muV(jj);
    PRMSE(jj,1)=sqrt(mean((PMat(jj,1:Nt-NumInPerResp)-YTest(jj,:)).^2));
end


PredictMat(r,:,:)=PMat;
PredictRMSE(r,:)=PRMSE;
    
        
figure(1);
for k=1:L_S
subplot(L_S,1,k);    
plot(IV,data(k,:),'-b');
hold on
    IVV=r*TDaysYr+1:r*TDaysYr+NumPredict;
    plot(IVV,reshape(PredictMat(r,k,:),NumPredict,1),'.-r');   
hold off
xlabel("time (days)");
ylabel("Stock Price");
title("Forecast");
%legend(["Observed" "Forecast"]);
end



pause;


end 


% delete(mypool);











</pre></body></html>]text/x-matlabUUTF-8    ( 7 N ` v �?@s�                           �