bplist00�_WebMainResource�	
^WebResourceURL_WebResourceFrameName_WebResourceData_WebResourceMIMEType_WebResourceTextEncodingName_�https://collab.its.virginia.edu/access/content/attachment/baf13f30-073e-4377-85b6-08f9a53082f3/_anon_/e68ac13b-6396-4159-9750-1a419a4bf54e/S7995_DL_MNIST_042519.mPO�<html><head></head><body><pre style="word-wrap: break-word; white-space: pre-wrap;">%

%

C=randperm(10000);
NTrainD=9000;
XTrain=TrainD0to9Full(:,:,:,C(1:NTrainD));
XValidation=TrainD0to9Full(:,:,:,C(NTrainD+1:10000));
YTrain=TrainLabels_Full(C(1:NTrainD));
YValidation=TrainLabels_Full(C(NTrainD+1:10000));
% XTrain, YTrain, XValidation, YValidation are defined

layers = [
    imageInputLayer([28 28 1])
    
    convolution2dLayer(3,8,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(10)
    softmaxLayer
    classificationLayer];


options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',4, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{XValidation,YValidation}, ...
    'ValidationFrequency',30, ...
    'Verbose',false, ...
    'Plots','training-progress');


% options = trainingOptions('sgdm', ...
%     'InitialLearnRate',0.01, ...
%     'MaxEpochs',4, ...
%     'Shuffle','every-epoch', ...
%     'Verbose',false, ...
%     'Plots','training-progress');

net = trainNetwork(XTrain,YTrain,layers,options);


YClsfy=classify(net,XValidation);
accuracy=sum(YClsfy == YValidation)/numel(YValidation)


YPred = predict(net,XValidation);




</pre></body></html>]text/x-matlabUUTF-8    ( 7 N ` v �9:-;                           A