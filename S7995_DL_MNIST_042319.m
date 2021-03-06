bplist00�_WebMainResource�	
^WebResourceURL_WebResourceFrameName_WebResourceData_WebResourceMIMEType_WebResourceTextEncodingName_�https://collab.its.virginia.edu/access/content/attachment/baf13f30-073e-4377-85b6-08f9a53082f3/_anon_/ac3a1130-04d6-4843-9747-d3e4ff713dc9/S7995_DL_MNIST_042319.mPO�<html><head></head><body><pre style="word-wrap: break-word; white-space: pre-wrap;">%

[XTrain,YTrain] = digitTrain4DArrayData;
% digitTrain4DArray has 5000 images from MNIST
% Here, 4000 will be used to train and 1000 for validation

idx = randperm(size(XTrain,4),1000);
XValidation = XTrain(:,:,:,idx);
XTrain(:,:,:,idx) = [];
YValidation = YTrain(idx);
YTrain(idx) = [];

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


% options = trainingOptions('sgdm', ...
%     'InitialLearnRate',0.01, ...
%     'MaxEpochs',4, ...
%     'Shuffle','every-epoch', ...
%     'ValidationData',{XValidation,YValidation}, ...
%     'ValidationFrequency',30, ...
%     'Verbose',false, ...
%     'Plots','training-progress');


options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',4, ...
    'Shuffle','every-epoch', ...
    'Verbose',false, ...
    'Plots','training-progress');

net = trainNetwork(XTrain,YTrain,layers,options);


YClsfy=classify(net,XValidation);
accuracy=sum(YClsfy == YValidation)/numel(YValidation)


YPred = predict(net,XValidation);




</pre></body></html>]text/x-matlabUUTF-8    ( 7 N ` v �9:;I                           O