bplist00�_WebMainResource�	
^WebResourceURL_WebResourceFrameName_WebResourceData_WebResourceMIMEType_WebResourceTextEncodingName_�https://collab.its.virginia.edu/access/content/attachment/baf13f30-073e-4377-85b6-08f9a53082f3/_anon_/a31a00c3-a5ee-444c-9c16-8a41e5e878d0/S7995_PCA_042719.mPOt<html><head></head><body><pre style="word-wrap: break-word; white-space: pre-wrap;">%




[coeff,score,latent]=pca(train0to9N);
MnN=mean(train0to9N,1);



TRvec=zeros(28,28,6,10000);
 


for k1=1:10000

TRvec(:,:,1,k1)=reshape(MnN+score(k1,1)*coeff(:,1)',28,28);
TRvec(:,:,2,k1)=reshape(MnN+score(k1,1:3)*coeff(:,1:3)',28,28);
TRvec(:,:,3,k1)=reshape(MnN+score(k1,1:10)*coeff(:,1:10)',28,28);
TRvec(:,:,4,k1)=reshape(MnN+score(k1,1:20)*coeff(:,1:20)',28,28);
TRvec(:,:,5,k1)=reshape(MnN+score(k1,1:50)*coeff(:,1:50)',28,28);
TRvec(:,:,6,k1)=reshape(MnN+score(k1,1:100)*coeff(:,1:100)',28,28);

end










</pre></body></html>^text/x-objcsrcUUTF-8    ( 7 N ` v �45��                           �