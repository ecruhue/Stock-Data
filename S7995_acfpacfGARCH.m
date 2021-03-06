bplist00�_WebMainResource�	
^WebResourceURL_WebResourceFrameName_WebResourceData_WebResourceMIMEType_WebResourceTextEncodingName_�https://collab.its.virginia.edu/access/content/attachment/baf13f30-073e-4377-85b6-08f9a53082f3/_anon_/b2429f54-47fa-411b-9e54-18175201e8a1/S7995_acfpacfGARCH.mPO><html><head></head><body><pre style="word-wrap: break-word; white-space: pre-wrap;">function S7995_acfpacfGARCH(X,m)



[d n]=size(X);
  if n==1
      n=d;
  end 
Y=X-mean(X);


subplot(2,2,1);
autocorr(Y,m);
title('sample acf - Data Not Squared','fontsize',14);
subplot(2,2,2);
parcorr(Y,m);
title('sample pacf - Data Not Squared','fontsize',14);

subplot(2,2,1)
axis([0 m -1 1])
subplot(2,2,2)
axis([0 m -1 1])
subplot(2,2,1)
set(gca,'fontsize',14)
subplot(2,2,2)
set(gca,'fontsize',14)

Z=Y.^2;

subplot(2,2,3);
autocorr(Z,m);
title('sample acf - Data Squared','fontsize',14);
subplot(2,2,4);
parcorr(Z,m);
title('sample pacf - Data Squared','fontsize',14);

subplot(2,2,3)
axis([0 m -1 1])
subplot(2,2,4)
axis([0 m -1 1])
subplot(2,2,3)
set(gca,'fontsize',14)
subplot(2,2,4)
set(gca,'fontsize',14)


% ac=xcorr(Y);
% %  two-sided and not divided by sample size yet
% r=ac(n:2*n-1)/n;
% alp=ones(m);
%  for k=1:m
%    a=levinson(r(1:k+2),k+1);
%    alp(k)=-a(k+1);
%  end
% cc=2/(n^.5);
% subplot(2,1,1);
% R=ones(m+1,1);
% R(2:m+1)=r(2:m+1)/r(1);
% plot(0:m,R,'*',1:m,cc*ones(1,m),'--',1:m,-cc*ones(1,m),'--','MarkerSize',12);
% % AX=axis;
% % MM=max([abs(AX(3)) abs(AX(4))]);
% axis([0 m+1 -1 1]);
% title('Sample ACF (lags \geq 0)','fontsize',14);
% xlabel('lag','fontsize', 14);
% A=ones(m+1,1);
% A(2:m+1)=alp(1:m);
% subplot(2,1,2);
% plot(0:m,A,'*',1:m,cc*ones(1,m),'--',1:m,-cc*ones(1,m),'--','MarkerSize',12);
% % AX=axis;
% % MM=max([abs(AX(3)) abs(AX(4))]);
% axis([0 m+1 -1 1]);
% title('Sample PACF (lags \geq 0)','fontsize',14);
% xlabel('lag','fontsize', 14);
</pre></body></html>]text/x-matlabUUTF-8    ( 7 N ` v �67y�                           �