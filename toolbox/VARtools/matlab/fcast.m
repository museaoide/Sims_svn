function  yhat = fcast(y0,By,Bx,xdata,horiz)
%function  yhat = fcast(y0,By,Bx,xdata,horiz)
%yhat: (horiz+lags) x nvar
%By: equations x variables x lags
%Bx: equations x nx
%xdata: (lags+horiz) x nx
%y0: lags x nvar
% produces horiz-step-ahead forecast
[lags,nvar]=size(y0);
nx=size(Bx,2);
yhat=zeros(horiz+lags,nvar);
yhat(1:lags,:)=y0;
% By is equations x variables x lags
Bmat=permute(By,[3,2,1]);
Bmat=flipdim(Bmat,1);
Bmat=reshape(Bmat,lags*nvar,nvar);
Bx=Bx;
for it=1:horiz
    ydata=yhat(it:it+lags-1,:);
    yhat(lags+it,:)=sum(Bmat.*repmat(ydata(:),1,nvar))+xdata(lags+it,:)*Bx';
end
