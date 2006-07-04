function response=impulsdt(A,nstep)
%function response=impulsdt(A,nstep)
% A should be nvar x nvar x nlag, from A0x(t)=A1x(t-1)+...+e(t).  Note that nlag is one
% more than the number lags in the reduced form var.
%
% Note that there is a routine impulsdtrf that more conveniently computes reduced form
% impulse responses.  
%
% If it is necessary to use this routine instead,
% starting with a reduced form coefficient matrix B(neqn,nvar,lags) like By from rfvar.m, 
% your input to this routine for non-orthoganalized responses is cat(3,eye(nvar), B).
% If you have a square root W of sigam, the r.f. covariance matrix (so W'*W=sigma), then
% orthogonalized impulse responses are found by using as input to this routine 
% reshape(W'\reshape(cat(3,eye(nvar),B),nvar,nvar*nlag),[nvar nvar nlag])
% 
%
[nvar,j,nlag]=size(A);
response=zeros(nvar,nvar,nstep);
A0=A(:,:,1);
if cond(A0)>1e200
   disp 'singular A0'
   response=[];
else
   A0i=inv(A0);
   Aplus=A(:,:,2:nlag);
   Aplus=A0\Aplus(:,:);
   Aplus=reshape(Aplus,[nvar,nvar,nlag-1]);
   response(:,:,1)=A0i;
   for it=1:nstep
      for ilag=1:min(nlag-1,it-1)
         response(:,:,it)=response(:,:,it)+Aplus(:,:,ilag)*response(:,:,it-ilag);
      end
   end
end
