numgrad <- function(fcn,x,...){
  delta <- 1e-6
  n <- length(x)
### we tolerate x's that may be n x 1, 1 x n, or R vectors (with no dim)
### but note that g comes out as n x 1 vector regardless
  tvec <- delta*diag(n)
  g <- matrix(0,n,1)
  f0 <- fcn(x,...)
  badg <- FALSE
  for (i in 1:n){
    scale <- 1
    tvecv <- tvec[,i]
    if(is.null(dim(x))){
      tvecv <- as.vector(tvecv)
    }else{
      dim(tvecv) <- dim(x)
    }
    g0 <- (fcn(x+scale*tvecv,...) - f0)/(scale*delta)
    if (abs(g0)< 1e15){
      g[i] <- g0
    }else{
      cat("bad gradient ------------------------\n")
      g[i] <- 0
      badg <- TRUE
    }
  }
  return(list(g=g,badg=badg))
}
###-------------------------------------------------------------
###     if g0 > 0
###        sided=2;
###        g1 = -(eval([fcn '(x-scale*tvec(:,i)''' tailstr]) - f0) ...
###           /(scale*delta);
###        if g1<0
###           scale = scale/10;
###        else
###           break
###        end
###     else
###        sided=1;
###        break
###     end
###  end
###  if sided==1
###     g(i)=g0;
###  else
###     if (g0<1e20)
###        if (g1>-1e20)
###           g(i)=(g0+g1)/2;
###        else
###           g(i)=0;
###           badg=1;
###           disp( ['Banging against wall, parameter ' int2str(i)] );
###        end
###     else
###        if g1>-1e20
###           if g1<0
###              g(i)=0;
###              badg=1;
###              disp( ['Banging against wall, parameter ' int2str(i)] );
###           else
###              g(i)=g1;
###           end
###        else
###           g(i)=0;
###           badg=1;
###           disp(['Valley around parameter ' int2str(i)])
###        end
###     end
###  end
###end
###save g.dat g x f0
###eval(['save g g x f0 ' stailstr]);
