#Internal Functions

# log of posterior density
lpost=function(Y,crit,b){
  p=diff(pnorm(c(-Inf,crit,Inf)))
  sum(Y*log(p),na.rm=T)-sum(crit^2/2*b)
}

# Maximum likelihood estimate of criteria
crit.ml=function(y){
  I=length(y)
  N=sum(y)
  Fy=c(0,cumsum(y)/N)
  return(qnorm(Fy)[2:I])}

# External Functions

# MCMC chain 
doMCMC=function(y,M=100000,b=1,sdTune=.05){
  I=length(y)
  crit=matrix(nrow=M,ncol=I-1,0)
  crit[1,]=crit.ml(y)
  count=0
  
  for (m in 2:M){
    crit[m,]=crit[m-1,]
    cand=rnorm(I-1,crit[m,],sdTune)
    if (sum(diff(cand)<0)==0){
      lpost.cur=lpost(y,crit[m,],b)
      lpost.cand=lpost(y,cand,b)
      prob=min(exp(lpost.cand-lpost.cur),1)
      flip=rbinom(1,1,prob)
      if (flip){
        crit[m,]=cand
        count=count+1
      }
    }
  }
  return(list(crit=crit,accept=count/M))
}

# compute Bayes factor
stochdomBF=function(y1,y2,sdTune=c(.05,.05),M=100000,b=1){
  ml1=crit.ml(y1)
  ml2=crit.ml(y2)
  I=length(y1)
  out1=doMCMC(y1,M=M,b=b,sdTune=sdTune[1])
  out2=doMCMC(y2,M=M,b=b,sdTune=sdTune[2])
  diffMat=(out2$crit-out1$crit)>0
  a=apply(diffMat,1,mean)
  a1=mean(a==0)
  a2=mean(a==1)
  out=c((a1+a2)*I/2,a1*I,a2*I,out1$accept,out2$accept)
  names(out)=c("2Side","1>2","2>1","accept1","accept2")
  return(out)
}
