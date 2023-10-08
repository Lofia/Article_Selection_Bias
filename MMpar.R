####################
### MM Algorithm ###
####################
MM=function(x,N,alpha,beta){
  n=length(x)
  
  est=function(c,theta){
    a=rep(1,n);a[1]=a[1]+n*alpha
    b=(N-n)/(1-c)*p;b[n]=b[n]+n*beta
    
    M=matrix(nrow=n,ncol=n)
    for(i in 1:n) for(j in i:n) M[i,j]=sum(a[i:j])/sum(b[i:j])
    
    v=rep(Inf,n)
    v[1]=min(M[1,]);v[n]=max(M[,n])
    for(k in 2:(n-1)) v[k]=max(apply(M[1:k,k:n],1,min))
    #cat(sum(p*v))
    #cat('\n')
    return(list(v=v,v_tilde=unique(v),index=match(unique(v),v),M=M))
  }
  
  bisection_c=function(left,right,theta,iter=1){ # f increasing
    mid=(left+right)/2
    fmid=mid-sum(p*est(mid,theta)$v)
    if(abs(left-right)<0.0001) return(list(c=mid,iter=iter,fmid=fmid))
    iter=iter+1
    #cat(c(left,right,mid,fmid))
    #cat('\n')
    if(fmid>=0) return(bisection_c(left,mid,theta,iter))
    if(fmid<0) return(bisection_c(mid,right,theta,iter))
  }
  
  bisection_theta=function(left,right,v,iter=1){ # f increasing
    mid=(left+right)/2
    pp=c(exp(-mid*x),0)
    p=pp[1:n]-pp[2:(n+1)]
    ppp=c(x*exp(-mid*x),0)
    dp=ppp[2:(n+1)]-ppp[1:n]
    #cat('\n',length(p),length(dp),length(v),'\n')
    fmid=(N-n)/(1-sum(p*v))*sum(dp*v)-n/mid+sum(x)
    if(abs(left-right)<0.0001) return(list(theta=mid,iter=iter,fmid=fmid))
    iter=iter+1
    #cat(c(left,right,mid,fmid))
    #cat('\n')
    if(fmid>=0) return(bisection_theta(left,mid,v,iter))
    if(fmid<0) return(bisection_theta(mid,right,v,iter))
  }
  
  theta0=1/mean(x)
  loop=1
  repeat{
    pp=c(exp(-theta0*x),0)
    p=pp[1:n]-pp[2:(n+1)]
    c=bisection_c(0,1,theta0)
    v=est(c$c,theta0)
    theta=bisection_theta(0,5/mean(x),v$v)
    
    llf=sum(log(v$v))-theta$theta*sum(x)+n*log(theta$theta)+(N-n)*log(1-sum(p*v$v))
    # cat(loop,'\n',
    #     'c:',c$c,c$iter,c$fmid,'\n',
    #     'v:',v$v_tilde,'\n',
    #     'theta:',theta$theta,theta$iter,theta$fmid,'\n',
    #     'llf:',llf,'\n')
    if(abs(theta$theta-theta0)<0.001) return(list(v=v$v_tilde,theta=theta$theta,llf=llf))
    theta0=theta$theta #;v0=v
    loop=loop+1
  }
}


###################
### Power of MM ###
###################
library(pbapply)
# load('critical_value_for_power.Rda')
M=100
bb=seq(0,5,0.5)
n=c(20)
betas=matrix(nrow=length(bb),ncol=length(n))
for(j in 1:length(n)){
  for(i in 1:length(bb)){
    b=bb[i]
    cat('b =',b,'n =',n[j],'\n')
    temp=pbsapply(1:M,function(nouse){
      x=sort(rgamma(n[j],shape=b+1,rate=1))
      MME=MM(x,60,0.03,0.3)
      cri=quantile(pbsapply(1:100,function(nouse){
        return(MM(sort(rexp(n[j],MME$theta)),60,0.03,0.3)$llf)
      }),0.05)
      if(MME$llf<cri) return(1)
      return(0)})
    betas[i,j]=sum(temp)/M
  }
}
save(betas,file="MM_power.Rda")
