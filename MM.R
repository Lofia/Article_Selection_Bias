x=sort(rexp(20,2))
# if theta<rate v_tilde is a constant ???
MM=function(x,N,alpha,beta){
  n=length(x)
  pp=c(exp(-theta*x),0)
  p=pp[1:n]-pp[2:(n+1)]
  
  est=function(c,theta){
    a=rep(1,n);a[1]=a[1]+n*alpha
    b=(N-n)/(1-c)*p;b[n]=b[n]+n*beta
    
    M=matrix(nrow=n,ncol=n)
    for(i in 1:n){
      for(j in i:n) M[i,j]=sum(a[i:j])/sum(b[i:j])
    }
    #print(M)
    v=rep(Inf,n)
    v[1]=min(M[1,]);v[n]=max(M[,n])
    for(k in 2:(n-1)) v[k]=max(apply(M[1:k,k:n],1,min))
    #cat(sum(p*v))
    #cat('\n')
    return(list(v=v,v_tilde=unique(v),index=match(unique(v),v)))
  }
  
  bisection_c=function(left,right,theta,n=1){ # f increasing
    mid=(left+right)/2
    if(abs(left-right)<0.001) return(list(c=mid,iter=n))
    n=n+1
    fmid=mid-sum(p*est(mid,theta)$v)
    #cat(c(left,right,mid,fmid))
    #cat('\n')
    if(fmid>0.001) return(bisection_c(left,mid,theta,n))
    if(fmid<0.001) return(bisection_c(mid,right,theta,n))
  }
  
  bisection_theta=function(left,right,v){ # f increasing
    mid=(left+right)/2
    fmid=
      if(fmid>0.001) return(bisection_theta(left,mid,v))
    if(fmid<0.001) return(bisection_theta(mid,right,v))
    return(mid)
  }
  
  theta0=1/mean(x)
  repeat{
    c=bisection_c(0,1,theta0) ###gai
    v=est(c,theta0)
    theta=bisection_theta(v)
    
    if(0) return(list(v,theta))
    theta0=theta #;v0=v
  }
}


# g=1;an=alpha*n
# Fx=function(x) return(1-exp(-100*x))
# est3=function(x,g,an,beta,i,n){ # estimate each omega
#   outs=rep(NA,n-i+1)
#   if(i>1) an=0 # no need to add an
#   if(i==n) return(list(i0=n,out=(n-i+1+an)/(g*(1-Fx(x[i]))+beta)/n))
#   outs[1:(n-i)]=sapply(i:(n-1),function(j){return((j-i+1+an)/(g*(Fx(x[j+1])-Fx(x[i]))))}) # j changes
#   outs[n-i+1]=(n-i+1+an)/(g*(1-Fx(x[i]))+beta)
#   outn=min(outs)
#   i0=which(outs==outn)+i-1 # i0 is the index of n1+n2+n3+...+n_i
#   return(list(i0=i0,out=outn/n)) # out is the estimated omega_k
# }
# 
# pmle3=function(x,n,an,beta){ # estimate each omega_hats
#   omegas=index=vector()
#   iout=0;i=1;sum=0
#   while(i<=n){
#     iout=iout+1
#     result=est3(x,g,an,beta,i,n) # estimate omegas[k]=omega_k and index i0
#     i0=result[[1]];omegas[iout]=result[[2]]
#     
#     index[iout]=i0 # index[i] = n1+n2+n3+...+n_i
#     i=i0+1 # continue with the next group
#     #cat('omegas:',omegas,'index:',index,'i0:',i0,'i:',i,'iout:',iout,'\n')
#   }
#   return(list(omegas=omegas,index=index,iout=iout))
# }
