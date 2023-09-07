MM=function(x,N,alpha,beta){
  n=length(x)
  est=function(c,theta){
    
    return(v_tilde,index)
  }
  
  bisection_c=function(left,right,theta){ # f increasing
    mid=(left+right)/2
    fmid=
    if(fmid>0.001) return(bisection_c(left,mid,theta))
    if(fmid<0.001) return(bisection_c(mid,right,theta))
    return(mid)
  }
  
  bisection_theta=function(left,right,v){ # f increasing
    mid=(left+right)/2
    fmid=
      if(fmid>0.001) return(bisection_theta(left,mid,v))
    if(fmid<0.001) return(bisection_theta(mid,right,v))
    return(mid)
  }
  
  theta0=0
  repeat{
    c=bisection_c(0,1,theta0)
    v=est(c,theta0)
    theta=bisection_theta(v)
    
    if(0) return(list(v,theta))
    theta0=theta #;v0=v
  }
}