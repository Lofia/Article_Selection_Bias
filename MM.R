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

# x=sort(rexp(300,5))
# 1/mean(x)
# MM(x,350,0.03,0.3)


#########################
### real MM Algorithm ###
#########################
GMM=function(x,N,alpha,beta){
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


#############################
## 3D Critical Value Table ## x x x
#############################
library(data.table)
B=CJ(n=c(10,20,50),alpha=c(0.03,0.1,0.3),theta=c(0.5,1,2))
B$value=rep(Inf,nrow(B))
library(pbapply)
for(i in 1:nrow(B)){
  B[i,'value']=quantile(pbsapply(1:100,function(nouse){
    cat('n=',as.numeric(B[i,'n']),' alpha=',as.numeric(B[i,'alpha']),' theta=',as.numeric(B[i,'theta']),sep='')
    return(MM(sort(rexp(as.numeric(B[i,'n']),as.numeric(B[i,'theta']))),60,as.numeric(B[i,'alpha']),0.3)$llf)
  }),0.05)
}
save(B,file="critical_value_table.Rda")
B$value=round(B$value,2)
# library(plotly)
fig=plot_ly(B,x=~n,y=~alpha,z=~theta,text=~value,mode='text',type='scatter3d') %>% 
  layout(#title='3D Critical Value Table',
    scene=list(xaxis=list(title='sample size n',type="category"),
               yaxis=list(title='constant alpha',type="category"),
               zaxis=list(title='parameter theta',type="category")))
fig


#########################################
## Critical Values for Computing Power ## x x x
#########################################
library(data.table)
Bp=CJ(n=c(10,20,50),alpha=c(0.03),theta=c(1))
Bp$value=rep(Inf,nrow(Bp))
library(pbapply)
for(i in 1:nrow(Bp)){
  cat('n =',Bp$n[i],'\n')
  Bp[i,'value']=quantile(pbsapply(1:10000,function(nouse){
    return(MM(sort(rexp(as.numeric(Bp[i,'n']),as.numeric(Bp[i,'theta']))),60,as.numeric(Bp[i,'alpha']),0.3)$llf)
  }),0.05)
}
save(Bp,file="critical_value_for_power.Rda")


###################
### Power of MM ###
###################
library(pbapply)
# load('critical_value_for_power.Rda')
M=50
bb=seq(0,5,0.5)
n=c(10,20,50)
betas=matrix(nrow=length(bb),ncol=length(n))
for(j in 1:length(n)){
  for(i in 1:length(bb)){
    b=bb[i]
    cat('b =',b,'n =',n[j],'\n')
    temp=pbsapply(1:M,function(nouse){
      x=sort(rgamma(n[j],shape=b+1,rate=1))
      MME=MM(x,60,0.03,0.3)
      cri=quantile(pbsapply(1:50,function(nouse){
        return(MM(sort(rexp(n[j],MME$theta)),60,0.03,0.3)$llf)
      }),0.05)
      if(MME$llf<cri) return(1)
      return(0)})
    betas[i,j]=sum(temp)/M
  }
}
save(betas,file="MM_power.Rda")


#######################
### Power of AUMPUT ###
#######################
h=function(alpha=0.05,n=10,theta=1){#simulate the critical value at
  #significance level alpha, sample size n and parameter theta
  N=100000
  N_sampels=sapply(1:N,function(no){
    x=rexp(n,theta)
    return(mean(log(x))-log(mean(x)))
  })
  return(quantile(N_sampels,1-alpha))
}
k=function(beta,alpha=0.05,n=10,theta=1){#simulate the power at beta, with
  #significance level alpha, sample size n and parameter theta
  N=100000
  N_sampels=sapply(1:N,function(no){
    x=rgamma(n,shape=beta+1,scale=1/theta)
    return(mean(log(x))-log(mean(x)))
  })
  c=h(alpha,n)
  return(sum(N_sampels>c)/N)
}

betas2=matrix(nrow=length(bb),ncol=length(n))
for(i in 1:length(bb)){#change beta at b
  for(j in 1:length(n)){#change n
    cat('b =',bb[i],'n =',n[j],'\n')
    betas2[i,j]=k(beta=bb[i],n=n[j])
  }
}
save(betas2,file="AUMPUT_power.Rda")


######################
### Power Curves 1 ###
######################
load('MM_power.Rda')
load('AUMPUT_power.Rda')
# par(bg = "white")
windowsFonts(A = windowsFont("Times New Roman"))
plot(1,type="n",ylab='power',xlab='b',xlim=c(0, 5),ylim=c(0,1),family="A",
     main='v(x)=x^b, N=60, alpha=0.03, beta=0.3')
# rect(par("usr")[1], par("usr")[3],
#      par("usr")[2], par("usr")[4],
#      col = "aliceblue")
for(i in 1:ncol(betas)){
  lines(bb,betas[,i],pch=18,col=heat.colors(5)[i],type="b",lty=1)
  lines(bb,betas2[,i],pch=19,col=heat.colors(5)[i],type="b",lty=2)
}
abline(h=0.05,lty=2)
legend(2.5,0.5,legend=c(as.vector(outer(c("  MM  ","AUMPUT"),n,paste,sep="  n="))),
       col=rep(heat.colors(5)[1:3],each=2),lty=1:2,pch=18:19,text.font=10,bg="#f7f7f7")



######################
### Power Curves 2 ###
######################
# wl(x) = 0 or 1 if 0 < x < b or b < x < 1,
# w2(x)= 0.1 or 1 if 0 < x < b or b < x < 1,
# w3(x)=0 or 1/2 or 1 if 0 < x < b or b < x < 1/2 or 1/2 < x < 1

# ?gamlss.family
library(gamlss)
library(gamlss.tr)
test1<-trun.r(par=c(1), family="EXP", type="left")
rr<-test1(100000,mu=1)
hist(rr,100)
