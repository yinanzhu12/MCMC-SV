library(MCMCpack)
library(truncnorm)

ln_p=function(y,h,sigma2,mu){
  return=-1.5*log(h)-y^2/(2*h)-(log(h)-mu)^2/(2*sigma2)
}

updateh=function(y,h,lhleft,lhright,alpha,delta,sigma_nu2){
  sigma2=sigma_nu2/(1+delta^2)
  mu=(delta*(lhleft+lhright)+(1-delta)*alpha)/(1+delta^2)
  s=(1-2*exp(sigma2))/(1-exp(sigma2))+0.5
  r=(s-1)*exp(mu+sigma2/2)+y^2/2
  qmode=r/(1+s)
  ln_c=log(1.1)+ln_p(y,qmode,sigma2,mu)-log(dinvgamma(qmode,scale=r,shape=s))
  hn=rinvgamma(1,shape=s,scale=r)
  accept1=exp(ln_p(y,hn,sigma2,mu)-ln_c)/dinvgamma(hn,scale=r,shape=s)
  roll=runif(1)
  tr=1
  rej=0
  rep=0
  while(roll>accept1){
    hn=rinvgamma(1,shape=s,scale=r)
    accept1=exp(ln_p(y,hn,sigma2,mu)-ln_c)/dinvgamma(hn,scale=r,shape=s)
    roll=runif(1)
    tr=tr+1
    rej=rej+1
  }
  if(accept1>=1){
  accept2=exp(ln_p(y,hn,sigma2,mu)-ln_p(y,h,sigma2,mu)+log(dinvgamma(h,scale=r,shape=s))-log(dinvgamma(hn,scale=r,shape=s)))
  roll=runif(1)
  if(roll>accept2){
      rep=1
      hn=h
  }
  }
  return=list(newh=hn,rep=rep,reject=rej,trial=tr)
}

updateh_rw=function(y,h,lhleft,lhright,alpha,delta,sigma_nu2,sd_rw){
  sigma2=sigma_nu2/(1+delta^2)
  mu=(delta*(lhleft+lhright)+(1-delta)*alpha)/(1+delta^2)
  ln_hn=rnorm(1,mean=log(h),sd=sd_rw)
  hn=exp(ln_hn)
  accept=exp(ln_p(y,hn,sigma2,mu)-ln_p(y,h,sigma2,mu))*hn/h
  roll=runif(1)
  rep=0
  if(roll>accept){
    hn=h
    rep=1
  }
  return=list(newh=hn,rep=rep,reject=0,trial=1)
}

updateh_rej=function(y,lhleft,lhright,alpha,delta,sigma_nu2){
  sigma2=sigma_nu2/(1+delta^2)
  mu=(delta*(lhleft+lhright)+(1-delta)*alpha)/(1+delta^2)
  mup=mu+sigma2/2*(y^2*exp(-mu)-1)
  ln_hn=rnorm(1,mean=mup,sd=sqrt(sigma2))
  accept=y^2*(exp(-mu)*(1+mu-ln_hn)-exp(-ln_hn))/2
  accept=exp(accept)
  roll=runif(1)
  rej=0
  tr=1
  while(roll>accept){
    ln_hn=rnorm(1,mean=mup,sd=sqrt(sigma2))
    accept=y^2*(exp(-mu)*(1+mu-ln_hn)-exp(-ln_hn))/2
    accept=exp(accept)
    roll=runif(1)
    tr=tr+1
    rej=rej+1
  }
  return=list(newh=exp(ln_hn),rep=0,reject=rej,trial=tr)
}

sampler=function(y,nu_0,s_0,delta_0,sigma_delta2,alpha_0,sigma_alpha2,t,burnin,sd_rw){
  n=length(y)
  delta=1
  alpha=0
  sigma_nu2=0.1
  h0=sd(y)^2
  h=rep(h0,n)
  #h=y^2+0.0000001
  ln_h=log(h)
  h_sample=matrix(nrow=n,ncol=t)
  alpha_sample=rep(0,t)
  delta_sample=rep(0,t)
  sigma_nu2_sample=rep(0,t)
  trial=0
  rej=0
  upd=0
  rep=0
  for(ite in 1:(t+burnin)){
    if(ite%%500==0){
      print(ite)
    }
    for(i in 2:(n-1)){
      r=updateh(y[i],h[i],ln_h[i-1],ln_h[i+1],alpha,delta,sigma_nu2)
      #r=updateh_rw(y[i],h[i],ln_h[i-1],ln_h[i+1],alpha,delta,sigma_nu2,sd_rw)
      #r=updateh_rej(y[i],ln_h[i-1],ln_h[i+1],alpha,delta,sigma_nu2)
      h[i]=r$newh
      upd=upd+1
      rep=rep+r$rep
      rej=rej+r$reject
      trial=trial+r$trial
      ln_h[i]=log(h[i])
    }
    ln_h[1]=rnorm(1,mean=alpha+delta*ln_h[2],sd=sqrt(sigma_nu2))
    h[1]=exp(ln_h[1])
    ln_h[n]=rnorm(1,mean=alpha+delta*ln_h[n-1],sd=sqrt(sigma_nu2))
    h[n]=exp(ln_h[n])
    s1=sum(ln_h)
    s2=sum(ln_h^2)
    s3=sum(ln_h[1:(n-1)]*ln_h[2:n])
    
    sigma_nu_scale=(s_0+(n-1)*alpha^2+(1+delta^2)*s2-ln_h[1]^2-delta^2*ln_h[n]^2-2*alpha*((1-delta)*s1-ln_h[1]+delta*ln_h[n])-2*delta*s3)/2
    sigma_nu2=rinvgamma(1,shape=(nu_0+n-1)/2,scale=sigma_nu_scale)
    
    deltamean=(sigma_nu2*delta_0+sigma_delta2*(s3-alpha*(s1-ln_h[n])))/(sigma_nu2+sigma_delta2*(s2-ln_h[n]^2))
    deltasd=sqrt(sigma_nu2*sigma_delta2/(sigma_nu2+sigma_delta2*(s2-ln_h[n]^2)))
    delta=rtruncnorm(1,a=-1,b=1,mean=deltamean,sd=deltasd)
    
    alphamean=(sigma_alpha2*((1-delta)*s1-ln_h[1]+delta*ln_h[n])+sigma_nu2*alpha_0)/(sigma_nu2+(n-1)*sigma_alpha2)
    alphasd=sqrt(sigma_nu2*sigma_alpha2/(sigma_nu2+(n-1)*sigma_alpha2))
    alpha=rnorm(1,mean=alphamean,sd=alphasd)
    #print(delta)
    
    if(ite>burnin){
    alpha_sample[ite-burnin]=alpha
    delta_sample[ite-burnin]=delta
    sigma_nu2_sample[ite-burnin]=sigma_nu2
    h_sample[,ite-burnin]=h
  }
  }
  return=list(h_sample=h_sample,alpha_sample=alpha_sample,delta_sample=delta_sample,sigma_nu2_sample=sigma_nu2_sample,rejectionrate=rej/trial,repeatrate=rep/upd)
}

simulate=function(delta,sigma_nu2,n){
  ln_h=rep(0,n)
  y=rep(0,n)
  for(i in 2:n)ln_h[i]=rnorm(1,mean=delta*ln_h[i-1],sd=sqrt(sigma_nu2))
  for(i in 1:n)y[i]=rnorm(1,sd=exp(ln_h[i]/2))
  return=list(y=y,ln_h=ln_h)
}