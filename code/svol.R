library(stats)
source('~/GitHub/SDS383D-course-work/Project/functions.R')

x=0.01
sigma2=0.05
mu=-3
h=seq(from=1,to=10,by=0.01)
s=(1-2*exp(sigma2))/(1-exp(sigma2))+0.5
r=(s-1)*exp(mu+sigma2/2)+x^2/2
qmode=r/(1+s)
c=p(x,qmode,sigma2,mu)/dinvgamma(qmode,scale=r,shape=s)
plot(p(x,h,sigma2,mu)~h,type='l',col="black",ylab="")
lines(c*dinvgamma(h,shape=s,scale=r)~h,col="blue")
lines(1.1*c*dinvgamma(h,shape=s,scale=r)~h,col="red")
lines(1.2*c*dinvgamma(h,shape=s,scale=r)~h,col="green")
legend(x=6.8,y=0.25,legend=c("p","c=1","c=1.1","c=1.2"),fill=c("black","blue","red","green"))


sp <- read.csv("~/GitHub/SDS383D-course-work/Project/SP.csv")
n=length(sp$Close)
logchange=log(sp$Close[2:n])-log(sp$Close[1:(n-1)])
plot(logchange,type='l',xlab='',ylab='log of price')

proc.time()
sample=sampler(mockdata$y,1,1,0,10,0,10,8000,4000,0.05)
proc.time()
plot(sample$delta_sample,type='l',ylab="",main="delta")
hist(sample$delta_sample,xlab='',ylab="",breaks=20,main="histogram of delta")
plot(sample$alpha_sample,type='l',ylab="",main="alpha")
hist(sample$alpha_sample,xlab='',ylab="",breaks=20,main="histogram of alpha")
plot(sample$sigma_nu2_sample,type='l',ylab="",main="sigma_nu2")
hist(sample$sigma_nu2_sample,xlab='',ylab="",main="histogram of sigma_nu^2")
plot(log(sample$h_sample[100,]),type='l',ylab="log(h)",xlab="iteration",main="evolution of h")
sample$rejectionrate
logmean=log(rowMeans(sample$h_sample))
plot(logmean,type='l',ylab="log(h)",main="log of posterior mean of h")

dat=sample
a=dat$alpha_sample
d=dat$delta_sample
s=dat$sigma_nu2_sample
h=dat$h_sample
hmean=rowMeans(h)
sqrt(sum((hmean-exp(mockdata$ln_h))^2)/1000)

logmean=log(rowMeans(h))
plot(logmean,type='l',ylab="log(h)",main="log of posterior mean of h")

mean(a)
sd(a)
mean(d)
sd(d)
mean(s)
sd(s)

hist(d,xlab='',ylab="",breaks=20,main="histogram of delta")
hist(a,xlab='',ylab="",breaks=20,main="histogram of alpha")
hist(s,xlab='',ylab="",breaks=20,main="histogram of sigma_nu^2")
cov(a,d)
cov(a,s)
cov(s,d)
acf(a,lag.max=8000,main="correlagram of alpha")
acf(d,lag.max=8000,main="correlagram of delta")
acf(s,lag.max=8000,main="correlagram of sigma_nu^2")

mockdata=simulate(0.95,0.06,1000)
