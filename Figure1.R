#S0 and I0 are the initial values of susceptibles and infectives
S0<-200
I0<-1
#gamma and beta are the removal and the pairwise infection paramteres
beta<-1
gamma<-110
#inttime is the time step for the discretization of the process
inttime<-0.0001
#Simulate the classical stochastic SIR model
M<-SIRsimexp(S0,I0,beta,gamma,inttime)
#Determine the absorption time and the discretized time vector
num<-ncol(M)
Time<-c(0:(num-1))
Time<-Time*inttime
#Extract the process I^1(t) from the matrix
N<-vector(length=num)
for (i in c(1:num)){
  N[i]<-M[2,i]
}
#Plot the sample path of the process
plot(Time,N,type="s",xlab="t",ylab="",main=expression(paste(alpha,"=1")),col="red")
#nu is the order of stability
nu<-0.7
#Simulate a stable subordinator
gamma2<-(cos(pi*nu/2))^(1/nu)
X<-vector(length=num)
X[1]<-0
i<-2
while(i<=num){
  Z<-rstable(1,nu,1,gamma2,0,pm=1)
  Z<-(inttime^(1/nu))*Z
  X[i]<-X[i-1]+Z
  i<-i+1
}
#Plot the sample path of the time-changed process
plot(X,N,type="s",xlab=expression(paste(sigma,"(t)")),ylab="",main=expression(paste(alpha,"=0.7")),col="blue")
