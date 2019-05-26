#Ntraj=Number of trajectories
Ntraj=10000
#S0 and I0 are the initial values of susceptibles and infectives
S0=500
I0=1
#gamma and beta are the removal and the pairwise infection paramteres
gamma=150
beta=1
#inttime is the time step for the discretization of the process
inttime=0.01
#Initialize two vectors
Timeofout=vector(length=Ntraj)
Outvalue=vector(length=Ntraj)
#Simulate Ntraj classical SIR processes
for (i in c(1:Ntraj)){
  print(i)
  M=SIRsimexp(S0,I0,beta,gamma,inttime)
  D=dim(M)
  Timeofout[i]=(D[2]-1)*inttime
  Outvalue[i]=M[1,D[2]]
}
#Save the vectors
Timeofoutexp<-Timeofout
Outvalueexp<-Outvalue
#nu is the order of stability
nu=0.70
#Initilize the vectors
Timeofout=vector(length=Ntraj)
Outvalue=vector(length=Ntraj)
#Simulate Ntraj time-changed SIR processes
for (i in c(1:Ntraj)){
  print(i)
  M=SIRsimfrac(S0,I0,beta,gamma,nu,inttime)
  D=dim(M)
  Timeofout[i]=(D[2]-1)*inttime
  Outvalue[i]=M[1,D[2]]
}
#Save the vectors
TimeofoutTC<-Timeofout
OutvalueTC<-Outvalue

#Produce the histograms for the distribution of $S_\infty$ where $S_\infty$<50
hist(Outvalueexp[which(Outvalueexp<50)],xlim=c(0,50),breaks=50, xlab="", ylab="", main=expression(paste(alpha,"=1")),col="red",prob=TRUE)
hist(Outvalue3[which(Outvalue3<50)],xlim=c(0,50),breaks=50, xlab="", ylab="", main=expression(paste(alpha,"=0.7")),col="blue",prob=TRUE)

#Calculate the intensity of the epidemics
Iexp<-(S0-Outvalueexp)/S0
IvalueTC<-(S0-OutvalueTC)/S0

#Define the function Pi for the classical and the time-changed process
Piexp<-function(t){
  z<-length(which(Iexp<t))/10000
  return(z)
}
PiTC<-function(t){
  z<-length(which(IvalueTC<t))/10000
  return(z)
}
#Discretize the function Pi
Piexpvec<-vector(length=1001)
for (i in c(1:1001)){
  Piexpvec[i]<-Piexp(0.001*(i-1))
}
PiTCvec<-vector(length=1001)
for (i in c(1:1001)){
  PiTCvec[i]<-PiTC(0.001*(i-1))
}
#Discretized time vector
Time<-vector(length=1001)
for (i in c(1:1001)){
  Time[i]<-0.001*(i-1)
}
#Plot of the function Pi
plot(Time, Piexpvec, type="l",xlab="",ylab="",main=expression(paste(alpha,"=1")),col="red")
plot(Time, PiTCvec, type="l",xlab="",ylab="",main=expression(paste(alpha,"=0.7")),col="blue")
#Definition of the Tail functions
Tailexp<-function(t){
  z<-length(which(Timeofoutexp>t))/10000
  return(z)
}
TailTC<-function(t){
  z<-length(which(TimeofoutTC>t))/10000
  return(z)
}
#Discretization of the Tail funtions
Tailexpvec<-vector(length=1001)
for (i in c(1:1001)){
  Tailexpvec[i]<-Tailexp(0.1*(i-1))
}
TailTCvec<-vector(length=1001)
for (i in c(1:1001)){
  TailTCvec[i]<-TailTC(0.1*(i-1))
}
#Discretized time vector
Time<-vector(length=1001)
for (i in c(1:1001)){
  Time[i]<-0.01*(i-1)
}
#Plot of the Tail functions
plot(Time, Tailexpvec, type="l",xlim=c(0.9,1),xlab="",ylab="",col="red")
par(new=TRUE)
plot(Time, TailTCvec, type="l",xlim=c(0.9,1),xlab="",ylab="",col="blue")

