#Function for the simulation of the classic stochastic SIR model:
#In input:
#S0 is the initial value of the suscpetible individuals and must be an integer;
#I0 is the initial value of the infective individuals and must be an integer;
#beta is the pairwise infection parameter;
#gamma is the removal parameter;
#inttime is the time interval for the discretization of the process.
#In output: the function gives a matrix M in output whose first row is the process S, the second row is I and the third row is R.
SIRsimexp<-function(S0,I0,beta,gamma,inttime){
  #Simulation of the first inter-jump time
  rate<-I0*(beta*S0+gamma)
  IntEv<-rexp(1, rate=rate)
  #Initialization of the Calendar vector (i. e. the vector that contains all the jump times)
  Calendar<-vector(length=1)
  Calendar[1]<-IntEv[1]
  i<-2
  #Initialization of the processes S, I and R
  S<-vector(length=1)
  I<-vector(length=1)
  R<-vector(length=1)
  S[1]<-S0
  I[1]<-I0
  R[1]<-0
  j<-1
  #Starting the simulation: the stop condition is I=0.
  while(I[i-1]!=0){
    S<-c(S,0)
    I<-c(I,0)
    R<-c(R,0)
    if((i-1)*inttime<Calendar[j]){
      #If a new event has not already occured, fix the values.
      S[i]<-S[i-1]
      I[i]<-I[i-1]
      R[i]<-R[i-1]
    }else{
      #If a new event has occured
      j<-j+1
      if(S[i-1]==0){
        #If there are no susceptibles, the only event that could happen is the decrease of infectives
        S[i]=0
        I[i]=I[i-1]-1
        R[i]=R[i-1]+1
      }else{
        #If there are susceptibles, decide what event has to happen
        U<-runif(1)
        rate<-I[i-1]*(beta*S[i-1]+gamma)
        prob<-gamma*I[i-1]/rate
        if(U<prob){
          S[i]<-S[i-1]
          I[i]<-I[i-1]-1
          R[i]<-R[i-1]+1
        }else{
          S[i]<-S[i-1]-1
          I[i]<-I[i-1]+1
          R[i]<-R[i-1]
        }
      }
      #If there are still infectives, simulate a new inter-jump time and update the calendar
      if(I[i]!=0){
        rate<-I[i]*(beta*S[i]+gamma)
        IntEv<-rexp(1,rate = rate)
        Calendar<-c(Calendar,0)
        Calendar[j]<-Calendar[j-1]+IntEv[1]
      }
    }
    i<-i+1
  }
  #Create the output matrix
  numcol<-i-1
  M<-matrix(nrow=3,ncol=numcol)
  for( j in c(1:numcol)){
    M[1,j]<-S[j]
    M[2,j]<-I[j]
    M[3,j]<-R[j]
  }
  return(M)
}
