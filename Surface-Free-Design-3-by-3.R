########################################
###### This code can be used to implement the design proposed in
##### `A Surface-Free Design for Phase~I Dual-Agent Combination Trials' by
##### Mozgunov, Gasparini & Jaki (2020)
##### in the case of 3 doses of each agent (Section 3)
########################################



############## Specifying the Model ##############
########################################
library("rjags")
model1.string <-"
model {
s[1] ~ dbin(p[1], n[1])  
p[1] <- 1-theta[1]
s[2] ~ dbin(p[2], n[2])  
p[2] <- 1-theta[1]*theta[2]
s[3] ~ dbin(p[3], n[3])  
p[3] <- 1-theta[1]*theta[2]*theta[3]
s[4] ~ dbin(p[4], n[4])  
p[4] <- 1-theta[1]*theta[4]
s[5] ~ dbin(p[5], n[5])  
p[5] <- 1-theta[1]*theta[2]*theta[4]
s[6] ~ dbin(p[6], n[6])  
p[6] <- 1-theta[1]*theta[2]*theta[3]*theta[4]
s[7] ~ dbin(p[7], n[7])  
p[7] <- 1-theta[1]*theta[4]*theta[5]
s[8] ~ dbin(p[8], n[8])  
p[8] <- 1-theta[1]*theta[2]*theta[4]*theta[5]
s[9] ~ dbin(p[9], n[9])  
p[9] <- 1-theta[1]*theta[2]*theta[3]*theta[4]*theta[5]

theta[1] ~ dbeta(a[1],b[1])T(0,0.999999)
theta[2] ~ dbeta(a[2],b[2])T(0,0.999999)
theta[3] ~ dbeta(a[3],b[3])T(0,0.999999)
theta[4] ~ dbeta(a[4],b[4])T(0,0.999999)
theta[5] ~ dbeta(a[5],b[5])T(0,0.999999)
}
"
model1.spec<-textConnection(model1.string) # the truncation is needed for the computational purposes

# Function to Compute Mean Prior Point Estimate of Connections Using the Monotherapy Data
compute.prior.means.SFD<-function(p1,p2){
  t.prior<-mat.or.vec(length(p1)+length(p2)-1,1)
  t.prior[1]<-1-p1[1]-p2[1]+p1[1]*p2[1]
  for (i in 2:(length(p1))){
    t.prior[i]<-(1-p1[i])/(1-p1[i-1])
  }
  for (i in 1:(length(p2)-1)){
    t.prior[length(p2)+i]<-(1-p2[i+1])/(1-p2[i])
  }
  return(t.prior)
}

# Computing the Mean Prior Point Estimates
pA<-c(0.05,0.10,0.20)
pB<-c(0.10,0.20,0.30)
t.prior<-compute.prior.means.SFD(p1=B,p2=pA)


# Defining the Strength of Prior
c.prior<-rep(4,5)   # stregth of Beta prior distributions on the connections
a.prior<-t.prior*c.prior # Finding the first parameters of Beta distribution
b.prior<-(1-t.prior)*c.prior # Finding the second parameters of Beta distribution


target<-0.30         # Target Toxicity Level
cohort<-3            # Cohort Size
start.dose<-1        # Starting Dose (corresponds to the combination (1,1))
n<-12                # Number of cohorts in the trial
nsims<-2000          # Number of simulations  
iterations<-5000     # Number of samples in MCMC
safety<-T            # Apply Safety Constraint?


# Specifying the true combiantion-toxicity scenario
true<-c(0.020,0.050,0.120,
        0.100,0.200,0.300,
        0.150,0.300,0.500)


# Creating matrices to store the results

result.original<-mat.or.vec(nsims,length(true))
result.experiment<-mat.or.vec(nsims,length(true))
result.toxicity<-mat.or.vec(nsims,1)
experiment<-array(0,dim=c(n,9,nsims))
counter<-0


### Running the Design

set.seed(100) 
for (z in 1:nsims){
  datad<-mat.or.vec(1,length(true))
  doses.exp<-mat.or.vec(1,length(true))
  doses.tox<-mat.or.vec(1,length(true))
  
  outcome<-sum(rbinom(cohort,1,true[start.dose]))
  current.dose<-start.dose
  doses.exp[current.dose]<-doses.exp[current.dose]+cohort
  doses.tox[current.dose]<-doses.tox[current.dose]+outcome
  datan<-c(doses.exp)
  experiment[1,current.dose,z]<-cohort
  datas<-c(doses.tox)
  nextdose<-start.dose
  model1.spec<-textConnection(model1.string)
  mydata <- list(s = datas,n = datan, a=a.prior,b=b.prior)
  jags <- jags.model(model1.spec,data =mydata,n.chains=2,n.adapt=iterations,quiet=TRUE)
  update(jags, iterations,progress.bar="none")
  tt<-jags.samples(jags,c('theta'),iterations,progress.bar="none")
  t<-cbind(c(tt$theta[1,,]),c(tt$theta[2,,]),c(tt$theta[3,,]),c(tt$theta[4,,]),c(tt$theta[5,,]))
  t.mean<-colMeans(t)
  p1<-1-t.mean[1]
  p2<-1-t.mean[1]*t.mean[2]
  p3<-1-t.mean[1]*t.mean[2]*t.mean[3]
  p4<-1-t.mean[1]*t.mean[4]
  p5 <- 1-t.mean[1]*t.mean[2]*t.mean[4]
  p6 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[4]
  p7 <- 1-t.mean[1]*t.mean[4]*t.mean[5]
  p8 <- 1-t.mean[1]*t.mean[2]*t.mean[4]*t.mean[5]
  p9 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[4]*t.mean[5]
  if(nextdose==1){
    p3<-p5<-p6<-p7<-p8<-p9<-1}
  if(nextdose==2){
    p6<-p7<-p8<-p9<-1}
  if(nextdose==3){
    p7<-p8<-p9<-1
  }
  if(nextdose==4){
    p6<-p8<-p9<-1}
  if(nextdose==5){
    p9<-1}
  if(nextdose==7){
    p9<-1}
  p<-c(p1,p2,p3,p4,p5,p6,p7,p8,p9)
  nextdose<-which(abs(p-target)==min(abs(p-target)))
  
  for (i in 2:n){
    outcome<-sum(rbinom(cohort,1,true[nextdose]))
    current.dose<-nextdose
    doses.exp[current.dose]<-doses.exp[current.dose]+cohort
    doses.tox[current.dose]<-doses.tox[current.dose]+outcome
    datan<-c(doses.exp)
    experiment[i,current.dose,z]<-cohort
    datas<-c(doses.tox)
    model1.spec<-textConnection(model1.string)
    mydata <- list(s = datas,n = datan, a=a.prior,b=b.prior)
    jags <- jags.model(model1.spec,data =mydata,n.chains=2,n.adapt=iterations,quiet=TRUE)
    update(jags, iterations,progress.bar="none")
    tt<-jags.samples(jags,c('theta'),iterations,progress.bar="none")
    t<-cbind(c(tt$theta[1,,]),c(tt$theta[2,,]),c(tt$theta[3,,]),c(tt$theta[4,,]),c(tt$theta[5,,]))
    t.mean<-colMeans(t)
    p1<-1-t.mean[1]
    p2<-1-t.mean[1]*t.mean[2]
    p3<-1-t.mean[1]*t.mean[2]*t.mean[3]
    p4<-1-t.mean[1]*t.mean[4]
    p5 <- 1-t.mean[1]*t.mean[2]*t.mean[4]
    p6 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[4]
    p7 <- 1-t.mean[1]*t.mean[4]*t.mean[5]
    p8 <- 1-t.mean[1]*t.mean[2]*t.mean[4]*t.mean[5]
    p9 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[4]*t.mean[5]
    if(nextdose==1){
      p3<-p5<-p6<-p7<-p8<-p9<-1}
    if(nextdose==2){
      p6<-p7<-p8<-p9<-1}
    if(nextdose==3){
      p7<-p8<-p9<-1
    }
    if(nextdose==4){
      p6<-p8<-p9<-1}
    if(nextdose==5){
      p9<-1}
    if(nextdose==7){
      p9<-1}
    
    if(safety){
      if(length(which(1-t[,1]>target))/length(t[,1])>0.7){
        p<-rep(1,9)
      }else{
        if(length(which(1-t[,1]*t[,2]>target))/length(t[,1])>0.7){
          p2<-1
        }
        if(length(which(1-t[,1]*t[,2]*t[,3]>target))/length(t[,1])>0.7){
          p3<-1
        }    
        if(length(which(1-t[,1]*t[,4]>target))/length(t[,1])>0.7){
          p4<-1
        } 
        if(length(which(1-t[,1]*t[,2]*t[,4]>target))/length(t[,1])>0.7){
          p5<-1
        } 
        if(length(which(1-t[,1]*t[,2]*t[,3]*t[,4]>target))/length(t[,1])>0.7){
          p6<-1
        } 
        if(length(which(1-t[,1]*t[,4]*t[,5]>target))/length(t[,1])>0.7){
          p7<-1
        }
        if(length(which(1-t[,1]*t[,2]*t[,4]*t[,5]>target))/length(t[,1])>0.7){
          p8<-1
        }
        if(length(which(1-t[,1]*t[,2]*t[,3]*t[,4]*t[,5]>target))/length(t[,1])>0.7){
          p9<-1
        }
      }
    }
    p<-c(p1,p2,p3,p4,p5,p6,p7,p8,p9)
    

      nextdose<-which(abs(p-target)==min(abs(p-target)))
    
    if(safety){
      if(length(which(1-t[,1]>target))/length(t[,1])>0.7){
        break
      }
    }
  }
  if(safety){
    if(length(which(1-t[,1]>0.2))/length(t[,1])>0.7){
      result.original[z,nextdose]<-0
      cat("termination","\n")
      counter<-counter+1
    }else{
  result.original[z,nextdose]<-1 
    }}
  cat("recommendation is",nextdose,"\n")
  result.original[z,nextdose]<-1 
  result.experiment[z,]<-datan
  result.toxicity[z]<-sum(datas)
  cat(z,"out of",nsims,"\n")
  cat("The current proportion of correct selections is",(colSums(result.original)[6]+colSums(result.original)[8])/z,"\n")
}


# Storing Proportion of Selections
y<-selections<-colMeans(result.original)
selections[6] # Selection of the first MTC
selections[8] # Selection of the second MTC
selections[6]+selections[8] # The proportion of correct selections


# Computing the allocation probabilities
mozg<-experiment[,,1]/cohort
for (w in 2:nsims){
  mozg<-mozg+experiment[,,w]/cohort
}
prop<-mozg/nsims

y1<-prop[,1]+prop[,2]+prop[,4]  # Allocation probabilities for the low combinations
y2<-prop[,3]+prop[,5]+prop[,7]  # Allocation probabilities for the medium combinations
y3<-prop[,6]+prop[,8]           # Allocation probabilities for the MTCs
y4<-prop[,9]                    # Allocation probabilities for the high combinations



# Plotting the allocation probabilities figure
par(mfrow=c(1,1),mai=c(0.8,0.8,0.3,0.2))
plot(1:13,c(y1,sum(c(selections[1],selections[2],selections[4]))),type="l",ylim=c(0,1),col="green",main="Allocation probability",xlab="cohort",ylab="Probability")
lines(1:13,c(y2,sum(c(selections[3],selections[5],selections[7]))),col="black")
lines(1:13,c(y3,sum(c(selections[6],selections[8]))),col="blue")
lines(1:13,c(y4,selections[9]),col="red")
abline(h=0.6784,lty=2) # Performance of the Benchmark
legend(8, 0.95, c("Low combinations" ,"Medium combinations","MTCs","High combinations"),pch=15,ncol = 1,col=c("green","black","blue","red"),cex=1.3)

