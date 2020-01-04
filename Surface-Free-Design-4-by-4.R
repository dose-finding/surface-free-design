########################################
###### This code can be used to implement the design proposed in
##### `A Surface-Free Design for Phase~I Dual-Agent Combination Trials' by
##### Mozgunov, Gasparini & Jaki (2020)
##### in the case of 4 doses of each agent (Section 4 and Supplementary Materials)
########################################

library("rjags") # Uploading JAGS library

# Defining the model with respect to the parametrisation (6)
model1.string <-"
model {
s[1] ~ dbin(p[1], n[1])  
p[1] <- 1-theta[1]
s[2] ~ dbin(p[2], n[2])  
p[2] <- 1-theta[1]*theta[2]
s[3] ~ dbin(p[3], n[3])  
p[3] <- 1-theta[1]*theta[2]*theta[3]
s[4] ~ dbin(p[4], n[4])  
p[4] <- 1-theta[1]*theta[2]*theta[3]*theta[4]
s[5] ~ dbin(p[5], n[5])  
p[5] <- 1-theta[1]*theta[5]
s[6] ~ dbin(p[6], n[6])  
p[6] <- 1-theta[1]*theta[2]*theta[5]
s[7] ~ dbin(p[7], n[7])  
p[7] <- 1-theta[1]*theta[2]*theta[3]*theta[5]
s[8] ~ dbin(p[8], n[8])  
p[8] <- 1-theta[1]*theta[2]*theta[3]*theta[4]*theta[5]
s[9] ~ dbin(p[9], n[9])  
p[9] <- 1-theta[1]*theta[5]*theta[6]
s[10] ~ dbin(p[10], n[10])  
p[10] <- 1-theta[1]*theta[2]*theta[5]*theta[6]
s[11] ~ dbin(p[11], n[11])  
p[11] <- 1-theta[1]*theta[2]*theta[3]*theta[5]*theta[6]
s[12] ~ dbin(p[12], n[12])  
p[12] <- 1-theta[1]*theta[2]*theta[3]*theta[4]*theta[5]*theta[6]
s[13] ~ dbin(p[13], n[13])  
p[13] <- 1-theta[1]*theta[5]*theta[6]*theta[7]
s[14] ~ dbin(p[14], n[14])  
p[14] <- 1-theta[1]*theta[2]*theta[5]*theta[6]*theta[7]
s[15] ~ dbin(p[15], n[15])  
p[15] <- 1-theta[1]*theta[2]*theta[3]*theta[5]*theta[6]*theta[7]
s[16] ~ dbin(p[16], n[16])  
p[16] <- 1-theta[1]*theta[2]*theta[3]*theta[4]*theta[5]*theta[6]*theta[7]

theta[1] ~ dbeta(a[1],b[1])T(0,0.999999)
theta[2] ~ dbeta(a[2],b[2])T(0,0.999999)
theta[3] ~ dbeta(a[3],b[3])T(0,0.999999)
theta[4] ~ dbeta(a[4],b[4])T(0,0.999999)
theta[5] ~ dbeta(a[5],b[5])T(0,0.999999)
theta[6] ~ dbeta(a[6],b[6])T(0,0.999999)
theta[7] ~ dbeta(a[7],b[7])T(0,0.999999)
}
"
model1.spec<-textConnection(model1.string) # the truncation above is needed for the computational purposes
# END of the Model Specification



### Defining the parameters of the Beta prior distribution for the connections
a.prior<-rep(3.81,7)  # prior values of a for beta distributions
b.prior<-rep(0.19,7)  # prior values of b for beta distributions

### Parameters of the Trial
target<-0.20    # Target Toxicity Level
cohort<-1       # Cohort Size
start.dose<-1   # Starting dose (corresponds to the lowest combination)
n<-50           # Number of Cohorts in the Trial
nsims<-2000     # Number of Simulations
iterations<-3000# Number of samples used by the MCMC
safety<-T       # Apply Safety Constraint (Default=TRUE)
randomize<-F    # Apply Randomisation Rule based on the point estimates? (Default=F) Note: The randomisation is discussed in Supplementary Materials

### Define the scenarios
### The scenarios are as in Supplementrary Materials 
### The scenarios number from the main body of the paper is in brackets

# true<-c(0.04,0.08,0.12,0.16,0.10,0.14,0.18,0.22,0.16,0.20,0.24,0.28,0.22,0.26,0.30,0.34)          #Scenario A
# true<-c(0.02,0.04,0.06,0.08,0.05,0.07,0.09,0.11,0.08,0.10,0.12,0.14,0.11,0.13,0.15,0.17)          #Scenario B (#1)
# true<-c(0.10,0.20,0.30,0.40,0.25,0.35,0.45,0.55,0.40,0.50,0.60,0.70,0.55,0.65,0.75,0.85)          #Scenario C (#2)
# true<-c(0.44,0.48,0.52,0.56,0.50,0.54,0.58,0.52,0.56,0.60,0.64,0.68,0.62,0.66,0.70,0.74)          #Scenario D (#3)
# true<-c(0.08,0.18,0.28,0.29,0.09,0.19,0.29,0.30,0.10,0.20,0.30,0.31,0.11,0.21,0.31,0.41)          #Scenario E
# true<-c(0.12,0.13,0.14,0.15,0.16,0.18,0.20,0.22,0.44,0.45,0.46,0.47,0.50,0.52,0.54,0.55)          #Scenario F (#4)
# true<-c(0.01,0.02,0.03,0.04,0.04,0.10,0.15,0.20,0.06,0.15,0.30,0.45,0.10,0.30,0.50,0.80)          #Scenario G (#5)
true<-c(0.06,0.08,0.10,0.15,0.10,0.12,0.30,0.45,0.15,0.30,0.50,0.60,0.50,0.55,0.60,0.70)/3*2      #Scenario H (#6)
# true<-c(0.08,0.12,0.16,0.18,0.10,0.15,0.30,0.45,0.12,0.30,0.50,0.55,0.30,0.50,0.55,0.60)/3*2      #Scenario I (#7)
# true<-c(0.05,0.10,0.15,0.30,0.10,0.15,0.30,0.45,0.15,0.30,0.45,0.50,0.30,0.45,0.50,0.60)/3*2      #Scenario J
# true<-c(0.15,0.30,0.45,0.50,0.30,0.45,0.50,0.60,0.45,0.55,0.60,0.70,0.55,0.60,0.70,0.80)/3*2      #Scenario K
# true<-c(0.02,0.07,0.10,0.15,0.07,0.10,0.15,0.30,0.10,0.15,0.30,0.45,0.15,0.30,0.45,0.55)/3*2      #Scenario L
# true<-c(0.30,0.45,0.60,0.70,0.45,0.55,0.65,0.75,0.50,0.60,0.70,0.80,0.60,0.70,0.80,0.90)/3*2      #Scenario M
# true<-c(0.01,0.02,0.08,0.10,0.03,0.05,0.10,0.13,0.07,0.09,0.12,0.15,0.09,0.12,0.15,0.30)/3*2      #Scenario N
# true<-c(0.05,0.08,0.10,0.13,0.09,0.12,0.15,0.30,0.15,0.30,0.45,0.50,0.45,0.50,0.60,0.70)/3*2      #Scenario O (#8)
# true<-c(0.07,0.10,0.15,0.30,0.15,0.30,0.45,0.52,0.30,0.50,0.60,0.65,0.45,0.60,0.70,0.75)/3*2      #Scenario P (#9)
# true<-c(0.02,0.10,0.15,0.50,0.05,0.12,0.30,0.55,0.08,0.15,0.45,0.60,0.15,0.45,0.60,0.80)/3*2      #Scenario Q
# true<-c(0.10,0.12,0.30,0.40,0.15,0.30,0.37,0.43,0.30,0.37,0.42,0.47,0.37,0.42,0.47,0.52)/3*2      #Scenario R
# true<-c(0.01,0.03,0.06,0.08,0.04,0.07,0.12,0.16,0.08,0.10,0.15,0.30,0.10,0.15,0.30,0.50)/3*2      #Scenario S
# true<-c(0.06,0.10,0.15,0.30,0.10,0.30,0.50,0.70,0.50,0.60,0.70,0.80,0.60,0.70,0.80,0.90)/3*2      #Scenario T
# true<-c(0.05,0.12,0.20,0.30,0.10,0.20,0.30,0.40,0.30,0.42,0.52,0.62,0.45,0.50,0.60,0.70)/3*2      #Scenario U (#10)
# true<-c(0.12,0.20,0.30,0.40,0.20,0.30,0.40,0.50,0.42,0.52,0.62,0.70,0.52,0.62,0.70,0.80)/3*2      #Scenario V
# true<-c(0.04,0.06,0.08,0.20,0.10,0.20,0.30,0.50,0.30,0.42,0.52,0.70,0.42,0.52,0.70,0.80)/3*2      #Scenario W



# Defining matrices to store the results
result.original<-mat.or.vec(nsims,length(true))
result.experiment<-mat.or.vec(nsims,length(true))
result.toxicity<-mat.or.vec(nsims,1)


# Starting the Simulations
counter<-0
for (z in 1:nsims){
  datad<-mat.or.vec(1,length(true))
  doses.exp<-mat.or.vec(1,length(true))
  doses.tox<-mat.or.vec(1,length(true))
  
  outcome<-sum(rbinom(cohort,1,true[start.dose]))
  current.dose<-start.dose
  doses.exp[current.dose]<-doses.exp[current.dose]+cohort
  doses.tox[current.dose]<-doses.tox[current.dose]+outcome
  datan<-c(doses.exp)
  datas<-c(doses.tox)
  
  model1.spec<-textConnection(model1.string)
  mydata <- list(s = datas,n = datan, a=a.prior,b=b.prior)
  jags <- jags.model(model1.spec,data =mydata,n.chains=1,n.adapt=iterations,quiet=TRUE)
  update(jags, iterations,progress.bar="none")
  tt<-jags.samples(jags,c('theta'),iterations,progress.bar="none")
  t<-cbind(c(tt$theta[1,,]),c(tt$theta[2,,]),c(tt$theta[3,,]),c(tt$theta[4,,]),c(tt$theta[5,,]),c(tt$theta[6,,]),c(tt$theta[7,,]))
  t.mean<-colMeans(t)
  p1<-1-t.mean[1]
  p2<-1-t.mean[1]*t.mean[2]
  p3<-1-t.mean[1]*t.mean[2]*t.mean[3]
  p4<-1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[4]
  p5 <- 1-t.mean[1]*t.mean[5]
  p6 <- 1-t.mean[1]*t.mean[2]*t.mean[5]
  p7 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[5]
  p8 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[4]*t.mean[5]
  p9 <- 1-t.mean[1]*t.mean[5]*t.mean[6]
  p10 <- 1-t.mean[1]*t.mean[2]*t.mean[5]*t.mean[6]
  p11 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[5]*t.mean[6]
  p12 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[4]*t.mean[5]*t.mean[6]
  p13 <- 1-t.mean[1]*t.mean[5]*t.mean[6]*t.mean[7]
  p14 <- 1-t.mean[1]*t.mean[2]*t.mean[5]*t.mean[6]*t.mean[7]
  p15 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[5]*t.mean[6]*t.mean[7]
  p16 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[4]*t.mean[5]*t.mean[6]*t.mean[7]
  p<-c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16)
  
  
  
  
  if(randomize){
    A<-which(abs(p-target)<=0.05)
    if(length(A)>1){
      nextdose<-sample(A, 1, prob = rep(1/length(A),length(A)),replace = TRUE)
    }
    else{
      nextdose<-which(abs(p-target)==min(abs(p-target)))
    }
  }
  else{
    nextdose<-which(abs(p-target)==min(abs(p-target)))
  }  
  
  
  for (i in 2:n){
    outcome<-sum(rbinom(cohort,1,true[nextdose]))
    current.dose<-nextdose
    doses.exp[current.dose]<-doses.exp[current.dose]+cohort
    doses.tox[current.dose]<-doses.tox[current.dose]+outcome
    datan<-c(doses.exp)
    datas<-c(doses.tox)
    model1.spec<-textConnection(model1.string)
    mydata <- list(s = datas,n = datan, a=a.prior,b=b.prior)
    jags <- jags.model(model1.spec,data =mydata,n.chains=1,n.adapt=iterations,quiet=TRUE)
    update(jags, iterations,progress.bar="none")
    tt<-jags.samples(jags,c('theta'),iterations,progress.bar="none")
    t<-cbind(c(tt$theta[1,,]),c(tt$theta[2,,]),c(tt$theta[3,,]),c(tt$theta[4,,]),c(tt$theta[5,,]),c(tt$theta[6,,]),c(tt$theta[7,,]))
    t.mean<-colMeans(t)
    p1<-1-t.mean[1]
    p2<-1-t.mean[1]*t.mean[2]
    p3<-1-t.mean[1]*t.mean[2]*t.mean[3]
    p4<-1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[4]
    p5 <- 1-t.mean[1]*t.mean[5]
    p6 <- 1-t.mean[1]*t.mean[2]*t.mean[5]
    p7 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[5]
    p8 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[4]*t.mean[5]
    p9 <- 1-t.mean[1]*t.mean[5]*t.mean[6]
    p10 <- 1-t.mean[1]*t.mean[2]*t.mean[5]*t.mean[6]
    p11 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[5]*t.mean[6]
    p12 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[4]*t.mean[5]*t.mean[6]
    p13 <- 1-t.mean[1]*t.mean[5]*t.mean[6]*t.mean[7]
    p14 <- 1-t.mean[1]*t.mean[2]*t.mean[5]*t.mean[6]*t.mean[7]
    p15 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[5]*t.mean[6]*t.mean[7]
    p16 <- 1-t.mean[1]*t.mean[2]*t.mean[3]*t.mean[4]*t.mean[5]*t.mean[6]*t.mean[7]
    p<-c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16)
    
    if(randomize){
      if(i==n){
        nextdose<-which(abs(p-target)<=0.05)
        if(length(nextdose)==0){
          nextdose<-which(abs(p-target)==min(abs(p-target)))
        }
      }else{
        A<-which(abs(p-target)<=0.05)
        if(length(A)>1){
          nextdose<-sample(A, 1, prob = rep(1/length(A),length(A)),replace = TRUE)
        }
        else{
          nextdose<-which(abs(p-target)==min(abs(p-target)))
        }
      }
    }
    else{
      nextdose<-which(abs(p-target)==min(abs(p-target)))
    }
    
    if(safety){
      if(length(which(1-t[,1]>0.2))/length(t[,1])>0.7){
        break
      }
    }
    # cat("Cohort",i,"\n")
  }
  if(safety){
    if(length(which(1-t[,1]>0.2))/length(t[,1])>0.7){
      result.original[z,nextdose]<-0
      cat("termination","\n")
      counter<-counter+1
    }else{
      result.original[z,nextdose]<-1 
    }}
  result.original[z,nextdose]<-1 
  result.experiment[z,]<-datan
  result.toxicity[z]<-sum(datas)
  cat(z,"out of",nsims,"\n")
}
y<-result.original

# END of Simulations

# Printing the Results

recommendations.scenario.H<-y  # Selection Proportions
terminations.scenario.H<-counter/nsims  # Proportion of Terminated Trials
experimentation.scenario.H<-colMeans(result.experiment)   # Experimentational Proportions
toxicities.scenario.H<-mean(result.toxicity)    # Average number of toxicity responses
mean.patients.scenario.H<-sum(colMeans(result.experiment))   #Average Number of Patients in the Trial

