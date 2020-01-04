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


pA<-c(0.05,0.10,0.20)
pB<-c(0.10,0.20,0.30)
t.prior<-compute.prior.means.SFD(p1=pB,p2=pA)
t.prior
