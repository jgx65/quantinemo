# function
plot.fitness.stabilizing<-function(Zopt=0, omega=1, min=-10, max=10)
{
  phenotype<-sort(runif(200,min,max))
  plot(phenotype, exp((phenotype - Zopt)^2/((-2)*omega^2)),
      ylab="fitness", ylim=c(0,1), type="l", lwd=2)
      
  lines(c(0,0),c(0.07,1.02))  
  text(0 ,0.02 ,expression(paste("Z"[Opt])))
  
  arrows(0, exp((omega - Zopt)^2/((-2)*omega^2)), omega,exp((omega - Zopt)^2/((-2)*omega^2)), code=3, length=0.1)
  text(omega/2, exp((omega - Zopt)^2/((-2)*omega^2))-0.05,expression(omega))
}

# how to call the function
plot.fitness.stabilizing(Zopt=0, omega=2, min=-8, max=8)# function



plot.fitness.directional<-function(r=1, r_max=0, s=1, min=-10, max=10)
{
  phenotype<-sort(runif(200,min,max))
  plot(phenotype, ((1+s*exp(r*(r_max-phenotype)))^(-1/s)),
      ylab="fitness", ylim=c(0,1), type="l", lwd=2)
      
  lines(c(0,0),c(0.07,1.02))  
  text(0 ,0.02 ,expression(paste("P"[rMax])))
  
  lines(c(0,0.7),c(0.5,0.5)) 
  x<-sort(runif(100))/2+0.5
  lines(0.05+x/2, 0.491+sqrt((1-x*x))/10)
  text(0.7, 0.55,"r")
  
  arrows(-2, 0.65, -2, 0.15, length=0.1)
  arrows(-1.8, 0.7, 2, 0.9, length=0.1)
  text(-2,0.69,"s")
  
}

# how to call the function
plot.fitness.directional(r=1, r_max=0, s=1, min=-4, max=4)


plot.fitness.directional<-function(r=1, r_max=0, s=1, min=-10, max=10, range.min=-10, range.max=10)
{
  phenotype<-sort(runif(2000,range.min,range.max))
  plot(phenotype, min+((max-min)/(1+(s*exp(r*(r_max-phenotype)))^(1/s))),
      ylab="factor", xlab="population density (N/K)", ylim=c(0,1), type="l", lwd=2)
    
}

# how to call the function
plot.fitness.directional(r=1e4, r_max=0.5, s=1, min=0, max=1, range.min=0, range.max=1.1)



## fitness landscape
Zopt1<- -2
Zopt2<-  2
omega<- 0.5
nbPoints<-10
min=-4
max=4
  phenotype<-0:nbPoints/nbPoints*(max-min)+min
w<-(exp((phenotype - Zopt1)^2/((-2)*omega^2)) + exp((phenotype - Zopt2)^2/((-2)*omega^2)))
plot(phenotype, w,
      ylab="fitness", ylim=c(0,1), type="l", lwd=2)
