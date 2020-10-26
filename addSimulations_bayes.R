
library(bandit)
library(abind)
source('PlayBandit_Regret.R')

#full Bayesian simulation (p varies!):


firstRun=FALSE
if (!firstRun) {
  load("xMC_bayes.rda") 
  TS = gsub(" |:", "-", as.character(Sys.time()))
  save(x1MC, x2MC, p, file=paste0("xMC_bayes_",TS ,".rda"))
} else {
  x1MC=x2MC=RunPars=list()
  RunPars$x1MC = list()
  RunPars$x2MC = list()
  }
p=list()
p[[1]] =seq(from=0.05, to=0.0725, by=0.0025)#10 arms
p[[2]] =seq(from=0.05, length=10, by=4*0.0025)#also 10
# Case I

  NumSims=100
  maxIters = 1000
  #p[[1]] =seq(from=0.05, to=0.0725, by=0.0025)#10 arms
  Narms=length(p[[1]])

  ## 1a 
  (shape12 = BetaShape( mean(p[[1]]),var(p[[1]])))
  a1=shape12[1];b1=shape12[2]
  #a few small tests:
  #tmp=PlayBanditManyArms(c(a=a1,b=b1),maxIters = 40,Narms=10,alpha=1,beta=1)
  #tmp=PlayBanditMC(rep( list(c(a=a1,b=b1)),5),maxIters = 40,Narms=10,alpha=1,beta=1,NumCores =5)
  #View(RunPars$x1MC[[1]])
  
  t0=Sys.time()
  extra_x1MC = PlayBanditMC(rep( list(c(a=a1,b=b1)),NumSims),maxIters = maxIters, Narms=Narms, NumCores = NumSims,alpha=1,beta=1)
  print(Sys.time()-t0)
  tmp= attr(extra_x1MC,"params")
  if (firstRun) {
    x1MC[[1]] = extra_x1MC
    RunPars$x1MC[[1]] = tmp$ActualP
  } else {
    cat("dimensions before merge:", dim(x1MC[[1]]), "\n")
    x1MC[[1]] = abind(x1MC[[1]],extra_x1MC)
    RunPars$x1MC[[1]] = rbind(RunPars$x1MC[[1]], tmp$ActualP)
  }
  cat("1a: dimensions after merge:", dim(x1MC[[1]]), "\n")
  
if (1){
## 1b. Perfectly agreeing prior:
   
   (shape12 = BetaShape( mean(p[[1]]),var(p[[1]])))
   t0=Sys.time()
   extra_x1MC = PlayBanditMC(rep( list(c(a=a1,b=b1)),NumSims),maxIters = maxIters, Narms=Narms, 
                  NumCores = NumSims,alpha=a1,beta=b1)
   print(Sys.time()-t0)
   tmp= attr(extra_x1MC,"params")
   if (firstRun) {
     x1MC[[2]] = extra_x1MC 
     RunPars$x1MC[[2]] = tmp$ActualP
   } else {
     cat("dimensions before merge:", dim(x1MC[[2]]), "\n")
     x1MC[[2]] = abind(x1MC[[2]],extra_x1MC)
     RunPars$x1MC[[2]] = rbind(RunPars$x1MC[[2]], tmp$ActualP)
   }
   cat("1b: dimensions after merge:", dim(x1MC[[2]]), "\n")
   
# Case II   
   
  Narms=length(p[[2]])
  (shape12 = BetaShape( mean(p[[2]]),var(p[[2]])))
  a2=shape12[1];b2=shape12[2]
  ## 2a 
  t0=Sys.time()
  extra_x2MC = PlayBanditMC(rep( list(c(a=a2,b=b2)),NumSims),maxIters = maxIters, Narms=Narms, NumCores = NumSims,alpha=1,beta=1)
  print(Sys.time()-t0)
  tmp= attr(extra_x2MC,"params")
  if (firstRun) {
    x2MC[[1]] = extra_x2MC 
    RunPars$x2MC[[1]] = tmp$ActualP
  } else {
    cat("dimensions before merge:", dim(x2MC[[1]]), "\n")
    x2MC[[1]] = abind(x2MC[[1]],extra_x2MC)
    RunPars$x2MC[[1]] = rbind(RunPars$x2MC[[1]], tmp$ActualP)
  }
  cat("2a: dimensions after merge:", dim(x2MC[[1]]), "\n")
  
  ## 2b. Perfectly agreeing prior:
  
  (shape12 = BetaShape( mean(p[[2]]),var(p[[2]])))
  t0=Sys.time()
  extra_x2MC = PlayBanditMC(rep( list(c(a=a2,b=b2)),NumSims),maxIters = maxIters, Narms=Narms, 
                            NumCores = NumSims,alpha=a2,beta=b2)
  print(Sys.time()-t0)
  tmp= attr(extra_x2MC,"params")
  if (firstRun) {
    x2MC[[2]] = extra_x2MC 
    RunPars$x2MC[[2]] = tmp$ActualP
  } else {
    cat("dimensions before merge:", dim(x2MC[[2]]), "\n")
    x2MC[[2]] = abind(x2MC[[2]],extra_x2MC)
    RunPars$x2MC[[2]] = rbind(RunPars$x2MC[[2]], tmp$ActualP)
  }
  cat("2b: dimensions after merge:", dim(x2MC[[2]]), "\n")

  
  ## 2c. VERY disagreeing prior:
  
  
  (shape12 = BetaShape( 3*mean(p[[2]]),var(p[[2]])))  
  a3=shape12[1];b3=shape12[2]
  t0=Sys.time()
  extra_x2MC = PlayBanditMC(rep( list(c(a=a3,b=b3)),NumSims),maxIters = maxIters, Narms=Narms, 
                            NumCores = NumSims,alpha=a2,beta=b2)
  print(Sys.time()-t0)
  tmp= attr(extra_x2MC,"params")
  if (firstRun) {
    x2MC[[3]] = extra_x2MC 
    RunPars$x2MC[[3]] = tmp$ActualP
  } else {
    cat("dimensions before merge:", dim(x2MC[[3]]), "\n")
    x2MC[[3]] = abind(x2MC[[3]],extra_x2MC)
    RunPars$x2MC[[3]] = rbind(RunPars$x2MC[[3]], tmp$ActualP)
  }
  cat("2c: dimensions after merge:", dim(x2MC[[3]]), "\n")
}
  
## ------------------------------------------------------------------------
save(x1MC, x2MC, p, RunPars, file="xMC_bayes.rda")

