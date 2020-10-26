if (0){
  PlayBandit(p= c(0.05,0.04))
  PlayBandit(p= c(0.05,0.0495),verbose=1)
  rep( list(c(0.05,0.04)) ,5)
  PlayBanditMC(rep( list(c(0.05,0.04)) ,5), NumCores = 1)
  
  PlayBanditMC(rep( list(c(0.05,0.04)) ,5), NumCores = 5, Narms = 5)
  
  pList = list()
  for (p in  rev(c(1,2,3,5,7.5,10,15,20, 25)) ){
    p2 = 0.05*(1-p/100)
    pList[[as.character(p)]] = c(0.05, p2)
  }
  
  PlayBanditMC(pList[[1]], NumCores = 1,verbose=1)
  
  for (p in pList){
    pt=power.prop.test(p1 = p[1], p2 = p[2], power=0.95, sig.level = 0.5)
    cat("working on", p, ", max N:",format(2*round(pt$n)), "\n")
    res=PlayBanditMC(rep(list(p),each=500), NumCores = 100, maxIters = 4*round(pt$n))
    
    save(res,file=paste0("res-",paste0(10^4*p,collapse="-"), ".rda"))
  }
    
  pList = rep( pList,each=100)
  
  res=PlayBanditMC(pList, NumCores = 50)
}
#Officially there are only two stopping criteria:
#By default, we force the bandit to run for at least two weeks. After that, we keep track of two metrics.
# (i) The first is the probability that each variation beats the original. If we're 95% sure that a variation beats the original then Google Analytics declares that a winner has been found. Both the two-week minimum duration and the 95% confidence level can be adjusted by the user.
# (ii) The second metric that we monitor is is the "potential value remaining in the experiment", which is particularly useful when there are multiple arms. At any point in the experiment there is a "champion" arm believed to be the best. If the experiment ended "now", the champion is the arm you would choose. The "value remaining" in an experiment is the amount of increased conversion rate you could get by switching away from the champion.
# We have added a "maximum number of iteration" filter to prevent stragglers
PlayBanditManyArms = function(p= c(0.05,0.08),##<< lower and upper bound of true proportions
                      Nday = 100, ##total impressions per day
                      Ninit = 14*Nday/2, ## burn-in
                      power=0.95,
                      confLeveL=0.95,
                      maxIters=NULL,
                      Narms = 4,
                      maxVR = 0.01, ##<< maximum remaining  Value (0.01 = 1%)
                      verbose= 0
){
  require(bandit)
  if (!dir.exists("Logs/done/")) dir.create("Logs/done/", showWarnings = FALSE, recursive =TRUE)
               
  #p=runif(Narms, p[1], p[2]) #random choices
  p=seq(from=p[1], to=p[2], length=Narms)
  try({
    pid=Sys.getpid()
    logFile = paste0(pid,"-",paste0(10^4*p,collapse="-"), ".csv") 
    options(digits=4)
    #if (!dir.exists("Logs/done")) dir.create("Logs/done", recursive = TRUE)
    
    #sink(paste0("Logs/", logFile))
    #if (0) 
    if (missing(maxIters) | is.null(maxIters) ) {
      pt=power.prop.test(p1 = p[1], p2 = p[2], power=power, sig.level = 0.5)
      N = 2*round(pt$n)
      
      if (missing(maxIters)) maxIters = 2*N
      if (is.null(maxIters)) maxIters = 4*N
      
      if (verbose) {
        #print(pt)
        cat("power.prop.test$n:", N, "enforcing maxIters=", maxIters, "\n")
      }
    }
    #initial "burn-in:
    #x1 = rbinom(1,Ninit,p1)
    x = rbinom(Narms,Ninit,p)
    
    N=rep(Ninit,Narms)
    BestArm = best_binomial_bandit(x,N)
    
    BestArmSequence=vector()
    BestArmSequence[1] = BestArm[1]
    k=2
    #potential value remaining:
    pvR = 0.1
    while(max(BestArm)< confLeveL & pvR > maxVR & sum(N) < maxIters & abs(sum(BestArm)-1) < 0.01 ) {
      #Thompson sampling assigns visits to arms in proportion to the probability that each arm is optimal:
      
      N12 =round(Nday * BestArm)
      
      x = x + rbinom(Narms, N12,p)
      #x2 = x2 + rbinom(1, N12[2],p2)
      N = N + N12
      
      BestArm = as.numeric(best_binomial_bandit(x, N))
      #potential value remaining:
      vRBB = value_remaining_best_bandit(x, N)
      if (abs(sum(BestArm)-1)>0.01) BestArm = vRBB$postWin
      
      pvR = quantile(vRBB$vR,0.95)
      
      BestArmSequence[k] = BestArm[1]
      k=k+1
      ret=matrix(round(c(k,x, N, BestArm, pvR, p),4), nrow=1)
      #browser()
      #every 10000
      if (k %% 10^2 == 0){
        #print(ret)
        write.table(ret, paste0("Logs/", logFile), append=TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE)
      }
    }
    ret=matrix(round(c(k,x, N, BestArm, pvR, p),4), nrow=1)
    write.table(ret, paste0("Logs/", logFile), append=TRUE,row.names = FALSE, col.names = FALSE, quote=FALSE)
    #print(ret)
    
    BestArm=round(BestArm*1000)
    pvR=round(pvR*1000)
    ret=c(x,N,BestArm, pvR, p)
    names(ret) = c(paste0("x",1:Narms),paste0("N",1:Narms),paste0("BA",1:Narms), "pVR", paste0("p",1:Narms))
    
    if (verbose>1) return(list(BestArmSequence,ret))
    
    cmd = paste0("mv Logs/", logFile, " Logs/done/", logFile)
    try(system(cmd))
    return(ret)
  })#end of try wrapper
  return(NULL)
}


PlayBandit = function(p= c(0.05,0.08),##<<true proportions
                      Nday = 100, ##total impressions per day
                      Ninit = 14*Nday/2, ## burn-in
                      power=0.95,
                      confLeveL=0.95,
                      maxIters=NULL,
                      maxVR = 0.01, ##<< maximum remaining  Value (0.01 = 1%)
                      verbose= 0
){
  require(bandit)
  p1 = p[1]; p2 = p[2]
  try({
    pid=Sys.getpid()
    logFile = paste0(pid,"-",paste0(10^4*p,collapse="-"), ".csv") 
    options(digits=4)
    #if (!dir.exists("Logs/done")) dir.create("Logs/done", recursive = TRUE)
    
    #sink(paste0("Logs/", logFile))
    #if (0) 
    if (missing(maxIters) | is.null(maxIters) ) {
      pt=power.prop.test(p1 = p1, p2 = p2, power=power, sig.level = 0.5)
      N = 2*round(pt$n)
      
      if (missing(maxIters)) maxIters = 2*N
      if (is.null(maxIters)) maxIters = 4*N
      
      if (verbose) {
        #print(pt)
        cat("power.prop.test$n:", N, "enforcing maxIters=", maxIters, "\n")
      }
    }
    #initial "burn-in:
    x1 = rbinom(1,Ninit,p1)
    x2 = rbinom(1,Ninit,p2)
    
    N1=N2=Ninit
    BestArm = best_binomial_bandit(c(x1,x2),c(N1,N2))
    
    BestArmSequence=vector()
    BestArmSequence[1] = BestArm[1]
    k=2
    #potential value remaining:
    pvR = 0.1
    while(max(BestArm)< confLeveL & pvR > maxVR & N1+N2 < maxIters & abs(sum(BestArm)-1) < 0.01 ) {
      #Thompson sampling assigns visits to arms in proportion to the probability that each arm is optimal:
      
      N12 =round(Nday * BestArm)
      
      x1 = x1 + rbinom(1, N12[1],p1)
      x2 = x2 + rbinom(1, N12[2],p2)
      N1 = N1 + N12[1];N2 = N2 + N12[2]
      
      BestArm = as.numeric(best_binomial_bandit(c(x1,x2),c(N1,N2)))
      #potential value remaining:
      vRBB = value_remaining_best_bandit(c(x1,x2),c(N1,N2))
      if (abs(sum(BestArm)-1)>0.01) BestArm = vRBB$postWin
      
      pvR = quantile(vRBB$vR,0.95)
      
      BestArmSequence[k] = BestArm[1]
      k=k+1
      ret=matrix(round(c(k=k,x1=x1,x2=x2,N1=N1,N2=N2, bb1=BestArm[1], bb2=BestArm[2], pvR=pvR, p1=p1, p2=p2),4), nrow=1)
      #browser()
      #every 10000
      if (k %% 10^2 == 0){
        #print(ret)
        write.table(ret, paste0("Logs/", logFile), append=TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE)
      }
    }
    ret=matrix(round(c(k=k,x1=x1,x2=x2,N1=N1,N2=N2, bb1=BestArm[1], bb2=BestArm[2], pvR=pvR, p1=p1, p2=p2),4), nrow=1)
    write.table(ret, paste0("Logs/", logFile), append=TRUE,row.names = FALSE, col.names = FALSE, quote=FALSE)
    #print(ret)
    
    BestArm=round(BestArm*1000)
    pvR=round(pvR*1000)
    ret=c(x1=x1,x2=x2,N1=N1,N2=N2, bb1=BestArm[1], bb2=BestArm[2], pvR=pvR, p1=p1, p2=p2)
    
    if (verbose>1) return(list(BestArmSequence,ret))
    
    cmd = paste0("mv Logs/", logFile, " Logs/done/", logFile)
    try(system(cmd))
    return(ret)
  })#end of try wrapper
  cat("failed x/n:",c(x1,x2),c(N1,N2), best_binomial_bandit(c(x1,x2),c(N1,N2)), "\n")
  return(NULL)
}

PlayBanditMC = function(
   p= list(c(0.05,0.08)),##<<true proportions
   Ntotal = 5000,
   Nday = 100, ##total impressions per day
   Ninit = 14*Nday/2, ## burn-in
   Narms = 2,
   power=0.95,
   confLeveL=0.95,
   NumCores = 20,
   maxIters=NULL,
   verbose= 0
){
  if (NumCores>1){
    if (require(parallel))
      if (Narms ==2) {
        res = mclapply(p, PlayBandit, mc.cores=NumCores, Ninit = Ninit, Nday = Nday, power=power, maxIters=maxIters, verbose=verbose)    
      } else {
        res = mclapply(p, PlayBanditManyArms, mc.cores=NumCores, Ninit = Ninit, Nday = Nday, power=power, maxIters=maxIters, Narms=Narms, verbose=verbose)
      }
    
  } else {
    if (Narms ==2) {
      res = lapply(p, PlayBandit, Ninit = Ninit, Nday = Nday, power=power, maxIters=maxIters, verbose=verbose)
    } else {
      res = lapply(p, PlayBanditManyArms, Ninit = Ninit, Nday = Nday, power=power, maxIters=maxIters, Narms=Narms, verbose=verbose)
    }
  }
  try({res=do.call("rbind", res)})
  return(res)

}

 mysample=function(p,N) return(sum(sample(c(1,0),N, replace=TRUE, prob = c(p, 1-p))))

SimplePowerMC = function(
               NumCores = 50,
               Ntotal = 5000,
               p= c(0.05,0.08),##<<true proportions
               verbose= 0
){
  require(parallel)
  Nmin=power.prop.test(p1 = p[1], p2 = p[2], power=0.95)$n 
  
  CoinFlips1 = mclapply(1:Ntotal, mysample, mc.cores=NumCores, N=Nmin, p = p[1])
  CoinFlips2 = mclapply(1:Ntotal, mysample, mc.cores=NumCores, N=Nmin, p = p[2])
  
  res=cbind.data.frame(CoinFlips1,CoinFlips2)

  pooledP = rowSums(res[,1:2])/(Nmin + Nmin)
  SE = sqrt( pooledP * ( 1 - pooledP ) * ( (1/Nmin) + (1/Nmin) ) )
  res[,"zStat"] = (res[,2]-res[,1])/(SE*Nmin)
  
  return(list(CoinFlips=res, Nmin=Nmin,p=p))

}

SimpleDiagnostics = function(
  res, 
  p = c(0.05,0.04),
  verbose = 1, 
  maxIters=30000 #assuming that arm1 in fact is superior!
){
  pvRcol = grep("pvR", colnames(res))
  N=nrow(res)
  #res = cbind(res,N=rowSums(res[,3:4]))
  res = as.data.frame(res)
  res$N = rowSums(res[,3:4])
  
  Miters = mean(res[,"N"])
  if ((p[1]-p[2])/p[1] <= 0.010000001) {
    print("counting RV<1% as a correct decision")
    Acc = res[,"bb1"]/1000 >=0.95 | (res[,"bb1"]>res[,"bb2"] & res[,pvRcol] <=10) 
  } else {
    Acc = res[,"bb1"]/1000 >=0.95 | (res[,"bb1"]>res[,"bb2"] & res[,pvRcol] <=10)
  }
  Wrong = res[,"bb2"]/1000 >=0.95 | (res[,"bb1"]<=res[,"bb2"] & res[,pvRcol] <=10) 
  Undecided = res[,"bb1"]/1000 < 0.95 & res[,"bb2"]/1000 < 0.95 & res[,pvRcol] >10
  ELSE = !Acc & !Wrong & !Undecided
  
  head(res[Undecided,])
  #browser()
  # res$pVal1 = NA
  # res$pVal2 = NA
  # 
  # for (i in 1:N){
  #   res$pVal1[i] = 
  # }
  #   
  pt1=power.prop.test(p1 = p[1], p2 = p[2], power=0.95, sig.level = 0.5, alternative = "one.sided")
  pt2=power.prop.test(p1 = p[1], p2 = p[2], power=mean(Acc), sig.level = 0.5, alternative = "one.sided")
  
  
  
  if (verbose){
    k=sum(res[,"bb1"]/1000 >=0.95)
    cat("The first arm was diagnosed to be superior ", k, "times out of ", N, "\n")
    k=sum(res[,"bb2"]/1000 >=0.95)
    cat("The second arm was diagnosed to be superior ", k, "times out of ", N, "\n")
    k=sum(rowSums(res[,3:4])>=maxIters)
    cat("No early stopping was possible ", k, "times out of ", N, "\n")
    k=sum(res[,pvRcol] <=10)
    cat("The potential value remaining was below 1% ", k, "times out of ", N, "\n")
    
    cat("The first arm was correctly identified ", sum(Acc), "times out of ", N, "\n")
    
    cat("The mean number of iterations was ", Miters, "\n")
  
    Ptotal=rowSums(res[,c("bb1", "bb2")])
    cat("The total probability was less than 0.95", 100*mean(Ptotal<950), "% \n")
  
    #new type-I error:
    k=sum((res[,"bb2"]>res[,"bb1"] & res[,pvRcol] <=10) )
    
    cat("Assuming the first arm was superior, the 2nd arm was accepted as superior with value remaining less than 1%", k, "times out of ", N, "\n")
  }
  
  
  SummaryStats =c(Accuracy = sum(Acc)/N, Undecided = sum(Undecided)/N, WrongDecision = sum(Wrong)/N, MeanIteration = Miters, Nfreq1 =2*pt1$n,  Nfreq2 =2*pt2$n)
  
  invisible(SummaryStats)
}

SimpleDiagnosticsMA = function( #multiple arms
  res, 
  p = c(0.05,0.04),
  verbose = 1, 
  maxIters=30000 #assuming that arm1 in fact is superior!
){
  #names(res) = c(paste0("x",1:Narms),paste0("N",1:Narms),paste0("BA",1:Narms), "pVR", paste0("p",1:Narms))
  #Narms = length(p)
  Narms = (ncol(res)-1)/5
  if (missing(p)) p = res[1,paste0("p",1:Narms)]
    
  pvRcol = "pVR" # grep("pvR", colnames(res))
  
  N=nrow(res)
  #res = cbind(res,N=rowSums(res[,3:4]))
  res = as.data.frame(res)
  res$N = rowSums(res[,paste0("N",1:Narms)])
  
  Miters = mean(res[,"N"])
  # if ((p[1]-p[2])/p[1] <= 0.010000001) {
  #   print("counting RV<1% as a correct decision")
  #   Acc = res[,"BA1"]/1000 >=0.95 | (res[,"BA1"]>res[,"BA2"] & res[,pvRcol] <=10) 
  # } else {
  #   Acc = res[,"BA1"]/1000 >=0.95 | (res[,"BA1"]>res[,"BA2"] & res[,pvRcol] <=10)
  # }
  otherBAs = apply(res[,paste0("BA",2:Narms)],1,max)
  Acc = res[,"BA1"]/1000 >=0.95 | (res[,"BA1"]> otherBAs & res[,pvRcol] <=10)
  Wrong = res[,"BA2"]/1000 >=0.95 | (res[,"BA1"]<= otherBAs & res[,pvRcol] <=10) 
  Undecided = res[,"BA1"]/1000 < 0.95 & otherBAs/1000 < 0.95 & res[,pvRcol] >10
  ELSE = !Acc & !Wrong & !Undecided
  
  #head(res[Undecided,])
  
  SummaryStats =c(Accuracy = sum(Acc)/N, Undecided = sum(Undecided)/N, WrongDecision = sum(Wrong)/N, MeanIteration = Miters) #, Nfreq1 =2*pt1$n,  Nfreq2 =2*pt2$n)
  
  invisible(SummaryStats)
}


SummarizeBandit = function(
  res,##<< result from PlayBanditMC()
  p= c(0.05,0.08),##<<true proportions
  alpha = c(0.1,0.25,0.5,0.75,0.95),##<< 1-significance levels for power computation
  Nmin##<<minimum number of trial needed
){
  k = round(100*sum(res[,2]>res[,1])/nrow(res),3)
  cat(k, "percent correctly identified \n Selected quantiles:\n")
  
  NumTrials = rowSums(res[,c("N1","N2")])
  AvgNumTrials = round(mean(NumTrials))
  hist(NumTrials, xlab = "number of trials", main ="")
  cat("avg. num of trials for the bandits:", AvgNumTrials, "\n")
  abline(v=AvgNumTrials,col="darkgreen",lty=2)
  print(quantile(NumTrials, probs = c(0.25,0.5,0.75, 0.95)))
  
  if (missing(Nmin) ) {
    Nmin = rep(NA,length(alpha))
    names(Nmin) = as.character(alpha)
    for (a in alpha)
      Nmin[as.character(a)] =round(power.prop.test(p1 = p[1], p2 = p[2], power=0.95, sig.level = 1-a)$n) # 1752
  }
  
  if (!is.null(Nmin)){
    for (i in 1:length(Nmin)){
      #cat("frequentist power calculation demands 2*", Nmin[i], "trials to achieve a power of", 0.95,"at significance level ",names(Nmin)[i], "\n")
      abline(v=2*Nmin[i],col=2,lty=2)
      #mtext("Nmin", side=1, at =2*Nmin[i], line = -2, col = 2)
    }
  }
}


PlotBeta =  function (x, n, alpha = 1, beta = 1) 
{
  k <- length(x)
  ans <- numeric(k)
  for (i in (1:k)) {
    indx <- (1:k)[-i]
    f <- function(z) {
      r <- dbeta(z, x[i] + alpha, n[i] - x[i] + beta)
      for (j in indx) {
        r <- r * pbeta(z, x[j] + alpha, n[j] - x[j] + 
                         beta)
      }
      return(r)
    }
    #ans[i] = integrate(f, 0, 1)$value
    curve(f,0,1, xname = "z");grid()
  }

}

####################################

best_binomial_bandit_sim = function(x, n, alpha = 1, beta = 1, ndraws = 10000) return(prob_winner(sim_post(x,n,alpha, beta, ndraws)))

value_remaining_best_bandit = function (x, n, alpha = 1, beta = 1, ndraws = 10000) 
{
  post = bandit::sim_post(x, n, alpha, beta, ndraws)
  postWin = bandit::prob_winner(post)
  iMax = which.max(postWin)
  thetaMax = apply(post, 1, max)
  vR = (thetaMax - post[, iMax])/post[, iMax]
  return(list(vR=vR,postWin=as.numeric(postWin)))
}