if (0){
  
  x0=PlayBanditManyArms(p=c(0.04,0.045,0.05),maxIters = 300)
  matplot(x0[,7:9], type="l",ylab="posterior");grid()
  plot(x0[,"regret"], type="l",ylab="regret")
  plot(x0[,"pvR"], type="l",ylab="value remaining")
  
  #quick test:
  NumSims=5
  maxIters = 6
  p=c(0.04,0.045,0.05)
  x0MC = PlayBanditMC(rep( list(p),NumSims),maxIters = maxIters , NumCores = NumSims)
  #test: 
  tmp = array(as.numeric(unlist(x0MC)), dim=c(maxIters, 3*length(p)+2, NumSims))
  dimnames(tmp) =list(paste0("batch",1:maxIters),colnames(x0MC[[1]]),paste0("sim",1:NumSims))
  x0=x0MC[[1]]
  matplot(x0[,7:9], type="l",ylab="posterior");grid()
  plot(x0[,"regret"], type="l",ylab="regret")
  plot(x0[,"pvR"], type="l",ylab="value remaining")
  
  #serious simulation 1
  NumSims=100
  maxIters = 1000
  p=seq(from=0.05, to=0.0725, by=0.0025)#10 arms
  t0=Sys.time()
  x0MC = PlayBanditMC(rep( list(p),NumSims),maxIters = maxIters , NumCores = NumSims)
  print(Sys.time()-t0)
  #test: 
  # tmp = x0MC#backup
  # x0MC = array(as.numeric(unlist(x0MC)), dim=c(maxIters, 3*length(p)+2, NumSims))
  # dimnames(x0MC) =list(paste0("batch",1:maxIters),colnames(tmp[[1]]),paste0("sim",1:NumSims))
  x0=x0MC[,,1]
  matplot(x0[,7:9], type="l",ylab="posterior");grid()
  plot(x0[-(1:5),"regret"],type="l",ylab="regret")#, ylim=c(0,1)
  plot(x0[,"pvR"], type="l",ylab="value remaining")
  
  boxplot(t(x0MC[seq(5,500,by=5),"regret",]),pch=20,cex=0.5);grid()
  
  #serious simulation 2
  NumSims=100
  maxIters = 1000
  p=seq(from=0.05, length=10, by=4*0.0025)#10 arms
  t0=Sys.time()
  x2MC = PlayBanditMC(rep( list(p),NumSims),maxIters = maxIters , NumCores = NumSims)
  print(Sys.time()-t0)
  # tmp = x2MC#backup
  # x2MC = array(as.numeric(unlist(x2MC)), dim=c(maxIters, 3*length(p)+2, NumSims))
  # dimnames(x2MC) =list(paste0("batch",1:maxIters),colnames(tmp[[1]]),paste0("sim",1:NumSims))
  x0=x2MC[,,1]
  matplot(x0[,7:9], type="l",ylab="posterior");grid()
  plot(x0[-(1:5),"regret"],type="l",ylab="regret")#, ylim=c(0,1)
  plot(x0[,"pvR"], type="l",ylab="value remaining")
  boxplot(t(x2MC[seq(5,1000,by=10),"regret",]),pch=20,cex=0.5);grid()
  #par(mfrow=c(3,3));
  for (j in 2:10){
    a=paste0("bba",j)
    boxplot(t(x2MC[seq(5,1000,by=10),a,]),pch=20,cex=0.5,main = p[j]);grid()
  } 
        
  
}
#Officially there are only two stopping criteria:
#By default, we force the bandit to run for at least two weeks. After that, we keep track of two metrics.
# (i) The first is the probability that each variation beats the original. If we're 95% sure that a variation beats the original then Google Analytics declares that a winner has been found. Both the two-week minimum duration and the 95% confidence level can be adjusted by the user.
# (ii) The second metric that we monitor is is the "potential value remaining in the experiment", which is particularly useful when there are multiple arms. At any point in the experiment there is a "champion" arm believed to be the best. If the experiment ended "now", the champion is the arm you would choose. The "value remaining" in an experiment is the amount of increased conversion rate you could get by switching away from the champion.
# We have added a "maximum number of iteration" filter to prevent stragglers
PlayBanditManyArms = function(
  p=seq(from=0.05, to=0.0725, by=0.0025),##<< true proportions
  alpha=1, ##<< prior alpha
  beta=1, ##<< prior beta
  Nday = 50, ##total impressions per day
  Ninit = Nday, ##14*Nday/2, ## burn-in
  power=0.95,
  #confLeveL=0.95,
  maxIters=NULL,
  Narms,
  #maxVR = 0.01, ##<< maximum remaining  Value (0.01 = 1%)
  verbose= 0
){
  require(bandit)
  #if (!dir.exists("Logs/done/")) dir.create("Logs/done/", showWarnings = FALSE, recursive =TRUE)
  
  if (!is.null(names(p)) & all(names(p) %in% c("a","b"))) {#no fixed p, instead we sample from a beta
    if (verbose) cat("no fixed p, instead we sample from a beta with params ", p, "\n")
    p = sort(round(rbeta(Narms, p["a"], p["b"]),4))
  } 
    
  stopifnot(Narms == length(p))
  
  pMax = max(p)
  #browser()
  #p=runif(Narms, p[1], p[2]) #random choices
  
  try({
    #pid=Sys.getpid()
    #logFile = paste0(pid,"-",paste0(10^4*p,collapse="-"), ".csv") 
    options(digits=4)
    #if (!dir.exists("Logs/done")) dir.create("Logs/done", recursive = TRUE)
    
    #sink(paste0("Logs/", logFile))
    #if (0) 
    if (missing(maxIters) | is.null(maxIters) ) {
      pt=power.prop.test(p1 = mean(p), p2 = mean(p)+mean(diff(p)), power=power, sig.level = 0.05)
      N = 2*round(pt$n)
      
      if (missing(maxIters)) maxIters = 2*N
      if (is.null(maxIters)) maxIters = 4*N
      
      if (verbose>=0) {
        #print(pt)
        cat("power.prop.test$n:", N, "enforcing maxIters=", maxIters, "\n")
      }
    }
    #initial "burn-in:
    

    x=matrix(0,nrow =maxIters,ncol=3*Narms+2)
    colnames(x)=c(round(p,5),paste0("N", 1:Narms),paste0("bba", 1:Narms),"pvR","regret")
    
    x[1,1:Narms] = rbinom(Narms, Ninit,p)
    x[1,paste0("N", 1:Narms)]= rep(Ninit,Narms)
    BestArm = best_binomial_bandit(x[1,1:Narms],x[1,paste0("N", 1:Narms)], alpha=alpha, beta=beta)
    vRBB = value_remaining_best_bandit(x[1,1:Narms], x[1,paste0("N", 1:Narms)], alpha=alpha, beta=beta)
    if (abs(sum(BestArm)-1)>0.01) BestArm = vRBB$postWin
    x[1,paste0("bba", 1:Narms)]= BestArm
    x[1,"pvR"] = quantile(vRBB$vR,0.95)
    x[1,"regret"] = pMax*sum(Ninit) - sum(p*Ninit)
    
    #cat("dim(x):", dim(x), "\n")
    #print(head(x),2)
    
    #while(max(BestArm)< confLeveL & pvR > maxVR & sum(N) < maxIters & abs(sum(BestArm)-1) < 0.01 ) {
    for (i in 2:maxIters) {
      #Thompson sampling assigns visits to arms in proportion to the probability that each arm is optimal:
      N12 =round(Nday * BestArm)
      currentDraw=rbinom(Narms, N12,p)
      #cumulative sums:
      cumX=x[i-1,1:Narms]+currentDraw
      x[i,1:Narms] = cumX
      cumN = x[i-1,paste0("N", 1:Narms)]+N12
      x[i,paste0("N", 1:Narms)]= cumN
      BestArm = best_binomial_bandit(cumX,cumN, alpha=alpha, beta=beta)
      vRBB = value_remaining_best_bandit(cumX, cumN, alpha=alpha, beta=beta)
      if (abs(sum(BestArm)-1)>0.01) BestArm = vRBB$postWin
      x[i,paste0("bba", 1:Narms)]= BestArm
      x[i,"pvR"] = quantile(vRBB$vR,0.95)
      x[i,"regret"] = pMax*sum(N12) - sum(p*N12)
    }  
    
    
    return(x)
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
  p=list(seq(from=0.05, to=0.0725, by=0.0025)),##<< true proportions OR params to beta!
  Narms, ##<< number of arms
  alpha=1, ##<< prior alpha
  beta=1, ##<< prior beta
  Nday = 50, ##total impressions per day
  Ninit = Nday, ##14*Nday/2, ## burn-in
  power=0.95,
  #confLeveL=0.95,
  maxIters=NULL,
  #Narms = length(p),
  NumCores,
  verbose= 0
){
  if (is.list(p)){
    if (missing(Narms)) Narms= length(p[[1]])
    Nsims= length(p)  
    if (missing(NumCores)) NumCores = Nsims
  }
  
  if (NumCores>1){
    if (require(parallel))
      res = mclapply(p, PlayBanditManyArms, mc.cores=NumCores, alpha=alpha, beta=beta, Ninit = Ninit, Nday = Nday, 
                     power=power, maxIters=maxIters,Narms=Narms, verbose=verbose)
  } else {
      res = lapply(p, PlayBanditManyArms, alpha=alpha, beta=beta, Ninit = Ninit, Nday = Nday, power=power, 
                   Narms=Narms, maxIters=maxIters, verbose=verbose)
  }
  
  if (verbose>1) return(res)
  try({
    successSims=!sapply(res,is.null)
    NsimsEff = sum(successSims)
    tmp = res[successSims]#backup
    #browser()
    res = array(as.numeric(unlist(tmp)), dim=c(maxIters, 3*Narms+2, NsimsEff))
    dimnames(res) =list(paste0("batch",1:maxIters),colnames(tmp[[1]]),paste0("sim",1:NsimsEff))
  })
  ret = try({
    ActualP = matrix(0,nrow=NsimsEff,ncol=Narms)
    for (i in 1:NsimsEff) ActualP[i,]=as.numeric(colnames(tmp[[i]])[1:Narms])
    attr(res,"params") = list(betaPars=c(alpha=alpha, beta=beta),ActualP=ActualP)
  })
  
  stopifnot(nrow(ActualP) == dim(res)[[3]])
  
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

value_remaining_best_bandit = function#
(x, n, 
 alpha = 1, beta = 1, 
 ndraws = 10000,
 verbose=0
) 
{
  #posterior distribution the Bayesian probabilities for each arm 
  #matrix with ndraws rows and length(x) cols:
  post = bandit::sim_post(x, n, alpha, beta, ndraws)
  #vector of length(x): prob. being the winner
  postWin = bandit::prob_winner(post)
  iMax = which.max(postWin)
  thetaMax = apply(post, 1, max)
  vR = (thetaMax - post[, iMax])/post[, iMax]
  return(list(vR=vR,postWin=as.numeric(postWin)))
}

BetaShape= function(#get beta shape parameters from mean and variance
  m=0.1, ##<< mean
  v=0.1^2 ##<< variance
){
  a = -m*(m^2-m+v)/v
  b=(m-1)*(m^2-m+v)/v
  return(c(a,b))
}

modifyRegret = function(
  x2MC, 
  ActualP = RunPars$x2MC,
  Narms=10,
  l, 
  getPhat = TRUE, ##<< also compute the estimated arm probabilities at termination
  BestArmCutoff=0.95,
  ValRem=0.01, ##<< optional stopping
  verbose = FALSE
){
  
  stopifnot(nrow(ActualP[[l]]) == dim(x2MC[[l]])[3])
  Nsims = dim(x2MC[[l]])[3]
  maxIters=dim(x2MC[[l]])[1]
  
  regret=x2MC[[l]][,"regret",]
  if (getPhat){
    TestStats = matrix(NA,nrow=Nsims, ncol=4+Narms)
    colnames(TestStats) = c("Nterminal", "pMax","pTerminal", "iTerminal", paste0("pHat", 1:Narms))
  } else {
    TestStats = matrix(NA,nrow=Nsims, ncol=4)
    colnames(TestStats) = c("Nterminal", "pMax","pTerminal", "iTerminal")
  }
  
  
  #postBestArm = apply(x2MC[[l]][,paste0("bba",1:Narms),],1,which.max)
  #postBestArm = postBestArmVal = matrix(0,nrow=dim(x2MC[[l]])[1],ncol=dim(x2MC[[l]])[3])
  #sigh, I need a for loop, could not figure out how to work this with the apply family
  #for (i in 1:nrow(postBestArm)){#time
  for (j in 1:Nsims){#sim replication
    #postBestArm[i,j] = which.max(x2MC[[l]][i,paste0("bba",1:Narms),j])
    #postBestArmVal[i,j] = 
    #maxPost = apply(x2MC[[l]][,paste0("bba",1:Narms),j],1,which.max)
    maxPostVal = apply(x2MC[[l]][,paste0("bba",1:Narms),j],1,max,na.rm=TRUE)
    
    OS1 = match(TRUE, maxPostVal > BestArmCutoff, nomatch =maxIters)
    OS2 = match(TRUE, x2MC[[l]][,"pvR",] < ValRem, nomatch =maxIters)
    OS12 = pmin(OS1,OS2)
    #from here one we stop the test and reap reward of the beat arm only:
    if (OS12<maxIters) {
      bestArm = which.max(x2MC[[l]][OS12,paste0("bba",1:Narms),j])
      pTerminal = ActualP[[l]][j,bestArm]
      pMax=max(ActualP[[l]][j,],na.rm=TRUE)
      TotalN= diff(rowSums(x2MC[[l]][(OS12):maxIters,paste0("N",1:Narms),j]))
      regret[(OS12+1):maxIters,j] = (pMax-pTerminal)*TotalN
      Nterminal= sum(x2MC[[l]][OS12,paste0("N",1:Narms),j])
      TestStats[j,1:4] = c(Nterminal, pMax,pTerminal, OS12)
    }
    if (getPhat){#compute the estimated arm probabilities at termination
      browser()
      TestStats[j,paste0("pHat", 1:Narms)] = x2MC[[l]][OS12,1:Narms,j]/x2MC[[l]][OS12,paste0("N",1:Narms),j]
    }
  }
  
  return(list(regret=regret, TestStats=TestStats ))
}

plotRegret = function(x2MC, 
                      l, 
                      p,
                      qLines = c(0.1,0.9), ##<< add quantile lines
                      batches2plot=seq(5,1000,by=10),
                      post2plot= c(2,4,6:10), ##<< posteriors to plot
                      regret = NULL,#
                      mainRegret ="regret",
                      ylim = NULL,
                      extraDiagnostics = FALSE){
  #rownames(x2MC[[l]][,,1]) = 1:nrow(x2MC[[l]][,,1]) 
  
  
  
  if (extraDiagnostics){
    par(mfrow=c(1,3));
    x0=x2MC[[l]][,,1]
    matplot(x0[,paste0("bba", 1:Narms)], type="l",ylab="posterior");grid()
    plot(x0[-(1:5),"regret"],type="l",ylab="regret")#, ylim=c(0,1)
    plot(x0[,"pvR"], type="l",ylab="value remaining")
  } else {
    
    if (!is.null(post2plot)) {
      mfrow=c(2,4) 
      par(mfrow=mfrow)
    } #else mfrow=c(1,1) 
  
    if (is.null(regret)) regret=x2MC[[l]][,"regret",]
    N=nrow(regret)
     #, border=rgb(1,0.9,0.75,0.5)
    Xbxp=boxplot(t(regret[batches2plot,]),pch=20,cex=0.25,outcol=rgb(0.55,0,1,0.3), ylab="regret",
            outbg=rgb(0.55,0,1,0.3) , col="bisque", ylim=ylim, xlab = "num batches",
            main = mainRegret,names=(1:N)[batches2plot]);
    grid();
    #mtext(text = "regret",side=2,line=2)
    #mtext(text = "num batches",side=1,line=2)
    #browser()
    if (!is.null(qLines)) {
      q =apply(regret[batches2plot,],1,quantile,p=qLines) 
      for (i in 1:length(qLines)) lines(smooth.spline(1:length(batches2plot),q[i,]),col="brown",lwd=1.85)
    }
    
    if (!is.null(post2plot)){
      p=colnames(x2MC[[l]][,,1])#[1:Narms]
      for (j in post2plot){
        a=paste0("bba",j)
        boxplot(t(x2MC[[l]][batches2plot,a,]),pch=20,cex=0.35,outcol=rgb(0,0,1,0.4),
                outbg=rgb(0,0,1,0.4), col="bisque",main = p[j], ylab = "posterior",names=(1:N)[batches2plot]);grid()
        if (!is.null(qLines)) {
          q =apply(x2MC[[l]][batches2plot,a,],1,quantile,p=qLines) 
          for (i in 1:length(qLines)) lines(smooth.spline(1:length(batches2plot),q[i,]),col="brown",lwd=1.85)
        }
      }
    }
  }
}