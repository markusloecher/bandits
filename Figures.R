
########################### Figure 1 ##########################

baseR=TRUE
fname="./figures/RegretNoStopping.pdf"
if (baseR){
  #m=8;png(fname, m*500,m*350)
  m=1;pdf(fname, m*8,m*3.5)
  par(mfrow=c(2,3), mar=c(2,2,1,1),cex=0.7)
  for (l in 1:3){
    if (l==1) {
      par(mar=c(4,4,1,1)) 
    } else {
      par(mar=c(2,2,1,1))
    }
    plotRegret(regret=OrigRegret$x2MC[[l]], 
               post2plot = NULL, mainRegret="",ylim=c(0,5))
  }
  for (l in 1:3)
    plotRegret(regret=OrigRegret$x1MC[[l]], 
               post2plot = NULL, mainRegret="",ylim=c(0,2))
  dev.off()
} else {
  l=1
  GGregret = cbind.data.frame(OrigRegret$x2MC[[l]], l=l, EffectSize = "large")
  tmp = cbind.data.frame(OrigRegret$x1MC[[l]], l=l, EffectSize = "small")
  GGregret = rbind(GGregret,tmp)
  
  #for (OS in c(TRUE, FALSE))
  for (l in 2:3) {
    tmp = cbind.data.frame(OrigRegret$x2MC[[l]], l=l, EffectSize = "large")
    GGregret = rbind(GGregret,tmp)
    
    tmp = cbind.data.frame(OrigRegret$x1MC[[l]], l=l, EffectSize = "small")
    GGregret = rbind(GGregret,tmp)
    
  }
  
  
  GGregret$l = factor(GGregret$l)
  GGregret$EffectSize = factor(GGregret$EffectSize)
  
}

########################### Figure 2 ##########################

#baseR=TRUE
load("cumRegret2.rda")

g1 <- subset(cumRegret2, EffectSize=="large") %>%
  ggplot( aes(x=cumRegret2)) +
  geom_histogram( color="black", fill="white")  + labs(fill="") + 
  facet_wrap(. ~   l,ncol=3) + scale_x_continuous(limits=c(0,300)) #+ scale_x_sqrt(limits=c(0,500)) #+ xlim(0,500)#, scales="free")

#p + scale_x_log10() + scale_y_sqrt()
g1 = g1 + scale_y_continuous(limits=c(0,175)) + xlab("cumulative regret") #+ scale_y_sqrt()

fname="./figures/CumRegret1.pdf"
m=3/4;ggsave(fname,g1,width=m*8,height=m*4)

########################### Figure 3 ##########################

g1 <- subset(TerminalRegret, EffectSize=="large") %>%
  ggplot( aes(x=TerminalRegret)) +
  geom_histogram( color="black", fill="white")  + labs(fill="") + facet_wrap(. ~   l,ncol=3)

g1 = g1  + scale_y_sqrt() + xlab("terminal regret [%]")
fname="./figures/TerminalRegretLarge.pdf"
m=3/4;ggsave(fname,g1,width=m*8,height=m*4)

########################### Figure 4 ##########################

g1 <- subset(TerminalRegret, EffectSize=="small") %>%
  ggplot( aes(x=TerminalRegret)) +
  geom_histogram( color="black", fill="white")  + labs(fill="") + facet_wrap(. ~l,ncol=3)

g1 = g1  + scale_y_sqrt() + xlab("terminal regret [%]")#+ scale_x_log10()
fname="./figures/TerminalRegretSmall.pdf"
m=3/4;ggsave(fname,g1,width=m*8,height=m*4)

########################### Figure 5 ##########################
load("pHatLong.rda")
mpHat = aggregate(bias ~ EffectSize +  l , data = subset(pHatLong, arm=="arm10"), FUN = mean, na.rm=TRUE)
g1  <-   subset(pHatLong, arm=="arm10") %>%
  ggplot( aes_string(x=paste0("bias"))) +  geom_vline(xintercept=0)  + xlim(c(-100, 100)) +
  geom_histogram( color="black", fill="white")  + labs(fill="") + facet_wrap(. ~ EffectSize +  l,ncol=3)

g1 = g1 + geom_vline(aes_string(xintercept = "bias"), mpHat,col=2) + xlab("bias [%]") + scale_y_sqrt()

fname="./figures/EstimationArm10.pdf"
m=3/4;ggsave(fname,g1,width=m*8,height=m*6)



########################### Figure 6 ##########################
Fig6 = function(m=1){
  load("nullPVR.rda")
  fname="./figures/NullSimulation.pdf"
  
  pdf(fname, m*5,m*3.5)
  
  #par(mfrow=c(2,1))
  batches2plot=seq(5,1000,by=10)
  l=1#for (l in 1:2){
    Xbxp=boxplot(t(nullPVR[[l]][batches2plot,]),pch=20,cex=0.35,outcol=rgb(0.55,0,1,0.4), ylab="value remaining",
                 outbg=rgb(0.55,0,1,0.4) , col="bisque", main = c("","Uniform prior","Informed prior")[l],
                 names=NULL,log="y",xlab="10xbatch", ylim = c(0.005,3));
    grid()
    abline(h=0.01,col=2, lty=2)
    qLines = c(0.1,0.9)
    q =apply(nullPVR[[l]][batches2plot,],1,quantile,p=qLines) 
    for (i in 1:length(qLines)) lines(smooth.spline(1:length(batches2plot),q[i,]),col="brown",lwd=1.5)
  #}
dev.off()
}
Fig7 = function(){
  load("pList.rda")
  SummaryStats = as.data.frame(SummaryStats)
  SummaryStats$Delta = 0.05-unlist(pList)[seq(2,2*length(pList),by=2)][-10]
  SummaryStats$Narms=2
  SummaryStats$p2 = unlist(pList)[seq(2,2*length(pList),by=2)][-10]
  SS = SummaryStats
  
  #multiple arms:
  load("pListMA.rda")
  SummaryStats$Delta = 0.05-SummaryStats[,"p2"]
  SS=SS[,colnames(SummaryStats)]
  SS = rbind.data.frame(SS,SummaryStats)
  SS$Delta = 10*SS$Delta
  SS$Nfreq = NA#frequentist sample size
  SS$Deltak = SS$Delta/(as.numeric(as.character(SS$Narms))-1)
  
  for (i in 1:nrow(SS)){
    p2=seq(from=0.05, to=SS$p2[i], length=SS$Narms[i])[2]
    #print(p2)
    #p2 = 0.05 + (SS$p2[i]-0.05)/SS$Narms[i]
    #SS$Nfreq[i] =ceiling(power.prop.test(p1 = 0.05, p2 = p2, power=SS$Accuracy[i], sig.level = 1)$n)
  }
  
  SS$Narms = factor(SS$Narms)
  
  g1 = ggplot(SS, aes(x=Delta,y=Accuracy)) + geom_line(aes(lty=Narms,colour=Narms), size=2) +labs(x=expression(-10*Delta)) + theme_bw() # + geom_point(aes(colour=Narms), size=2)
  g1 = g1 + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  g2 = ggplot(SS, aes(x=Delta,y=MeanIteration)) + geom_line(aes(lty=Narms,colour=Narms), size=2) + scale_y_log10() +labs(x=expression(-10*Delta), y = "ASN") + theme_bw()
  g2 = g2 + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  g12=grid.arrange(g1+ theme(legend.position=c(0.8,0.4)), g2+ theme(legend.position="none"), ncol=1)
  
  #p12
  fname="./figures/Accuracy_ASN.pdf"
  ggsave(fname, g12, width=5,height=4, device=cairo_pdf)
  
}