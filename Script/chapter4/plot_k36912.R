#Usage: Rscript --vanilla ../Script/chapter3/plot_k36912.R 

suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

#setwd("~/Dropbox (ASU)/Indel_project/chapter4")

readFile = function(inD, pat){
  Files= list.files(inD,pattern=pat,full.names=T)
  index = as.numeric(str_extract(basename(Files),'[^.]+')) #rectify the order
  Files = Files[order(index)]
  return(Files)
}

doMat18 = function(Files,n,nvar){
  dMat = matrix(0,n,nvar)
  colnames(dMat)=c('A','C','G','T','s1','s2','s3','s4','s5','s6','omega','tau','r0*t','r1*t','r2*t','ext.I','ext.D','omegaZ')
  for (i in 1:n) {
    Jmp       = fromJSON(Files[i])
    go        = Jmp$gap.openning
    go.comb   = c(go[1]+go[4],go[2]+go[5],go[3]+go[6])
    rt        = -log(1-go.comb)
    ge        = Jmp$gap.extension
    dMat[i,]  = c(Jmp$nuc.freq, Jmp$sigmas, Jmp$omega, Jmp$branch.length, rt, ge, Jmp$omega.z) 
  }
  return(dMat)
}

################################
readM = function(inD,n,nvar){
  Files   = readFile(inD,'est')
  res.mat = doMat18(Files,n,nvar)
  return(res.mat)
}





###################################################################

#read true input
#n=100;m=4
inF  = "trueP.100.txt"
tP   = read.table(inF,header=T)
tMat = as.matrix(tP)
tz   = tMat[,18]
n    = dim(tMat)[1]
nvar = dim(tMat)[2]


targ = c("k3","k6","k9","k12")
m    = length(targ)
ouD1 = "Figure/wz_k36912.pdf"

dz = matrix(0,n,m)
for(i in 1:m){
  inD1   = paste0(targ[i],"/Results/Gse")
  dMat1  = readM(inD1,n,nvar)
  dz[,i] = dMat1[,18]  
}

#error perc
epc = abs(dz-tz)/tz
colnames(epc) = targ
print(colMeans(epc))

edf = as.data.frame(as.table(epc))
colnames(edf)[2:3] = c("k.value","error_perc")

edF_sorted = edf %>% mutate(k.value = fct_reorder(k.value,error_perc))
#ggplot(edF_sorted, aes(x=k.value,y=error_perc,color=k.value)) + geom_boxplot() + coord_flip()

gs =  ggplot(edF_sorted, aes(x=k.value,y=error_perc,color=k.value)) + geom_boxplot() + coord_flip()

pdf(ouD1)
#gs + geom_boxplot(outlier.alpha=0) + geom_point(size=2, alpha=0.6) 
gs + geom_jitter(size=2,alpha=0.25,width=0.2) + stat_summary(fun=mean,geom='point',size=5)
