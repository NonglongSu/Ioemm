# Apply the importance sampling Wi = g/f 
# g(x)--ziqi's phase_coati model; f(x)--juan's coatiM 

# A -- ancestor, B -- descendent
# insertion and deletion length distribution seperated.
# inserting nucleotides are ignored like deletions.

suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(dfoptim))
suppressPackageStartupMessages(library(profvis))

#setwd("~/Dropbox (ASU)/Indel_project/chapter4")

######################################
#Count the nucleotide freq
count_freq = function(input){
  nuc.count = 0
  for (i in input){
    dna       = readDNAStringSet(i)
    nuc.count = nuc.count + oligonucleotideFrequency(dna,width=1)
  }
  nuc.freq  = colSums(nuc.count)/sum(nuc.count)
  cat(sprintf("sum of nuc freq is: %.3f\n",sum(nuc.freq)))
  return(nuc.freq)
}


#-Log-likelihood of subs
LL_min  = function(theta){
  nuc_codon_freq = init_f(cod,DD,64)
  fw  = nuc_codon_freq[[1]]
  cfw = nuc_codon_freq[[2]]
  
  rmat = GTR(theta[1:6], fw)
  Rmat = MG94(rmat,theta[7],cod,codonstrs,syn,64)
  -sum(log(expm(Rmat)*cfw)*DD)
}

#Cal. the number of 16 indicator patterns
grabI = function(iz,pw,npat){
  ip = rep(0,npat)
  for (pat in 1:npat) {
    patm = t(t(iz)==pw[pat,])
    lenI = length(which(rowSums(patm) == ncol(patm)))
    ip[pat] = lenI
  }
  return(ip)
} 


##-Log-likelihood of indels
LL_min_Z = function(gv){
  g012 = gv[1]
  wz   = gv[2]
  #gap LL
  scoreG = sum(phasew*log(g012)) + sum(phasew[,1]*log(wz))
  
  #non-gap LL(must:kk1[1]=0,kk1[16]=1)
  kk1 = rowSums(t(t(pw)*nk1))
  kk2 = rowSums(t(t(pw)*nk2))
  
  #6x16 pat.mat
  # scoreM = 
  #   sum(iw[1:3,]*log(1-g012[1:3] + outer((1-wz)*g012[1:3],kk1))) + 
  #   sum(iw[4:6,]*log(1-g012[4:6] + outer((1-wz)*g012[4:6],kk2)))  
  
  tiw = t(iw)
  scoreM = sum(tiw[,1:3]*log(1-g012 + (1-wz)*g012*kk1)) + 
           sum(tiw[,4:6]*log(1-g012 + (1-wz)*g012*kk2))  
  
  # sum(iw[1,]*log(1-g012[1] + (1-wz)*g012[1]*kk1)) +
  # sum(iw[2,]*log(1-g012[2] + (1-wz)*g012[2]*kk1)) +
  # sum(iw[3,]*log(1-g012[3] + (1-wz)*g012[3]*kk1)) +
  # sum(iw[4,]*log(1-g012[4] + (1-wz)*g012[4]*kk2)) +
  # sum(iw[5,]*log(1-g012[5] + (1-wz)*g012[5]*kk2)) +
  # sum(iw[6,]*log(1-g012[6] + (1-wz)*g012[6]*kk2))
  
  -(scoreG + scoreM)
}






##############################################
#inD  = "k12/Gs/36"
#ouF  = "Results/Gse/36.est.json"
main = function(inD, ouF, omega_l){
  
  sub1 = "../Script/sources/"
  sub2 = "../Script/chapter3/"
  source(paste0(sub1,"codon_call.R"))
  source(paste0(sub1,"gtr.R"))
  source(paste0(sub1,"mg94.R"))
  source(paste0(sub1,"init_f.R"))
  source(paste0(sub1,"LL.R"))
  source(paste0(sub2,"phase_indel_prob_wz.R"))
  
  #construct codons and its degeneracy
  co.res    = codon_call()
  codons    = co.res[[1]]
  codonstrs = co.res[[2]]
  syn       = co.res[[3]]
  
  #nmkb method
  cod      <<- codons
  codonstrs<<- codonstrs
  syn      <<- syn
  
  
  #read input
  Files = list.files(inD, full.names=T)
  index = as.numeric(str_extract(basename(Files),'[^.]+'))
  Files = Files[order(index)]
  K     = length(Files)
  Alist = list()
  for (i in 1:K) {
    Alist[[i]] = readBStringSet(Files[i])
  }
  
  
  #Set up initial parameters
  #indel rates & omega_z
  set.seed(8088)
  omegal = as.numeric(omega_l)  #default: 12
  f0     = count_freq(Files)
  p0     = rep(0.1,7)
  g0     = 0.01
  e0     = rep(0.4,2)
  wz0    = 0.2
  
  rmat   = GTR(p0[1:6],f0)
  t0     = -sum(diag(rmat)*f0)
  
  
  #summary statistic
  #>seq likelihood
  llz  = rep(0,K)                    
  #>number of gaps 
  gavg = matrix(0,K,2)         
  #>gap array
  gArr = array(0,c(6,2,K))       
  #>codon array
  codArr= array(0,c(64,64,K))
  #>Indicator list of Non-gap edge.
  Iz   = list()
  
  #>phase_wz_matrix
  #pw = as.matrix(expand.grid(c(0,1),c(0,1),c(0,1),c(0,1)))  
  pwl = list()
  for(i in 1:(omegal/3)){
    pwl[[i]]=0:1
  }
  pw  = as.matrix(expand.grid(pwl))  
  pw  = unname(pw)
  pw  <<- pw
  npat= 2^(omegal/3)

  
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>Iteration through the black box 
  iter  = 1
  max.it= 20
  repeat{
    if(iter>max.it) {
      cat("Pass the", max.it, "iterations limits!")
      break}
    
    
    avg.gap = 1/(1-e0)
    if(any(avg.gap<1)){
      print("Average gap length less than 3!")
      break}
    
    #prellocate gtr/mg94 mat, codon freq.
    cf0   = sapply(seq(64), function(x){prod(f0[match(cod[x,],DNA_BASES)])})
    Rmat  = MG94(rmat,p0[7],cod,codonstrs,syn)
    Pmat  = log(expm(Rmat))
    
    #geometric gap length (default:0,1,2,3)
    omebin  = seq(0,omegal/3-1)   
    norm.k1 = (1-e0[1])*e0[1]^omebin
    norm.k2 = (1-e0[2])*e0[2]^omebin                    
    #norm.k1 = lk1/sum(lk1)
    #norm.k2 = lk2/sum(lk2)
    print(c(sum(norm.k1),sum(norm.k2)))
    
    #profvis({
      for(i in 1:K){#E step
        A = Alist[[i]]
        res1        = ziqi_prob_wz(A,g0,e0,Pmat,codonstrs,syn,f0,wz0,omegal,norm.k1,norm.k2)
        llz[i]      = res1[[1]]
        gArr[,,i]   = res1[[2]]
        gavg[i,]    = res1[[3]]
        Iz[[i]]     = res1[[4]]
        codArr[,,i] = res1[[5]]
        print(i)
      }
    #})
    
    
    ####################################
    #M-step: parameter estimates
    #pseudo weight
    Wi = rep(1/K,K)
    
    #nmkb for p0
    #>weighted codon matrix
    datw = matrix(0,64,64)
    for (j in 1:K) {
      dat.wei = Wi[j]*codArr[,,j]
      datw    = datw + dat.wei
    }
    DD <<- datw
    pb = nmkb(fn=LL_min, par=p0, lower=0.0, control=list(tol=1e-4,trace=F)) 
    if(pb$convergence != 0){
      cat("Warning: failed convergence!")
    }else{
      pnew = pb$par
    }
    
    ##cal. the tau
    nuc_f = init_f(cod,DD,64)
    fnew  = nuc_f[[1]]
    rmat  = GTR(p0[1:6],fnew)
    tnew  = -sum(diag(rmat)*fnew)
    
    
    ###############################
    #gap extension 
    w.avg.gap = colMeans(gavg)
    #w.avg.gap = c(mean(gavg[which(gavg[,1]>0),1]),mean(gavg[which(gavg[,2]>0),2]))
    enew      = 1-1/w.avg.gap
    
    #weighted gap phases
    phasew = matrix(0,6,2)
    for (j in 1:K) {
      phasew = phasew + gArr[,,j]*Wi[j]
    }
    
    
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    iw = matrix(0,6,npat) #based on pw patterns.
    
    #profvis({
      for (h in 1:K) {
        ipos1 = seq(1,nrow(Iz[[h]][,,1])) %% 3
        ipos2 = seq(1,nrow(Iz[[h]][,,2])) %% 3
        
        #insertion of phase 0,1,2
        Iz0 = Iz[[h]][which(ipos1==1),,1]
        Iz1 = Iz[[h]][which(ipos1==2),,1]
        Iz2 = Iz[[h]][which(ipos1==0),,1]
        
        #Deletion of phase 0,1,2
        Iz0D = Iz[[h]][which(ipos2==1),,2]
        Iz1D = Iz[[h]][which(ipos2==2),,2]
        Iz2D = Iz[[h]][which(ipos2==0),,2]
        
        iw[1,] = iw[1,] + grabI(Iz0,pw,npat)/K
        iw[2,] = iw[2,] + grabI(Iz1,pw,npat)/K
        iw[3,] = iw[3,] + grabI(Iz2,pw,npat)/K
        iw[4,] = iw[4,] + grabI(Iz0D,pw,npat)/K
        iw[5,] = iw[5,] + grabI(Iz1D,pw,npat)/K
        iw[6,] = iw[6,] + grabI(Iz2D,pw,npat)/K
      }
   #})
    
     
    #nmkb for g0
    phasew<<- phasew
    iw    <<- iw
    
    ############################
    # lk_geo = function(ext,omebin){
    #   lk1     = (1-ext[1])*ext[1]^omebin
    #   lk2     = (1-ext[2])*ext[2]^omebin
    #   return(rbind(lk1,lk2))
    # }
    #
    # norm.k1 = lk_geo(enew,omebin)[1,]
    # norm.k2 = lk_geo(enew,omebin)[2,]
    
    # nk1   <<- lk_geo(enew,omebin)[1,]
    # nk2   <<- lk_geo(enew,omebin)[2,]
    
    
    nk1   <<- norm.k1
    nk2   <<- norm.k2
    ###################################
    
    gv  = c(g0,wz0)
    gb  = nmkb(fn=LL_min_Z, par=gv, lower=0.0, upper=1.0, control=list(tol=1e-4,trace=T)) 
    if(gb$convergence != 0){
      cat("Warning: failed convergence!")
    }else{
      gv    = gb$par
      gnew  = gv[1]
      wznew = gv[2]
    }
    
    #print output
    print(pnew[1:6]/tnew)
    print(pnew[7])
    print(tnew)
    print(enew)
    print(-log(1-2*gnew)/2)
    print(wznew)
    
    
    #rmse tolerance
    p    = c(g0,e0,p0,wz0,t0)
    q    = c(gnew,enew,pnew,wznew,tnew)
    delta= (q-p)^2
    rmse = sqrt(mean(delta))
    cat(sprintf("iter:%i, rmse:%.6f\n",iter, rmse))
    if(rmse<=1e-4){
      break
    }else{
      iter=iter+1
      f0  = fnew
      p0  = pnew
      g0  = gnew
      e0  = enew
      wz0 = wznew
      t0  = tnew
    }
  }
  
  
##output Jsons
par_lst = list('nuc.freq'=fnew, 'sigmas'=pnew[1:6]/tnew, 'omega'=pnew[7], 'branch.length'=tnew, 
               'gap.openning'=gnew,'gap.extension'=enew, 'omega.z'=wznew)

sum_lst = list('zn.phases'=phasew[,1],'zs.phases'=phasew[,2],'avg.gap.size'=w.avg.gap)

par_est = toJSON(par_lst)
sum_stat= toJSON(sum_lst) 

dname = dirname(ouF)
fname = str_extract(basename(ouF),'[^.]+')
ouF2  = paste0(dname,'/',fname,'.sum.json')
write(par_est,ouF)
write(sum_stat,ouF2)
  
  
  
}

########################################
args = commandArgs(trailingOnly=T)
main(args[1],args[2],args[3])

