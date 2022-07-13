# Take the coatiM (viterbi algo) and weight score as input
# Apply the importance sampling Wi = g/f
# g(x)--ziqi's phase_coati model; f(x)--juan's coatiM
#
# A -- ancestor, B -- descendent
# insertion and deletion length distribution combined.
# inserting bases are ignored since they are treated as deletions.

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
  if(sum(nuc.freq)==1){
    return(nuc.freq)
  }else{
    print("The nucleotide frequency sums up not 1!")
  }
}

#Read sample
read_sample = function(D,K){
  Alist = list()
  for (i in 1:K) {
    Alist[[i]] = DNAStringSet(c(D$Seq1[i],D$Seq2[i]))
  }
  return(Alist)
}

#-Log-likelihood of
LL_min  = function(theta){
  nuc_codon_freq = init_f(cod,DD,64)
  fw  = nuc_codon_freq[[1]]
  cfw = nuc_codon_freq[[2]]
  
  rmat = GTR(theta[1:6], fw)
  Rmat = MG94(rmat,theta[7],cod,codonstrs,syn)
  -sum(log(expm(Rmat)*cfw)*DD)
}

#remove effect of scaling factor [ancestor/descendant]
test_pab = function(ab,f0){
  dnaAB    = DNAStringSet(gsub('-','',ab))
  ab.count = oligonucleotideFrequency(dnaAB,width=1)
  sumPab   = sum(ab.count[1,]*log(f0)) + sum(ab.count[2,]*log(f0))
  return(sumPab)
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


#Generate the 6x16 indicator-phase matrix
id_phase_m = function(Izi,Iposi,wii,pw,npat){
  iwi = matrix(0,6,npat) 
  
  ipos1 = which(Iposi[,1]!=0) %%3
  ipos2 = which(Iposi[,2]!=0) %%3
  
  #insertion of phase 0,1,2
  Iz0 = Izi[[1]][which(ipos1==1),]
  Iz1 = Izi[[1]][which(ipos1==2),]
  Iz2 = Izi[[1]][which(ipos1==0),]
  
  #Deletion of phase 0,1,2
  Iz0D = Izi[[2]][which(ipos2==1),]
  Iz1D = Izi[[2]][which(ipos2==2),]
  Iz2D = Izi[[2]][which(ipos2==0),]
  
  iwi[1,] =  grabI(Iz0,pw,npat) 
  iwi[2,] =  grabI(Iz1,pw,npat) 
  iwi[3,] =  grabI(Iz2,pw,npat) 
  iwi[4,] =  grabI(Iz0D,pw,npat)
  iwi[5,] =  grabI(Iz1D,pw,npat)
  iwi[6,] =  grabI(Iz2D,pw,npat)
  
  return(iwi*wii)
}

##-Log-likelihood of indels
LL_min_Z = function(gv){
  g012 = gv[1]
  wz   = gv[2]
  #gap LL
  scoreG = sum(Nw*log(g012)) + sum(Nw[,1]*log(wz))
  
  #non-gap LL
  kk1 = rowSums(t(t(pw)*nk1))
  kk2 = rowSums(t(t(pw)*nk2))
  
  #6x16 pat.mat
  tIw = t(Iw)
  
  # scoreM = 
  #   sum(Iw[1:3,]*log(1-g012[1:3] + outer((1-wz)*g012[1:3],kk1))) + 
  #   sum(Iw[4:6,]*log(1-g012[4:6] + outer((1-wz)*g012[4:6],kk2)))  
  
  scoreM = 
    sum(tIw[,1:3]*log(1-g012 + (1-wz)*g012*kk1)) + 
    sum(tIw[,4:6]*log(1-g012 + (1-wz)*g012*kk2)) 
  
  -(scoreG + scoreM)
}

#generate next lk
lk_geo = function(ext,omebin){
  lk1     = (1-ext[1])*ext[1]^omebin
  lk2     = (1-ext[2])*ext[2]^omebin
  return(rbind(lk1,lk2))
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>test

#Judge whether non-gap state would add the failed zn/zs
#Only use deletion method
znzs_treat = function(rs, rs1, L2k, i, syn, wid){
  
  rem = i %% 3
  #indicator of each position
  l     = length(wid)
  idcat = rep(0,l)   #syn
  
  if(rem==2){#phase1
    sub1 = substring(rs,i-1,i+1)
    for (j in seq(l)) {
      if((i+wid[j])>L2k){
        next
      }
      sub2= substring(rs,i+wid[j]-1,i+wid[j]+1)
      ref = paste0(rs1[i-1],rs1[i+wid[j]],rs1[i+wid[j]+1],collapse='')
      sec = syn[[ref]]
      if(all(!(c(sub1,sub2) %in% sec))){
        idcat[j]=1
      }
    }
  }else if(rem==0){#phase2
    sub1 = substring(rs,i-2,i)
    for (j in seq(l)) {
      if((i+wid[j])>L2k){
        next
      }
      sub2= substring(rs,i+wid[j]-2,i+wid[j])
      ref = paste0(rs1[i-2],rs1[i-1],rs1[i+wid[j]],collapse='')
      sec = syn[[ref]]
      if(all(!(c(sub1,sub2) %in% sec))){
        idcat[j]=1
      }
    }
  }else{#phase0
    return(idcat)
  }
  
  return(idcat)
}

#Generate the indicator matrix of non-gap edge of every position for two seqs separately
#non-syn:1 || syn:0
nogap_sum = function(rA, syn, omegal){
  rAs = str_split(rA,'')
  L2  = width(rA)     
  
  idm1 = matrix(NA,L2[1],omegal/3)
  idm2 = matrix(NA,L2[2],omegal/3)
  idm  = list(idm1,idm2)
  wid  = seq(3,omegal,3)
  
  for (k in seq(2)) {
    rs1 = rAs[[k]]  
    rs  = paste0(rAs[[k]],collapse="") 
    
    i=1 
    while (i <= L2[k]) {
      idm[[k]][i,] = znzs_treat(rs,rs1,L2[k],i,syn,wid) 
      i = i+1
    }
  }
  
  return(idm)
}

#non-gap prob
score_M  = function(M,lk1,lk2,g0,wz0){
  n1 = dim(M[[1]])[1]
  n2 = dim(M[[2]])[1]
  m  = dim(M[[1]])[2]
    
  if(n1%%3 ==1){
    M[[1]] = rbind(M[[1]],matrix(NA,2,m))
  }else if(n1%%3 == 2){
    M[[1]] = rbind(M[[1]],rep(NA,m))
  }
  
  if(n2%%3 ==1){
    M[[2]] = rbind(M[[2]],matrix(NA,2,m))
  }else if(n2%%3 == 2){
    M[[2]] = rbind(M[[2]],rep(NA,m))
  }
  
  k1     = t(t(M[[1]])*lk1)
  k2     = t(t(M[[2]])*lk2)
  scoreM = sum(log(1-g0 + (1-wz0)*g0*rowSums(k1)), na.rm=TRUE) +
           sum(log(1-g0 + (1-wz0)*g0*rowSums(k2)), na.rm=TRUE) 
  return(scoreM)
}



##############################################
# inD   = "k12/Gs_trim/1"
# ouD   = "k12/JsonD/1/"

main = function(inD,ouD,ouF,ssize,omega_l){
  
  sub1 = "../Script/sources/"
  sub2 = "../Script/chapter3/"
  source(paste0(sub1,"codon_call.R"))
  source(paste0(sub1,"gtr.R"))
  source(paste0(sub1,"mg94.R"))
  source(paste0(sub1,"init_f.R"))
  source(paste0(sub1,"LL.R"))
  source(paste0(sub2,"phase_indel_prob_wz_1.R"))
  
  #construct codons and its degeneracy
  co.res    = codon_call()
  codons    = co.res[[1]]
  codonstrs = co.res[[2]]
  syn       = co.res[[3]]
  
  #nmkb method
  cod      <<- codons
  codonstrs<<- codonstrs
  syn      <<- syn
  
  
  ####################################
  Files = list.files(inD, full.names=T)
  namev = str_extract(basename(Files),'[^.]+')
  
  if(namev[1]=='1'){#simulation
    index = as.numeric(namev) #rectify the order
    Files = Files[order(index)]
    tag   = as.numeric(basename(inD))
  }else{#90 species
    Fbname = list.files(inD, full.names=F)
    kept   = read.table(spList,header=F)[,1]
    index  = which(Fbname %in% kept)
    Files  = Files[index]
  }
  n = length(Files)
  
  #Initial parameters
  set.seed(8088)
  f0  = count_freq(Files)
  p0  = rep(0.1,7)
  
  g0  = 0.01
  e0  = rep(0.8,2)               #nuc extension
  wz0 = 0.2
  
  rmat= GTR(p0[1:6],f0)
  t0  = -sum(diag(rmat)*f0) 
  
  #Iterate through the black box
  max.it = 20
  K      = as.numeric(ssize)     #default:ssize=100
  
  #LL
  lljv   = matrix(0,n,max.it)
  llzv   = matrix(0,n,max.it)
  #avg gap len
  E           = matrix(0,n,2)
  colnames(E) = c('avg.gap.len.I','avg.gap.len.D')
  #weight 
  Wv = list()          
  
  
  #Prepare pattern matrix[16x4], number of patterns.
  omegal = as.numeric(omega_l)   #default:omega_l=12
  omebin = seq(0,omegal/3-1)  
  pwl = list()
  for(i in 1:(omegal/3)){
    pwl[[i]]=0:1
  }
  pw  = as.matrix(expand.grid(pwl))  
  pw  = unname(pw)
  pw  <<- pw  
  npat= 2^(omegal/3)
  
  
  
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>RUNNING
  iter=1
  tm = system.time({
    repeat{
      #Reset those summary stat in each iteration.
      Dw  = matrix(0,64,64)
      Iw  = matrix(0,6,npat)
      Nw  = matrix(0,6,2)
      
      if(iter>max.it){
        cat("Pass the", max.it, "iterations limits!")
        break}
      
      #pre-cal avg gap size
      avg.gap = 1/(1-e0)
      if(any(avg.gap<3)){
        print("Average gap length less than 3!")
        break
      }
      e3 = 1-1/(avg.gap/3)
      
      #prellocate gtr/mg94 mat, codon freq.
      cf0   = sapply(seq(64), function(x){prod(f0[match(cod[x,],DNA_BASES)])})
      Rmat  = MG94(rmat,p0[7],cod,codonstrs,syn)
      Pmat  = log(expm(Rmat))
      
      #geometric gap length (default:0,1,2,3): unnormalized
      lk1 = (1-e3[1])*e3[1]^omebin
      lk2 = (1-e3[2])*e3[2]^omebin 
      
      #tm1 =system.time({
        for (j in 1:n) {#E-step
          ##coati-sampler
          input = Files[j]
          ouJ   = paste0(ouD,j,'.json')
          cmd   = paste("bash ../Script/chapter3/coati_sampler.sh", input, ouJ, K,
                        g0,mean(e0),f0[1],f0[2],f0[3],f0[4],p0[1],p0[2],p0[3],p0[4],p0[5],p0[6],p0[7],t0,sep=' ')
          system(cmd)
          Input.tmp    = fromJSON(ouJ)
          Dat.tmp      = Input.tmp %>% tidyr::unpack(aln)
          
          aln     = Dat.tmp %>% dplyr::group_by(Seq1,Seq2)  
          aln2    = aln %>% dplyr::group_keys(Seq1,Seq2,weight,log_weight)
          
          gpsize  = group_size(aln)
          ngroups = n_groups(aln)
          print(ngroups)
          
          llj          = aln2$log_weight    #unnormalized
          lljv[j,iter] = sum(llj*gpsize)/K
          Alist        = read_sample(aln2,ngroups)
          
          #llj2 = exp(llj)/sum(exp(llj)*gpsize)
          #print(sum(llj2*gpsize))
          
          
          ##preset stats for each alignment
          Gz        = array(0,c(6,2,ngroups))
          codon_arr = array(0,c(64,64,ngroups))
          gl        = matrix(0,ngroups,2)
          llz       = rep(0,ngroups)
          Ipos      = list()
          Iz        = list()
          
          ##adding scalling factor
          scal.pab = test_pab(Alist[[1]],f0)
          
          ##pre-calculate the indicator M
          rA     = DNAStringSet(gsub('-','',Alist[[1]]))
          M.012  = nogap_sum(rA, syn, omegal)
          
    
          #profvis({
          for(i in 1:ngroups){
            A    = Alist[[i]]
            res1 = ziqi_prob_wz(A,g0,e3,Pmat,codonstrs,syn,f0,wz0,omegal,lk1,lk2)
            llz[i]         = res1[[1]] - scal.pab 
            Gz[,,i]        = res1[[2]]
            gl[i,]         = res1[[3]]
            codon_arr[,,i] = res1[[4]]
            Ipos[[i]]      = res1[[5]]
          } 
          #})
          
          #correct the indicator matrix, LL, indicator position.  
          for(i in 1:ngroups){
            Iz1 = M.012[[1]][Ipos[[i]][,1],]
            Iz2 = M.012[[2]][Ipos[[i]][,2],]
            Iz[[i]] = list(Iz1,Iz2) 
            scoreM  = score_M(Iz[[i]],lk1,lk2,g0,wz0)
            llz[i]  = llz[i] + scoreM
          }
          
          
          ############
          #cal. the weight
          llzv[j,iter] = sum(llz*gpsize)/K
          
          uwi          = llz-llj
          uwi.max      = max(llz-llj)
          uwi.dif      = uwi - uwi.max
          Wi           = exp(uwi.dif)/sum(exp(uwi.dif)*gpsize)
          print(sum(Wi*gpsize))
          Wv[[j]]      = Wi*gpsize
          
          
          ##weighted average gap length
          E[j,] = colSums(Wi*gl*gpsize) 
          
          ##weighted gap phases, non-gap indicators, codon matrix 
          N    = matrix(0,6,2)
          iw   = matrix(0,6,npat) 
          datw = matrix(0,64,64)    
          
          for (i in 1:ngroups) {
            #>
            N.wei = Wi[i]*Gz[,,i]
            N     = N + N.wei*gpsize[i]
            #>
            iw.wei= id_phase_m(Iz[[i]],Ipos[[i]],Wi[i],pw,npat) 
            iw    = iw + iw.wei*gpsize[i]
            #>
            dat.wei = Wi[i]*codon_arr[,,i]
            datw    = datw + dat.wei*gpsize[i]
          }
          
          Nw = Nw + N/n
          Iw = Iw + iw/n
          Dw = Dw + datw/n
        }
      #})
      
      
      
      
      
      #######################
      #M step: para estimates
      
      ##gap extension
      #w.avg.gap = c(mean(E[which(E[,1]!=0),1]),mean(E[which(E[,2]!=0),2]))
      w.avg.gap = c(mean(E[which(E[,1]>=1),1]),mean(E[which(E[,2]>=1),2]))
      enew1     = 1 - 1/(w.avg.gap*3)
      enew      = 1 - 1/w.avg.gap
      if(any(enew<0)){
        print("Warning: gap extension prob (unit of 3) < 0!")
        break
      }
      
      
      ##gap opening 
      # nnew = colMeans(N)      
      # mnew = colMeans(M)
      # gnew = nnew/(nnew+mnew)
      
      #nmkb for p0
      DD <<- Dw
      pb = nmkb(fn=LL_min, par=p0, lower=0, upper=1, control=list(tol=1e-5,trace=T,maxfeval=5000))  #change the tolerance
      if(pb$convergence != 0){
        cat("Warning: failed convergence!")
      }else{
        pnew = pb$par
      }
      
      ##cal. the tau
      nuc_f = init_f(cod,DD,64)     
      fnew  = nuc_f[[1]]
      rmat  = GTR(pnew[1:6],fnew)
      tnew  = -sum(diag(rmat)*fnew)
      
      ##nmkb for g0
      Nw <<- Nw
      Iw <<- Iw
      
      nk1   <<- lk_geo(enew,omebin)[1,]
      nk2   <<- lk_geo(enew,omebin)[2,]
      
      gv  = c(g0,wz0)
      gb  = nmkb(fn=LL_min_Z, par=gv, lower=0.0, upper=1.0, control=list(tol=1e-4,trace=T)) 
      if(gb$convergence != 0){
        cat("Warning: failed convergence!")
      }else{
        gv    = gb$par
        gnew  = gv[1]
        wznew = gv[2]
      }
      
      
      #print the results
      print(pnew[1:6]/tnew)
      print(pnew[7])
      print(tnew)
      print(enew)
      print(-log(1-2*gnew)/2)
      print(wznew)
      
      #record the pars
      # gv   = rbind(gv,gnew)
      # pv   = rbind(pv,pnew)
      # ev   = c(ev,enew)
      # tv   = c(tv,tnew)
      
      #rmse tolerance
      p    = c(g0,e0,p0,wz0,t0)
      q    = c(gnew,enew1,pnew,wznew,tnew)
      delta= (q-p)^2
      rmse = sqrt(mean(delta))
      cat(sprintf("iter:%i, rmse:%.6f\n",iter, rmse))
      if(rmse<=1e-4){#adjusted
        break
      }else{
        if(iter>10){
          #K=1e+3
          # if(tag==98){
          #   max.it=11
          # }
        }
        iter= iter+1
        f0  = fnew
        p0  = pnew
        g0  = gnew
        e0  = enew1
        wz0 = wznew
        t0  = tnew
      }
    }
  })
  
  cat(sprintf("Running loops: %i\n  Running time:%.3f mins", iter, tm[3]/60))
  #########################################################
  #plot ll
  # lljv = lljv[,which(colSums(lljv)!=0)]
  # llzv = llzv[,which(colSums(llzv)!=0)]
  # maxw = maxw[,which(colSums(maxw)!=0)]
  # 
  # set.seed(8088)
  # colorB  = c("#000000", "#E69F00", "#56B4E9", "#009E73",
  #             "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # x = 1:iter
  # plot(x,lljv[1,],type='o',col=colorB[1],pch='*',ylim=c(min(lljv),max(lljv)))
  # for (i in 2:n) {
  #   y = lljv[i,]
  #   cr= colorB[sample(1:8,1)]
  #   points(x,y,col=cr,pch='*')
  #   lines(x,y,col=cr,lty=1)
  # }
  # 
  # plot(x,llzv[1,],type='o',col=colorB[1],pch='*',ylim=c(min(llzv),max(llzv)))
  # for (i in 2:n) {
  #   y = llzv[i,]
  #   cr= colorB[sample(1:8,1)]
  #   points(x,y,col=cr,pch='*')
  #   lines(x,y,col=cr,lty=1)
  # }
  
  # #plot maximum weight
  # plot(x,maxw[1,],type='o',col=colorB[1],pch='*',ylim=c(min(maxw),max(maxw)))
  # for (i in 2:n) {
  #   y = maxw[i,]
  #   cr= colorB[sample(1:8,1)]
  #   points(x,y,col=cr,pch='*')
  #   lines(x,y,col=cr,lty=1)
  # }
  
  
  
  #Cal. the extension prob (T-G)/T
  #will the total length be weighted as well?
  
  ##pipe out the json file.
  ##confirm with gillespie simu:
  ##p0[1:6]/tnew vs Sigma/T
  ##p0[7] vs W
  ##Pi
  ##ext0 vs ext
  ##tnew vs brlen
  ##r1*brlen vs (g0[1]+g0[4])/2, (g0[2]+g0[5])/2, (g0[3]+g0[6])/2
  
  # ouF1 = "pise.json"
  # ouF2 = "weight_mat.txt" 
  
  par_lst = list('nuc.freq'=fnew, 'sigmas'=pnew[1:6]/tnew, 'omega'=pnew[7], 'branch.length'=tnew, 
                 'gap.openning'=gnew,'gap.extension'=enew,'omega.z'=wznew)
  sum_lst = list('zn.phases'=Nw[,1],'zs.phases'=Nw[,2],'avg.gap.size'=w.avg.gap)
  
  par_est = toJSON(par_lst)
  sum_stat= toJSON(sum_lst)
  #output the parameter estimates, summary stats, weight matrix
  #ouF = "Results/Pise/1.est.json"
  dname = dirname(ouF)
  fname = str_extract(basename(ouF),'[^.]+')
  ouF2  = paste0(dname,'/',fname,'.sum.json')
  write(par_est,ouF)
  write(sum_stat,ouF2)
  
  #dname2 = "Results/weightM"
  #ouF3   = paste0(dname2,'/',fname,'.txt')
  #lapply(Wv,write,ouF3,append=T,ncolumns=10000)
  
  
  
}


################################
args = commandArgs(trailingOnly=T)
main(args[1],args[2],args[3],args[4],args[5])
