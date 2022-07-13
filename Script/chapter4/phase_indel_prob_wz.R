#1: output the prob of ziqi's model
#2: output the summuary stats

#[gap.length>12] are unlikely to preserve the originial aa. treat 

#Cal. the inserted freq prob
prob_insert = function(As2,g1,f0){
  if(length(g1)==0){
    return (0)
  }else{
    pos1 = start(g1)
    pos2 = end(g1)
    insN = c()
    for(j in 1:length(g1)){
      insN = c(insN,As2[pos1[j]:pos2[j]])
    }
  }
  f.rank = match(insN,DNA_BASES)
  f.prob = log(prod(f0[f.rank]))
  return(f.prob)
}

#Cal. the prob of the end state
end_prob = function(lenA,g,g0,e0){
  if(lenA %in% end(g[[1]])){#ins>end
    endP = log(1-e0[1])
  }else if(lenA %in% end(g[[2]])){#del>end
    endP = log(1)
  }else{#match>end
    endP = log(1-g0)
  }
  return(endP)
}



#del-->subs-->ins.
#Judge whether the gap state is zn/zs. 
#record the one-after-gap-right-end pos.
znzs_cure = function(seqc, gz, k, syn){
  if(k==1){#ins 
    s =seqc[[2]]
    s3=seqc[[1]]
  }else{#del
    s =seqc[[1]]
    s3=seqc[[2]]
  }
  
  imat = matrix(0,3,2)
  ps1  = start(gz)
  ps2  = end(gz)
  
  rem = ps1 %% 3
  pos = c()      #record the 0,1-edge phase position
  
  #1A-- -AA AAA  2A-- -AA AAA   4AAA AA- --A  3AAA AAA AAA  >phase1
  #1AAA A-- -AA  2AAA AA- --A   4A-- -AA AAA  3A-- -A- --A
  
  #5AA- --A AAA   >phase2
  #5AAA AA- --A  
  
  #6--- NNN NNN  7--- NNN  8NNN NNN  9NNN --- end  10NNN NNN end >phase0 
  #6NNN --- NNN  7NNN NNN  8--- NNN  9NNN NNN end  10NNN --- end
  
  # for (i in 1:length(rem)){
  #   if(rem[i]==2){#phase1
  #     
  #   }else if(rem[i]==0){#phase2
  #     
  #   }else{#phase0
  #     
  #   }
  # }
  
  
  for (i in 1:length(rem)) {
    if(rem[i] == 2){#phase 1 
      if(s3[ps2[i]+2]=='-'){#3
        tmp.0=s[ps1[i]-1]
        tmp.1=s[ps2[i]+1]
        j = ps2[i]+2
        repeat{
          j     = j+1
          if(s3[j]!='-'){
            break
          }
        }
        tmp.2  = s[j]
        pos    = c(pos,ps2[i]+1)
        ref    = paste0(tmp.0,tmp.1,tmp.2,          collapse="")
        sub1   = paste0(tmp.0,s[ps1[i]],s[ps1[i]+1],collapse="")
        sub2   = paste0(s[ps2[i]],tmp.1,s[ps2[i]+2],collapse="")
        sub3   = paste0(s[j-2],s[j-1],tmp.2        ,collapse="")
        sec    = syn[[ref]]
        if(sub1 %in% sec || sub2 %in% sec || sub3 %in% sec){ 
          imat[2,2]=imat[2,2]+1 
        }else{
          imat[2,1]=imat[2,1]+1 
        }
      }else{
        if(s[ps2[i]+1]=='-'){#1
          j = ps2[i]+1
          repeat{
            j      = j+1
            tmp.1  = s[j]
            if(tmp.1!='-'){
              break
            }
          }
          tmp.0=s[ps1[i]-1]
          tmp.2=s[j+1]
          pos  = c(pos,j)
        }else if(s[ps2[i]+2]=='-'){#2,4
          tmp.0=s[ps1[i]-1]
          tmp.1=s[ps2[i]+1]
          j = ps2[i]+2
          repeat{
            j     = j+1
            tmp.2 = s[j]
            if(tmp.2!='-'){
              break
            }
          }
          pos  = c(pos,ps2[i]+1)
        }else if(s[ps1[i]-1]=='-'){#1
          j = ps1[i]-1
          repeat{
            j      = j-1
            tmp.0  = s[j]
            if(tmp.0!='-'){
              break
            }
          }
          tmp.1 = s[ps2[i]+1]
          tmp.2 = s[ps2[i]+2]
          pos   = c(pos,ps2[i]+1)
        }else{#normal
          tmp.0 = s[ps1[i]-1]
          tmp.1 = s[ps2[i]+1]
          tmp.2 = s[ps2[i]+2]
          pos   = c(pos,ps2[i]+1)
        }
        ref    = paste0(tmp.0,tmp.1,tmp.2,          collapse="")
        sub1   = paste0(tmp.0,s[ps1[i]],s[ps1[i]+1],collapse="")
        sub2   = paste0(s[ps2[i]],tmp.1,tmp.2,      collapse="")
        sec    = syn[[ref]]
        if(sub1 %in% sec || sub2 %in% sec){ 
          imat[2,2]=imat[2,2]+1 
        }else{
          imat[2,1]=imat[2,1]+1 
        }
      }
    }else if(rem[i] == 0){#phase2
      if(s3[ps1[i]-2]=='-'){#3
        j = ps1[i]-2
        repeat{
          j = j-1
          if(s3[j]!='-'){
            break
          }
        }
        tmp.0  = s[j]
        tmp.1  = s[ps1[i]-1]
        tmp.2  = s[ps2[i]+1]
        pos    = c(pos,ps2[i]+1)
        ref    = paste0(tmp.0,tmp.1,tmp.2,          collapse="")
        sub1   = paste0(tmp.0,s[j+1],s[j+2],        collapse="")
        sub2   = paste0(s[ps1[i]-2],tmp.1,s[ps1[i]],collapse="")
        sub3   = paste0(s[ps2[i]-1],s[ps2[i]],tmp.2,collapse="")
        sec    = syn[[ref]]
        if(sub1 %in% sec || sub2 %in% sec || sub3 %in% sec){ 
          imat[3,2]=imat[3,2]+1 
        }else{
          imat[3,1]=imat[3,1]+1 
        }
      }else{
        if(s[ps1[i]-2]=='-'){#2,4
          j = ps1[i]-2
          repeat{
            j      = j-1
            tmp.0  = s[j]
            if(tmp.0!='-'){
              break
            }
          }
          tmp.1 = s[ps1[i]-1]
          tmp.2 = s[ps2[i]+1]
          pos   = c(pos,ps2[i]+1)
        }else if(s[ps2[i]+1]=='-'){#5
          j = ps2[i]+1
          repeat{
            j      = j+1
            tmp.2  = s[j]
            if(tmp.2!='-'){
              break
            }
          }
          tmp.0=s[ps1[i]-2]
          tmp.1=s[ps1[i]-1]
          pos  = c(pos,j)
        }else if(s[ps1[i]-1]=='-'){#5
          j = ps1[i]-1
          repeat{
            j      = j-1
            tmp.1  = s[j]
            if(tmp.1!='-'){
              break
            }
          }
          tmp.0 = s[j-1]
          tmp.1 = s[j]
          tmp.2 = s[ps2[i]+1]
          pos   = c(pos,ps2[i]+1)
        }else{
          tmp.0 = s[ps1[i]-2]
          tmp.1 = s[ps1[i]-1]
          tmp.2 = s[ps2[i]+1]
          pos   = c(pos,ps2[i]+1)
        }
        ref    = paste0(tmp.0,tmp.1,tmp.2,          collapse="")
        sub1   = paste0(tmp.0,tmp.1,s[ps1[i]],      collapse="")
        sub2   = paste0(s[ps2[i]-1],s[ps2[i]],tmp.2,collapse="")
        sec    = syn[[ref]]
        if(sub1 %in% sec || sub2 %in% sec){ 
          imat[3,2]=imat[3,2]+1  
        }else{
          imat[3,1]=imat[3,1]+1  
        }
      }
    }else{#phase0
      if(is.na(s[ps2[i]+1])){#9,10
        pos = 0
      }else{
        if(s[ps2[i]+1]=='-'){#6 
          j = ps2[i]+1
          repeat{
            j    = j+1
            tmp  = s[j]
            if(tmp!='-'){
              break
            }
          }
          pos =c(pos,j)
        }else{#7,8
          pos =c(pos,ps2[i]+1) 
        }
      }
      imat[1,2]=imat[1,2]+1   
    }
  }
  res = list(imat,pos)
  return(res)
}


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
      if((i+wid[j])>(L2k+1)){
        next
      }
      sub2= substring(rs,i+wid[j]-1,i+wid[j]+1)
      ref = paste0(rs1[i-1],rs1[i+wid[j]],rs1[i+wid[j]+1],collapse='')
      sec = syn[[ref]]
      if(!(sub1 %in% sec) && !(sub2 %in% sec)){
        idcat[j]=1
      }
    }
  }else if(rem==0){#phase2
    sub1 = substring(rs,i-2,i)
    for (j in seq(l)) {
      if((i+wid[j])>(L2k+1)){
        next
      }
      sub2= substring(rs,i+wid[j]-2,i+wid[j])
      ref = paste0(rs1[i-2],rs1[i-1],rs1[i+wid[j]],collapse='')
      sec = syn[[ref]]
      if(!(sub1 %in% sec) && !(sub2 %in% sec)){
        idcat[j]=1
      }
    }
  }else{#phase0
    return(idcat)
  }
  return(idcat)
}

#Generate the phase matrix [zn,zs,long]
#Generate a list of non-gap-site-indicator-removed positions
#Count the number of indels of {3,6,9,12}
#Sum the length of indels
gap_sum = function(As,g,syn,omegal){
  Imat   = list()          
  sphase = matrix(0,6,2)   #short gap phase         
  lphase = matrix(0,3,2)   #long gap phase
  posEdt = list('ins.s'=0, 'del.s' =0,
                'ins.l'=0, 'del.l' =0)  
  
  #gap number; gap length
  lG = c(length(g[[1]]),length(g[[2]]))
  lT = c(sum(width(g[[1]])),sum(width(g[[2]])) )/3
  
  #short gap (znzs); long gap
  if(length(unlist(g))==0){#no gaps
    return(list(sphase,posEdt,c(0,0),c(0,0)))
  }else{
    for (k in seq(2)) {
      wid = width(g[[k]])
      if(length(wid)==0){next}
      
      if(any(wid>omegal)){
        id  = which(wid>omegal)
        pos = start(g[[k]][id])
        stp = end(g[[k]][id])
        rem = pos %% 3
        
        lphase[1,k]   = length(which(rem == 1))
        lphase[2,k]   = length(which(rem == 2))
        lphase[3,k]   = length(which(rem == 0))
        
        if(length(which(stp!=length(As[[1]]))) != 0){
          id2           = which(stp!=length(As[[1]]))
          posEdt[[k+2]] = stp[id2]+1
        }
      }
      
      gz = g[[k]][which(wid<=omegal)] 
      if(length(gz)==0){
        next
      }else{
        Imat[[k]] = znzs_cure(As,gz,k,syn) 
        if(k==1){
          k0 = k
        }else{
          k0 = k+2
        }
        sphase[k0:(k0+2),] = Imat[[k]][[1]]   
        posEdt[[k]]        = Imat[[k]][[2]]
      }
    }
  }
  
  sphase[,2]       = sphase[,2]+c(lphase)
  colnames(sphase) = c('zn','zs')
  rownames(sphase) = c('I0','I1','I2','D0','D1','D2')
  
  res  = list(sphase,posEdt,lG,lT)
  return(res)
}



#########################

#Generate the indicator matrix of non-gap edge of every position for two seqs separately
#non-syn:1 || syn:0
nogap_sum = function(As, rAs, syn, omegal){
  L1 = length(As[[1]])                              #align len
  L2 = c(length(rAs[[1]]),length(rAs[[2]]))         #each len
  
  idm = array(NA,c(L1,omegal/3,2))
  wid = seq(3,omegal,3)
  
  for (k in seq(2)) {
    s1 = As[[k]]
    rs = paste0(rAs[[k]],collapse="")
    rs1= rAs[[k]]
    
    i=j=1
    while (i<=L2[k]) {
      if(s1[j]=='-'){
        repeat{
          j=j+1
          if(s1[j]!='-'){
            break
          }
        }
      }
      if((i+2)<=L2[k]){#not in the last codon
        idm[j,,k] = znzs_treat(rs,rs1,L2[k],i,syn,wid)
      }
      i=i+1
      j=j+1
    }
  }
  
  #check
  if((L1==j-1) && (L2[2]==i-1)){
    return(idm)
  }else{
    print("Error:The indicator matrix setup is wrong!")
  }
}

#Remove within-seq gaps only
rmgap1 = function(A,g){
  # Av   = str_split(A,'')
  # g    = IRangesList(lapply(Av, function(x){IRanges(x=='-')}))
  rm.A1 = stri_sub_replace_all(A[[1]],from=sort(start(g[[1]])),to=sort(end(g[[1]])),replacement='')
  rm.A2 = stri_sub_replace_all(A[[2]],from=sort(start(g[[2]])),to=sort(end(g[[2]])),replacement='')
  rm.A  = DNAStringSet(c(rm.A1,rm.A2))
  if(all(width(rm.A)%%3 == 0)){
    return(rm.A)
  }else{
    print("Warning:removed-gap of each sequence is not multiple of three")
  }
}

#Trimming the non-gap edges 
edge_trim = function(M.012, posEdt, omechunk){
  
  ##put the non-gap-indicator-removed-site within one seq as NA
  for (k in seq(2)) {
    M.012[c(posEdt[[k]],posEdt[[k+2]]),,k]=NA
  }
  ##treat long-non-gap-site within the opposite seq as syn(0), keep the same seq as NA. 
  # M.012[posEdt[[3]],,2]=rep(0,omechunk)
  # M.012[posEdt[[4]],,1]=rep(0,omechunk)
  
  return(M.012)
}


#Remove all gap-positioned string
rmgap2 = function(A,As,g){
  g    = unlist(g)
  rm.A = stri_sub_replace_all(A,from=sort(start(g)),to=sort(end(g)),replacement='')
  rm.A = DNAStringSet(rm.A)
  if(all(width(rm.A)%%3 == 0)){
    return(rm.A)
  }else{
    print("Warning:removed-gap sequences are not multiple of three")
    break
  }
}

#Generate obs-codon matrix
countN = function(A, codonstrs){
  seqs = str_split(A,'') 
  len  = width(A)[1]
  #codon trans matrix  
  nmat = matrix(0,64,64)
  i=1
  while(i<len) {
    c1 = paste0(seqs[[1]][i:(i+2)], collapse = '')
    c2 = paste0(seqs[[2]][i:(i+2)], collapse = '')
    coor1 = which(codonstrs %in% c1)
    coor2 = which(codonstrs %in% c2)
    nmat[coor1,coor2] = nmat[coor1,coor2] + 1
    i=i+3
  }
  
  return(nmat)
}


#test A
# A = DNAStringSet(c("A---AAAAA","AAAAA---A")) #2
# A = DNAStringSet(c("AA---AAAA","AAAAA---A")) #5
# A = DNAStringSet(c("---AAAAAA","AAA---AAA")) #6
# A = DNAStringSet(c("AAAAAA","AAA---")) #10
# A = DNAStringSet(c("ATCGATCGA","A---A---A")) #3
###################################################################
ziqi_prob_wz = function(A,g0,e0,Pmat,codonstrs,syn,f0,wz0,omegal,norm.k1,norm.k2){
  
  
  As    = str_split(A,'')
  lenA  = length(As[[1]])
  g     = IRangesList(lapply(As, function(x){IRanges(x=='-')}))
  
  
  #inserted freq prob
  insP = prob_insert(As[[2]],g[[1]],f0)
  
  #end prob
  endP = end_prob(lenA,g,g0,e0)
  
  ###############################################
  #summary the gap indicator matrix
  #default omegal:12 
  sumG  = gap_sum(As,g,syn,omegal)
  N.012 = sumG[[1]]
  posEdt= sumG[[2]]
  num.g = sumG[[3]] 
  len.g = sumG[[4]] 
  
  #average gap length; ext prob
  avg.g = len.g/num.g
  if(all(len.g==0)){
    avg.g = c(0,0)
  }else if(len.g[1]==0){
    avg.g[1] = 0
  }else if(len.g[2]==0){
    avg.g[2] = 0
  }
  scoreE = (len.g[1]-num.g[1])*log(e0[1])+num.g[1]*log(1-e0[1]) + (len.g[2]-num.g[2])*log(e0[2])+num.g[2]*log(1-e0[2])
  
  
  #Generate the non-gap indicator matrix & post-editing
  rA     = rmgap1(A,g)             
  rAs    = str_split(rA,'')
  M.012  = nogap_sum(As, rAs, syn, omegal)
  M.trim = edge_trim(M.012, posEdt, omegal/3)          
  
  
  ##############################################
  #gap prob
  #scoreG = sum(N.012[,1]*log(g0*wz0) + N.012[,2]*log(g0))
  scoreG = sum(N.012*log(g0)) + sum(N.012[,1]*log(wz0))
  
  #non-gap prob
  k1 = t(t(M.trim[,,1])*norm.k1)
  k2 = t(t(M.trim[,,2])*norm.k2)
  scoreM = sum(log(1-g0 + (1-wz0)*g0*rowSums(k1)), na.rm=TRUE) +
           sum(log(1-g0 + (1-wz0)*g0*rowSums(k2)), na.rm=TRUE) 
  
  ##################################################################
  #substituion prob
  rA2    = rmgap2(A,As,g)
  dat    = countN(rA2,codonstrs)
  scoreP = sum(Pmat*dat)
  
  #sum LL
  score_ziqi = scoreE + scoreG + scoreM + scoreP + insP + endP
  
  #summary stat
  res = list(score_ziqi,N.012,avg.g,M.trim,dat) 
  return(res)
}




