### October 6, 2020
# functions from Ian's GitHub
# added eta parameter
# added gamma and eta function options

#library(mvtnorm)
#library(Rcpp)
#library(RcppArmadillo)
#library(mvtnorm)
#sourceCpp("~/Box/sleep_code/LSE_fast.cpp")
#sourceCpp("~/sleep_code/LSE_fast.cpp")

Mod2Orig = function(tmat,d0,format="%m/%d/%Y %H:%M:%S",tz="EST"){
  tmat_out = matrix(NA,nrow=nrow(tmat),ncol=2)
  t0=as.POSIXct(paste(d0,"00:00:00",sep=" "),tz=tz,format)
  for(i in 1:nrow(tmat)){
    tmat_out[i,1]=strftime(as.POSIXct(as.numeric(t0)+60*60*(tmat[i,1]+24*(tmat[i,3]-1)),tz=tz,origin="1970-01-01"),tz="EST",format)
    tmat_out[i,2]=strftime(as.POSIXct(as.numeric(t0)+60*60*(tmat[i,2]+24*(tmat[i,3]-1)),tz=tz,origin="1970-01-01"),tz="EST",format)
  }
  return(data.frame(t0=tmat_out[,1],t1=tmat_out[,2],stringsAsFactors=F))
}

Orig2Mod = function(tmat,anchor_hr,format="%m/%d/%Y %H:%M:%S",tz="EST"){
  tmat_out = matrix(NA,nrow=nrow(tmat),ncol=3)
  for(i in 1:nrow(tmat)){
    anchor_cur=anchor_hr*60*60+as.numeric(as.POSIXct(strftime(as.POSIXct(tmat[i,1],tz=tz,format=format,origin="1970-01-01"),tz=tz,format="%m/%d/%Y"),tz=tz,format="%m/%d/%Y"))
    if(as.numeric(as.POSIXct(tmat[i,1],tz=tz,format=format,origin="1970-01-01"))<anchor_cur){
      anchor_cur=anchor_hr*60*60+as.numeric(as.POSIXct(strftime(as.POSIXct(as.numeric(as.POSIXct(tmat[i,1],tz=tz,format=format,origin="1970-01-01"))-24*60*60,tz=tz,origin="1970-01-01"),tz=tz,format="%m/%d/%Y"),tz=tz,format="%m/%d/%Y"))
    }
    if(i==1){anchor0=anchor_cur}
    tmat_out[i,1]=anchor_hr+(as.numeric(as.POSIXct(tmat[i,1],tz=tz,format=format,origin="1970-01-01"))-anchor_cur)/(60*60)
    tmat_out[i,2]=anchor_hr+(as.numeric(as.POSIXct(tmat[i,2],tz=tz,format=format,origin="1970-01-01"))-anchor_cur)/(60*60)
    tmat_out[i,3]=round((anchor_cur-anchor0)/(24*60*60))+1
  }
  return(tmat_out)
}


## estimating model parameters
trap_quad_points = function(numpts){
  maxpt=qnorm(1-.5/numpts)
  pts1=seq(from=-maxpt,to=maxpt,length.out=numpts)
  wts1 = dnorm(pts1)
  dx=2*maxpt/(numpts-1)
  ptinds=as.matrix(expand.grid(rep(list(1:numpts),2)))
  pts=matrix(pts1[ptinds],nrow(ptinds),2)
  wts=apply(ptinds,1,function(xx) prod(wts1[xx])*dx*dx)
  return(list('pts'=pts,'wts'=wts))
}

d_w_cond_x = function(t_init,wt,xs,xw,lambda_s,lambda_w,mu_s,mu_w){
  if(xs>=xw || lambda_s <=0 || lambda_w <= 0 || mu_s >= mu_w){return(0)}
  if(t_init<xs){
    denom = 1-exp(-lambda_w*(xs-t_init))+exp(-lambda_s*(xs-t_init))-exp(-lambda_s*(xw-t_init))+exp(-lambda_w*(xw-t_init))
  }else if(t_init<xw){
    denom = 1-exp(-lambda_s*(xw-t_init))+exp(-lambda_w*(xw-t_init)) 
  }else{
    denom = 1
  }
  if(t_init+wt>xs && t_init+wt<xw){
    numer= lambda_s*exp(-lambda_s*wt)
  }else{
    numer= lambda_w*exp(-lambda_w*wt)
  }
  return(numer/denom)
}

#par order : lambda_s,lambda_w,mu_s,mu_w,sd_s,sd_w
QuadratureIntegral = function(mat,cpar,d){
  if(length(as.numeric(mat))==0){return(1)}
  if(is.numeric(mat)){
    mat=matrix(mat,nrow=1)
  }
  n_i=nrow(mat)
  quad_out=trap_quad_points(d)
  lik_wcondx = matrix(NA,nrow=d^2,ncol=n_i)
  for(j in 1:n_i){
    for(i in 1:d^2){
      lik_wcondx[i,j]=d_w_cond_x(t_init=mat[j,1],wt=mat[j,2]-mat[j,1],xs=quad_out$pts[i,1]*cpar[5]+cpar[3],xw=quad_out$pts[i,2]*cpar[6]+cpar[4],cpar[1],cpar[2],cpar[3],cpar[4])
    }
  }
  quad_scaled_w=quad_out$wts/(sum(quad_out$wts))
  logintegrand_p1=sum(rowSums(log(lik_wcondx))*quad_scaled_w)
  return(logintegrand_p1)
}

#par order : lambda_s,lambda_w,mu_s,mu_w,sd_s,sd_w
LogLikelihood =function(dat,cpar,d,gamma,eta){
  if(cpar[5]<=0 || cpar[6]<=0 || cpar[4]<=cpar[3] || cpar[1] < 0 || cpar[2] <0 || cpar[2]<cpar[1]){
    return(-Inf)
  }
  labels=unique(dat[,3])
  ls_ids = list()
  for(i in 1:length(labels)){
    ls_ids[[i]]=which(dat[,3]==i)
  }
  liktot=0
  for(i in 1:length(labels)){
    mat=dat[ls_ids[[i]],1:2]
    liktot=liktot+QuadratureIntegral(mat,cpar,d=5)
  }
  #Bayesian prior on sleep duration
  POP_AVG_DUR=7.4
  POP_SD_DUR=1.2
  liktot=liktot+(eta*length(labels)^gamma)*log(dnorm(abs(cpar[4]-cpar[3]),mean=POP_AVG_DUR,sd=POP_SD_DUR))
  #Bayesian prior on sd of time to sleep
  liktot=liktot+(eta*length(labels)^gamma)*log(dnorm(cpar[5],mean=1.2,sd=.5))
  #Bayesian prior on sd of time to wake
  liktot=liktot+(eta*length(labels)^gamma)*log(dnorm(cpar[6],mean=1.2,sd=.5))
  return(liktot)
}

InitialParameters = function(mat_mod,anchor_t=14){
  MIN_LAMBDA = .00001 # smallest value lambda can be
  itrvl_len_v = seq(6,9,.5)
  out_ls=list()
  ratio_v = rep(NA,length(itrvl_len_v))
  # find mu_s0, mu_w0
  for(j in 1:length(itrvl_len_v)){
    itrvl_len=itrvl_len_v[j]
    start_vals=seq(0,24-itrvl_len,.25)
    frac_vals = rep(NA,length(start_vals))
    for(i in 1:length(start_vals)){
      frac_vals[i]=length(intersect(which(mat_mod[,1]>anchor_t+start_vals[i]),which(mat_mod[,1]<anchor_t+start_vals[i]+itrvl_len)))/nrow(mat_mod)
    }
    mu_s0=start_vals[order(frac_vals)[1]]+anchor_t
    mu_w0=mu_s0+itrvl_len
    # find rate during average sleep and average waking interval
    ndays= length(unique(mat_mod[,3]))
    lambda_s0=max(c((min(frac_vals)*nrow(mat_mod)/ndays)/itrvl_len,MIN_LAMBDA))
    lambda_w0=max(c(((1-min(frac_vals))*nrow(mat_mod)/ndays)/(24-itrvl_len),MIN_LAMBDA))
    out_ls[[j]]=list(mu_s0,mu_w0,lambda_s0,lambda_w0)
    ratio_v[j]=lambda_s0/lambda_w0
  }
  return(unlist(out_ls[[order(ratio_v)[1]]]))
}

#par order : lambda_s,lambda_w,mu_s,mu_w,sd_s,sd_w
GridSearchInitPars = function(mat_mod,mu_s0,mu_w0,lambda_s0,lambda_w0,d){
  sd_s_v = c(.25,.5,1)
  sd_w_v = c(.25,.5,1)
  minval=-Inf
  for(sd_s in sd_s_v){
    for(sd_w in sd_w_v){
      par_v=c(lambda_s0,lambda_w0,mu_s0,mu_w0,sd_s,sd_w)
      curval=LogLikelihood(mat_mod,par_v,d,gamma,eta)
      if(curval>minval){
        cur_par=par_v
        minval=curval
      }
    }
  }
  return(par_v)
}

#par order : lambda_s,lambda_w,mu_s,mu_w,sd_s,sd_w
GetIndSleepEstimates =function(mat_mod,cpar){
  labels=1:max(mat_mod[,3])      ### this was labels=unique(mat_mod[,3])
  ls_ids = list()
  for(i in 1:length(labels)){
    ls_ids[[i]]=which(mat_mod[,3]==i)
  }
  xmat = matrix(NA,nrow=length(labels),ncol=3)
  for(i in 1:length(labels)){
    mat=mat_mod[ls_ids[[i]],1:2]    
    if(length(mat)==2) {               ###### added this
      mat=matrix(mat,nrow=1)        
    } else if (length(mat)==0) {       ###### added this
      xmat[i,]=c(cpar[3],cpar[4],labels[i])
      next
    }
    g3=function(par_v){
      tot=0
      for(j in 1:nrow(mat)){
        tot=tot-log(d_w_cond_x(t_init=mat[j,1],wt=mat[j,2]-mat[j,1],xs=par_v[1],xw=par_v[2],cpar[1],cpar[2],cpar[3],cpar[4]))
      }
      tot=tot-dnorm(par_v[1],cpar[3],cpar[5],log=TRUE)-dnorm(par_v[2],cpar[4],cpar[6],log=TRUE)
      ##### PENALIZE SLEEP DURATION
      IND_SD_DUR=1
      ############ SET GAMMA HERE #################
      gamma=0 #gamma can range from 0 to 1. 1/0 gives largest/smallest penalties
      #eta=0
      #############################################
      tot=tot-(nrow(mat)^gamma)*dnorm(par_v[2]-par_v[1],mean=cpar[4]-cpar[3],IND_SD_DUR,log=TRUE)
      ##### 
      return(tot)
    }
    optim.out3=optim(par=cpar[3:4],g3,control=list(maxit=1000))
    xmat[i,]=c(optim.out3$par,labels[i])
  }
  return(xmat)
}

#par order : lambda_s,lambda_w,mu_s,mu_w,sd_s,sd_w
NumMaxLogLik = function(outmat_mod,gamma,eta,anchor_t=14,d=3){
  iparinit=InitialParameters(outmat_mod,anchor_t)
  ipar=GridSearchInitPars(outmat_mod,mu_s0=iparinit[1],mu_w0=iparinit[2],lambda_s0=iparinit[3],lambda_w0=iparinit[4],d)
  #  ipar = c(1,3,mu_s+24+2,mu_w+24+2,2,.5)
  g1=function(cpar){
    return(-LogLikelihood(outmat_mod,cpar,d,gamma,eta))
  }
  optim.out1=optim(par=ipar,g1,control=list(maxit=1000))
  return(optim.out1$par)
}










