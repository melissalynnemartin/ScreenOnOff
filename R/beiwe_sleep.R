# 10/6/2020
# FIXED - times were off by an hour

setwd("~")
source("~/sleep_code/sleep_functions_2020-10-06.R")
library(data.table)
#library(zoo)
#library(dplyr)
library(parallel)

#path = "TedSleep/data/3e1xnphl/power_state"

#sub.paths = list.files(path="TedSleep/data",full.names=TRUE)
#paths = paste0(sub.paths,"/power_state")

#setwd("~/TedSleep/data")

beiwe_sleep = function(id,gamma,eta) {
  
  # list files in given path
  files = list.files(path=paste0(id,"/power_state"), full.names=TRUE)
  
  # create a list where each element is a power state csv file from path
  data = list()
  for (i in 1:length(files)) {
    data[[i]] = read.csv(files[i], header=TRUE)
  }
  # collapse the list to make a single data frame
  df = rbindlist(data)
  
  ##### why add time variable?
  #  ref.time = as.POSIXct("1970-01-01")
  #  df$time = ref.time + df$timestamp/1000 ##################### CHECK THIS
  
  # only grab locked/unlocked or screen on/off events
  if (ncol(df)==4) {
    df.use = df[which(df$event=="Locked" | df$event=="Unlocked"),]
  } else if (ncol(df)==3) {
    df.use = df[which(df$event=="Screen turned off" | df$event=="Screen turned on")]
  }
  
  # define outmat_orig data frame with screen off start time t0 and screen on start time t1
  outmat_orig = data.frame(t0=vector(length=ceiling(nrow(df.use)/2)+1), t1=vector(length=ceiling(nrow(df.use)/2)+1))
  
  counter = 0
  
  # grab screen off/locked intervals
  for (i in 1:(nrow(df.use)-1)) {
    if ((df.use[i,"event"]=="Locked" & df.use[i+1,"event"]=="Unlocked") | 
        (df.use[i,"event"]=="Screen turned off" & df.use[i+1,"event"]=="Screen turned on")) {
      counter=counter+1
      outmat_orig$t0[counter] = as.POSIXct(df.use$timestamp[i]/1000,origin="1970-01-01")
      outmat_orig$t1[counter] = as.POSIXct(df.use$timestamp[i+1]/1000,origin="1970-01-01")
    } else next
  }
  
  
  outmat_orig = outmat_orig[which(outmat_orig$t0 != 0),]
  
  outmat_orig$t0 = as.POSIXct(outmat_orig$t0, origin="1970-01-01")
  #outmat_orig$t0 = as.POSIXct(format(outmat_orig$t0), tz="UTC")
  #attr(outmat_orig$t0,"tzone") = "America/New_York"
  
  outmat_orig$t1 = as.POSIXct(outmat_orig$t1, origin="1970-01-01")
  #outmat_orig$t1 = as.POSIXct(format(outmat_orig$t1), tz="UTC")
  #attr(outmat_orig$t1,"tzone") = "America/New_York"
  
  outmat_mod = Orig2Mod(outmat_orig, anchor_hr=14)
  
  anchor_t = 14
  mle.out=NumMaxLogLik(outmat_mod,anchor_t,gamma=gamma,eta=eta)
  
  xest=GetIndSleepEstimates(outmat_mod,mle.out)
  
  d0=format(outmat_orig[1,1], "%m/%d/%Y") #### what exactly should this be? for id 1, why is it two days later?
  xest_orig=Mod2Orig(xest,d0,format="%m/%d/%Y %H:%M:%S")
  names(xest_orig)=c("bedtime","wake-up time")
  
  bedtime_est = as.POSIXct(xest_orig[,1],format="%m/%d/%Y %H:%M:%S")
  waketime_est = as.POSIXct(xest_orig[,2],format="%m/%d/%Y %H:%M:%S")
  bed_date = as.Date(substr(xest_orig[,1],1,10),format="%m/%d/%Y")
  wake_date = as.Date(substr(xest_orig[,2],1,10),format="%m/%d/%Y")
  date_night = ifelse(bed_date==wake_date,bed_date-1,bed_date)
  date_night = as.Date(date_night,origin="1970-01-01")
  #df = data.frame(method="proposed",date_night=date_night,bed_time=bedtime_est,wake_time=waketime_est)
  
  sleep_t_h=floor(mle.out[3]%%24)
  sleep_t_m=floor((mle.out[3]%%24-floor(mle.out[3]%%24))*60)
  if(sleep_t_m<10){
    sleep_t=paste(sleep_t_h,":0",sleep_t_m,sep="")
  }else{
    sleep_t=paste(sleep_t_h,":",sleep_t_m,sep="")
  }
  wake_t_h=floor(mle.out[4]%%24)
  wake_t_m=floor((mle.out[4]%%24-floor(mle.out[4]%%24))*60)
  if(wake_t_m<10){
    wake_t=paste(wake_t_h,":0",wake_t_m,sep="")
  }else{
    wake_t=paste(wake_t_h,":",wake_t_m,sep="")
  }
  cat(paste(" Avg. time to sleep = ",sleep_t," (+/- ",round(mle.out[5],1)," hour)\n",sep="")
      ,(paste("Avg. time to wake  = ",wake_t," (+/- ",round(mle.out[6],1)," hour)\n",sep=""))
      #,(paste("Correlation between time to sleep and time to wake = ",round(mle.out[5],2),"\n",sep=""))
      ,(paste("Rate (per hour) of frequency of phone use while asleep = ", round(mle.out[1],5),"\n",sep=""))
      ,(paste("Rate (per hour) of frequency of phone use while awake = ", round(mle.out[2],5),"\n",sep="")))
  
  ## maybe change format of mu_s and mu_w
  param_est = data.frame(id=id, mu_s=mle.out[3],mu_w=mle.out[4],sd_s=mle.out[5],sd_w=mle.out[6],
                         lambda_s=mle.out[1],lambda_w=mle.out[2])
  sleep_est = data.frame(id=id,date_night=date_night,beiwe_bedtime=bedtime_est,beiwe_waketime=waketime_est)
  
  return(list(param_est=param_est,sleep_est=sleep_est))
}

#mclapply(paths,beiwe_sleep,mc.cores=13)
