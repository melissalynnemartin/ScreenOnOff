##### 10/6/20
source("~/sleep_code/beiwe_sleep.R")

setwd("/storage/xia_mobile/beiwe_data/raw_data")
### ^ this is where everything is stored, but for now you can try the code on the data here:
# setwd("~/Box/TedSleep/data")

ids = list.files()

do.beiwe = function(id) {
  tryCatch({
    beiwe = beiwe_sleep(id,gamma=0,eta=0.1)
    save(beiwe,file=paste0("~/sleep_code/results_2020-10-18/beiwe/",id,"_beiwe.RData"))
  },error=function(e){})
}

mclapply(ids,do.beiwe,mc.cores=24)