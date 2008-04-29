`missing_at_random` <-
function(data,class){
   # AREA UNDER THE UN-CERTAINTY CURVE FOR EACH TWO MODELS
   auuc         <- mean(1-apply(class,1,max,na.rm=TRUE))
   auuc_sd      <- sd(  1-apply(class,1,max,na.rm=TRUE))
   out          <- matrix(NA,ncol(class),2,dimnames=list(1:ncol(class),c("AUUC","AUUC sdev")))
   out[1,1]     <- auuc
   out[1,2]     <- auuc_sd
   return(list(auuc = auuc, sd = auuc_sd, out = out))
}

