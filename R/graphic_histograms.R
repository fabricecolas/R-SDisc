`graphic_histograms` <-
function(cdata,fileName="2007-11-24_IMG_distributions.ps",xlim=NULL){
        local_data <- cdata[["data"]][,cdata[["canalysis_variables"]]]
        local_data <- as.data.frame(local_data)
	par('mai'=c(0.5,0.5,0.5,1),'mfrow'=c(4,5))
        #
        hist_list  <- list()
        hist_stats <- list()
        hist_counts <- c()
        models <- cdata[["model"]]
        for(i in colnames(local_data)){
           hist_list[[i]]  <- hist(local_data[,i],plot=FALSE)
           local_stat_txt <- ""
           hist_stats[[i]] <- list(txt="",xlim=list(),ylim=list())
           for(j in names(models)){
               local_model <- models[[j]][["model"]][[i]]
               if(models[[j]][["type"]] == "lm")
                  local_stat_txt <- paste2(local_stat_txt,local_model[["formula"]],"\n")
               else
                  local_stat_txt <- paste2(local_stat_txt,sprintf("%1.2f",local_model[[1]]),", ")
           }
           local_stat_txt <- paste2(local_stat_txt,length(table(local_data[,i]))," distinct values")
           hist_stats[[i]][["txt"]] <- local_stat_txt
           hist_stats[[i]][["xlim"]] <- c(min(local_data[,i],na.rm=TRUE),max(local_data[,i],na.rm=TRUE))
#           hist_stats[[i]][["xlim"]] <- c(quantile(data.matrix(local_data[,i]),probs=c(0.05),na.rm=TRUE),
#                quantile(data.matrix(local_data[,i]),probs=c(0.95),na.rm=TRUE))
           hist_stats[[i]][["ylim"]] <- range(hist_list[[i]][["counts"]])
        }
        #
	for(i in 1:ncol(local_data)){
            par('lab'=c(3,3,3),  # number of ticks
                'mar'=c(2.5,2.5,1.5,1) # margins bottom,left,top,right
                )
            plot(hist_list[[i]],
                main=dimnames(local_data)[[2]][i],
                xlab="",ylab="",xlim=hist_stats[[i]][["xlim"]],ylim=hist_stats[[i]][["ylim"]]
                #,ylog=TRUE
                ) #!!!!
            text(x=mean(hist_stats[[i]][["xlim"]]),y=hist_stats[[i]][["ylim"]][2]*0.85,hist_stats[[i]][["txt"]])
            }
		
}

