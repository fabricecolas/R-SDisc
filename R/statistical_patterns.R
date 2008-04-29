`statistical_patterns` <-
function(data,fun_list=list(avg=mean,median=median)){
        data_out <- list()
        for(fun in names(fun_list))
                data_out[[fun]] <- compute_pattern(data,fun=fun_list[[fun]])
        return(data_out)
}

