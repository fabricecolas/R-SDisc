`analysis` <-
function(x){
   cat('\n Plot cdata for EDA\n')
   plot(attr(x,"cdata"))
   cat('\n Perform cluster analysis\n')
   x <- doModeling(x)
   cat("\nSave best models as CSV files\n")
   write.cresult(x)
   cat('\n Plot cresult\n')
   plot(x)
   return(x)
}

`funStatsAuuc` <-
function(class){
   # AREA UNDER THE UN-CERTAINTY CURVE FOR EACH TWO MODELS
   auuc         <- mean(1-apply(class,1,max,na.rm=TRUE))
   auuc_sd      <- sd(  1-apply(class,1,max,na.rm=TRUE))
   out          <- matrix(NA,ncol(class),2,dimnames=list(1:ncol(class),c("AUUC","AUUC sdev")))
   out[1,1]     <- auuc
   out[1,2]     <- auuc_sd
   return(list(auuc = auuc, sd = auuc_sd, out = out))
}

`stabilityAnalysis` <-
function(x,q,rseed=6013,nnoise=10,nrep=10,ps=TRUE, ncolors=9,minmax=c(0.5,1)){
   settings <- strsplit(q,',')[[1]]
   cdata <- print(attr(x, "cdata"))
   r <- list()
   v <- matrix(NA,nnoise,nrep*(nrep-1)/2,dimnames=list(sprintf("%4f",(1/2)^(1:nnoise)),list()))
   # NOISE LEVELS
   for(sd_l in (1/2)^(1:nnoise)){ 
      l <- sprintf("%4f",sd_l)
      r[[l]] <- list()
      for(ridx in 1:nrep){
         s <- rseed+ridx
         r[[l]][[s]] <- list()
         set.seed(s)
         Y <- cdata+matrix(rnorm(nrow(cdata)*ncol(cdata),sd=sd_l,mean=0),nrow(cdata),ncol(cdata))
         r[[l]][[s]] <- clusterModelBasedEM(Y,modelName=settings[1],G=settings[2],rseed=as.numeric(settings[3]))
         r[[l]][[s]][["Y"]] <- Y
      }
      tmp <- c()
      for(ridx1 in 1:nrep){
         s1 <- rseed+ridx1
         for(ridx2 in 1:nrep){
            s2 <- rseed+ridx2
            if(ridx2>ridx1)
               tmp <- c(tmp, V.score(r[[l]][[s1]]$labelling, r[[l]][[s2]]$labelling))
         }
      }
      v[l,] <- sort(tmp)
   }
   breaks <- seq(minmax[1],minmax[2],(minmax[2]-minmax[1])/ncolors)
   if(ps){
      postscript(paste2(attr(attr(x,'cdata'),'figDir'),"/stabilityAnalysis_",gsub(',','_',q),".ps"))
      image(v,col=brewer.pal(ncolors,"Greys"),breaks=breaks,axes=FALSE,ylab='Association level V, quantiles',cex.lab=1.7)
      contour(v, add=TRUE,labcex=1.3)
      axis(1, at=seq(0,1,0.11),labels=sprintf("%.1f%%",100*as.numeric(row.names(v))),las=2,cex.axis=1.7)
      axis(2, at=seq(0,1,0.2),labels=sprintf("%d%%",100*seq(0,1,0.2)),cex.axis=1.7)
      title(paste(settings,sep="",collapse=" "))
      graphics.off()
   }
   v <- structure(v,class='stabilityAnalysis')
   return(v)
}

`compareCresult` <-
function(x,y,probs=seq(0,1,0.1)){
   m <- data.frame(attr(x,'cfun_params'),V=NA)
   row.names(m) <- apply(as.matrix(m[,1:3]),1,paste,sep='',collapse=',') 
   for(i in names(x)){
      tmp <- table(data.frame(x=attr(x[[i]],'model')$labelling, y=attr(y[[i]],'model')$labelling))
      tmp_summary <- summary(tmp)
      n <- as.numeric(tmp_summary["n.cases"])
      X2 <- as.numeric(tmp_summary["statistic"])
      k <- min(dim(tmp))
      m[i,'V'] <- sqrt(X2/(n*(k-1)))
   }
   m <- na.omit(m)
   m_x <- xtabs(formula=V~modelName+G,aggregate(data.frame(V=m[,'V']),by=list(modelName=m[,'modelName'],G=m[,'G']),mean,na.rm=TRUE))
   return(list(summary=m_x, m=m))
}

`doPattern` <-
function(data,fun=mean){
   data                 <- data.matrix(aggregate(data,fun,by=list(class=data[,"class"])))[,-1]
   row.names(data)      <- 1:nrow(data)	
   if(!is.na(match("class", dimnames(data)[[2]])))
      data              <- data[,-match("class", dimnames(data)[[2]])]
   return(data)
}


V.score <- function(x, y=NULL){
   if(length(x) == length(y))
      x <- as.table(ftable(x,y))
   if(is.table(x)){
      xSummary <- summary(x)
      n <- as.numeric(xSummary["n.cases"]) 
      X2 <- as.numeric(xSummary["statistic"]) 
      k <- min(dim(x)) 
      V <- sqrt(X2/(n*(k-1)))
      return(V)
   }
   else
      cat('x must be a table\n')
}

`compareModel` <-
function(x, y=NULL, fmt=c('%d','%.2f')){
   if((class(x) == 'cresult') & (is.null(y))){
      m <- list()
      mPairs <- t(combn(getBestModel(x),m=2))
      cNames <- gsub(",","_",paste("(",apply(mPairs,1,paste,collapse=")x("),")",sep=""))
      for(i in 1:nrow(mPairs)){
         m1 <- x[[which(mPairs[i,1] == names(x))]]
         m2 <- x[[which(mPairs[i,2] == names(x))]]
         if( !is.na(attr(m1,"model")$loglik) && !is.na(attr(m2,"model")$loglik))
            m[[cNames[i]]] <- compareModel(m1,m2)
      }
      return(m)
   }
   else if(class(x) == 'cmodel' & class(y) == 'cmodel'){
      xLabel <- attr(x,"model")$labelling
      yLabel <- attr(y,"model")$labelling
      xG <- attr(x,"model")$G
      yG <- attr(y,"model")$G
      xyTab <- as.table(ftable(xLabel, yLabel))
      # CHI2_STATS
      chi2 <- chisq.test(xyTab, simulate.p.value=TRUE)
      # CRAMER'S V 
      V <- V.score(xLabel, yLabel)  
      # ROW AND COLUMN ORDERING
      rowOrder <- hclust(dist(xyTab[1:xG, 1:yG]) )$order
      colOrder <- hclust(dist(t(xyTab[1:xG, 1:yG])))$order
      xyTab <- xyTab[rowOrder, colOrder]
      xyTab[xyTab == 0] <- NA
      # BIND STATS
      xyTab <- rbind(xyTab, CHI2=NA)
      xyTab[nrow(xyTab),1:2] <- c(chi2$p.value, chi2$statistic)
      xyTab <- rbind(xyTab, V=NA)
      xyTab[nrow(xyTab),1] <- V
      return(xyTab)
   }
}

`clusterModelBasedEM` <-
function(data, G=3, modelName="EII", rseed=6013) {
   G <- as.numeric(G)
   modelName <- as.character(modelName)
   # PREPARE A RANDOM-UNIFORM MATRIX OF CLUSTER MEMBERSHIPS WITH G COLS AND N
   # ROWS (MIN IS 0 AND MAX 1)
   if(G >= 2){
      set.seed(rseed)
      cmembership_matrix <- matrix(runif(G*nrow(data)), nrow(data), G)
      cmembership_matrix <- cmembership_matrix/apply(cmembership_matrix, 1, sum)
   }
   else
      cmembership_matrix <- matrix(1,nrow(data),G)
   # ESTIMATE AN INITIAL MODEL WHICH WE PROVIDE AS STARTING POINT TO EM
   mstep_est <- mstep(data=data, modelName=modelName, z=cmembership_matrix,
                  warn=FALSE)
   # AND THEN, STARTS EM WITH THE E-STEP
   model <- em(data=data, modelName=mstep_est$modelName,
               parameters=mstep_est$parameters, warn=FALSE)
   model[["labelling"]] <- map(model[["z"]],warn=FALSE)
   # CALCULATE THE BIC SCORE
   model[["BIC"]] <- bic(modelName=modelName, loglik=model$loglik, n=nrow(data)
                           , d=ncol(data), G=G)
   return(model)
}

`funStatsGeneralization` <-
function(data,class
   , fun_classifier=list(model=model_naive_bayes, predict=predict_naive_bayes)
   , K=7, fun_eval=splitTraintestStrata, title = NULL) {
   class                <- map(class)
   ldata                <- cbind(print(data),class=class)
   # INITIALIZATION
   indexes <- testfold  <- list()
   models <- predictions<- contingency.tables <- perfmeasures <- list()
   G                    <- length(table(class))
   if(!(G == 0)){
      # GENERATE THE DIFFERENT STRATA FROM THE K-FOLD CROSS VALIDATION
      out_cv <- fun_eval(ldata,K)
      indexes <- out_cv[["indexes"]]
      K <- out_cv[["K"]]
      # PROCEED TO EVALUATION FOR EACH K-FOLD
      for(k in 1:K){
         # INIT THE K-TH TEST FOLD
         testfold[[k]]          <- list()
         for(g in 1:G)
            testfold[[k]]       <- append(testfold[[k]],row.names(ldata[ldata[,"class"] == g,][indexes[[g]][k,],]))
         testfold[[k]]          <- unlist(testfold[[k]])
         # DEFINE TRAIN AND TEST SETS
         cdata_train            <- setCdata(data = ldata[-pmatch(testfold[[k]], row.names(ldata)),]
                                                , settings = attr(data,"settings"))
         cdata_test             <- setCdata(data = ldata[ pmatch(testfold[[k]], row.names(ldata)),]
                                                , tdata = attr(cdata_train,"tdata")
                                                , settings = attr(data,"settings"))
         trainset               <- cbind(print(cdata_train),class=cdata_train[,"class"])
         testset                <- cbind(print(cdata_test),class=cdata_test[,"class"])
         # TRAIN DIFFERENT MODELS, NAMELY NAIVE BAYES, SVM AND 1-NN RK: SVM
         # USES G(G-1) CLASSIFIERS, NAIVE BAYES AND KNN ARE SINGLE MULTI-CLASS
         models[[k]]            <- fun_classifier[["model"]](as.formula("class ~ ."),trainset,testset)
         predictions[[k]]       <- fun_classifier[["predict"]](models[[k]],testset)
         contingency.tables[[k]]<- table(predictions[[k]],testset[,"class"])
         perfmeasures[[k]]      <- sum(diag(contingency.tables[[k]]))/sum(contingency.tables[[k]])
      }
      # SUMMARY STATISTICS
      perfs		        <- unlist(perfmeasures)
      out                       <-t( matrix( c( mean(perfs),
                                                1.96*sd(perfs,na.rm=TRUE),K)
                                             , 3
                                             , G
                                             , dimnames=list(list(title,"CI95%","K-folds"),as.list(1:G)),byrow=FALSE))
      out[2:G,]                 <- NA
      return(list(      models = models
                        , predictions = predictions
                        , contingency.tables = contingency.tables
                        , perfmeasures = perfmeasures
                        , avg=mean(perfs,na.rm=TRUE)
                        , sd=1.96*sd(perfs,na.rm=TRUE)
                        , out=out))
   }
   else return(NULL)
}

`fullproba_ftable` <-
function(z1, z2){
   m_out                <- matrix(NA,ncol(z1),ncol(z2))
   row_name_vector      <- unique(c(row.names(z1),row.names(z2)))
   m                    <- matrix(NA, length(row_name_vector), ncol(z1)+ncol(z2),dimnames=list(row_name_vector,list()))
   m[row.names(z1),1:ncol(z1)] <- z1
   m[row.names(z2),(ncol(z1)+1):(ncol(z1)+ncol(z2))] <- z2
   for(i in 1:ncol(z1))
      for(j in 1:ncol(z2))
         m_out[i,j] <- mean(m[,i] * m[,j+ncol(z1)])
   return(m_out)
}

`genCdataSettings` <-
function(data, grBreak=6, asCSV=FALSE){
   # DEFINE THE DEFAULT VALUES OF OUR CONF MATRIX
   cNames <- list('oddGroup','inCAnalysis','tFun','vParGroup','vParY','vHeatmapY')
   dValues <- c("factor_1",TRUE,"tAvg tSigma","var_group_1","NA","NA")
   x <- matrix(dValues,ncol(data),length(dValues),byrow=TRUE,
                  dimnames=list(colnames(data),cNames))
   x[,'vHeatmapY'] <- 1:nrow(x)
   xSeq <- 1+as.integer(nrow(x)/grBreak)
   xGroup <- matrix(c('oddGroup',NA,NA),nrow(x),3,byrow=TRUE)
   for(i in 0:(grBreak-1)){
      xGroup[(i*xSeq+1):min((i+1)*xSeq,nrow(x)),2] <- i+1
      xGroup[(i*xSeq+1):min((i+1)*xSeq,nrow(x)),3] <- 1:(min((i+1)*xSeq,nrow(x))-(i*xSeq))
   }
   x[,'oddGroup'] <- x[,'vParGroup'] <- apply(xGroup[,1:2],1,paste,collapse='_')
   x[,'vParY'] <- xGroup[,3]
   if(asCSV){
      fSettings <- paste2(today(),'_settings.csv')
      write.csv(x,file=fSettings)
      return(fSettings)
   }
   else
      return(x)
}

`getCdataPrcomp` <- function(x, cumproportion=0.95, scale=TRUE, center=FALSE){
   # RETRIEVE PREVIOUS PARAMETERS
   s1 <- as.matrix(attr(x,"settings"))
   # PRINCIPAL COMPONENTS, 95%
   pcomp <- prcomp(x[,getCvar(s1)], scale=scale, center=center)
   vars <- pcomp[["sdev"]]^2
   vars <- vars/sum(vars)
   first_comps <- (cumsum(vars) < cumproportion)
   pcomp <- predict(pcomp)[,first_comps]
   df2 <- cbind(x, pcomp)
   # SETTINGS
   s2 <- genCdataSettings(df2)
   s2[row.names(s1),] <- s1[,colnames(s2)]
   s2[,"inCAnalysis"] <- "FALSE"
   new_names <- row.names(na.omit(s2[s2[,"vParY"] == "NA",]))
   s2[new_names,"inCAnalysis"] <- "TRUE"
   s2[new_names,"oddGroup"] <- new_names
   s2[new_names,"vParGroup"] <- "PCA"
   s2[new_names,"vParY"] <- 10*(1:length(new_names)/length(new_names))
   s2[new_names,"tFun"] <- "tAvg tSigma"
   max_heatmap <- max(as.numeric(s2[,"vHeatmapY"]),na.rm=TRUE)
   s2[new_names,"vHeatmapY"] <- as.character((max_heatmap+1):(max_heatmap+length(new_names)))
   x2 <- setCdata(data=df2, settings=s2, prefix=paste2(attr(x,"prefix"),"_PRCOMP"), data_o=attr(x,'data_o'))
   return(x2)
}

`getCdataPrincomp` <- function(x, cumproportion=0.95, scale=TRUE, center=FALSE){
   # RETRIEVE PREVIOUS PARAMETERS
   s1 <- as.matrix(attr(x,"settings"))
   # PRINCIPAL COMPONENTS, 95%
   pcomp_m <- princomp(x[,getCvar(s1)], scale=scale, center=center)
   vars <- pcomp_m[["sdev"]]^2
   vars <- vars/sum(vars)
   first_comps <- names(vars)[cumsum(vars) < cumproportion]
   pcomp <- pcomp_m[["scores"]][,first_comps]
   df2 <- cbind(x, pcomp)
   # SETTINGS
   s2 <- genCdataSettings(df2)
   s2[row.names(s1),] <- s1[,colnames(s2)]
   s2[,"inCAnalysis"] <- "FALSE"
   new_names <- row.names(na.omit(s2[s2[,"vParY"] == "NA",]))
   s2[new_names,"inCAnalysis"] <- "TRUE"
   s2[new_names,"oddGroup"] <- new_names
   s2[new_names,"vParGroup"] <- "PCA"
   s2[new_names,"vParY"] <- 10*(1:length(new_names)/length(new_names))
   s2[new_names,"tFun"] <- "tAvg tSigma"
   max_heatmap <- max(as.numeric(s2[,"vHeatmapY"]),na.rm=TRUE)
   s2[new_names,"vHeatmapY"] <- as.character((max_heatmap+1):(max_heatmap+length(new_names)))
   x2 <- setCdata(data=df2, settings=s2, prefix=paste2(attr(x,"prefix"),"_PRINCOMP"), data_o=attr(x,'data_o'))
   attr(x2,'princomp') <- pcomp_m
   return(x2)
}

`getStatsFun` <-
function(fun_name="oddratios",...){
   # CHEMO-INF: CHI2TEST STATS AND RESIDUALS 
   if(fun_name == "chi2test")
      return(function(data,class) { return(funStatsChi2test(data,class,...)) })
   # CHEMO-INF: JOINT DISTRIBUTION
   if(fun_name == "jointdistrib")
      return(function(data,class) { return(funStatsJointDistrib(data,class,...)) })
   # ODD RATIOS CHARACTERIZING EACH CLUSTER
   if(fun_name == "oddratios")
      return(function(data,class) { return(funStatsOddratios(data,class,...)) })
   # AREA UNDER THE UNCERTAINTY CURVE
   if(fun_name == "auuc")
      return(function(data,class) { return(funStatsAuuc(class)) })
   # GENERALIZATION ESTIMATES OF MACHINE LEARNING ALGORITHMS
   if(fun_name == "gen_naive_bayes")
      return(function(data,class) { return(
         funStatsGeneralization( data, class, title = "naive Bayes",
            fun_classifier=list(model = model_naive_bayes,predict =
            predict_naive_bayes), ...)) })
   if(fun_name == "gen_knn")
      return(function(data,class) { return(
         funStatsGeneralization( data, class, title = "1 nearest neighbour",
            fun_classifier=list(model = model_knn,predict = predict_knn), ...))})
   if(fun_name == "gen_svm")
      return(function(data,class) { return(
         funStatsGeneralization( data, class, title = "Support Vector Machine",
            fun_classifier=list(model = model_svm,predict = predict_svm), ...))})
else
      return(function(data,class,fun=fun_name) { return(call(eval(fun), data, class,...)) })
}

`getBestModel` <-
function(x, n=NULL){
   bTable <- attr(attr(x,"bicanalysis"),'data')[, c("modelName","G","rseed","relativeBic")]
   if(is.null(n))
      idx <- 1:attr(x,"nTopModels")
   else if(n < 1)
      idx <- 1:length(x)
   else
      idx <- n
   bestBic <- bTable[na.omit(order(bTable[,"relativeBic"], decreasing=FALSE)[idx]),]
   bestModel  <- apply(bestBic[,c("modelName","G","rseed")],1,paste,collapse=",",sep="")
   return(as.character(bestModel))
}

`getCdataCC` <-
function(data,settings)
{
   var_names         <- getCvar(settings)
   data_mat          <- data[row.names(na.omit(data[,var_names])),]
   return(data_mat)
}


`getBicStats`          <- function(x, fun='mean', vAgg="BIC", vX="G",
   vY="modelName", fmt='%.2f'){
   rList <- list()
   # DEFINE FUNCTIONS 
   getRankArray <- function(m){return(apply(-m,c(2,3),rank))}
   fReshape <- function(m, pX=vX, pY=vY, pZ=vAgg){
      m <- reshape(m[,c(pX,pY,pZ)],idvar=pX, timevar=pY, direction='wide')
      dimnames(m) <- list(m[,1], gsub(paste2(pZ,'.'),'',colnames(m)))
      return(m[,-1])
   }
   fCI <- function(x,pFmt=fmt){
      x <- quantile(x,probs=c(0.025,0.975))
      x <- paste(sprintf(x, fmt=pFmt), collapse=', ')
      x <- gsub('$',')',gsub('^',' (',x))
      return(x)
   }
   fFun <- function(x, f=fun, pFmt=fmt){ return(sprintf(eval(call(f,x)),fmt=pFmt)) }
   # MAKE THE ARRAY OF BIC SCORES
   vRSeed <- unique(x[,'rseed'])
   tmp <- fReshape(x[x$rseed == vRSeed[1],])
   if(length(vRSeed)>1){
      for(s in 2:length(vRSeed))
         tmp <- abind(tmp, fReshape(x[x$rseed == vRSeed[1],]), along=3)
      dimnames(tmp)[[3]] <- vRSeed
      if(fun == 'rank')
         tmp <- abind(apply(getRankArray(-tmp), c(1,2), mean, na.rm=TRUE), 
                      apply(getRankArray(-tmp),c(1,2),fCI, pFmt=fmt),along=3) 
      else
         tmp <- abind(apply(tmp,c(1,2),fFun,f=fun,pFmt=fmt), apply(tmp, c(1,2), fCI, pFmt=fmt), along=3)
      tmp <- apply(tmp,c(1,2),paste2)
   }
   return(tmp)
}

`setBicanalysis` <- function(x){
   bicTable <- cbind(attr(x, "cfun_params"),BIC=NA,relativeBic=NA)
   for(i in 1:length(x))
      bicTable[i,"BIC"] <- attr(x[[i]],"model")[["BIC"]]
   bicTable[,"relativeBic"] <- 100*(bicTable[,"BIC"]/max(bicTable[,"BIC"],na.rm=TRUE)-1)
   attr(x,'bicanalysis') <- structure(table(bicTable[,1:2]),data=bicTable,class='bicanalysis')
   return(x)
}


`initPlot` <- 
function(cdata, figIdx, plotNbr=1, type="H"){
   # CREATE FIGURE DIRECTORY
   # NUMBER OF SUBPLOT IN A FIGURE
   plotI <- plotJ <- as.integer(sqrt(plotNbr))
   if(plotI^2 < plotNbr)
      plotJ <- plotJ+1 
   # FIGURE OUTPUT
   maiParam <- attr(cdata,'mai')
   fDir <- attr(cdata,'figDir')
   fFig <- paste2(fDir,'/oddGroup_',sprintf('%03d',figIdx),'-',type,'.pdf')
   if(type=='H')
      maiParam <- maiParam+c(0,0,0,0)
   else if(type=='BB')
      maiParam <- maiParam+c(0.5,0.2,0.5,0)
   else{
      maiParam <- maiParam+c(0.1,0.5,0.2,0)
      fFig <- paste2(fDir,'/',type,'.pdf')
      }
   print(fFig)
   pdf(fFig)
   par(mfrow = c(plotI,plotJ), mai=maiParam)
   figIdx <- figIdx+1
   return(figIdx)
}

`closePlot` <-
function(){
   dev.off()
}

`plot.cdata` <- 
function(x, ...) {
   xSettings <- as.matrix(na.omit(attr(x,"settings")[,c("vParGroup","vParY")]))
   splitGroup <- split(row.names(xSettings),xSettings[,"vParGroup"])
   figIdxBB <- figIdxH <- 1
   for(i in 1:length(splitGroup)){
      # PRODUCE BOXPLOT BY SPLIT GROUP
      figIdxBB <- initPlot(x, figIdxBB, plotNbr=1, type='BB')
      boxplot(as.data.frame(x[,splitGroup[[i]]]),'las'=2,cex=0.7, main=paste2(" Boxplot ",names(splitGroup)[i]))
      closePlot()
      figIdxH <- initPlot(x, figIdxH, plotNbr=length(splitGroup[[i]]), type='H')
      for(varName in splitGroup[[i]]){
         # PREPARE THE TEXT FOR EACH HISTOGRAM
         vIdx <- which(names(attr(x,"tdata")) == varName)
         xlabText <- paste2(varName,' ', length(table(x[,varName])), " values\n")
         if(length(vIdx) == 1){
            for(tn in names(attr(x,"tdata")[[vIdx]])){
               if(tn == "tAdjust")
                  xlabText <-  paste2(xlabText, tolower(gsub("transform_","",tn))
                     , "=", attr(x,"tdata")[[vIdx]][[tn]][["print"]], " ")
               else
                  xlabText <-  paste2(xlabText, tolower(gsub("transform_","",tn))
                     , "=", sprintf("%.2f",attr(x,"tdata")[[vIdx]][[tn]])," ")
               }
            }
         # PRODUCE HISTOGRAM BY SPLIT GROUP
         hist(x[,varName], main=NULL,xlab=xlabText, ylab="", cex=0.5, las=2,'new'=TRUE, border=NULL)
      }
      closePlot()
   }
}

`plot.cresult` <- 
function(x, query = NULL, ...) {
   q <- query
   if(is.null(q))   
      q <- 1:length(x)
   plotNbr <- length(attr(x,"fun_plot"))
   figIdx <- 1
   for(i in q){
      if(!is.na(attr(x[[i]],"model")$loglik)){
         figIdx <- initPlot(attr(x,'cdata'), figIdx, plotNbr=plotNbr, type=paste2('MM-',gsub(',','_',names(x)[i])))
         for(fun_plot in unlist(attr(x,"fun_plot")))
            try(fun_plot(x[[i]]))
         closePlot()
      }
   }
}

`doModeling` <-
function(x) {
   cdata        <- print(attr(x,"cdata"))
   for(i in 1:length(x)){
      cat("\n",names(x)[[i]])
      attr(x[[i]],"model") <- x[[i]](x[[i]],cdata)
      if(!is.na(attr(x[[i]],"model")$loglik)){
         cdata_l <- attr(x,"cdata")[,order(attr(attr(x,"cdata"),"settings")[,"vHeatmapY"],na.last=NA)] 
         cdata_l <- cbind(cdata_l,class=attr(x[[i]],"model")[["labelling"]])
         cat("-> Patterns")
         # COMPUTE PATTERNS
         attr(x[[i]],"pattern") <- NULL
         for(p in names(attr(x,"fun_pattern")))
            attr(x[[i]],"pattern")[[p]] <- doPattern(cdata_l,attr(x,"fun_pattern")[[p]])
         # FOR ORDERING, MAKE A DENDROGRAM FROM THE CENTER PATTERN (NUMBER 1)
         if(!is.null(attr(x[[i]],"pattern"))){
            cat("-> Dendros")
            avg_p <- attr(x[[i]],"pattern")[[1]]
            attr(x[[i]],"dendro_cluster") <- hclust(dist(avg_p))
            attr(x[[i]],"dendro_var") <- hclust(dist(t(avg_p)))
            # REORDER THE PATTERNS
            cat("-> Ordering")
            for(p in names(attr(x,"fun_pattern")))
               attr(x[[i]],"pattern")[[p]] <- attr(x[[i]],"pattern")[[p]][attr(x[[i]],"dendro_cluster")$order,]
            attr(x[[i]],"pattern")[["class_count"]] <- table(cdata_l[,"class"])[attr(x[[i]],"dendro_cluster")$order]
            # SELECT DIVERGING COLORS FROM RColorBrewer
            attr(x[[i]],"cluster_colors") <- brewer.pal(nrow(avg_p),"Set1")
         }
         # STATS
         fun_stats <- attr(x,"fun_stats")
         cat("-> Stats")
         for(n in names(fun_stats))
            attr(x[[i]],"stats")[[n]] <- fun_stats[[n]](data=attr(x,"cdata"), class=attr(x[[i]],"model")[["z"]])
      }
   }
   x <- setBicanalysis(x)
   fSave <- paste2(attr(attr(x,"cdata"),'baseDir'),"/IMAGE.RData")
   cat("\nSave modeling into ",fSave,'\n')
   save(list='x', file=fSave)
   return(x)
}

`model_knn` <-
function(formula, trainset, testset, k = 1, l = 0, prob = FALSE, use.all = TRUE){
   return(knn(    trainset[,-which(colnames(trainset) == "class")]
                , testset[ ,-which(colnames(testset) == "class")] 
                , trainset[,"class"]
                , k = k
                , l = l
                , prob = prob
                , use.all = use.all))
}

`model_naive_bayes` <-
function(formula, trainset, testset, laplace = 1){
   return(naiveBayes(formula, data=trainset, laplace = 1))
}

`model_svm` <-
function(formula, trainset, testset, kernel="linear",type="C-classification"){
   return(svm(formula, data=trainset,kernel=kernel,type=type))
}

`paste2` <-
function(...){
	return(paste(...,collapse="",sep=""))
}

`predict_knn` <-
function(model, testset){
   return(model)
}

`predict_naive_bayes` <-
function(model,testset){
   return(map(predict(model,testset[,-match("class",colnames(testset))],type="raw")))
}

`predict_svm` <-
function(model, testset){
   return(predict(model,testset[,-match("class",colnames(testset))]))
}

`write.cresult` <-
function(x, query = NULL){
   if(is.null(query))
      query             <- getBestModel(x)
   for(q in query){
      fCSV        <- paste2(attr(attr(x,'cdata'),'tabDir'),'/MM_',gsub(",","_",q),".csv") 
      print(fCSV)
      write.csv2(print(x[[q]]), file=fCSV)
   }
}

`patternMean` <- function(x){
   return(mean(x,na.rm=TRUE))
}

`patternMedian` <- function(x){
   return(median(x,na.rm=TRUE))
}

`patternLowquant` <- function(x){
   return(quantile(x,probs = 0.025,na.rm=TRUE))
}

`patternUpquant` <- function(x){
   return(quantile(x,probs = 0.975,na.rm=TRUE))
}

`patternMax` <- function(x){
   return(max(x,na.rm=TRUE))
}

`patternSd` <- function(x){
   return(sd(x,na.rm=TRUE))
}

`setCresult` <-
function(cdata=NULL, cfun=clusterModelBasedEM, cfun_settings=list(modelName=c("EII",
   "VII"), G=3:5, rseed=6013:6015), fun_pattern=list(mean=patternMean,
   median=patternMedian, lowquant=patternLowquant, upquant=patternUpquant) ,
   fun_plot=list(plotParcoord=getPlotFun(type="plotParcoord"),
   plot_legend=getPlotFun(type="plot_legend"),
   plot_image=getPlotFun(type="plot_image"),
   plot_dendro_cluster=getPlotFun(type = "plot_dendro_cluster"),
   plot_dendro_var=getPlotFun(type="plot_dendro_var")),
   fun_stats=list(oddratios=getStatsFun(fun_name="oddratios")),
   nTopModels=5){
   # COMPUTES THE COMBINATION OF THE CLUSTERING FUNCTION PARAMETERS 
   param_set_df         <- cfun_settings[[1]]
   i <- 2
   while(i <= length(cfun_settings)){
      param_set_df      <- merge(param_set_df,cfun_settings[[i]])
      colnames(param_set_df)<-names(cfun_settings)[1:i]
      i <- i+1
   }
   # CONSTRUCT A LIST OF 'CMODEL' DATA STRUCTURES
   cresult_out          <- list()
   for(i in 1:nrow(param_set_df)){
      id_name           <- paste(as.matrix(param_set_df[i, names(cfun_settings)]),collapse=",")
      fun_call          <- function(cmodel,cdata){return(do.call(cfun,append(list(data=cdata),attr(cmodel,"cfun_settings"))))}
      cresult_out[[id_name]] <- structure(fun_call,
         cfun_settings=param_set_df[i,], model=NULL, pattern=list(),
         dendro_cluster=NULL, dendro_var=NULL, cluster_colors=NULL,
         stats=list(), class="cmodel")
   }
   # PRE-EVAL FUN_PLOT TO RETRIEVE THE DIFFERENT PLOTTING FUNCTIONS
   fun_plot2 <- list()
   for(lfun in unlist(fun_plot))
      fun_plot2 <- c(fun_plot2, lfun(cdata))
   # CONSTRUCT A 'CRESULT' DATA STRUCTURE WITH CMODELS IN IT
   cresult_out <- structure(cresult_out, cdata=cdata,
      prefix=attr(cdata,"prefix"), cfun_params=param_set_df,
      fun_plot=fun_plot2, fun_stats=fun_stats, fun_pattern=fun_pattern,
      nTopModels=nTopModels, bicanalysis=NULL, rinfo=sessionInfo(),
      class = "cresult")
   return(cresult_out)
}

`setCdata` <-
function(data, data_o = NULL, tdata = NULL, settings, init_fun = list(getCdataCC), prefix =""){
   tdata_list <- list()
   settings <- settings[row.names(settings)[!is.na(settings[,"oddGroup"])],]
   # BACKUP DATA INTO DATA_O 
   if(is.null(data_o) )
      data_o <- data
   # DO INIT FUNCTIONS OF THE DATA
   for(fun in init_fun)
      data <- fun(data,settings)
   # TRANSFORM THE DATA
   data <- data[,row.names(settings)]
   t_result <- getCdata(data,settings,tdata)
   data_out <- data.matrix(t_result[["data"]])
   tdata_out <- t_result[["tdata"]]
   baseDir <- paste2(today(),"_",prefix)
   lIdx <- 1
   dirTmp <- paste2(baseDir,'-',letters[lIdx])
   while(!dir.create(dirTmp, showWarnings=FALSE)){
      lIdx <- lIdx+1
      dirTmp <- paste2(baseDir,'-',letters[lIdx])
      if(lIdx > length(letters))
         break
   }
   baseDir <- dirTmp
   figDir <- paste2(baseDir,'/figures')
   tabDir <- paste2(baseDir,'/tables')
   dir.create(figDir, showWarnings=FALSE)
   dir.create(tabDir, showWarnings=FALSE)
   cdata_out <- structure(data_out, data_o  = data_o  , tdata   = tdata_out ,
            settings= settings , xlim    = c(-3,3) , ylim    = c(0,50) , mai
            = c(0.6,0.3,0.05,0.05) , prefix = prefix, baseDir = baseDir ,
            figDir  = figDir , tabDir  = tabDir , rinfo = sessionInfo(), class
            = "cdata")
   return(cdata_out)
}

`plotParcoord` <- 
function(xlim = c(-3,3)
      , cex = NULL
      , title = NULL
      , T_s = NULL
      , pattern = NULL){
   return(function(cmodel){
      T_s <- as.matrix(T_s)
      min2 <- function(x){if(min(x)<0){return(min(x))}else{return(0)}}
      max2 <- function(x){if(max(x)>10){return(max(x))}else{return(10)}}
      ylim<- c(min2(as.numeric(T_s[,"vParY"])),max2(as.numeric(T_s[,"vParY"])))
      plot(x=0 ,new=TRUE ,ann=FALSE ,pch=18,col="white",axes=FALSE ,xlim=xlim , ylim=ylim)
      title(main = paste2("(",paste(as.matrix(attr(cmodel, "cfun_settings")),collapse=","),") ",title), cex=cex)
      for(x_name in pattern){
         pattern <- attr(cmodel,"pattern")[[x_name]][,row.names(T_s)]
         lwd <- 3
         lty <- "solid"
         if(x_name != "mean") {
            lwd <- 1
            lty <- "dashed"
         }
         for(g in 1:nrow(pattern)){
            D_i_s <- cbind(X=as.numeric(pattern[g,]),Y=as.numeric(T_s[,"vParY"]))
            D_i_s <- D_i_s[sort.list(D_i_s[,"Y"]),]
            gap <- 1.05*abs(median(D_i_s[1:(nrow(D_i_s)-1),"Y"] - D_i_s[2:nrow(D_i_s),"Y"]))
            for(l in 1:(nrow(D_i_s)-1))
               if(!(D_i_s[l+1,"Y"] - D_i_s[l,"Y"] > gap))
                  arrows(D_i_s[l,1], D_i_s[l,2], D_i_s[l+1,1], D_i_s[l+1,2], col=attr(cmodel,"cluster_colors")[g], length=0, lwd=lwd, lty=lty)
               # ELSE, AS (is_white_gap == TRUE) THEN DO NOT DRAW ANY ARROW...
         }
      }
      fv_idx <- 1:nrow(T_s)
      if(nrow(T_s) > 20)
         fv_idx <- as.integer(quantile(1:nrow(T_s)))
      axis(2,at=as.numeric(T_s[fv_idx,"vParY"]),labels=colnames(pattern)[fv_idx],las=2,tick=FALSE)
      axis(1,at=seq(from = xlim[1] , to = xlim[2], by = (xlim[2] - xlim[1])/4),cex=cex)
   })
} 

`getPlotFun` <- function(type="plotParcoord", title=NULL, xlim=c(-3, 3), zlim=c(-2,2),
   xy=c(-2.2, 0), cex=0.7, pattern="mean",
   color_gradient=rev(brewer.pal(9,"RdBu")),range_fv=NULL){
#   color_gradient=colorRampPalette(rev(brewer.pal(11,"RdBu")))(8)){
   return(function(cdata){
      return(getPlotFun2(type=type, cdata=cdata, xlim=xlim, xy=xy, cex=cex,
         title=title, pattern=pattern, color_gradient=color_gradient,range_fv=range_fv, zlim=zlim))
      })
}

`getPlotFun2` <- function(type, cdata, title, xlim, xy, cex, pattern,
   color_gradient,range_fv,zlim){
   if(type == "plotParcoord"){
      cfun_parcoord     <- c()
      for(v in unique(as.character(attr(cdata,"settings")[,"vParGroup"]))){
         if(!is.na(v)){
            T_s            <- na.omit(attr(cdata,"settings")[attr(cdata,"settings")[,"vParGroup"] == v,])
            lfun           <- eval(call("plotParcoord",xlim = xlim , cex = cex , title = v, T_s = T_s , pattern = pattern))
            cfun_parcoord  <- c(cfun_parcoord, lfun)
         }
      }
      return(cfun_parcoord)
   }
   if(type == "plot_series"){ 
      if(is.null(range_fv))
         range_fv <- 1:ncol(cdata)
      range_fv <- colnames(cdata)[range_fv]
      return(
         function(cmodel, apattern = pattern, acex = cex, axlim=xlim,rangeFeatVector=range_fv){
               x <- attr(cmodel,"pattern")[[apattern]][,rangeFeatVector]
               plot(ts(t(x)), plot.type='single', col=attr(cmodel,"cluster_colors"), lwd=2, axes=FALSE, 
                  ylim=axlim, xlab='',ylab='', main=paste(as.matrix(attr(cmodel, "cfun_settings")),sep=',',collapse=','))
               fv_idx <- 1:ncol(x)
               if(ncol(x) > 20)
                  fv_idx <- as.integer(quantile(fv_idx))
               axis(1, at=(1L:ncol(x))[fv_idx], labels = colnames(x)[fv_idx], las = 2, line = -0.5, tick = TRUE, cex.axis = acex)
               axis(2, at=as.integer(seq(axlim[1],axlim[2],length.out=4)),las=2)
               d <- table(attr(cmodel,'model')$labelling)[row.names(x)]
               d <- apply(cbind(names(d),' (',d,')'),1,paste2)
               legend("bottomleft",legend=d,col=attr(cmodel,"cluster_colors"),lwd=2, bty='n')
            })
            }
   if(type == "plot_image")            
      return(
         function(cmodel, apattern = pattern, acolgrad = color_gradient, atitle = title, acex = cex, azlim=zlim){
               x <- attr(cmodel,"pattern")[[apattern]]
               image(1L:nrow(x), 1L:ncol(x), x, xlim = 0.5 + c(0, nrow(x)),
                  ylim = 0.5 + c(0, ncol(x)), axes = FALSE, xlab = "", ylab = "",
                  col=acolgrad, main=atitle,add=FALSE, zlim=azlim)
               fv_idx <- 1:ncol(x)
               if(ncol(x) > 20)
                  fv_idx <- as.integer(quantile(fv_idx))
               axis(2, at=(1L:ncol(x))[fv_idx], labels = colnames(x)[fv_idx], las = 2, line = -0.5, tick = 0, cex.axis = acex)
               axis(1, 1L:nrow(x), labels = row.names(x),las = 1, line = -0.5, tick = 0, cex.axis = acex)
            })
   if(type == "plot_legend")           
      return(
         function(cmodel, axy = xy){
            pattern              <- attr(cmodel, "pattern")[["class_count"]]
            for(g in 1:nrow(pattern)){
               g_name      <- attr(cmodel,"cluster_colors")[g]
               y1 <- y0    <- axy[2] + 5 * g / nrow(pattern)
               x0          <- axy[1]
               x1          <- 1.2 * axy[1]
               arrows(x0, y0, x1, y1, col=attr(cmodel,"cluster_colors")[g], length=0, lwd=3)
               text(x0+0.15, y1, labels=paste2(names(pattern)[g], " (",pattern[g],")"), pos=4)
            }
         })
   if(type == "plot_dendro_cluster" | type == "plot_dendro_var")   
        return( 
           function(cmodel
                ,  acex = cex
                , atitle = title ){
                     if(type == "plot_dendro_var") 
                        x <- attr(cmodel,"dendro_var")
                     else
                        x <- attr(cmodel,"dendro_cluster")
                     plot(x, axes = FALSE, yaxs = "i", main=atitle, ylab=NULL,cex=acex)
           })
}

`setTdata` <-
function(testimate=NULL,var=NULL){
  tdata_out <- structure(
       testimate
     , var   = var
     , class = "tdata")
  return(tdata_out)
}

`funStatsOddratios` <-
function(data,class,fun_midthreshold=median){
   #
   nclass <- ncol(class)
   class <- map(class, warn=FALSE)
   ldata <- cbind(data,class=class)
   sgroup <- as.data.frame(na.omit(attr(data,"settings")[,c("vParGroup","oddGroup")]))
   s <- matrix(0,nclass, length(unique(sgroup[,"oddGroup"])),
                dimnames=list(1:nclass, sort(unique(sgroup[,"oddGroup"]))))
   class_vector <- 1:nclass
   gr_vector <- sort(unique(as.character(sgroup[,"oddGroup"])))
   nrow2 <- function(x){y <- nrow(x) ; if(is.null(y)){return(0)}else{return(y)}}
   sum2 <- function(x){if(!is.null(dim(x))){return(apply(x,1,sum))}else{return(x)}}
   # FOR EACH HIERARCHICAL SUBSET OF OUTCOMES, MAKE SUM SCORES
   for(gr in gr_vector){
      m <- cbind(SScore=sum2(ldata[,row.names(sgroup[sgroup[,"oddGroup"]==gr,])]),class=class)
      # FOR EACH PATTERN, DERIVE ITS STATISTIC
      for(g in class_vector){
         middleScore <- fun_midthreshold(m[,"SScore"],na.rm=TRUE)
         m11 <- nrow2(m[m[,"class"]==g & m[,"SScore"] >= middleScore,])
         m12 <- nrow2(m[m[,"class"]==g & m[,"SScore"] < middleScore,])
         m21 <- nrow2(m[m[,"class"]!=g & m[,"SScore"] >= middleScore,])
         m22 <- nrow2(m[m[,"class"]!=g & m[,"SScore"] < middleScore,])
         mat <- matrix(c(m11,m12,m21,m22),2,2,byrow=TRUE)
         # LOG OF THE ODD-RATIO ALSO NAMED CROSS PRODUCT
         s[g,match(gr,gr_vector)] <- log(mat[1,1] * mat[2,2] / (mat[1,2] * mat[2,1]))
         # dimnames(s)[[2]][match(gr,gr_vector)] <- paste(gr,sprintf("%.1f",middleScore),sep="_",collapse="")
      }
   }
   return(s)
}

`funStatsChi2test` <- 
function(data, class, class_col='Class', probs=seq(0,1,0.4)[2:3]){
   ldata <- cbind(data, modelClass=map(class))
   ldata <- cbind(ldata, target=as.character(attr(data,'data_o')[row.names(ldata),class_col]))
   s_test <- chisq.test(xtabs( ~.,ldata[,c('modelClass',"target")]), simulate.p.value=TRUE)
   s <- cbind(s_test[["residuals"]]^2, chi2test=NA)
   s[1:2,"chi2test"] <- c(sprintf('%.1e',s_test[["p.value"]]), sprintf('%.1f',s_test[["statistic"]]))
   colnames(s)[match("chi2test",colnames(s))] <- paste2(class_col, " (residual)")
   return(s)
}

`funStatsJointDistrib` <- 
function(data, class, class_col='Class'){
   ldata <- cbind(data, modelClass=map(class))
   ldata <- cbind(ldata, target=as.character(attr(data,'data_o')[row.names(ldata),class_col]))
   s <- xtabs( ~.,ldata[,c("modelClass","target")])
   return(s)
}

`splitTraintestStrata` <-
function(data,K=10,perc=0.7,rseed = 6013){
   indexes <- list()
   G <- length(table(data[,"class"]))
   for(g in 1:G){
      strata_idx <- nrow(data[data[,"class"] == g,])
      indexes[[g]] <- matrix(0,0,strata_idx * (1-perc)) 
      set.seed(rseed)
      for(k in 1:K)
            indexes[[g]] <- rbind(indexes[[g]],sample(1:strata_idx)[1:(strata_idx * (1-perc))])
   }
   return(list(indexes=indexes,K=K))
}

`today` <-
function(){
   tmp_date <- format(Sys.time(), "%Y-%m-%d")
   return(tmp_date)
}

`tAbsMax` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun         <- list(function(x){return(max(abs(x),na.rm=TRUE))}
                      , function(x){return(do.call('/',x))})
   return(doTdata(data=vdata,tfun=tfun,tdata=tdata))
}

`doTdata` <-
function(data,tfun=list(function(x){return(0)},function(x){return(do.call('+',x))}),tdata){
   if(is.null(attr(tdata,"tfun"))) 
      testimate         <- tfun[[1]](data)
   else 
      testimate         <- as.numeric(tdata)
   vdata                <- tfun[[2]](list(data,testimate))
   tdata_out            <- setTdata(testimate=testimate,var=attr(tdata,"var"))
   return(list(data=vdata,tdata=tdata_out))
}


`tAvg` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun         <- list(function(x){return(mean(x,na.rm=TRUE))}
                      , function(x){return(do.call('-',x))})
   return(doTdata(data=vdata,tfun=tfun,tdata=tdata))
}

`tL1` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun      <- list(function(x){return(x/sqrt(sum(x,na.rm=TRUE)))}
                      , function(x){return(do.call('/',x))})
   return(doTdata(data=vdata,tfun=tfun,tdata=tdata))
}

`tL2` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun      <- list(function(x){return(x/sqrt(sum(x^2,na.rm=TRUE)))}
                      , function(x){return(do.call('/',x))})
   return(doTdata(data=vdata,tfun=tfun,tdata=tdata))
}

`tMax` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun         <- list(function(x){return(max(x,na.rm=TRUE))}
                      , function(x){return(do.call('-', x))})
   return(doTdata(data=vdata,tfun=tfun,tdata=tdata))
}

`tMedian` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun         <- list(function(x){return(median(x,na.rm=TRUE))}
                      , function(x){return(do.call('-',x))})
   return(doTdata(data=vdata,tfun=tfun,tdata=tdata))
}

`tMin` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun       <- list(function(x){return(min(x,na.rm=TRUE))}
                      , function(x){return(do.call('-',x))})
   return(doTdata(data=vdata,tfun=tfun,tdata=tdata))
}

`tSigma` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun         <- list(function(x){sd_est <- sd(x,na.rm=TRUE) ; if(sd_est>0) return(sd_est) else return(1)}
                      , function(x){return(do.call('/', x))})
   return(doTdata(data=vdata,tfun=tfun,tdata=tdata))
}


`tAdjust` <-
function(data, tdata, tformula){
   var          <- strsplit(tformula,"~")[[1]][1]
   data         <- as.data.frame(data)
   if(length(tdata) == 0){
      print(tformula)
      tdata_lm          <- lm(as.formula(tformula),data=data)
      tdata_lm[["print"]]<- tformula
      tdata             <- setTdata(tdata_lm,var=attr(tdata,'var'))
      # AS THE VARS WHICH ARE NOT INVOLVED IN THE CANALYSIS MAY HAVE MISSING
      # VALUES, TO PRESERVE THE NA'S WE REBUILD A MATRIX DATA STRUCTURE WHICH
      # SERVES TO INDEX PROPERLY THE RESIDUALS. AND WE RETURN A SIMPLE NUMERIC 
      # VECTOR WHICH HAS SOME NA'S IN IT.
      vdata             <- matrix(NA,nrow(data),1,dimnames=list(row.names(data),var))
      tresiduals        <- residuals(tdata)
      vdata[names(tresiduals),] <- tresiduals
      vdata             <- as.numeric(vdata)
   }
   else
      vdata             <- data[,var]-predict(tdata,data)
   return(list(data=vdata,tdata=tdata))
}

`getCdata` <-
function(data,settings,tdata=NULL){
   # BEFORE TRANSFORMING THE DATA, CONVERT IT TO NUMERIC
   data <- data.matrix(data)
   if(is.null(tdata)) tdata <- list()
   settings <- as.matrix(settings)
   # FOR EACH VAR
   for(var in row.names(settings)){
      tdata_idx                 <- which(names(tdata) == var)
      if(length(tdata_idx) == 0) 
         tdata[[var]]  <- list()
      # USE (MULTIPLE)-SPACES AS SEPARATOR BETWEEN TRANSFORMATION FUNCTIONS
      tFun             <- strsplit(settings[var,"tFun"],"[ ]+")[[1]]
      for(lfun in tFun){
         if(is.na(lfun)) 
            break ;
         # IN CASE 'LFUN' AS AN INSIDE PARAMETER, RETRIEVE IT (1) REPLACE ALL
         # ('\") BY SPACES AND THEN (2) SPLIT
         lfun                   <- gsub(" ","",lfun,extended=FALSE)
         lfun                   <- gsub("\"","'",lfun,extended=FALSE)
         lfun                   <- strsplit(gsub("')","",gsub("('"," ",lfun,extended=FALSE),extended=FALSE),"[ ]+")[[1]]
         if(length(tdata_idx) == 0)
            tdata_var_lfun      <- setTdata(testimate=NULL,var=var)
         else
            tdata_var_lfun      <- tdata[[var]][[lfun[1]]]
         arg_list               <- list(data,tdata=tdata_var_lfun)
         # IF LFUN HAS SOME ARGUMENTS, WE APPEND THEM TO THE ARGUMENT LIST
         if(!is.na(lfun[2]))
            arg_list            <- append(arg_list,lfun[2])
         # PROCEED TO THE TRANSFORM AND OUTPUT A DATA VECTOR AND A TDATA
         var_result             <- do.call(lfun[1],args=arg_list)
         # UPDATE THE MAIN DATA FRAME WITH THE VECTOR
         data[,var]             <- var_result[["data"]]
         # KEEP RECORD OF EACH INDIVIDUAL TRANSFORMATION FOR THE 'VAR'
         tdata[[var]][[lfun[1]]] <- var_result[["tdata"]]
      }
   }
   return(list(data=data,tdata=tdata))
}

`doChi2Testing` <- 
function(m,tonominal=FALSE,simulate.p.value = TRUE, B = 10000){
   m2 <- m
   if(tonominal)
      m2 <- cbind(m,sup_med=as.numeric(m[,1])>median(as.numeric(m[,1])))[,3:2]
   a <- xtabs(formula=~.,data=m2)
   if(nrow(a) > 2) range_row <- 1:nrow(a)
   else range_row <- 1
   apvalues <- matrix(NA,max(range_row),ncol(a))
   for(i in range_row){
      for(j in 1:ncol(a)){
         a2by2 <- matrix(c(a[i,j],sum(a[i,][-j]),sum(a[,j][-i]),sum(a[-i,-j])),2,2,byrow=TRUE)
         # FOR TESTING: cat("\n",i,j) print(a2by2)
         chi2_test <- chisq.test(a2by2,simulate.p.value = simulate.p.value, B = B)
         apvalues[i,j] <- sign(chi2_test$residuals[1,1]) * chi2_test$p.value
      }
   }
   if(nrow(apvalues) > 1)
      row.names(apvalues) <- row.names(a)
   else if(nrow(apvalues) == 1 && tonominal) 
      row.names(apvalues)[1] <- paste2('> med ',colnames(m)[1])
   else 
      row.names(apvalues)[1] <- row.names(a)[1]
   return(apvalues)
}



`plot_patternAndPvalues` <-
function(a, title="", ylim=c(-2,38), ylim2=c(15,0), xlab=""){
   # P-VALUES
   plot(x=0, ann=FALSE, pch=18, col="white", axes=FALSE, xlim=c(1,nrow(a)), ylim=ylim2)
   lines(list(x=1:nrow(a),y=-log(a[,'pvalues'])), col='grey', lty=1, lwd=1)
   lines(list(x=c(5,nrow(a)), y=-log(c(0.05,0.05))),col='grey',lwd=1)
   lines(list(x=c(5,nrow(a)), y=-log(c(0.01,0.01))),col='grey',lwd=1)
   axis(4, at=-log(c(1,0.05,0.01)), las='1', labels=c('1   ',' .05',' .01'))
   # IMAGE
   b <- matrix(0,nrow(a),1)
   b[which(abs(a[,'pvalues'])<0.05),] <- 0.5
   b[which(abs(a[,'pvalues'])<0.01),] <- 1
   image(x=1:nrow(a),y=-log(c(1,0.05)),z=b,col=brewer.pal(5,"Blues"), breaks=c(-0.1,0.1,0.3,0.6,0.8,1),add=TRUE)
   par(new=TRUE)
   plot(x=0, new=TRUE, ann=FALSE, pch=18, col="white", axes=FALSE, xlim=c(1,nrow(a)), ylim=ylim,xlab=xlab)
   # SCORES
   lines(list(x=1:nrow(a), y=a[,'scores']), lwd=2)
   # RANGES OF FREQUENCIES: HIGHLY DISCRIMINATIVES
   topFeats <- as.numeric(which(abs(a[,'pvalues'])<0.05))
   k <- j <- topFeats[1]
   for(i in topFeats){
      if(i>j+1) k <- c(k,j,i)
      j<- i
   }
   if(!is.integer(length(k)/2))
   k <- c(k,k[length(k)])
   rangeFeats <- matrix(k,length(k)/2,2,byrow=TRUE)
   # SCORE AXES  
   axis(1, at=c(1,unique(rangeFeats),nrow(a)), labels=row.names(a)[c(1,unique(rangeFeats),nrow(a))],las='2')
   seqY <- seq(ylim[1],ylim[2],length.out=4)
   axis(2, at=seqY, labels=sprintf('%1.1e',seqY), las=2)
   title(title)
}


`bootstrap_pvalues` <- 
function(x, B_max=1000,rseed=6013){
   set.seed(rseed)
   master_bootstrap <- sample(1:nrow(x),B_max*nrow(x),replace=TRUE)
   x2 <- matrix(NA,B_max,ncol(x)-1,dimnames=list(list(),colnames(x[,-ncol(x)])))
   for(b in 1:B_max){
      idx <- master_bootstrap[(1+(b-1)*(nrow(x)-1)):(b*(nrow(x)-1))]
      cat(".",b)
      for(f in 1:(ncol(x)-1)){
         x2[b,f] <- chisq.test(xtabs(formula=~.,data=x[idx,c(f,ncol(x))]),simulate.p.value=TRUE, B = 1000)$p.value
         }

   }
   return(x2)
}


`toBinary` <- 
function(x){
   x_bin <- x
   for(i in 1:ncol(x))
      x_bin[,i] <- as.numeric(x[,i])>median(as.numeric(x[,i]))
   return(x_bin)
}



`plot.test.feats` <- 
function(x, range_fv=NULL, p_quantile=0.95, ylim=NA, ylim2=c(15,0), xlab="", text="", ...){
   postscript(paste2(attr(x,'prefix'),'.ps'))
   par('mfrow'=c(3,2),mai=c('0.5','0.5','0.5','1'))
   for(i in 1:length(attr(x,'fv'))){
      if(length(range_fv)>0)
         a <- cbind(scores=attr(x,'sc_med')[[i]][range_fv], pvalues=apply(attr(x,'x2')[[i]][,range_fv],2,quantile,probs=p_quantile))
      else
         a <- cbind(scores=attr(x,'sc_med')[[i]], pvalues=apply(attr(x,'x2')[[i]],2,quantile,probs=p_quantile))
      if(length(ylim) != 2)
         ylim <- as.integer(range(a[,'scores']))
      title <- paste2(text, ' (#',i,')')
      plot_patternAndPvalues(a, title=title, ylim2=ylim2, ylim=ylim, xlab=xlab)
   }
   graphics.off()
}

`test.feats` <- 
function(x, q=NA, B=4){
   if(class(x) == "cresult"){
      # RETRIEVE APPROPRIATE CLASS/SUBTYPE LABELS 
      if(!is.na(q) && length(which(q %in% names(x))) == 0)
         l <- as.character(attr(attr(x,'cdata'),'data_o')[row.names(attr(x,'cdata')),q])
      else{
         if(is.na(q))
            q <- getBestModel(x,1)
         l <- as.character(attr(x[[q]],'model')$labelling)
      }
      m <- data.frame(attr(attr(x,'cdata'),'data_o')[,colnames(attr(x,'cdata'))],class=l)
      p <- attr(x, 'prefix')
   }
   else{
      p <- today()
      m <- x 
      l <- as.character(x[,'class'])
   }
   # INIT DATA
   prefix <- paste2(p,'_discriminative_features_',gsub(",","_",q))
   mBin <- data.frame(toBinary(m[,-ncol(m)]),class=l)
   x2 <- fv <- sc_med <- list()
   # 
   for(i in sort(unique(l))){
      mTmp <- data.frame(mBin[,-ncol(mBin)], class=(l==i))
      x2[[i]] <- bootstrap_pvalues(mTmp, B_max=B)
      fv[[i]] <- names(which(apply(x2[[i]], 2, quantile, probs=0.5)<0.05))
      sc_med[[i]] <- apply(m[(l==i), -which(colnames(m) %in% c('class',q))],2,median)
   }
   obj_out <- structure(fv, x2=x2, fv=fv, sc_med=sc_med, prefix=prefix, class="test.feats")
   if((class(x) == "cresult") && (length(which(q %in% names(x)))>0)){
      attr(x[[q]],'test.feats') <- obj_out
      return(x)
   }
   else
      return(obj_out)
}

`bestNFeatures` <- 
function(x, N=10){
   # FEATURE SELECTION BY RANKING OF THE P-VALUES
   cList <- names(attr(x,'x2'))
   m <- attr(x,'x2')[[1]]
   pTable <- rTable <- fSelFreq <- matrix(NA, ncol(m),length(cList), dimnames=list(colnames(m),cList)) 
   maxP <- matrix(NA,length(cList),1,dimnames=list(cList,'maxP'))
   # NFREQ CONTROLS THE QUANTITY OF FEATURE SELECTED BY CLASS
   nPerClass <- as.integer(N/length(cList))
   fList <- NULL
   for(i in cList){
      pTable[,i] <- abs(attr(x,'sc_med')[[i]])
      rTable[,i] <- order(attr(x,'sc_med')[[i]])
      fSelFreq[,i] <- order(attr(x,'sc_med')[[i]])<=nPerClass
      fList <- unique(c(fList,names(which(fSelFreq[,i]))))
   }
   j <- 1
   while(length(fList)<N){
      i <- 1+j%%length(cList)
      fSelFreq[,i] <- order(attr(x,'sc_med')[[i]]) <= (nPerClass+as.integer(1+j/length(cList)))
      fList <- unique(c(fList,names(which(fSelFreq[,i]))))
      j <- j+1
   }
   for(i in cList)
      maxP[i]<- max(pTable[which(fSelFreq[,i]),i])
   attr(x,'bestNFeatures') <- structure(fList, N=N, pTable=pTable, rTable=rTable, maxP=maxP,class='bestNFeatures')
   return(x)
}


`varFeatAnalysis` <- 
function(m, a=NULL, maxFeat=NA, postscript=TRUE,clusterModel='VVI', clusterG=6,
   clusterRseed=6013:6063){
   nFeatMax <- min(ncol(m)-1,maxFeat,na.rm=TRUE)
   iterMax <- as.integer((2*(log(nFeatMax))/log(2)-4))
   iterSeq <- as.integer(2^(2+(1:iterMax)/2))
   x <- NULL
   if(is.null(a))
      a <- test.feats(m)
   for(i in iterSeq){
      cat('\n',i,' features...')
      aTmp <- bestNFeatures(a,N=i)
      # PREPARE SETTINGS
      settings <- genCdataSettings(m)
      settings[,'inCAnalysis'] <- FALSE 
      settings[as.character(attr(aTmp,'bestNFeatures')),'inCAnalysis'] <- TRUE
      # PREPARE CDATA
      cdata <- setCdata(m,settings=settings, prefix=paste2('With_',i,'_Features'))
      # PREPARE CRESULT
      x[[as.character(i)]] <- setCresult(cdata=cdata,
         fun_stats=list(), fun_plot=list(), fun_pattern=NULL, nTopModels=5,
         cfun=clusterModelBasedEM, cfun_settings=list(modelName=clusterModel,
         G=clusterG, rseed=clusterRseed))
      x[[as.character(i)]] <- doModeling(x[[as.character(i)]])
   }
   # V, SCALE UP
   combi <- t(combn(names(x),2))
   compList <- list()
   v<- matrix(NA,length(x),length(x),dimnames=list(names(x),names(x)))
   for(i in 1:nrow(combi)){
       idx <- paste(combi[i,],collapse=",")
       compList[[idx]] <- compareCresult(x[[combi[i,1]]],x[[combi[i,2]]])
       v[combi[i,2],combi[i,1]] <- v[combi[i,1],combi[i,2]] <- compList[[idx]]$summary[1]
    }
    print(v)
    # GRAPHIC OUTPUT
   fileName <- paste2(today(),'varFeatAnalysis')
   if(postscript){
      postscript(file=paste2(fileName,'.ps'))
      image(1L:nrow(v), 1L:ncol(v), v, xlim = 0.5 + c(0, nrow(v)),ylim = 0.5 +
         c(0, ncol(v)),    axes = FALSE, xlab = "", ylab = "",
         col=brewer.pal(9,"Greys"), main='scaling up
         (V)',zlim=range(v,na.rm=TRUE))
      fv_idx <- 1:ncol(v)
      if(ncol(v) > 20) 
         fv_idx <- as.integer(quantile(fv_idx))
         axis(2, at=(1L:ncol(v))[fv_idx], labels = colnames(v)[fv_idx], las = 2, line = -0.5, tick = 0)
         axis(1, 1L:nrow(v), labels = row.names(v),las = 1, line = -0.5, tick = 0)
         graphics.off()
   }
   # SAVE
   cat("\nSave 'cresult' into ", fileName)
   save(list = "x", file = paste2(fileName,'.RData'))
   return(x)
}

`getCvar` <- 
function(x){
   varSel <- which(x[,"inCAnalysis"] == "TRUE" | x[,"inCAnalysis"] == " TRUE")
   return(row.names(x)[varSel])
}

`print.bicanalysis` <-
function(x, ...){
   return(x[1:nrow(x),])
}

`print.cdata` <-
function(x, ...){
   return(x[1:nrow(x),getCvar(attr(x,"settings"))])
}

`print.cmodel` <-
function(x, ...){
   xModel <- attr(x, "model")
   m <- cbind(xModel[["z"]], class=map(xModel[["z"]]))
   colnames(m)[1:(ncol(m)-1)] <- 1:(ncol(m)-1)
   return(m)
}

`print.cresult` <- 
function(x, m1=NULL, m2=NULL, type='term', label='tab:jointdistrib', ...) {
   doCaption <- function(o, str){
      str <- gsub('x','and',gsub('[()]', ' ', gsub('_',',',str)))
      str <- paste2("The level of association of models", str,"is $V=", sprintf('%.1f', 100*o["V",1]))
      str <- paste2(str, "\\%$ and the $\\chi^2$ test, on which $V$ is based, has its $p=")
      str <- paste2(str, sprintf('%.4f', o['CHI2',1]), '$ ($\\chi^2=',sprintf('%.1f',o['CHI2',2]),'$).')
      return(str)
   }
   if(is.null(m1) & is.null(m2)){
      d <- compareModel(x)
      for(i in names(d))
         print(xtable(d[[i]][1:(nrow(d[[i]])-2),],caption=doCaption(d[[i]],i),digits=0,label=paste2(label,i), type=type))
   }
   else if((m1 %in% names(x)) & (m2 %in% names(x))){
      m12 <- compareModel(x[[m1]], x[[m2]])
      str <- paste2('(',m1,')x(',m2,')')
      if(type == 'latex' | type == 'html')
         print(xtable(m12[1:(nrow(m12)-2),],caption=doCaption(m12,str),digits=0,label=label), type=type)
      else
         return(m12)
   }
}

`summary.bicanalysis` <- 
function(object, fun='mean', bic='relativeBic', type='latex', caption=NULL, label='tab:bic', fmt='%.2f', ...){
   x.summary <- getBicStats(attr(object,'data'), fun=fun, vAgg=bic, fmt=fmt)
   if(tolower(type)=='latex' | tolower(type) == 'html'){
      if(is.null(caption)) 
         caption <- paste2(fun,' of ', gsub('_',' ',bic))
      print(xtable(x.summary, caption=caption, label=label,type=type))
   }
   else
      print(as.matrix(x.summary))
   return(x.summary)
}


`summary.cdata` <- 
function(object, shift='', ...){
   cat(shift, '\'',attr(object,'prefix'),'\' dataset for subtype discovery analysis\n')
   cat(shift, attr(object,'rinfo')[['R.version']]$version.string,',',attr(object,'rinfo')[['R.version']]$platform,'\n')
   cat(shift, 'SubtypeDiscovery ', attr(object,'rinfo')$otherPkgs[['SubtypeDiscovery']][['Version']],'\n')
   cat(shift, '\t number of variables in the cluster analysis/originally: ', ncol(print(object)),'/', ncol(attr(object,'data_o')),'\n')
   cat(shift, '\t number of complete and incomplete cases: ', nrow(print(object)), '/', nrow(attr(object,'data_o'))-nrow(print(object)), '\n')
   cat(shift, '\t xlim: ',attr(object,'xlim'), '\n')
   cat(shift, '\t ylim: ',attr(object,'ylim'), '\n')
   cat(shift, '\t par() \'mai\' parameter: ',attr(object,'mai'),'(to define the base margin of the plotting regions)\n')
   cat(shift, '\t base, figure and table directories:\n\t\t', attr(object,'baseDir'),'\n\t\t', attr(object,'figDir'),'\n\t\t', attr(object,'tabDir'), '\n')
}

`summary.cmodel` <- 
function(object, type='term', label='tab:stats', caption=NULL, ...){
   objStats <- attr(object,"stats")
   m <- attr(object,'cfun_settings')
   mCaption <- ''
   summaryStr <- function(str){
      cList <- strsplit(strsplit(gsub('and[ ]*','',str),' ')[[1]],'')
      str <- ''
      for(i in 1:length(cList))
         str <- paste2(str,cList[[i]][1],'.')
      return(toupper(str))
   }
   if(length(objStats) >=1){
      mCaption <- paste2('Statistics of model ', m[[1]], ',', m[[2]], ',', m[[3]],' (')
      x <- matrix(0,0,nrow(objStats[[1]]))
      i <- 1
      while(i <= length(objStats)){
         x <- rbind(x, t(objStats[[i]]))
         mCaption <- paste2(mCaption, names(objStats)[[i]],', ')
         i <- i+1
      }
      mCaption <- gsub(', )',')',paste2(mCaption,').'))
      row.names(x) <- sapply(row.names(x),summaryStr)
   }
   else
      x <- 'no stats'
   if(tolower(type)=='latex' | tolower(type) == 'html')
      print(xtable(x, caption=mCaption, label=label, type=type))
   else 
      print(x)
}

`summary.cresult` <- 
function(object, type='term', ...){
   cat('\'',attr(object,'prefix'),'\' subtype discovery analysis summary\n')
   cat('------- data ----------\n')
   summary(attr(object,'cdata'), shift=' ')
   cat('----- settings --------\n')
   cat(' ',attr(object,'rinfo')[['R.version']]$version.string,',',attr(object,'rinfo')[['R.version']]$platform,'\n')
   cat(' SubtypeDiscovery ', attr(object,'rinfo')$otherPkgs[['SubtypeDiscovery']][['Version']],')\n')
   cat(rep(' ',3),'model based cluster analysis\n')
   p <- as.matrix(attr(object,'cfun_params'))
   cat(rep(' ',6),'covariance models: ', sort(unique(p[,1])),'\n')
   cat(rep(' ',6),'number of clusters: ', sort(unique(p[,2])),'\n')
   cat(rep(' ',6),'random initialization numbers: ',sort(unique(p[,3])),'\n')
   cat(rep(' ',6),'BIC table of relative difference with respect to most likely model\n')
   x.summary <- summary(attr(object,'bicanalysis'),type=type,bic='relativeBic')
   x.summary <- summary(attr(object,'bicanalysis'),type=type,bic='BIC')
   x.summary <- summary(attr(object,'bicanalysis'),type=type,bic='BIC',fun='rank')
   vBModels <- getBestModel(object)
   cat(rep(' ',3),'model ranking, top: \n')
   for(i in 1:length(vBModels))
      cat('\t\t', i,' ', vBModels[i],'\n')
   if(is.null(q))
      q <- vBModels[1]
}


