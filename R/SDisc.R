`cat2` <- 
function(str, ...){
   cat('\n',paste2(str, ...),'\n')
}

`SDStability.default` <-
function(x, q, rseed=6013, nnoise=10, nrep=10){
   settings <- strsplit(q,',')[[1]]
   df <- print(SDData(x))
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
         Y <- df+matrix(rnorm(nrow(df)*ncol(df),sd=sd_l,mean=0),nrow(df),ncol(df))
         #r[[l]][[s]] <- modelBasedEM(Y,modelName=settings[1],G=settings[2],rseed=as.numeric(settings[3]))
         #r[[l]][[s]][["Y"]] <- Y
      }
      tmp <- c()
      for(ridx1 in 1:nrep){
         s1 <- rseed+ridx1
         for(ridx2 in 1:nrep){
            s2 <- rseed+ridx2
            if(ridx2>ridx1)
               tmp <- c(tmp, agreementScores(r[[l]][[s1]]$label, r[[l]][[s2]]$label))[["V"]]
         }
      }
      v[l,] <- sort(tmp)
   }
   v <- structure(v,class='SDStability')
   return(v)
}

`plot.SDStability` <- 
function(x, ncolors=9, minmax=c(0.5,1), ...){
   breaks <- seq(minmax[1],minmax[2],(minmax[2]-minmax[1])/ncolors)
   plotInit(plotNbr=1, type="SDStability")
   image(x, col=brewer.pal(ncolors,"Greys"), breaks=breaks, axes=FALSE, ylab='Association level V, quantiles',
      cex.lab=1.7)
   contour(x, add=TRUE,labcex=1.3)
   axis(1, at=seq(0,1,0.11), las=2, cex.axis=1.7, labels=sprintf("%.1f%%",100*as.numeric(row.names(x))))
   axis(2, at=seq(0,1,0.2), cex.axis=1.7, labels=sprintf("%d%%",100*seq(0,1,0.2)))
   title(paste(SDCModelSettings(x),sep="",collapse=" "))
   plotClose()
}

`predict.SDData` <-
function(object, newdata, prefix='Newdata', subset=NULL, ...){
   df <- dataOrig <- newdata <- SDDataDimnames(newdata)
   if(class(newdata) == 'SDisc' || class(newdata) == 'SDData')
      df <- dataOrig <- SDDataOrig(newdata)
   else{
      settings <- SDDataSettings(object) 
      for(fun in SDDataInitFun(object))
         df <- fun(df, settings)
      df <- df[,row.names(settings)]
      if(!is.null(subset))
         df <- df[which(subset %in% row.names(df)),]
      tRes <- list()
      newdf <- df
      for(i in row.names(settings[!is.na(settings[,'tFun']),])){
         res <- SDTData(i, tFun=settings[i,'tFun'], data=newdf, TData=SDTData(object)[[i]])# SDTDataElmt(tEst[[i]]) 
         tRes[[i]] <- res[["tdata"]]
         newdf <- res[["data"]]
      }
      res <- structure(data.matrix(newdf), dataOrig=dataOrig, SDTData=tRes, settings=settings, prefix=prefix,
         class='SDData')
      res <- SDDataSetup(res, prefix=prefix)
      return(res)
   }
}

`predict.SDisc` <-
function(object, newdata, ...){
      #data <- predict(SDData(object), newdata)
   if(class(object) =='SDisc' && class(newdata) =='SDData'){
      data <- SDData(SDDataOrig(newdata), settings=SDDataSettings(newdata), TData=SDTData(object),
         prefix=paste2(SDPrefix(newdata), '-Data-PredictedBy-', SDPrefix(object)), subset=SDDataSubset(newdata)) 
      df <- print(data)
      m <- bestModel(object, 1)
      mName <- strsplit(m, ',')[[1]]
      mOrig <- attr(object[[m]],'model')
      if(!is.null(mOrig$parameters$variance$Sigma))
         mOrig$parameters$variance$Sigma <- SparseM::as.matrix(mOrig$parameters$variance$Sigma)
      denseSigma <- array(NA, dim=c(ncol(df), ncol(df),as.numeric(mName[[2]])))
      for(gId in 1:as.numeric(mName[[2]]))
         denseSigma[,,gId] <- SparseM::as.matrix(mOrig$parameters$variance$sigma[[gId]])
      mPred <- estep(modelName=mOrig$modelName, data=df, parameters=mOrig$parameters)
      mPred[["label"]] <- map(mPred[["z"]],warn=FALSE)
      res <- SDisc(data, cFunSettings=list(modelName=mName[[1]], G=mName[[2]], rseed=mName[[3]]), nTopModels=1,
               nnodes=1)
      attr(res[[1]],'model') <- mPred
      return(res)
   }
   else if(class(object) == 'SDisc' && class(newdata) == 'SDisc')
      predict(object, newdata=SDData(newdata))
   else
      cat2('You must provide two \'SDisc\' objects or an \'SDisc\' and a \'SDData\'. ')
}

`SDPattern` <-
function(m, fun='mean'){
   m <- data.matrix(aggregate(m, function(y){return(do.call(fun, args=list(y, na.rm=TRUE)))},
      by=list(class=m[,"class"])))[,-1]
   row.names(m) <- 1:nrow(m)	
   if(!is.na(match("class", dimnames(m)[[2]])))
      m <- m[,-match("class", dimnames(m)[[2]])]
   return(m)
}

`agreementScores` <- 
function(x, y=NULL){
   if(dim(x)[[2]] == 2 && is.null(y))
      y <- x[,2]
   if(length(x) == length(y)){
      res <- compareMatchedClasses(x, y, method="rowmax", verbose=FALSE)
      xtab <- as.table(ftable(x,y))
      xSummary <- summary(xtab)
      n <- as.numeric(xSummary["n.cases"]) 
      X2 <- as.numeric(xSummary["statistic"]) 
      k <- min(dim(xtab)) 
      res[["V"]] <- sqrt(X2/(n*(k-1)))
      return(res)
   }
}

`compareModel` <-
function(x, y=NULL, fmt=c('%d','%.2f')){
   if(class(x) == 'SDisc' && class(y) == 'SDisc')
      return(compareModel(x[[bestModel(x,1)]], y[[bestModel(y,1)]]))
   if(class(x) == 'SDisc' && is.null(y)){
      mPairs <- t(combn(bestModel(x),m=2))
      cNames <- gsub(",","_",paste("(",apply(mPairs,1,paste,collapse=")x("),")",sep=""))
      m1 <- x[[which(mPairs[1,1] == names(x))]]
      m2 <- x[[which(mPairs[1,2] == names(x))]]
      return(compareModel(m1,m2))
   }
   else if(class(x) == 'SDCModel' & class(y) == 'SDCModel'){
      m1 <- attr(x,"model")
      m2 <- attr(y,"model")
      if( !is.na(m1$loglik) && !is.na(m2$loglik)){
         res <- agreementScores(m1$label, m2$label) 
         m <- as.table(ftable(factor(m1$label, levels=1:m1$G), factor(m2$label, levels=1:m2$G)))
         m <- m[hclust(dist(m))$order, hclust(dist(t(m)))$order]
         #m <- m[,as.numeric(matchClasses(m))]
         m[m == 0] <- NA
         res[["xtab"]] <- m  
         res[["chi2"]] <- chisq.test(ftable(m1$label, m2$label), simulate.p.value=TRUE)
         return(res)
      }
      else{
         cat2('One of the two \'SDCModels\' provided to \'compareModel\' fail to optimize (NA likelihood). ')
         return(NULL)
      }
   }
   else{
      cat2('You must provide one or two \'SDisc\' objects, or two \'SDCModel\'.')
      return(NULL)
   }
}


`modelBasedEM` <-
function(x, data) {
   p <- as.matrix(SDCModelSettings(x))
   G <- as.numeric(p[,'G'])
   rseed <- as.numeric(p[,'rseed'])
   modelName <- as.character(p[,'modelName'])
   cat(modelName,',', G, ',', rseed, ' ',sep="")
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
   mstep_est <- mstep(data=data, modelName=modelName, z=cmembership_matrix, warn=FALSE)
   # AND THEN, STARTS EM WITH THE E-STEP
   model <- em(data=data, modelName=mstep_est$modelName, parameters=mstep_est$parameters, warn=FALSE)
   model[["label"]] <- map(model[["z"]],warn=FALSE)
   model[["z"]] <- NULL 
   # CALCULATE THE BIC SCORE
   model[["BIC"]] <- bic(modelName=modelName, loglik=model$loglik, n=nrow(data) , d=ncol(data), G=G)
   if(!is.na(model$loglik)){
      if(!is.null(model$parameters$variance$Sigma))
         model$parameters$variance$Sigma <- as.matrix.csr(model$parameters$variance$Sigma)
      sparseSigma <- list()
      for(gId in 1:G)
         sparseSigma[[gId]] <- as.matrix.csr(model$parameters$variance$sigma[,,gId])
      model$parameters$variance$sigma <- sparseSigma
      data <- cbind(data, class=model[["label"]]) 
      attr(x,'pattern') <- list(mean=t(model$parameter$mean), clCount=table(data[,'class']))
      attr(x, 'GColors') <- brewer.pal(length(unique(data[,'class'])),"Set1")
   }
   attr(x,'model') <- model
   x <- isModeled(x, set=TRUE)
   return(x)
}

`bestModel` <-
function(x, n=NULL, modelName=NULL, G=NULL){
   if(is.null(n))
      n <- attr(x,"nTopModels")
   m <- print(attr(x,"bicTable"), n=n, modelName=modelName, G=G, latex=FALSE)
   if(n < 1)
      mList <- unique(apply(m[,c("modelName","G")], 1, paste, collapse=",", sep=""))
   else
      mList <- apply(m[,c("modelName","G","rseed")], 1, paste, collapse=",", sep="")
   return(as.character(mList))
}

`SDDataDimnames` <-
function(data, settings){
      if(is.null(row.names(data)))
         row.names(data) <- 1:nrow(data)
      if(is.null(colnames(data)))
         colnames(data) <- paste('v', 1:ncol(data), sep='')
   return(data)
}

`SDDataCC` <-
function(data, settings){
   vNames <- SDDataInCAnalysis(settings)
   tFunNames <- unique(unlist(strsplit(as.character(settings[,'tFun']),'[[:punct:] ]+')))
   vNames <- unique(c(vNames, tFunNames[tFunNames %in% colnames(data)]))
   return(data[row.names(na.omit(data[,vNames])),])
}

`SDDataInitFun` <- 
function(x){
   if(class(x) == 'SDisc')
      return(SDDataInitFun(SDData(x)))
   else if(class(x) == 'SDData')
      return(attr(x, 'initFun'))
   else
      cat2('You must provide either an \'SDisc\' or \'SDData\' object to \'SDDataInitFun\'')
}

`SDDataSetup` <- 
function(x, prefix=NULL, xlim=c(-3,3), ylim=c(0,50), mai=c(0.6,0.3,0.05,0.05)){
   if(class(x) == 'SDData'){
      x <- SDPrefix(x, prefix) # do basedir, figdir and tabdir
      x <- SDDataXlim(x, xlim)
      x <- SDDataYlim(x, ylim)
      x <- SDDataMai(x, mai)
      x <- SDDataRInfo(x) 
      return(x)
   }
   else
      cat2('You must provide an \'SDData\' to \'SDDataSetup\'.')
}

`bicTable.default` <- 
function(x, ...){
   if(class(x) == 'SDisc')
      return(attr(x, 'bicTable'))
   else
      cat2('x must of class \'SDisc\'.')
}

`bicStats` <- 
function(x, fun='mean', vAgg="BIC", vX="G", vY="modelName", fmt='%.2f'){
   rList <- list()
   # DEFINE FUNCTIONS 
   getRankArray <- function(m){return(apply(-m,c(2,3),rank))}
   fReshape <- function(m, pX=vX, pY=vY, pZ=vAgg){
      m <- reshape(m[,c(pX,pY,pZ)],idvar=pX, timevar=pY, direction='wide')
      dimnames(m) <- list(m[,1], gsub(paste2(pZ,'.'),'',colnames(m)))
      return(m[,-1])
   }
   fCI <- function(x,pFmt=fmt){
      x <- quantile(x,probs=c(0.025,0.975), na.rm=TRUE)
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
         tmp <- abind(tmp, fReshape(x[x$rseed == vRSeed[s],]), along=3)
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

`plotInit` <- 
function(x=NULL, figIdx=1, plotNbr=1, type="H", latex=FALSE, legend='', lab=NULL){
   plotI <- plotJ <- 1
   while(plotI * plotJ < plotNbr){
      if(plotI == plotJ)
         plotI <- plotI+1
      else
         plotJ <- plotJ+1
   }
   maiParam <- c(0.6,0.3,0.05,0.4) 
   fDir <- '.'
   fLeg <- paste2(SDPrefix(x), ', ')
   if(!is.null(x)){
      maiParam <- attr(x,'mai')
      fDir <- SDFigDir(x)
   }
   width <- '8cm'
   fLabel <- paste2('oddGroup-', sprintf('%03d', figIdx), '-', type)
   if(type=='H'){
      maiParam <- maiParam+c(0,0,0,0)
      fLeg <- paste2(fLeg, '\\textbf{histograms} of the variables of the factor \\textbf{', legend, '}.')
   }
   else if(type=='LM'){
      maiParam <- maiParam+c(0.5,0.2,0.5,0)
      fLeg <- paste2(fLeg, legend)
      width <- '8cm'
   }
   else if(type=='BB'){
      maiParam <- maiParam+c(0.5,0.2,0.5,0)
      fLeg <- paste2(fLeg, '\\textbf{boxplots} of the variables of the factor \\textbf{', legend, '}.')
   }
   else{
      maiParam <- maiParam+c(0.1,0.5,0.2,0)
      fLabel <- type
      fLeg <- paste2(fLeg, ' visual representation of \\textbf{model} ', gsub('-',',', gsub('MM-','',type)), '.')
      width <- '16cm'
   }
   fPDF <- paste2(fDir, '/', fLabel,'.pdf')
   if(latex){
      if(is.null(lab) && (type=='H' || type=='B'))
         lab <- paste2('oddGroup-',sprintf('%03d',figIdx),'-',gsub('-',',',type))
      cat2('\\begin{figure}\\begin{center}')
      cat2('\\href{',fPDF,'}{\\includegraphics[width=',width,']{',fPDF,'}}')
      cat2('\\caption{\\label{fig:',lab,'}',fLeg,'}')
      cat2('\\end{center}\\end{figure}')
   }
   else 
      cat2(fPDF)
   pdf(fPDF)
   par(mfrow = c(plotI,plotJ), mai=maiParam)
   figIdx <- figIdx+1
   return(figIdx)
}

`plotClose` <-
function(){ graphics.off() }

`plot.SDData` <- 
function(x, q=NULL, est=1, zlim=c(-2,2), latex=FALSE, ...) {
   xSettings <- as.matrix(na.omit(SDDataSettings(x)[,c("vParGroup","vParY")]))
   splitGroup <- split(row.names(xSettings),xSettings[,"vParGroup"])
   figIdxBB <- figIdxH <- 1
   if(is.null(q)){
      for(i in 1:length(splitGroup)){
         # PRODUCE BOXPLOT BY SPLIT GROUP
         figIdxBB <- plotInit(x, figIdxBB, plotNbr=1, type='BB', latex=latex, legend=names(splitGroup)[i])
         boxplot(as.data.frame(x[,splitGroup[[i]]]),'las'=2,cex=0.7, main=paste2(" Boxplot ",names(splitGroup)[i]))
         plotClose()
         figIdxH <- plotInit(x, figIdxH, plotNbr=length(splitGroup[[i]]), type='H', latex=latex,
            legend=names(splitGroup)[i])
         for(v in splitGroup[[i]]){
            lab <- paste2(v,'\n', length(table(x[,v])), " values\n")
            hist(x[,v], main=NULL, xlab=NULL, ylab="", cex.lab=0.75, las=2,'new'=TRUE, col='grey',border=0)
            mtext(lab, cex=0.75, side=1)
         }
         plotClose()
      }
   }
   else if(!is.null(q)){
      res <- list()
      cNames <- colnames(print(x))
      m <- matrix(NA, length(cNames), length(q),dimnames=list(cNames, q))
      for(i in q){
         for(j in cNames){
            val <- summary(SDTData(x)[[j]], q=i, digits=9)
            if(length(val) > 1)
               val <- strsplit(val[grep(i, colnames(val))], "[ )(;]+")[[1]][est]
            m[j,i] <- as.numeric(val)
            if(est==3)
               m[j,i] <- -log(m[j,i])
         }
      }
      m <- t(m)
      cex <- 0.7
      figIdxBB <- plotInit(x, 1, plotNbr=1, type='LM', latex=latex, legend=paste2(q))
      colGradient <- rev(brewer.pal(9,"RdBu"))
      breaks <- seq(zlim[1],zlim[2],length.out=10)
      if(est==3){
         colGradient <- brewer.pal(3,"Blues")
         breaks <- c(0, -log(0.05),-log(0.01), 1000)
      }
      image(1L:nrow(m), 1L:ncol(m), m, xlim=0.5+c(0,nrow(m)), ylim=0.5+c(0,ncol(m)), axes= FALSE, xlab= "",
         ylab="", col=colGradient, main=paste(q, collapse=', '), breaks=breaks, add=FALSE, ...) 
      axis(2, at=(1L:ncol(m)), labels=colnames(m), las=2, line=-0.5, tick=0, cex.axis=cex)
      axis(1, 1L:nrow(m), labels=row.names(m), las=1, line=-0.5, tick=0, cex.axis=cex)
      plotClose();
   }
}

`plot.SDisc` <- 
function(x, q=NULL, type=c('plotParcoord', 'plotLegend', 'plotPC1', 'plotPC2', 'plotDendroCluster', 'plotDendroVar'),
   latex=FALSE, title=NULL, xlim=c(-3, 3), zlim=c(-2,2), xy=c(-2.2, 0), pattern="mean", cex=0.7,
   colGrad=rev(brewer.pal(9,"RdBu")), rangeFV=NULL, lab=NULL, ...) {
   if(is.numeric(q))
      q <- which(names(x) %in% bestModel(x, n=q))
   else if(is.null(q))
      q <- which(names(x) %in% bestModel(x, n=1))
   else
      q <- which(names(x) %in% q)
   m <- x[[q]]
   if(!is.na(attr(x[[q]],"model")$loglik)){
      plotNbr <- length(type)
      if(('plotLegend' %in% type))
         plotNbr <- plotNbr-1
      if(('plotParcoord' %in% type))
         plotNbr <- plotNbr+length(na.omit(unique(SDDataSettings(x)[,'vParGroup'])))-1
      figIdx <- plotInit(SDData(x), figIdx=1, plotNbr=plotNbr, type=paste2('MM-',names(x)[q]), lab=lab,
         latex=latex)
      for(plotFun in type)
         do.call(plotFun, args=list(x=m, data=SDData(x), title=title, xlim=xlim, zlim=zlim, xy=xy, cex=cex,
            pattern=pattern, colGrad=colGrad, rangeFV=rangeFV))
      plotClose()
   }
}

`plotSeries` <-
function(x, xlim=c(-2,2), cex=0.7, pattern='mean', colGrad=rev(brewer.pal(9,"RdBu")), rangeFV=NULL, ...){
   if(is.null(rangeFV))
      rangeFV <- 1:ncol(x)
   rangeFV <- colnames(x)[rangeFV]
   p <- SDCModelPattern(x, fun=pattern)[,rangeFV]
   plot(ts(t(p)), plot.type='single', col=SDCModelColors(x), lwd=2, axes=FALSE, 
      ylim=xlim, xlab='',ylab='', main=paste(SDCModelSettings(x),sep=',',collapse=','))
   fv_idx <- 1:ncol(p)
   if(ncol(p) > 20)
      fv_idx <- as.integer(quantile(fv_idx))
   axis(1, at=(1L:ncol(p))[fv_idx], labels = colnames(p)[fv_idx], las = 2, line = -0.5, tick = TRUE, cex.axis = cex)
   axis(2, at=as.integer(seq(xlim[1],xlim[2],length.out=4)),las=2)
   d <- table(SDCModel(x)$label)[row.names(p)]
   d <- apply(cbind(names(d),' (',d,')'),1,paste2)
   legend("bottomleft",legend=d,col=SDCModelColors(x),lwd=2, bty='n')
}

`plotImage` <- 
function(x, data, title=NULL, zlim=c(-2,2), pattern='mean', cex=0.7, colGrad=rev(brewer.pal(9,"RdBu")),
orderPattern=TRUE, ...){
   muP <- SDCModelPattern(x, data=data, fun=pattern)
   if(orderPattern){
      cOrder <- hclust(dist(muP))$order
      row.names(muP) <- as.character(1:nrow(muP))
      muP <- muP[cOrder,]
   }
   image(1L:nrow(muP), 1L:ncol(muP), muP, xlim=0.5+c(0,nrow(muP)), ylim=0.5+c(0,ncol(muP)), axes= FALSE, xlab=
      "", ylab="", col=colGrad, main=title,add=FALSE, zlim=zlim)
   fId <- 1:ncol(muP)
   if(ncol(muP) > 20)
      fId <- as.integer(quantile(fId))
   axis(2, at=(1L:ncol(muP))[fId], labels=colnames(muP)[fId], las=2, line=-0.5, tick=0, cex.axis=cex)
   axis(1, 1L:nrow(muP), labels=row.names(muP), las=1, line=-0.5, tick=0, cex.axis=cex)
}

`plotLegend` <- 
function(x, xy=c(-2.2, 0), cex=0.7, colGrad=rev(brewer.pal(9,"RdBu")), ...){
   p  <- SDCModelPattern(x, fun="clCount")
   for(g in 1:nrow(p)){
      g_name <- SDCModelColors(x)[g]
      y1 <- y0 <- xy[2] + 5 * g / nrow(p)
      x0 <- xy[1]
      x1 <- 1.2 * xy[1]
      arrows(x0, y0, x1, y1, col=SDCModelColors(x)[g], length=0, lwd=3)
      text(x0+0.15, y1, labels=paste2(names(p)[g], " (",p[g],")"), pos=4)
   }
}

`plotDendroCluster` <- 
function(x, title=NULL, cex=0.7, pattern='mean', ...){
   mu <- SDCModelPattern(x, fun=pattern)
   plot(hclust(dist(mu)), axes = FALSE, yaxs = "i", main=title, ylab=NULL, cex=cex)
}

`plotDendroVar` <- 
function(x, title=NULL, cex=0.7, pattern='mean', ...){
   mu <- SDCModelPattern(x, fun=pattern)
   plot(hclust(dist(t(mu))), axes=FALSE, yaxs="i", main=title, ylab=NULL, cex=cex)
}


`plotParcoord` <- 
function(x, data, title=NULL, xlim=c(-3, 3), pattern='mean', cex=0.7, colGrad=rev(brewer.pal(9,"RdBu")), ...){
   set <- as.matrix(na.omit(SDDataSettings(data)))
   for(v in unique(set[,'vParGroup'])){
      T_s <- as.matrix(set[set[,"vParGroup"] == v,])
      min2 <- function(x){if(min(x)<0){return(min(x))}else{return(0)}}
      max2 <- function(x){if(max(x)>10){return(max(x))}else{return(10)}}
      ylim<- c(min2(as.numeric(T_s[,"vParY"])),max2(as.numeric(T_s[,"vParY"])))
      plot(x=0 ,new=TRUE ,ann=FALSE ,pch=18,col="white",axes=FALSE ,xlim=xlim , ylim=ylim)
      title(main = paste2("(",paste(as.matrix(SDCModelSettings(x)),collapse=","),") ", v), cex=cex)
      for(pName in pattern){
         p <- SDCModelPattern(x, data=data, fun=pName)[,row.names(T_s)]
         lwd <- 3
         lty <- "solid"
         if(pName != "mean") {
            lwd <- 1
            lty <- "dashed"
         }
         for(g in 1:nrow(p)){
            D_i_s <- cbind(X=as.numeric(p[g,]),Y=as.numeric(T_s[,"vParY"]))
            D_i_s <- D_i_s[sort.list(D_i_s[,"Y"]),]
            gap <- 1.05*abs(median(D_i_s[1:(nrow(D_i_s)-1),"Y"] - D_i_s[2:nrow(D_i_s),"Y"]))
            for(l in 1:(nrow(D_i_s)-1))
               if(!(D_i_s[l+1,"Y"] - D_i_s[l,"Y"] > gap))
                  arrows(D_i_s[l,1], D_i_s[l,2], D_i_s[l+1,1], D_i_s[l+1,2], col=SDCModelColors(x)[g], length=0, lwd=lwd, lty=lty)
               # ELSE, AS (is_white_gap == TRUE) THEN DO NOT DRAW ANY ARROW...
         }
      }
      idFV <- 1:nrow(T_s)
      if(nrow(T_s) > 20)
         idFV <- as.integer(quantile(1:nrow(T_s)))
      axis(2, at=as.numeric(T_s[idFV,"vParY"]), labels=colnames(p)[idFV], las=2, tick=FALSE)
      axis(1, at=seq(from=xlim[1], to=xlim[2], by=(xlim[2]-xlim[1])/4), cex=cex)
   }
} 

`SDisc.default` <-
function(x, cfun='modelBasedEM', cFunSettings=list(modelName=c("EII", "VII"), G=3:5, rseed=6013:6015), nTopModels=5, nnodes=1, ...){
   # TRANSFORM THE DATA
   cat2('Prepare the data')
   data <- SDData(x, ...)
   # CALCULATE THE ENSEMBLE OF MODELS
   mSet <- expand.grid(cFunSettings)
   row.names(mSet) <- apply(as.matrix(mSet[, names(cFunSettings)]), 1, paste, collapse=",")
   # TEST IF THE CALCULATION SHOULD BE RE-DONE 
   fIMG <- paste2(SDBaseDir(data),"/IMAGE.RData")
   x <- list()
   if(file.exists(fIMG)){
      cat2('Load and test for consistency: ', fIMG)
      load(fIMG)
      eqData <- (digest(SDDataOrig(data)) == digest(SDDataOrig(x)))
      eqSettings <- (digest(SDDataSettings(data)) == digest(SDDataSettings(x)))
      eqCfun <- (digest(cfun) == digest(SDCFun(x))) 
      eqCfunSettings <- (digest(mSet) == digest(SDSettings(x))) 
      if(isModeled(x) && eqData && eqSettings && eqCfun && eqCfunSettings)
         return(x)
      else
         x <- list()
   }
   # ELSE, CONTINUE
   for(i in row.names(mSet))
      x[[i]] <- SDCModel(i, settings=mSet[i,])
   df <- print(data)
   fun <- function(y, yData=df){ return(do.call(cfun, args=list(y, data=yData))) }
   cat2('Modeling for clusters ')
   startTime <- Sys.time()
   x <- SDLapply(x, fun, nnodes=nnodes)
   stopTime <- Sys.time()
   x <- structure(x, SDData=data, prefix=SDPrefix(data), cFunSettings=mSet, nTopModels=nTopModels,
      bicTable=NULL, startTime=startTime, stopTime=stopTime, rinfo=sessionInfo(), cfun=cfun, class = "SDisc")
   cat2('Collect BICs (likelihood) of the models ')
   bicTable <- cbind(mSet, BIC=NA, relativeBic=NA)
   for(i in 1:length(x))
      bicTable[names(x)[i],"BIC"] <- attr(x[[i]],"model")[["BIC"]]
   bicTable[,"relativeBic"] <- 100*(bicTable[,"BIC"]/max(bicTable[,"BIC"],na.rm=TRUE)-1)
   attr(x,'bicTable') <- structure(table(bicTable[,1:2]), prefix=SDPrefix(x), data=bicTable, class='bicTable')
   save.SDisc(x)
   write.SDisc(x)
   return(x)
}

`paste2` <-
function(...){ return(paste(...,collapse="",sep="")) }

`write.SDisc` <-
function(x, q = NULL){
   cat2("Save best models as CSV files ")
   for(qName in bestModel(x,n=q)){
      fCSV        <- paste2(SDTabDir(x),'/MM-',qName,".csv") 
      cat2(fCSV)
      write.csv2(print(x[[qName]], data=SDData(x)), file=fCSV)
   }
}

`lowQuant` <- 
function(x, ...){ return(quantile(x, probs=0.025, ...)) }

`upQuant` <- 
function(x, ...){ return(quantile(x, probs=0.975, ...)) }

`SDDataSettings` <- 
function(x, asCSV=FALSE, inCAnalysis=NULL, latex=FALSE){ 
   if((class(x) == 'SDData' || class(x) == 'SDisc') & !latex)
      return(SDDataAttr(x, attrName='settings')) 
   if((class(x) == 'SDData' || class(x) == 'SDisc') & latex)
      texTable(SDDataAttr(x, attrName='settings'), type='latex', cap='SDDataSettings')
   else{
      cNames <- list('oddGroup','inCAnalysis','tFun','vParGroup','vParY','vHeatmapY')
      x <- SDDataDimnames(x)
      v <- c("NA",TRUE,"mean sd","varGroup1","NA","NA")
      x <- matrix(v,ncol(x),length(v),byrow=TRUE, dimnames=list(colnames(x),cNames))
      x[,'vHeatmapY'] <- x[,'vParY']<- 1:nrow(x)
      x[,'oddGroup'] <- row.names(x)
      if(length(inCAnalysis)>0){
         x[,'inCAnalysis'] <- 'FALSE' 
         x[inCAnalysis,'inCAnalysis'] <- 'TRUE'
         x[x[,'inCAnalysis'] == FALSE,'tFun'] <- NA
      }
      if(!latex && asCSV == TRUE)
         write.csv(x,file='settings.csv')
      else if(!latex && asCSV == FALSE)
         return(x)
      else if(!latex && is.character(asCSV))
         write.csv(x,file=asCSV)
      else if(latex)
         texTable(x, type='latex', cap='SDDataSettings')
   }
}

`SDSettings` <- 
function(x){ 
   if(class(x) == 'SDisc')
      return(attr(x, 'cFunSettings') ) 
   else
      cat2('A \'SDisc\' object must be provided.')
}

`SDCFun` <-
function(x){
   if(class(x) == 'SDisc')
      return(attr(x,'cfun'))
   else
      cat2('A \'SDisc\' object must be provided.')
}

`SDRInfo` <-
function(x){
   if(class(x) == 'SDisc' || class(x) == 'SDData')
      return(attr(x,'rinfo'))
   else
      cat2('A \'SDisc\' object must be provided.')
}

`SDCModel` <- 
function(x, n=1, settings=NULL){
   if(class(x) == 'SDCModel')
      return(attr(x, 'model'))
   else if(class(x) == 'SDisc')
      return(attr(x[[bestModel(x,n=n)]],'SDCModel'))
   else
      return(structure(x, cFunSettings=settings, model=NULL, pattern=list(), GColors=NULL, stats=list(),
         isModeled=FALSE, class="SDCModel"))
}

`SDCModelAttr` <- 
function(x, attrName){
   if(class(x) == 'SDCModel')
      return(attr(x, attrName))
   else
      cat2('A \'SDCModel\' object must be provided.')
}

`SDCModelPattern` <- 
function(x, data=NULL, fun='mean'){ 
   if(is.null(data)){
      p <- SDCModelAttr(x, attrName='pattern')
      if(fun %in% names(p))
         return(p[[fun]]) 
      else
         cat2('Data must be provided to compute this pattern ')
   }
   else{
      df <- cbind(print(data, allNumVars=TRUE), class=print(x,data=data)[,'class'])
      return(SDPattern(df, fun))
   }
}

`SDCModelColors` <- 
function(x){ return(SDCModelAttr(x, attrName='GColors')) }

`SDCModelSettings` <- 
function(x){ return(SDCModelAttr(x, attrName='cFunSettings')) }

`SDPrefix` <- 
function(x, value=NULL){ 
   if(is.null(value))
      return(SDDataAttr(x, attrName='prefix')) 
   else{
      x <- SDDataAttr(x, attrName='prefix', value=value)
      x <- SDBaseDir(x, SDPrefix(x))
      x <- SDFigDir(x, paste2(SDPrefix(x),'/figures'))
      x <- SDTabDir(x, paste2(SDPrefix(x),'/tables'))
   }
   return(x)
}

`SDBaseDir` <- 
function(x, value=NULL){ 
   if(!is.null(value))
      dir.create(value, showWarnings=FALSE)
   return(SDDataAttr(x, attrName='baseDir', value=value)) }

`SDFigDir` <- 
function(x, value=NULL){ 
   if(!is.null(value))
      dir.create(value, showWarnings=FALSE)
   return(SDDataAttr(x, attrName='figDir', value=value)) 
}

`SDTabDir` <- 
function(x, value=NULL){ 
   if(!is.null(value))
      dir.create(value, showWarnings=FALSE)
   return(SDDataAttr(x, attrName='tabDir', value=value)) 
}

`SDDataSubset` <- 
function(x, value=NULL){ return(SDDataAttr(x, attrName='subset', value=value)) }

`SDDataOrig` <- 
function(x, value=NULL){ return(SDDataAttr(x, attrName='dataOrig', value=value)) }

`SDDataXY` <- 
function(x, value=NULL){ return(SDDataAttr(x, attrName='xy', value=value)) }

`SDDataXlim` <- 
function(x, value=NULL){ return(SDDataAttr(x, attrName='xlim', value=value)) }

`SDDataYlim` <- 
function(x, value=NULL){ return(SDDataAttr(x, attrName='ylim', value=value)) }

`SDDataMai` <- 
function(x, value=NULL){ return(SDDataAttr(x, attrName='mai', value=value)) }

`SDDataRInfo` <- 
function(x, value=NULL){ return(SDDataAttr(x, attrName='rinfo', value=sessionInfo())) }

`SDTDataElmt` <- 
function(x, data=NULL, var=NULL, fModel=NULL){
   if(class(x) == 'SDTData')
      return(attr(x,'elmt'))
   else if(class(x) == 'SDData')
      cat2('There are ', length(SDTData(x)), ' elements in the \'SDTData\'. Select one.')
   else if(class(x) == 'SDisc')
      return(SDTDataElmt(SDData(x)))
   else if(!is.null(data) && !is.null(var)){
      f <- strsplit(x,"[()']+")[[1]]
      if(is.null(fModel)){
         if(!is.na(f[2]))
            fArgs <- list(formula=f[2], data=as.data.frame(data[, all.vars(as.formula(f[2]))]))
         if(f[1] %in% c('mean', 'min', 'median','max', 'sd'))
            fArgs <- list(data[,var], na.rm=TRUE)
         fModel <- do.call(f[1], args=fArgs)
      }
      if(!is.null(fModel)){
         fArgs <- list(data[,var], fModel)
         if(f[1] %in% c('mean', 'min', 'median'))
            fRes <- do.call('-', fArgs)
         if(f[1] %in% c('max', 'sd'))
            fRes <- do.call(function(a, b, lambda=quantile(abs(a), probs=0.01)/1000){
               return(a/(lambda+b))}, fArgs)
         if(!is.na(f[2])){
            fArgs <- list(formula=f[2], data=as.data.frame(data[,all.vars(as.formula(f[2]))]))
            fRes <- predict(fModel, newdata=fArgs[['data']])-data[,var]
         }
         names(fRes) <- row.names(data)
      }
      return(structure(x, vName=var, fName=f[1], fArgs=fArgs, fModel=fModel, fRes=fRes, class="SDTDataElmt"))
   }
   else
      return(NULL)
}

`SDTData` <- 
function(x, tFun=NULL, data=NULL, TData=NULL){ 
   if(class(x) == 'SDisc')
      return(SDTData(SDData(x))) 
   else if(class(x) == 'SDData')
      return(SDDataAttr(x, attrName='SDTData') ) 
   else{
      tFunElmts <- strsplit(as.character(tFun), "[, ]+")[[1]]
      elmts <- list()
      newdata <- data
      for(f in tFunElmts){
         fModel <- attr(SDTDataElmt(TData)[[f]], 'fModel')
         elmts[[f]] <- res <- SDTDataElmt(f, data=newdata, var=x, fModel=fModel)
         newdata[names(attr(res,'fRes')), x] <- attr(res,'fRes')
      }
      return(list(data=newdata, tdata=structure(x, vName=x, tFun=tFun, elmts=elmts, class='SDTData')))
   }
}

`SDData.default` <-
function(x, prefix, dataOrig=NULL, TData=NULL, settings=NULL, initFun=list(SDDataCC), subset=NULL, ...){
   if(class(x) == 'SDisc')
      return(attr(x,'SDData'))
   else if(class(x) == 'SDData'){
      # x <- SDPrefix(x, prefix) # do basedir, figdir and tabdir
      if(!is.null(subset)) # extract a subset of the data
         return(structure(print(x)[which(subset %in% row.names(print(x))),], dataOrig=SDDataOrig(x),
            initFun=initFun, TData=SDTData(x), settings=SDDataSettings(x), xlim=SDDataXlim(x),
            ylim=SDDataYlim(x), mai=SDDataMai(x), prefix=SDPrefix(x), baseDir=SDBaseDir(x), figDir=SDFigDir(x),
            tabDir=SDTabDir(x), rinfo=sessionInfo(), subset=subset, class='SDData'))
      else # or simply returns the original object 
         return(x)
   }
   else{
      x <- SDDataDimnames(x)
      if(is.null(settings))
         settings <- SDDataSettings(x)
      else if(length(settings) == 1)
         if(file.exists(settings))
            settings <- read.csv(settings, row.names=1, sep=';')
      if(is.null(dataOrig))
         dataOrig <- x
      for(fun in initFun)
         x <- fun(x, settings)
      x <- x[,row.names(settings)]
      if(!is.null(subset))
         x <- x[which(subset %in% row.names(x)),]
      tRes <- list()
      newdata <- x
      for(i in row.names(settings[!is.na(settings[,'tFun']),])){
         res <- SDTData(i, tFun=settings[i,'tFun'], data=newdata)
         tRes[[i]] <- res[["tdata"]]
         newdata <- res[["data"]]
      }
      res <- structure(data.matrix(newdata), dataOrig=dataOrig, initFun=initFun, SDTData=tRes, settings=settings,
         prefix=prefix, class='SDData')
      res <- SDDataSetup(res, prefix=prefix)
      return(res)
   }
}


`SDDataAttr` <- 
function(x, attrName, value=NULL){
   if(class(x) == 'SDisc' && !is.null(value)){
      attr(x, 'SDData') <- SDDataAttr(SDData(x), attrName=attrName, value=value)
      return(x)
   }
   else if(class(x) == 'SDData' && !is.null(value)){
      attr(x, attrName) <- value
      return(x)
   }
   else if(class(x) == 'SDData' || class(x) == 'bicTable')
      return(attr(x, attrName))
   else if(class(x) == 'SDisc')
      return(SDDataAttr(SDData(x), attrName=attrName, value=value))
   else
      cat2('An \'SDisc\', \'SDData\' or \'bicTable\' object must be provided.')
}

`isModeled` <- 
function(x, set=NULL){
   if(class(x) == 'SDisc'){
      res <- TRUE
      for(i in 1:length(x))
         res <- res && isModeled(x[[i]])
      return(res)
   }
   else if(class(x) == 'SDCModel' && is.null(set))
      return(attr(x,'isModeled'))
   else if(class(x) == 'SDCModel' && !is.null(set)){
      attr(x,'isModeled') <- TRUE
      return(x)
   }
   else
      return(NULL)
}

`oddRatios` <-
function(x, class, fun=median, lambda=0.01){
   df <- cbind(x, class=map(class, warn=FALSE))
   set <- as.data.frame(na.omit(SDDataSettings(x)[,c("vParGroup","oddGroup")]))
   vGr <- sort(unique(as.character(set[,"oddGroup"])))
   m <- matrix(0,ncol(class), length(vGr),dimnames=list(1:ncol(class), vGr))
   nrow2 <- function(y){y <- nrow(y) ; if(is.null(y)){return(0)}else{return(y)}}
   sum2 <- function(y){if(!is.null(dim(y))){return(apply(y,1,sum))}else{return(y)}}
   for(gr in vGr){
      mat <- cbind(SScore=sum2(df[,row.names(set[set[,"oddGroup"]==gr,])]),class=map(class, warn=FALSE))
      for(g in 1:ncol(class)){
         mu <- fun(mat[,"SScore"], na.rm=TRUE)
         inG <- mat[,"class"]==g
         supMu <- mat[,"SScore"] >= mu
         m11 <- nrow2(mat[ inG &  supMu,])
         m12 <- nrow2(mat[ inG & !supMu,])
         m21 <- nrow2(mat[!inG &  supMu,])
         m22 <- nrow2(mat[!inG & !supMu,])
         p <- matrix(c(m11,m12,m21,m22),2,2,byrow=TRUE)
         m[g, match(gr,vGr)] <- log(abs(lambda+p[1,1]*p[2,2]) / abs(lambda+p[1,2]*p[2,1]))
      }
   }
   legend <- paste2(SDPrefix(x), ', \\textbf{oddratios} for the main factors ')
   res <- structure(m, legend=legend, class='SDStats')
   return(res)
}

`oddRatiosB` <-
function(x, class, fun=median, lambda=0.01){
   df <- cbind(x, class=map(class, warn=FALSE))
   set <- as.data.frame(na.omit(SDDataSettings(x)[,c("vParGroup","oddGroup")]))
   vGr <- sort(unique(as.character(set[,"oddGroup"])))
   m <- matrix(0,ncol(class), length(vGr),dimnames=list(1:ncol(class), vGr))
   sum2 <- function(y){if(!is.null(dim(y))){return(apply(y,1,sum))}else{return(y)}}
   for(gr in vGr){
      mat <- cbind(SScore=sum2(df[,row.names(set[set[,"oddGroup"]==gr,])]),class=map(class, warn=FALSE))
      for(g in 1:ncol(class)){
         mu <- fun(mat[,"SScore"], na.rm=TRUE)
         inG <- mat[,"class"]==g
         supMu <- mat[,"SScore"] >= mu
         m11 <- sum(mat[ inG &  supMu,] - mu*mat[ inG &  supMu,])
         m12 <- sum(mat[ inG & !supMu,] - mu*mat[ inG & !supMu,])
         m21 <- sum(mat[!inG &  supMu,] - mu*mat[!inG &  supMu,])
         m22 <- sum(mat[!inG & !supMu,] - mu*mat[!inG & !supMu,])
         p <- matrix(c(m11,m12,m21,m22), 2, 2, byrow=TRUE)
         m[g,match(gr,vGr)] <- log(abs(lambda+p[1,1]*p[2,2])/abs(lambda+p[1,2]*p[2,1]))
      }
   }
   legend <- paste2(SDPrefix(x), ', (Bayesian) \\textbf{oddratios} for the main factors ')
   res <- structure(m, legend=legend, class='SDStats')
   return(res)
}

`chi2test` <- 
function(x, class, target='Class'){
   m <- cbind(x, modelClass=map(class))
   test <- list()
   for(tName in target){
      m2 <- cbind(modelClass=m[, 'modelClass'], target=SDDataOrig(x)[row.names(m),tName])
      test[[tName]] <- chisq.test(xtabs( ~.,m2[,c('modelClass',"target")]), simulate.p.value=TRUE)
   }
   if(length(target) == 1){
      test <- test[[1]]
      legend <- paste2('For \\textbf{', target, '}: $p_{\\chi^2}=',
         sprintf('%0.3f',test[["p.value"]]), '$ $(\\chi^2=',
         sprintf('%.1f',test[["statistic"]]),')$')
      res <- structure(test[["residuals"]], test=test, legend=legend, class='SDStats')
   }
   else if(length(target) > 1){
      tab <- matrix(NA, length(target), 2, dimnames=list(target, list('p-value','X2-sum')))
      for(tName in target)
         tab[tName,] <- c(test[[tName]][["p.value"]], test[[tName]][["statistic"]])
      signTargets <- row.names(tab)[which(tab[,'p-value'] < 0.05)]
      if(length(signTargets) >= 1)
         legend <- paste2(SDPrefix(x), ', variable(s) \\textbf{', paste(signTargets, collapse=' ', sep=''), '}
            show a significant association ($\\chi^2$) with the identified subtypes ')
      else
         legend <- paste2('No variable does show a significant association with
            the identified subtypes.')
      res <- structure(t(tab), test=test, legend=legend, class='SDStats')
   }
   else
      cat2('You must provide at least one \'target\' for chi2-testing')
   return(res)
}

`jointDistrib` <- 
function(x, class, target='Class'){
   df <- cbind(x, modelClass=map(class))
   df <- cbind(df, target=as.character(SDDataOrig(x)[row.names(df),target]))
   m <- xtabs( ~.,df[,c("modelClass","target")])
   test <- chisq.test(m, simulate.p.value=TRUE)
   p <- sprintf('%.3f', test[["p.value"]])
   chi2 <- sprintf('%3f', test[["statistic"]])
   legend <- paste2('Joint distribution of the \\textbf{', target,'} and the different subtypes
      ($p_{\\chi^2}=',p,', \\chi^2=', chi2,'$)')
   res <- structure(m, test=test, legend=legend, class='SDStats')
   return(res)
}

`today` <-
function(){
   tmp_date <- format(Sys.time(), "%Y-%m-%d")
   return(tmp_date)
}

`SDLapply` <- 
function(x, fun, nnodes=1, ...){
   if(nnodes > 1){
      cat2('nnodes:',nnodes)
      cl <- startSDCluster(nnodes)
      set.seed(6013)
      sList <- lapply(split(sample(names(x)),rep(1:nnodes, length.out=length(x))), function(i) x[i])
      res <- docall(c, clusterApply(cl, sList, lapply,fun, ...))
      #res <- parLapply(cl, x, fun, ...) 
      stopCluster(cl)
   }
   else 
      res <- lapply(x, fun, ...) 
   return(res)
}

`startSDCluster` <- 
function(nNodes){
   cl <- makeCluster(rep('localhost',nNodes), type = "SOCK")
   clusterCall(cl, library, 'SDisc',character.only=TRUE)
   return(cl)
}

`save.SDisc` <- 
function(x){
   fSave <- paste2(SDBaseDir(x),"/IMAGE.RData")
   cat2("Save modeling into ",fSave)
   save(list='x', file=fSave)
}

`SDDataInCAnalysis` <- 
function(x){
   if(class(x) == 'SDData')
      x <- SDDataSettings(x)
   varSel <- which(x[,"inCAnalysis"] == "TRUE" | x[,"inCAnalysis"] == " TRUE")
   return(row.names(x)[varSel])
}

`print.bicTable` <-
function(x, n=NULL, modelName=NULL, G=NULL, latex=FALSE, lab='bic5', ...){
   m <- attr(x,'data')
   m <- m[order(m[,"relativeBic"], decreasing=FALSE, na.last=NA),]
   if(is.null(n))
      m <- m[1:nrow(m),]
   else if(n < 1)
      m <- m[which(m[,'relativeBic'] < n*100),]
   else
      m <- m[1:n,]
   if(!is.null(G))
      m <- m[m$G == G, ]
   if(!is.null(modelName))
      m <- m[m$modelName == modelName, ]
   if(latex)
      texTable(m, type='latex', cap=paste2('Top ranking \\textbf{',SDPrefix(x),'} models.'), lab=lab)
   else
      return(m)
}

`print.SDStats` <-
function(x, ...){
   return(t(x[1:nrow(x),]))
}

`print.SDData` <-
function(x, rseed=NULL, range=1:3, allNumVars=FALSE, latex=FALSE, ...){
   if(allNumVars)
      m <- as.matrix(x[1:nrow(x),!is.na(SDDataSettings(x)[,'tFun'])])
   else
      m <- as.matrix(x[1:nrow(x),SDDataInCAnalysis(x)])
   if(is.null(rseed)){
      if(latex)
         texTable(m, cap=paste2('\\textbf{',SDPrefix(x),'}, extract of the \\textbf{transformed} data matrix.'),
            lab=paste2('SDData',SDPrefix(x)))
      else
         return(m)
   }
   else if(!is.null(rseed)){
      set.seed(rseed) ; rNames <- row.names(m)[sample(1:nrow(m))[range]]
      set.seed(rseed) ; cNames <- colnames(m)[sample(1:ncol(m))[range]]
      m <- m[rNames, cNames] 
      mOrig <- SDDataOrig(x)[rNames, cNames]
      if(latex){
         texTable(mOrig, cap=paste2('\\textbf{',SDPrefix(x),'}, extract of the \\textbf{original} data matrix.'),
            lab=paste2('SDData',SDPrefix(x)))
         texTable(m, cap=paste2('\\textbf{',SDPrefix(x),'}, extract of the \\textbf{transformed} data matrix.'),
            lab=paste2('SDData',SDPrefix(x)))
      }
      else
         return(list(data=m, dataOrig=mOrig))
   }
}

`print.SDCModel` <-
function(x, data, ...){
   model <- em(data=print(data), modelName=SDCModel(x)$modelName, parameters=SDCModel(x)$parameters, warn=FALSE)
   z <- model[["z"]]
   z <- cbind(z, class=map(z))
   colnames(z)[1:(ncol(z)-1)] <- 1:(ncol(z)-1)
   return(z)
}

`print.SDisc` <- 
function(x, y=NULL, m1=1, m2=2, latex=FALSE, lab='jointdistrib', ...) {
   doCaption <- function(o, cap){
      cap <- paste2(SDPrefix(x), ", the \\textbf{comparison of model} ", gsub('x','and', gsub('[()]', ' ', cap)), "
         exhibits a random index ", sprintf('%.1f',100*o$rand)," (a $\\kappa=", sprintf('%.1f', 100*o$kappa),"$ ,
         and a relative degree of association $V=", sprintf('%.1f', 100*o$V), "\\%$ with $p_{\\chi^2}=",
         sprintf('%.4f', o$chi2$p.value), '$, $\\chi^2=',sprintf('%.1f',o$chi2$statistic),'$).')
      return(cap)
   }
   if(class(x) == 'SDisc'){
      if(class(y) == 'SDisc'){
         res <- compareModel(x, y)
         cap <- paste2('(',SDPrefix(x), ' ', bestModel(x,1), ')x(', SDPrefix(y),' ',bestModel(y,1), ')')
      }
      else{
         if(is.numeric(m1) && is.numeric(m2)){
            m1 <- bestModel(x, m1)[m1]
            m2 <- bestModel(x, m2)[m2]
         }
         res <- compareModel(x[[m1]], x[[m2]])
         cap <- paste2('(', m1,')x(', m2,')')
      }
      if(latex && !is.null(res))
         texTable(res[["xtab"]], cap=doCaption(res,cap), digits=0, lab=lab)
      else
         return(res)
   }
   else
      cat2('You must provide one or two \'SDisc\' objects for comparison.')
}

`summary.bicTable` <- 
function(object, fun='min', bic='relativeBic', latex=FALSE, lab='bic', fmt='%.2f', ...){
   df <- attr(object,'data')
   x.summary <- bicStats(df, fun=fun, vAgg=bic, fmt=fmt)
   idBestM <- which(df[,'relativeBic'] == 0)
   bestName <- paste(as.matrix(df[idBestM,1:3]), collapse=',', sep=' ')
   if(latex)
      texTable(x.summary, cap=paste2('\\textbf{',SDPrefix(object),'}, model ', bestName, ' shows the \\textbf{highest BIC} score
         over: the repeated random starts, type of model and number of component. '), lab=lab) 
   else
      return(as.matrix(x.summary))
}

`naPattern` <- 
function(x, latex=FALSE, ...){
   if(class(x) == 'SDData' & latex == FALSE)
      return(row.names(SDDataOrig(x)[!(row.names(SDDataOrig(x)) %in% row.names(x)),]))
   else if(class(x) == 'SDData' & latex == TRUE){
      `rateMissingData` <- function(x){
         rate <- length(naPattern(x))/nrow(SDDataOrig(x))*100
         return(sprintf('%.2f', rate)) }
      rowId <- row.names(SDDataOrig(x)[!(row.names(SDDataOrig(x)) %in% row.names(x)),])
      rate <- rateMissingData(x)
      m <- t(apply(is.na(SDDataOrig(x)[rowId, ]), 1, table))[,c('TRUE','FALSE')]
      colnames(m) <- c('isNA','isNotMissing')
      m <- cbind(m, naRate=100*m[,'isNA']/apply(m,1,sum))
      texTable(m, lab=paste2('missing',SDPrefix(x)), cap=paste2( SDPrefix(x), ', index of the observations 
         presenting \\textbf{missing values} along with the number of missings and non-missings; the observations with a
         missing value represent ', rate,'\\% of the available observations.'), ...)
   }
   else
      cat2('A \'SDData\' object must be provided to the function.')
}

`summary.SDTDataElmt` <-
function(object, digits=2, ...){
   fModel <- attr(object,'fModel')
   fName <- attr(object,'fName')
   n <- digits
   if(class(fModel) == 'lm'){
      m <- summary(fModel)
      k <- m$coefficients
      cNames <- c(paste(row.names(k), c('(SE; Pr(>|t|))')),'$R^2$ (adj-$R^2$; N)')
      res <-  matrix(NA, 1, length(cNames), dimnames=list(summary(fModel)$call[[2]],cNames)) 
      for(i in 1:nrow(k))
         res[,cNames[i]] <- sprintf(paste2('%.',n,'f (%.',n,'f; %.1e)'),k[i,1],k[i,2],k[i,4]) 
      res[,ncol(res)] <- sprintf(paste2('%.',n,'f (%.',n,'f; %d)'),m$r.squared,m$adj.r.squared,length(m$residuals))
   }
   else
      res <- matrix(sprintf(paste2('%.', n,'e'),as.numeric(fModel)),1,1,dimnames=list(attr(object,'vName'), fName))
   return(res)
}

`summary.SDTData` <- 
function(object, q=NULL, digits=2, ...){
   res <- NULL 
   tList <- SDTDataElmt(object)
   for(tElmt in tList){
      tRes <- summary(tElmt, digits=digits)
      if(is.null(res))
         res <- tRes
      else{
         cNames <- colnames(tRes)[which(!(colnames(tRes) %in% colnames(res)))]
         res <- cbind(res, matrix(NA, nrow(res), length(cNames), dimnames=list(row.names(res), cNames))) 
         res[row.names(tRes), colnames(tRes)] <- tRes 
      }
   }
   return(res)
}

`summary.SDData` <- 
function(object, q=NULL, latex=FALSE, digits=3, ...){
   if(!is.null(q)){
      lName <- lapply(SDTData(object), function(y){return(names(SDTDataElmt(y)))})
      lName <- lapply(lName, function(y, p=q){return(grep(p, y))})
      cNames <- NULL
      for(i in names(lName))
         if(length(lName[[i]]) > 0)
            cNames <- c(cNames, i)
   }
   else
      cNames <- names(SDTData(object))
   resSummary <- list()
   for(vName in cNames){
      vRes <- summary(SDTData(object)[[vName]], q)
      if(!is.null(vRes))
         resSummary[[vName]] <- vRes
   }
   rName <- unique(unlist(lapply(resSummary, row.names)))
   res <- matrix(NA, length(rName), 0, dimnames=list(rName, list()))
   for(vName in names(resSummary)){
      mSum <- resSummary[[vName]]
      cNames <- unique(c(colnames(mSum), colnames(res)))
      if(length(cNames) != length(colnames(res)))
         res <- cbind(res, matrix(NA, nrow(res), length(cNames)-length(colnames(res)),
            dimnames=list(row.names(res), cNames[which(!(cNames %in% colnames(res)))])))
      res[row.names(mSum), colnames(mSum)] <- mSum
   }
   if(latex)
      texTable(res, cap=paste2('\\textbf{',SDPrefix(object), '} summary of the different data treatments operated on ', 
         'the data.'), lab=paste2('tab:',SDPrefix(object),'DataSummary'), align=paste2(rep('r',ncol(res)+1)), ...)
   else
      return(data.frame(res))
}

`texTable` <- 
function(x, cap='', lab='', digits=NULL, sanitize=TRUE, align=NULL, oddColor='blueLines', longtab=FALSE, ...){
   f <- NULL
   if(!sanitize)
      f <- function(y){y}
   if(!is.null(oddColor))
      cat2(paste2("\\rowcolors{2}{white}{", oddColor,"}"))
   if(longtab)
      print(xtable(x, caption=cap, label=paste2('tab:',lab), digits=digits, align=align, type='latex'),
         table.placement="h!", tabular.environment='longtable', floating=FALSE, sanitize.text.function=f, ...)
   else
      print(xtable(x, caption=cap, label=paste2('tab:',lab), digits=digits, align=align, type='latex'),
         table.placement="h!", sanitize.text.function=f, ...)
}

`summary.SDCModel` <- 
function(object, data, type='oddRatiosB', latex=FALSE, lab='', shortStr=FALSE, ...){
   if(is.null(data) || class(data) != 'SDData')
      cat2('To summary an SDCModel, an SDData-instance must be provided as second parameter (data).')
   else{
      # class <- SDCModel(object)[["z"]]
      class <- print(object, data)
      m <- do.call(type, args=list(data, class=class[,-ncol(class)], ...))
      if(shortStr){
         summaryStr <- function(str){
            cList <- strsplit(strsplit(gsub('and[ ]*','',str),' ')[[1]],'')
            str <- ''
            for(i in 1:length(cList))
               str <- paste2(str,cList[[i]][1],'.')
            return(toupper(str))
         }
         row.names(m) <- sapply(row.names(m), summaryStr)
      }
      set <- SDCModelSettings(object)
      if(latex){
         df <- print(m)
         caption <- paste2(attr(m,'legend'), paste2(' in model ', set[[1]], ',', set[[2]], ',', set[[3]],'.'))
         cell.format <- matrix(rep("\\color{grey}",nrow(df)*ncol(df)), nrow(df), ncol(df))
         if(type=='oddRatiosB' || type=='oddRatios'){
            cell.format[df <  -2] <- "\\color{blue3}"
            cell.format[df < -10] <- "\\color{blue4}"
            cell.format[df >   2] <- "\\color{red3}"
            cell.format[df >  10] <- "\\color{red4}"
            digits <- 2
         }
         else if(type == 'jointDistrib'){
            val <- quantile(abs(df),probs=0.8)
            cell.format[df < -val] <- "\\color{blue3}"
            cell.format[df >  val] <- "\\color{red3}"
            digits <- 0
         }
         else if(type == 'chi2test'){
            cell.format <- matrix(rep("",nrow(df)*ncol(df)), nrow(df), ncol(df))
            cell.format[df[,1]<0.05,1] <- "\\color{red3}"
            cell.format[df[,1]<0.01,1] <- "\\color{red4}"
            digits <- 3
         }
         df <- apply(df, c(1,2), sprintf, fmt=paste2('%.', digits,'f'))
         df <- apply(abind(cell.format, df, along=3), c(1,2), paste, collapse=' ')
         texTable(df, cap=caption, lab=lab, sanitize=FALSE, align=rep('r',ncol(df)+1))
      }
      else 
         return(m)
   }
}

`summary.SDisc` <- 
function(object, q=1, ...){
   if(is.null(q)){
      cat2('- DATA SUMMARY -')
      print(summary(SDData(object)))
      cat2('- DATA TREATMENT SETTINGS-')
      print(SDDataSettings(object))
      cat2('- \'',SDPrefix(object),'\'-SDisc SUMMARY -')
      cat2(' ',SDRInfo(object)[['R.version']]$version.string,',',SDRInfo(object)[['R.version']]$platform)
      cat2(' SDisc ', SDRInfo(object)$otherPkgs[['SubtypeDiscovery']][['Version']],')')
      cat2('- MIXTURE MODELS -')
      print(t(xtabs(formula=~., SDSettings(object)[,1:2])))
      cat2('- BIC TABLE -')
      print(as.data.frame(summary(bicTable(object))))
      cat2('- BEST MODEL -')
      print(data.frame("bestModels"=bestModel(object)))
   }
   else
      return(summary(object[[bestModel(object, n=q)]], data=SDData(object), ...))
}

#
# expDesign:
#
expDesign <- function(x, gr1=NULL, gr2=NULL, tex=FALSE, exclude=NULL){
   if(is.null(gr1) & is.null(gr2))
      gr1 <- colnames(x)
   if(!is.null(gr1) & is.null(gr2))
      m <- t(combn(gr1, 2))
   if(!is.null(gr1) & !is.null(gr2))
      m <- as.matrix(expand.grid(gr1,gr2))
   res <- list()
   for(i in 1:nrow(m)){
      p <- chisq.test(x[,m[i,1]], x[,m[i,2]], simulate.p.value=TRUE, B=2000)$p.value
      res[[i]] <- table(x[,c(m[i,1], m[i,2])], exclude=exclude)
      if(is.null(exclude)){
         row.names(res[[i]])[is.na(row.names(res[[i]]))] <- 'NA'
         colnames(res[[i]])[is.na(colnames(res[[i]]))] <- 'NA'
      }
      names(res)[i] <- paste2(m[i,1],'-', m[i,2])
      # new
      res[[i]] <- res[[i]][hclust(dist(res[[i]]))$order, hclust(dist(t(res[[i]])))$order]
      if(tex){
         cap <- paste2('Joint distribution of \\textbf{', m[i,1], '} and \\textbf{',m[i,2],'}, $p_{\\chi^2}=',
            signif(p, digits=-1),'$.')
         print(xtable(res[[i]], cap=cap), table.placement="!ht",tabular.environment="longtable",
            latex.environments=c("center", "footnotesize"), floating=FALSE)
      }
      else
         print(res[[i]])
   }
   return(invisible(res))
}

`SDiscReportHead` <- function(LO='Report', packages=list('Sweave','amsmath', 'underscore', 'setspace', 'pdflscape',
   'multirow', 'glossaries'=c('toc', 'acronym', 'xindy'), 'multicol', inputenc='latin1', babel='english', 'pdfpages',
   caption=c('small','bf'), 'graphicx','fancyhdr', 'lastpage','longtable','color',xcolor='table',
   geometry=c('left=1.25cm','top=2cm','right=1.25cm','bottom=2cm'), hyperref=c('colorlinks=true', 'citecolor=blueDoc',
   'filecolor=blueDoc', 'linkcolor=blueDoc', 'urlcolor=blueDoc'),'lscape','sectsty','colortbl','wrapfig','array'),
   author=list(name='MyName',email='MyEmail',address='address')){
   cat('\\documentclass{article}')
   for(pId in 1:length(packages)){
      if(names(packages)[pId] != ''){
         cat(paste2('\n\\usepackage[', paste(packages[[pId]], collapse=',')))
         cat(paste2(']{',names(packages)[pId],'}'))
      }
      else
         cat(paste2('\n\\usepackage{',packages[[pId]],'}'))
   }

   ### SDISC STYLE FILE
   cat('\n\\definecolor{blueDoc}{rgb}{0.2156863,0.4941176,0.7215686}')
   cat('\n\\definecolor{blueLines}{rgb}{0.8705882,0.9215686,0.9686275}')
   cat('\n\\definecolor{blue4}{rgb}{0.1294118,0.4000000,0.6745098}')
   cat('\n\\definecolor{blue3}{rgb}{0.2627451,0.5764706,0.7647059}')
   cat('\n\\definecolor{blue2}{rgb}{0.5725490,0.7725490,0.8705882}')
   cat('\n\\definecolor{blue1}{rgb}{0.8196078,0.8980392,0.9764706}')
   cat('%%\\definecolor{white}{rgb}{0.9686275,0.9686275,0.9686275}')
   cat('\n\\definecolor{grey}{rgb}{0.8,0.8,0.8}')
   cat('\n\\definecolor{grey2}{rgb}{0.3215686,0.3215686,0.3215686}')
   cat('\n\\definecolor{red1}{rgb}{0.9921569,0.8588235,0.7803922}')
   cat('\n\\definecolor{red2}{rgb}{0.9568627,0.6470588,0.5098039}')
   cat('\n\\definecolor{red3}{rgb}{0.8392157,0.3764706,0.3019608}')
   cat('\n\\definecolor{red4}{rgb}{0.69803922,0.09411765,0.16862745}')
   cat('\n\\fancyhead{}')
   cat('\n\\fancyfoot{}')
   cat('\n\\renewcommand{\\rmdefault}{phv}')
   cat('\n\\renewcommand{\\sectionmark}[1]{\\markright{#1}{}}')
   cat('\n\\renewcommand{\\headrulewidth}{3pt}')
   cat('\n\\renewcommand{\\footrulewidth}{3pt}')
   cat('\n\\newcommand{\\captionfonts}{\\sffamily\\small\\it}')
   cat('\n\\fancyhead[RO]{\\color{grey2}{\\sffamily\\nouppercase\\leftmark}}')
   cat('\n\\fancyfoot[CO]{ }')
   cat('\n\\fancyfoot[RO]{\\color{grey2}{\\sffamily\\thepage/\\pageref{LastPage}}}')
   if(all(c('name','email','address') %in% names(author))){
      cat(paste2('\n\\fancyfoot[LO]{\\color{grey2}\\sffamily\\href{mailto:',author$email,'}{',author$name,'},',author$address,'}}'))
      cat(paste2('\n\\author{\\href{mailto:',author$email,'}{',author$name,'}'))
   }
   else{
      if(!is.null(author))
         cat(paste2('\n\\fancyfoot[LO]{\\color{grey2}\\sffamily\\href{mailto:',author[[1]]$email,'}{',author[[1]]$name,'}, ',author[[1]]$address,'}'))
      cat('\n\\author{')
      for(a in author)
         cat(paste2('\\href{mailto:',a$email,'}{',a$name,'}\\\\'))
      cat('}')
   }
   cat('\n\\makeglossaries\n \n\\makeatletter \n\\def\\thickhrulefill{\\leavevmode \\leaders \\hrule height
      1pt\\hfill \\kern \\z@} \n\\def\\maketitle{% \n\\null \n\\thispagestyle{empty}% \n\\vskip 1cm
      \n\\begin{flushright} \n\\normalfont\\vspace{7cm}\\fontsize{40}{48}\\selectfont\\sffamily
      \\bfseries\\color{grey2}\\@title\\par \n\\end{flushright} \n\\vfil \n\\begin{flushright} \n\\LARGE \\strut
      \\sffamily\\bfseries\\color{grey2}\\@author \\par \n\\end{flushright} \n\\par \n\\vfil \n\\vfil \n\\null
      \n\\cleardoublepage \n} \n\\makeatother \n\\allsectionsfont{\\usefont{OT1}{phv}{bc}{n}\\selectfont}')
}

`plotPC1` <- function(x, data, ncomp=5, ...){ plotPC(xx=x, xdat=data, xNComp=ncomp, xNPC=1, ...) }
`plotPC2` <- function(x, data, ncomp=5, ...){ plotPC(xx=x, xdat=data, xNComp=ncomp, xNPC=2, ...) }
`plotPC3` <- function(x, data, ncomp=5, ...){ plotPC(xx=x, xdat=data, xNComp=ncomp, xNPC=3, ...) }
`plotPC4` <- function(x, data, ncomp=5, ...){ plotPC(xx=x, xdat=data, xNComp=ncomp, xNPC=4, ...) }
`plotPC5` <- function(x, data, ncomp=5, ...){ plotPC(xx=x, xdat=data, xNComp=ncomp, xNPC=5, ...) }

`plotPC` <- function(xx, xdat, xNComp, xNPC, ...) {
   y <- print(SDData(xdat))
   G <- attr(xx,'model')$label 
   GSet <- levels(factor(G))
   PC <- data.frame(matrix(NA, dim(y)[[1]], length(GSet), dimnames=list(row.names(y),GSet))) 
   varExpl <- data.frame(matrix(NA, xNPC, length(GSet), dimnames=list(list(),GSet)))
   RelVar <- data.frame(matrix(NA, length(GSet), dim(y)[[2]], dimnames=list(GSet, colnames(y))))
   names(PC) <- paste("PC", GSet, sep="")
   for(i in 1:length(GSet)){
      yG <- y[as.character(G) == GSet[i],]
      if(ncol(yG)<=nrow(yG)){
         ySvd <- svd(t(yG), nu=xNComp, nv=xNComp)
         varExpl[1:xNPC, i] <- (ySvd$d[1:xNPC])^2/sum(ySvd$d^2)
         PC[as.character(G) == GSet[i],i] <- ySvd$v[,xNPC]
         tmp <- data.frame(cbind(y, PC=PC[,i]))
         RelVar[i,] <- sapply(colnames(y), function(ww) anova(lm(as.formula(paste2("PC~", ww)), data=tmp))[1,5])
      }
   }
   ### sum up the relative contribution of each variable to 1
   RelVar <- as.matrix(t(RelVar / t(apply(RelVar,1,sum, na.rm=TRUE)))) 
   ### scale down each component to the percentage of variability explained
   RelVar <- RelVar %*% diag(as.numeric(varExpl[xNPC,])) 
   ### go for the barplot
   barplot(RelVar, col=colorRampPalette(brewer.pal(9, 'Pastel1'))(nrow(RelVar)), border=FALSE, axes=FALSE, axisname=FALSE,
      legend.text=row.names(RelVar), args.legend=list(border=FALSE, box.col='transparent'), main=paste2('PC',xNPC,' ~ Predictor'))
   axis(2, at=seq(0,max(varExpl[xNPC,], na.rm=T, length.out=4)), labels=round(seq(0,max(varExpl[xNPC,], na.rm=T,
      length.out=4)),d=2), las=2)
   axis(1, at=1:ncol(RelVar), labels=apply(cbind(GSet,t(round(100*varExpl[xNPC,],d=1))), 1, function(w)paste2(w[1] , ' (',w[2],'%)')), las=2)
}


SDData <- function(x, ...) UseMethod("SDData")
SDisc <- function(x, ...) UseMethod("SDisc")
bicTable <- function(x, ...) UseMethod("bicTable")
SDStability <- function(x, ...) UseMethod("SDStability")

