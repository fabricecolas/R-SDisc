\name{mclust2_BIC}
\alias{mclust2_BIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mclust2_BIC(data, G = NULL, modelNames = NULL, prior = NULL, control = emControl(), initialization = list(hcPairs = NULL, subset = NULL, noise = NULL), Vinv = NULL, warn = FALSE, x = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{G}{ ~~Describe \code{G} here~~ }
  \item{modelNames}{ ~~Describe \code{modelNames} here~~ }
  \item{prior}{ ~~Describe \code{prior} here~~ }
  \item{control}{ ~~Describe \code{control} here~~ }
  \item{initialization}{ ~~Describe \code{initialization} here~~ }
  \item{Vinv}{ ~~Describe \code{Vinv} here~~ }
  \item{warn}{ ~~Describe \code{warn} here~~ }
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, G = NULL, modelNames = NULL, prior = NULL, control = emControl(), 
    initialization = list(hcPairs = NULL, subset = NULL, noise = NULL), 
    Vinv = NULL, warn = FALSE, x = NULL, ...) 
{
    #
    # CONTROL OF THE PARAMETERS
    # 
    if (!is.null(x)) {
        if (!missing(prior) || !missing(control) || !missing(initialization) || !missing(Vinv)) 
            stop("only G and modelNames may be specified as arguments when x is supplied")
        #
        # SET LOCAL SCOPE VARIABLES WITH ATTRIBUTES OF THE SAME NAME FROM OBJECT x
        #
        prior           <- attr(x, "prior")
        control         <- attr(x, "control")
        initialization  <- attr(x, "initialization")
        Vinv            <- attr(x, "Vinv")
        warn            <- attr(x, "warn")
    }
    #
    # DATA MATRIX CONTROL (DIMENSION)
    #
    dimData             <- dim(data)
    oneD                <- is.null(dimData) || length(dimData[dimData > 1]) == 1
    if (!oneD && length(dimData) != 2) 
        stop("data must be a vector or a matrix")
    # SET LOCAL SCOPE DIMENSION VARIABLES
    if (oneD) {
        data            <- drop(as.matrix(data))
        n               <- length(data)
        d               <- 1
    }
    else {
        data            <- as.matrix(data)
        n               <- nrow(data)
        d               <- ncol(data)
    }
    #
    # SET DEFAULT MODEL PARAMETERS IF NONE ARE PROVIDED (X IS NULL)
    #
    if (is.null(x)) {
        if (is.null(modelNames)) {
            # IF DIMENSION IS OF LENGTH 1, ONLY TWO AVAILABLE MODELS
            if (d == 1) {
                modelNames <- c("E", "V")
            }
            # MODELS FOR OTHER DIMENSIONS
            else {
                modelNames <- .Mclust$emModelNames
                if (n <= d) {
                  m             <- match(c("EEE", "EEV", "VEV", "VVV"), .Mclust$emModelNames,nomatch = 0)
                  modelNames    <- modelNames[-m]
                }
            }
        }
        # NUMBER OF CLUSTERS TO LOOK FOR, MAY BE A VECTOR OF NUMBERS
        if (is.null(G)) {
            # IF THE MODELING WAS CHOSEN WITHOUT POISSON NOISE, THEN DEFAULT G
            # RANGE FROM 1 TO 9, ELSE 0 IS ALLOWED FOR NOISE MODELING ONLY (SEE
            # DOCUMENTATION)
            G           <- if (is.null(initialization$noise)) 
                1:9
            else 0:9
        }
        else {
            G           <- sort(as.numeric(G))
        }
        # ORDERED CLUSTER NUMBERS AND MODELS
        Gall            <- G
        Mall            <- modelNames
    }
    #
    # CONTROL PROVIDED PARAMETERS
    #
    else {
        Glabels         <- dimnames(x)[[1]]
        Mlabels         <- dimnames(x)[[2]]
        if (is.null(G)) 
            G           <- Glabels
        if (is.null(modelNames)) 
            modelNames  <- Mlabels
        # 
        # ARE THE PROVIDED PARAMETERS POSSIBLE, MATCH CONTROL
        #
        Gmatch                  <- match(as.character(G), Glabels, nomatch = 0)
        Mmatch                  <- match(modelNames, Mlabels, nomatch = 0)
        if (all(Gmatch) && all(Mmatch)) {
            attr(x, "G")                <- as.numeric(G)
            attr(x, "modelNames")       <- modelNames
            attr(x, "returnCodes")      <- attr(x, "returnCodes")[as.character(G), modelNames, drop = FALSE]
            return(x[as.character(G), modelNames, drop = FALSE])
        }
        Gall            <- sort(as.numeric(unique(c(as.character(G), Glabels))))
        Mall            <- unique(c(modelNames, Mlabels))
    }
    if (any(as.numeric(G)) < 0) {
        if (is.null(initialization$noise)) {
            stop("G must be positive")
        }
        else {
            stop("G must be nonnegative")
        }
    }
    #
    # IF DIM OF THE DATA IS ONE, THEN A MATCHING ON THE PROVIDED MODELS TO
    # COMPUTE ON IS DONE E, V
    #
    if (d == 1 && any(nchar(modelNames) > 1)) {
        Emodel          <- any(sapply(modelNames, function(x) charmatch("E", x, nomatch = 0)[1]) == 1)
        Vmodel          <- any(sapply(modelNames, function(x) charmatch("V", x, nomatch = 0)[1]) == 1)
        modelNames      <- c("E", "V")[c(Emodel, Vmodel)]
    }
    #
    # INITIALIZE VARIABLES
    #
    l                   <- length(Gall)
    m                   <- length(Mall)
    EMPTY               <- -.Machine$double.xmax
    BIC  <- AIC  <- SIG <- RET <- matrix(EMPTY, nrow = l, ncol = m, dimnames = list(as.character(Gall), as.character(Mall)))
    if (!is.null(x)) {
        AIC[dimnames(x)[[1]], dimnames(x)[[2]]]         <- x
        BIC[dimnames(x)[[1]], dimnames(x)[[2]]]         <- x
        SIG[dimnames(x)[[1]], dimnames(x)[[2]]]         <- x
        RET[dimnames(x)[[1]], dimnames(x)[[2]]]         <- attr(x, "returnCodes")
        AIC                                             <- AIC[as.character(G), modelNames, drop = FALSE]
        BIC                                             <- BIC[as.character(G), modelNames, drop = FALSE]
        RET                                             <- RET[as.character(G), modelNames, drop = FALSE]
        SIG                                             <- SIG[as.character(G), modelNames, drop = FALSE]
    }
    G                   <- as.numeric(G)
    Glabels             <- as.character(G)
    Gout                <- G
    stats				<- NULL
    #
    # LIST TO STORE THE RESULTS OF THE MULTIVARIATE NORMAL FITS
    #
    outputList          <- NULL
    #
    if (is.null(initialization$noise)) {
    	#
    	# HERE!!!
    	#
        if (G[1] == 1) {
            for (mdl in modelNames[BIC["1", ] == EMPTY]) {
                out             <- mvn(modelName = mdl, data = data, prior = prior)
                out.scores      <- mclust2_bicaic(modelName = mdl, loglik = out$loglik, n = n, d = d, G = 1, equalPro = FALSE)
                BIC["1", mdl]   <- out.scores[1] 
                AIC["1", mdl]   <- out.scores[2] 
                RET["1", mdl]   <- attr(out, "returnCode")
                outputList    	<- c(outputList,list(out))                
            }
            if (l == 1) {
                BIC[BIC == EMPTY] <- NA
#                return(structure(BIC, BIC=BIC, AIC=AIC,  SIG=SIG, NSignificant = NSignificant, G = G, 
                return(structure(BIC, BIC=BIC, AIC=AIC,  SIG=SIG, G = G, 
                				modelNames = modelNames, prior = prior, control = control, 
                                initialization = list(hcPairs = initialization$hcPairs, subset = initialization$subset), 
                                warn = warn, n = n, d = d, oneD = oneD, returnCodes = RET, out = outputList, 
                                class = "mclust2_BIC"))
            }
            G                   <- G[-1]
            Glabels             <- Glabels[-1]
        }
        if (is.null(initialization$subset)) {
            if (is.null(initialization$hcPairs)) {
                if (d != 1) {
                  if (n > d) {
                    hcPairs     <- hc(modelName = "VVV", data = data)
                  }
                  else {
                    hcPairs     <- hc(modelName = "EII", data = data)
                  }
                }
                else {
                  hcPairs       <- NULL
                }
            }
            else hcPairs        <- initialization$hcPairs
            if (d > 1 || !is.null(hcPairs)) 
                clss            <- hclass(hcPairs, G)
            #
            # AFTER FITTING A SINGLE MULTIVARIATE MODEL ON THE DATA (1-CLUSTER, PREVIOUS LOOP),
            # CALCULATE BY EM ALGORITHM, EVENTUALY INITIALIZED BY HIERARCHICAL CLUSTERING, THE 
            # STATISTIC FOR "G" CLUSTERS IN THE DATA.
            #
            for (g in Glabels) {
                if (d > 1 || !is.null(hcPairs)) {
                  z             <- unmap(clss[, g])
                }
                else {
                  z             <- unmap(qclass(data, as.numeric(g)))
                }
                for (modelName in modelNames[BIC[g, ] == EMPTY]) {
                  out           <- me(modelName = modelName, data = data, z = z, prior = prior, control = control, warn = warn)
                  out.scores    <- mclust2_bicaic(modelName = modelName, loglik = out$loglik, n = n, d = d, G = as.numeric(g), equalPro = control$equalPro)
                  BIC[g, modelName] <- out.scores[1] 
                  AIC[g, modelName] <- out.scores[2] 
                  RET[g, modelName] <- attr(out, "returnCode")
                  if(RET[g, modelName] == 0){
                  	stats <- NULL
               		# stats 		<- mclust2_statistics(mclust2_mergeSibpairs(matrix(map(out$z),
                  	#							dimnames=list(dimnames(out$z)[[1]], NULL))))
             	    # SIG[g, modelName]	<- stats["Significant?","T"]
                  }
                  else{
                  	stats				<- NULL
                  	SIG[g, modelName]	<- 0
                  	}
                    
                  outputList    <- c(outputList,list(out))
                }
            }
        }
        else {
            if (is.logical(initialization$subset)) 
                initialization$subset <- (1:n)[initialization$subset]
            if (is.null(initialization$hcPairs)) {
                if (d != 1) {
                  if (n > d) {
                    hcPairs     <- hc(modelName = "VVV", data = data[initialization$subset, ])
                  }
                  else {
                    hcPairs     <- hc(modelName = "EII", data = data[initialization$subset, ])
                  }
                }
                else {
                  hcPairs <- NULL
                }
            }
            else hcPairs <- initialization$hcPairs
            if (d > 1 || !is.null(hcPairs)) 
                clss <- hclass(hcPairs, G)
            for (g in Glabels) {
                if (d > 1 || !is.null(hcPairs)) {
                  z <- unmap(clss[, g])
                }
                else {
                  z <- unmap(qclass(data[initialization$subset], as.numeric(g)))
                }
                dimnames(z) <- list(as.character(initialization$subset), NULL)
                for (modelName in modelNames[!is.na(BIC[g, ])]) {
                  ms 			<- mstep(modelName = modelName, z = z, data = as.matrix(data)[initialization$subset, ], 
                  						prior = prior, control = control, warn = warn)
                  es 			<- do.call("estep", c(list(data = data, warn = warn), ms))
                  out 			<- me(modelName = modelName, data = data, z = es$z, prior = prior, control = control, warn = warn)
                  out.scores    <- mclust2_bicaic(modelName = modelName, loglik = out$loglik, n = n, d = d, G = as.numeric(g), 
                  						equalPro = control$equalPro)
                  BIC[g, modelName] <- out.scores[1] 
                  AIC[g, modelName] <- out.scores[2] 
                  RET[g, modelName] <- attr(out, "returnCode")
                  if(RET[g, modelName] == 0){
                  		stats 	<- NULL
#               		 stats 		<- mclust2_statistics(mclust2_mergeSibpairs(data,matrix(map(out$z),
#                  								dimnames=list(dimnames(out$z)[[1]], NULL))))
             	     # SIG[g, modelName]	<- stats["Significant?","T"]
                  }
                  else{
                  	stats			<- NULL
                  	SIG[g, modelName]	<- 0
                  	}
                  outputList    <- c(outputList,list(out))
                }
            }
        }
    }
    structure(BIC, BIC=BIC, AIC=AIC, SIG=SIG, G = Gout, modelNames = modelNames, prior = prior, control = control, 
        initialization = list(hcPairs = hcPairs, subset = initialization$subset, noise = initialization$noise), 
        Vinv = Vinv, warn = warn, n = n, d = d, oneD = oneD, 
        returnCodes = RET, out = outputList,class = "mclust2_BIC")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
