`summary_mclust2_BIC` <-
function (object, data, G = NULL, modelNames = NULL, ...) 
{
    dimData             <- dim(data)
    oneD                <- is.null(dimData) || length(dimData[dimData > 1]) == 
        1
    if (!oneD && length(dimData) != 2) 
        stop("data must be a vector or a matrix")
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
    initialization      <- attr(object, "initialization")
    hcPairs             <- initialization$hcPairs
    subset              <- initialization$subset
    prior               <- attr(object, "prior")
    control             <- attr(object, "control")
    warn                <- attr(object, "warn")
    oldClass(object)    <- NULL
    attr(object, "prior")               <- attr(object, "warn") <- NULL
    attr(object, "modelNames")          <- attr(object, "oneD") <- NULL
    attr(object, "initialization")      <- attr(object, "control") <- NULL
    d                   <-      if (is.null(dim(data))) 1
                                else ncol(data)
    if (is.null(G)) 
        G               <- dimnames(object)[[1]]
    if (is.null(modelNames)) 
        modelNames      <- dimnames(object)[[2]]
    bestBICs            <- pickBIC(object[as.character(G), modelNames, drop = FALSE], k = 3)
    if(is.na(bestBICs)[1])
    	return(NULL)
    temp                <- unlist(strsplit(names(bestBICs)[1], ","))
    bestModel           <- temp[1]
    G                   <- as.numeric(temp[2])
    if (G == 1) {
        out             <- mvn(modelName = bestModel, data = data, prior = prior)
        ans             <- c(list(bic = bestBICs, classification = rep(1, n), uncertainty = rep(0, n)), out)
        orderedNames    <- c("modelName", "n", "d", "G", "bic", "loglik", "parameters", "z", "classification", "uncertainty")
        return(structure(ans[orderedNames], bestBICvalues = bestBICs, 
                prior = prior, control = control, initialization = initialization, class = "summary_mclust2_BIC"))
    }
    if (is.null(subset)) {
        if (d > 1 || !is.null(hcPairs)) {
            z           <- unmap(hclass(hcPairs, G))
        }
        else {
            z           <- unmap(qclass(data, G))
        }
        out             <- me(modelName = bestModel, data = data, z = z, prior = prior, control = control, warn = warn)
    }
    else {
        if (d > 1 || !is.null(hcPairs)) {
            z           <- unmap(hclass(hcPairs, G))
        }
        else {
            z           <- unmap(qclass(data[subset], G))
        }
        ms              <- mstep(modelName = bestModel, prior = prior, z = z, data = as.matrix(data)[subset, ], control = control, warn = warn)
        es              <- do.call("estep", c(list(data = data), ms))
        out             <- me(modelName = bestModel, data = data, z = es$z, prior = prior, control = control, warn = warn)
    }
    obsNames <- if (is.null(dim(data))) {
        names(data)
    }
    else {
        dimnames(data)[[1]]
    }
    classification      <- map(out$z)
    uncertainty         <- 1 - apply(out$z, 1, max)
    names(classification) <- names(uncertainty) <- obsNames
    ans                 <- c(list(bic = as.vector(bestBICs[1]), classification = classification, uncertainty = uncertainty), out)
    orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", "parameters", "z", "classification", "uncertainty")
    structure(ans[orderedNames], bestBICvalues = bestBICs, prior = prior, 
        control = control, initialization = initialization, class = "summary_mclust2_BIC")
}

