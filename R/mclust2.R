`mclust2` <-
function (data, G = NULL, modelNames = NULL, prior = NULL, control = emControl(), initialization = NULL, warn = FALSE, ...) 
{
    mc                  <- match.call(expand.dots = FALSE)
    mc[[1]]             <- as.name("mclust2_BIC")
    Bic                 <- eval(mc, parent.frame())
    G                   <- attr(Bic, "G")
    modelNames          <- attr(Bic, "modelNames")
    Sumry               <- summary_mclust2_BIC(Bic, data, G = G, modelNames = modelNames)
    if (!(length(G) == 1)) {
        bestG           <- length(unique(Sumry$cl))
        if (bestG == max(G)) 
            warning("optimal number of clusters occurs at max choice")
        else if (bestG == min(G)) 
            warning("optimal number of clusters occurs at min choice")
    }
    attr(Bic, "n")      <- attr(Bic, "warn") <- NULL
    attr(Bic, "initialization") <- attr(Bic, "control") <- NULL
    attr(Bic, "d")      <- attr(Bic, "returnCodes") <- attr(Bic, "class") <- NULL
    oldClass(Sumry)     <- NULL
    Sumry$bic           <- Sumry$bic[1]
    ans                 <- c(list(out = Bic), Sumry)
    orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", "parameters", "z", "classification", "uncertainty","out")
    structure(ans[orderedNames], class = "mclust2_mclust2_Mclust")
}

