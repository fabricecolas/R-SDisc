`mclust2_bicaic` <-
function (modelName, loglik, n, d, G, noise = FALSE, equalPro = FALSE,...)
{
   modelName    <- switch(EXPR = modelName, XII = "EII", XXI = "EEI",XXX = "EEE", modelName)
   if (G == 0) {
      if (!noise)
      stop("undefined model")
      nparams <- 1
   }
   else {
         nparams <- nVarParams(modelName, d, G) + G * d
      if (!equalPro)
         nparams <- nparams + (G - 1)
      if (noise)
         nparams <- nparams + 2
   }
   #
   # RETURN VECTOR WITH LOG BIC AND AIC SCORES
   #
   return(c(2 * loglik - nparams * logb(n), 2 * loglik - 2 * nparams))
}

