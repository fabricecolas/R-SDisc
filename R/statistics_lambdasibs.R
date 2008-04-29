`statistics_lambdasibs` <-
function(data,class)
{
   class <- map(class)
   b                    <- cbind(data[["orig"]],class = class)
   b                    <- reshape(b[,c("family","member","class")],timevar="member",idvar="family",direction="wide")
   b                    <- b[,-1]
   # MAKE A DATA STRUCTURE FROM PREVIOUSLY RESHAPED 'B' SPARSE MATRIX 
   # . IN THE 1ST COLUMN THE PROBANT, 
   # . IN THE 2ND ITS SIBLING 
   # . ALL THE NA'S ARE REMOVED
   for(j in 1:2)
      for(i in 1:nrow(b))
         while(is.na(b[i,j]))
            b[i,j:(ncol(b)-1)] <- b[i,(j+1):ncol(b)]
   b <- b[,1:2]
   class_set            <- levels(factor(as.matrix(class)))
   sibship_size		<- dim(b)[[1]]
   rowNames		<- c("Counts","P_cl","var_cl","P_cl_s1","P_cl_s2","P_cl_s1_s2",
                                          "var_cl_s1_s2","cov_P_cl_P_cl_s1_s2","lambdaSib","var_lambda",
                                          "se_lambda","ci_lambda","Significant?")
   s                    <- as.data.frame(matrix(0,length(rowNames),length(class_set)+1))
   # SET DIM NAMES FOR THE RESULT TABLE 'S'
   dimnames(s)[[1]]     <- rowNames
   dimnames(s)[[2]]     <- c(class_set,"Total")
   # COUNT THE NUMBER OF INDIVIDUALS IN EACH CLASS
   s["Counts",names(table(class))] <- table(class)
   s["Counts","Total"]  <- sum(s["Counts",names(table(class))])
   # P_cl          := CLASS PREVALENCE
   s["P_cl",]           <- s["Counts",]/s["Counts","Total"]
   # var_cl        := VAR( CLASS PREVALENCE ESTIMATE)
   s["var_cl",]         <- s["P_cl",]*(1-s["P_cl",])/s["Counts","Total"]
   # P_cl_s1       := CLASS PREVALENCE FOR SIB_1
   s["P_cl_s1",]        <- c(table(factor(b[,1],levels=levels(as.factor(class)))),sibship_size)/sibship_size
   # P_cl_s2       := CLASS PREVALENCE FOR SIB_2
   s["P_cl_s2",]        <- c(table(factor(b[,2],levels=levels(as.factor(class)))),sibship_size)/sibship_size
   # P_cl_s1_s2    := P ( s_i^1 AND s_i^2 )
   for(i in as.numeric(class_set))
            s["P_cl_s1_s2",as.character(i)] <- dim(b[b[,1] == i & b[,2] == i,])[[1]]/sibship_size
   # var_cl_s1_s2  := VAR (P ( s_i^1 AND s_i^2 ))
   s["var_cl_s1_s2",] <- s["P_cl_s1_s2",]*(1-s["P_cl_s1_s2",])/sibship_size
   # lambdasib     := ESTIMATION OF THE RISK RATIO
   s["lambdaSib",]    <- s["P_cl_s1_s2",]/s["P_cl",]^2
   # "cov_P_cl_P_cl_s1_s2",
   s["cov_P_cl_P_cl_s1_s2",] <- s["P_cl_s1_s2",]*(1-s["P_cl",])/s["Counts","Total"]
   # var_lambda
   s["var_lambda",]   <- (1/s["P_cl",]^4)*(s["var_cl_s1_s2",]
                                    -4*s["cov_P_cl_P_cl_s1_s2",]*s["P_cl_s1_s2",]/s["P_cl",]
                                    +4*s["var_cl",]*(s["P_cl_s1_s2",]/s["P_cl",])^2)
   # se_lambda
   s["se_lambda",]   <- sqrt(s["var_lambda",]) 
   # ci_lambda
   s["ci_lambda",]   <- 1.96*s["se_lambda",]
   s["Significant?",]<- (s["lambdaSib",]-s["ci_lambda",])>1
   s["Significant?","Total"]<- sum(s["Significant?",1:length(class_set)])
   #
   return(list(out=t(s[,1:length(class_set)])))
}

