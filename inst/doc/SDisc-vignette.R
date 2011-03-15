###################################################
### chunk number 1: dir
###################################################
#line 23 "SDisc-vignette"
options(continue="   ")
options(width=85)
dir0 <- dir(recursive=TRUE)


###################################################
### chunk number 2: installRSDisc eval=FALSE
###################################################
## #line 69 "SDisc-vignette"
## install.packages('SDisc', dep=TRUE)


###################################################
### chunk number 3: loadRSDisc
###################################################
#line 72 "SDisc-vignette"
library(SDisc)


###################################################
### chunk number 4: normdep
###################################################
#line 96 "SDisc-vignette"
set.seed(6015)
Eps <- runif(50)
time <- sample(1:5, 50, replace=TRUE)
normdep.df <- matrix(c(rnorm(50), time, 2*time+Eps), 50, 3, dimnames=list(list(),c("norm","t","v")))
normdep.df[sample(1:nrow(normdep.df))[1:5]] <- NA


###################################################
### chunk number 5: SDDataSettings
###################################################
#line 105 "SDisc-vignette"
normdep.set <- SDDataSettings(normdep.df)
normdep.set[,'tFun'] <- c('mean sd','', 'lm(v~t)')


###################################################
### chunk number 6: SDData
###################################################
#line 109 "SDisc-vignette"
normdep <- SDData(normdep.df, settings=normdep.set, prefix='normdep')


###################################################
### chunk number 7: 
###################################################
#line 112 "SDisc-vignette"
SDDataSettings(normdep, latex=TRUE)


###################################################
### chunk number 8: 
###################################################
#line 117 "SDisc-vignette"
naPattern(normdep, latex=TRUE)
print(normdep, rseed=6013, latex=TRUE)
plot(normdep, latex=TRUE)
summary(normdep, q='mean|sd', latex=TRUE)
summary(normdep, q='lm', latex=TRUE)


###################################################
### chunk number 9: 
###################################################
#line 127 "SDisc-vignette"
set.seed(6016)
Eps <- runif(30)
time <- sample(1:5, 30, replace=TRUE)
normdep.df2 <- matrix(c(rnorm(30), time, 2*time+Eps), 30, 3, dimnames=list(list(), c("norm","t","v")))
normdep2 <- predict(normdep, newdata=normdep.df2, prefix='normdep2')


###################################################
### chunk number 10: 
###################################################
#line 134 "SDisc-vignette"
summary(normdep2, q='lm', latex=TRUE, sanitize=FALSE)
summary(normdep2, q='mean|sd', latex=TRUE)


###################################################
### chunk number 11: 
###################################################
#line 148 "SDisc-vignette"
class(normdep)
normdep <- SDisc(normdep, settings=normdep.set, prefix='normdep', 
   cFunSettings=list(modelName=c("EII", "VII", "VEI","VVI"), G=3:5, rseed=6013:6023))
class(normdep)
class(SDData(normdep))
class(SDDataSettings(normdep))


###################################################
### chunk number 12:  eval=FALSE
###################################################
## #line 157 "SDisc-vignette"
## normdep <- SDisc(normdep)


###################################################
### chunk number 13: SDiscBictable
###################################################
#line 163 "SDisc-vignette"
summary(bicTable(normdep), latex=TRUE)
print(bicTable(normdep), modelName='VII', G=4, latex=TRUE)


###################################################
### chunk number 14: SDiscCompare
###################################################
#line 173 "SDisc-vignette"
print(normdep, latex=TRUE)
print(normdep,  m1=1, m2=bestModel(normdep, modelName='VII', G=4)[1], latex=TRUE)


###################################################
### chunk number 15: SDiscCharacterize
###################################################
#line 179 "SDisc-vignette"
plot(normdep, latex=TRUE)
summary(normdep, q=1, latex=TRUE)


###################################################
### chunk number 16: 
###################################################
#line 192 "SDisc-vignette"
sessionInfo()


###################################################
### chunk number 17: cleanup
###################################################
#line 195 "SDisc-vignette"
dir1 <- dir(recursive=TRUE)
dir1 <- dir1[grep('pdf', dir1, invert=T)]
file.remove(dir1[which(!(dir1 %in% dir0))])
file.remove(c('normdep2/figures','normdep2/tables'))
file.remove(c('normdep2'))


