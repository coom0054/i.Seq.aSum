### V1&V2: using permutations to calculate p-values
### for sequentially selected final "grouping" model like the aSum test;
### V2 differs from V1 and the aSum test in its SNP selection:
###  assign each SNP to a positive or negative or NULL group, 
###  as indicated by SNPstatus=1 or -1 or 0, while V1 only assign to
###   positive or negative group with SNPstatus=1 or -1.
#### A key: use the Score test to compare the models to SPEED up since
####        aymptotically the score test and LRT are the same!

library(RcppEigen)
library(inline)
library(compiler)
library(glmnet)
library(rrBLUP)
library(hglm)
library(nlme)

####INPUT:   Y: binary outcome; = 0 or 1
####         Xg: SNPs;  nObs by nSNP
####         Xe: SNPs;  nObs by nEnviro
####         true.G and true.G and true.GE; vectors of the correct allocation
####         null.lr; a vector of Likelihood Ratio Test statistics to compare our interaction model to
####OUTPUT: a list of test stats, p-values and SNPstatus and EnviroStatus for two versions of
####        the test: V1 does NOT do SNP selection while V2 does; both 
####        share a common reg coef for the two non-null groups
####        and allows a separate intercept term.
############Sequential allocation (-1,1 or 0) marginal model:
al3.M <- function(Y,Xfix,X,link=gaussian()){
	nPred <- ncol(X); n <- length(Y)
	signs1 <- signs2 <- signs3 <- rep(1,nPred)
	ones <- rep(1,length(Y))
	res1 <- glm.fit(cbind(ones,Xfix,X%*%signs1),Y,family=link)$dev
	###seq check whether to move each predictor to another group:
	for (k in 1:nPred){
		signs2 <- signs3 <- signs1   #reset temp signs
		signs2[k] <- -1    #switch kth predictor to negative group
		signs3[k] <- 0    #switch kth predictor to negative group
		res2 <- glm.fit(cbind(ones,Xfix,X%*%signs2),Y,family=link)$dev
		res3 <- glm.fit(cbind(ones,Xfix,X%*%signs3),Y,family=link)$dev
 		if (res2 <= res1 && res2 <= res3) {
    		signs1 <- signs2
    		res1 <- res2
    	}
    	if (res3 <= res1 && res3 <= res2) {
    		signs1 <- signs3
    		res1 <- res3
    	}
	}
	list(U=res1,signs=signs1)
}	

allocate3.M <- cmpfun(al3.M)	


############Sequential allocation (-1,1 or 0) marginal model with separate betas:
al3.Ib <- function(Y,Xfix,X,link=gaussian()){
	nPred <- ncol(X); n <- length(Y)
	signs1 <- signs2 <- signs3 <- rep(0,nPred)
	ones <- rep(1,length(Y))
	res1 <- glm.fit(cbind(ones,Xfix),Y,family=link)$dev
	###seq check whether to move each predictor to another group:
	for (k in 1:nPred){
		signs2 <- signs3 <- signs1   #reset temp signs
		signs2[k] <- -1    #switch kth predictor to negative group
		signs3[k] <- 1    #switch kth predictor to negative group
		if (sum(signs2 == 1) == 0)
			Xa <- cbind(ones,Xfix,X%*%(signs2 == -1))
		if (sum(signs2 == 1) > 0)
			Xa <- cbind(ones,Xfix,X%*%(signs2 == -1),X%*%(signs2 == 1))
		res2 <- glm.fit(Xa,Y,family=link)$dev
		if (sum(signs3 == -1) == 0)
			Xa <- cbind(ones,Xfix,X%*%(signs3 == 1))
		if (sum(signs3 == -1) > 0)
			Xa <- cbind(ones,Xfix,X%*%(signs3 == -1),X%*%(signs3 == 1))
		res3 <- glm.fit(Xa,Y,family=link)$dev
 		if (res2 <= res1 && res2 <= res3) {
    		signs1 <- signs2
    		res1 <- res2
    	}
    	if (res3 <= res1 && res3 <= res2) {
    		signs1 <- signs3
    		res1 <- res3
    	}
	}
	list(U=res1,signs=signs1)
}	

allocate3.Ibeta <- cmpfun(al3.Ib)	

############Sequential allocation (-1,1 or 0) joint model:
al3.J <- function(Y,Xfix,Xg,Xe,link=gaussian()){
	nSNP <- ncol(Xg); nE <- ncol(Xe); n <- length(Y)
	gsigns1 <- gsigns2 <- rep(1,nSNP)
	esigns1 <- esigns2 <- rep(1,nE)
	gscore <- Xg%*%gsigns1
	escore <- Xe%*%esigns1 
	ones <- rep(1,length(Y))
	res1 <- logLik(glm(Y~cbind(ones,Xfix,gscore,escore),family=link))
	###seq check whether to move each predictor to another group:
	for (k in 1:nSNP){
		gsigns2 <- gsigns3 <- gsigns1   #reset temp signs
		gsigns2[k] <- -1    #switch kth predictor to negative group
		gsigns3[k] <- 0    #switch kth predictor to negative group
  		gscore1 <- Xg%*%gsigns2
  		gscore2 <- Xg%*%gsigns3
  		res2 <- logLik(glm(Y~cbind(ones,Xfix,gscore1,escore),family=link))
  		res3 <- logLik(glm(Y~cbind(ones,Xfix,gscore2,escore),family=link))
  		#cat("Current Gsigns: ",gsigns1,res1," -1: ",gsigns2,res2," 0: ",gsigns3,res3,"\n")
  		if (res2 >= res1 && res2 >= res3) {
    		gsigns1 <- gsigns2
    		res1 <- res2
    	}
    	if (res3 >= res1 && res3 >= res2) {
    		gsigns1 <- gsigns3
    		res1 <- res3
    	}
	}
	gscore <- Xg%*%gsigns1
	for (k in 1:nE){
		esigns2 <- esigns3 <- esigns1   #reset temp signs
		esigns2[k] <- -1    #switch kth predictor to negative group
		esigns3[k] <- 0    #switch kth predictor to negative group
  		escore1 <- Xe%*%esigns2
  		escore2 <- Xe%*%esigns3
  		res2 <- logLik(glm(Y~cbind(ones,Xfix,gscore,escore1),family=link))
  		res3 <- logLik(glm(Y~cbind(ones,Xfix,gscore,escore2),family=link))
  		#cat("Current Esigns: ",esigns1,res1," -1: ",esigns2,res2," 0: ",esigns3,res3,"\n")
  		if (res2 >= res1 && res2 >= res3) {
    		esigns1 <- esigns2
    		res1 <- res2
    	}
    	if (res3 >= res1 && res3 >= res2) {
    		esigns1 <- esigns3
    		res1 <- res3
    	}
	}
	list(U=res1,gsigns=gsigns1,esigns=esigns1)
}

allocate3.J <- cmpfun(al3.J)	

#########################################################################
### FUNCTION FOR FINDING ALLOCATION OF GE INTERACTION
############Sequential allocation (-1,1 or 0) 
#########################################################################
al3.reI <- function(Y,Xg,Xe){
	n <- length(Y)
	Xge = NULL
	for (k in 1:length(Y)){
		Xge = rbind(Xge,kronecker(Xg[k,],Xe[k,]))
	}
	nPred <- ncol(Xge)
	signs1 <- signs2 <- signs3 <- rep(0,nPred)
	ones <- rep(1,n)
	X <- cbind(ones,Xe)
	res1 <- -mixed.solve(Y,Z=Xg,X=X,method="ML")$LL
	###seq check whether to move each predictor to another group:
	for (k in 1:nPred){
		print(k)
		signs2 <- signs3 <- signs1   #reset temp signs
		signs2[k] <- -1    #switch kth predictor to negative group
		signs3[k] <- 1    #switch kth predictor to negative group
		X <- cbind(ones,Xe,Xge%*%signs2)
		res2 <- -mixed.solve(Y,Z=Xg,X=X,method="ML")$LL
		X <- cbind(ones,Xe,Xge%*%signs3)
		res3 <- -mixed.solve(Y,Z=Xg,X=X,method="ML")$LL
		cat(signs1,res1,res2,res3,"\n")
 		if (res2 <= res1 && res2 <= res3) {
    		signs1 <- signs2
    		res1 <- res2
    	}
    	if (res3 <= res1 && res3 <= res2) {
    		signs1 <- signs3
    		res1 <- res3
    	}
	}
	list(U=res1,signs=signs1)
}

allocate3.reI <- cmpfun(al3.reI)

#########################################################################
### FUNCTION FOR FINDING ALLOCATION OF GE INTERACTION
############Sequential allocation (-1,1 or 0) with different betas for each group
#########################################################################
al3.reIbeta <- function(Y,Xg,Xe){
	n <- length(Y)
	Xge = NULL
	for (k in 1:length(Y)){
		Xge = rbind(Xge,kronecker(Xg[k,],Xe[k,]))
	}
	nPred <- ncol(Xge)
	signs1 <- signs2 <- signs3 <- rep(0,nPred)
	ones <- rep(1,n)
	X <- cbind(ones,Xe)
	res1 <- -mixed.solve(Y,Z=Xg,X=X,method="ML")$LL
	###seq check whether to move each predictor to another group:
	for (k in 1:nPred){
		print(k)
		signs2 <- signs3 <- signs1   #reset temp signs
		signs2[k] <- -1    #switch kth predictor to negative group
		signs3[k] <- 1    #switch kth predictor to negative group
		if (sum(signs2 == 1) == 0)
			X <- cbind(ones,Xe,Xge%*%(signs2 == -1))
		if (sum(signs2 == 1) > 0)
			X <- cbind(ones,Xe,Xge%*%(signs2 == -1),Xge%*%(signs2 == 1))
		res2 <- -mixed.solve(Y,Z=Xg,X=X,method="ML")$LL
		if (sum(signs3 == -1) == 0)
			X <- cbind(ones,Xe,Xge%*%(signs3 == 1))
		if (sum(signs3 == -1) > 0)
			X <- cbind(ones,Xe,Xge%*%(signs3 == -1),Xge%*%(signs3 == 1))
		res3 <- -mixed.solve(Y,Z=Xg,X=X,method="ML")$LL
		cat(signs1,res1,res2,res3,"\n")
 		if (res2 <= res1 && res2 <= res3) {
    		signs1 <- signs2
    		res1 <- res2
    	}
    	if (res3 <= res1 && res3 <= res2) {
    		signs1 <- signs3
    		res1 <- res3
    	}
	}
	list(U=res1,signs=signs1)
}

allocate3.reIbeta <- cmpfun(al3.reIbeta)

#########################################################################
### FUNCTION FOR FINDING ALLOCATION OF GE INTERACTION using random effect with R package "nlme"
############Sequential allocation (-1,1 or 0) 
#########################################################################
al3.lmeI <- function(Y,Xg,Xe,Xfix){
	n <- length(Y)
	dummyId <- factor(rep(1,n))
	Z.block <- list(dummyId=pdIdent(~-1+Xg))
	Xge = NULL
	for (k in 1:length(Y)){
		Xge = rbind(Xge,kronecker(Xg[k,],Xe[k,]))
	} 
	nPred <- ncol(Xge)
	signs1 <- signs2 <- signs3 <- rep(0,nPred)
	X <- cbind(Xfix,Xe)
	res1 <- -lme(Y~1+X,random=Z.block,method="ML")$logLik
	###seq check whether to move each predictor to another group:
	for (k in 1:nPred){
		#print(k)
		signs2 <- signs3 <- signs1   #reset temp signs
		signs2[k] <- -1    #switch kth predictor to negative group
		signs3[k] <- 1    #switch kth predictor to negative group
		X <- cbind(Xfix,Xe,Xge%*%signs2)
		res2 <- -lme(Y~1+X,random=Z.block,method="ML")$logLik
		X <- cbind(Xfix,Xe,Xge%*%signs3)
		res3 <- -lme(Y~1+X,random=Z.block,method="ML")$logLik
		#cat(signs1,res1,res2,res3,"\n")
 		if (res2 <= res1 && res2 <= res3) {
    		signs1 <- signs2
    		res1 <- res2
    	}
    	if (res3 <= res1 && res3 <= res2) {
    		signs1 <- signs3
    		res1 <- res3
    	}
	}
	X <- cbind(Xfix,Xe,Xge%*%signs1)
	final.fit <- -lme(Y~1+X,random=Z.block,method="ML")$logLik
	X <- cbind(Xfix,Xe)
	null.fit <- -lme(Y~1+X,random=Z.block,method="ML")$logLik 
	list(LRT=2*(null.fit-final.fit),signs=signs1)
}

allocate3.lmeI <- cmpfun(al3.lmeI)

#########################################################################
### FUNCTION FOR FINDING ALLOCATION OF GE INTERACTION using ridge penalty with R package "penalized"
############Sequential allocation (-1,1 or 0) 
#########################################################################
al3.ridI <- function(Y,Xg,Xe){
	n <- length(Y)
	Xge = NULL
	for (k in 1:length(Y)){
		Xge = rbind(Xge,kronecker(Xg[k,],Xe[k,]))
	}
	nPred <- ncol(Xge)
	signs1 <- signs2 <- signs3 <- rep(0,nPred)
	ones <- rep(1,n)
	X <- cbind(ones,Xe)
	fit <- optL2(Y,penalized=Xg,unpenalized=X,approximate=TRUE)
	res1 <- -fit$cvl
	###seq check whether to move each predictor to another group:
	for (k in 1:nPred){
		print(k)
		signs2 <- signs3 <- signs1   #reset temp signs
		signs2[k] <- -1    #switch kth predictor to negative group
		signs3[k] <- 1    #switch kth predictor to negative group
		X <- cbind(ones,Xe,Xge%*%signs2)
		res2 <- -cvl(Y,penalized=Xg,unpenalized=X,lambda2 = fit$lambda,approximate=TRUE)$cvl
		X <- cbind(ones,Xe,Xge%*%signs3)
		res3 <- -cvl(Y,penalized=Xg,unpenalized=X,lambda2 = fit$lambda,approximate=TRUE)$cvl
		#cat(signs1,res1,res2,res3,"\n")
 		if (res2 <= res1 && res2 <= res3) {
    		signs1 <- signs2
    		res1 <- res2
    	}
    	if (res3 <= res1 && res3 <= res2) {
    		signs1 <- signs3
    		res1 <- res3
    	}
	}
	list(U=res1,signs=signs1)
}

allocate3.ridgeI <- cmpfun(al3.ridI)

#########################################################################
### FUNCTION FOR FINDING ALLOCATION OF GE INTERACTION using random effect with R package "hglm"
############Sequential allocation (-1,1 or 0) 
#########################################################################
al3.reI2 <- function(Y,Xg,Xe,Xfix){
	n <- length(Y)
	Xge = NULL
	for (k in 1:length(Y)){
		Xge = rbind(Xge,kronecker(Xg[k,],Xe[k,]))
	}
	nPred <- ncol(Xge)
	signs1 <- signs2 <- signs3 <- rep(0,nPred)
	X <- cbind(Xfix,Xe)
	res1 <- hglm(X=X,y=Y,Z=Xg,calc.like=T)
	if (!is.null(res1$likelihood)){
		res1 <- -res1$likelihood$pbvh	
		###seq check whether to move each predictor to another group:
		for (k in 1:nPred){
			#print(k)
			signs2 <- signs3 <- signs1   #reset temp signs
			signs2[k] <- -1    #switch kth predictor to negative group
			signs3[k] <- 1    #switch kth predictor to negative group
			X <- cbind(Xfix,Xe,Xge%*%signs2)
			res2 <- -hglm(X=X,y=Y,Z=Xg,calc.like=T)$likelihood$pbvh
			X <- cbind(Xfix,Xe,Xge%*%signs3)
			res3 <- -hglm(X=X,y=Y,Z=Xg,calc.like=T)$likelihood$pbvh
			#cat(signs1,res1,res2,res3,"\n")
	 		if (res2 <= res1 && res2 <= res3) {
	    		signs1 <- signs2
	    		res1 <- res2
	    	}
	    	if (res3 <= res1 && res3 <= res2) {
	    		signs1 <- signs3
	    		res1 <- res3
	    	}
		}
	}
	X <- cbind(Xfix,Xe,Xge%*%signs1)
	final.fit <- hglm(X=X,y=Y,Z=Xg,calc.like=T)
	X <- cbind(Xfix,Xe)
	null.fit <- hglm(X=X,y=Y,Z=Xg,calc.like=T)
	if(is.null(final.fit$likelihood) | is.null(null.fit$likelihood)) lrtre <- NA
	if(!is.null(final.fit$likelihood) & !is.null(null.fit$likelihood)) lrtre <- -2*(null.fit$likelihood$pbvh-final.fit$likelihood$pbvh) 
	list(LRT=lrtre,signs=signs1)
}

allocate3.reI2 <- cmpfun(al3.reI2)


#########################################################################
### FUNCTION FOR FINDING ALLOCATION OF GE INTERACTION using ridge penalty with R package "glmnet"
############Sequential allocation (-1,1 or 0) 
#########################################################################
al3.ridI2 <- function(Y,Xg,Xe,Xfix){
	n <- length(Y)
	Xge = NULL
	for (k in 1:length(Y)){
		Xge = rbind(Xge,kronecker(Xg[k,],Xe[k,]))
	}
	nPred <- ncol(Xge)
	signs1 <- signs2 <- signs3 <- rep(0,nPred)
	X <- cbind(Xfix,Xe)
	cvla <- cv.glmnet(cbind(X,Xg),Y,penalty.factor=c(rep(0,ncol(X)),rep(1,ncol(Xg))),alpha=0)
	lam <- cvla$lambda.min
	res1 <- deviance(glmnet(cbind(X,Xg),Y,penalty.factor=c(rep(0,ncol(X)),rep(1,ncol(Xg))),alpha=0,lambda=lam))
	###seq check whether to move each predictor to another group:
	for (k in 1:nPred){
		#print(k)
		signs2 <- signs3 <- signs1   #reset temp signs
		signs2[k] <- -1    #switch kth predictor to negative group
		signs3[k] <- 1    #switch kth predictor to negative group
		X <- cbind(Xfix,Xe,Xge%*%signs2)
		res2 <- deviance(glmnet(cbind(X,Xg),Y,penalty.factor=c(rep(0,ncol(X)),rep(1,ncol(Xg))),alpha=0,lambda=lam))
		X <- cbind(Xfix,Xe,Xge%*%signs3)
		res3 <- deviance(glmnet(cbind(X,Xg),Y,penalty.factor=c(rep(0,ncol(X)),rep(1,ncol(Xg))),alpha=0,lambda=lam))
		#cat(signs1,res1,res2,res3,"\n")
 		if (res2 <= res1 && res2 <= res3) {
    		signs1 <- signs2
    		res1 <- res2
    	}
    	if (res3 <= res1 && res3 <= res2) {
    		signs1 <- signs3
    		res1 <- res3
    	}
	}
	X <- cbind(Xfix,Xe,Xge%*%signs1)
	final.fit <- glmnet(cbind(X,Xg),Y,penalty.factor=c(rep(0,ncol(X)),rep(1,ncol(Xg))),alpha=0,lambda=lam)
	X <- cbind(Xfix,Xe)
	null.fit <- glmnet(cbind(X,Xg),Y,penalty.factor=c(rep(0,ncol(X)),rep(1,ncol(Xg))),alpha=0,lambda=lam) 
	list(LRT=deviance(null.fit)-deviance(final.fit),signs=signs1)
}

allocate3.ridgeI2 <- cmpfun(al3.ridI2)


#########################################################################
### FUNCTION FOR FINDING ALLOCATION OF GE INTERACTION using ridge penalty with R package "glmnet"
############Sequential allocation (-1,1 or 0) with separate betas for -1 and 1 
#########################################################################
al3.ridI2b <- function(Y,Xg,Xe,Xfix){
	n <- length(Y)
	Xge = NULL
	for (k in 1:length(Y)){
		Xge = rbind(Xge,kronecker(Xg[k,],Xe[k,]))
	}
	nPred <- ncol(Xge)
	signs1 <- signs2 <- signs3 <- rep(0,nPred)
	X <- cbind(Xfix,Xe)
	cvla <- cv.glmnet(cbind(X,Xg),Y,penalty.factor=c(rep(0,ncol(X)),rep(1,ncol(Xg))),alpha=0)
	lam <- cvla$lambda.min
	res1 <- deviance(glmnet(cbind(X,Xg),Y,penalty.factor=c(rep(0,ncol(X)),rep(1,ncol(Xg))),alpha=0,lambda=lam))
	###seq check whether to move each predictor to another group:
	for (k in 1:nPred){
		#print(k)
		signs2 <- signs3 <- signs1   #reset temp signs
		signs2[k] <- -1    #switch kth predictor to negative group
		signs3[k] <- 1    #switch kth predictor to negative group
		if (sum(signs2 == 1) == 0)
			X <- cbind(Xfix,Xe,Xge%*%(signs2 == -1))
		if (sum(signs2 == 1) > 0)
			X <- cbind(Xfix,Xe,Xge%*%(signs2 == -1),Xge%*%(signs2 == 1))
		res2 <- deviance(glmnet(cbind(X,Xg),Y,penalty.factor=c(rep(0,ncol(X)),rep(1,ncol(Xg))),alpha=0,lambda=lam))
		if (sum(signs3 == -1) == 0)
			X <- cbind(Xfix,Xe,Xge%*%(signs3 == 1))
		if (sum(signs3 == -1) > 0)
			X <- cbind(Xfix,Xe,Xge%*%(signs3 == -1),Xge%*%(signs3 == 1))
		res3 <- deviance(glmnet(cbind(X,Xg),Y,penalty.factor=c(rep(0,ncol(X)),rep(1,ncol(Xg))),alpha=0,lambda=lam))
		#cat(signs1,res1,res2,res3,"\n")
 		if (res2 <= res1 && res2 <= res3) {
    		signs1 <- signs2
    		res1 <- res2
    	}
    	if (res3 <= res1 && res3 <= res2) {
    		signs1 <- signs3
    		res1 <- res3
    	}
	}
	if (sum(signs1 == 1) == 0 & sum(signs1 == -1) > 0)
			X <- cbind(Xfix,Xe,Xge%*%(signs1 == -1))
	if (sum(signs1 == 1) > 0 & sum(signs1 == -1) == 0)
			X <- cbind(Xfix,Xe,Xge%*%(signs1 == 1))
	if (sum(signs1 == 1) > 0 & sum(signs1 == -1) > 0)
			X <- cbind(Xfix,Xe,Xge%*%(signs1 == -1),Xge%*%(signs1 == 1))
	final.fit <- glmnet(cbind(X,Xg),Y,penalty.factor=c(rep(0,ncol(X)),rep(1,ncol(Xg))),alpha=0,lambda=lam)
	X <- cbind(Xfix,Xe)
	null.fit <- glmnet(cbind(X,Xg),Y,penalty.factor=c(rep(0,ncol(X)),rep(1,ncol(Xg))),alpha=0,lambda=lam) 
	list(LRT=deviance(null.fit)-deviance(final.fit),signs=signs1)
}

allocate3.ridgeI2b <- cmpfun(al3.ridI2b)

#########LIKELIHOOD RATIO TEST
lrt <- function (obj1, obj2) {
    L0 <- logLik(obj1)
    L1 <- logLik(obj2)
    L01 <- as.vector(- 2 * (L0 - L1))
    df <- attr(L1, "df") - attr(L0, "df")
    return(c(pchisq(L01, df, lower.tail = FALSE),L01))
}
 
