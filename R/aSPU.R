###SumSqU and SumSqUw tests for interactions between GENE AND ENVIRONMENT
### BC, 4/6/2015


### Copyright: Wei Pan, 12/17/08
###            weip@biostat.umn.edu
###            http://www.biostat.umn.edu/~weip/

#input: sdat: a data-frame with components named
#              1) trait: a vector of 0's and 1's for the disease indicator;
#              2) Ro: genotypes for gene 1
#              3) Rt: genotypes for gene 2
#       B: # of simulations used to estimate the p-value for UminP; default=1000
# output: pvalues for 
#         1) main-effects + (2-way) interactions model: SumSqU, SumSqUw, score, (simulation-based) UminP, (Bonferonni-adjusted) UminP
#         2) main-effects model: SumSqU, SumSqUw, score, (simulation-based) UminP, (Bonferonni-adjusted) UminP
library(mvtnorm)
library(nlme)

PowerUniv <- function(U,V){
      n <- dim(V)[1]
      x <- as.numeric(max(abs(U)))
      TER <- as.numeric(1-pmvnorm(lower=c(rep(-x,n)),upper=c(rep(x,n)),mean=c(rep(0,n)),sigma=V))
      
      return(TER)
}

#######MODIFIED FOR QUANTITATIVE TRAIT AND GENE-ENVIRONMENT INTERACTION#######

aSPU_ge<-function(Y,Xg,Xe,Xfix,B=1000){
	n <- length(Y)
	p <- ncol(Xg)
	q <- ncol(Xe)
	##INTERACTIONS
	Xge = NULL
	for (k in 1:n){Xge = rbind(Xge,kronecker(Xg[k,],Xe[k,]))}
		
	tdat1<-data.frame(trait=Y,Xg,Xe,Xfix)
	fit1<-glm(trait~.,data=tdat1)
	pis<-fitted.values(fit1)
	Us <- matrix(0, nrow=n, ncol=p*q)
	for(i in 1:(p*q)){
	        tdat2<-data.frame(X1=Xge[,i],Xg,Xe,Xfix)
       		fit2<-glm(X1~.,data=tdat2)
      		X1mus <- fitted.values(fit2)
       		Us[, i]<-(Y-pis)*(Xge[,i] - X1mus)
    	}
    
    	U <- apply(Us,2,sum)
    	CovS <- matrix(0,nrow=p*q,ncol=p*q)
    	for (i in 1:n){CovS <- CovS + Us[i,]%*%t(Us[i,])}
     	
	V <- CovS

        ## SPU test statistics with gamma = 1,...,8, and infinity
        res <- c(abs(sapply(1:8,function(x) sum(U^x))),max(abs(U)))

        ## simulate B null score vectors
        permmat <- mvrnorm(B,rep(0,length(U)),V)

        ## apply SPU tests to null score vectors
        permstat <- cbind(abs(sapply(1:8,function(y) apply(permmat,1,function(x) sum(x^y)))), apply(permmat,1, function(z) max(abs(z))))

        p_res <- sapply(1:9, function(x) (sum(res[x] < permstat[,x])+1)/(B+1))
        aspu <- min(p_res)

        nullpval <- t(sapply(1:B, function(y) sapply(1:9, function(x) (sum(permstat[y,x] < permstat[,x])+1)/(B+1))))
        null_aspu <- apply(nullpval,1,min)

        p_aspu <- (sum(aspu > null_aspu)+1)/(B+1)
	
	return( list( stat = c(res,aspu) , pval = c(p_res,p_aspu) ))
}

### Exact ridge SPU test
ridgeSPU <- function(Y,Xg,Xe,Xfix,B=1000){
	## Has to be assigned globally because of the stupid groupedData function ##
	n <- length(Y)
	assign("dummyId",factor(rep(1,n)),envir=globalenv())
	assign("X",cbind(Xfix,Xe),envir=globalenv())
	dat <- groupedData(Y~1+X|dummyId)

	Z.block <- list(dummyId=pdIdent(~0+Xg))
	
	## run model and obtain variance terms
	res1 <- try(lme(Y~1+X,random=Z.block,data=dat,method="ML"))
	if (inherits(res1,"try-error"))
		{res1 <- lme(Y~1+X,random=Z.block,method="ML",data=dat,control=list(opt="optim"))}
	vars <- VarCorr(res1)
	sig_g <- as.numeric(vars[1])
	sig_e <- res1$sigma^2

	## build estimated sigma under H_0
	Sigma <- sig_g*Xg%*%t(Xg)+sig_e*diag(1,n)
	sinv <- solve(Sigma)

	tilX <- cbind(rep(1,n),X)
	alphahat <- solve(t(tilX)%*%sinv%*%tilX)%*%t(tilX)%*%sinv%*%Y

	Xge = NULL
	Xge <- t(sapply(1:n,function(x) rbind(Xge,kronecker(Xg[x,],Xe[x,])))) 

	p <- ncol(Xg)
	q <- ncol(Xe)
	Us <- Xgeb <- matrix(0, nrow=n, ncol=p*q)
	for(i in 1:(p*q)){
		X1 <- Xge[,i]
		assign("X1",X1,envir=globalenv())
	        dat <- groupedData(X1~1+X|dummyId)
       		fit2 <- try(lme(X1~1+X,random=Z.block,data=dat,method="ML"))
      		if (inherits(fit2,"try-error"))
                	{fit2 <- lme(X1~1+X,random=Z.block,data=dat,method="ML",control=list(opt="optim"))}
		X1mus <- fitted.values(fit2)
       		Xgeb[,i] <- (Xge[,i] - X1mus)
		Us[,i] <- Xgeb[,i]*sinv%*%(Y-tilX%*%alphahat)
    	}

	U <- colSums(Us)
	
	CovS <- matrix(0,nrow=p*q,ncol=p*q)
    	for (i in 1:n){CovS <- CovS + Us[i,]%*%t(Us[i,])}
	V <- CovS

	## SPU test statistics with gamma = 1,...,8, and infinity
	res <- c(abs(sapply(1:8,function(x) sum(U^x))),max(abs(U)))

	## simulate B null score vectors
	permmat <- mvrnorm(B,rep(0,length(U)),V)

	## apply SPU tests to null score vectors
	permstat <- cbind(abs(sapply(1:8,function(y) apply(permmat,1,function(x) sum(x^y)))), apply(permmat,1, function(z) max(abs(z))))

	p_res <- sapply(1:9, function(x) (sum(res[x] < permstat[,x])+1)/(B+1))
	aspu <- min(p_res)
	
	nullpval <- t(sapply(1:B, function(y) sapply(1:9, function(x) (sum(permstat[y,x] < permstat[,x])+1)/(B+1))))
	null_aspu <- apply(nullpval,1,min)

	p_aspu <- (sum(aspu > null_aspu)+1)/(B+1)	

	return( list( stat = c(res,aspu) , pval = c(p_res,p_aspu) ))

}
