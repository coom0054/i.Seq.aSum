setwd("/share/bulk/gedi/coom0054/GE/Independent/Output")

library(parallel)
library(iSKAT)
library(glmnet)
source("../simgenenv.R")
source("../seqSum_ind.R")
source("../aSPU.R")

hapmap.path = "/share/bulk/gedi/coom0054/GE/HapMap/"

n.mz = 1000

step = as.numeric(Sys.getenv("PBS_ARRAYID"))
situation = "main"
inttype = "equal"

b0 <- 0
sig_err <- 10
sig_fam <- 1 
gene <- "ADH1B"        # 10 SNPS
# causal <- c("rs1042026","rs1789882")
causal <- c("rs1042026","rs1789882","rs17033","rs2075633")   
# gene <- "ALDH2"          # 14 SNPS
# causal <- c("rs4767939","rs886205")  ## two different LD blocks
# gene <- "SLC6A4"          # 19 SNPS
# causal <- c("rs140701","rs6354")
# gene <- "ADH1C"          # 25 SNPS
# causal <- c("rs698","rs1693424")    ##rs698 from literature
# gene <- "DRD2"          # 74 SNPS
# causal <- c("rs4581480","rs4936274")
gene <- "GABRA2"          # 99 SNPS
causal <- c("rs534459","rs505474","rs540363","rs519270")  ##SNPs from two different LD blocks

if (situation == "nomain"){
        be <- rep(0.5,1)
        h2 <- rep(0,length(causal))
}
if (situation == "main"){
        be <- rep(0.5,1)
        h2 <- rep(0.2,length(causal))
}
if (inttype == "equal"){
        beta <- 0.05*step
        if (beta < 0){add="-"}
        else {add = ""}
        bge <- c(rep(beta,length(causal)/2),rep(beta,length(causal)/2))
}
if (inttype == "opp"){
        beta <- 0.05*step
        add <- ""
        bge <- c(rep(beta,length(causal)/2),rep(-beta,length(causal)/2))
}

if (beta == 0){nsim = 10000}
if (beta != 0){nsim = 1000}

title <- "block"

sfile = paste("stat",title,"_",situation,inttype,"_",add,beta,"_",gene,".txt",sep="")
pfile = paste("pval",title,"_",situation,inttype,"_",add,beta,"_",gene,".txt",sep="")
tfile = paste("time",title,"_",situation,inttype,"_",add,beta,"_",gene,".txt",sep="")

calc.stats <- function(Y,Xg,Xe,Xcov){
	ones <- rep(1,n)
	Xge <- NULL
	for (id in 1:n){
		Xge = rbind(Xge,kronecker(Xg[id,],Xe[id,]))
	}

	## min test 
        t <- proc.time()[1]
        single_int <- function(snp,env){
                G <- Xg[,snp]
                E <- Xe[,env]
                altlk <- lm(Y ~ Xcov + G + E + G*E)
		return(coef(summary(altlk))[4+ncol(Xcov),4])
        }
	res <- c(sapply(1:ncol(Xg), function(x) sapply(1:ncol(Xe), function(y) single_int(x,y))))
        min_test <- min(res)*length(res)
        time.min <- proc.time()[1]-t

	t <- proc.time()[1]
	#######FULL MODEL######
	alt <- glm(Y ~ Xcov + Xg + Xe + Xge)
	null <- glm(Y ~ Xcov + Xg + Xe)
	lrt <- lrt(null,alt)[2]
	p.lrt <- lrt(null,alt)[1]
	time.lrt <- proc.time()[1]-t
	
	t <-  proc.time()[1]
	######SPU TESTS#####
	SPU <- aSPU_ge(Y,Xg,Xe,Xcov,B=1000)
	time.spu <- proc.time()[1]-t

	t <-  proc.time()[1]
	######SPU TESTS#####
	rSPU <- ridgeSPU(Y,Xg,Xe,Xcov,B=1000)
	time.rspu <- proc.time()[1]-t
		
	t <- proc.time()[1]
	#########iSKAT#########
	if (length(unique(Xcov)) == 1){Xcov=NULL}
	iSKAT <- iSKAT(Z=Xg,Y=matrix(Y,nrow=n,ncol=1),E=Xe,X=Xcov,out_type="C")
	time.skat <- proc.time()[1]-t

	t <- proc.time()[1]
	#######SSI#######
	int2 <- allocate3.M(Y,cbind(Xcov,Xg,Xe),Xge)
	#print(int2$signs)
	alt <- glm(Y ~ Xcov + Xg + Xe + Xge%*%int2$signs)
	null <- glm(Y ~ Xcov + Xg + Xe)
	SSI <- lrt(null,alt)[2]
	p.ssi <- NA
	time.ssi <- proc.time()[1]-t

	t <- proc.time()[1]
	########SS.MARGINAL#########
	snp1 <- allocate3.M(Y,Xcov,Xg)
	env1 <- allocate3.M(Y,Xcov,Xe)
	alt <- glm(Y ~ Xcov + Xg%*%snp1$signs + Xe%*%env1$signs +  Xg%*%snp1$signs*Xe%*%env1$signs)
	null <- glm(Y ~ Xcov + Xg%*%snp1$signs + Xe%*%env1$signs)
	SSM <- lrt(null,alt)[2]
	p.ssm <- 1-pchisq(SSM,df=1)
	time.ssm <- proc.time()[1]-t
		
	t <- proc.time()[1]
	########SS.JOINT#########
	both1 <- allocate3.J(Y,Xcov,Xg,Xe)
	alt <- glm(Y ~ Xcov + Xg%*%both1$gsigns + Xe%*%both1$esigns +  Xg%*%both1$gsigns*Xe%*%both1$esigns)
	null <- glm(Y ~ Xcov + Xg%*%both1$gsigns + Xe%*%both1$esigns)
	SSJ <- lrt(null,alt)[2]
	p.ssj <- 1-pchisq(SSJ,df=1)
	time.ssj <- proc.time()[1]-t
		
	t <- proc.time()[1]
	########SS.GENE#########
	snp2 <- allocate3.M(Y,cbind(Xcov,Xe),Xg)
	gscore <- Xg%*%snp2$signs
	ge <- apply(Xe,2,function(x) gscore*x)
	alt <- glm(Y ~ Xcov + Xe + gscore + ge)
	null <- glm(Y ~ Xcov + Xe + gscore)
	SSG <- lrt(null,alt)[2]
	p.ssg <- 1-pchisq(SSG,df=ncol(Xe))
	time.ssg <- proc.time()[1]-t
					

	stats <- c(min_test,lrt,SPU$stat,rSPU$stat,iSKAT$pvalue,SSI,SSM,SSJ,SSG)
	times <- c(time.min,time.lrt,rep(time.spu,10),rep(time.rspu,10),time.skat,time.ssi,time.ssm,time.ssj,time.ssg)
	pvals <- c(min_test,p.lrt,SPU$pval,rSPU$pval,iSKAT$pvalue,p.ssi,p.ssm,p.ssj,p.ssg)
	tests <- c("min_test","LRT","Sum","SSU",paste("SPU",3:8,sep=""),"UminP","aSPU",
			"ridge.Sum","ridge.SSU",paste("ridge.SPU",3:8,sep=""),"ridge.UminP","ridge.aSPU","iSKAT","SSI","SSM","SSJ","SSG")
	
	return(list(stat=stats,time=times,p=pvals,method=tests))
}

Oneiteration<-function(i){ #one simulation to compare different methods
set.seed(i)
print(i)

Xg.mz <- simgen(gene,hapmap.path,n.person=n.mz)
maf <- Xg.mz$MAF[which(dimnames(Xg.mz$X)[[2]] %in% causal)]
Xe.mz <- as.matrix(simenv(n.person=n.mz)[,1])

Y.mz <- simphe_ge(Xg.mz$X,Xe.mz,causal=causal,h2=h2,MAF=maf,b0=b0,s2e=sig_err,s2g=sig_fam,
						betae=be,bge=bge,nfam=n.mz,type="MZ")

ind.id <- seq(from=1,to=n.mz*2,by=2)

Xg <- as.matrix(Xg.mz$X[ind.id,])
Xe <- as.matrix(Xe.mz[ind.id,])
Y <- Y.mz$Y[ind.id]

n <- length(Y)

male <- rbinom(n,1,0.48)
Xcov <- cbind(male)

Y <- rnorm(n,Y+male*0.5,sd=0.01)

assign("Xg",as.matrix(Xg),envir=globalenv())
assign("Xe",as.matrix(Xe),envir=globalenv())
assign("Xcov",as.matrix(Xcov),envir=globalenv())
assign("Xfix",as.matrix(Xcov),envir=globalenv())
assign("Y",Y,envir=globalenv())
assign("n",length(Y),envir=globalenv())

nolasso <- calc.stats(Y,Xg,Xe,Xcov)

##############################################################################################################
##### We can also shrink down the number of interactions by using Lasso to shrink Xg #########################
lasso <- cv.glmnet(cbind(Xcov,Xe,Xg),Y,penalty.factor=c(rep(0,ncol(Xcov)+ncol(Xe)),rep(1,ncol(Xg))),alpha=1)
lambda <- which(lasso$lambda==lasso$lambda.min)
fit <- lasso$glmnet.fit
keep <- as.numeric(coef(lasso,s="lambda.min"))[-(1:(ncol(Xcov)+ncol(Xe)+1))]>0
#print(as.numeric(keep))
Xgnew <- Xg[,keep]

while (sum(keep==FALSE) == ncol(Xg)){
                lambda <- lambda+1
                keep <- as.numeric(coef(fit)[,lambda])[-(1:(ncol(Xe)+2))]>0
                #print(keep)
}

lasso <- list(stat=rep(NA,12),A=NULL,time=rep(NA,12),p=rep(NA,12),method=NULL)
if (sum(keep==FALSE) < ncol(Xg)){
	Xgnew <- Xg[,keep]
	lasso <- calc.stats(Y,Xgnew,Xe,Xcov)
}

if (i == 1){
	write.table(t(c("ID",nolasso$method)),file=pfile,append=F,sep=",",quote=F,row.names=F,col.names=F)
	write.table(t(c("ID",nolasso$method)),file=sfile,append=F,sep=",",quote=F,row.names=F,col.names=F)
	write.table(t(c("ID",lasso$method)),file=paste("lasso",sfile,sep="_"),append=F,sep=",",quote=F,row.names=F,col.names=F)
	write.table(t(c("ID",lasso$method)),file=paste("lasso",pfile,sep="_"),append=F,sep=",",quote=F,row.names=F,col.names=F)
	write.table(t(c(i,nolasso$p)),file=pfile,append=T,sep=",",quote=F,row.names=F,col.names=F)
	write.table(t(c(i,nolasso$stat)),file=sfile,append=T,sep=",",quote=F,row.names=F,col.names=F)
	write.table(t(c(i,lasso$p)),file=paste("lasso",pfile,sep="_"),append=T,sep=",",quote=F,row.names=F,col.names=F)
	write.table(t(c(i,lasso$stat)),file=paste("lasso",sfile,sep="_"),append=T,sep=",",quote=F,row.names=F,col.names=F)
	if (beta == 0){
		write.table(t(c("ID",nolasso$method)),file=tfile,append=F,sep=",",quote=F,row.names=F,col.names=F)
		write.table(t(c(i,nolasso$time)),file=tfile,append=T,sep=",",quote=F,row.names=F,col.names=F)
	}
}
if (i > 1){
	write.table(t(c(i,nolasso$p)),file=pfile,append=T,sep=",",quote=F,row.names=F,col.names=F)
	write.table(t(c(i,nolasso$stat)),file=sfile,append=T,sep=",",quote=F,row.names=F,col.names=F)
	write.table(t(c(i,lasso$p)),file=paste("lasso",pfile,sep="_"),append=T,sep=",",quote=F,row.names=F,col.names=F)
	write.table(t(c(i,lasso$stat)),file=paste("lasso",sfile,sep="_"),append=T,sep=",",quote=F,row.names=F,col.names=F)
	if (beta == 0){
		write.table(t(c(i,nolasso$time)),file=tfile,append=T,sep=",",quote=F,row.names=F,col.names=F)
	}
}

}

if (!file.exists(sfile)){
        Oneiteration(1)
        output<-mclapply(2:nsim,FUN=Oneiteration, mc.cores=5)
}
if (file.exists(sfile)){
        res <- read.csv(sfile,stringsAsFactors=F)
        do.sim <- (1:nsim)[! 1:nsim %in% res$ID]
        output<-mclapply(do.sim,FUN=Oneiteration, mc.cores=20)
}

