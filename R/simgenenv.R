#####USING HAPMAP DATA####
## gene is the gene from HapMap data (there should be a text file with that data, e.g. gene="ADH1B")
## path is where to find that data (e.g. path = "/Users/brandoncoombes/Documents/Grad School/PhD Work/GE/HapMap_Genes/")
## mz is indicator of whether we are simulating a MZ twin (this skips simulating a whole family)
## n.person is the total number of each family member
simgen<-function(gene,path,mz=FALSE,n.person=100){
	hapmap <- read.csv(paste(path,gene,".txt",sep=""))
	## make the string into a genotype matrix
	geno <- do.call(cbind,strsplit(as.character(hapmap$CEU.genotype)," "))
	hapX <- matrix(NA,nrow=nrow(geno),ncol=ncol(geno))
	for (i in 1:ncol(geno)){
		## impute missing with MAF
		miss <- ifelse(grepl("N",geno[,i]),sum(rbinom(2,1,hapmap$reference.allele.frequency[i])),0)
		## Make the Minor Allele Count matrix
		hapX[,i] <- 2 - miss - grepl(hapmap$reference.allele[i],substr(geno[,i],1,1))-grepl(hapmap$reference.allele[i],substr(geno[,i],2,2))
	}
	dimnames(hapX)[[2]]<-hapmap$marker.id
    nhap <- nrow(hapX)
    k <- ncol(hapX)
    
## Generate genotype of Mother
	mo <- hapX[sample(1:nhap,n.person,replace=T),]

	switch <- which(colSums(mo)/(2*nrow(mo))>0.5)
	if (length(switch) > 0) 
		mo[,switch] <- sapply(switch,function(x) 2-mo[,x]) ##reverse minor allele so MAF < 0.5	

	## we can skip the below if MZ because they are identical to their twin    
    if (mz == TRUE){
    		xmat<-matrix(NA, nrow=2*n.person,ncol=k)
    	xmat[seq(1,2*n.person,2),]=mo
    	xmat[seq(2,2*n.person,2),]=mo 

    	dimnames(xmat)[[2]]<-hapmap$marker.id
    		return(list(MAF=1-hapmap$reference.allele.frequency,X=xmat))
    	}
    
## Generate genotype of Father
	fa <- hapX[sample(1:nhap,n.person,replace=T),]

	switch <- which(colSums(fa)/(2*nrow(fa))>0.5)
        if (length(switch) > 0)	
		fa[,switch] <- sapply(switch,function(x) 2-fa[,x]) ##reverse minor allele so MAF < 0.5
	
ofs1 <- ofs2 <- matrix(NA,n.person,k)
## offspring 1    
    for (snp in 1:k){      
      ofs1[which(mo[,snp]==0 & fa[,snp]==0),snp]=0
      ofs1[which((fa[,snp]==0 & mo[,snp]==2) | (fa[,snp]==2 & mo[,snp]==0)) ,snp]=1
      ofs1[which(fa[,snp]==2 & mo[,snp]==2),snp]=2
      ofs1[which((fa[,snp]==0 & mo[,snp]==1) | (fa[,snp]==1 & mo[,snp]==0)),snp] = sample(c(0,1),size=length(which((fa[,snp]==0 & mo[,snp]==1) | (fa[,snp]==1 & mo[,snp]==0))),prob=c(0.5,0.5),replace=T)
      ofs1[which((fa[,snp]==1 & mo[,snp]==2) | (fa[,snp]==2 & mo[,snp]==1)),snp] = sample(c(1,2),size=length(which((fa[,snp]==1 & mo[,snp]==2) | (fa[,snp]==2 & mo[,snp]==1))),prob=c(0.5,0.5),replace=T)
      ofs1[which(fa[,snp]==1 & mo[,snp]==1),snp] = sample(c(0,1,2),size=length(which(fa[,snp]==1 & mo[,snp]==1)),prob=c(0.25,0.5,0.25),replace=T)
    }
## offspring 2
    for (snp in 1:k){      
      ofs2[which(mo[,snp]==0 & fa[,snp]==0),snp]=0
      ofs2[which((fa[,snp]==0 & mo[,snp]==2) | (fa[,snp]==2 & mo[,snp]==0)) ,snp]=1
      ofs2[which(fa[,snp]==2 & mo[,snp]==2),snp]=2
      ofs2[which((fa[,snp]==0 & mo[,snp]==1) | (fa[,snp]==1 & mo[,snp]==0)),snp] = sample(c(0,1),size=length(which((fa[,snp]==0 & mo[,snp]==1) | (fa[,snp]==1 & mo[,snp]==0))),prob=c(0.5,0.5),replace=T)
      ofs2[which((fa[,snp]==1 & mo[,snp]==2) | (fa[,snp]==2 & mo[,snp]==1)),snp] = sample(c(1,2),size=length(which((fa[,snp]==1 & mo[,snp]==2) | (fa[,snp]==2 & mo[,snp]==1))),prob=c(0.5,0.5),replace=T)
      ofs2[which(fa[,snp]==1 & mo[,snp]==1),snp] = sample(c(0,1,2),size=length(which(fa[,snp]==1 & mo[,snp]==1)),prob=c(0.25,0.5,0.25),replace=T)
    }
## putting snps together
    xmat<-matrix(NA, nrow=2*n.person,ncol=k)
    xmat[seq(1,2*n.person,2),]=ofs1
    xmat[seq(2,2*n.person,2),]=ofs2

    dimnames(xmat)[[2]]<-hapmap$marker.id  
    return(list(MAF=1-hapmap$reference.allele.frequency,X=xmat))
}


## rho is correlation for the environments
## sym is indicates which correlation matrix is used (1 = symmetric, 0 = AR1)
## r is vector of within family correlation for each environment 
## n.person is the total number of each family member

simenv<-function(rho = 0.5,
			 	 sym = 1,
			 	 r = c(0,0,0,0),
                 n.person=100
                 ){
  library(mvtnorm)
  	n.env <- length(r)
    if (sym == 1){
    		Ve <- matrix(rho,nrow=n.env,ncol=n.env)
    		diag(Ve) <- 1
    } 
    else {
    		Ve <- diag(n.env)
    		for (i in 1:n.env){for (j in 1:n.env) Ve[i,j] <- rho^abs(i-j)}
    }
	
## offspring 1    
    ofs1 <- rmvnorm(n.person,rep(0,n.env),Ve)
## offspring 2
	ofs2 <- matrix(NA,n.person,n.env)
	for (i in 1:n.env) ofs2[,i] <- rnorm(n.person,ofs1[,i],1-r[i])

## putting snps together
    xmat<-matrix(NA, nrow=2*n.person,ncol=n.env)
    xmat[seq(1,2*n.person,2),]=ofs1
    xmat[seq(2,2*n.person,2),]=ofs2  
    return(xmat)
}


### This is a situation where genes and environments interact
#############################################################################################################
simphe_ge<-function(gendat=gen,		# gene matrix
	        envdat=env,				# environment matrix
		causal,					# causal snps marker.ids; causal SNPs are NOT randomly selected; it is the snps we choose
      	        h2=0,						# h2 is vector of heritability for gene betas (same length as causal) 
	        MAF=0.3,					# MAF is vector of the MAF (same length as causal)
                b0=0,						# intercept for the model
                s2e=5,						# s2e is the random error variance
                s2g=5,						# s2g is the family error variance
                betae=c(0,0,0,0),			# betae is a vector of enviromental betas (same length as number of environments)
                bge=0,
                nfam=100,					# number of families
                type="DZ"){				# type of siblings (MZ, DZ , or adopted)
   library(mvtnorm)
   library(Matrix)
   if(type=="DZ") phi<-matrix(c(1,0.5,0.5,1),nrow=2,byrow=T) else
   if(type=="MZ") phi<-matrix(c(1,1,1,1),nrow=2,byrow=T) else
   if (type=="ADOPT") phi<-diag(2)
 
   q <- length(betae)  ## number of environments
   
### first generate Random effect term, which is polygenic effect + environment effect
### construction of covariance matrix
	Vmat<-s2g*phi-s2g*diag(2)+s2e*diag(2)
 
 ##MAKE ALL INTERACTIONS Matrix###
	Xge <- NULL
	for (k in 1:(2*nfam)){
		Xge <- rbind(Xge,kronecker(gendat[k,],envdat[k,]))
	}
   
   ## simulate correlated errors between twins
   e<-rmvnorm(nfam,rep(0,2),Vmat)
   RE<-array(t(e[,1:2]))
     
## additive model
   h2<-as.matrix(h2)
   MAF<-as.matrix(MAF)
   b<-c(sqrt(h2/(2*MAF*(1-MAF))))
   nsnp<-dim(gendat)[2]
   nenv<-dim(envdat)[2]
   betage <- rep(0,nsnp*nenv)
   betag<-rep(0,nsnp)
   b0.index <- which(dimnames(gendat)[[2]] %in% causal)
   betag[b0.index]<-b
   betag<-as.matrix(betag) 

   for (k in 1:length(causal)){betage[(nenv*(b0.index[k]-1)+1):(nenv*(b0.index[k]-1)+nenv)] <- bge[(nenv*(k-1)+1):(nenv*(k-1)+nenv)]}   
   mu=as.matrix(gendat)%*%betag+as.matrix(envdat)%*%betae+as.matrix(Xge)%*%betage+b0
   Y=as.matrix(RE)+mu 
   return(list(Y=Y,truege=betage))     
}

### This is a situation where the environment regulates the gene to create an interaction
######################################################################################
simphe_envreg<-function(gendat=gen,		# gene matrix
			     envdat=env,				# environment matrix
			     causal,					# causal snps marker.ids; causal SNPs are NOT randomly selected; it is the snps we choose
			     h2=0,						# h2 is vector of heritability for gene betas (same length as causal) 
			     MAF=0.3,					# MAF is vector of the MAF (same length as causal)
                b0=0,						# intercept for the model
                s2e=5,						# s2e is the random error variance
                s2g=5,						# s2g is the family error variance
                betae=c(0,0,0,0),			# betae is a vector of enviromental betas (same length as number of environments)
                nfam=100,					# number of families
                type="DZ"){				# type of siblings (MZ, DZ , or adopted)
   library(Matrix)
   if(type=="DZ") phi<-matrix(c(1,0.5,0.5,1),nrow=2,byrow=T) else
   if(type=="MZ") phi<-matrix(c(1,1,1,1),nrow=2,byrow=T) else
   if (type=="ADOPT") phi<-diag(2)
 
   q <- length(betae)  ## number of environments
   
### first generate Random effect term, which is polygenic effect + environment effect
### construction of covariance matrix
	Vmat<-s2g*phi
   ## simulate correlated errors between twins
   e<-rmvnorm(nfam,rep(0,2),Vmat)
   RE<-array(t(e[,1:2]))
     
## additive model
   h2<-as.matrix(h2)
   MAF<-as.matrix(MAF)
   b<-c(sqrt(h2/(2*MAF*(1-MAF))))
   nsnp<-dim(gendat)[2]
   nenv<-dim(envdat)[2]
   b0.index <- which(dimnames(gendat)[[2]] %in% causal)
   Xg <- gendat[,b0.index]
   reg <- rowSums(envdat > 2)*2+rowSums(envdat > 1)+rowSums(envdat > 0)*0.5+rowSums(envdat > -1)*0.25
   mu=(Xg*bg)%*%b+as.matrix(envdat)%*%betae+b0
   Y=as.matrix(RE)+mu  
   return(list(Y=Y))     
}
   
  

	



