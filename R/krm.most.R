krm.most = function(formula, data, regression.type=c("logistic","linear"), 
    kern.type=c("rbf","mi","mm","prop"), n.rho=10, range.rho=0.99, 
    n.mc=2000, 
    seq.file.name=NULL, formula.kern=NULL,
    inference.method=c("parametric.bootstrap", "perturbation", "LGL2008"),
    verbose=FALSE) 
{
    
    # match input parameters
    if (length(regression.type)!=1) stop("regression.type has to be specified")
    regression.type <- match.arg(regression.type)    
    kern.type <- match.arg(kern.type)    
    inference.method <- match.arg(inference.method)
    if (inference.method=="LGL2008") {
        if (kern.type=="mi") {
            cat("Profile HMM MI kernel is currently not implemented with LGL2008\n")
            return (NA) # the way rho is defined for LGL2008 does not apply to MI kernel
        }
        n.rho=500 # according to LGL2008 paper
    }
    
    # check input parameters
    if (!is.null(seq.file.name) & !kern.type %in% c("mi","mm","prop")) stop("choose a kernel from mi, mm, prop for sequence kernels")
    if (!is.null(formula.kern) & !kern.type %in% c("rbf")) stop("choose a kernel from rbf for Euclidean kernels")
    if (is.null(formula.kern) & is.null(seq.file.name)) stop("either formula.kern or seq.file.name has to be supplied")
    if (!is.null(seq.file.name)) {stopifnot (is.character(seq.file.name)); if(verbose) myprint(seq.file.name) }
    if (!is.null(formula.kern)) {stopifnot (is(formula.kern,"formula")); if(verbose) {cat("formula.kern: "); print(formula.kern)}}
    
    # save rng state before set.seed in order to restore before exiting this function
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {        
        set.seed(1)
        save.seed <- get(".Random.seed", .GlobalEnv)
    }                        
    
    # define rhos
    if (!is.null(seq.file.name)) {
        l2dist = -log(getSeqKernel(seq.file.name, kern.type, tau=1, call.C=T))
    } else if (!is.null(formula.kern)) {
        Z=model.matrix(formula.kern, data)
        if (colnames(Z)[1]=="(Intercept)") Z=Z[,-1,drop=FALSE]
        l2dist = -log(getK(Z,kernel=kern.type,para=1))        
    }
    diag(l2dist)=NA
    if (inference.method!="LGL2008") {
        median.log.K = median(l2dist, na.rm=T)
        if (median.log.K==0) {print("l2dist has mean and median 0"); return (rep(NA,4))}    
        rho.alpha=(1-range.rho)/2
        rhos=sqrt(-median.log.K/log(seq(rho.alpha,1-rho.alpha, length=n.rho))) # try to space rho uniformly on the correlation scale
        
    } else {
        # implementing Liu et al 2008
        maxk = max(l2dist, na.rm=T)
        mink = min(l2dist, na.rm=T)
        if (mink==0) {
            l2dist[l2dist==0]=NA
            mink = min(l2dist, na.rm=T)
        }
        rhos=sqrt(seq(mink/10, maxk*10, length=n.rho))
    }
    if (verbose) {myprint(summary(c(l2dist))); myprint(summary(c(rhos)))}
    
    # terms used in parametric boostrap and perturbation
    X=model.matrix(formula, data)
    y=model.frame(formula, data)[[1]]
    if (verbose==2) {cat("str(y):\n"); str(y)}
    n=nrow(data)
    if (regression.type=="logistic") {    
        fit=glm(formula, data, family="binomial") # null model fit
        # terms derived from the model and don't depend on kernel, many are not needed for parametric bootstrap, but computed anyway
        beta.h=coef(fit)
        mu.h=drop(expit(X %*% beta.h))
        D.h = c(mu.h*(1-mu.h)) 
        XD.5 <- X * D.h^.5
        . <- eigen(crossprod(XD.5))
        V.beta.h <- crossprod(t(.$vectors) * .$values^(-0.5))  # solve(t(X) %*% D.h %*% X)
        V.eta.h  <- crossprod(t(X %*% .$vectors) * .$values^(-0.5)) # X %*% V.beta.h %*% t(X)
        V.mu.h   <- DXD(D.h,V.eta.h,D.h) # D.h %*% V.eta.h %*% (D.h)
        P.h1=diag(D.h) - V.mu.h # this is the variance of Y - mu^hat to be used in simulation mvn
        extra.kurtosis = (mu.h*(1-mu.h)^4 + (1-mu.h)*mu.h^4) - 3 * D.h**2 # needed for mean estimate
        A.h=diag(n) - DXD(D.h,V.eta.h,rep(1,n))# D.h %*% V.eta.h, needed for variance estimate        
            
    } else if (regression.type=="linear") {    
        fit=glm(formula, data, family="gaussian") # null model fit
        # terms derived from the model and don't depend on kernel, many are not needed for parametric bootstrap, but computed anyway
        beta.h=coef(fit)
        mu.h=drop(X %*% beta.h)
        noise.sd = summary(fit)$dispersion ** .5
        P = diag(n) - X %*% solve(crossprod(X)) %*% t(X)        
    
    }     
    resid.=y-mu.h
    if(verbose) {myprint(dim(X)); myprint(beta.h)}
    
    # compute kernel for rho=1
    if (!is.null(seq.file.name)) {
        K.1 = getSeqKernel(seq.file.name, kern.type, tau = 1, call.C = T)
    } else if (!is.null(formula.kern)) {
        K.1 = getK(Z, kernel=kern.type, para=1) 
    }
    
    # do inference
    Q.rho.stats = sapply(rhos, simplify = "array", function(rho) {
        if (verbose==2) myprint(rho)
        
        # compute K_rho
        K = K.1 ^ (rho^-2)
    
        # parametric bootstrap-based inference
        if (inference.method=="parametric.bootstrap"){        
            test.stats.obs = krm.score.test (formula, data, K, regression.type) 
            
            test.stats.rep = sapply(1:n.mc, simplify = "array", function(j){
                set.seed(j+1e6)                         
                if (regression.type=="logistic") {
                    new.y <- rbern(n, mu.h)
                } else if (regression.type=="linear") {
                    new.y <- rnorm(n, mean=mu.h, sd=noise.sd) 
                }
                new.dat = data.frame(new.y, X)
                names(new.dat)[1] = as.character(formula)[2]
                krm.score.test (formula, new.dat, K, regression.type, verbose=verbose) 
            })
            cbind(test.stats.obs,test.stats.rep)        
        
        # perturbation-based inference 
        } else if (inference.method=="perturbation") {
            require(MASS)
            if (regression.type=="logistic") {            
                m1=tr(P.h1 %*% K)
                W.h=crossprod(A.h,symprod(K,A.h)) # t(A.h) %*% K %*% A.h
                v2 = varQ (W.h, variance=D.h, extra.kurtosis=extra.kurtosis, do.C=TRUE)
                test.stats.obs = perturbation.test(resid., K, m1, v2)  
                if (verbose) {
                    test.stats.obs.2 = krm.score.test (formula, data, K, regression.type, verbose=FALSE)           
                    stopifnot(all(test.stats.obs.2==test.stats.obs, na.rm=TRUE)) # sanity check
                }
    
                test.stats.rep = sapply(1:n.mc, simplify = "array", function(j){
                    set.seed(j+1e6) 
                    new.resid <- mvrnorm(1,mu=rep(0,n),Sigma=P.h1) 
                    perturbation.test(new.resid, K, m1, v2)                        
                })
            
            } else if (regression.type=="linear") {
                PK <- P %*% K            
                m=tr(PK)        
                V.Q.norm = 2*tr(crossprod(PK)) # 2*tr( P %*% K %*% P %*% K )
                test.stats.obs = perturbation.test(resid./noise.sd, K, m, V.Q.norm)  
                if (verbose) {
                    test.stats.obs.2 = krm.score.test (formula, data, K, regression.type)
                    stopifnot(all(test.stats.obs.2==test.stats.obs, na.rm=TRUE)) # sanity check
                }
                               
                test.stats.rep = sapply(1:n.mc, simplify = "array", function(j){
                    set.seed(j+1e6) 
                    new.resid <- mvrnorm(1,mu=rep(0,n),Sigma=P) 
                    perturbation.test(new.resid, K, m, V.Q.norm)                        
                })
                
            }
            cbind(test.stats.obs,test.stats.rep)        
            
        } else if (inference.method=="LGL2008") {
            if (regression.type=="logistic") {            
                Ph1K <- P.h1 %*% K
                m1=tr(Ph1K)
                V.Q.norm.h = 2*tr(crossprod(Ph1K)) # 2*tr( P.h1 %*% K %*% P.h1 %*% K )
                Q = txSy(resid.,K,resid.) # drop(t(new.resid) %*% K %*% new.resid)    
                (Q-m1)/sqrt(V.Q.norm.h)
                
            } else stop("LGL2008 is logistic regression only")
        
        } else stop ("something wrong in krm.most")
                
        
    }) 
    ## if inference.method is not LGL2008, Q.rho.stats is 3-dimensional array: # of test stats considered  by   1+n.mc   by   # of rhos
    ## if inference.method is LGL2008, Q.rho.stats is vector of length n.rho
    assign(".Random.seed", save.seed, .GlobalEnv) # restore rng state 
    if (verbose==2) {cat("str(Q.rho.stats):\n"); str(Q.rho.stats)}
        
    if (inference.method!="LGL2008") {
        # maximize over # of rhos
        p.chi.sup.1 <- drop(apply(Q.rho.stats["chiI",,,drop=FALSE],1:2,max))
        p.chi.sup.2 <- drop(apply(Q.rho.stats["chiII",,,drop=FALSE],1:2,max))
        Q.norm.sup.1 <- drop(apply(Q.rho.stats["normI",,,drop=FALSE],1:2,max))
        Q.norm.sup.2 <- drop(apply(Q.rho.stats["normII",,,drop=FALSE],1:2,max))    
        if (verbose==2) {cat("str(p.chi.sup.1):\n"); str(p.chi.sup.1)}    
            
        c(
            "chiI"=mean(p.chi.sup.1[-1] > p.chi.sup.1[1]), 
            "chiII"=mean(p.chi.sup.2[-1] > p.chi.sup.2[1]), 
            "normI"=mean(Q.norm.sup.1[-1] > Q.norm.sup.1[1]), 
            "normII"=mean(Q.norm.sup.2[-1] > Q.norm.sup.2[1])
        )
    } else {
        M = max(Q.rho.stats)
        W = sum(abs(diff(Q.rho.stats)))
        pnorm(-M) + W*exp(-0.5*(M^2))/sqrt(8*pi)
    }
}


#
perturbation.test = function (new.resid, K, m, v) {    
    Q = txSy(new.resid,K,new.resid) # drop(t(new.resid) %*% K %*% new.resid)    
    s=v/(2*m); k=2*m^2/v
    c("chiI"=pchisq(Q/s, df=k), "chiII"=NA, "normI"=pnorm((Q-m)/sqrt(v)), "normII"=NA)    
}


# keep getK here so that when this file is source, krm.most can run without having to source another file
# note this function exists in both svmw package and aucm package. be sure to change both if some changes are needed.
#calculates the kernel matrix between X and itself and returns a n by n matrix. Alternatively, it calculates the kernel matrix between X and X2 and returns a n by n2 matrix. 
#{X}{covariate matrix with dimension n by d. Note this is not the paired difference of covariate matrix.}
#{kernel}{string specifying type of kernel:
#    polynomial or p (1 + <x,y>)^para,
#    rbf or r exp(-para*||x-y||^2),
#    linear or l <x,y>,
#    no default.
#{para}{parameter of the kernel fucntion}
#{X2}{optional second covariate matrix with dimension n2 by d}
# return a kernel matrix
getK=function(X,kernel,para=NULL,X2=NULL){
    kernel=substr(kernel,1,1)
    if (kernel=="r" | kernel=="e") {
        if (!is.null(X2)) {
            aux = X[rep(1:nrow(X),nrow(X2)),,drop=F] - X2[rep(1:nrow(X2),each=nrow(X)),,drop=F]
            dist.mat = matrix(rowSums(aux^2), nrow=nrow(X))
#            aux=X2[rep(1:nrow(X2),nrow(X)),] - X[rep(1:nrow(X),each=nrow(X2)),]
#            dist.mat = matrix(rowSums(aux^2), nrow=nrow(X2))
        } else {
            dist.mat = as.matrix(dist(X))^2
        }
    }
    
    if (is.null(X2)) X2=X
    switch(kernel, 
        p=(tcrossprod(X,X2)+1)^para, # polynomial
        r=exp(-para*dist.mat), # rbf
        e=dist.mat, # Euclidean distance
        l=tcrossprod(X,X2), # linear
        stop(kernel %+% " kernel not supported")
    )
}
