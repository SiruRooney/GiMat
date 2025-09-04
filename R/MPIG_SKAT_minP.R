#'@title Resampling p-values of associations function.
#'@description Sample p-values from different Gimat with different kernels
#'@param n.pheno The number of phenotypes
#'@param n The number of individuals
#'@param res_model The residuals from null models.
#'@param obj.list A list of outputs of MPiG with different kernels.
#'@param Z1 The gene-by-gene or gene-by-env interactions
#'@param X Covariates.
#'@param resample The sample sizes.
#'@import SKAT
#'@import stats
#'@return p-values.

MPiG_resampling_pvalues<-function(n.pheno, n,res_model,obj.list, Z1, X, resample = 200)
{
  n.obj<-length(obj.list)

  m<-NCOL(Z1)
  #R<-R1<-A<-Lambda<-wt<-vector("list", n.obj)
  R<-R1<-vector("list", n.obj)

  for (i in 1:n.obj) {
    R[[i]]<-obj.list[[i]]$R
    R1[[i]]<-mat.sqrt(R[[i]])
  }

  #V.item=cov(as.matrix(obj.res$res_model))
  V.item=cov(as.matrix(res_model))
  V.item.inv=Get_GenInverse(V.item)
  V.half<-mat.sqrt(V.item.inv)
  ############################################################
  Qs<-matrix(0,ncol=n.obj,nrow=resample)
  for (j in 1:resample){
    res.1<-scale(as.matrix(rnorm(n*n.pheno),ncol=n.pheno),center=T, scale=F)

    res.1<-as.vector(t(res.1))
    res.Vhalf<-rep(0, n*n.pheno)
    for (id1 in 1:n){
      idx<-((id1-1)*n.pheno+1):(id1*n.pheno)
      res.Vhalf[idx]<-V.half%*%res.1[idx]
    }

    res.mat.2=scale(matrix(res.Vhalf,byrow=TRUE,ncol=n.pheno),center=T,scale=F)
    Q.temp1<-t(res.mat.2)%*%Z1
    trs<-vector("list",n.obj)
    for(indx in 1:n.obj){
      trs[[indx]]<-matrix(0, ncol=m, nrow=n.pheno)
    }
    # tr1 <- tr2 <- tr3 <- tr4 <- matrix(0, ncol = m, nrow = n.pheno)
    for (indx in 1:n.obj){
      for (i in 1:m){
        trs[[indx]][, i]<-R1[[indx]]%*%Q.temp1[,i]
      }
    }
    for (indx in 1:n.obj) {
      Qs[j, indx]<-sum(trs[[indx]]^2)
    }
    #  m1<-m*n.pheno
  }

  pv<-matrix(0,ncol=n.obj,nrow=resample)
  for (i in 1:n.obj) pv[[i]]<-array()
  Lambda<-vector("list",n.obj)


  for (j in 1:n.obj) {
    Lambda[[j]]<- Get_Lambda(obj.list[[j]]$eig_chi)
    for (i in 1:resample)pv[i,j]<-SKAT:::Get_PValue.Lambda(Lambda[[j]],Qs[i, j])$p.value
  }
  return(list(p.values = pv, Sigma_Ps = R))
}



#'@title min P-values of all kernels
#'@param n.pheno The number of phenotypes
#'@param n The number of individuals
#'@param res_model The residuals from null models.
#'@param obj.list A list of outputs of MPiG with different kernels.
#'@param Z1 The gene-by-gene or gene-by-env interactions
#'@param X Covariates.
#'@param resample The resampling times.
#'@return P-values among all p-values from different kernels.
#'@import copula
#'@export
#'@import stats
#'@import SKAT
MPiG_minP_base<-function(n.pheno,n,res_model,obj.list,Z1,X,resample=200){

  Q <-MPiG_resampling_pvalues(n.pheno,n,res_model,obj.list, Z1, X, resample = 200)$p.values
  n.obj <- length(obj.list)

  pv <- array()
  for (i in 1:n.obj) pv[i] <- obj.list[[i]]$p.value
  minP <- min(pv)
  p.minp = n.obj * minP
  ind <-adj.Lambda(eigen(cov(Q))$values)
  R.mat <- cor(Q, method = "kendall")
  obj.t.Copula =tCopula(R.mat[lower.tri(R.mat)],dim=nrow(R.mat),dispstr="un")
  obj.norm.Copula = tCopula(R.mat[lower.tri(R.mat)],dim=nrow(R.mat),dispstr="un")
  p.minp <- 1 - pCopula(1 - minP, obj.t.Copula)

  if(minP==0){
    p.minp=minP
  }else{

  if ((p.minp/minP > n.obj) || (p.minp/minP < 1)){
    p.minp <- 1 - pCopula(1 - minP, obj.norm.Copula)
  }
  if ((p.minp/minP > n.obj) || (p.minp/minP < 1))
    p.minp <- n.obj * minP
  }
  return(list(p.value = p.minp))

}

adj.Lambda<-function (lambda1)
{
    IDX1 <- which(lambda1 >= 0)
    IDX2 <- which(lambda1 > mean(lambda1[IDX1])/1e+05)
    if (length(IDX2) == 0) {
        return(0)
    }
    else {
        lambda <- lambda1[IDX2]
        return(lambda)
    }
}

