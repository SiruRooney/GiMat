#'@title Create weights for rare variants.
#'@param maf Minor allele frequencies.
#'@param weights_beta Two parameters in Beta distribution: Beta(a1,a2)
#'@return Output the weights from Beta distribution.
#'@import stats
Beta_weights<-function(maf,weights_beta) {weights=dbeta(maf,weights_beta[1],weights_beta[2]);return(weights)}


#'@title Impute function
#'@description Impute the missing values for genotype.
#'@param Z Genotype matrix with N times p.
#'@param impute.method Imputed method consists of "random", "fixed" and "bestguess".
#'@return Output the complete genotype matrix without missing values.
#'@import stats

Impute<-function(Z, impute.method){

  p<-dim(Z)[2]

  if(impute.method =="random"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-rbinom(length(IDX),2,maf1)
      }
    }
  } else if(impute.method =="fixed"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-2 * maf1
      }
    }
  } else if(impute.method =="bestguess") {

    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-round(2 * maf1)
      }
    }

  } else {
    stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
  }

  return(Z)
}


#'@title Gene-by-gene interaction.
#'@description Create gene-by-gene interaction.
#'@param G1 The genotype matrix (N times p1) from the first phenotype.
#'@param G2 The genotype matrix (N times p2) from the second phenotype.
#'@param D The environmental factor. By defaul, D=NULL.
#'@param nrow_G The number of individuals.
#'@param iG_with The variable used for indicating gene-by-gene interactions or gene-by-environment interactions. iG_with="G" refers to gene-by-gene interactions and iG_with="D" refers to gene-by-env interactions.
#'@param impute.method Imputed method consists of "random", "fixed" and "bestguess".
#'@param missing_cutoff The threshold for the missing rate. If the missing ratio is less than the threshold, impute method is utilized.
#'@param max_MAF The threshold of maximum MAF.
#'@return The merged matrix includes the interactions, G1 and G2.
#'@import stats
#'@export
MPiG_MAIN_Check_Z<- function(G1,G2,D,nrow_G,iG_with,impute.method,missing_cutoff,max_MAF){
  Z=cbind(G1,G2)
  m <- ncol(Z)
  ncol_G1=NCOL(G1)
  ncol_G2<-NCOL(G2)
  id_mis<-union(which(is.na(G1)), which(G1== 9))
  if(length(id_mis) > 0){
    G1[id_mis] <- NA
  }


  in_include_g1<-apply(G1,2,function(gg){

    missing.ratio<-length(which(is.na(gg))/nrow_G);
    sd1 <- sd(gg, na.rm=TRUE)
    maf1 <- min(mean(gg, na.rm=TRUE)/2, 1- mean(gg, na.rm=TRUE)/2)


    if((missing.ratio < missing_cutoff) &&( sd1 > 0) && (maf1<max_MAF)){
      in_include_snp <- TRUE
    }else{
      in_include_snp <- FALSE
    }

    return(in_include_snp)
  })

  G11=G1[,in_include_g1]


  if(!is.null(G2)){

    id_mis<-union(which(is.na(G2)), which(G2== 9))
    if(length(id_mis) > 0){
      G2[id_mis] <- NA
    }


    in_include_g2<-apply(G2,2,function(gg){

      missing.ratio<-length(which(is.na(gg))/nrow_G);
      sd1 <- sd(gg, na.rm=TRUE)
      maf1 <- min(mean(gg, na.rm=TRUE)/2, 1- mean(gg, na.rm=TRUE)/2)

      if((missing.ratio < missing_cutoff) && (sd1 > 0) && (maf1<max_MAF)){
        in_include_snp <- TRUE
      }else{
        in_include_snp <- FALSE
      }

      return(in_include_snp)
    })
    G21<-G2[,in_include_g2]
    id_mis2<-which(is.na(G21));


    if(length(id_mis2)>0){
      G2_imp=Impute(G21,impute.method)
    }else{
      G2_imp=G21
    }

  }

  id_mis1<-which(is.na(G11));
  G1_imp<-ifelse(length(id_mis1)>0,Impute(G11,impute.method),G11)
  if(length(id_mis1)>0){
    G1_imp=Impute(G11,impute.method)
  }else{
    G1_imp=G11
  }

  G2_imp=as.matrix(G2_imp)
  G1_imp=as.matrix(G1_imp)


  if(iG_with=="G"){

    V<-apply(G1_imp,2,function(g1,g2){
      return(drop(g1)*g2)
    },G2_imp)

    V=matrix(V,nrow=nrow_G,ncol=ncol_G1*ncol_G2)

    re<-list(V=V,Z=cbind(G1_imp,G2_imp))
  }else{

    V<-apply(D,2,function(d,g1){return(drop(d)*g1)},G1_imp)

    re<-list(V=V,Z=cbind(D,G1_imp))
  }
  return(re)
}

#'@title Inverse matrix calculation.
#'@param V.item The inputted matrix.
#'@return Output the inversed matrix.

Get_GenInverse<-function(V.item){
  temp <- eigen(V.item)
  idx <- which(abs(temp$values) > 10^-5)

  if(length(idx)>1){
    V.item.inv <- temp$vectors[, idx] %*% diag(1/temp$values[idx]) %*%t(temp$vectors[, idx])
  }else{
    V.item.inv <- temp$vectors[, idx] %*% as.matrix(1/temp$values[idx])%*%t(temp$vectors[, idx])
  }
  return(V.item.inv)
}


mat.sqrt<-function (A)
{
  ei <- eigen(A)
  d <- ei$values
  d <- (d + abs(d))/2
  d2 <- sqrt(d)
  ans <- ei$vectors %*% diag(d2) %*% t(ei$vectors)
  return(ans)
}

#'@title Null model for GiMat association tests.
#'@param Y_mat The dichotomous phenotype matrix (N times K).
#'@param X Covariates, such as sex, age, PCs.
#'@param G1 The genotype matrix (N times p1) from the first phenotype.
#'@param G2 The genotype matrix (N times p2) from the second phenotype.
#'@param D The environmental factor. By defaul, D=NULL.
#'@param iG_with The variable used for indicating gene-by-gene interactions or gene-by-environment interactions. iG_with="G" refers to gene-by-gene interactions and iG_with="D" refers to gene-by-env interactions.
#'@import stats
#'@return A list consisting of residual, variance and fitted values and the combined matrix W<-cbind(rep(1,nrow_G),X,G1,G2) or W<-cbind(rep(1,nrow_G),X,G1,D).

MPiG_SKAT_Null<-function(Y_mat,X,G1,G2=NULL,D=NULL,iG_with="G"){
  nrow_G=NROW(G1)
  n.pheno<-NCOL(Y_mat)


  if(iG_with=="G"){
    W<-cbind(rep(1,nrow_G),X,G1,G2)
  }else{
    W<-cbind(rep(1,nrow_G),X,G1,D)
  }

  if(is.null(cbind(G1,G2,D))){
    Z=X
  }else{
    Z=cbind(X,G1,G2,D)
  }
  model=list()
  res_model=NULL#matrix(1,nrow=nrow_G,ncol=n.pheno)
  phi_model=NULL#matrix(1,nrow=nrow_G,ncol=n.pheno)
  mu_model=NULL#matrix(1,nrow=nrow_G,ncol=n.pheno)

    for (i in 1:n.pheno){
      #model[[i]]<-glm(Y_mat[,i]~X+Z-1,family=binomial(logit))
      model[[i]]<-glm(Y_mat[,i]~Z,family=binomial(logit))
      res_model<-cbind(res_model,model[[i]]$y-model[[i]]$fitted.values)
      phi_model<-cbind(phi_model,model[[i]]$fitted.values*(1-model[[i]]$fitted.values))
      mu_model<-cbind(mu_model,model[[i]]$fitted.values)
    }

  out=list(W=W,model=model,res_model=res_model,phi_model=phi_model,mu_model=mu_model)
  return(out)

}

#'@title GiMat association tests.
#'@param Y_mat The dichotomous phenotype matrix (N times K).
#'@param X Covariates, such as sex, age, PCs.
#'@param G1 The genotype matrix (N times p1) from the first phenotype.
#'@param G2 The genotype matrix (N times p2) from the second phenotype.
#'@param D The environmental factor. By defaul, D=NULL.
#'@param W The explainable variables in the null model. It contains intercept, covariates and genotype matrix (if the main genetic effects exist) and environment matrix (if main environmental effects exist).
#'@param model The complete output of the null model.
#'@param res_model The matrix in which each column contains the residual. res_model can be obtained from the output of null model.
#'@param phi_model The matrix in which each column contains the variance. phi_model can be obtained from the output of null model.
#'@param mu_model The matrix in which each column contains the mean. mu_model can be obtained from the output of null model.
#'@param Sigma_p The kernel matrix. HOM is a diagonal matrix, Het is a identical matrix. PC is the correlation matrix.
#'@param iG_with The variable used for indicating gene-by-gene interactions or gene-by-environment interactions. iG_with="G" refers to gene-by-gene interactions and iG_with="D" refers to gene-by-env interactions.
#'@param Is.Common To indicate whether gene-by-gene interactions are limited to common variants. For now, Is.Common=TRUE by default.
#'@param weights_beta A vector of length 2. Two parameters in Beta distribution: Beta(a1,a2)
#'@param weights_G Weights for genetic variants.
#'@param weights_V Weights for gene-by-gene or gene-by-env interactions.
#'@param impute.method Imputed method consists of "random", "fixed" and "bestguess".
#'@param missing_cutoff The threshold for the missing rate. If the missing ratio is less than the threshold, impute method is utilized.
#'@param max_MAF The threshold of maximum MAF.
#'@author Siru Wang
#'@export
#'@return Output the statistics, the corresponding p-values.
#'@import SKAT
#'@import copula
#'@import CompQuadForm
#'@import stats

MPiG_SKAT_Fast_robust<-function(Y_mat,X,G1,G2=NULL,D=NULL,
                                W,model,res_model,phi_model,mu_model,Sigma_p,
                                iG_with="G",
                                Is.Common=TRUE,
                                weights_beta=c(1,25),
                                weights_G = NULL,weights_V=NULL,impute.method = "fixed",
                                missing_cutoff=0.15,max_MAF=NULL){


  #Genotype.kernal function
  ####################################################################################################
  #Z=as.matrix(cbind(G1,G2));
  if(!is.matrix(G1)){
    G1=as.matrix(G1)
  }

  if(!is.matrix(G2)){
    G2=as.matrix(G2)
  }

  ncol_G1=NCOL(G1);ncol_G2=NCOL(G2);
  if (!is.null(D)) {D=as.matrix(D);ncol_D=NCOL(D)}
  nrow_G=NROW(G1)

  if(iG_with=="G"){
    Z.out<-MPiG_MAIN_Check_Z(as.matrix(G1),G2=as.matrix(G2),D=NULL,nrow_G,iG_with,impute.method,missing_cutoff,max_MAF)

  }else{
    Z.out<-MPiG_MAIN_Check_Z(as.matrix(G1),G2=NULL,D=as.matrix(D),nrow_G,iG_with,impute.method,missing_cutoff,max_MAF)

  }

  Z<-Z.out$Z
  V<-Z.out$V

  if(isTRUE(Is.Common)){
    weight_V=NULL;weight_G=NULL
    #  MAF_for_betaweights_V <- colMeans(as.matrix(V), na.rm=T)/2
    #  MAF_for_betaweights_V <- pmin(MAF_for_betaweights_V, 1-MAF_for_betaweights_V)

    #  weight_V=Beta_weights(MAF_for_betaweights_V,weights_beta)
    #  weight_G=NULL
  }else{
    MAF_for_betaweights_G <- colMeans(as.matrix(Z), na.rm=T)/2
    MAF_for_betaweights_G <- pmin(MAF_for_betaweights_G, 1-MAF_for_betaweights_G)
    weight_G<-Beta_weights(MAF_for_betaweights_G,weights_beta)
    weight_V=weight_G
  }


  V_w <- NULL
  if(is.null(weights_V)==FALSE){
    V_w <- t(t(V) * (weights_V))
  }else{
    V_w<-V
  }

  # V_w <- colSums(t(V) * (weight_V))
  # V_w=as.matrix(V_w)

  X<- as.matrix(X)
  n.pheno<-NCOL(Y_mat)
  p_Z <- NCOL(Z)
  q_X <- NCOL(X) + 1
  # +1 is for intercept,

    M1<-apply(phi_model,2,function(phi,w){
      return(as.vector(phi)*w)
    },W)

    M1_temp<-tapply(M1,rep(1:n.pheno,each=NROW(W)*NCOL(W)),function(m1){return(m1)})


    P_list<-lapply(M1_temp,function(m1,w){

      m1_mat=matrix(m1,nrow=NROW(W),ncol=NCOL(W))
      wdw=t(w)%*%m1_mat
      #wdw_p=wdw+temp
      #p=diag(n+1)-w%*%solve(wdw)%*%t(m1_mat)
      p=diag(nrow_G)-w%*%Get_GenInverse(wdw)%*%t(m1_mat)
      return(p)

    },W)

    #IkP<-Matrix::bdiag(P_list)
    #IkPt<-Matrix::bdiag(t(P_list))
    IkP=matrix(0,nrow=n.pheno*nrow_G,ncol=n.pheno*nrow_G)
    for(i in 1:n.pheno){
      IkP[1:nrow_G+(i-1)*nrow_G,1:nrow_G+(i-1)*nrow_G]=P_list[[i]]
    }
    IkPt=t(IkP)

  #  print("test which part is o-o-m1")

  Sigma_y=cov(as.matrix(res_model))
  #  print("test which part is o-o-m2")
  Sigma_y_inv=Get_GenInverse(Sigma_y)
  #  print("test which part is o-o-m3")
  R1=mat.sqrt(Sigma_p)
  #  print("test which part is o-o-m4")

  V_list=lapply(as.list(1:n.pheno),function(ii){return(V_w)})
  IV=matrix(0,nrow=n.pheno*NROW(V_w),ncol=n.pheno*NCOL(V_w))
  for(i in 1:n.pheno){
    IV[(i-1)*NROW(V_w)+(1:NROW(V_w)),(i-1)*NCOL(V_w)+(1:NCOL(V_w))]=V_list[[i]]
  }
  IVt=t(IV)
  #  print("test which part is o-o-m5")


  #R1Siginv_I=(R1%*%Sigma_y_inv)%x%diag(nrow_G)
  R1Siginv_I=kronecker(R1%*%Sigma_y_inv,diag(nrow_G))
  #  print("test which part is o-o-m6")
  SiginvR1_I=(Sigma_y_inv%*%R1)%x%diag(nrow_G)
  #  print("test which part is o-o-m7")

  SigI=Sigma_y%x%diag(nrow_G)
  #  print("test which part is o-o-m8")

  eig_chi=IVt%*%R1Siginv_I%*%IkP%*%SigI%*%IkPt%*%SiginvR1_I%*%IV
  #  print("test which part is o-o-m9")
  #########################################################################################
  diff_Ymu=matrix(c(Y_mat-mu_model),ncol=1)
  Q=t(diff_Ymu)%*%((Sigma_y_inv%*%Sigma_p%*%Sigma_y_inv)%x%(V%*%t(V)))%*%diff_Ymu
  lambda_skat <- SKAT::Get_Lambda(eig_chi)
  #type I error simulation acc=10^{-6}
  out<-davies(Q,lambda_skat,lim = 10^(5),acc=10^(-6), rep(1,length(lambda_skat)))
  #power simulation acc=10
  #out<-davies(Q,lambda_skat,lim = 10^(5),acc=10^(-10), rep(1,length(lambda_skat)))
  pvalue <- out$Qq

  return(list(Q = Q,p.value =pvalue,R = Sigma_p,eig_chi=eig_chi))
}
