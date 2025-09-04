
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GiMat

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17055370.svg)](https://doi.org/10.5281/zenodo.17055370)
<!-- badges: end -->

This is an R package for testing the multiple phenotype associatin with gene-by-gene (G$\times$G) interactions, which accounts for the  relationship among the causal interaction effects on each of the phenotypes.

## System requirements
GiMat could be installed with ease on versions of R > 4.0.2.


## Installation

Before installing R package of GiMat, the following packages from CRAN are required:

``` r
install.packages("CompQuadForm")
install.packages("SKAT")
install.packages("copula")
install.packages("stats")
```
You can install the GitHub version of GiMat:
``` r
# install.packages("devtools")
devtools::install_github("SiruRooney/GiMat")
```
## Example
Here's an example workflow using GiMat:


### Data Preparation

First, prepare multiple dichotomous phenotypes and the common variant sets that are used for testing gene-by-gene interaction associations with these phenotypes. `GiMat` R package provides an example.  

``` r
library(GiMat)
#> Warning: replacing previous import 'copula::profile' by 'stats::profile' when
#> loading 'GiMat'
#> Warning: replacing previous import 'copula::logLik' by 'stats::logLik' when
#> loading 'GiMat'
#> Warning: replacing previous import 'copula::confint' by 'stats::confint' when
#> loading 'GiMat'
#> Warning: replacing previous import 'copula::coef' by 'stats::coef' when loading
#> 'GiMat'
```

```r
summary(Covariates)
#     sex              age       
# Min.   :0.0000   Min.   :38.58  
# 1st Qu.:0.0000   1st Qu.:45.99  
# Median :0.0000   Median :48.02  
# Mean   :0.4794   Mean   :48.03  
# 3rd Qu.:1.0000   3rd Qu.:50.05  
# Max.   :1.0000   Max.   :57.82 
apply(cbind(pheno1,pheno2),2,table)
#    pheno1          pheno2      
#0    4015            4011
#1    985             989

table(geno1)
#geno1
#   0    1    2 
# 6029 3310  661
table(geno2)
#geno2
#    0     1     2 
# 18525  5294  1181 
```

Then we construct the different kernels for modelling the phenotype covariates.
``` r
Y_mat=cbind(pheno1,pheno2)
# Heterogeneous (Het) kernel
Sigma_het=diag(2)
# Homogeneous (Hom) kernel
Sigma_hom=matrix(1,nrow=2,ncol=2)
# Phenotype covariance (PhC) kernel
Sigma_phc=binary_cov_matrix(Y_mat)
# Gene association with multi-trait (GAM) kernel
Sigma_gam=binary_cov_matrix(Y_mat)^2
```

### Running

In this example, we assume that the genetic variants, `geno1` and `geno2`, are not associated with both phenotypes but the gene-by-gene interactions are associated with phenotypes.

We conduct the null model under the null hypothesis that the gene-by-gene interactions have no effects on the phenotypes.


#### Inputted data for `MPiG_SKAT_Null` function

- `Y_mat`: Phenotype matrix. Each column refers to one dichotomous phenotype.

- `X`: Covariates including age, sex, PCs.

- `G1`: The first genotype matrix. Each column refers to one variant.

- `G2`: The second genotype matrix. Each column refers to one variant. By default, `G2` is NULL.

- `D`: The environmental exposure. If we test the gene-by-env interactions associated with the environmental exposure, we input the individual-level environmental exposure and set `G2` as NULL. By default, `D` is NULL.

- `iG_with`: An indicator which represents the gene-by-gene interactions or gene-by-env interactions. `iG_with` is set to "G" for gene-by-gene and "D" for gene-by-env. By default, `iG_with` is set to "G".


``` r
out<-MPiG_SKAT_Null(Y_mat,X=Covariates,G1=NULL,G2=NULL,D=NULL,iG_with="G")
```

Then we conduct the gene-by-gene interaction association tests under different kernel assumptions. 

#### Inputted data for `MPiG_SKAT_Fast_robust` function

There are additional options in `MPiG_SKAT_Fast_robust` function except some overlapping inputted options appearing in `MPiG_SKAT_Null` function. 

- `W`: The explainable variables in the null model. It contains intercept, covariates and genotype matrix (if the main genetic effects exist) and environment matrix (if main environmental effects exist).

- `model`: The complete output from the null model.

- `res_model`: Residuals under the null model.

- `phi_model`: Variance under the null model.

- `mu_model`: The fitted values (means) under the null model.

- `Sigma_p`: Phenotype covariance kernel.

- `Is.Common`: To indicate the variants are common or rare. If variants are common, `Is.Common` is TRUE. The current version is restricted to common variant interactions.

- `weights_beta`, `weights_G` and `weights_V`: The weights for allelic effect, variants and interactions. These options are applicable for gene-by-env interaction association tests. The current version is restricted to common variant interactions.

- `impute.method`: The imputed method for imputting missing values of the variants. By default, `impute.method` is set to "fixed".

- `missing_cutoff`: The threshold for variants having missing values. By default,`missing_cutoff=0.15`.

- `max_MAF`: The maximum of MAF. By default, `max_MAF=0.5`.

``` r
#Under heterogeneous (Het) kernel assumption
Gimat.het=MPiG_SKAT_Fast_robust(Y_mat,X=Covariates,G1=geno1,G2=geno2,D=NULL,W=out$W,model=out$model,res_model=out$res_model,phi_model=out$phi_model,mu_model=out$mu_model,Sigma_p=Sigma_het,iG_with="G",Is.Common=TRUE,weights_beta=c(1,25),weights_G=NULL,weights_V=NULL,impute.method="fixed",missing_cutoff=0.15,max_MAF=0.5)

#Under homogeneous (Hom) kernel assumption
Gimat.hom=MPiG_SKAT_Fast_robust(Y_mat,X=Covariates,G1=geno1,G2=geno2,D=NULL,W=out$W,model=out$model,res_model=out$res_model,phi_model=out$phi_model,mu_model=out$mu_model,Sigma_p=Sigma_hom,iG_with="G",Is.Common=TRUE,weights_beta=c(1,25),weights_G=NULL,weights_V=NULL,impute.method="fixed",missing_cutoff=0.15,max_MAF=0.5)

#Under the phenotype covariance (phc) kernel assumption
Gimat.phc=MPiG_SKAT_Fast_robust(Y_mat,X=Covariates,G1=geno1,G2=geno2,D=NULL,W=out$W,model=out$model,res_model=out$res_model,phi_model=out$phi_model,mu_model=out$mu_model,Sigma_p=Sigma_phc,iG_with="G",Is.Common=TRUE,weights_beta=c(1,25),weights_G=NULL,weights_V=NULL,impute.method="fixed",missing_cutoff=0.15,max_MAF=0.5)

#Under the gene association with multi-trait (GAM) kernel
Gimat.gam=MPiG_SKAT_Fast_robust(Y_mat,X=Covariates,G1=geno1,G2=geno2,D=NULL,W=out$W,model=out$model,res_model=out$res_model,phi_model=out$phi_model,mu_model=out$mu_model,Sigma_p=Sigma_gam,iG_with="G",Is.Common=TRUE,weights_beta=c(1,25),weights_G=NULL,weights_V=NULL,impute.method="fixed",missing_cutoff=0.15,max_MAF=0.5)

#Minimum P value‐based omnibus kernel tests

obj.list=list(Gimat.het,Gimat.hom,Gimat.phc,Gimat.gam)
#Extract the gene-by-gene or gene-by-env interactions
out.iG=MPiG_MAIN_Check_Z(G1=geno1,G2=geno2,D=NULL,nrow_G=NROW(geno1),iG_with="G",impute.method="fixed",missing_cutoff=0.15,max_MAF=0.5)

iG=out.iG$V
Gimat.minP=MPiG_minP_base(n.pheno=NCOL(Y_mat),n=NROW(Y_mat),res_model=out$res_model,obj.list=obj.list,Z1=iG,X=covariates, resample = 200)

```
## Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Siru Wang](https://github.com/SiruRooney)


