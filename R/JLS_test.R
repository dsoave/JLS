#' Joint Location Scale (JLS) Test
#'
#' This function performs the Joint Location Scale (JLS) test (Soave et al. 2015) to simultaneously test for mean and variance differences between groups.  The JLS test uses Fisher's combined p-value method to combine evidence from the individual locaiton (regression t-test) and scale (Levene's test of homogeneity of variances) tests.
#' @param y a qunatitative outcome variable
#' @param x.loc a design matrix (or vector) for the location model
#' @param x.scale a design matrix (or vector) for the scale model
#' @keywords JLS
#' @export
#' @author David Soave
#' @import quantreg
#' @details No missing data are allowed - function will return an "error".  Absolute residuals, used in Levene's test (1960), are estimated using least absolute deviation (LAD) regression. LAD residuals correspond to deviations from group medians in the presence of a single categorical covariate. Outcome (phenotype) must be quantitative and covariate (genotype) must be discrete (categorical).
#' @return p_L the location test (regression t-test) p-value
#' @return p_S the scale test (Levene's test) p-value
#' @return p_JLS the JLS test (Fisher's combined method) p-value
#' @references Soave, D., Corvol, H., Panjwani, N., Gong, J., Li, W., Boelle, P.Y., Durie, P.R., Paterson, A.D., Rommens, J.M., Strug, L.J., and Sun, L. (2015). A Joint Location-Scale Test Improves Power to Detect Associated SNPs, Gene Sets, and Pathways. American journal of human genetics 97, 125-138.
#' @examples
#' #################################################################################
#' ## Example simulating data from model [i] (Soave et al. 2015 AJHG)
#' #################################################################################
#'
#' n<-2000  ## total sample size
#' pA<-0.3  ## MAF
#' pE1<-0.3  ## frequency of exposure E1
#'
#' ## Genotypes (XG)
#' genocount<-rmultinom(1,size=n,prob=c(pA*pA, 2*pA*(1-pA), (1-pA)*(1-pA)))
#' XG<-c(rep(0, genocount[1]), rep(1, genocount[2]), rep(2,genocount[3]))
#' XG<-sample(XG,size=length(XG),replace=FALSE)
#'
#' ## Exposures (E1)
#' E1<-rbinom(n,1,prob=pE1)
#'
#' ## Phenotype (y)
#' y<-0.01*XG+0.3*E1+0.5*XG*E1+rnorm(n,0,1)
#'
#'# Additive model (will work with dosages)
#'JLS_test(y,XG,XG)
#'# Genotypic model (will not work with dosages --> factor() will create many groups)
#'JLS_test(y,factor(XG),factor(XG))
#'
#'X2=round(cbind(XG==1,XG==2)) #convert to 2 columns
#'
#'# Genotypic model --> same result as results above using JLS_test(y,factor(X),factor(X))
#'# This is how genotype probabilities will be analyzed (using a 2 column design matrix)
#'JLS_test(y,X2,X2)



#Note that p_S is obtained using a stage 1 median regression (rq function, tau=0.5) where group medians are chosen to be the larger of two middle values when the group size is even.
JLS_test <-function(y,x.loc,x.scale){

  ## check if there is missing data
  data <- cbind(y, x.loc,x.scale)
  if(sum(is.na(data)) > 0)  stop("missing value(s) not allowed")


  ## Obtain the p-values from the individual location and scale test
  p_L <- anova(lm(y~x.loc))[1,5]

  lm1 <- rq(y~x.scale,tau=.5)
  d1 <- abs(y-predict(lm1))
  p_S <- anova(lm(d1~x.scale))[1,5]

  ## JLS test statistic (and corresponding p-value) using Fisher's combined p-value method
  t_JLS <- -2*log(p_L)-2*log(p_S)
  p_JLS <- 1-pchisq(t_JLS,4)

  ## return the location, scale and JLS p-values
  return(cbind(p_L,p_S,p_JLS))
}


