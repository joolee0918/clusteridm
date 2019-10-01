#===========================================
# Full history of Proband : Pairwise + Auxiliary
#=========================================

pair.logL <-	function( par, Y.fam, X.fam,  Y.proband, X.proband, Y.R, X.R, Y.S,
                        outdata.F0 = NULL, outdata.proband = NULL, cut, R.fR, A.fR, cut.R, R.fS, A.fS, cut.S = NULL,  G, Age, Cal, design, full, no.death){


  #if(design == 1){
   # res <- loglikFD1_pch(par, cut, outdata.F0, outdata.proband, Age, Cal, lam03, full, gauleg.f, utils::combn)

  theta <- ifelse(lam03$Age.f < 14, exp(par[2]), 0)
  theta <- ifelse(lam03$Age.f <16 & lam03$Age.f >=14 , exp(par[3]), theta)
  theta <- ifelse( lam03$Age.f >=16 , exp(par[4]), theta)

  lam12 <- lam03$rate*theta
  nr <- length(cut.R)
  LAM03.R <-  lapply(1:nr, function(j) sapply(1:length(R.fR[[j]]), function(i) lam03[lam03$Year.f==R.fR[[j]][i] & lam03$Age.f==A.fR[[j]][i], ]$rate))
  LAM12.R <-  lapply(1:nr, function(j) sapply(1:length(R.fR[[j]]), function(i) lam12[lam03$Year.f==R.fR[[j]][i] & lam03$Age.f==A.fR[[j]][i]]))

  if(!is.null(cut.S)){
    ns <- length(cut.S)
    LAM03.S <-  lapply(1:ns, function(j) sapply(1:length(R.fS[[j]]), function(i) lam03[lam03$Year.f==R.fS[[j]][i] & lam03$Age.f==A.fS[[j]][i], ]$rate))
    LAM12.S <-  lapply(1:ns, function(j) sapply(1:length(R.fS[[j]]), function(i) lam12[lam03$Year.f==R.fS[[j]][i] & lam03$Age.f==A.fS[[j]][i]]))
  }

  if(is.null(G)) res <- loglikFD2_pch(par, theta, Y.fam, X.fam, as.matrix(Y.proband), as.matrix(X.proband),  Age, Cal, cut, lam03, gauleg.f, utils::combn)
  else res <- loglikFD2_pch_gene(par, theta, Y.fam, X.fam, as.matrix(Y.proband), as.matrix(X.proband),  Age, Cal, cut, lam03, gauleg.f, utils::combn)


  auxtmp1 <- 0

  if(is.null(G)) auxtmp1 <- loglikR_pch( par,  cut, Y.R, X.R,  LAM03.R, LAM12.R, cut.R, gauleg.f)
  else auxtmp1 <- loglikR_pch_gene( par,  cut, Y.R, X.R,  LAM03.R, LAM12.R, cut.R, gauleg.f)

  auxtmp2 <- 0
  if(!is.null(LAM03.S) & !is.null(cut.S)) {
    if(is.null(G)) auxtmp2 <- loglikS_pch(par,  cut, Y.S, LAM03.S, LAM12.S, cut.S, gauleg.f)
    else auxtmp2 <- loglikS_pch_gene(par,  cut, Y.S, LAM03.S, LAM12.S, cut.S, gauleg.f)
  }
  res = res + auxtmp1 + auxtmp2

  return(-res)
}



loglikFD2_pch_R <- function(par, Y_F, X_F,  Y_proband, X_proband,
                 Age, Cal,   cut_F, lam03, fgau, combn){

  theta <- ifelse(lam03$Age.f < 14, exp(par[2]), 0)
  theta <- ifelse(lam03$Age.f <16 & lam03$Age.f >=14 , exp(par[3]), theta)
  theta <- ifelse( lam03$Age.f >=16 , exp(par[4]), theta)

  res <- loglikFD2_pch(par, theta, Y_F, X_F, Y_proband, X_proband,  Age, Cal, cut_F, lam03, fgau, combn)
  return(res)
}
                                                
loglikFD2_pch_gene_R <- function(par, Y_F, X_F,  Y_proband, X_proband,
                 Age, Cal,   cut_F, lam03, fgau, combn){

  theta <- ifelse(lam03$Age.f < 14, exp(par[2]), 0)
  theta <- ifelse(lam03$Age.f <16 & lam03$Age.f >=14 , exp(par[3]), theta)
  theta <- ifelse( lam03$Age.f >=16 , exp(par[4]), theta)

  
  res <- loglikFD2_pch_gene(par, theta, Y_F, X_F, Y_proband, X_proband,  Age, Cal, cut_F, lam03, fgau, combn)
  return(res)
}

                                                

loglikR_pch_R <- function(par, cut_F, Y_R, X_R,  R.fR, A.fR, cutR, fgau){

  theta <- ifelse(lam03$Age.f < 14, exp(par[2]), 0)
  theta <- ifelse(lam03$Age.f <16 & lam03$Age.f >=14 , exp(par[3]), theta)
  theta <- ifelse( lam03$Age.f >=16 , exp(par[4]), theta)

  lam12 <- lam03$rate*theta
  LAM03.R <- LAM12.R <- list()
  LAM03.R[[1]] <- sapply(1:length(R.fR), function(i) lam03[lam03$Year.f==R.fR[i] & lam03$Age.f==A.fR[i], ]$rate)
  LAM12.R[[1]] <- sapply(1:length(R.fR), function(i) lam12[lam03$Year.f==R.fR[i] & lam03$Age.f==A.fR[i]])


  res <- loglikR_pch(par,  cut_F, Y_R, X_R,  LAM03.R[1], LAM12.R[1], cutR, fgau)

  return(res)
}

loglikR_pch_gene_R <- function(par, cut_F, Y_R, X_R,  R.fR, A.fR, cutR, fgau){

  theta <- ifelse(lam03$Age.f < 14, exp(par[2]), 0)
  theta <- ifelse(lam03$Age.f <16 & lam03$Age.f >=14 , exp(par[3]), theta)
  theta <- ifelse( lam03$Age.f >=16 , exp(par[4]), theta)

  lam12 <- lam03$rate*theta
  LAM03.R <- LAM12.R <- list()
  LAM03.R[[1]] <- sapply(1:length(R.fR), function(i) lam03[lam03$Year.f==R.fR[i] & lam03$Age.f==A.fR[i], ]$rate)
  LAM12.R[[1]] <- sapply(1:length(R.fR), function(i) lam12[lam03$Year.f==R.fR[i] & lam03$Age.f==A.fR[i]])


  res <- loglikR_pch_gene(par, cut_F, Y_R, X_R,  LAM03.R[1], LAM12.R[1], cutR, fgau)

  return(res)
}

loglikS_pch_R <- function(par, cut_F, Y_S, X_S,  R.fS, A.fS, cutS, fgau){

  theta <- ifelse(lam03$Age.f < 14, exp(par[2]), 0)
  theta <- ifelse(lam03$Age.f <16 & lam03$Age.f >=14 , exp(par[3]), theta)
  theta <- ifelse( lam03$Age.f >=16 , exp(par[4]), theta)

  lam12 <- lam03$rate*theta

  LAM03.S <- LAM12.S <- list()
  LAM03.S[[1]] <- sapply(1:length(R.fS), function(i) lam03[lam03$Year.f==R.fS[i] & lam03$Age.f==A.fS[i], ]$rate)
  LAM12.S[[1]] <- sapply(1:length(R.fS), function(i) lam12[lam03$Year.f==R.fS[i] & lam03$Age.f==A.fS[i]])


  res <- loglikS_pch( par, cut_F,  Y_S, LAM03.S[1], LAM12.S[1], cutS, fgau)


  return(res)
}

loglikS_pch_gene_R <- function(par, cut_F, Y_S, X_S,  R.fS, A.fS, cutS, fgau){

  theta <- ifelse(lam03$Age.f < 14, exp(par[2]), 0)
  theta <- ifelse(lam03$Age.f <16 & lam03$Age.f >=14 , exp(par[3]), theta)
  theta <- ifelse( lam03$Age.f >=16 , exp(par[4]), theta)

  lam12 <- lam03$rate*theta
  LAM03.S <- LAM12.S <- list()
  LAM03.S[[1]] <- sapply(1:length(R.fS), function(i) lam03[lam03$Year.f==R.fS[i] & lam03$Age.f==A.fS[i], ]$rate)
  LAM12.S[[1]] <- sapply(1:length(R.fS), function(i) lam12[lam03$Year.f==R.fS[i] & lam03$Age.f==A.fS[i]])


  res <- loglikS_pch_gene( par,  cut_F,  Y_S, LAM03.S[1], LAM12.S[1], cutS, fgau)


  return(res)
}



