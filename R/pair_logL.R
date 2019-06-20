#===========================================
# Full history of Proband : Pairwise + Auxiliary
#=========================================

pair.logL <-	function( par, Y.fam, X.fam,  Y.proband, X.proband, Y.R, X.R, Y.S,
                        outdata.F0 = NULL, outdata.proband = NULL, cut, LAM03.R, cut.R, LAM03.S = NULL, cut.S = NULL,  Age, Cal, design, full, no.death){


  if(design == 1){
    if(no.death == TRUE) res <- NloglikFD1(par, outdata.F0, outdata.proband, Age, Cal, lam03, full, gauleg.f, msm::dpexp, msm::ppexp, utils::combn)
    else res <- loglikFD1_pch(par, cut, outdata.F0, outdata.proband, Age, Cal, lam03, full, gauleg.f, utils::combn)



  }else{
    if(no.death == TRUE) res <- NloglikFD2(par, Y.fam, X.fam, as.matrix(Y.proband), as.matrix(X.proband),  Age, Cal, cut, lam03, gauleg.f, utils::combn)

    else res <- loglikFD2_pch(par, Y.fam, X.fam, as.matrix(Y.proband), as.matrix(X.proband),  Age, Cal, cut, lam03, gauleg.f, utils::combn)



  }

  auxtmp1 <- 0
  if(no.death == FALSE) auxtmp1 <- loglikR_pch( par, cut, Y.R, X.R,  LAM03.R, cut.R, gauleg.f)
  auxtmp2 <- 0
  if(!is.null(LAM03.S) & !is.null(cut.S)) {
    if(no.death == TRUE) fitter <- get("NloglikS")
    else fitter <- get("loglikS_pch")

    auxtmp2 <- fitter(par, cut, Y.S, LAM03.S, cut.S, gauleg.f)

  }
  res = res + auxtmp1 + auxtmp2

  print(c(c(exp(par[1]), par[2], exp(par[-c(1:2)])), res))

  return(-res)
}


