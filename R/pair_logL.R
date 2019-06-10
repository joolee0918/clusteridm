#===========================================
# Full history of Proband : Pairwise + Auxiliary
#=========================================
pair.logL <-	function(outdata.F, outdata.proband, outdata.R, outdata.S = NULL, par, Cr_R, LAM03.R, cut.R, LAM03.S = NULL, cut.S = NULL,  Age, Cal, design, no.death){

  if(design == 1){
    if(no.death == TRUE) fitter <- get("NloglikFD1")
    else fitter <- get("loglikFD1")

    res <- fitter(par, outdata.F, outdata.proband[, Cr_R], Age, Cal, lam03, gauleg.f, msm::dpexp, msm::ppexp, utils::combn)

  }else{
    if(no.death == TRUE) fitter <- get("NloglikFD2")
    else fitter <- get("loglikFD2")

    res <- fitter(par, outdata.F, outdata.proband[, Cr_R], Age, Cal, lam03, gauleg.f, msm::dpexp, msm::ppexp, utils::combn, eha::ppch)

  }

  auxtmp1 <- 0
  if(no.death == FALSE) auxtmp1 <- loglikR(outdata.R, par, LAM03.R, cut.R, gauleg.f, msm::dpexp, msm::ppexp)

  auxtmp2 <- 0
  if(!is.null(LAM03.S) & !is.null(cut.S)) {
    if(no.death == TRUE) fitter <- get("NloglikS")
    else fitter <- get("loglikS")

    auxtmp2 <- fitter(outdata.S, par, LAM03.S, cut.S, gauleg.f, msm::dpexp, msm::ppexp)

  }
  res = res + auxtmp1 + auxtmp2

  print(c(exp(par), res))

  return(-res)
}


