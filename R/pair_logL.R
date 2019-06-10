#===========================================
# Full history of Proband : Pairwise + Auxiliary
#=========================================
pair.logL <-	function(outdata.F, outdata.proband, outdata.R, outdata.S = NULL, par, Cr_R, LAM03.R, cut.R, LAM03.S = NULL, cut.S = NULL,  Age, Cal, design){

  if(design == 1){
    res <- loglikFD1(par, outdata.F, outdata.proband[, Cr_R], Age, Cal, lam03, gauleg.f, msm::dpexp, msm::ppexp, utils::combn)

  }else{
    res <- loglikFD2(par, outdata.F, outdata.proband[, Cr_R], Age, Cal, lam03, gauleg.f, msm::dpexp, msm::ppexp, utils::combn, eha::ppch)

  }

  auxtmp1 <- loglikR(outdata.R, par, LAM03.R, cut.R, gauleg.f, msm::dpexp, msm::ppexp)

  auxtmp2 <- 0
  if(!is.null(LAM03.S) & !is.null(cut.S)) auxtmp2 <- loglikS(outdata.S, par, LAM03.S, cut.S, gauleg.f, msm::dpexp, msm::ppexp)

  res = res + auxtmp1 + auxtmp2

  print(c(exp(par), res))

  return(-res)
}


