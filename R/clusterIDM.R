

#' @importFrom numDeriv grad
#' @importFrom msm dpexp ppexp
#' @importFrom utils combn
#' @importFrom eha ppch
#'
#' @export
clusterIDM <- function(outdata.F, outdata.proband, outdata.R, outdata.S = NULL, Cr_R, Ar_R, Cr_S =NULL, lam03, Age = NULL, Cal = NULL, fam.id, birth, design, no.death = FALSE){

if(is.null(Age)) Age = c(0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90)
if(is.null(Cal)) Cal = seq(1920, 2011, by=5)

nr <- nrow(outdata.R)
A.R <- Ar_R - outdata.R[, birth]
C.B0 <- lapply(1:nr, function(i) outdata.R[i, birth]+ Age)
C <- lapply(1:nr, function(i) sort(c(Cal, C.B0[[i]])))

Cf.F <- lapply(1:nr, function(i) unique(C[[i]][outdata.R[i, birth] <= C[[i]] & C[[i]] <= Ar_R]))
Af.F <- lapply(1:nr, function(i) Cf.F[[i]] - outdata.R[i, birth])

R.f <- lapply(1:nr, function(i) sapply(1:length(Cf.F[[i]]), function(k) findInterval(Cf.F[[i]][k], Cal)))
A.f <- lapply(1:nr, function(i) sapply(1:length(Af.F[[i]]), function(k) findInterval(Af.F[[i]][k], Age)))

LAM03.R <-  lapply(1:nr, function(j) sapply(1:length(R.f[[j]]), function(i) lam03[lam03$Year.f==R.f[[j]][i] & lam03$Age.f==A.f[[j]][i], ]$rate))
cut.R <- Af.F

LAM03.S <- cut.S <- NULL
if(!is.null(outdata.S)){
ns <- nrow(outdata.S)
A.S <- Cr_S - outdata.S[, birth]
C.B0 <- lapply(1:nr, function(i) outdata.S[i, birth]+ Age)
C <- lapply(1:nr, function(i) sort(c(Cal, C.B0[[i]])))

Cf.F <- lapply(1:nr, function(i) unique(C[[i]][outdata.S[i, birth] <= C[[i]] & C[[i]] <= Cr_S]))
Af.F <- lapply(1:nr, function(i) Cf.F[[i]] - outdata.S[i, birth])

R.f <- lapply(1:nr, function(i) sapply(1:length(Cf.F[[i]]), function(k) findInterval(Cf.F[[i]][k], Cal)))
A.f <- lapply(1:nr, function(i) sapply(1:length(Af.F[[i]]), function(k) findInterval(Af.F[[i]][k], Age)))

cut.S <- Af.F
LAM03.S <-  lapply(1:ns, function(j) sapply(1:length(R.f[[j]]), function(i) lam03[lam03$Year.f==R.f[[j]][i] & lam03$Age.f==A.f[[j]][i], ]$rate))
}

outdata.F0 <- split(outdata.F, as.factor(outdata.F[, fam.id]))

if(no.death == TRUE) par<- c(log(lam01), log(rho))
else par<- c(log(lam01), log(theta), log(rho))


pairle <- optim(par, pair.logL, outdata.F = outdata.F0,  outdata.proband = outdata.proband, outdata.R = outdata.R, outdata.S = outdata.S,
                Cr_R = Cr_R, LAM03.R = LAM03.R, cut.R = cut.R, LAM03.S = LAM03.S, cut.S = cut.S, Age = Age, Cal = Cal, design = design, no.death = no.death, method = "BFGS", hessian = TRUE)

parameter.pair <- exp(pairle$par)

if(design == 1){
  if(no.death == TRUE) {

    score_i <- sapply(1:nf, function(i) numDeriv::grad(func=NloglikFD1,  x=pairle$par, outdata_F = outdata.F0[i], outdata_proband=outdata.proband[i,Cr_R],
                                                       Age = Age, Cal = Cal, lam03 = lam03,
                                                       fgau = gauleg.f, fdpexp = msm::dpexp, fppexp = msm::ppexp, combn = utils::combn))


  }else {

    score_i <- sapply(1:nf, function(i) numDeriv::grad(func=loglikFD1,  x=pairle$par, outdata_F = outdata.F0[i], outdata_proband=outdata.proband[i,Cr_R],
                                                       Age = Age, Cal = Cal, lam03 = lam03,
                                                      fgau = gauleg.f, fdpexp = msm::dpexp, fppexp = msm::ppexp, combn = utils::combn))

    }
  }else{
  if(no.death == TRUE) {

    score_i <- sapply(1:nf, function(i) numDeriv::grad(func=NloglikFD2,  x=pairle$par, outdata_F = outdata.F0[i], outdata_proband=outdata.proband[i,Cr_R],
                                                       Age = Age, Cal = Cal, lam03 = lam03,
                                                       fgau = gauleg.f, fdpexp = msm::dpexp, fppexp = msm::ppexp, combn = utils::combn, ppch = eha::ppch))


  }else {

    score_i <- sapply(1:nf, function(i) numDeriv::grad(func=loglikFD2,  x=pairle$par, outdata_F = outdata.F0[i], outdata_proband=outdata.proband[i,Cr_R],
                                                       Age = Age, Cal = Cal, lam03 = lam03,
                                                       fgau = gauleg.f, fdpexp = msm::dpexp, fppexp = msm::ppexp, combn = utils::combn, ppch = eha::ppch))

  }
}
score_r <- sapply(1:nr, function(i) numDeriv::grad(loglikR, x=pairle$par,  outdata_R = outdata.R[i,], LAM03R = LAM03.R[i], cutR = cut.R[i], fgau = gauleg.f, fdpexp = msm::dpexp, fppexp = msm::ppexp))

if(!is.null(outdata.S)) {
  if(no.death == TRUE) {
  score_s <- sapply(1:ns, function(i) numDeriv::grad(NloglikS, x=pairle$par,  outdata_S = outdata.S[i,], LAM03S = LAM03.S[i], cutS = cut.S[i], fgau = gauleg.f, fdpexp = msm::dpexp, fppexp = msm::ppexp))
  }else{
    score_s <- sapply(1:ns, function(i) numDeriv::grad(loglikS, x=pairle$par,  outdata_S = outdata.S[i,], LAM03S = LAM03.S[i], cutS = cut.S[i], fgau = gauleg.f, fdpexp = msm::dpexp, fppexp = msm::ppexp))
  }
}
B <- score_i%*%t(score_i) + score_r%*%t(score_r)
if(!is.null(outdata.S)) B <- B +  score_s%*%t(score_s)

A <- pairle$hessian
var <- solve(A)%*%B%*%t(solve(A))
sd.pair <- sqrt(diag(var))

result <- list()
result$estimate <- parameter.pair
result$var <- var
return(result)
}


