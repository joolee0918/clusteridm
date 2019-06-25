

#' @importFrom numDeriv grad
#' @importFrom msm dpexp ppexp
#' @importFrom utils combn
#' @importFrom eha ppch
#' @importFrom lubridate ymd ydate year
#'
#' @export
clusterIDM <- function(fam.formula, R.formula, S.formula,
                       outdata.fam, outdata.R, outdata.S = NULL,
                       fam.id, fam.rel, recruit.age.fam,
                       first.visit.age.R, R.id, recruit.age.S =NULL,
                       birth, cut=0, lam03, Age = NULL, Cal = NULL, design, full = FALSE, no.death = FALSE, init = NULL){

  new.fam.formula <- update.formula(fam.formula,  paste("~.+", paste(fam.id, fam.rel, recruit.age.fam, birth, sep="+")))
  new.R.formula <- update.formula(R.formula,  paste("~.+", paste(first.visit.age.R, birth, sep="+")))

  #outdata.fam[, birth] <- lubridate::year(outdata.fam[, birth]) + lubridate::yday(outdata.fam[,birth])/365
  #outdata.R[, birth] <- lubridate::year(outdata.R[, birth]) + lubridate::yday(outdata.R[,birth])/365
  outdata.proband <- outdata.R[outdata.R[,R.id] %in% unique(outdata.fam[, fam.id]), ]
  outdata.R <- outdata.R[!outdata.R[,R.id] %in% unique(outdata.fam[, fam.id]), ]

  m.fam <- model.frame(new.fam.formula, data = outdata.fam)
  m.proband <- model.frame(new.R.formula, data = outdata.proband)
  m.R <- model.frame(new.R.formula, data = outdata.R)

  Y.fam <- as.matrix(model.extract(m.fam, "response"))
  Y.proband <- as.matrix(model.extract(m.proband, "response"))
  Y.R <- as.matrix(model.extract(m.R, "response"))

  data.fam <- model.matrix(new.fam.formula, outdata.fam)[,-1]
  colnames(data.fam) <- c("fid", "rel.id", "exam.age", "B")

  data.proband <- model.matrix(new.R.formula, outdata.proband)[,-1]
  data.R <- model.matrix(new.R.formula, outdata.R)[,-1]
  colnames(data.proband) <- c("exam.age", "B")
  colnames(data.R) <- c("exam.age", "B")

  Y.S <- data.S <- NULL
  if(!is.null(S.formula)) {
    new.S.formula <- update.formula(S.formula,  paste("~.+", paste(recruit.age.S, birth, sep="+")))
    m.S <- model.frame(new.S.formula, outdata.S)
    Y.S <- as.matrix(model.extract(m.S, "response"))
    data.S <- model.matrix(new.S.formula, outdata.S)[,-1]
    colnames(data.R) <- c("exam.age", "B")
  }

 if(is.null(Age)) Age = c(0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90)
 if(is.null(Cal)) Cal = seq(1890, 2015, by=5)

  lam03$Year.f <- rep(seq(1, length(Cal)), each=length(Age))
  lam03$Age.f <- rep(seq(1, length(Age)), length(Cal))

  nr <- nrow(outdata.R)
  C.B0 <- lapply(1:nr, function(i) data.R[i, "B"]+ Age)
  C <- lapply(1:nr, function(i) sort(c(Cal, C.B0[[i]])))

  Cf.F <- lapply(1:nr, function(i) unique(C[[i]][data.R[i, "B"] <= C[[i]] & C[[i]] <= Y.R[i,2] + data.R[i, "B"] ]))
  Af.F <- lapply(1:nr, function(i) Cf.F[[i]] - data.R[i, "B"])

  R.f <- lapply(1:nr, function(i) sapply(1:length(Cf.F[[i]]), function(k) findInterval(Cf.F[[i]][k], Cal)))
  A.f <- lapply(1:nr, function(i) sapply(1:length(Af.F[[i]]), function(k) findInterval(Af.F[[i]][k], Age)))

  LAM03.R <-  lapply(1:nr, function(j) sapply(1:length(R.f[[j]]), function(i) lam03[lam03$Year.f==R.f[[j]][i] & lam03$Age.f==A.f[[j]][i], ]$rate))
  cut.R <- lapply(1:nr, function(j) Af.F[[j]][-1])

  LAM03.S <- cut.S <- NULL
  if(!is.null(outdata.S)){
  ns <- nrow(outdata.S)
  C.B0 <- lapply(1:nr, function(i) data.S[i, birth]+ Age)
  C <- lapply(1:nr, function(i) sort(c(Cal, C.B0[[i]])))

  Cf.F <- lapply(1:nr, function(i) unique(C[[i]][data.S[i, birth] <= C[[i]] & C[[i]] <= Y.S[i, 1] + data.S[i, birth]]))
  Af.F <- lapply(1:nr, function(i) Cf.F[[i]] - data.S[i, birth])

  R.f <- lapply(1:nr, function(i) sapply(1:length(Cf.F[[i]]), function(k) findInterval(Cf.F[[i]][k], Cal)))
  A.f <- lapply(1:nr, function(i) sapply(1:length(Af.F[[i]]), function(k) findInterval(Af.F[[i]][k], Age)))

  LAM03.S <-  lapply(1:ns, function(j) sapply(1:length(R.f[[j]]), function(i) lam03[lam03$Year.f==R.f[[j]][i] & lam03$Age.f==A.f[[j]][i], ]$rate))
  cut.S <- lapply(1:ns, function(j) Af.F[[j]][-1])
}



Y.fam<- split(as.data.frame(Y.fam), as.factor(data.fam[, "fid"]))
Y.fam <- lapply(1:length(Y.fam), function(i) as.matrix(Y.fam[[i]]))
X.fam<- split(as.data.frame(data.fam), as.factor(data.fam[, "fid"]))
outdata.F0 <- split(outdata.fam, as.factor(outdata.fam[, fam.id]))

if(is.null(init)){
  if(no.death == TRUE) par <- c(rho, log(lam01))
  else par<- c( rho, log(theta), log(lam01))

}else{
  par <- init
}


pairle <- optim(par, pair.logL, Y.fam = Y.fam, X.fam = X.fam,  Y.proband = Y.proband, X.proband = data.proband, Y.R = Y.R, X.R = data.R, Y.S = Y.S,
                outdata.F0 = outdata.F0, outdata.proband = outdata.proband[, first.visit.age.R], cut = cut, LAM03.R = LAM03.R, cut.R = cut.R, LAM03.S = LAM03.S, cut.S = cut.S, Age = Age, Cal = Cal, design = design, full = full, no.death = no.death, method = "BFGS", hessian = TRUE)

parameter.pair <- exp(pairle$par)

if(design == 1){
  if(no.death == TRUE) {

    score_i <- sapply(1:nf, function(i) numDeriv::grad(func=NloglikFD1,  x=pairle$par, outdata_F = outdata.F0[i], outdata_proband=outdata.proband[i,first.visit.age.R],
                                                       Age = Age, Cal = Cal, lam03 = lam03, full = full,
                                                       fgau = gauleg.f,  combn = utils::combn))


  }else {

    score_i <- sapply(1:nf, function(i) numDeriv::grad(func=loglikFD1_pch,  x=pairle$par, outdata_F = outdata.F0[i], outdata_proband=outdata.proband[i,first.visit.age.R],
                                                       Age = Age, Cal = Cal, lam03 = lam03, full = full,
                                                      fgau = gauleg.f, combn = utils::combn))

    }
  }else{
  if(no.death == TRUE) {

    score_i <- sapply(1:nf, function(i) numDeriv::grad(func=NloglikFD2,  x=pairle$par, outdata_F = outdata.F0[i], outdata_proband=outdata.proband[i,first.visit.age.R],
                                                       Age = Age, Cal = Cal, lam03 = lam03,
                                                       fgau = gauleg.f, combn = utils::combn))


  }else {
    score_i <- sapply(1:nf, function(i) numDeriv::grad(func=loglikFD2_pch,  x=pairle$par, Y_F = Y.fam[i], X_F = X.fam[i], Y_proband = t(as.matrix(Y.proband[i,])), X_proband = t(as.matrix(data.proband[i,])),
                                                       Age = Age, Cal = Cal, cut_F = cut, lam03 = lam03, fgau = gauleg.f, combn = utils::combn))

  }
  }
score_r <- matrix(0, nrow=nr, ncol=2)
if(no.death == FALSE) score_r <- sapply(1:nr, function(i) numDeriv::grad(loglikR_pch, x=pairle$par,  cut_F = cut, Y_R = t(as.matrix(Y.R[i,])), X_R = t(as.matrix(data.R[i,])), LAM03R = LAM03.R[i], cutR = cut.R[i], fgau = gauleg.f))
else score_r <- sapply(1:nr, function(i) numDeriv::grad(NloglikR_pch, x=pairle$par,  cut_F = cut, Y_R = t(as.matrix(Y.R[i,])), X_R = t(as.matrix(data.R[i,])), fgau = gauleg.f))


if(!is.null(outdata.S)) {
  if(no.death == TRUE) {
  score_s <- sapply(1:ns, function(i) numDeriv::grad(NloglikS_pch, x=pairle$par,  outdata_S = t(as.matrix(outdata.S[i,])), fgau = gauleg.f))
  }else{
    score_s <- sapply(1:ns, function(i) numDeriv::grad(loglikS_pch, x=pairle$par,  cut_F = cut, Y_S = t(as.matrix(Y.S[i,])) ,LAM03S = LAM03.S[i], cutS = cut.S[i], fgau = gauleg.f))
  }
}
B <- score_i%*%t(score_i)
B <- B + score_r%*%t(score_r)

if(!is.null(outdata.S)){
   B <- B +  score_s%*%t(score_s)
}

A <- pairle$hessian
var <- solve(A)%*%B%*%t(solve(A))
sd.pair <- sqrt(diag(var))

result <- list()
result$estimate <- parameter.pair
result$var <- var
return(result)
}


