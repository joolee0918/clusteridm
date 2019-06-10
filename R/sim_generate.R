#==================================================
# Generate latent illness-death model
#
# Biased sampling i) tracing proband ii) alive non-proband only data observed; m=2, 4, 6 with 1/3 prob
# Last Update : April 09, 2018
#==================================================

#' @importFrom timereg pc.hazard
#' @importFrom copula cCopula claytonCopula iTau
# ----------------------------------------------------------
# Marginal function of Exponential distribution
# ----------------------------------------------------------

marg.f <- function(lam, kap, beta, xx, tt){
  #p <- length(beta)
  #xx <- matrix(xx, ncol=p)
  #beta <- matrix(beta, ncol=1)
  #tmp <- as.numeric(xx %*% beta)
  St <- exp(-lam*tt)
  ft <- lam*exp(-lam*tt)
  Ft <- 1-St
  out <- data.frame(St = St, Ft = Ft, ft = ft)
  return(out)
}


dClayton <- function(u1, u2, rho){
  res <- (1+rho)*(u1*u2)^(-rho-1)*(u1^(-rho) + u2^(-rho)-1)^(-1/rho-2)
  return(res)
}
pClayton <- function(u1, u2, rho){
  res <- (u1^(-rho) + u2^(-rho)-1)^(-1/rho)
  return(res)
}

hf <- function(u, del, k, mi, rho){
  term <- rep(0, mi)
  for(j in 1:mi){
    term[j] <- u[j]^(-rho)
  }
  if(k>1) {
    term1=1
    for(j in 1:(k-1)) {tmp <- (1+j*rho);
    term1 = term1*tmp}
  } else {
    term1=1}

  term2 <- rep(0, ncol=mi)
  for(j in 1:mi){
    term2[j] <- u[j]^(del[j])
  }
  res = term1*(prod(term2))^(-rho-1)*(sum(term)-mi+1)^(-1/rho-k)
  return(res)
}

h1 <- function(u0, u1, rho){


  res = u0^(-rho-1)*(u1^(-rho) + u0^(-rho) -1) ^(-1/rho-1)
  return(res)
}

ff <- function(u1, u2, j1, j2, lam01, n.rho, rho, cut, LAM03, LAM12, T0, A1) {
  dClayton(h1(pexp(T0, rate=lam01, lower.tail=F), pexp(u1, rate=lam01, lower.tail=F), rho=rho) , h1(pexp(T0, rate=lam01, lower.tail=F), pexp(u2, rate=lam01, lower.tail=F), rho=rho),  n.rho)*dClayton(pexp(T0, rate=lam01, lower.tail=F), pexp(u1, rate=lam01, lower.tail=F), rho)*dClayton(pexp(T0, rate=lam01, lower.tail=F), pexp(u2, rate=lam01, lower.tail=F), rho)*dexp(u1, rate=lam01)*dexp(u2, rate=lam01)*ppch(u1, levels=LAM03[[j1+1]], cuts=cut[[j1+1]][-1], lower.tail=F)*ppch(A1[j1], levels=LAM12[[j1+1]], cuts=cut[[j1+1]][-1], lower.tail=F)/ppch(u1, levels=LAM12[[j1+1]], cuts=cut[[j1+1]][-1], lower.tail=F)*ppch(u2, levels=LAM03[[j2+1]], cuts=cut[[j2+1]][-1], lower.tail=F)*ppch(A1[j2], levels=LAM12[[j2+1]], cuts=cut[[j2+1]][-1], lower.tail=F)/ppch(u2, levels=LAM12[[j2+1]], cuts=cut[[j2+1]][-1], lower.tail=F)
}


#==================================================
# Generate data
#==================================================

## Generate family data
#' @export
generate_con <- function(id, mi, lam01, lam12, lam03, ktau, proband, sen, u0) {

  Age = c(0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90)
  Cal = seq(1920, 2011, by=5)

  rho <- iTau(claytonCopula(), ktau)
  C2 <- claytonCopula(rho, dim = mi)

  mulT <- rep(0, mi-1)
  if(sen == T) {
    u <- rgamma(mi-1, shape = 1/0.2, scale=0.2)
  }else{
    u <- rep(1, mi-1)
  }

  multU <- copula::cCopula(c(pexp(proband$X, rate=u0*lam01, lower.tail=F), 1-runif(mi-1, 0, 1)), copula=C2, inverse=T)[-1]

  multT <- qexp(multU, rate=u*lam01, lower.tail=FALSE)
  tmp <- data.frame(fam.id = rep(id, mi), mi=rep(mi, mi), mem=seq(1, mi), rel.id=c(rep(1,2), rep(2, mi-2)))
  tmp$proband <- ifelse(proband$mem==tmp$mem, 1, 0)
  tmp <- tmp[order(tmp$proband, decreasing=T), ]

  if(tmp[tmp$proband==1,]$rel.id==1) {

    B1 <- proband$B + runif(1, 0, 10)
    B2 <- proband$B + runif(mi-2, 20, 30)
    B <- c(B1, B2)

  } else {

    B1 <- proband$B  - runif(2, 20, 30)
    B2 <- proband$B  + runif(mi-3, 0, 10)
    B <- c(B1, B2)
  }


  t = lapply(1:(mi-1), function(i) B[i] + Age)
  C = lapply(1:(mi-1), function(i) sort(c(Cal, t[[i]])))
  Cf = lapply(1:(mi-1), function(i) unique(C[[i]][B[i] <=C[[i]] & C[[i]] <= 2010.5]))
  Af = lapply(1:(mi-1), function(i) Cf[[i]] - B[i])

  R = lapply(1:(mi-1), function(j) sapply(1:length(Cf[[j]]), function(k) findInterval(Cf[[j]][k], Cal)))
  A =  lapply(1:(mi-1), function(j) sapply(1:length(Af[[j]]), function(k) findInterval(Af[[j]][k], Age)))

  LAM03 <- lapply(1:(mi-1), function(j) sapply(1:length(R[[j]]), function(i) lam03[lam03$Year.f==R[[j]][i] & lam03$Age.f==A[[j]][i], ]$rate))
  LAM12 <- lapply(1:(mi-1), function(j) sapply(1:length(R[[j]]), function(i) lam12[lam12$Year.f==R[[j]][i] & lam12$Age.f==A[[j]][i], ]$rate))

  multD = t(sapply(1:(mi-1), function(i) timereg::pc.hazard(cbind(breaks=c(Af[[i]], 2010.5-B[i]) , rate=c(0, u[i]*LAM03[[i]])), rr=1, cum.hazard=F)))


  multX <- rep(0, mi-1)
  multY <- rep(0, mi-1)
  multdel1 <- rep(0, mi-1)
  multdel2 <- rep(0, mi-1)
  multdel3 <- rep(0, mi-1)

  for(i in 1:(mi-1)){
    multX[i] <- min(multT[i], multD[i,]$time)
    multdel1[i] <- multD[i,]$status
    multdel2[i] <- ifelse(multT[i] < multD[i, ]$time, 1,0)
  }
  for(i in 1:(mi-1)){

    if(multdel2[i]==1) {
      tmpy = timereg::pc.hazard(cbind(breaks=c(Af[[i]], 2010.5-B[i]) , rate=c(0, u[i]*LAM12[[i]])), rr=1,  entry=multX[i], cum.hazard=F)
      multY[i] = tmpy$time
      multdel3[i] = tmpy$status
    } else {
      multY[i] = 0
      multdel3[i] = 0
    }
  }
  multX = c(proband$X, multX)
  multY = c(proband$Y, multY)
  B = c(proband$B, B)
  exam.age = c(2010.5 - B)
  del1 = c(proband$del1, multdel1)
  del2 = c(proband$del2, multdel2)
  del3 = c(proband$del3, multdel3)
  u = c(u0, u)

  tmp <- cbind(tmp, B = B, X = multX, Y = multY, exam.age=exam.age, del1=del1, del2=del2, del3=del3, u = u)
  return(tmp)
}


#' @export
generate_mar <- function(id, B0, R0, lam01, lam12, lam03, sen){

  v = runif(1, 0, 1)
  if(sen == T) {
    u <- rgamma(1, shape = 1/0.2, scale=0.2)
  }else{
    u <- 1
  }
  multT <- qexp(v, rate=u*lam01, lower.tail=F)
  Age = c(0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90)
  Cal = seq(1920, 2011, by=5)

  C.B0 = B0 + Age
  C = sort(c(Cal, C.B0))

  Cf.R = unique(C[B0 <=C & C <= R0])
  Af.R = Cf.R - B0

  R = sapply(1:length(Cf.R), function(k) findInterval(Cf.R[k], Cal))
  A = sapply(1:length(Af.R), function(k) findInterval(Af.R[k], Age))

  F0 = 2010.5
  Cf.F = unique(C[B0 <=C & C <= F0])
  Af.F = Cf.F - B0

  R.F = sapply(1:length(Cf.F), function(k) findInterval(Cf.F[k], Cal))
  A.F = sapply(1:length(Af.F), function(k) findInterval(Af.F[k], Age))

  LAM03 <-  sapply(1:length(R), function(i) lam03[lam03$Year.f==R[i] & lam03$Age.f==A[i], ]$rate)
  LAM12 <-  sapply(1:length(R.F), function(i) lam12[lam12$Year.f==R.F[i] & lam12$Age.f==A.F[i], ]$rate)

  multD = timereg::pc.hazard(cbind(breaks=c(Af.R, R0-B0) , rate=c(0, u*LAM03)), rr=1, cum.hazard=F)

  multX <- min(multT, multD$time)
  multdel1 <- multD$status
  multdel2 <- ifelse(multT < multD$time, 1,0)

  if(multdel2==1){
    tmp <-  timereg::pc.hazard(cbind(breaks=c(Af.F, F0-B0) , rate=c(0, u*LAM12)), rr=1,  entry=multX, cum.hazard=F)
    multY <- tmp$time
    multdel3 <- tmp$status
    multdel1 <- 0
  }else {
    multY <- 0; multdel3 <- 0
  }

  res <- data.frame(id=id, B=B0, X=multX, Y=multY, del1=multdel1, del2=multdel2,  del3=multdel3, u = u )
  return(res)
}

## Generage current status data
#' @export
generate_cur <- function(id, B0, lam01, lam12, lam03, sen){

  v = runif(1, 0, 1)
  if(sen == T) {
    u <- rgamma(1, shape = 1/0.2, scale=0.2)
  }else{
    u <- 1
  }

  multT <- qexp(v, rate=u*lam01, lower.tail=F)
  Age = c(0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90)
  Cal = seq(1920, 2011, by=5)
  C.B0 = B0 + Age
  C = sort(c(Cal, C.B0))

  Cf.R = unique(C[B0 <=C & C <= 2000.5])
  Af.R = Cf.R - B0

  R = sapply(1:length(Cf.R), function(k) findInterval(Cf.R[k], Cal))
  A = sapply(1:length(Af.R), function(k) findInterval(Af.R[k], Age))

  LAM03 <-  sapply(1:length(R), function(i) lam03[lam03$Year.f==R[i] & lam03$Age.f==A[i], ]$rate)
  LAM12 <-  sapply(1:length(R), function(i) lam12[lam12$Year.f==R[i] & lam12$Age.f==A[i], ]$rate)

  multD = timereg::pc.hazard(cbind(breaks=c(Af.R, 2000.5-B0) , rate=c(0,u*LAM03)), rr=1, cum.hazard=F)

  multX <- min(multT, multD$time)
  multdel1 <- multD$status
  multdel2 <- ifelse(multT < multD$time, 1,0)

  if(multdel2==1){
    tmp <-  timereg::pc.hazard(cbind(breaks=c(Af.R, 2000.5-B0) , rate=c(0, u*LAM12)), rr=1,  entry=multX, cum.hazard=F)
    multY <- tmp$time
    multdel3 <- tmp$status
  }
  else {multY <- 0; multdel3 <- 0}



  res <- data.frame(id=id, B=B0, X=multX, Y=multY, del1=multdel1, del2=multdel2,  del3=multdel3 )
  return(res)
}

