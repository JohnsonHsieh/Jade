
#
#
###########################################
#' Estimating detected species relative abundance
#' 
#' \code{DetAbu} Estimating detected species relative abundance
#' @param x a vector of species abundance frequency
#' @param zero reserve zero frequency or not. Default is \code{FALSE}.
#' @return a numerical vector
DetAbu <- function(x, zero=FALSE){
  x <- unlist(x)
  n <- sum(x)  
  f1 <- sum(x==1)
  f2 <- sum(x==2)
  f3 <- sum(x==3)
  
  if(f1>0 & f2>0){
    A1 <- f1 / n * ((n-1)*f1 / ((n-1)*f1 + 2*f2))
  }else if(f1>1 & f2==0){
    A1 <- (f1-1) / n * ((n-1)*f1 / ((n-1)*f1 + 2))
  }else{
    return(x/n)
  }

  if(f2>0 & f3>0){
    A2 <- f2 / choose(n, 2) * ((n-2)*f2 / ((n-2)*f2 + 3*f3))^2  
  }else if(f2>1 & f3==0){
    A2 <- (f2-1) / choose(n, 2) * ((n-2)*f2 / ((n-2)*f2 + 3))^2  
  }else{
    A2 <- 0
  }
  
  if(zero==FALSE) x <- x[x>0]
  q.solve <- function(q){
    e <- A1 / sum(x/n*exp(-q*x))
    out <- sum((x/n * (1 - e * exp(-q*x)))^2) - sum(choose(x,2)/choose(n,2)) + A2
    abs(out)
  }
  #q <- tryCatch(uniroot(q.solve, lower=0, upper=1)$root, error = function(e) {1})
  q <- tryCatch(optimize(q.solve, c(0,1))$min, error = function(e) {1})
  e <- A1 / sum(x/n*exp(-q*x))
  o <- x/n * (1 - e * exp(-q*x))
  o
}

#
#
###########################################
#' Estimating undetected species relative abundance
#' 
#' \code{UndAbu} Estimating undetected species relative abundance
#' @param x a vector of species abundance frequency
#' @return a numerical vector
UndAbu <- function(x){
  x <- unlist(x)
  n <- sum(x)  
  f1 <- sum(x==1)
  f2 <- sum(x==2)
  f3 <- sum(x==3)
  f4 <- max(sum(x == 4), 1)
  f0.hat <- ceiling(ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2))  #estimation of unseen species via Chao1
  if(f0.hat < 0.5){
    return(NULL)
  } 
 
  if(f1>0 & f2>0){
    A1 <- f1 / n * ((n-1)*f1 / ((n-1)*f1 + 2*f2))    
  }else if(f1>1 & f2==0){
      A1 <- (f1-1) / n * ((n-1)*f1 / ((n-1)*f1 + 2))
  }else{
    return(NULL)
  }
  
  if(f2>0 & f3>0){
    A2 <- f2 / choose(n, 2) * ((n-2)*f2 / ((n-2)*f2 + 3*f3))^2  
  }else if(f2>1 & f3==0){
    A2 <- (f2-1) / choose(n, 2) * ((n-2)*f2 / ((n-2)*f2 + 3))^2  
  }else{
    A2 <- 0
  }

  
  R <- ifelse(A2>0, A1^2/A2, 1)
  j <- 1:f0.hat
  f.solve <- function(x){ 
    out <- sum(x^j)^2 / sum((x^j)^2) - R
    abs(out)
  }
  b <-  tryCatch(optimize(f.solve, lower=(R-1)/(R+1), upper=1, tol=1e-5)$min, error = function(e) {(R-1)/(R+1)})
  a <- A1 / sum(b^j)
  p <- a * b^j
  if(f0.hat == 1) p <- A1
  p
}


#
#
###########################################
#' Estimating detected species incidence probability
#' 
#' \code{DetInc} Estimating detected species incidence probability
#' @param x a vector of species incidence frequency
#' @param zero reserve zero frequency or not. Default is \code{FALSE}.
#' @return a numerical vector
DetInc <- function(y, zero=FALSE){
  y <- unlist(y)
  nT <- max(y)  
  y <- y[-1]
  Q1 <- sum(y==1)
  Q2 <- sum(y==2)
  Q3 <- sum(y==3)
  
  if(Q1>0 & Q2>0){
    A1 <- Q1 / nT * ((nT-1)*Q1 / ((nT-1)*Q1 + 2*Q2))
  }else if(Q1>1 & Q2==0){
    A1 <- (Q1-1) / nT * ((nT-1)*Q1 / ((nT-1)*Q1 + 2))
  }else{
    return(y/nT)
  }
  
  if(Q2>0 & Q3>0){
    A2 <- Q2 / choose(nT, 2) * ((nT-2)*Q2 / ((nT-2)*Q2 + 3*Q3))^2
  }else if(Q2>1 & Q3==0){
    A2 <- (Q2-1) / choose(nT, 2) * ((nT-2)*Q2 / ((nT-2)*Q2 + 3))^2
  }else{
    A2 <- 0
  }
  
  
  if(zero==FALSE) y <- y[y>0]
  q.solve <- function(q){
    e <- A1 / sum(y/T*exp(-q*y))
    out <- sum((y/nT * (1 - e * exp(-q*y)))^2) - sum(choose(y,2)/choose(nT,2)) + A2
    abs(out)
  }
  #q <- tryCatch(uniroot(q.solve, lower=0, upper=1)$root, error = function(e) {1})
  q <- tryCatch(optimize(q.solve, c(0,1))$min, error = function(e) {1})
  e <- A1 / sum(y/nT*exp(-q*y))
  o <- y/nT * (1 - e * exp(-q*y))
  o
}


#
#
###########################################
#' Estimating undetected species incidence probability
#' 
#' \code{UndInc} Estimating undetected species incidence probability
#' @param x a vector of species incidence frequency.
#' @return a numerical vector
UndInc <- function(y){
  y <- unlist(y)
  nT <- max(y)  
  y <- y[-1]
  Q1 <- sum(y==1)
  Q2 <- sum(y==2)
  Q3 <- sum(y==3)
  Q4 <- max(sum(y == 4), 1)
  Q0.hat <- ceiling(ifelse(Q2 == 0, (nT - 1) / nT * Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2))  #estimation of unseen species via Chao2
  if(Q0.hat < 0.5){
    return(NULL)
  } 
  
  if(Q1>0 & Q2>0){
    A1 <- Q1 / nT * ((nT-1)*Q1 / ((nT-1)*Q1 + 2*Q2))
  }else if(Q1>1 & Q2==0){
    A1 <- (Q1-1) / nT * ((nT-1)*Q1 / ((nT-1)*Q1 + 2))
  }else{
    return(NULL)
  }
  
  if(Q2>0 & Q3>0){
    A2 <- Q2 / choose(nT, 2) * ((nT-2)*Q2 / ((nT-2)*Q2 + 3*Q3))^2
  }else if(Q2>1 & Q3==0){
    A2 <- (Q2-1) / choose(nT, 2) * ((nT-2)*Q2 / ((nT-2)*Q2 + 3))^2
  }else{
    A2 <- 0
  }

  R <- ifelse(A2>0, A1^2/A2, 1)
  j <- 1:Q0.hat
  f.solve <- function(x){ 
    out <- sum(x^j)^2 / sum((x^j)^2) - R
    abs(out)
  }
  b <-  tryCatch(optimize(f.solve, lower=(R-1)/(R+1), upper=1, tol=1e-5)$min, error = function(e) {(R-1)/(R+1)})
  a <- A1 / sum(b^j)
  p <- a * b^j
  if(Q0.hat ==1) p <- A1
  p
}



#
#
###########################################
#' Estimating species rank dstribution for abundance/incidence based data
#' 
#' \code{SpecDist} Estimating species rank dstribution for abundance/incidence based data
#' \code{datatype} the data type of input data. That is individual-based abundance data (\code{datatype = "abundance"}) or sample-based incidence data (\code{datatype = "incidence"}).
#' @param x a a vector of species abundance or incidence frequency. If \code{datatype = "incidence"}, then the input format of first entry should be total number of sampling units, and followed by species incidence frequency.
#' @return a \code{data.frame} object of RAD/RID
#' 
SpecDist <- function(x, datatype="abundance"){
  TYPE <- c("abundance", "incidence")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  
  if(datatype=="abundance"){  
    obs <- DetAbu(x, zero=FALSE)
    und <- UndAbu(x)
    
  }else if(datatype=="incidence"){
    obs <- DetInc(x, zero=FALSE)
    und <- UndInc(x)
  }
  if(is.null(und)){
    out <- data.frame("probability"=obs, "method"="detected")
  }else{
    out <-  rbind(data.frame("probability"=obs, "method"="detected"),
                  data.frame("probability"=und, "method"="undetected"))
  }       
  out[sort(-out[,1]),]
}


#
#
###########################################
# Examples
#


data(ExAbun)
SpecDist(ExAbun, "abundance")


data(ExInci)
SpecDist(ExInci, "abundance")
