#' @title Ball Correlation
#' @param x a numeric matirx contain explanatory variables
#' @param y a numeric matirx or vector
#' @param weight when weight = TRUE, Weighted Ball Covariance is used instead of Standard Ball Covariance(Default weight = FALSE)
#' @param R the number of replication. When R = 0, the function return Ball Covariance or Weighted Ball Covariance.
#' @param seed the random seed
#' @noRd
#' @useDynLib bcorsis
#'
BI <- function(x, y, weight = FALSE, R = 0, seed = 2015){
	x <- as.matrix(x)
	y <- as.matrix(y)
	dim_x <- dim(x)
	dim_y <- dim(y)
	n <- dim_x[1]
	p <- dim_x[2]
	q <- dim(y)[2]
	RCT <- numeric(1)
	Dx <- numeric(n*n)
	Dy <- numeric(n*n)
	dst <- TRUE
	Dx <- .C("distance", as.double(t(x)), as.double(t(Dx)), as.integer(n),as.integer(p))
	x <- matrix(Dx[[2]],n,n)
	Dy <- .C("distance", as.double(t(y)), as.double(t(Dy)), as.integer(n),as.integer(q))
	y <- matrix(Dy[[2]],n,n)
	RCT<-.C("BI", as.double(t(x)), as.double(t(y)), as.integer(n), as.integer(p), as.integer(q), as.integer(dst), HRC=as.double(RCT), as.integer(seed), as.integer(R), as.integer(weight))
	return(RCT$HRC)
}

#' @title Ball Correlation for univariate random variable
#' @inheritParams BI
#' @param x a numeric vector
#' @param y a numeric vector
#' @noRd
#' @useDynLib bcorsis
#'
UBI <- function(x, y, R = 0, weight = FALSE, seed = 2015){
  if(is.vector(x) & is.vector(y)) {
    n <- length(x)
    RCT <- numeric(1)
    RCT <- .C("UBI", as.double(x), as.double(y), as.integer(n), HRC=as.double(RCT), as.integer(seed), as.integer(R), as.integer(weight))
    return(RCT$HRC)
  }
  else {
    stop("x and y must be vector!")
  }
}


#' Standardized Ball Correlation for Multivariate data
#'
#' Calculate Ball Correlation statistic for Multivariate data
#'
#' @aliases bcor
#' @param x a numeric matirx
#' @param y a numeric matirx
#' @param fast Fast version for univariate case
#' @export
#' @examples
#' n <- 100
#' x <- rnorm(n)
#' y <- rnorm(n)
#' bcor(x, y)
#'
bcor <- function(x, y, fast = FALSE){
  if(fast) {
    return(sqrt(UBI(x = y, y = x, R = 0)/sqrt(UBI(x = x, y = x, R = 0)*UBI(x = y, y = y, R = 0))))
  } else {
    return(sqrt(BI(x = y, y = x, R = 0)/sqrt(BI(x = x, y = x, R = 0)*BI(x = y, y = y, R = 0))))
  }
}



