#' @title Ball Correlation Sure Independence Screening For Survival data
#' @description Utilize extension of Ball Correlation in survival to select candidate variables related to survival status.
#' @inheritParams bcorsis
#' @param Y a numeric matirx(first column should be event time, second column should be survival status) or Surv object
#' @param standized allows the user to standardize the covariate
#' @import stats
#' @import utils
#' @import survival
#' @export
#' @examples
#' data(survdat)
#' re <- bcorsis.surv(Y = survdat[, c(1, 2)], X = survdat[, -c(1, 2)], candidate = "large")
#' re$candidate.set
#'
bcorsis.surv <- function(Y, X, candidate = "small", standized = TRUE){

  n <- dim(X)[1]; p <- dim(X)[2]
  n <- as.numeric(n); p <- as.numeric(p)
  ids <- 1:p
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  colnames(X) <- paste0("X", 1:p)
  colnames(Y) <- paste0("Y", 1:ncol(Y))
  if(any(apply(Y,2,anyNA))) {stop("NA appear in matrix Y")}
  if(any(apply(X,2,anyNA))) {stop("NA appear in matrix X")}

  # decide candicate size
  d_logn <- round(n/log(n))
  d <- n
  if(is.integer(candidate)) {
    final_d <- candidate
  } else {
    if(candidate == "small"){
      final_d <- d_logn
    } else {
      final_d <- d
    }
  }

  # prepare for screening
  time <- Y[,1]
  delta <- Y[,2]
  ord.t <- sort(time)
  ix <- order(time)
  ord.delta <- delta[ix]
  xo <- X[ix,]
  if(standized) {
    xo <- apply(xo, 2, scale)
  }

  # BCor Screening(survival)
  fitc <- survfit(Surv(time,1-delta)~1)
  Sc <- fitc$surv
  if(length(unique(ord.t)) != n) {
    rep_num <- as.data.frame(table(ord.t))[, "Freq"]
    Sc <- mapply(function(x, y) {
      rep(x, y)
    }, Sc, rep_num, SIMPLIFY = FALSE)
    Sc <- unlist(Sc)
  }
  # num_non_zero = n - sum()
  rcory_result <- apply(xo,2,function(x){
    SRCT(x = x,t = ord.t,delta = ord.delta,Sc = Sc,n = n)
  })
  rcory_result <- unlist(rcory_result)
  max_ids <- order(rcory_result,decreasing=T)
  chooseids <- max_ids[1:final_d]
  Xhavepickout <- ids[chooseids]

  return(list(candidate.set=Xhavepickout,
              candidate.data=list(X=X[,Xhavepickout],Y=Y),
              candidate.size=length(Xhavepickout)))
}

#' @title Ball Correlation in survival.
#' @param x ordered covariate
#' @param t ordered survival event time
#' @param delta ordered survival event status
#' @param Sc Survfit object
#' @param n Sample size
#' @useDynLib bcorsis
#' @noRd
#'
SRCT <- function(x, t, delta, Sc, n){
  RCT <- numeric(1)
  RCT<-.C("SRCT", as.double(t(x)), as.double(t(t)), as.double(t(delta)),
          as.double(t(Sc)), as.integer(n), RC=as.double(RCT))
  return(RCT$RC)
}


#' @title Ball Correlation in survival analysis
#' @description Calculate Ball Correlation Statistic in survival
#' @param time Time
#' @param status Status
#' @param x single explanatory variable
#' @return Survival Ball Correlation
#' @import stats
#' @import utils
#' @import survival
#' @export
#' @examples
#' data(survdat)
#' bcor.surv(time = survdat[["time"]], status = survdat[["status"]], x = survdat[,3])
#'
bcor.surv <- function(time, status, x) {
  n <- length(time)
  ord.t <- sort(time)
  ix <- order(time)
  ord.delta <- status[ix]
  x <- scale(x[ix])
  #
  fitc <- survfit(Surv(time,1-status)~1)
  Sc <- fitc$surv
  if(length(unique(ord.t)) != n) {
    rep_num <- as.data.frame(table(ord.t))[, "Freq"]
    Sc <- mapply(function(x, y) {
      rep(x, y)
    }, Sc, rep_num, SIMPLIFY = FALSE)
    Sc <- unlist(Sc)
  }
  #
  SRCT(x = x, t = ord.t, delta = ord.delta, Sc = Sc, n = n)
}



