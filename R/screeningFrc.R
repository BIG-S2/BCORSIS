#' @title Ball Correlation Sure Independence Screening
#' @description Ball Correlation Sure Independence Screening For Regression and Classification task(Include: BcorSIS, BCor-IISIS, Bcor-Interaction, Bcor-pvalue method)
#' @param Y a numeric matirx included univariate or multivariate response variable
#' @param X a numeric matirx included explanatory variables
#' @param candidate size of candidate set, only available to iterative BCor-SIS(I-BCor-SIS).
#' if candidate = "large", n variables are selected where n equal to sample size.
#' if candidate = "small", n/log(n) variables are selected.
#' if candidata is integer, candidata variables are selected.
#' @param method Method for screening procedure, include: "bcorsis", "bcorsis-pvalue",
#' "ibcorsis-lm", "ibcorsis-gam" and "interaction".
#' "bcorsis" means Standard Sure Independence Screening Procedure;
#' Another way for bcorsis procedure with p-value, is "bcorsis-pvalue".
#' "ibcorsis-lm" and "ibcorsis-gam" carry out Iterative BCor-SIS procedure with ordinary linear
#' regression and generalized additive models, respectively.
#' Option "interaction" is design for detect Variables with potential Interaction.
#' @param parms Parameters list only available for Iterative BCor-SIS. d1 is the
#' size of initial set, d2 is the variable set size added in each iteration.
#' df is degree freedom of basis in generalized additive models, so it only
#' play a role when method = "ibcorsis-gam".
#' @param R the number of replication. A argument only available when method == "bcorsis-pvalue"
#' @param parallel If parallel = TRUE, snowfall parallel framework will be used in features screening.
#' @param ncore If ncore = x, we use x threads to screening features.
#' @import stats
#' @import utils
#' @import gam
#' @import snowfall
#' @export
#' @examples
#' ###############  Standard Method  ###############
#' data(regdat)
#' X <- regdat$X
#' Y <- regdat$Y
#' re <- bcorsis(X = X, Y = Y, candidate = "large")
#' re$candidate.set
#' ###############  Iteration Method  ###############
#' data(interactiondat)
#' X <- interactiondat$X
#' Y <- interactiondat$Y
#' re<- bcorsis(X = X, Y = Y, candidate = "large", method = "interaction")
#' re$candidate.set
#' ###############  BCorSIS procedure with p-value  ###############
#' data(regdat)
#' X <- regdat$X
#' Y <- regdat$Y
#' re <- bcorsis(X = X, Y = Y, candidate = "large", method = "bcorsis-pvalue", R = 10)
#' re$candidate.set
#' ###############  Iterative Method  ###############
#' data(regdat)
#' re <- bcorsis(X = X, Y = Y, candidate = "large",
#'             method = "ibcorsis-gam", parms = list(d1 = 5, d2 = 5, df = 4))
#' re$candidate.set
#'
bcorsis <- function(Y, X, candidate = c("large"),
                    method = "bcorsis",
                    parms = list(d1 = 5, d2 = 5, df = 3),
                    R = 99, parallel = FALSE, ncore = 2)
{
  if(parallel) {
    sfInit(parallel = TRUE, cpus = ncore)
    sfLibrary(bcorsis)
  }
  n <- dim(X)[1]
  p <- dim(X)[2]
  ids <- 1:p
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  colnames(X) = paste0("X", 1:p)
  colnames(Y) = paste0("Y", 1:ncol(Y))
  if(any(apply(Y, 2, anyNA))) {stop("NA appear in matrix Y")}
  if(any(apply(X, 2, anyNA))) {stop("NA appear in matrix X")}

  # decide candicate size
  d_logn <- round(n/log(n))
  d <- n
  d1 <- parms$d1
  d2 <- parms$d2
  df <- parms$df
  if(is.integer(candidate)) {
     final_d <- candidate
  } else {
    if(candidate == "small"){
      final_d <- d_logn
    } else {
      final_d <- d
    }
  }

  rcory_result <- c()

  if(method != "interaction") {
    if(method == "bcorsis-pvalue") {
      rcory_result = lapply(as.list(1:length(ids)),function(id, x = X, y = Y){
        xo <- x[,id]
        BI(y, xo, R = R)
      })
      rcory_result <- 1 - unlist(rcory_result)
    } else {
      if(ncol(Y) > 1) {
        if(parallel) {
          sfExport(list = c("X", "Y"), local = TRUE)
          rcory_result <- sfLapply(as.list(1:length(ids)), function(id, x = X, y = Y){
            xo <- x[,id]
            bcor(y, xo)
          })
        } else {
          rcory_result <- lapply(as.list(1:length(ids)),function(id,x = X, y = Y){
            xo <- x[,id]
            bcor(y, xo)
          })
        }
      } else {
        rcory_result <- lapply(as.list(1:length(ids)),function(id, x = X, y = Y[,1]){
          xo <- x[,id]
          bcor(y, xo, fast = TRUE)
        })
      }
      rcory_result <- unlist(rcory_result)
    }
    max_ids <- order(rcory_result, decreasing = T)

    method_sub <- head(unlist(strsplit(method, "-")), n = 1)
    method_sub1 <- tail(unlist(strsplit(method,"-")), n = 1)
    if(method_sub == "ibcorsis"){
      chooseids <- max_ids[1:d1]
      Xhavepickout <- ids[chooseids]
      Xlastpickout <- ids[chooseids]  # flag the pick out variables, used to remove the effect of selected variables

      ids <- ids[-chooseids]
      #iteration round
      if(method_sub1 == 'lm'){
        while(length(Xhavepickout) < final_d)
        {
          # lm fit for X
          Xnew <- lm(X[,ids]~X[,Xhavepickout])$resid
          # lm fit for Y
          Y <- lm(Y~X[,Xlastpickout])$resid

          # BCor-screening
          if(parallel) {
            sfExport(list = c("Xnew", "Y"))
            rcory_result <- sfLapply(as.list(1:length(ids)),function(id,x=Xnew,y=Y){
              xo <- x[,id]
              bcor(y, xo)
            })
          } else {
            rcory_result <- lapply(as.list(1:length(ids)),function(id,x=Xnew,y=Y){
              xo <- x[,id]
              bcor(y,xo)
            })
          }
          rcory_result <- unlist(rcory_result)
          max_ids <- order(rcory_result,decreasing=T)

          # prepare for next iteration
          chooseids <- max_ids[1:d2]
          Xhavepickout <- c(Xhavepickout,ids[chooseids])
          Xlastpickout <- ids[chooseids]
          ids <- ids[-chooseids]
        }
      }
      if(method_sub1 == 'gam'){
        while(length(Xhavepickout) < final_d)
        {
          # gam fit for X
          lastpickout_formula <- paste0(' + s(',colnames(X)[Xlastpickout],collapse = paste0(",df = ",df,")"))
          lastpickout_formula <- paste0(lastpickout_formula,paste0(",df = ",df,")"),collapse = "")
          lastpickout_dat <- X[,Xlastpickout]
          Xnew <- sapply(ids,function(x){
            formula_one <- paste0(colnames(X)[x],"~",lastpickout_formula)
            formula_one <- as.formula(formula_one)
            dat <- as.data.frame(cbind(X[,x],lastpickout_dat))
            colnames(dat)[1] <- colnames(X)[x]
            # dat <- as.data.frame(dat)
            # colnames(dat) <- paste0("X",c(x,Xhavepickout))
            gam(formula_one,data = dat)$residuals
          })

          # gam fit for Y
          dat <- data.frame(Y,lastpickout_dat)
          names(dat)[1] <- c("Y")
          formula_Y <- as.formula(paste("Y~", lastpickout_formula))
          Y <- gam(formula = formula_Y,data = dat)$residuals

          # BCor-screening
          rcory_result <- lapply(as.list(1:length(ids)),function(id, x = Xnew, y = Y){
            xo <- x[,id]
            bcor(y, xo)
          })
          rcory_result <- unlist(rcory_result)
          max_ids <- order(rcory_result,decreasing=T)

          # prepare for next iteration
          chooseids <- max_ids[1:d2] #
          Xhavepickout <- c(Xhavepickout,ids[chooseids])
          Xlastpickout <- ids[chooseids]
          ids <- ids[-chooseids]
        }
      }
    } else {
      chooseids <- max_ids[1:final_d]
      Xhavepickout <- ids[chooseids]
    }
  } else {
    bcorValue <- apply(X, 2, bcor, y = Y)
    bcor2Value <- apply((X)^2, 2, bcor, y = Y)
    max_ids <- order(bcorValue, decreasing = TRUE)
    chooseids <- max_ids[1:final_d]
    Xhavepickout <- ids[chooseids]
    max_ids <- order(bcor2Value, decreasing = TRUE)
    chooseids <- max_ids[1:final_d]
    Xhavepickout <- unique(c(ids[chooseids], Xhavepickout))
  }

  return(list(candidate.set = Xhavepickout,
              candidate.data = list(X=X[,Xhavepickout],Y=Y),
              candidate.size = length(Xhavepickout)))
}


#' @param n Sample Size
#' @param p Dimension of Variables
#' @param rho Correlatin between Each Pairs of Variables
#' @param seed Random Seed
#' @noRd
#' @return A List contain Y(response) and X(covariate)
simulation.regression <- function(n = 100, p = 200, rho = 0, seed = 2015) {
  set.seed(seed)
  n <- n
  p <- p
  error <- rnorm(n,0,1)
  vsubset <- c(1,2,3,4)
  if(rho==0) {
    X <- matrix(data = rnorm(n*p), nrow = n, ncol = p)
  } else {

  }
  beta <- matrix(c(5,5,5,-15),4,1)
  Y <- as.vector(cbind(X[,1],X[,2],X[,3],X[,4])%*%beta) + error
  return(list("Y" = Y, "X" = X))
}



#' @param n Sample Size
#' @param p Dimension of Variables
#' @param rho Correlatin between Each Pairs of Variables
#' @param seed Random Seed
#' @noRd
#' @return A List contain Y(response) and X(covariate)
#'
simulation.interaction <- function(n = 100, p = 100, rho = 0, seed = 2015)
{
  set.seed(seed)
  error <- rnorm(n)
  if(rho==0) {
    X <- matrix(data = rnorm(n*p), nrow = n, ncol = p)
  } else {

  }
  Y <- X[,1]*X[,2] + error
  return(list("Y" = Y, "X" = X))
}
