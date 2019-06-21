#' Title
#'
#' @param X
#' @param Y
#' @param type
#' @param sepx
#' @param reg
#' @param group
#' @param lambda
#' @param lambda1
#' @param lambda2
#'
#' @return
#' @export
#'
#' @examples
train_VariSel <- function(X =NULL, Y, type, sepx = NULL, regressors, group, lambda = NULL,
                          lambda1=NULL, lambda2 = NULL, Sigma_12inv = NULL, type_S12_inv ="fixed"){

  if(!is.null(X) & !is.null(regressors)) stop("either regressors and group must be supplied or X (the design matirx already full)")
  if(!is.null(sepx) & !is.null(group)) stop("either sepx or group must be supplied see the help of the function")
  if (is.null(X)){
    if(is.null(group) | is.null(regressors)) stop('Error if the design matrix X is not supplied both group and regressors must be')
  X <- model.matrix(~regressors:group -1)
  colnmaes(X) <- gsub('regressors', "", gsub("group", "", colnames(X)))
  sepx <- ":"
  }
  if(grepl("fused", type) & !is.null(lambda1)){
    lambda = list(lambda1 = lambda1, lambda2 = lambda2)
  }
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  if(ncol(Y) !=1 & is.null(Sigma_12inv)) Sigma_12inv <- get_S12inv(Y, X, type)
   mod <- type_to_varisel(X, Y, type,
                          Sigma_12inv = Sigma_12inv, sep= sepx)
   mod$estime(lambda = lambda)
   return(mod)
}


## caret

#' Title
#'
#' @param Y
#' @param X
#' @param type
#'
#' @return
#' @export
#'
#' @examples
get_S12inv <- function(Y, X, type){

    res<- lm(Y~X)$residuals
  return(solve(chol(cor(res))))
}

# tuned_grid
