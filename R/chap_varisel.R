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
train_VariSel <- function(X, Y, type, sepx, reg, group, lambda = NULL,
                          lambda1=NULL, lambda2 = NULL){
  if(!is.null(X) & !is.null(reg)) stop("either reg and group must be supplied or X (the design matirx already full)")
  if(!is.null(sepx) & !is.null(group)) stop("either sepx or group must be supplied see the help of the function")
  if (is.null(X)){
    if(is.null(group) | is.null(reg)) stop('Error if the design matrix X is not supplied both group and reg must be')
  X <- model.matrix(~reg:group -1)
  sepx <- ":"
  }
  if(!is.null(lambda1)){
    lambda = list(lambda1 = lambda1, lambda2 = lambda2)
  }
   mod <- type_to_varisel(X, Y, type,
                          Sigma_12inv = Sigma_12inv, sep= sepx)
   mod$estime(lambda = lambda)
}


## caret

# tuned_grid
