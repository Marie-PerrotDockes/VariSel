#' Title
#'
#' @param new_x
#' @param new_group
#' @param lambda
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
mod_fus2_univ <- R6Class("mod_fus2_univ",
  inherit = mod_lasso,
  public = list(
    initialize = function(X, Y, group,  a = 1){
      X  <- as.matrix(X)
      p  <- ncol(X)
      X1 <- model.matrix(~group + group:X - 1)
      X  <- cbind(X1, X)
      b  <- (3 * p - a * p + 2) / (2 * p )
      penalty.factor <- c(0, 0, rep(b, (ncol(X1) - 2)), rep(a, p))
      super$initialize(X, Y, penalty.factor = penalty.factor)
    },
    predict = function(new_x, new_group, lambda = NULL, ...){
      new_x  <- as.matrix(new_x)
      p  <- ncol(new_x)
      new_x1 <- model.matrix(~group + group:new_x - 1)
      new_x  <- cbind(new_x1, new_x)
      super$predict(new_x = new_x, lambda = lambda)
    }
  ))
#' Title
#'
#' @param X
#' @param Y
#' @param a
#'
#' @return
#' @export
#'
#' @examples
mod_fus2_resp <- R6Class("mod_fus2_resp",
 inherit = mod_lasso,
 public = list(
  initialize = function(X, Y, a = 1){
    X <- as.matrix(X)
    q <- ncol(Y)
    if (q != 2) stop("Y must have two columns to use fus2resp")
    if (is.null(colnames(Y))) colnames(Y) <- paste0("rep", 1:q)
    self$trait <- as.list(colnames(Y))
    X <- bdiag(rep(list(X), q)) %>% as.matrix()
    group <- rep(colnames(Y), each = nrow(Y))
    Y <- Y %>% as.matrix() %>%  as.numeric()
    p <- ncol(X)
    X1  <- model.matrix(~group + group:X - 1)
    X <- cbind(X1, X)
    b <- (3 * p - a * p + 2) / (2 * p)
    penalty.factor <- c(0, 0, rep(b, (ncol(X1) - 2)), rep(a, p))
    super$initialize(X, Y, penalty.factor = penalty.factor, univ = TRUE)
  },
  predict = function(new_x, names_y, lambda = NULL, ...){
    if(!is.matrix(new_x)) new_x <- as.matrix(new_x)
    group <- rep(names_y, each = nrow(new_x))
    new_x <- bdiag(rep(list(new_x), 2)) %>% as.matrix()
    new_x1 <- model.matrix(~group + group:new_x - 1)
    new_x  <- cbind(new_x1, new_x)
    super$predict(new_x = new_x, lambda = lambda)
  }
))
