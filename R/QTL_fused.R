#' Description of the function
#'
#' @param  response a vector response variable
#' @param  regressors a quantitative matrix of regressor
#' @param  group a vector with two levels. (The group of the ANCOVA)
#' @param  a the parameters that indicate how much the coefficients will be fused
#' @param  lambda if the user wants to use it owns values of lambdas
#' @return The coefficients of the fused lasso ANCOVA for the different value of lambda
#' @examples
#' B <- c(1, -1, 1.5, 1.5, rep(0, 6), 2, 0, 2, 0)
#'group <- c(rep('M1', 10), rep('M2', 10))
#'regressors <- matrix(rnorm(6*20), ncol = 6)
#'X  <- model.matrix(~group + group:regressors - 1)
#'y <- X%*%B + rnorm(20)
#'y <- scale(y)
#'mod <- fl2(y, regressors, group)
#'colors <- c(rep("grey",2), rep('green',2),rep('black', 6), rep(c("orange","blue"), 2), 'darkgreen', rep('yellow',3), rep('purple',2))
#'matplot(mod$lambda ,t(mod$beta),type='l',col=colors)
#'#' @import R6 Matrix gglasso tidyverse glmnet stabs magrittr viridis stringr FusedLasso
#' @export
mod_fused <- R6Class("mod_fused",
  inherit = VariSel,
  public = list(
  group = NULL,
  graphe = NULL,
  estime = function( lambda2 = NULL, lambda1  = 0.002,
    ratio = 1e-3, nlambda2 = 100){
    self$mod <-  private$tb %>%
      mutate(Model = map(Data,
        ~fused_lasso(X = private$x, response = ., G = self$graphe,
          lambda2 = lambda2, lambda1  = lambda1,
          ratio = ratio, nlambda2 = nlambda2)
                              ))
    self$res <- self$mod  %>%
      mutate( Beta = map(Model, function(mod){
          a <- mod$beta
          rownames(a) <- colnames(private$x)
          return(a)
      }),
        Intercept = map(Model, ~.$Intercept),
        Lambda1 = map(Model, ~.$lambda1),
        Lambda2 = map(Model, ~.$lambda2),
        Df = map(Beta, ~.x %>%
          as.matrix() %>%
          as.data.frame() %>%
          mutate(group = self$group) %>%
          group_by(group) %>%
          summarise_all(sum) %>%
          summarise_all(~sum(. != 0)) %>%
          select(-group))
        )
  }
))

#' Description of the function
#'
#' @export
mod_fused_univ <- R6Class("mod_fused_univ", inherit = mod_fused,
  private = list(r = NULL),
  public  = list(
    initialize = function(x, y, sep = "\\."){
      super$initialize(x, y)
      self$group <- get_group(private$name_x, sep = sep)
      self$graphe <- getGraphe(self$group)
    }
 )
)

#' Description of the function
#' @export
mod_fused_multi <- R6Class("mod_fused_multi", inherit = mod_fused,
  private = list(r = NULL),
  public  = list(
    initialize = function(x, y,  Sigma_12inv = diag(1, ncol(as.data.frame(y))),
                          univ = FALSE){
      super$initialize(x, y, univ = univ, Sigma_12inv = Sigma_12inv)
    }
  ))



#' Description of the function
#'
#' @export
mod_fused_multi_both <- R6Class("mod_fused_multi_both",
  inherit = mod_fused_multi,
  public = list(
    initialize = function(x, y, univ = FALSE,
                          Sigma_12inv = diag(1, ncol(as.data.frame(y)))){
      super$initialize(x, y, univ = FALSE, Sigma_12inv = Sigma_12inv)
      self$group <- get_group_both(private$name_x, r = private$r)
      self$graphe <- getGraphe(self$group)
    }
 ))

#' Description of the function
#'
#' @export
mod_fused_multi_regr <- R6Class("mod_fused_multi_regr",
  inherit = mod_fused_multi,
  public = list(
    initialize = function(x, y, univ = FALSE,
                          Sigma_12inv = diag(1, ncol(as.data.frame(y)))){
      super$initialize(x, y, univ = FALSE, Sigma_12inv = Sigma_12inv)
      self$group <- get_group_marker(private$name_x, r = private$r)
      self$graphe <- getGraphe(self$group)
    }
))
