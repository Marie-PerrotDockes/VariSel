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
#'#' @import R6 Matrix gglasso tidyverse glmnet stabs magrittr viridis stringr
#' @export
mod_group <- R6::R6Class("mod_group",
  inherit = VariSel,
  public = list(
    group = NULL,
    estime = function(lambda = NULL){
      if(!is.null(lambda)) {
        self$mod <- private$tb %>%
         mutate( Lambda = lambda,
          Model = map2(Data, Lambda,
           ~grp_lasso(private$x, .x,
            self$group, lambda = .y)
        ))
      }else{
        self$mod <- private$tb %>%
          mutate(Model = map(Data,
            ~grp_lasso(private$x, .,
              self$group)
              ))
      }
      self$res <- self$mod  %>%
        mutate( Beta = map(Model, ~.$beta),
          Intercept = map(Model, ~.$b0),
          Group = map(Model, ~.$group),
          Lambda = map(Model, ~.$lambda),
          Df = map2(Beta, Group, ~.x %>%
            as.data.frame() %>%
            mutate(group = .y) %>%
            group_by(group) %>%
            summarise_all(sum) %>%
            summarise_all(~sum(. != 0)) %>%
            select(-group))
        )
      super$estime()
    },
    sel_cv = function( s = "lambda.min"){
      self$cv <-  private$tb %>%
        mutate(
          Model = map(Data,
            ~cv_grp_lasso(private$x, ., self$group)
          ),
          Beta = map(Model, ~coef(., s)[-1, ]),
          Beta = map(Beta, function(beta){
            beta %>% as.data.frame() %>%
              rownames_to_column() %>%
              mutate(order = c(match(colnames(private$x), rowname))) %>%
              arrange(order) %>% select(-order) %>%
              filter(. != 0) %>% dplyr::rename(., value = .)
              }
            ))
    },
    sel_stab = function( nb.cores = 7, B = 500, PFER = 1){
      self$stab <-  private$tb %>%
      mutate(
        Model = map(Data,
          ~grp_lasso_st(private$x, .,
          self$group, nb.cores = nb.cores, B = B, PFER = PFER)),
        Frequencies = map(Model, ~.$phat),
        Selected    = map(Model, function(mod){
          t <- rep(FALSE, ncol(private$x))
          t[mod$selected] <- TRUE
          t
        })
    )
  },
    get_coef = function(){
    super$get_coef()
    self$coef <- self$coef %>%
      separate(Marker, sep =self$sepx, into = c ("Reg", "group"))
  }
))

#' Description of the function
#'
#' @export
mod_group_univ <- R6::R6Class("mod_group_univ", inherit = mod_group,
  private = list(r = NULL),
  public = list(
    initialize = function(x, y,
                          Sigma_12inv = diag(1, ncol(as.data.frame(y))),
                          sepx = "\\."){
      super$initialize(x, y,Sigma_12inv = Sigma_12inv, sepx = sepx)
      self$group <- get_group(private$name_x, sep = sepx)
    }
))

#' Description of the function
#' @export
mod_group_multi <- R6::R6Class("mod_group_multi", inherit = mod_group,
  private = list(r = NULL),
  public = list(
    initialize = function(x, y,
                  Sigma_12inv = diag(1, ncol(as.data.frame(y))),
                  univ = FALSE, sepx ="\\." ){
      super$initialize(x, y, univ = univ,
                       Sigma_12inv = Sigma_12inv,
                       sepx = sepx)
    },
    sel_cv = function(s = "lambda.min"){
      super$sel_cv(s = s)
      self$cv <- self$cv %>% mutate(Beta = map(Beta ~ . %>%
        separate(rowname, sep = self$sepy,
                 into = c("Trait", "rowname"), fill = "left")))
    }
))



#' Description of the function
#'
#' @export
mod_group_multi_both <- R6::R6Class("mod_group_multi_both",
  inherit = mod_group_multi,
  public = list(
    initialize = function(x, y, univ = FALSE,
                  Sigma_12inv = diag(1, ncol(as.data.frame(y))), sepx =  "\\."){
      super$initialize(x, y, univ = FALSE,
        Sigma_12inv = diag(1, ncol(as.data.frame(y))),
        sepx = sepx)
      self$group <- get_group_both(private$name_x, sep = sepx, r = private$r)
    }
  ))

#' Description of the function
#'
#' @export
mod_group_multi_regr <- R6::R6Class("mod_group_multi_regr",
  inherit = mod_group_multi,
  public = list(
    initialize = function(x, y, univ = FALSE,
                  Sigma_12inv = diag(1, ncol(as.data.frame(y))), sepx = "\\."){
      super$initialize(x, y, univ = FALSE,
        Sigma_12inv = diag(1, ncol(as.data.frame(y))),
        sepx = sepx)
      self$group <- get_group_marker(private$name_x, sep = sepx, r = private$r)
    }
))
