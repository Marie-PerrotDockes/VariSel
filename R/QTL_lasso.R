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
mod_lasso <- R6::R6Class("mod_lasso",
  inherit = VariSel,
  public = list(
    penalty.factor = NULL,
    initialize = function(x, y, univ = TRUE,
                          Sigma_12inv = diag(1, ncol(as.data.frame(y))),
                          sepy ="__", penalty.factor = NULL){
      super$initialize(x, y, univ,  Sigma_12inv, sepy)
      if (is.null(penalty.factor)){
        self$penalty.factor <- rep(1, ncol(private$x))
      } else {
        self$penalty.factor <- penalty.factor
          if (length(self$penalty.factor) != ncol(private$x))
            stop("weight must be of length ncol(X)")
      }
    },

    estime = function(lambda = NULL, standardize = FALSE, intercept = TRUE){
      if(!is.null(lambda)){
        self$mod <-  private$tb %>%
          mutate( Lambda = lambda,
            Model = purrr::map2(Data, Lambda,
            ~ glmnet(private$x, .x,
              lambda = .y,
              standardize = standardize,
              intercept = intercept,
              penalty.factor = self$penalty.factor)))
      }else{
        self$mod <-  private$tb %>%
          mutate(Model = purrr::map(Data,
              ~ glmnet(private$x, .,
                standardize = standardize,
                intercept = intercept,
                penalty.factor = self$penalty.factor)))
      }
      self$res <- self$mod  %>%
        mutate( Beta = map(Model, ~.$beta %>% as.matrix()),
          Intercept = map(Model, ~.$a0),
          Lambda = map(Model, ~.$lambda),
          Df = map(Beta, ~ as.matrix(.) %>%
            as.tibble() %>%
            summarise_all( ~sum(. != 0)))) %>%
          select(-Model)
      super$estime()
  },

    sel_cv = function(s = "lambda.min"){
      browser()
      if (is.null(private$cv)){
        self$cv <- private$tb
        mutate(Model = purrr::map(Data,
          ~ cv.glmnet(private$x, .,
                      penalty.factor = self$penalty.factor)))
      }
      private$s <- s
      self$cv <- self$cv %>%
        mutate(Beta = map(Model, ~coef(., private$s)[-1, ] %>%
          as.data.frame() %>%
          rownames_to_column() %>%
          dplyr::rename(., value = .) %>%
          filter(value != 0)
        ))
    },
  sel_stabs = function( nb.cores = 7, B = 500, PFER = 1){
    nc <- round(ncol(private$x) / max(group))
    q <- max(5, min(2 * round(nrow(private$x) / log(nc) / 10) * 10, nc))
    my_stab <- partial(stabsel, x = private$x,   q = q, PFER = 1, B = B,
      fitfun = glmnet.lasso,
      args.fitfun = list( penalty.factor = self$penalty.factor),
      mc.cores = nb.cores, mc.preschedule = TRUE)
    self$stab <-  private$tb %>%
      mutate(Model = map(Data,
                         ~my_stab(y = .)),
             Frequencies = map(Model, ~.$phat),
             Selected = map(Model, function(mod){
               t <- rep(FALSE, ncol(private$x))
               t[mod$selected] <- TRUE
               t
             })
      )
  },
  plot_path = function( type ="first", nb = 6){
    if(is.null(self$coef)) self$get_coef()
    res <- self$coef %>%
      mutate(ret  =  paste(Marker,  Trait, sep =" on "))

    sel <- res %>%  pull(ret) %>%
      unique()
    sel <- sel[1:nb]
    res_sel <- res %>% filter(ret %in% sel) %>%
      mutate(col = factor(ret, levels = unique(ret)))
    res_other <- res %>% filter(!ret %in% sel)
    ggplot(data = res_other,aes(x = Lambda, y = value, group = ret))+
      geom_line(color = "lightgray") +
      geom_line(data = res_sel,aes(color =col)) +
      scale_x_log10() +
      labs ( color = " ", y = "value of the coefficients", title = "Regularization Path", x = "Lambda")


  }
))
