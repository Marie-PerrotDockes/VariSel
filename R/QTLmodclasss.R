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
#' @import R6 Matrix gglasso tidyverse glmnet stabs magrittr viridis stringr
#' @export
VariSel <- R6Class(
  "VariSel",
  private = list(
    x = NULL,
    y = NULL,
    intercept = NULL,
    r = NULL,
    Data = NULL,
    tresh = NULL,
    df = NULL,
    s = NULL,
    n = NULL,
    name_x =  NULL,
    tb = NULL
  ),
  public = list(
    res = NULL,
    sepy = NULL,
    pe = NULL,
    mod = NULL,
    cv = NULL,
    stab = NULL,
    univ = NULL,
    name_y = NULL,
    trait = NULL,
    initialize = function(x,
                          y,
                          univ = TRUE,
                          Sigma_12inv = diag(1, ncol(as.data.frame(y))),
                          sepy = "__") {
      self$univ <- univ
      private$x <- as.matrix(x)
      private$y <- as.data.frame(y)
      if (nrow(private$x) != nrow(private$y))
        stop("X and Y must have the same number of rows!")
      private$r <- ncol(private$y)
      private$n <- nrow(private$x)
      private$name_x <- colnames(private$x)
      self$sepy <- sepy
      if (is.null(self$trait))  self$trait <- as.list(colnames(private$y))
      if (!univ) {
        private$y <-
          as.numeric(as.matrix(private$y) %*% Sigma_12inv) %>%
          as.data.frame()
        colnames(private$y) <- "allresponses"
        Eg <- expand.grid(colnames(x), colnames(y))
        private$x <-
          kronecker(Matrix::t(Sigma_12inv), private$x)
        colnames(private$x) <-
          paste(Eg[, 2], Eg[, 1], sep = sepy)
      }
      self$name_y <- colnames(private$y)
      private$tb <- private$y %>%
        gather(key = "Trait") %>%
        group_by(Trait) %>%
        nest(.key = Data) %>%
        mutate(Data = map(Data, ~ .$value))
    },

    estime = function() {
      self$res <- self$res %>%
        mutate(
          Yhat = map2(
            Beta,
            Intercept,
            ~ cbind(1, private$x) %*%
              rbind(.y, as.matrix(.x))
          ),
          Ehat = map2(Data, Yhat, ~ .x - .y),
          Ehat = map(Ehat, ~ as.data.frame(.)),
          MSE = map(Ehat, ~ summarise_all(., ~ sum(. ^ 2))),
          # BIC = map2(MSE, Df, ~ private$n * log(.x / private$n) +
          #              log(private$n) * .y)
          BIC = pmap(list(MSE, Df, Yhat), function(x ,y, z)
            log(x /length(z)) + log(length(z)) * y)
        )
    },

    predict = function(new_x, lambda = NULL, ...) {
      if (is.null(self$mod))
        self$estime(lambda = lambda)
      if (!self$univ)
        new_x <- bdiag(rep(list(as.matrix(new_x)), private$r))
      self$res <- self$res %>%
        mutate(New_pred = map2(
          Beta,
          Intercept,
          ~ cbind(1, as.matrix(new_x)) %*%
            rbind(.y, as.matrix(.x))
        ))

    },
    get_beta = function(){
      if(is.null(self$res$Beta)){self$estime()}
      do.call(rbind.data.frame, as.matrrix(self$res$Beta)) %>%
        rownames_to_column() %>%
        separate(rowname, sep =self$sepy,
                 into = c("Trait","Marker"))
    },
    plot_error = function(print = TRUE) {
      if (is.null(self$res$MSE))
        self$predict()
      self$res$MSE %>%
        map(function(x) {
          x %>% gather() %>%
            mutate(key = as.numeric(gsub("s", "", key))) %>%
            ggplot(aes(x = key, y = value)) +
            geom_point() +
            geom_line() + theme_bw()
        }) %>%
        gridExtra::grid.arrange(grobs = .)
    },
    plot_BIC = function(print = TRUE) {
      if (is.null(self$res$BIC))
        self$predict()
      self$res$BIC %>%
        map(function(x) {
          x %>% gather() %>%
            mutate(key = as.numeric(gsub("s", "", key))) %>%
            ggplot(aes(x = key, y = value)) +
            geom_point() +
            geom_line() + theme_bw()
        }) %>%
        gridExtra::grid.arrange(grobs = .)
    },
    plot_coef = function(tresh, sel = colnames(private$y)) {
      if (is.null(self$mod))
        self$estime()
      if (is.null(private$Data) ||
          private$tresh != tresh) {
        private$Data <- self$res$Beta %>% map( ~ .[, tresh]) %>%
          set_names(colnames(private$y)) %>%
          bind_cols(id = colnames(private$x))  %>%
          set_names(c(colnames(private$y), "id"))
        private$tresh <- tresh
      }
      private$Data %>%
        filter_at(vars(sel), any_vars(. != 0)) %>%
        gather(key, value,-id) %>%
        ggplot(aes(x = id, y = key, fill = value)) +
        scale_fill_viridis() +
        geom_tile() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90)) + ylab("") + xlab("")
    },
    ROC = function(b) {
      self$res$Beta %>%
        map(function(x) {
          x  %>%
            as.matrix() %>%
            as.data.frame() %>%
            gather() %>%
            group_by(key) %>%
            summarise_all(funs(TPR(., b = b), FPR(., b = b)))
        })
    },
    plot_anime = function(name_pos = NULL, num_max = 30) {
      if (is.null(name_pos)) {
        name_pos <- 1:ncol(private$x)
      }
      if (is.null(self$res))
        self$estime()
      self$res %>%
        mutate(
          Beta = map(
            Beta,
            ~ rowid_to_column(.) %>%
              gather(key, value,-rowid) %>% as.tibble()
          ),
          Beta = map2(
            Beta,
            Trait,
            ~ mutate(
              .x,
              key = as.numeric(gsub("s", "", key)),
              Trait = .y,
              pos = name_pos[rowid]
            )
          )
        ) %>%
        select(Beta) %>%
        summarise_all( ~ list(bind_rows(.))) %>%
        extract2(1) %>% extract2(1) %>%
        filter(key < num_max & value != 0) %>%
        ggplot(aes(x = pos, y = Trait, color = value)) + geom_point() +
        transition_time(key) + theme_bw() +
        scale_color_gradientn(colours = brewer.pal(10, "Spectral"))
    },
    plot_cv = function(sel = colnames(private$y), s = "lambda.min") {
      if (is.null(self$cv) || private$s != s) {
        self$sel_cv(s)
        private$s <- s
      }
      self$cv %>% pull(Beta) %>% map_dfr( ~ .) %>%
        ggplot(aes(x = Trait, y = Marker, fill = value)) + geom_tile() +
        scale_fill_viridis() +
        theme_bw() + theme(axis.text.x = element_text(angle = 90))
    }
  )
)
