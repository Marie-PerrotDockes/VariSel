
#' Title
#'
#' @param graphe
#' @param lambda2
#' @param lambda1
#' @param ratio
#' @param nlambda2
#' @param fixe_lambda1
#'
#' @return fused lasso model
#' @export
fused_lasso <- function(X, response, G, lambda2 = NULL, lambda1  = 0,
                        ratio = 1e-3, nlambda2 = 100){
if (is.null(lambda2)) {
  lambdas <- fusedlassoMaxLambdas(X, response,
                                  family = "gaussian", graph = G)
  lambda2.max <- lambdas$maxLambda2

  lambda2     <- 10 ^ seq(log10(lambda2.max),
                          log10(ratio * lambda2.max), len = nlambda2)
  lambda1 <- (lambdas$maxLambda1 / lambdas$maxLambda2) * lambda2
}
weight_l1 <- rep(1, ncol(X))
mod<- FusedLasso::fusedlasso(X, response, graph = G,
         addIntercept = FALSE, family = "gaussian", lambda2 = lambda2,
         lambda1 = lambda1, wLambda1 = weight_l1, accuracy = 1e-6)
return(mod)
}



#' Title
#'
#' @param graphe
#' @param lambda2
#' @param lambda1
#' @param ratio
#' @param nlambda2
#' @param fixe_lambda1
#'
#' @return fused lasso model
#' @export
cv_fused_lasso <- function(X, response, G, lambda2 = NULL, lambda1  = 0,
                        ratio = 1e-3, nlambda2 = 100){
  browser()
  if (is.null(lambda2)) {
    lambdas <- fusedlassoMaxLambdas(X, response,
                                    family = "gaussian", graph = G)
    lambda2.max <- lambdas$maxLambda2

    lambda2     <- 10 ^ seq(log10(lambda2.max),
                            log10(ratio * lambda2.max), len = nlambda2)
    lambda1 <- (lambdas$maxLambda1 / lambdas$maxLambda2) * lambda2
  }
  dat <- X %>%  as.data.frame() %>%
    mutate(response = response) %>% vfold_cv()
  weight_l1 <- rep(1, ncol(X))
  my_fus <- partial(fusedlasso,  graph = G,
              addIntercept = FALSE, family = "gaussian", lambda2 = lambda2,
              lambda1 = lambda1, wLambda1 = weight_l1, accuracy = 1e-6)
 de<- dat %>%
    mutate(reg = map(splits, function(x){
      my_fus(y = as.tibble(x) %>% pull(response),
             X = as.tibble(x) %>% select(-response) %>%
               as.matrix() )
     }),
     mse = map2(splits, reg)
    )
  return(mod)
}

mse_fus <- function(x, mod){
  y <- as.tibble(x) %>% pull(response)
  X <- as.tibble(x) %>% select(-response) %>%
    as.matrix()

  mse <- (y- X %*%  mod$beta)^2 %>%
    colSums() / length(y)
}


