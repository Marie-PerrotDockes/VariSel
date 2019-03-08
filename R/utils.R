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
#' @export
trans_beta <- function(beta, name){
  beta %>% as.data.frame() %>%
    rownames_to_column() %>%
    separate(rowname, c("Trait", "Marker"), "_", extra = "merge") %>%
    group_by(Trait) %>%
    nest(.key = Beta) %>%
    mutate(Beta = map(Beta, function(x){
      x[match(name, x$Marker), ] %>%
        select (-Marker)
    }
    )) %>%
    select(-Trait)
}



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
#' @export
ord_beta <- function(beta, name){
  beta %>% as.data.frame() %>%
    rownames_to_column() %>%
    separate(rowname, c("Trait", "Marker"), "_", extra = "merge") %>%
    mutate(ord = match(Marker, name)) %>% arrange(ord) %>%
    select(-Trait, -Marker, -ord)
}



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
#' @export
univ_y <- function(y){
  y %>%
   gather(key = "Trait") %>%
    group_by(Trait) %>%
    nest(.key = Data) %>%
    mutate(Data = map(Data, ~.$value))

}

#' Description of the function
#'
#' @param  b a vector with non null values
#' @param  b_hat an estimation of this vector
#' @return The proportion of non nulle values of B found non null in B_hat
#' @examples
#' B <- c(1, -1, 1.5, 1.5, rep(0, 6), 2, 0, 2, 0)
#' B_hat <- c(0, 0, 1, 1, rep(0,4) , rep(1 , 0))
#' TPR(b_hat = B_hat , b =B)
#' @export
TPR <- function(b_hat, b){
  TP <- sum(b_hat != 0 & b != 0)
  P <- sum(b != 0)
  TP / P
}


#' Description of the function
#'
#' @param  b a vector with non null values
#' @param  b_hat an estimation of this vector
#' @return The proportion of null values of B found non null in B_hat
#' @examples
#' B <- c(1, -1, 1.5, 1.5, rep(0, 6), 2, 0, 2, 0)
#' B_hat <- c(0, 0, 1, 1, rep(0,4) , rep(1 , 0))
#' FPR(b_hat = B_hat , b =B)
#' @export
FPR <- function(b_hat, b){
  FN <- sum(b_hat != 0 & b == 0)
  N  <- sum(b == 0)
  FN / N
}



#' Title
#'
#' @param group
#'
#' @return toc
#' @export
getGraphe <- function(group){
  G      <- list(conn = list(), weight = list())
  G[[1]] <- lapply(seq_along(group), function(i){
    as.integer(setdiff(which(group == group[i]), i))
  })
  G[[2]] <- lapply(seq_along(group), function(i){
    rep(1, (length(which(group == group[i])) - 1))
  })
  return(G)

}

#' @export
sel_ols <- function(b, y, x){
"touc"
  }
