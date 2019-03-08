require(testthat)
require(tidyverse)
require(gglasso)
theme_set(theme_minimal())
options("tibble.max_extra_cols" = 5)
context("test-grp_lasso")


mt <- mtcars
names(mt)[c(3,5,7,9)] <- paste("a",names(mt)[c(3,5,7,9)],sep=".")
names(mt)[-c(1:3,5,7,9)] <- paste("B",names(mt)[-c(1:3,5,7,9)],sep=".")
group_real <- rep(2,9)
group_real[(c(3,5,7,9)-2)] <- 1
X <- mt[,-c(1:2)]
Y <- mt[,1:2]
test_that("group from get_grp", {

  group <- get_group(names(X), sep = "\\.")
   expect_equal(group, group_real)
})



test_that("univ coef eq gglasso", {
  mod <- QTLmod_group_univ$new(X, Y)
  mod$estime()
  mod$res$Beta
  ord <- order(group_real)
  group_ord <- group_real[ord]
  X_ord <- X[, ord]
  # Eg <- expand.grid(colnames(X), colnames(Y))
  # colnames(X_ord) <- paste(Eg[,2],Eg[,1], sep ="_")[ord]
  # y <- as.numeric(as.matrix(Y))
 tac <-    Y %>%
    gather(key = "Trait") %>%
    group_by(Trait) %>%
    nest(.key = Data) %>%
    mutate(Data = map(Data, ~.$value),
           Beta = map(Data, ~ gglasso(x = as.matrix(X_ord),
                                      y = ., group = group_ord)$beta %>%
             as.data.frame() %>%
             rownames_to_column())) %>%
   unnest(Beta, .drop = TRUE) %>%
   gather(-Trait, -rowname, key= key, value = value) %>%
   transmute(value = paste(round(value,5),Trait,rowname, key,sep="_")) %>%
   pull(value)

  touc <- mod$res %>%
    mutate(Beta = map(Beta,~ as.data.frame(.) %>%
                        rownames_to_column())) %>%
    unnest(Beta,.drop=TRUE) %>% gather(-Trait, -rowname, key= key, value = value)  %>%
    transmute(value = paste(round(value,5),Trait,rowname, key,sep="_")) %>%
    pull(value)

  expect_setequal(tac, touc )
})



test_that("multi_both coef eq gglasso", {
Sigma_12inv <- diag(1, ncol(as.data.frame(Y)))
  mod <- QTLmod_group_multi_both$new(X, Y,Sigma_12inv = Sigma_12inv)
  mod$estime()
  mod$res$Beta

  Eg <- expand.grid(colnames(X), colnames(Y))

  y <- as.numeric(as.matrix(Y) %*% Sigma_12inv)
  X2 <- kronecker(Matrix::t(Sigma_12inv), as.matrix(X))
  colnames(X2) <- paste(Eg[,2],Eg[,1], sep = "__")
  group_real2 <- rep(group_real, ncol(Y))
  ord <- order(group_real2)
  group_ord <- group_real2[ord]
  X_ord <- X2[, ord]
  # Eg <- expand.grid(colnames(X), colnames(Y))
  # colnames(X_ord) <- paste(Eg[,2],Eg[,1], sep ="_")[ord]
  # y <- as.numeric(as.matrix(Y))
  tac <-    gglasso(x = X_ord, y = y, group = group_ord)$beta %>%
                        as.data.frame() %>%
                        rownames_to_column() %>%
    separate(rowname, into = c("Trait", "Marker"), sep="__") %>% gather(-Trait, -Marker, key= key, value = value) %>%
    transmute(value = paste(round(value,5),Trait,Marker, key,sep="_")) %>% pull(value)


  touc <- mod$res %>% mutate(Beta = map(Beta,~ as.data.frame(.) %>%
                                          rownames_to_column()) ) %>%
    unnest(Beta,.drop=TRUE) %>%  select(-Trait) %>%
    separate(rowname, into = c("Trait", "Marker"), sep="__") %>%
    gather(-Trait, -Marker, key= key, value = value) %>%
    transmute(value = paste(round(value,5),Trait,Marker, key,sep="_")) %>% pull(value)
  expect_setequal(tac, touc )
})




test_that("multi_marker coef eq gglasso", {
  Sigma_12inv <- diag(1, ncol(as.data.frame(Y)))
  mod <- QTLmod_group_multi_marker$new(X, Y,Sigma_12inv = Sigma_12inv)
  mod$estime()
  mod$res$Beta

  group_real2 <- 0:(ncol(Y)-1) %>%
    map(~ group_real + max(group_real) * .) %>%
    unlist()
  Eg <- expand.grid(colnames(X), colnames(Y))

  y <- as.numeric(as.matrix(Y) %*% Sigma_12inv)
  X2 <- kronecker(Matrix::t(Sigma_12inv), as.matrix(X))
  colnames(X2) <- paste(Eg[,2],Eg[,1], sep = "__")
  ord <- order(group_real2)
  group_ord <- group_real2[ord]
  X_ord <- X2[, ord]
  # Eg <- expand.grid(colnames(X), colnames(Y))
  # colnames(X_ord) <- paste(Eg[,2],Eg[,1], sep ="_")[ord]
  # y <- as.numeric(as.matrix(Y))
  tac <-    gglasso(x = X_ord, y = y, group = group_ord)$beta %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    separate(rowname, into = c("Trait", "Marker"), sep="__") %>% gather(-Trait, -Marker, key= key, value = value) %>%
    transmute(value = paste(round(value,5),Trait,Marker, key,sep="_")) %>% pull(value)


  touc <- mod$res %>% mutate(Beta = map(Beta,~ as.data.frame(.) %>%
                                          rownames_to_column()) ) %>%
    unnest(Beta,.drop=TRUE) %>%  select(-Trait) %>%
    separate(rowname, into = c("Trait", "Marker"), sep="__") %>%
    gather(-Trait, -Marker, key= key, value = value) %>%
    transmute(value = paste(round(value,5),Trait,Marker, key,sep="_")) %>% pull(value)
  expect_setequal(tac, touc )
})

