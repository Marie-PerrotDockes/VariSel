
#' Title
#'
#' @param X
#' @param Y
#' @param type
#' @param Sigma_12inv
#' @param group
#' @param a
#'
#' @return
#' @export
type_to_varisel <- function(X,Y, type, Sigma_12inv = diag(1, ncol(as.data.frame(Y))), penalty.factor){
  mod <- switch(type,
                "group_univ" =  mod_group_univ$new(X,Y),
                "group_multi_both" = mod_group_multi_both$new(X,Y, Sigma_12inv = Sigma_12inv),
                "group_multi_marker" = mod_group_multi_marker$new(X,Y, Sigma_12inv = Sigma_12inv),
                "fused_univ" =  mod_group_univ$new(X,Y),
                "fused_multi_both" = mod_group_multi_both$new(X,Y, Sigma_12inv = Sigma_12inv),
                "fused_multi_regr" = mod_group_multi_marker$new(X,Y, Sigma_12inv = Sigma_12inv),
                "lasso_univ" = mod_lasso$new(X,Y, univ = TRUE),
                "lasso_multi" = mod_lasso$new(X,Y, Sigma_12inv = Sigma_12inv, univ = FALSE),
                "fus2mod_univ" = mod_lasso$new(X, Y, univ = TRUE,
                                               penalty.factor = penalty.factor),
                "fus2resp" = mod_lasso$new(X, Y, univ = TRUE,
                                               penalty.factor = penalty.factor))
  return(mod)
}



#' Title
#'
#' @param rsamp a splits object with a collumn: response the response vector
#' @param type type for varisel object
#'
#' @return
#' @export
rsamples_to_mse<- function(rsamp, type, resp, Sigma_12inv = diag(1, length(resp)), penalty.factor = NULL){
 Y_tr <- rsamp %>% analysis() %>% select(resp)
 X_tr <- rsamp %>% analysis() %>% select(-resp)
 Y_t <- rsamp %>% assessment() %>% select(resp)
 X_t <- rsamp %>% assessment() %>% select(-resp)
 mod <- type_to_varisel(X_tr, Y_tr, type, Sigma_12inv, penalty.factor)
 mod$predict(new_x = X_t)
 if(!mod$univ){
   y <- (list(as.numeric(as.matrix(Y_t))))
   trait <- rep(colnames(Y_t), each = nrow(Y_t))
   a<- mod$res %>%  mutate(New_y  = y,
                           MSE = map2(New_pred, New_y, ~ (.x- .y)^2 %>%
                                        as.matrix() %>%
                                        as.data.frame() %>%
                                        mutate(Trait = trait ) %>%
                                        gather(key, value, -Trait) %>% group_by(key, Trait) %>%
                                        summarise(value = mean(value))))
   a %>% select(MSE) %>% unnest(MSE) %>%  mutate(type = type)

    } else {
   y <- Y_t %>%
     gather(key = "Trait") %>%
     group_by(Trait) %>%
     nest(.key = Data) %>%
     mutate(Data = map(Data, ~.$value)) %>%
     pull(Data)
 a<- mod$res %>%  mutate(New_y  = y,
                      MSE = map2(New_pred, New_y, ~ (.x- .y)^2 %>%
                                   as.data.frame() %>%
                                   gather() %>% group_by(key) %>%
                                   summarise(value = mean(value))))
 a %>% select(MSE,Trait) %>% unnest(MSE) %>%  mutate(type = type)
}

}


#' Title
#'
#' @param X
#' @param Y
#' @param type
#' @param Sigma_12inv
#' @param group
#' @param a
#'
#' @return
#' @export
bt_error <- function(X,Y, type, Sigma_12inv = diag(1, ncol(as.data.frame(Y))), group= NULL, a = 1, times = 10){
  if(type == "fus2mod_univ"){
    X <- as.matrix(X)
    p <- ncol(X)
    X1  <- model.matrix(~group + group:X - 1)
    X <- cbind(X1, X)
    b <- (3 * p - a * p + 2) / (2 * p )
    penalty.factor <- c(0, 0, rep(b, (ncol(X1) - 2)), rep(a, p))
  } else{
    if(type == "fus2resp"){
      X <- as.matrix(X)
      q <- ncol(Y)
      if(q !=2) stop("Y must have two columns to use fus2resp")
      if(is.null(colnames(Y))) colnames(Y) <- paste0("rep",1:q)
      X <- bdiag(rep(list(X), q)) %>% as.matrix()
      group <- rep(colnames(Y), each = nrow(Y))
      Y <- Y %>% as.matrix() %>%  as.numeric()
      p <- ncol(X)
      X1  <- model.matrix(~group + group:X - 1)
      X <- cbind(X1, X)
      b <- (3 * p - a * p + 2) / (2 * p )
      penalty.factor <- c(0, 0, rep(b, (ncol(X1) - 2)), rep(a, p))
    }else{
      penalty.factor <- NULL
    }
  }
  Y <- as.data.frame(Y)
  if(is.null(colnames(Y))) colnmaes(Y) <- paste('rep', 1:ncol(Y), sep='_')
  resp <- colnames(Y)
  cbind.data.frame(X, Y) %>%
    bootstraps(times) %>%
    mutate(MSE = map(splits, ~rsamples_to_mse(., type, resp, Sigma_12inv, penalty.factor))) %>%
    unnest(MSE)

}


#' Title
#'
#' @param X
#' @param Y
#' @param types
#' @param Sigma_12inv
#' @param group
#' @param a
#' @param times
#'
#' @return
#' @export
compar_type <- function(X, Y, types, Sigma_12inv = diag(1, ncol(as.data.frame(Y))),
                        group= NULL, a = 1, times = 10){
types %>% as.list() %>% tibble() %>%
  transmute(compar = map(.,~bt_error(X, Y, ., Sigma_12inv, group, a, times))) %>%
  unnest()
}

#' Title
#'
#' @param X
#' @param Y
#' @param type
#' @param Sigma_12inv
#'
#' @return
#' @export
chapeau <- function(X,Y, type, Sigma_12inv = diag(1, ncol(as.data.frame(Y))),
                    group= NULL, a = 1){
  if(type == "fus2mod_univ"){
    X <- as.matrix(X)
    p <- ncol(X)
    X1  <- model.matrix(~group + group:X - 1)
    X <- cbind(X1, X)
    b <- (3 * p - a * p + 2) / (2 * p )
    penalty.factor <- c(0, 0, rep(b, (ncol(X1) - 2)), rep(a, p))
  } else{
    if(type == "fus2resp"){
      X <- as.matrix(X)
      q <- ncol(Y)
      if(q !=2) stop("Y must have two columns to use fus2resp")
      if(is.null(colnames(Y))) colnames(Y) <- paste0("rep",1:q)
      X <- bdiag(rep(list(X), q)) %>% as.matrix()
      group <- rep(colnames(Y), each = nrow(Y))
      Y <- Y %>% as.matrix() %>%  as.numeric()
      p <- ncol(X)
      X1  <- model.matrix(~group + group:X - 1)
      X <- cbind(X1, X)
      b <- (3 * p - a * p + 2) / (2 * p )
      penalty.factor <- c(0, 0, rep(b, (ncol(X1) - 2)), rep(a, p))
    }else{
      penalty.factor <- NULL
    }
  }
  mod <- type_to_varisel(X,Y, type, Sigma_12inv, penalty.factor)

  univ <- mod$univ

  mod$estime()
 res<- mod$res %>%
    mutate(Beta = map(Beta, ~as.data.frame(t(.)) %>%
                        rowid_to_column(var="num_lambda")),
                     Lambda = map(Lambda,~as_tibble(.) %>% dplyr::rename(Lambda = value))) %>%
    select(Lambda, Beta, Trait) %>%
    unnest(Lambda, Beta) %>%
    gather(-Trait, -Lambda, -num_lambda, key=Marker, value = value) %>% filter(value!=0) %>% arrange(num_lambda)
 if(type == "fus2mod_univ"){
   res <- res %>% separate(Marker, sep =":X", into = c("Group", "Marker")) %>%
     mutate(Group = gsub('group', "", Group))
 }
 if(!univ) {
  res <-  res %>% select(-Trait) %>% separate(Marker, sep =mod$sepy, into = c("Trait", "Marker"))
 }

 return(res)
}


cv.chapeau <- function(X,Y, univ, type, Sigma_12inv = diag(1, ncol(as.data.frame(Y)))){
"d"
}

predict.chapeau <- function(X, Y, univ, type, Sigma_12inv = diag(1, ncol(as.data.frame(Y)))){

"d"
}
