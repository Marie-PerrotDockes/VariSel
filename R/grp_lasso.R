#' Description of the function
#' @export
get_group <-function(name_x, sep = "\\."){
  group <- name_x %>%
    str_split(sep) %>%
    map_chr(~paste(.[-length(.)], collapse = ".")) %>%
    as.factor() %>%
    as.integer()
}

#' This function give the group for r Traits if the group are just among the Marker and not among the Trait
#' @export
get_group_marker <- function(name_x, sep = "\\.", r){
  group <- get_group(name_x, sep)
  0:(r-1) %>%
    map(~ group + max(group) * .) %>%
    unlist()
}
#' Description of the function
#'
#' @param  X a matrix with colnames indicating group vaiables
#' @param  sep the separators between the group (for exemlpe the marker) and the specificity of the columns (for exemple the allele)
#' @param  r the number of Trait
#' @return The group for r Trait if the columns have to be group if they have the same Marker, it is the same group for all Trait
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
get_group_both <-function(name_x, sep = "\\.", r){
  rep(get_group(name_x, sep), r)
}
#' Description of the function
#' @export
grp_lasso <- function(X,y, group, 		dfmax = as.integer(max(group)) + 1,
                      pmax = min(dfmax * 1.2, as.integer(max(group)))){
  # Y <- as.matrix(Y)
  # r <- ncol(Y)
  # if(is.null(colnames(Y))) colnames(Y) <- 1:r
  ord <- order(group)
  group_ord <- group[ord]
  X_ord <- X[, ord]
  # Eg <- expand.grid(colnames(X), colnames(Y))
  # colnames(X_ord) <- paste(Eg[,2],Eg[,1], sep ="_")[ord]
  # y <- as.numeric(as.matrix(Y))
  mod <- gglasso(x = X_ord, y = y, group = group, dfmax = dfmax, pmax = pmax)
  a <- mod$beta
  mod$beta <- a[match(colnames(X),rownames(a)),]
  mod$group <- group
  mod
  }

#' Description of the function
#' @export
cv_grp_lasso <- function(X,y, group,s ="lambda.min"){
  # Y <- as.matrix(Y)
  # r <- ncol(Y)
  # if(is.null(colnames(Y))) colnames(Y) <- 1:r
  ord       <- order(group)
  group_ord <- group[ord]
  X_ord     <- X[,ord]
  mod       <- cv.gglasso(x = X_ord, y = y, group = group)
  # coef_min  <- coef(mod,"lambda.min") %>% as.data.frame() %>%
  #   rownames_to_column() %>%
  #   mutate(order = c(1,  match(colnames(X),rowname)+1)) %>%
  #   arrange(order) %>% select(-order)
  # coef_1se  <- coef(mod,"lambda.1se") %>% as.data.frame() %>%
  #   rownames_to_column() %>%
  #   mutate(order = c(1,  match(colnames(X),rowname)+1)) %>%
  #   arrange(order) %>% select(-order)
  # list(coef_min = coef_min, coef_1se = coef_1se, cvm = mod$cvm, cvsd = mod$cvsd)
  mod
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
grp_lasso_st <- function(X,y, group, nb.cores = 7, B = 500, PFER = 1){
  nc <- round(ncol(X) / max(group))
  q <- max(5, min(2 * round(nrow(X)/log(nc)/10)*10, nc))
  st <-  stabsel(X, y,  q=q, PFER = 1, B = B,
                 fitfun=gglasso_forstab,
                 args.fitfun=list( group = group),
                 mc.cores=nb.cores, mc.preschedule=TRUE)
  st
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
gglasso_forstab <- function (x, y, group, q)
{
  fit <- grp_lasso(X = x, y = y, group = group, pmax = q)
  selected <- which(fit$beta[,ncol(fit$beta)] !=0 )
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  names(ret) <- colnames(x)
  cf <- fit$beta[, ]
  sequence <- as.matrix(cf != 0)
  return(list(selected = ret, path = sequence))
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
gglasso_st_tot <- function(X,Y, group = NULL, sep ="\\.", nb.cores = 7, B = 500, PFER = 1, type_group){
  if (type_group == "Trait") {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    q <- max(5, min(2 * round(nrow(X)/log(ncol(X))/10)*10, ncol(X)))
    return(stabsel(X, Y, q = q, PFER = PFER, B = 500,
                   fitfun=glmnet.multivar,
                   args.fitfun=list( family="mgaussian",
                                    standardize=FALSE),
                   mc.cores=nb.cores, mc.preschedule=TRUE))
  } else{
    r <- ncol(Y)
    if (is.null(group))
      group <- switch (type_group,
                       both = get_group_both(X, sep, r),
                       marker = get_group_both(X, sep, r)
      )

    return(grp_lasso_st(X,Y, group, nb.cores, B, PFER))
  }

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
summary.stabsel <- function(st, cutoff = 0.85, map, mrk2lg){
    select <- stabsel(st, cutoff = cutoff)$selected

    # Trait  <- select %>% names() %>%
    #   str_split('_',2, simplify = FALSE) %>%
    #   map_dfc(~c(.[1],
    #          str_split("\\.", .[2]) %>%
    #            map_dfr(~c(paste(.[-length(.)], collapse = "."), .[length(.)]))))

    if (length(select) != 0){

    Tr  <- select %>% names() %>%
         str_split('_',2, simplify = FALSE)
    trait <- Tr %>%  map_chr(~ if (length(.) > 1 ){ .[1]
      } else {"all"})

    Mr <- Tr %>%
      map_chr(~.[length(.)]) %>%
      str_split( .,'\\.')

    marker <- Mr %>%
      map_chr(~.[1])

    allele <- Mr %>%
      map_chr(~.[2])

    mp <- suppressWarnings( as.numeric(map[,2]))
    names(mp) <- map$V1
    cbind.data.frame(trait, marker, allele, linkage_group = mrk2lg[marker], position =mp[marker] , Probability = st$max[select], idx = 0)
    } else {
      cbind.data.frame(trait = 0, marker = 0, allele = 0, linkage_group = 0, position =0 , probability = 0, idx = 0)

}
}

# summary.stabsel(st, map = map, mrk2lg)

#' Give the result for a group lasso
#'
#' @param  X a matrix with the value of different allele for different Marker
#' @param  Y a vector or a matrix with the values of different Trait
#' @param  type_group a character indicating the way of creating the group :
#' "both" means that the variable are selected for all a Marker for all the Trait.
#' "marker" means that the variable are selected for a Marker but for one Trait.
#' "trait" means that a variable are selected for all the Trait but not for all the marker but just for one allele
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
grpLassoQTL <- function(X, Y, marker, map,
                        PFER=1, B=500, cutoff=0.85,  nb.cores, sep ="\\.", mrk2lg, type_group="both"){
  Y   <- as.matrix(Y)
  st  <- gglasso_st_tot(X, Y, group = NULL, sep = sep, nb.cores = nb.cores, B = 500, PFER = PFER, type_group = type_group)
  summary(st, cutoff = cutoff, map, mrk2lg)
}

# md <-grp_lasso(X, Y)


