
#' Title
#'
#' @param X a design matrix of explicatives variables for the linear model
#' @param Y a matrix of responses
#' @param type the type of models to fit. See details for a description of the different types.
#' @param sepx if the type use fused lasso or group lasso on the columns of X it is the character separating the group value to the identifier values.
#'  For instance if the names of the columns of X are gene.allele  and if we want to group variable that belong to the same gene use sepx=".".
#' @param regressors if X is null and we are interested in the association between regressors and reponses in different group this is a matrix of regressors.
#' @param group this is a charcater vector containing the different groups. It side must be the number of row of the matrix regressors.
#' It indiciate for each rows in which group it is.
#' @param lambda default to null. A numeric vector containing values for the regularisation parameter for lasso and group-lasso penalties.
#' @param lambda1 default to null. A numeric vector containing values for the lambda1 parameters (sparsity parmater for the fused lasso penalties).
#' @param lambda2 default to null. A numeric vector containing values for the lambda2 parameters (fusion parmater for the fused lasso penalties).
#' @param Sigma_12inv a matrix of the square root of the inverse of the covariance matrix of the residuals. ( Use to remoove the dependance that may exist among the responses)
#' @param type_S12_inv if Sigma_12inv is null it can be estimated using different type. See details for the differents types available.
#'
#' @return a VariSel object.
#' @export
#'
#' @examples
train_VariSel <- function(X = NULL, Y, type, sepx = NULL, regressors =NULL, group=NULL, lambda = NULL,
                          lambda1=NULL, lambda2 = NULL, Sigma_12inv = NULL, type_S12_inv = "emp",p = NULL, q =NULL){

  if(!is.null(X) & !is.null(regressors)) stop("either regressors and group must be supplied or X (the design matirx already full)")
  if(!is.null(sepx) & !is.null(group)) stop("either sepx or group must be supplied see the help of the function")
  if (is.null(X)){
    if(is.null(group) & is.null(regressors)) stop('Error if the design matrix X is not supplied at least group or regressors must be')
    if(is.null(regressors)) regressors <- matrix(1,nrow= length(group))
    if(!is.matrix(regressors)) regressors <- as.matrix(regressors)
    X <- model.matrix(~regressors:group -1)
    colnames(X) <- gsub('regressors', "", gsub("group", "", colnames(X)))
    sepx <- ":"
  }
  if(!is.matrix(X)) X <- as.matrix(X)
  if(grepl("fused", type) & !is.null(lambda1)){
    lambda = list(lambda1 = lambda1, lambda2 = lambda2)
  }
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  if(ncol(Y) !=1 & is.null(Sigma_12inv)) Sigma_12inv <- get_S12inv(Y, X, type_S12_inv, p = p, q = q)
   mod <- type_to_varisel(X, Y, type,
                          Sigma_12inv = Sigma_12inv, sep= sepx)
   mod$estime(lambda = lambda)
   return(mod)
}


## caret

#' Title
#'
#' @param Y
#' @param X
#' @param type
#'
#' @return
#' @importFrom MASS lm.ridge
#' @export
#'
#' @examples
get_S12inv <- function(Y, X, type, p = NULL, q = NULL){

  if(nrow(X) > ncol(X)){
    res <- lm(as.matrix(Y)~X -1)$residuals
  } else{
    res <- sapply(1:ncol(Y),function(i){Y[,i]-X%*%lm.ridge(as.matrix(Y[,i])~X -1)$coef})
  }
  if(type =="emp"){
   S12_inv <-  solve(chol(cor(res)))
  }
  if (type =="Block"){
    S12_inv<- Sigma_estimation(res, inv_12 = TRUE)$S_inv_12
  }
  if(! type %in% c("emp", "Block")){
    S12_inv <-  whitening(res, type, p , q)
  }
  return(S12_inv)
}

# tuned_grid


#' Title
#'
#' @param mod
#'
#' @return
#' @export
#'
#' @examples
coef.VariSel<- function(mod){
mod$get_coef()

  return(mod$coef)
}


##
#' Title
#'
#' @param mod
#'
#' @return
#' @export
#'
#' @examples
plot.VariSel <- function(mod, type ="first", nb = 6){

  mod$plot_path(type= type, nb = nb)

    }


#' Title
#'
#' @param mods
#'
#' @return
#' @export
#'
#' @examples
compar_path <- function(mods ){

 p<-  tibble(mods) %>%
    mutate(type = map_chr(mods, ~class(.)[1]),
           coef = map(mods, ~coef(.))) %>%
    separate(type, sep ="_", into =c("Mod", "Type", "Univ","Group")) %>%
    select(-mods) %>%
    unnest() %>%
    mutate(col  =  paste(Reg,  Trait, sep =" on "),
           col = factor(col, levels = unique(col))) %>%
    ggplot(aes(x = Lambda, y = value,linetype = group,group = paste(col, group),color =col))+
    geom_line() +
    scale_x_log10() +
    labs ( color = " ", y = "value of the coefficients", title = "Regularization Path", x = "Lambda") +
    facet_grid(Type~Group)

 print(p)

}
