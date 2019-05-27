train_VariSel <- function(X, Y, type, sepx, reg, group){
  if(!is.null(X) & !is.null(reg)) stop("either reg and group must be supplied or X (the design matirx already full)")
  if(!is.null(sepx) & !is.null(group)) stop("either sepx or group must be supplied see the help of the function")
  if (is.null(X)){
    if(is.null(group) | is.null(reg)) stop('Error if the design matrix X is not supplied both group and reg must be')
  X <- model.matrix(~reg:group -1)
  sepx <- ":"
  }
   mod <- type_to_varisel(X, Y, type,
                          Sigma_12inv = Sigma_12inv, sep= sepx)
   mod$estime()
}


## caret

# tuned_grid
