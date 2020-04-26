type_to__S12inv <- function(X, Y, type){

  if(nrow(X) > ncol(X)){residus <- lm(as.matrix(Y)~X -1)$residuals
  } else{
    residus <- sapply(1:q,function(i){Y[,i]-X%*%MASS::lm.ridge(as.matrix(Y[,i])~X -1)$coef})
  }


}


#' Title
#'
#' @param X design matrix
#' @param Y response matrix
#' @param type type of models
#' @param Sigma_12inv the inverse of the square root of the covariance matrix of Y
#' @param group an optional vector with 2 modalities use only if type = "fus2mod_univ",
#'  the coefficients will be abble to be distinct for the two groups but encourage to be fused
#'
#' @import tidyr dplyr tibble
#' @importFrom R6 R6Class
#' @return
#' @export
type_to_varisel <- function(X, Y, type,
                            Sigma_12inv = diag(1, ncol(as.data.frame(Y))),
                            group,  a = 1,  sep = "\\."){
  mod <- switch(type,
    "group_univ"         =  mod_group_univ$new(X, Y, sepx = sep),

    "group_multi_both"   = mod_group_multi_both$new(X, Y,
                            Sigma_12inv = Sigma_12inv, sepx = sep),

    "group_multi_regr" = mod_group_multi_regr$new(X, Y,
                            Sigma_12inv = Sigma_12inv, sepx = sep),

    "fused_univ"         = mod_fused_univ$new(X, Y, sepx = sep),

    "fused_multi_both"   = mod_fused_multi_both$new(X, Y,
                            Sigma_12inv = Sigma_12inv, sepx = sep),

    "fused_multi_regr"   = mod_fused_multi_regr$new(X, Y,
                            Sigma_12inv = Sigma_12inv, sepx = sep),

    "lasso_univ"         = mod_lasso$new(X, Y, univ = TRUE),

    "lasso_multi"        = mod_lasso$new(X, Y,
                            Sigma_12inv = Sigma_12inv,
                            univ = FALSE),

    "fus2mod_univ"       = mod_fus2_univ$new(X, Y, group = group,  a = a),

    "fus2resp"           = mod_fus2_resp$new(X, Y,  a = a))
  return(mod)
}



#' Title
#'
#' @param rsamp a splits object with a collumn: response the response vector
#' @param type type of penalty to use.
#'
#' @return
#' @export
rsamples_to_mse <- function(rsamp, type, resp,
                           Sigma_12inv = diag(1, length(resp)),
                           lambda = NULL, grp = NULL,  sep = "\\."){
 Y_tr <- rsamp %>% analysis() %>% dplyr::select(resp)
 X_tr <- rsamp %>% analysis() %>% dplyr::select(-c(resp, grp))
 Y_t <- rsamp %>% assessment() %>% dplyr::select(resp)
 X_t <- rsamp %>% assessment() %>% dplyr::select(-c(resp, grp))
 if(!is.null(grp)){
   group_tr <- rsamp %>% analysis() %>% dplyr::select(grp)
   group_t <- rsamp %>% assessment() %>% dplyr::select(grp)
 } else{
   group_tr <- NULL
   new_group <- NULL
 }
 mod <- type_to_varisel(X_tr, Y_tr, type,
                        Sigma_12inv, group = group_tr,
                        sep = sep)
 mod$predict(new_x = X_t, lambda = lambda, new_group = new_group_t,
             names_y = resp)
 if (!mod$univ | type =="fus2resp"){
   y     <- list(as.numeric(as.matrix(Y_t)))
   trait <- rep(colnames(Y_t), each = nrow(Y_t))
   a     <- mod$res %>%
     mutate(New_y  = y,
            MSE = map2(New_pred, New_y, ~ (.x - .y) ^ 2 %>%
                        as.matrix() %>%
                        as.data.frame() %>%
                         mutate(Trait = trait) %>%
                        gather(key, value, -Trait, factor_key = FALSE) %>%
                        group_by(key, Trait) %>%
                        summarise(value = mean(value)))
            ) %>%
     ungroup() %>%
     dplyr::select(MSE) %>%
     unnest(MSE) %>%
     mutate(type = type)
    } else {
   y <- Y_t %>%
     gather(key = "Trait") %>%
     group_by(Trait) %>%
     nest() %>%
     mutate(Data = map(data, ~.$value)) %>%
     pull(Data)
   a <- mod$res %>%
    mutate(New_y  = y,
           MSE = map2(New_pred, New_y, ~ (.x - .y) ^ 2 %>%
                       as.data.frame() %>%
                       gather(factor_key = FALSE) %>%
                       group_by(key) %>%
                       summarise(value = mean(value)))
           ) %>%
    dplyr::select(MSE, Trait) %>%
    unnest(MSE) %>%
    mutate(type = type)
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
bt_error <- function(X, Y, type,
                     Sigma_12inv = diag(1, ncol(as.data.frame(Y))),
                     group = NULL, a = 1, times = 10,
                      sep = "\\."){
  if(!is.data.frame(Y))  Y <- as.data.frame(Y)
  if (is.null(colnames(Y))) colnames(Y) <- paste("rep", 1:ncol(Y), sep = "_")
  resp <- colnames(Y)
  mod_tot <- type_to_varisel(X, Y, type, Sigma_12inv, group, sep = sep)
  mod_tot$estime()

  if(type == "fus2mod_univ"){
    res_MSE <- cbind.data.frame(X, Y, group) %>%
      bootstraps(times) %>%
      mutate(MSE = future_map(splits,
        ~rsamples_to_mse(., type, resp,
            Sigma_12inv, grp = group,
            lambda = mod_tot$res$Lambda))) %>%
      unnest(MSE)
    } else{
      res_MSE <- cbind.data.frame(X, Y) %>%
        bootstraps(times) %>%
        mutate(MSE = future_map(splits,
          ~rsamples_to_mse(., type, resp,
            Sigma_12inv,
            lambda = mod_tot$res$Lambda))) %>%
        unnest(MSE)
    }
list(res_MSE, mod_tot)
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
compar_type <- function(X = NULL, Y, types,
                Sigma_12inv = NULL,
                group= NULL, a = 1, times = 10,  sep = "\\.",
                 regressors = NULL,
                type_S12_inv ="emp",
                p = NULL, q = NULL){
  if (is.null(X)){
    if(is.null(group) & is.null(regressors)) stop('Error if the design matrix X is not supplied at least one of group or regressors must be')

    if(is.null(regressors)) regressors <- matrix(1,nrow= length(group))
    if(!is.matrix(regressors)) regressors <- as.matrix(regressors)
    X <- model.matrix(~regressors:group -1)
    colnames(X) <- gsub('regressors', "", gsub("group", "", colnames(X)))
    sep <- ":"
  }
  if(!is.matrix(X)) X <- as.matrix(X)
    if(ncol(Y) !=1 & is.null(Sigma_12inv)) Sigma_12inv <- get_S12inv(Y, X, type_S12_inv, p = p, q =q)
  result <- types %>% as.list() %>% tibble() %>%
   transmute(compar = map(., ~bt_error(X, Y, .,
                              Sigma_12inv = Sigma_12inv, group = group, a = a ,
                              times = times,  sep = sep)))
  Models <- result %>% mutate(Models = map(compar, ~extract2(.,2))) %>%
    dplyr::select(-compar) %>% mutate(type = types)
  res_MSE <- result %>% mutate(MSE = map(compar, ~extract2(.,1))) %>%
    dplyr::select(-compar) %>%
   unnest()
  all_res <- Models %>%
    mutate(Trait = map(Models, ~.$trait),
           BIC = map(Models,~if(length(.$res$BIC) != 1){
             .$res$BIC
           }else{
             rep(.$res$BIC, length(.$trait))
           }),
           Lambda = map(Models,~if(length(.$res$Lambda) != 1){
             .$res$Lambda
           }else{
             rep(.$res$Lambda, length(.$trait))
           }),
           Name = map(Models,~if(length(.$name_y) != 1){
             .$name_y
           }else{
             rep(.$name_y, length(.$trait))
           })
           ) %>%
    dplyr::select(-Models) %>%
    unnest() %>%
    mutate(BIC = map(BIC, ~gather(., key = "key", value ="BIC")),
           Lambda = map(Lambda, ~ if(is.list(.)){
             as_tibble(.)
           } else{
             enframe(., name = NULL, value ="lambda1")
           })
    ) %>%
    unnest(Trait) %>%
    unnest() %>%
    left_join(res_MSE, by = c("key", "type", "Trait")) %>%
    rename(MSE_boot = value)


  return(list(res_MSE= res_MSE, all_res= all_res, Models = Models))
}


#' Title
#'
#' @param ct
#'
#' @return
#' @export
#'
#' @examples
plot_ct <- function(ct, criterion ="MSE"){
  if(criterion == "MSE"){
    p <- ct$all_res %>% group_by(key, type) %>% filter(type!="fus2resp") %>%
    # summarise(MSE = mean(MSE_boot), BIC = mean(BIC), lambda1 = mean(lambda1)) %>%
    ggplot(aes(x = lambda1, color = type, fill = type, y = MSE_boot, shape = type)) +
    geom_smooth() + theme_bw() + scale_x_log10()+
    labs(y = "MSE", title ="Bootstrap MSE", x = "Regularization Path")
  }
  if(criterion == "BIC"){
    p <- ct$all_res %>% group_by(key, type) %>% filter(type!="fus2resp") %>%
      # summarise(MSE = mean(MSE_boot), BIC = mean(BIC), lambda1 = mean(lambda1)) %>%
      ggplot(aes(x = lambda1, color = type, fill = type, y = BIC, shape = type)) +
      geom_smooth() + theme_bw() +  scale_x_log10()+
      labs(y = "BIC", title ="BIC", x = "Lambda")
  }
p
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
chapeau <- function(X, Y, type,
            Sigma_12inv = diag(1, ncol(as.data.frame(Y))),
            group= NULL, a = 1){

  mod  <- type_to_varisel(X, Y, type, Sigma_12inv, penalty.factor)
  mod$estime()
  univ <- mod$univ
  res  <- mod$res %>%
    mutate(Beta = map(Beta, ~as.data.frame(t(as.matrix(.))) %>%
                         rowid_to_column(var = "num_lambda")),
           Lambda = map(Lambda, ~as_tibble(.) %>%
                          dplyr::rename(Lambda = value))) %>%
    dplyr::select(Lambda, Beta, Trait) %>%
    unnest(Lambda, Beta) %>%
    gather(-Trait, -Lambda, -num_lambda, key = Marker, value = value) %>%
    filter(value != 0) %>% arrange(num_lambda)
 if (type == "fus2mod_univ"){
   res <- res %>% separate(Marker, sep = ":X", into = c("Group", "Marker")) %>%
     mutate(Group = gsub("group", "", Group))
 }
 if (!univ) {
  res <-  res %>% dplyr::select(-Trait) %>%
    separate(Marker, sep = mod$sepy, into = c("Trait", "Marker"))
 }

 return(res)
}


cv.chapeau <- function(X, Y, univ, type,
               Sigma_12inv = diag(1, ncol(as.data.frame(Y)))){
"d"
}

predict.chapeau <- function(X, Y, univ, type,
                    Sigma_12inv = diag(1, ncol(as.data.frame(Y)))){

"d"
}



#' Title
#'
#' @param ct
#' @param criterion
#' @param sepy
#'
#' @return
#' @export
#'
#' @examples
get_best_models <- function(ct, criterion = "MSE_boot", sepy="_"){
  if(criterion == "MSE_boot"){
Result <- ct$all_res %>% group_by(key, Name, type) %>%
    dplyr::select(-id) %>% summarise_if(is.numeric, mean) %>% ungroup() %>%
    group_by(Name, type) %>%
    slice(which.min(!!parse_expr(criterion))) %>%
    mutate(num = case_when(
      grepl("s", key)~ as.numeric(gsub("s","",key)) + 2,
      grepl("V", key)~ as.numeric(gsub("V","",key)) + 1
    )) %>% group_by(type) %>% nest() %>%
    left_join(ct$Models, by ="type") %>%
    mutate(data = map2(data,Models, ~.x %>%
         mutate(Beta=.y$res$Beta ,
          Beta = map2(num, Beta, ~.y %>%
            as.matrix() %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            dplyr::select(coef =.x,rowname) %>%
            separate(rowname, sep = "__",
                     into = c("Trait", "Marker"),
                     fill = "left"))) %>%
         unnest(Beta,.drop = FALSE))
                       ) %>%
    unnest(data) %>%
    mutate(Trait = case_when(
      is.na(Trait) ~ Name,
     ( !is.na(Trait)) ~Trait
      ))
  }


  if(criterion == "MSE_boot_1se"){
    Result <- ct$all_res %>% group_by(key, Name, type) %>%
      dplyr::select(-id) %>% summarise(MSE_boot_mean = mean(MSE_boot), MSE_boot_se = sd(MSE_boot),
                                MSE_boot_sum = MSE_boot_mean + MSE_boot_se) %>% ungroup() %>%
      dplyr::group_by(Name, type) %>%
      mutate( MSE_min = min(MSE_boot_sum)) %>%
      filter(MSE_boot_mean < min(MSE_boot_sum)) %>%
      slice(which.max(MSE_boot_mean)) %>%
      mutate(num = case_when(
        grepl("s", key)~ as.numeric(gsub("s","",key)) + 2,
        grepl("V", key)~ as.numeric(gsub("V","",key)) + 1
      )) %>% group_by(type) %>% nest() %>%
      left_join(ct$Models, by ="type") %>%
      mutate(data = map2(data,Models, ~.x %>%
                           mutate(Beta=.y$res$Beta ,
                                  Beta = map2(num, Beta, ~.y %>%
                                                as.matrix() %>%
                                                as.data.frame() %>%
                                                rownames_to_column() %>%
                                                dplyr::select(coef =.x,rowname) %>%
                                                separate(rowname, sep = "__",
                                                         into = c("Trait", "Marker"),
                                                         fill = "left"))) %>%
                           unnest(Beta,.drop = FALSE))
      ) %>%
      unnest(data) %>%
      mutate(Trait = case_when(
        is.na(Trait) ~ Name,
        ( !is.na(Trait)) ~Trait
      ))
  }
  return(Result)
    }


#' Title
#'
#' @param bmd
#'
#' @return
#' @export
#'
#' @examples
plot_md <- function(bmd, types =NULL){
  if(is.null(types)) types <- unique(bmd$type)

  if (length(types) == 1){
    dat <- bmd %>% filter(coef!=0 & type !="fus2resp" & type %in% types)
    if(length(unique(dat$Marker)) > length(unique(dat$Trait))){
        ggplot(dat, aes(y = Trait, x = Marker, fill =coef)) +
        geom_tile() + theme(axis.text.x = element_text(angle =90))+
       scale_fill_viridis()
    }else{
      ggplot(dat, aes(x = Trait, y = Marker, fill =coef)) +
        geom_tile() + theme(axis.text.x = element_text(angle =90))+
        scale_fill_viridis()
      }
  }else{
  bmd %>% filter(coef!=0 & type !="fus2resp" & type %in% types) %>%
    ggplot(aes(y = type, x = Marker, fill =coef)) +
    geom_tile() + theme(axis.text.x = element_text(angle =90))+
    facet_wrap(~Trait, scale = "free") +
    scale_fill_viridis()
  }
}


