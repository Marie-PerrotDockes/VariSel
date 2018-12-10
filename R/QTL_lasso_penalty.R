#' Description of the function
#' @export
QTLmod_lasso_adaptive <- R6Class("QTLmod_lasso",
                        inherit = QTLmod,

                        public = list(
                          initialize = function(x, y){
                            super$initialize(x, y)
                          },
                          estime = function(lambda = NULL, standardize = FALSE, intercept = TRUE){

                            self$mod <-  private$y %>%
                              gather() %>%
                              group_by(key) %>%
                              nest(.key = Data) %>%
                              mutate(Data = map(Data, ~.$value)) %>%
                              mutate(Model = purrr::map(Data,
                                                        ~ cv.glmnet(private$x, .,
                                                                 lambda = lambda,
                                                                 standardize = standardize,
                                                                 intercept = intercept))) %>%
                              mutate(Penalty_factor = map(Model ~ (coef(.)[-1] + 1 / sqrt(nrow(private$x))^(-1)))) %>%
                              mutate(Model = map2(Model, Penalty_factor,
                                                  ~ glmnet(private$x, .x,
                                                           penalty.factor = .y,
                                                               lambda = lambda,
                                                                standardize = standardize,
                                                                intercept = intercept)))
                            self$res <- self$mod  %>%
                              mutate( Beta = map(Model, ~.$beta), Intercept = map(Model,~.$a0)) %>%
                              select(-Model)


                            private$df <-self$res %>%
                              select(Beta) %>%
                              map( function(x){
                                x %>% as.matrix() %>%
                                  as.data.frame() %>%
                                  summarise_all(~sum(.!=0))
                              })
                          },

                          sel_cv = function(s = "lambda.min"){
                            if(is.null(self$cv) || (s != private$s)){

                              private$s <- s
                              self$cv <- private$y %>%
                                gather() %>%
                                group_by(key) %>%
                                nest(.key = Data) %>%
                                mutate(Data = map(Data, ~.$value)) %>%
                                mutate(Model = map(Data,~cv.glmnet(private$x, .)))%>%
                                mutate(Penalty_factor = map(Model, ~ (coef(.)[-1] + 1 / sqrt(nrow(private$x))^(-1)))) %>%
                                mutate(Model = map2(Data, Penalty_factor,~cv.glmnet(private$x, .x,
                                                                 penalty.factor = .y
                                                             )))

                            }

                            self$cv <- self$cv %>%
                              mutate(Beta = map(Model, ~coef(.,"lambda.min")[-1])) %>%
                              select(Beta, Data) %>% mutate(Beta = set_names(Beta, colnames(private$y)),
                                                            Data = set_names(Data, colnames(private$y))) %>%

                              summarise_all(~list(bind_rows(.)))

                          }
                        ))
