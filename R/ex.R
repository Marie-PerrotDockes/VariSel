# require(tidyverse)
# require(magrittr)
# require(gganimate)
# require(R6)
# miaou <- QTLmod_group_univ$new(X,Y)
# miaou$plot_coef()
# miaou$plot_anime()
# miaou$sel_cv()
#
#  miaou <- QTLmod_group_univ$new(X,Y)
#  miaou$plot_cv()
#  miaou <- QTLmod_group_multi_both$new(X,Y)
#  miaou$plot_cv()
#
# miaou <- QTLmod_group_multi_marker$new(X,Y)
#  miaou$plot_cv()
#  miaou <- QTLmod_lasso$new(X,Y)
#  miaou$plot_cv()
#  miaou <-mod_lasso$new(X,Y,univ=FALSE)
 # miaou$sel_cv()
# types <- c("group_univ" ,
# "group_multi_both",
# "group_multi_marker",
# "fused_univ",
# "fused_multi_both",
# "fused_multi_regr",
# "lasso_univ",
# "lasso_multi")
# test <- lapply(types, function(type){
#   print(type)
#   chapeau(X,Y,type)
# })
#
# group <- factor(c(rep("a", round(nrow(X)/2)),rep("b", nrow(X) - round(nrow(X)/2))))
# a<- chapeau(X,Y,"fus2mod_univ", group = group)
#
# require(rsample)
# require(tidyverse)
# iris2 <- iris %>% filter(Species %in% levels(Species)[1:2])
# rsamp <- iris2 %>%  select(-Species) %>%
#   bootstraps(5) %>% pull(splits) %>% first()
# resp <-c("Sepal.Length" ,"Sepal.Width")
#
# chapeau(iris2[,1:2], iris2[,3:4], type= "lasso_multi", resp, group = iris2$Species)
#
# X<-X[,700:800]
#
# ct <- compar_type(X, scale(Y), types = c(  "fus2resp", "group_univ" ,"group_multi_both" ,"group_multi_marker" ,
#                               "fused_univ" ,"fused_multi_both","fused_multi_regr",
#                               "lasso_univ","lasso_multi" ))
# p <- ct %>% mutate(key = as.numeric(gsub('s','',key))) %>%
#   ggplot(aes(x = key, color = type,fill = type, y = value))+geom_smooth() + theme_bw() +
# labs(y = "MSE", title ="Bootstrap MSE", x = "Regularization Path")
# ggsave(p , file = "ex_bootmse.pdf")
# group <- iris$Species
# X <- as.matrix(iris[,1:2])
# X1  <- model.matrix(~group + group:X - 1)
# X1 %>% as.tibble() %>% gather() %>%
#   mutate(Group = str_extract_all(key,paste(grp_value,collapse="|"))) %>% pull(Group)
#
# paste0('{',paste(paste0("group",unique(group)), collapse =","),'}')
