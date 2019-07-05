# VariSel_mod <- list(type = "both",
#                     library = "VariSel",
#                     loop = NULL)
#
# prm <- data.frame(parameter = c("type", "lambda", "lambda1", "lambda2"),
#                   class = c("character",rep("numeric", 2)),
#                   label = c("Type", "Lambda", "Lambda1", "Lambda2"))
#
# VariSel_mod$parameters <- prm
#
# VariSelGrid <- function(X, Y, Sigma_12inv,  search = "grid") {
#   y <-  as.numeric(Y %*% Sigma_12inv)
#   x <-    kronecker(Matrix::t(Sigma_12inv), X)
#
#
#   library(kernlab)
#   ## This produces low, middle and high values for sigma
#   ## (i.e. a vector with 3 elements).
#   sigmas <- kernlab::sigest(as.matrix(x), na.action = na.omit, scaled = TRUE)
#   ## To use grid search:
#   if(search == "grid") {
#     out <- expand.grid(sigma = mean(as.vector(sigmas[-2])),
#                        C = 2 ^((1:len) - 3))
#   } else {
#     ## For random search, define ranges for the parameters then
#     ## generate random values for them
#     rng <- extendrange(log(sigmas), f = .75)
#     out <- data.frame(sigma = exp(runif(len, min = rng[1], max = rng[2])),
#                       C = 2^runif(len, min = -5, max = 8))
#   }
#   out
# }
# filter <- ksvm(type~.,data=spamtrain,kernel="rbfdot",
#                kpar=list(sigma=c(0.05,0.1)),C=c(5,2),cross=3)
# filter
#
#
# get_lambda <- function(x, y, type){
#   if(grepl("lasso", type)){
#   mysd <- function(y) sqrt(sum((y-mean(y))^2)/length(y))
#   sx <- scale(x,scale=apply(x, 2, mysd))
#   lambda.max <- max(abs(t(sx) %*% y))/length(y)
#   lambda.min.ratio <- ifelse(nrow(sx)< ncol(sx), 0.01,0.0001)
#   nlambda <- 100
#   lambdapath <- round(exp(seq(log(lambda_max),
#                               log(lambda_max*lambda.min.ratio),
#                               length.out = nlambda)), digits = 10)
#   }
# }
