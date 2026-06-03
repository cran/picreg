## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse   = FALSE,
  comment    = "##",
  fig.width  = 5,
  fig.height = 3,
  fig.align  = "center"
)

# Optional truncation of long printed outputs, used selectively below.
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
  if (!is.null(n <- options$out.lines)) {
    x <- xfun::split_lines(x)
    if (length(x) > n) {
      x <- c(head(x, n), "...")
    }
    x <- paste(x, collapse = "\n")
  }
  hook_output(x, options)
})

set.seed(42)

## ----eval = FALSE-------------------------------------------------------------
# # install.packages("remotes")
# remotes::install_github("VcMaxouuu/picreg")

## -----------------------------------------------------------------------------
library(picreg)

## ----gen-data-----------------------------------------------------------------
data(QuickStartExample)
X <- QuickStartExample$X
y <- QuickStartExample$y

## ----fit-basic----------------------------------------------------------------
fit <- pic(X, y)

## ----print-fit----------------------------------------------------------------
print(fit)

## ----fit-fields---------------------------------------------------------------
# Number of selected features and their identifiers
fit$df
fit$selected

# Selected regularisation parameter
fit$lambda

# Family and penalty descriptors
fit$family
fit$penalty

## ----coef, out.lines = 12-----------------------------------------------------
coef(fit)                          # original scale
# coef(fit, standardized = TRUE)   # standardised scale; same as fit$beta

## ----predict, out.lines = 6---------------------------------------------------
predict(fit, newx = X[1:5,])                  # response
predict(fit, newx = X[1:5,], type = "link")   # linear predictor

## ----plot-fit, fig.height = 3-------------------------------------------------
plot(fit)

## ----plot-pdb, fig.height = 3.2-----------------------------------------------
plot(fit$lambda_pdb)

## ----pdb-summary--------------------------------------------------------------
pdb_summary(fit)

## ----manual-lambda------------------------------------------------------------
fit_lasso <- pic(X, y, penalty = "lasso", lambda = fit$lambda)
fit_scad  <- pic(X, y, penalty = "scad",  lambda = fit$lambda)
fit_mcp   <- pic(X, y, penalty = "mcp",   lambda = fit$lambda)

data.frame(
  penalty = c("lasso", "scad", "mcp"),
  df      = c(fit_lasso$df, fit_scad$df, fit_mcp$df),
  lambda  = c(fit_lasso$lambda, fit_scad$lambda, fit_mcp$lambda)
)

## ----binomial-ex--------------------------------------------------------------
data(BinomialExample)
X <- BinomialExample$X
y <- BinomialExample$y

fit_binom <- pic(X, y, family = 'binomial')
predict(fit_binom, newx = X[1:5, ], type = 'class')

## ----cox-example--------------------------------------------------------------
data(CoxExample)
X_cox <- CoxExample$X
y_cox <- CoxExample$y
print(y_cox[1:7, ])

fit_cox <- pic(X_cox, y_cox, family = "cox")
print(fit_cox)

## ----cox-baseline, fig.height = 4.5-------------------------------------------
plot_baseline(fit_cox)

## ----cox-survival, fig.height = 3.2-------------------------------------------
sf <- predict_survival_function(fit_cox, newx = X_cox[1:3, ])
plot_survival_curves(sf)

## ----cox-feature-effect, fig.height = 3.5-------------------------------------
fx <- feature_effects_on_survival(fit_cox, idx = "gene_1")
plot_survival_curves(fx)

## ----cox-feature-effect-custom, fig.height = 3.5------------------------------
fx <- feature_effects_on_survival(fit_cox,
                                  idx    = "gene_3",
                                  values = c(-2, -1, 0, 1, 2))
plot_survival_curves(fx)

## ----assess-gaussian----------------------------------------------------------
data(QuickStartExample)
X <- QuickStartExample$X
y <- QuickStartExample$y
fit <- pic(X, y)

assess(fit, newx = X, newy = y)

## ----assess-cox---------------------------------------------------------------
data(CoxExample)
fit_cox <- pic(CoxExample$X, CoxExample$y, family = "cox")
assess(fit_cox, newx = CoxExample$X, newy = CoxExample$y)

## ----assess-relax-------------------------------------------------------------
fit_relax <- pic(X, y, relax = TRUE)
assess(fit_relax, newx = X, newy = y)

## ----assess-with-truth--------------------------------------------------------
true_active <- paste0("gene_", 1:5)
assess(fit, newx = X, newy = y, true_features = true_active)

## ----phase-transition, eval=FALSE---------------------------------------------
# pt <- phase_transition(
#   n        = c(50, 100),
#   p        = c(100, 100),
#   type     = "gaussian",
#   s_max    = 20,
#   m        = 100,
#   penalty  = c("lasso", "scad"),
#   parallel = TRUE
# )
# plot(pt)

## ----phase-transition-img, echo = FALSE, out.width = "85%", fig.align = "center"----
knitr::include_graphics("phase_transition_pic.png")

## ----eval=FALSE---------------------------------------------------------------
# plot(pt, metric = "fdr")

## ----phase-transition-fdr-img, echo = FALSE, out.width = "85%", fig.align = "center"----
knitr::include_graphics("phase_transition_pic_fdr.png")

## ----pdb-asymptotic, fig.height = 3.5-----------------------------------------
as_ <- pdb_asymptotic(
  n_grid = c(20, 50, 100, 200),
  p      = 200,
  type   = "binomial",
  n_simu = 2000L
)

plot(as_)

## ----compare-glmnet-data, eval = requireNamespace("glmnet", quietly = TRUE), include = FALSE----
library(glmnet)

## ----compare-glmnet-sim-------------------------------------------------------
set.seed(1)
n <- 100; p <- 100; s <- 5
X <- matrix(rnorm(n * p), n, p)
beta <- numeric(p)
beta[1:s] <- 3
y <- as.numeric(X %*% beta + rnorm(n))
true_support <- seq_len(s)

## ----compare-glmnet-fit, eval = requireNamespace("glmnet", quietly = TRUE), echo = FALSE----
# Hidden: the fits themselves and the small helpers that pull the
# active set out of each fitted object.
fit_pic <- pic(X, y, family = "gaussian", penalty = "lasso")
sel_pic <- fit_pic$selected

fit_cv <- cv.glmnet(X, y, family = "gaussian", alpha = 1)
get_support <- function(fit, lambda) {
  sel <- which(as.numeric(stats::coef(fit, s = lambda)) != 0)
  sel <- sel[sel > 1L] - 1L
  sel
}
sel_min <- get_support(fit_cv, "lambda.min")
sel_1se <- get_support(fit_cv, "lambda.1se")

metrics <- function(sel, true_support) {
  data.frame(
    selected       = length(sel),
    exact_recovery = setequal(sel, true_support)
  )
}

## ----compare-glmnet-table, eval = requireNamespace("glmnet", quietly = TRUE)----
rbind(
  cbind(method = "picreg (lasso)",         metrics(sel_pic, true_support)),
  cbind(method = "cv.glmnet (lambda.min)", metrics(sel_min, true_support)),
  cbind(method = "cv.glmnet (lambda.1se)", metrics(sel_1se, true_support))
)

## ----compare-phase-transition, eval = FALSE, echo = FALSE---------------------
# set.seed(2)
# m      <- 100L                    # Monte Carlo replications per s
# s_max  <- 12L
# n      <- 100L; p <- 100L
# 
# recovery <- function(method, s) {
#   hits <- replicate(m, {
#     X <- matrix(rnorm(n * p), n, p)
#     b <- numeric(p)
#     true_idx <- if (s > 0L) sample.int(p, s) else integer(0)
#     if (s > 0L) b[true_idx] <- 3
# 
#     y <- as.numeric(X %*% b + rnorm(n))
# 
#     sel <- switch(
#       method,
# 
#       picreg = {
#         pic(X, y, family = "gaussian", penalty = "lasso")$selected
#       },
# 
#       glmnet_min = {
#         f <- cv.glmnet(
#           X, y,
#           family = "gaussian",
#           alpha = 1
#         )
#         idx <- which(as.numeric(stats::coef(f, s = "lambda.min")) != 0)
#         idx[idx > 1L] - 1L
#       },
# 
#       glmnet_1se = {
#         f <- cv.glmnet(
#           X, y,
#           family = "gaussian",
#           alpha = 1
#         )
#         idx <- which(as.numeric(stats::coef(f, s = "lambda.1se")) != 0)
#         idx[idx > 1L] - 1L
#       }
#     )
# 
#     setequal(sort(sel), sort(true_idx))
#   })
# 
#   mean(hits)
# }
# 
# s_grid <- 0:s_max
# 
# r_pic  <- vapply(s_grid, recovery, numeric(1L), method = "picreg")
# r_glm  <- vapply(s_grid, recovery, numeric(1L), method = "glmnet_min")
# r_1se  <- vapply(s_grid, recovery, numeric(1L), method = "glmnet_1se")
# 
# plot(
#   s_grid, r_pic,
#   type = "l",
#   lwd = 1.6,
#   lty = 1,
#   ylim = c(0, 1),
#   xlab = "s",
#   ylab = "PESR",
#   main = "picreg vs cv.glmnet - gaussian, lasso",
#   font.main = 1,
#   bty = "l"
# )
# 
# lines(
#   s_grid, r_glm,
#   type = "l",
#   lty = 2,
#   lwd = 1.6
# )
# 
# lines(
#   s_grid, r_1se,
#   type = "l",
#   lty = 3,
#   lwd = 1.6
# )
# 
# legend(
#   "topright",
#   legend = c("picreg", "cv.glmnet (lambda.min)", "cv.glmnet (lambda.1se)"),
#   lty = c(1, 2, 3),
#   lwd = 1.6,
#   bty = "n"
# )

## ----phase-transition-comp-img, echo = FALSE, out.width = "85%", fig.align = "center"----
knitr::include_graphics("phase_transition_comparaison.png")

## ----eval = FALSE-------------------------------------------------------------
# X <- matrix(rnorm(1e4 * 200), 1e4, 200)
# Y <- rnorm(1e4)

## ----eval = FALSE-------------------------------------------------------------
# system.time(pic(X, Y))

## ----echo = FALSE-------------------------------------------------------------
structure(c(1.238, 0.104, 1.404, 0, 0), class = "proc_time",
          .Names = c("user.self", "sys.self", "elapsed",
                     "user.child", "sys.child"))

## ----eval = FALSE-------------------------------------------------------------
# system.time(pic(X, Y, lambda_method = "analytical"))

## ----echo = FALSE-------------------------------------------------------------
structure(c(0.188, 0.023, 0.212, 0, 0), class = "proc_time",
          .Names = c("user.self", "sys.self", "elapsed",
                     "user.child", "sys.child"))

## ----eval = FALSE-------------------------------------------------------------
# system.time(cv.glmnet(X, Y))

## ----echo = FALSE-------------------------------------------------------------
structure(c(0.998, 0.052, 1.054, 0, 0), class = "proc_time",
          .Names = c("user.self", "sys.self", "elapsed",
                     "user.child", "sys.child"))

## ----eval = FALSE-------------------------------------------------------------
# X <- matrix(rnorm(1e4 * 200), 1e4, 200)
# Y <- rbinom(1e4, size = 1, prob = 0.5)

## ----eval = FALSE-------------------------------------------------------------
# system.time(pic(X, Y, family = "binomial"))

## ----echo = FALSE-------------------------------------------------------------
structure(c(1.805, 0.115, 1.941, 0, 0), class = "proc_time",
          .Names = c("user.self", "sys.self", "elapsed",
                     "user.child", "sys.child"))

## ----eval = FALSE-------------------------------------------------------------
# system.time(pic(X, Y, family = "binomial", lambda_method = "analytical"))

## ----echo = FALSE-------------------------------------------------------------
structure(c(0.189, 0.023, 0.213, 0, 0), class = "proc_time",
          .Names = c("user.self", "sys.self", "elapsed",
                     "user.child", "sys.child"))

## ----eval = FALSE-------------------------------------------------------------
# system.time(cv.glmnet(X, Y, family = "binomial"))

## ----echo = FALSE-------------------------------------------------------------
structure(c(3.191, 0.084, 3.377, 0, 0), class = "proc_time",
          .Names = c("user.self", "sys.self", "elapsed",
                     "user.child", "sys.child"))

