#' Extract coefficients for autotune_lasso objects
#'
#' This function extracts the regression coefficients estimated by an autotune
#' lasso model, using the stored \code{"autotune_lasso"} object
#'
#'
#' @param object An \code{"autotune_lasso"} object (for \code{coef.autotune_lasso})
#'   or \code{"autotune_lasso_path"} object (for \code{coef.autotune_lasso_path}).
#' @param intercept Logical; include the intercept term? If \code{FALSE}, the intercept
#'   row is omitted.
#' @param path Logical; if \code{TRUE}, return the full coefficient path across 
#' different lambdas.
#' @param sparse Logical; controls the output format of \code{coef.autotune_lasso()} when
#'   \code{full = FALSE}. If \code{TRUE}, returns a sparse \code{dgCMatrix} (glmnet-style);
#'   if \code{FALSE}, returns a base numeric vector.
#' @param a0 Intercept of the fitted autotune regression. Used only when \code{intercept = TRUE}.
#' @param ... Other arguments to coef.
#'
#' @return
#' \itemize{
#'   \item \code{coef.autotune_lasso(full = FALSE, sparse = TRUE)} returns a sparse \code{dgCMatrix}
#'     of dimension (\code{nvars} + 1) \eqn{\times 1} (or \code{nvars}  \eqn{\times 1} if \code{intercept = FALSE}).
#'   \item \code{coef.autotune_lasso(full = FALSE, sparse = FALSE)} returns a named numeric vector.
#'   \item \code{coef.autotune_lasso(full = TRUE)} calls \code{coef.autotune_lasso_path()} which return a sparse
#'     \code{dgCMatrix} of dimension (\code{nvars}  + 1) \eqn{\times} (k + 1) (or \code{nvars}  \eqn{\times} (k + 1) if
#'     \code{intercept = FALSE}), with columns named \code{lambda_1, ..., lambda_k, final_lambda}.
#' }
#'
#' @examples
#' set.seed(10)
#' n = 300
#' p = 500
#' s = 10
#' beta = c(rep(1, s), rep(0, p - s))
#' x = matrix(rnorm(n * p), ncol = p)
#' # Maunal sigma allocation
#' # y = x %*% beta + rnorm(n, sd = 1)
#' # Dynamic sigma allocation with snr specified
#' snr = 2
#' y = x %*% beta + rnorm(n, sd = sqrt(var(x%*%beta)/snr))
#' fit <- autotune_lasso(x, y)
#' coef(fit)                    # final solution
#' coef(fit, full = TRUE)        # full path (sparse)
#'
#' @importFrom Matrix sparseMatrix Matrix t rBind
#' @name coef.autotune_lasso
NULL