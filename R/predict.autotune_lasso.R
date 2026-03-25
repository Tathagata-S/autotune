#' Predict from a autotune_lasso fit
#'
#' @param object An object of class \code{autotune_lasso}.
#' @param newx Numeric matrix of new \code{x} (rows = observations,
#'   columns = predictors) at which predictions are required. This 
#'   argument is not used for \code{type="coefficients")}.
#' @param type Character; one of \code{"response"}, \code{"coefficients"}.
#' @param ... Further arguments passed from the generic.
#'
#' @return
#' If \code{type = "response"}, returns the prediction at the x provided for prediction.
#'
#' If \code{type = "coefficients"}, returns the fitted coefficients.
#'
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
#' predict(fit, newx = x[1:5, , drop = FALSE])
#' predict(fit, type = "coefficients")
#' 
#' predict(fit, type = "coefficients", path = TRUE)
#'
#' @method predict autotune_lasso
#' @export
predict.autotune_lasso <- function(object,
                              newx = NULL,
                              type = c("response", "coefficients"),
                              path = NULL,
                              ...) {
  type <- match.arg(type)
  
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' must be installed to use predict.autotune_lasso().")
  }
  
  # ---- coefficients / nonzero types -----------------------------------------
  if (type == "coefficients") {
    if (is.null(path))
    {
      return(coef(object, intercept = TRUE, path = FALSE, sparse = TRUE))
    }
    if (path) 
    {
      return(coef(object, intercept = TRUE, path = TRUE, sparse = TRUE))
    }
  }
  
  # ---- response type --------------------------------------------------------
  # If no newx provided, return stored fitted values if you keep them
  if (is.null(newx)) {
    stop("newx is NULL and it needs new values for x at which predictions are to be made.")
  }
  
  newx <- as.matrix(newx)
  if (!is.null(object$nvars) && ncol(newx) != object$nvars) {
    stop("newx must have ", object$nvars, " columns.")
  }
  
  m <- nrow(newx)
  yhat <- newx %*% object$beta + rep(object$a0, m)
  return(yhat)
}
