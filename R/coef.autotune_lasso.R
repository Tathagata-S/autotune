#' @rdname coef.autotune_lasso
#' @method coef autotune_lasso
#' @export

coef.autotune_lasso <- function(object, intercept = TRUE, path = FALSE, sparse = TRUE, ...) {
  
  if(path){
    coef(object$CD.path.details, intercept = intercept, a0 = object$a0)
  }
  
  if (sparse && !requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' must be installed to return sparse coefficients.")
  }
  
  beta_hat <- numeric(object$nvars + 1)
  
  beta_hat <- c(object$a0, object$beta)
  varnames <- c("(Intercept)", paste0("V", seq_len(object$nvars)))
  names(beta_hat) <- varnames
  
  if (!intercept) {
    beta_hat <- beta_hat[-1]
    varnames <- varnames[-1]
  }
  
  if (!sparse) return(beta_hat)
  
  nz <- which(beta_hat != 0)
  
  sparse_beta <- Matrix::sparseMatrix(
    i = nz,
    j = rep.int(1L, length(nz)),
    x = beta_hat[nz],
    dims = c(length(beta_hat), 1L),
    dimnames = list(varnames, NULL)
  )
  
  sparse_beta
}
