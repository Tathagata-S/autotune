#' @rdname coef.autotune_lasso
#' @method coef autotune_lasso_path
#' @export
coef.autotune_lasso_path <- function(object,
                                intercept = TRUE,
                                a0 = NULL, ...) {
  
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' must be installed to return sparse coefficients.")
  }
  
  beta_mat <- object$beta_matrix  # (k+1) x p
  if (!is.matrix(beta_mat)) beta_mat <- as.matrix(beta_mat)
  
  beta_sp <- Matrix::Matrix(t(object$beta_mat), sparse = TRUE)  # p x (k+1)
  
  k_plus_1 <- ncol(beta_sp)
  k <- k_plus_1 - 1
  p <- nrow(beta_sp)
  
  soln_names <- c(paste0("lambda_", seq_len(k)), "final_lambda")
  colnames(beta_sp) <- soln_names
  varnames <- paste0("V", seq_len(p))
  rownames(beta_sp) <- varnames
  if(intercept) {
    beta_sp <- rbind(rep(a0,k_plus_1), beta_sp)
    rownames(beta_sp) <- c("(Intercept)", varnames)
  }
  beta_sp
}