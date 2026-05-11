#' Plots for autotune_lasso objects
#'
#' @param object fitted \code{"autotune_lasso"} object 
#' @param max_preds Integer input stating maximum number of predictors to include in the plot. 
#'   Default NULL (then max_preds is calculated internally).
#' @param cumulative Logical input. If TRUE, plots cumulative R-squared; if FALSE, plots 
#'   adjusted R-squared. 
#' @param ... Other graphical parameters to plot
#'
#' @importFrom graphics abline grid legend lines par points text
#' @importFrom stats coef lm
#' @importFrom utils tail
#' @return 
#' Returns a data frame with columns:
#' \itemize{
#'   \item r_squared: Cumulative R-squared values
#'   \item adj_r_squared: Adjusted R-squared values
#'   \item has_nonzero_beta: logical vector whose \eqn{i^{th}} element 
#'   indicates whether the \eqn{i^{th}} ranked predictor has a
#'    non-zero coefficient or not.
#' }
#'
#' @examples
#' 
#' # Fit autotune lasso
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
#' 
#' # Basic diagnostic plot
#' plot(fit)
#' 
#' # Plot adjusted R-squared for first 15 predictors
#' plot(fit, max_preds = 15, cumulative = FALSE)
#' 
#' 
#' @seealso \code{\link{autotune_lasso}}
#'
#' @method plot autotune_lasso
#' @export
plot.autotune_lasso <- function(object,
                                max_preds = NULL,
                                cumulative = TRUE,...) {
  
  if (!inherits(object, "autotune_lasso")) {
    stop("Object must be of class 'autotune_lasso'")
  }

  dots <- list(...)
  cex.main <- if("cex.main" %in% names(dots)) dots$cex.main else 1.4
  cex.lab <- if("cex.lab" %in% names(dots)) dots$cex.lab else 1.6
  cex.axis <- if("cex.axis" %in% names(dots)) dots$cex.axis else 1.6
  dots$cex.main <- NULL
  dots$cex.lab <- NULL
  dots$cex.axis <- NULL

  sorted_predictors <- object$CD.path.details$sorted_predictors
  x <- object$x
  y <- object$y
  n_vars <- object$nvars
  n_obs <- object$nobs
  nnz_beta_count <- sum(object$beta != 0)

  if (is.null(x) || is.null(y)) {
    stop("Original data (x and y) not found in object. Please update your autotune_lasso function to store the original data.")
  }

  adj <- ifelse(n_obs > 10, 10, 1)
  if (is.null(max_preds)) {
    max_preds <- min(n_vars, n_obs - adj, 2 * nnz_beta_count)
  } else {
    max_preds <- min(max_preds, n_vars, n_obs - adj)
  }

  if (max_preds < 1) {
    stop("max_preds must be at least 1")
  }

  r_squared <- numeric(max_preds)
  adj_r_squared <- numeric(max_preds)

  cat(paste("Computing multiple R-squared of",max_preds, "nested linear models"))

  for (i in 1:max_preds) {
    selected_vars <- sorted_predictors[1:i]

    if (is.matrix(x) || is.data.frame(x)) {
      x_subset <- x[, selected_vars, drop = FALSE]
    } else {
      stop("x must be a matrix or data frame")
    }

    fit <- lm(y ~ x_subset)
    df_used <- i + 1

    r_sq <- summary(fit)$r.squared
    r_squared[i] <- r_sq

    if (n_obs > df_used) {
      adj_r_squared[i] <- 1 - (1 - r_sq) * (n_obs - 1) / (n_obs - df_used)
    } else {
      adj_r_squared[i] <- NA
    }
  }

  predictor_indices <- sorted_predictors[1:max_preds]
  has_nonzero_beta <- object$beta[predictor_indices] != 0

  n_predictors <- 1:max_preds

  if (cumulative) {
    y_values <- r_squared
    y_label <- "Cumulative R-squared"
    main_title <- "R-squared vs Number of Predictors (Autotune Lasso Order)"
  } else {
    y_values <- adj_r_squared
    y_label <- "Adjusted R-squared"
    main_title <- "Adjusted R-squared vs Number of Predictors (Autotune Lasso Order)"
  }

 
  point_colors <- ifelse(has_nonzero_beta, "steelblue", "red")

  old_par <- par(mar = c(5.1, 4.5, 4.1, 2.1))  # Default is c(5.1, 4.1, 4.1, 2.1)
  on.exit(par(old_par))

  plot(n_predictors, y_values,
       type = "n",
       xlab = "Number of Predictors",
       ylab = y_label,
       # main = main_title,
       xlim = c(1, max_preds),
       ylim = c(min(y_values, na.rm = TRUE) * 0.95,
                max(y_values, na.rm = TRUE) * 1.05),
       cex.main = cex.main,
       cex.lab = cex.lab,
       cex.axis = cex.axis,
       ...)

  grid(col = "lightgray", lty = "dotted")

  for (i in 1:(max_preds - 1)) {
    if (!is.na(y_values[i]) && !is.na(y_values[i + 1])) {
      lines(c(i, i + 1), c(y_values[i], y_values[i + 1]),
            col = point_colors[i], lwd = 2)
    }
  }

  # Add points on top
  points(n_predictors, y_values, pch = 16, col = point_colors, cex = 1.2)

  # Add a reference line at R² = 0 if relevant
  if (min(y_values, na.rm = TRUE) < 0.1) {
    abline(h = 0, col = "gray", lty = 3)
  }


  legend("bottomright",
         legend = c("Non-zero coefficients", "Zero coefficients"),
         col = c("steelblue", "red"),
         pch = 16,
         bg = "white",
         cex = 1.5)


  final_r2 <- tail(y_values[!is.na(y_values)], 1)
  text(max_preds * 0.8, final_r2,
       paste("Final:", round(final_r2, 3)),
       pos = 3, col = "darkblue", cex = 0.9)

  invisible(data.frame(
    n_predictors = n_predictors,
    predictor_index = predictor_indices,
    has_nonzero_beta = has_nonzero_beta,
    r_squared = r_squared,
    adj_r_squared = adj_r_squared
  ))
}
