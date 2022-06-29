#' print a \code{spatgeom} object
#'
#' Print method for objects of class \code{spatgeom}.
#'
#' @param x an object of class \code{spatgeom}
#' @param return_table  if \code{TRUE}, returns a data frame with the
#'   estimated values. Otherwise, print the data frame in console. Defaults to
#'   \code{FALSE}
#' @param ... further arguments passed to the \code{plot} function
#'
#' @return Print the estimate given by \code{\link{alphastats}}.
#' @export
#'
#' @examples
#' n <- 30
#' a <- -1
#' b <- 1
#' theta <- runif(n, 0, 2 * pi)
#' r <- (sqrt(runif(n))) * (0.5) + 0.5
#' X1 <- r * cos(theta)
#' X2 <- runif(n, a, b)
#' Y <- data.frame(Y = r * sin(theta))
#' X <- data.frame(X1, X2)
#'
#' estimation <- alphastats(y = Y, x = X)
#'
#' print(estimation)
#'
#' @export
#'
print.spatgeom <- function(x, return_table = FALSE, ...) {
  out <- lapply(
    X = x$results,
    FUN = function(xx) {
      cbind(
        variable_name = xx$variable_name,
        mean_n = xx$mean_n,
        intensity = xx$intensity,
        xx$data_frame_triangles
      )
    }
  )

  out <- do.call(rbind, out)
  out <- as.data.frame(out)


  if (return_table == TRUE) {
    return(out)
  }

  variable_name <- mean_n <-
    intensity <- geom_corr <- alpha <- NULL
  o <- dplyr::group_by(.data = out, variable_name)
  o <- dplyr::summarise(
    .data = o,
    mean_n = min(mean_n),
    intensity = dplyr::first(cut(intensity, breaks = 2)),
    alpha = dplyr::first(cut(alpha, breaks = 2)),
    geom_corr = dplyr::first(cut(geom_corr, breaks = 2))
  )

  cat("\nCall:\n", deparse(x[["call"]]), "\n", sep = "")
  cat("\nNumber of variables:", ncol(x[["x"]]), "\n")
  cat("\nNumber of observations:", nrow(x[["y"]]), "\n")
  print(o)
}
