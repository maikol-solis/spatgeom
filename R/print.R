#' print \code{topsa} objects
#'
#' Print method for objects of class \code{topsa}.
#'
#' @param topsaObj an object of class \code{topsa}
#' @param only.return.table  if \code{TRUE}, returns a data frame with the
#'   estimated values. Otherwise, print the data frame in console. Defaults to
#'   \code{FALSE}
#' @param ... further arguments passed to the \code{plot} function
#'
#' @return Print the threshold used, the box area, manifold embedding area, geometric
#' correlation index and symmetric sensitivity index for and object of class
#' \code{topsa}.
#' @export
#'
#' @examples
#'
#' ishigami.fun <- function(X) {
#' A <- 7
#' B <- 0.1
#' sin(X[, 1]) + A * sin(X[, 2])^2 + B * X[, 3]^4 * sin(X[, 1])
#' }
#'
#' X <- matrix(runif(3*50, -pi, pi), ncol = 3)
#' Y <- ishigami.fun(X)
#'
#' estimation <- topsa(Ydat = Y, Xdat = X)
#'
#' print(estimation)
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
