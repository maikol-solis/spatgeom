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
#' @return Print the estimate given by \code{\link{spatgeom}}.
#'
#' @examples
#'
#' xy <- donut_data(n = 30, a = -1, b = 1, theta = 2 * pi)
#'
#' # Basic print — shows alpha and geom_survival ranges per variable:
#' estimation <- spatgeom(y = xy[, 1], x = xy[, -1])
#' print(estimation)
#'
#' \donttest{
#' # Print a spatgeom object that includes hypothesis testing results:
#' est_ht <- spatgeom(
#'   y = xy[, 1], x = xy[, -1],
#'   hypothesis_testing = TRUE, method = "MAD"
#' )
#' print(est_ht)
#'
#' # Return the underlying data frame instead of printing:
#' tbl <- print(est_ht, return_table = TRUE)
#' head(tbl)
#' }
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
        xx$geom_indices
      )
    }
  )

  out <- do.call(rbind, out)
  out <- as.data.frame(out)

  if (return_table == TRUE) {
    return(out)
  }

  variable_name <- mean_n <-
    intensity <- geom_survival <- alpha <- NULL
  o <- dplyr::group_by(.data = out, variable_name)
  o <- dplyr::summarise(
    .data = o,
    mean_n = min(mean_n),
    intensity = min(intensity),
    alpha = dplyr::first(cut(alpha, breaks = 2)),
    geom_survival = dplyr::first(cut(geom_survival, breaks = 2))
  )

  cat("\nCall:\n", deparse(x[["call"]]), "\n", sep = "")
  cat("\nNumber of variables:", ncol(x[["x"]]), "\n")
  cat("\nNumber of observations:", nrow(x[["y"]]), "\n")
  print(o)
}
