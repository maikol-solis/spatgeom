#' @title Donut example
#'
#' @description Generate data points with the shape of a donut.
#' @param n Number of data points.
#' @param a Lower bound of the second variable.
#' @param b Upper bound of the second variable.
#' @param theta Angle of the donut.
#'
#'
#' @return A data frame with three variables. Variable 'y' is the response,
#'   variable 'x1' makes the donut shape with 'y', and 'x2' is a uniform random
#'   variable between a and b. '
#' @examples
#' xy <- donut_data(n = 30, a = -1, b = 1, theta = 2 * pi)
#' @export
donut_data <- function(n, a, b, theta) {
  theta <- stats::runif(n, 0, theta)
  r <- (sqrt(stats::runif(n))) * (0.5) + 0.5
  return(data.frame(
    y = r * sin(theta),
    x1 = r * cos(theta), x2 = stats::runif(n, a, b)
  ))
}


#' @title Linear example
#' @description Generate data points with a linear relationship.
#' @param n Number of data points.
#' @param a,b Lower and upper bound of the uniform distribution.
#' @return A data frame with three variables. Variable 'y = 0.6 * x1 + 0.3 * x2
#'   + 0.1 * x3' is the response, and 'x1', 'x2' and 'x3' are uniform random
#'   variables between a and b.
#' @examples
#' xy <- linear_data(n = 30, a = -1, b = 1)
#' @export


linear_data <- function(n = 100, a = -3, b = 3) {
  x1 <- stats::runif(n, a, b)
  x2 <- stats::runif(n, a, b)
  xnoise <- matrix(stats::runif(n * 1, a, b), nrow = n)
  x <- data.frame(x1, x2, X3 = xnoise)
  y <- data.frame(y = 0.6 * x1 + 0.3 * x2 + 0.1 * xnoise)
  return(cbind(y, x))
}
