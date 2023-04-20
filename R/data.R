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
    x1 = r * cos(theta), x2 = runif(n, a, b)
  ))
}
