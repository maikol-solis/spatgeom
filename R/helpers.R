

#' @export
estimate_symmetric_reflection <- function(polygon) {
  # true_polygon_coords <- sf::st_coordinates(polygon)
  #
  # polygon <-
  #   polygon - c(min(true_polygon_coords[, "X"]), min(true_polygon_coords[, "Y"]))
  #
  # polygon_coords <- sf::st_coordinates(polygon)
  #
  affine_transformation <-  matrix(c(1, 0, 0, -1), 2, 2)
  cntrd = sf::st_centroid(polygon)

  polygon_reflected <-
    (polygon - cntrd) * affine_transformation + cntrd
  # polygon_reflected <- polygon_reflected +
  #   c(0, min(polygon_coords[, "Y"]) + max(polygon_coords[, "Y"]))
  #   #c(0, 2 * mean(polygon_coords[, "Y"]))
  #   #c(0, 2 * diff(range(polygon_coords[, "Y"])) / 2)
  #
  # polygon_reflected + c(min(true_polygon_coords[, "X"]), min(true_polygon_coords[, "Y"]))
}


#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom magrittr %>%
st_segment <- function(x) {
  sf::st_cast(x, "LINESTRING") %>%
    sf::st_coordinates() %>%
    dplyr::as_tibble() %>%
    dplyr::rename_all(tolower) %>%
    dplyr::group_by(l1) %>%
    dplyr::mutate(x1 = dplyr::lead(x), y1 = dplyr::lead(y)) %>%
    stats::na.omit() %>%
    dplyr::mutate(
      geometry = purrr::pmap(list(x, x1, y, y1), ~ sf::st_linestring(matrix(
        c(..1, ..2, ..3, ..4), 2
      ))),
      geometry = sf::st_as_sfc(geometry)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(geometry)
}


RotMat <- function(angle) {
  matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), nrow = 2, ncol =
           2)
}


num_deriv <- function(y, x) {
  if (length(x) != length(y)) {
    stop('x and y vectors must have equal length')
  }

  n <- length(x)

  # Initialize a vector of length n to enter the derivative approximations
  fdx <- vector(length = n)

  fdx[1] <- (y[2] - y[1]) / (x[2] - x[1])

  # Iterate through the values using the forward differencing method
  for (i in 2:(n - 1)) {
    fdx[i] <- (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
  }

  # For the last value, since we are unable to perform the forward differencing method
  # as only the first n values are known, we use the backward differencing approach
  # instead. Note this will essentially give the same value as the last iteration
  # in the forward differencing method, but it is used as an approximation as we
  # don't have any more information
  fdx[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])

  return(fdx)

}
