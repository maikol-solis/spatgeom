

estimate_symmetric_reflection <- function(polygon) {
  affine_transformation <- matrix(c(1, 0, 0, -1), 2, 2)
  cntrd <- sf::st_centroid(polygon)
  polygon_reflected <- (polygon - cntrd) * affine_transformation + cntrd
  ## RETURN
  polygon_reflected
}


st_segment <- function(x) {
  l1 <- x1 <- y <- y1 <- geometry <- NULL
  segment <- sf::st_cast(x, "LINESTRING")
  segment <- sf::st_coordinates(segment)
  segment <- dplyr::as_tibble(segment)
  segment <- dplyr::rename_with(segment, .fn = tolower)
  segment <- dplyr::group_by(segment, l1)
  segment <- dplyr::mutate(segment, x1 = dplyr::lead(x), y1 = dplyr::lead(y))
  segment <- stats::na.omit(segment)
  segment <- dplyr::mutate(segment,
    geometry = purrr::pmap(
      list(x, x1, y, y1),
      ~ sf::st_linestring(matrix(c(..1, ..2, ..3, ..4), 2))
    ),
    geometry = sf::st_as_sfc(geometry)
  )
  segment <- dplyr::ungroup(segment)
  segment <- dplyr::select(segment, geometry)
  ## RETURN
  segment
}


rot_mat <- function(angle) {
  matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)),
    nrow = 2, ncol = 2
  )
}


num_deriv <- function(y, x) {
  if (length(x) != length(y)) {
    stop("x and y vectors must have equal length")
  }

  n <- length(x)

  # Initialize a vector of length n to enter the derivative approximations
  fdx <- vector(length = n)

  fdx[1] <- (y[2] - y[1]) / (x[2] - x[1])

  # Iterate through the values using the forward differencing method
  for (i in 2:(n - 1)) {
    fdx[i] <- (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
  }

  # For the last value, since we are unable to perform the forward differencing
  # method as only the first n values are known, we use the backward
  # differencing approach instead. Note this will essentially give the same
  # value as the last iteration in the forward differencing method, but it is
  # used as an approximation as we don't have any more information
  fdx[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])

  return(fdx)
}
