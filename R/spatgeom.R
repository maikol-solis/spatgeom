#' @title Geometric Spatial Point Pattern Analysis
#'
#' @description Function to estimate the geometric correlation between
#'   variables.
#'
#' @param x numeric matrix or data.frame of covariables.
#' @param y numeric vector of responses in a model.
#' @param scale boolean to make the estimations with scaled variables. Default
#'   \code{FALSE}.
#' @param nalphas a single number for the number of alphas generated between the
#'   minimum and maximum edge distance on the Delanauy triangulation.
#' @param envelope boolean to determine if the Monte-Carlo is estimated. Default
#'   \code{FALSE}.
#' @param mc_cores an integer to determine how many parallel process should be
#'   run. Default \code{mc_core=1}.
#'
#'
#' @return A list of class \code{spatgeom} with  the following elements:
#'
#' \describe{
#' \item{\strong{call}}{The function call.}
#'
#' \item{\strong{x}}{\code{x} input.}
#'
#' \item{\strong{y}}{\code{y} output.}
#'
#' \item{\strong{results}}{A list of size \code{ncol(x)} corresponding to each
#' column of \code{x}. Each element of the list has:
#' \describe{
#'
#' \item{\strong{triangles}}{a data frame of class \code{sfc} (see
#' [`sf::st_sf()`])with columns \code{geometry}, \code{segments},
#' \code{max_length} and \code{alpha}. The data.frame contains the whole
#' Delanauy triangulation for the corresponding column of \code{x} and \code{y}.
#' The \code{segments} column are the segments of each individual triangle and
#' \code{max_length} is the maximum length of them.}
#'
#' \item{\strong{geom_indices}}{a data frame with columns \code{alpha}
#'  and \code{geom_corr}. The \code{alpha} column is a numeric vector of size
#'  \code{nalphas} from the minimum to the maximum distance between points
#'  estimated in the data. The \code{geom_corr} column is the value \code{1 -
#'  (alpha shape Area)/(containing box Area).}}
#'
#' \item{\strong{intensity}}{the intensity estimated for the corresponding
#' column of \code{x} and \code{y}.}
#'
#' \item{\strong{mean_n}}{the mean number of points in the point process.}
#'
#' \item{\strong{envelope_data}}{a data frame in tidy format with 40 runs of a
#' CSR process, if \code{envelope=TRUE}, The CSR is created by generating
#' \emph{n} uniform points in the plane, where \emph{n} is drawn from Poisson
#' distribution with parameter \code{mean_n}.
#'
#' }}}}
#'
#' @references
#'
#' Hernández, A.J., Solís, M. Geometric goodness of fit measure to detect
#' patterns in data point clouds. Comput Stat (2022).
#' https://doi.org/10.1007/s00180-022-01244-1
#'
#' @examples
#' xy <- donut_data(n = 30, a = -1, b = 1, theta = 2 * pi)
#' estimation <- spatgeom(y = xy[, 1], x = xy[, -1])
#'
#' # If you want to estimate the envelope, you can use the envelope argument to
#' # TRUE. This will take a while to run.
#' \dontrun{
#' estimation_with_envelope <- spatgeom(
#'   y = xy[, 1], x = xy[, -1],
#'   envelope = TRUE
#' )
#' }
#' @export


spatgeom <- function(x, y,
                     scale = FALSE,
                     nalphas = 100,
                     envelope = FALSE,
                     mc_cores = 1) {
  if (missing(y)) {
    message("Running with only x")
  } else {
    message("Running with x and y")
    spatgeom_xy(x, y,
      scale = scale,
      nalphas = nalphas,
      envelope = envelope,
      mc_cores = mc_cores
    )
  }
}



spatgeom_xy <- function(x, y,
                        scale = FALSE,
                        nalphas = 100,
                        envelope = FALSE,
                        mc_cores = 2) {
  x <- as.data.frame(x)
  y <- as.data.frame(y)

  ans <- list()
  ans[["call"]] <- match.call()
  ans[["x"]] <- x
  ans[["y"]] <- y

  message("Index estimation")

  out_list <- parallel::mclapply(
    mc.cores = mc_cores,
    X = seq_len(ncol(x)),
    FUN = function(i) {
      message(paste0("Estimating R2 Geom for variable = ", i))
      estimate_curves(
        x1 = x[, i],
        x2 = y[, 1],
        scale = scale,
        nalphas = nalphas
      )
    }
  )

  out_list <- lapply(
    X = seq_len(ncol(x)),
    FUN = function(i) {
      append(out_list[[i]], list(variable_name = colnames(x)[i]))
    }
  )


  if (envelope == TRUE) {
    out_list <- estimate_envelope(
      triangles_list = out_list,
      x = x,
      y = y,
      scale = scale,
      nalphas = nalphas,
      mc_cores = mc_cores
    )
  }

  ans[["results"]] <- out_list
  class(ans) <- "spatgeom"
  return(ans)
}




spatgeom_x <- function(x, ...) {

}




estimate_curves <- function(x1, x2, scale, nalphas, intensity = NULL) {
  if (scale) {
    pts <- sf::st_cast(
      sf::st_sfc(
        sf::st_multipoint(scales::rescale(cbind(x1, x2)))
      ), "POINT"
    )
  } else {
    pts <- sf::st_cast(
      sf::st_sfc(
        sf::st_multipoint(cbind(x1, x2))
      ), "POINT"
    )
  }
  bb <- sf::st_make_grid(pts, n = 1)

  if (is.null(intensity)) {
    intensity <- length(pts) / sf::st_area(bb)
  }

  v2 <- sf::st_triangulate(sf::st_geometrycollection(pts))

  polygons <- sf::st_collection_extract(sf::st_sfc(v2))
  linestrings <- sf::st_cast(polygons, "LINESTRING")
  linestrings_splitted <- lwgeom::st_split(linestrings, pts)
  max_length <-
    sapply(linestrings_splitted, function(x) {
      max(sf::st_length(sf::st_cast(sf::st_sfc(x))))
    })

  triangles <-
    sf::st_sf(
      list(
        geometry = polygons,
        segments = linestrings_splitted,
        max_length = max_length,
        alpha = max_length / 2
      )
    )
  triangles <- triangles[order(triangles$alpha), ]

  geom_corr <- NULL

  out <- lapply(
    X = seq_along(triangles$alpha),
    FUN = function(k) {
      alpha_shape <- triangles[1:k, ]

      if (nrow(alpha_shape) > 0) {
        poly_union <- sf::st_union(alpha_shape$geometry)

        ## Geometric R2 index
        geom_corr <-
          1 - sf::st_area(poly_union) / sf::st_area(bb)
      } else {
        geom_corr <- 1
      }
      return(list(geom_corr = geom_corr))
    }
  )

  geom_corr <- sapply(out, function(x) {
    x$geom_corr
  })

  geom_indices <- data.frame(
    alpha = triangles$alpha,
    geom_corr = geom_corr
  )

  return(
    list(
      triangles = triangles,
      geom_indices = geom_indices,
      intensity = intensity,
      mean_n = sf::st_area(bb) * intensity
    )
  )
}

estimate_envelope <- function(triangles_list,
                              x,
                              y,
                              scale,
                              nalphas,
                              mc_cores = 2) {
  for (i in seq_len(ncol(x))) {
    message(paste0("Estimating envelope for variable = ", i))
    envelope_data <-
      data.frame(
        y = numeric(),
        x = numeric(),
        nsim = numeric()
      )

    envelope_data <- parallel::mclapply(
      mc.cores = mc_cores,
      X = seq_len(40),
      FUN = function(k) {
        n <- stats::rpois(1, lambda = triangles_list[[i]]$mean_n)
        x <- stats::runif(n, min = min(x[, i]), max = max(x[, i]))
        y <- stats::runif(n, min = min(y[, 1]), max = max(y[, 1]))
        enve <-
          estimate_curves(
            x1 = x,
            x2 = y,
            scale = scale,
            nalphas = nalphas,
            intensity = triangles_list[[i]]$intensity
          )
        enve_approx <-
          stats::approx(
            x = enve$geom_indices$alpha,
            y = enve$geom_indices$geom_corr,
            xout = triangles_list[[i]]$geom_indices$alpha
          )
        data.frame(enve_approx, nsim = k)
      }
    )
    envelope_data <- do.call("rbind", envelope_data)
    triangles_list[[i]]$envelope_data <- envelope_data
  }
  return(triangles_list)
}
