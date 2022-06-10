#' Geometric Sensitivity Analysis
#' @param y numeric vector of responses in a model.
#' @param x numeric matrix or data.frame of covariables.
#' @return A list of class \code{topsa} with the following elements:
#'
#' \describe{
#' \item{\strong{call}}{The function call.}
#' \item{\strong{Xdat}}{\code{X} input.}
#' \item{\strong{Ydat}}{\code{Y} output.}
#' \item{\strong{dimension}}{dimension to estimate the homology order.}
#' \item{\strong{threshold}}{cutoff level for the radius or area.}
#' \item{\strong{results}}{A list for each variable with:
#' \describe{
#' \item{\strong{threshold}}{threshold used to limit the radius or area.}
#' \item{\strong{Manifold_Area}}{geometrical area of the estimated manifold.}
#' \item{\strong{Box.Area}}{geometrical area of the estimated manifold.}
#' \item{\strong{Geometric.R2}}{geometric correlation between each
#' \code{x} and \code{y}.}
#' \item{\strong{Geometric.Sensitivity}}{symmetric sensitivity
#' index of each estimated manifold.}
#' \item{\strong{manifold_plot}}{a \code{sf}
#' object with the estimated manifold.}
#' } } }
#' @examples
#'
#' ishigami.fun <- function(X) {
#'   A <- 7
#'   B <- 0.1
#'   sin(X[, 1]) + A * sin(X[, 2])^2 + B * X[, 3]^4 * sin(X[, 1])
#' }
#'
#' X <- matrix(runif(3 * 50, -pi, pi), ncol = 3)
#' Y <- ishigami.fun(X)
#' estimation <- topsa(y = Y, x = X, method = "Alpha")
#' @importFrom methods is
#' @importFrom magrittr %>%
#' @export

alphastats <- function(y,
                       x,
                       scale = FALSE,
                       nalphas = 100,
                       envelope = TRUE) {
  x <- as.data.frame(x)
  y <- as.data.frame(y)

  ANS <- list()
  ANS[["call"]] <- match.call()
  ANS[["x"]] <- x
  ANS[["y"]] <- y

  # Xr <- matrix()
  # Yr <- matrix()
  # l <- lapply(seq_along(x), function(k) {
  #   scales::rescale(cbind(x[, k], y[, 1]))
  # })
  #
  # lx <- lapply(l, function(x)
  #   x[, 1])
  #
  # Xr <- as.data.frame(do.call("cbind", lx))
  # Yr <- as.data.frame(sapply(y, scales::rescale))
  # ANS[['Xr']] <- Xr
  # ANS[['Yr']] <- Yr
  # ANS[["Xr"]] <- as.data.frame(lapply(x, scales::rescale))
  # ANS[["Yr"]] <- as.data.frame(lapply(y, scales::rescale))
  # ANS[["angle"]] <- angle
  #


  # if (length(threshold.radius) == 1) {
  #   threshold.radius <- rep(threshold.radius, ncol(x))
  # } else if (length(threshold.radius) < ncol(x)) {
  #   stop("Please provide a numeric threshold vector of size 1 or ncol(x)")
  # }




  message("Index estimation")



  out_list <- parallel::mclapply(
    mc.cores = 6,
    X = seq_len(ncol(x)),
    FUN = function(i) {
      message(paste0("Estimating R2 Geom for variable = ", i))
      estimate_curves(
        x = x[, i],
        y = y[, 1],
        scale = scale,
        nalphas = nalphas
      )
    }
  )

  if (envelope == TRUE) {
    for (i in seq_len(ncol(x))) {
      message(paste0("Estimating envelope for variable = ", i))
      envelope_data <-
        data.frame(
          y = numeric(),
          x = numeric(),
          nsim = numeric()
        )

      envelope_data <- parallel::mclapply(
        mc.cores = 6,
        X = seq_len(40),
        FUN = function(k) {
          n <- rpois(1, lambda = out_list[[i]]$mean_n)
          x <- runif(n, min = min(x[, i]), max = max(x[, i]))
          y <- runif(n, min = min(y[, 1]), max = max(y[, 1]))
          enve <-
            estimate_curves(
              x = x,
              y = y,
              scale = scale,
              nalphas = nalphas,
              intensity = out_list[[i]]$intensity
            )
          enve_approx <-
            approx(
              enve$data_frame_triangles$alpha,
              y = enve$data_frame_triangles$geom_corr,
              xout = out_list[[i]]$data_frame_triangles$alpha
            )
          data.frame(enve_approx, nsim = k)
        }
      )
      envelope_data <- do.call("rbind", envelope_data)
      out_list[[i]]$envelope_data <- envelope_data
    }
  }

  ANS[["results"]] <- out_list
  class(ANS) <- "geomsensitivity"
  return(ANS)
}



estimate_curves <- function(x, y, scale, nalphas, intensity = NULL) {
  if (scale) {
    pts <-
      sf::st_cast(sf::st_sfc(sf::st_multipoint(scales::rescale(cbind(
        x, y
      )))), "POINT")
  } else {
    pts <-
      sf::st_cast(sf::st_sfc(sf::st_multipoint(cbind(x, y))), "POINT")
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
    st_sf(
      list(
        geometry = polygons,
        segments = linestrings_splitted,
        max_length = max_length,
        alpha = max_length / 2
      )
    )
  triangles <- triangles[order(triangles$alpha), ]

  geom_corr <- geom_sens <- geom_sens2 <- NULL
  # d <- dist(cbind(X, Y))
  d_min <- min(triangles$alpha)
  d_max <- max(triangles$alpha)
  alpha_seq <- seq(d_min, d_max * 1.1, length.out = nalphas)

  out <- lapply(
    X = alpha_seq,
    FUN = function(s) {
      alpha_shape <- subset(triangles, alpha <= s)

      # alpha_shape <- triangles %>%
      #   dplyr::filter(max_length < 2 * alpha_seq[s])
      if (nrow(alpha_shape) > 0) {
        poly_union <- sf::st_union(alpha_shape$geometry)
        poly_reflection <- estimate_symmetric_reflection(poly_union)


        poly_sym_difference <-
          sf::st_sym_difference(poly_union, poly_reflection)
        poly_sym_difference_bb <-
          sf::st_sym_difference(bb, poly_sym_difference)


        geom_corr <-
          1 - sf::st_area(poly_union) / sf::st_area(bb)
        geom_sens <-
          sf::st_area(poly_sym_difference) / (2 * sf::st_area(poly_union))
        geom_sens2 <- sf::st_area(poly_sym_difference_bb) /
          sf::st_area(bb)
      } else {
        geom_corr <- 1
        geom_sens <- 1
      }
      return(list(geom_corr = geom_corr, geom_sens = geom_sens))
    }
  )

  geom_corr <- sapply(out, function(x) {
    x$geom_corr
  })
  geom_sens <- sapply(out, function(x) {
    x$geom_sens
  })

  data_frame_triangles <- data.frame(
    alpha = alpha_seq,
    geom_corr,
    geom_sens
  )

  return(
    list(
      triangles = triangles,
      data_frame_triangles = data_frame_triangles,
      intensity = intensity,
      mean_n = sf::st_area(bb) * intensity
    )
  )
}
