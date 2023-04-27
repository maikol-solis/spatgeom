#' @title Geometric Spatial Point Pattern Analysis
#'
#' @description Function to estimate the geometric correlation between
#'   variables.
#'
#' @param x numeric matrix or data.frame of covariables.
#' @param y numeric vector of responses in a model.
#' @param scale_pts boolean to make the estimations with scaled variables.
#'   Default \code{FALSE}.
#' @param nalphas number of alphas generated for creating the geometric measure
#'   of fit index. Default 100.
#' @param envelope boolean to determine if the Monte-Carlo should be estimated.
#'   Default \code{FALSE}.
#' @param domain_type character with the type of domain to use. It can be either
#'   "bounding-box" or "convex-hull". Default "bounding-box".
#' @param mc_cores integer with the number of parallel process to run (if
#'   available). Default \code{1}.
#' @param hypothesis_testing logical. If \code{TRUE}, performs a CSR hypothesis
#'   test using a global envelope approach after estimating the geometric
#'   survival curves. Automatically sets \code{envelope = TRUE}. Default
#'   \code{FALSE}.
#' @param significance_level a numeric significance level passed to
#'   \code{\link{csr_test}}. Default \code{0.05}.
#' @param r numeric scaling parameter for the theoretical CSR curve used in
#'   \code{\link{csr_test}}: \eqn{\exp(-\lambda \pi (\alpha r)^2)}. Default
#'   \code{0.5}.
#' @param method character, one of \code{"MAD"} (Maximum Absolute Deviation)
#'   or \code{"DCLF"} (Diggle-Cressie-Loosmore-Ford). Specifies the global
#'   envelope test statistic passed to \code{\link{csr_test}}. Default
#'   \code{"MAD"}.
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
#' \item{\strong{geom_indices}}{a data frame with columns \code{alpha} and
#'  \code{geom_survival}. The \code{alpha} column is a numeric vector of size
#'  \code{nalphas} from the minimum to the maximum distance between points
#'  estimated in the data. The \code{geom_survival} column is the value \code{1
#'  - (alpha shape Area)/(containing box Area).}}
#'
#' \item{\strong{intensity}}{the intensity estimated for the corresponding
#' column of \code{x} and \code{y}.}
#'
#' \item{\strong{mean_n}}{the mean number of points in the point process.}
#'
#' \item{\strong{envelope_data}}{a data frame in tidy format with 40 runs of a
#' CSR process, if \code{envelope=TRUE}, The CSR is created by generating
#' \emph{n} uniform points in the plane, where \emph{n} is drawn from Poisson
#' distribution with parameter \code{mean_n}.}
#'
#' }}
#'
#' \item{\strong{hypothesis_testing_results}}{Only present when
#' \code{hypothesis_testing = TRUE}. A list returned by
#' \code{\link{csr_test}} with elements \code{hypothesis_testing_df} (a data
#' frame of test statistics and confidence bands for all variables) and
#' \code{details} (per-variable summary tables of the test results).}
#' }
#'
#' @references
#'
#' Hernández, A.J., Solís, M. Geometric goodness of fit measure to detect
#' patterns in data point clouds. Comput Stat (2022).
#' https://doi.org/10.1007/s00180-022-01244-1
#'
#' @examples
#'
#' xy <- donut_data(n = 30, a = -1, b = 1, theta = 2 * pi)
#'
#' # Basic usage: estimate the geometric survival curves only.
#' estimation <- spatgeom(y = xy[, 1], x = xy[, -1])
#' print(estimation)
#'
#' \donttest{
#' # With Monte Carlo envelope (takes a few seconds):
#' est_env <- spatgeom(y = xy[, 1], x = xy[, -1], envelope = TRUE)
#' plot_curve(est_env, type = "curve")
#'
#' # With integrated CSR hypothesis testing using the MAD statistic:
#' est_ht <- spatgeom(
#'   y = xy[, 1], x = xy[, -1],
#'   hypothesis_testing = TRUE,
#'   method = "MAD"
#' )
#' print(est_ht)
#' plot_curve(est_ht, type = "curve")
#'
#' # Inspect the per-variable test results:
#' est_ht$hypothesis_testing_results$details
#' }
#' @export

spatgeom <- function(
  x,
  y,
  scale_pts = FALSE,
  nalphas = 100,
  envelope = FALSE,
  domain_type = c("bounding-box", "convex-hull"),
  hypothesis_testing = FALSE,
  significance_level = 0.05,
  mc_cores = 1,
  r = 0.5,
  method = c("MAD", "DCLF")
) {
  if (missing(y)) {
    stop("The argument y is required when running spatgeom.")
  } else {
    message("Running spatgeom with both x and y.")
    domain_type <- domain_type[1]
    # If hypothesis testing is requested, force envelope simulation.
    if (hypothesis_testing && !envelope) {
      envelope <- TRUE
      message(
        "Hypothesis testing requires envelope simulation.",
        " Setting envelope = TRUE."
      )
    }

    # Compute the basic spatgeom object (without hypothesis testing)
    spatgeom_obj <- spatgeom_xy(
      x,
      y,
      scale_pts = scale_pts,
      nalphas = nalphas,
      envelope = envelope,
      domain_type = domain_type,
      mc_cores = mc_cores
    )

    # If requested, run the CSR hypothesis test and add its output to the
    # object.
    if (hypothesis_testing) {
      spatgeom_obj$hypothesis_testing_results <- csr_test(
        spatgeom_obj,
        significance_level = significance_level,
        r = r,
        method = method
      )
    }

    return(spatgeom_obj)
  }
}

# --- Internal function that computes the geometric survival curves ---
spatgeom_xy <- function(
  x,
  y,
  scale_pts = FALSE,
  nalphas = 100,
  envelope = FALSE,
  domain_type = c("bounding-box", "convex-hull"),
  mc_cores = 2
) {
  x <- as.data.frame(x)
  y <- as.data.frame(y)
  domain_type <- domain_type[1]

  ans <- list()
  ans[["call"]] <- match.call()
  ans[["x"]] <- x
  ans[["y"]] <- y

  message("Estimating geometric survival of empty space...")
  out_list <- parallel::mclapply(
    X = seq_len(ncol(x)),
    FUN = function(i) {
      message(paste0("Processing variable ", i, " (", colnames(x)[i], ")"))
      estimate_curves(
        x1 = x[, i],
        x2 = y[, 1],
        scale_pts = scale_pts,
        nalphas = nalphas,
        domain_type = domain_type
      )
    },
    mc.cores = mc_cores
  )

  # Attach the variable names.
  out_list <- lapply(
    X = seq_len(ncol(x)),
    FUN = function(i) {
      append(out_list[[i]], list(variable_name = colnames(x)[i]))
    }
  )

  if (envelope == TRUE) {
    out_list <- estimate_envelope(
      spatgeom_obj = out_list,
      x = x,
      y = y,
      scale_pts = scale_pts,
      nalphas = nalphas,
      domain_type = domain_type,
      mc_cores = mc_cores
    )
  }

  ans[["results"]] <- out_list
  class(ans) <- "spatgeom"
  return(ans)
}

# --- Estimate the geometric survival curves using alpha-shapes ---
estimate_curves <- function(
  x1,
  x2,
  scale_pts,
  nalphas,
  intensity = NULL,
  domain_type
) {
  # Scale points if requested.
  coords <- if (scale_pts) {
    scales::rescale(cbind(x1, x2))
  } else {
    cbind(x1, x2)
  }
  pts <- sf::st_multipoint(coords)
  pts <- sf::st_sfc(pts)
  pts <- sf::st_cast(pts, "POINT")

  # Step 1
  # Create the domain.
  if (domain_type == "convex-hull") {
    bb <- sf::st_convex_hull(sf::st_union(pts))
  } else if (domain_type == "bounding-box") {
    bb <- sf::st_make_grid(pts, n = 1)
  } else {
    stop("domain_type must be either 'convex-hull' or 'bounding-box'")
  }

  # Step 2
  # Estimate intensity if not provided
  if (is.null(intensity)) {
    intensity <- length(pts) / sf::st_area(bb)
  }

  # Step 3
  # Alpha-shape construction

  # Triangulation
  delaunauy <- sf::st_triangulate(sf::st_geometrycollection(pts))

  # Extract polygons and linestrings
  polygons <- sf::st_collection_extract(sf::st_sfc(delaunauy))
  linestrings <- sf::st_cast(polygons, "LINESTRING")
  linestrings_splitted <- lwgeom::st_split(linestrings, pts)

  # Vectorized computation of max_length
  max_length <- sapply(
    linestrings_splitted,
    function(x) {
      max(sf::st_length(sf::st_cast(sf::st_sfc(x))))
    }
  )

  # Construction of the data frame
  triangles <- data.frame(
    geometry = polygons,
    max_length = max_length,
    alpha = max_length / 2
  )
  # Order by alpha
  triangles <- triangles[order(triangles$alpha), ]

  ## Create a sequence indices with step = number of triangles / nalphas
  number_triangles <- nrow(triangles)
  step_seq <- ifelse(
    nalphas < number_triangles,
    ceiling(number_triangles / nalphas),
    1
  )
  idx_triangles <- unique(
    c(seq(1, number_triangles, by = step_seq), number_triangles)
  )

  # Precompute all individual geometries outside the loop
  individual_geometries <- lapply(seq_along(idx_triangles), function(k) {
    if (k == 1) {
      triangles[1, ]$geometry
    } else {
      sf::st_union(
        triangles[idx_triangles[k - 1]:(idx_triangles[k] - 1), ]$geometry
      )
    }
  })

  # Compute the cumulative union
  cumulative_union <- Reduce(
    sf::st_union,
    individual_geometries,
    accumulate = TRUE
  )

  # Step 4
  # Calculate the geometric indices

  # Calculate areas and print them
  areas <- sapply(cumulative_union, sf::st_area)
  geom_survival <- 1 - areas / sf::st_area(bb)

  # Construct the data frame with the results
  geom_indices <- data.frame(
    alpha = triangles$alpha[idx_triangles],
    geom_survival = geom_survival
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

# --- Estimate the envelope using Monte Carlo simulation ---
estimate_envelope <- function(
  spatgeom_obj,
  x,
  y,
  scale_pts,
  nalphas,
  domain_type,
  mc_cores = 2
) {
  for (i in seq_len(ncol(x))) {
    message(paste0("Estimating envelope for variable ", i))
    envelope_data <- parallel::mclapply(
      X = seq_len(40),
      FUN = function(k) {
        n <- stats::rpois(1, lambda = spatgeom_obj[[i]]$mean_n)
        x_rand <- stats::runif(n, min = min(x[, i]), max = max(x[, i]))
        y_rand <- stats::runif(n, min = min(y[, 1]), max = max(y[, 1]))
        enve <- estimate_curves(
          x1 = x_rand,
          x2 = y_rand,
          scale_pts = scale_pts,
          nalphas = nalphas,
          intensity = spatgeom_obj[[i]]$intensity,
          domain_type = domain_type
        )
        enve_approx <- stats::approx(
          x = enve$geom_indices$alpha,
          y = enve$geom_indices$geom_survival,
          xout = spatgeom_obj[[i]]$geom_indices$alpha
        )
        data.frame(enve_approx, nsim = k)
      },
      mc.cores = mc_cores
    )
    envelope_data <- do.call("rbind", envelope_data)
    spatgeom_obj[[i]]$envelope_data <- envelope_data
  }
  return(spatgeom_obj)
}
