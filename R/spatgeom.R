#' @title Geometric Spatial Point Pattern Analysis
#'
#' @description Function to estimate the geometric correlation between
#'   variables.
#'
#' @param x numeric matrix or data.frame. Either a matrix of covariables
#'   (paired with \code{y}), or a single \eqn{n \times p} matrix when \code{y}
#'   is omitted. In the latter case \code{x} must have at least 2 columns: if
#'   exactly 2, they are used directly as 2D coordinates; if more than 2, set
#'   \code{reduce = "pca"}, \code{"umap"}, or \code{"tsne"} to project to 2D
#'   first.
#' @param y numeric vector of responses. Optional: when omitted, \code{x} is
#'   treated as a single point cloud (see \code{reduce}).
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
#' @param reduce character; dimensionality reduction method to apply when
#'   \code{x} has more than 2 columns. One of \code{"none"} (default, keeps
#'   all columns and pairs each with \code{y} as separate variables),
#'   \code{"pca"} (Principal Component Analysis, no extra package required),
#'   \code{"umap"} (requires \pkg{uwot}), or \code{"tsne"} (requires
#'   \pkg{Rtsne}). When \code{reduce != "none"}, \code{x} is projected to 2
#'   dimensions via \code{\link{reduce_dim}} before the analysis.
#' @param reduce_args a named list of additional arguments forwarded to
#'   \code{\link{reduce_dim}} (and in turn to the underlying method). For
#'   example, \code{list(n_neighbors = 30)} when \code{reduce = "umap"}.
#'   Default \code{list()}.
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
#'
#' \item{\strong{reduce_embedding}}{Only present when \code{reduce != "none"}
#' and \code{x} had more than 2 columns. The \eqn{n \times 2} numeric matrix
#' produced by \code{\link{reduce_dim}}, whose columns were used as the 2D
#' point cloud.}
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
#'
#' # Wide matrix: reduce to 2D with PCA before running the analysis:
#' est_pca <- spatgeom(y = xy[, 1], x = xy[, -1], reduce = "pca")
#' plot_curve(est_pca, type = "curve")
#'
#' # Wide matrix: reduce with UMAP (requires the 'uwot' package):
#' if (requireNamespace("uwot", quietly = TRUE)) {
#'   est_umap <- spatgeom(y = xy[, 1], x = xy[, -1], reduce = "umap")
#'   plot_curve(est_umap, type = "curve")
#' }
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
  method = c("MAD", "DCLF"),
  reduce = c("none", "pca", "umap", "tsne"),
  reduce_args = list()
) {
  # --- Dimensionality reduction (applies regardless of whether y is present) --
  x <- as.data.frame(x)
  reduce <- match.arg(reduce)
  reduce_embedding <- NULL

  if (ncol(x) > 2 && reduce != "none") {
    reserved_reduce_args <- c("x", "method", "n_components")
    conflicting_reduce_args <- intersect(
      names(reduce_args),
      reserved_reduce_args
    )
    if (length(conflicting_reduce_args) > 0) {
      stop(
        "'reduce_args' cannot override reserved arguments: ",
        paste(conflicting_reduce_args, collapse = ", "),
        "."
      )
    }

    message(
      "Reducing x from ",
      ncol(x),
      " to 2 columns via method = '",
      reduce,
      "'."
    )
    reduce_embedding <- do.call(
      reduce_dim,
      c(list(x = as.matrix(x), method = reduce, n_components = 2L), reduce_args)
    )

    if (NCOL(reduce_embedding) != 2) {
      stop(
        "'reduce_dim' must return exactly 2 columns, but returned ",
        NCOL(reduce_embedding),
        "."
      )
    }

    x <- as.data.frame(reduce_embedding)
  }

  domain_type <- domain_type[1]

  if (hypothesis_testing && !envelope) {
    envelope <- TRUE
    message(
      "Hypothesis testing requires envelope simulation.",
      " Setting envelope = TRUE."
    )
  }

  if (missing(y)) {
    if (ncol(x) != 2) {
      stop(
        "'x' has ",
        ncol(x),
        " columns and 'y' is missing. ",
        "Set 'reduce' to \"pca\", \"umap\", or \"tsne\" to project 'x' to ",
        "2 columns, or supply a 'y' argument."
      )
    }

    x_2d <- x[, 1, drop = FALSE]
    y_2d <- x[, 2, drop = FALSE]

    spatgeom_obj <- spatgeom_xy(
      x = x_2d,
      y = y_2d,
      scale_pts = scale_pts,
      nalphas = nalphas,
      envelope = envelope,
      domain_type = domain_type,
      mc_cores = mc_cores
    )
  } else {
    spatgeom_obj <- spatgeom_xy(
      x,
      y,
      scale_pts = scale_pts,
      nalphas = nalphas,
      envelope = envelope,
      domain_type = domain_type,
      mc_cores = mc_cores
    )
  }

  if (!is.null(reduce_embedding)) {
    spatgeom_obj$reduce_embedding <- reduce_embedding
  }

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
