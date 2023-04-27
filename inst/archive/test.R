#' # --- Load required packages ---
#' library(sf)
#' library(lwgeom)
#' library(parallel)
#' library(dplyr)
#' library(tidyr)
#' library(scales)
#' library(stats)
#'
#' # --- Unified spatgeom function with integrated hypothesis testing ---
#' #' @title Geometric Spatial Point Pattern Analysis with CSR Hypothesis Testing
#' #'
#' #' @description Estimates the geometric survival (empty space) of a point process
#' #'   and (optionally) performs CSR hypothesis testing via a global envelope approach.
#' #'
#' #' @param x A numeric matrix or data.frame of covariables.
#' #' @param y A numeric vector or one‐column data.frame of responses.
#' #' @param scale_pts Logical. If TRUE the points are rescaled.
#' #' @param nalphas Integer. Number of alpha values (steps) for computing the geometric measure.
#' #' @param envelope Logical. If TRUE the envelope (Monte Carlo simulation) is computed.
#' #'   (Note: If `hypothesis_testing=TRUE` then `envelope` is forced to TRUE.)
#' #' @param domain_type Character. Either `"bounding-box"` or `"convex-hull"` (default is `"bounding-box"`).
#' #' @param hypothesis_testing Logical. If TRUE, the function performs a CSR hypothesis test.
#' #' @param significance_level Numeric. The significance level for the hypothesis test.
#' #' @param mc_cores Integer. The number of cores to use in parallel computations.
#' #' @param r Numeric. A parameter used in the theoretical function in the hypothesis test.
#' #' @param method Character. One of `"MAD"` or `"DCLF"`. Specifies the test statistic to use.
#' #' @param name_ex Character. A name string used to label output.
#' #'
#' #' @return An object of class `spatgeom` (a list) with elements including the geometric survival curves
#' #'         and (if requested) the hypothesis testing results.
#' #' @export
#' spatgeom <- function(
#'   x,
#'   y,
#'   scale_pts = FALSE,
#'   nalphas = 100,
#'   envelope = FALSE,
#'   domain_type = c("bounding-box", "convex-hull"),
#'   hypothesis_testing = FALSE,
#'   significance_level = 0.05,
#'   mc_cores = 1,
#'   r = 0.5,
#'   method = c("MAD", "DCLF"),
#'   name_ex = "ex1"
#' ) {
#'   if (missing(y)) {
#'     stop("The argument y is required when running spatgeom.")
#'   } else {
#'     message("Running spatgeom with both x and y.")
#'     domain_type <- domain_type[1]
#'     # If hypothesis testing is requested, force envelope simulation.
#'     if (hypothesis_testing && !envelope) {
#'       envelope <- TRUE
#'       message("Hypothesis testing requires envelope simulation. Setting envelope=TRUE.")
#'     }
#'
#'     # Compute the basic spatgeom object (without hypothesis testing)
#'     spatgeom_obj <- spatgeom_xy(
#'       x,
#'       y,
#'       scale_pts = scale_pts,
#'       nalphas = nalphas,
#'       envelope = envelope,
#'       domain_type = domain_type,
#'       significance_level = significance_level,
#'       mc_cores = mc_cores
#'     )
#'
#'     # If requested, run the CSR hypothesis test and add its output to the object.
#'     if (hypothesis_testing) {
#'       spatgeom_obj$hypothesis_testing_results <- csr_test(
#'         spatgeom_obj,
#'         significance_level = significance_level,
#'         r = r,
#'         method = method,
#'         name_ex = name_ex
#'       )
#'     }
#'
#'     return(spatgeom_obj)
#'   }
#' }
#'
#' # --- Internal function that computes the geometric survival curves ---
#' spatgeom_xy <- function(
#'   x,
#'   y,
#'   scale_pts = FALSE,
#'   nalphas = 100,
#'   envelope = FALSE,
#'   domain_type = c("bounding-box", "convex-hull"),
#'   significance_level = 0.05,
#'   mc_cores = 2
#' ) {
#'   x <- as.data.frame(x)
#'   y <- as.data.frame(y)
#'   domain_type <- domain_type[1]
#'
#'   ans <- list()
#'   ans[["call"]] <- match.call()
#'   ans[["x"]] <- x
#'   ans[["y"]] <- y
#'
#'   message("Estimating geometric survival of empty space...")
#'   out_list <- parallel::mclapply(
#'     X = seq_len(ncol(x)),
#'     FUN = function(i) {
#'       message(paste0("Processing variable ", i, " (", colnames(x)[i], ")"))
#'       estimate_curves(
#'         x1 = x[, i],
#'         x2 = y[, 1],
#'         scale_pts = scale_pts,
#'         nalphas = nalphas,
#'         domain_type = domain_type
#'       )
#'     },
#'     mc.cores = mc_cores
#'   )
#'
#'   # Attach the variable names.
#'   out_list <- lapply(
#'     X = seq_len(ncol(x)),
#'     FUN = function(i) {
#'       append(out_list[[i]], list(variable_name = colnames(x)[i]))
#'     }
#'   )
#'
#'   if (envelope == TRUE) {
#'     out_list <- estimate_envelope(
#'       spatgeom_obj = out_list,
#'       x = x,
#'       y = y,
#'       scale_pts = scale_pts,
#'       nalphas = nalphas,
#'       domain_type = domain_type,
#'       mc_cores = mc_cores
#'     )
#'   }
#'
#'   ans[["results"]] <- out_list
#'   class(ans) <- "spatgeom"
#'   return(ans)
#' }
#'
#' # --- Estimate the geometric survival curves using alpha-shapes ---
#' estimate_curves <- function(
#'   x1,
#'   x2,
#'   scale_pts,
#'   nalphas,
#'   intensity = NULL,
#'   domain_type
#' ) {
#'   # Scale points if requested.
#'   coords <- if (scale_pts) {
#'     scales::rescale(cbind(x1, x2))
#'   } else {
#'     cbind(x1, x2)
#'   }
#'   pts <- sf::st_multipoint(coords)
#'   pts <- sf::st_sfc(pts)
#'   pts <- sf::st_cast(pts, "POINT")
#'
#'   # Create the domain.
#'   if (domain_type == "convex-hull") {
#'     bb <- sf::st_convex_hull(sf::st_union(pts))
#'   } else if (domain_type == "bounding-box") {
#'     bb <- sf::st_make_grid(pts, n = 1)
#'   } else {
#'     stop("domain_type must be either 'convex-hull' or 'bounding-box'")
#'   }
#'
#'   # Estimate the intensity.
#'   if (is.null(intensity)) {
#'     intensity <- length(pts) / sf::st_area(bb)
#'   }
#'
#'   # Perform Delaunay triangulation.
#'   delaunay <- sf::st_triangulate(sf::st_geometrycollection(pts))
#'   polygons <- sf::st_collection_extract(sf::st_sfc(delaunay))
#'   linestrings <- sf::st_cast(polygons, "LINESTRING")
#'   linestrings_splitted <- lwgeom::st_split(linestrings, pts)
#'
#'   # Compute the maximum segment length for each linestring.
#'   max_length <- sapply(
#'     linestrings_splitted,
#'     function(x) {
#'       max(sf::st_length(sf::st_cast(sf::st_sfc(x))))
#'     }
#'   )
#'
#'   triangles <- data.frame(
#'     geometry = polygons,
#'     max_length = max_length,
#'     alpha = max_length / 2
#'   )
#'   triangles <- triangles[order(triangles$alpha), ]
#'
#'   number_triangles <- nrow(triangles)
#'   step_seq <- ifelse(nalphas < number_triangles, ceiling(number_triangles / nalphas), 1)
#'   idx_triangles <- unique(c(seq(1, number_triangles, by = step_seq), number_triangles))
#'
#'   individual_geometries <- lapply(seq_along(idx_triangles), function(k) {
#'     if (k == 1) {
#'       triangles[1, ]$geometry
#'     } else {
#'       sf::st_union(triangles[idx_triangles[k - 1]:(idx_triangles[k] - 1), ]$geometry)
#'     }
#'   })
#'
#'   cumulative_union <- Reduce(sf::st_union, individual_geometries, accumulate = TRUE)
#'
#'   areas <- sapply(cumulative_union, sf::st_area)
#'   geom_survival <- 1 - areas / sf::st_area(bb)
#'
#'   geom_indices <- data.frame(
#'     alpha = triangles$alpha[idx_triangles],
#'     geom_survival = geom_survival
#'   )
#'
#'   return(list(
#'     triangles = triangles,
#'     geom_indices = geom_indices,
#'     intensity = intensity,
#'     mean_n = sf::st_area(bb) * intensity
#'   ))
#' }
#'
#' # --- Estimate the envelope using Monte Carlo simulation ---
#' estimate_envelope <- function(
#'   spatgeom_obj,
#'   x,
#'   y,
#'   scale_pts,
#'   nalphas,
#'   domain_type,
#'   mc_cores = 2
#' ) {
#'   for (i in seq_len(ncol(x))) {
#'     message(paste0("Estimating envelope for variable ", i))
#'     envelope_data <- parallel::mclapply(
#'       X = seq_len(40),
#'       FUN = function(k) {
#'         n <- stats::rpois(1, lambda = spatgeom_obj[[i]]$mean_n)
#'         x_rand <- stats::runif(n, min = min(x[, i]), max = max(x[, i]))
#'         y_rand <- stats::runif(n, min = min(y[, 1]), max = max(y[, 1]))
#'         enve <- estimate_curves(
#'           x1 = x_rand,
#'           x2 = y_rand,
#'           scale_pts = scale_pts,
#'           nalphas = nalphas,
#'           intensity = spatgeom_obj[[i]]$intensity,
#'           domain_type = domain_type
#'         )
#'         enve_approx <- stats::approx(
#'           x = enve$geom_indices$alpha,
#'           y = enve$geom_indices$geom_survival,
#'           xout = spatgeom_obj[[i]]$geom_indices$alpha
#'         )
#'         data.frame(enve_approx, nsim = k)
#'       },
#'       mc.cores = mc_cores
#'     )
#'     envelope_data <- do.call("rbind", envelope_data)
#'     spatgeom_obj[[i]]$envelope_data <- envelope_data
#'   }
#'   return(spatgeom_obj)
#' }
#'
#' # --- CSR Hypothesis Testing Function (Global Envelope Test) ---
#' csr_test <- function(
#'   spatgeom_obj,
#'   significance_level = 0.05,
#'   r = 0.5,
#'   method = c("MAD", "DCLF"),
#'   name_ex = "ex1"
#' ) {
#'   # For safety, declare local variables
#'   nsim <- y <- sum_sq_dif <- sum_sq_dif2 <- abs_dif <- abs_dif2 <- x <- NULL
#'
#'   if (length(method) > 1 || !(method %in% c("MAD", "DCLF"))) {
#'     stop("Parameter method must be one of 'MAD' or 'DCLF'.")
#'   }
#'
#'   # We'll use the input x only to get variable names.
#'   x <- spatgeom_obj$x
#'
#'   # List to store per-variable results
#'   results_list <- list()
#'
#'   for (i in seq_len(ncol(x))) {
#'     # -------------------------------
#'     # 1. Get the envelope data and observed values
#'     # -------------------------------
#'     f_simul <- spatgeom_obj$results[[i]]$envelope_data
#'     # f_obs: observed geometric survival values (from the alpha–shape curve)
#'     f_obs <- spatgeom_obj$results[[i]]$geom_indices[, 2]
#'
#'     # To ensure matching keys, round the x values.
#'     f_obs_df <- data.frame(
#'       x = round(spatgeom_obj$results[[i]]$geom_indices[, 1], 8),
#'       y = f_obs
#'     )
#'
#'     f_obs_with_num <- f_obs_df %>%
#'       dplyr::group_by(x) %>%
#'       dplyr::summarise(y = mean(y, na.rm = TRUE)) %>%
#'       dplyr::ungroup() %>%
#'       dplyr::mutate(nsim = 0)
#'
#'     # Also round f_simul$x so that the later join succeeds.
#'     f_simul <- f_simul %>% dplyr::mutate(x = round(x, 8))
#'
#'     f_simul_and_obs <- rbind(f_simul, f_obs_with_num)
#'
#'     # -------------------------------
#'     # 2. Compute the grouped (observed) mean curve
#'     # -------------------------------
#'     f_mean <- f_simul_and_obs %>%
#'       dplyr::group_by(x) %>%
#'       dplyr::summarise(mean = mean(y, na.rm = TRUE)) %>%
#'       dplyr::ungroup() %>%
#'       dplyr::mutate(x = round(x, 8))
#'
#'     # -------------------------------
#'     # 3. Join with the original observed survival data and compute the theoretical curve
#'     # -------------------------------
#'     # Here, spatgeom_obj$results[[i]]$geom_indices has columns "alpha" and "geom_survival".
#'     # We use "alpha" as our x.
#'     f_dataframe <- dplyr::full_join(
#'       f_mean,
#'       spatgeom_obj$results[[i]]$geom_indices %>% dplyr::mutate(x = round(alpha, 8)),
#'       by = "x"
#'     )
#'
#'     intensity <- spatgeom_obj$results[[i]]$intensity
#'     f_dataframe <- f_dataframe %>%
#'       dplyr::mutate(theor = exp(-intensity * pi * (x * r)^2))
#'
#'     # -------------------------------
#'     # 4. Compute test statistics (either DCLF or MAD)
#'     # -------------------------------
#'     if (method == "DCLF") {
#'       # Use a squared-difference measure. Use the observed "geom_survival" column.
#'       delta_r <- f_mean$x[2] - f_mean$x[1]
#'       f_dataframe <- f_dataframe %>%
#'         dplyr::mutate(sum_sq_diff_mean = (geom_survival - mean)^2 * delta_r,
#'                       sum_sq_diff_theor = (geom_survival - theor)^2 * delta_r)
#'       t_obs_mean <- sum(f_dataframe$sum_sq_diff_mean, na.rm = TRUE)
#'       t_obs_theor <- sum(f_dataframe$sum_sq_diff_theor, na.rm = TRUE)
#'
#'       # Instead of using rep(), join the grouped data to f_simul by "x".
#'       f_simul <- dplyr::left_join(f_simul, f_mean, by = "x")
#'       f_simul <- dplyr::left_join(f_simul, f_dataframe %>% dplyr::select(x, theor), by = "x")
#'
#'       f_simul <- f_simul %>%
#'         dplyr::mutate(sum_sq_dif_mean = (mean - y)^2 * delta_r,
#'                       sum_sq_dif_theor = (theor - y)^2 * delta_r)
#'
#'       t_simul_mean <- f_simul %>%
#'         dplyr::group_by(nsim) %>%
#'         dplyr::summarise(maxim = sum(sum_sq_dif_mean, na.rm = TRUE))
#'       t_simul_theor <- f_simul %>%
#'         dplyr::group_by(nsim) %>%
#'         dplyr::summarise(maxim = sum(sum_sq_dif_theor, na.rm = TRUE))
#'
#'       tmax_mean <- max(t_simul_mean$maxim, na.rm = TRUE)
#'       tmax_theor <- max(t_simul_theor$maxim, na.rm = TRUE)
#'     } else {
#'       # MAD method: compare the absolute differences
#'       f_dataframe <- f_dataframe %>%
#'         dplyr::mutate(abs_dif_mean = abs(geom_survival - mean),
#'                       abs_dif_theor = abs(geom_survival - theor))
#'       t_obs_mean <- max(f_dataframe$abs_dif_mean, na.rm = TRUE)
#'       t_obs_theor <- max(f_dataframe$abs_dif_theor, na.rm = TRUE)
#'
#'       f_simul <- dplyr::left_join(f_simul, f_mean, by = "x")
#'       f_simul <- dplyr::left_join(f_simul, f_dataframe %>% dplyr::select(x, theor), by = "x")
#'
#'       f_simul <- f_simul %>%
#'         dplyr::mutate(abs_dif = abs(mean - y),
#'                       abs_dif2 = abs(theor - y))
#'
#'       t_simul_mean <- f_simul %>%
#'         dplyr::group_by(nsim) %>%
#'         dplyr::summarise(maxim = max(abs_dif_mean, na.rm = TRUE))
#'       t_simul_theor <- f_simul %>%
#'         dplyr::group_by(nsim) %>%
#'         dplyr::summarise(maxim = max(abs_dif_theor, na.rm = TRUE))
#'
#'       tmax_mean <- max(t_simul_mean$maxim, na.rm = TRUE)
#'       tmax_theor <- max(t_simul_theor$maxim, na.rm = TRUE)
#'     }
#'
#'     # -------------------------------
#'     # 5. Compute p-values using the ordered test statistics from the simulations
#'     # -------------------------------
#'     t_simul_ord <- sort(t_simul_mean$maxim, decreasing = TRUE)
#'     m <- length(t_simul_ord)
#'     k <- if (t_obs_mean < min(t_simul_ord)) m + 1 else which(t_obs_mean > t_simul_ord)[1]
#'
#'     t_simul_ord2 <- sort(t_simul_theor$maxim, decreasing = TRUE)
#'     m2 <- length(t_simul_ord2)
#'     k2 <- if (t_obs_theor < min(t_simul_ord2)) m2 + 1 else which(t_obs_theor > t_simul_ord2)[1]
#'
#'     results <- data.frame(
#'       Feature = c("H_0", "Variable", "Type of Test", "Number of Monte Carlo Simulations",
#'                   "T_max", "t_obs_mean", "t_obs_mean > T_max", "p-value",
#'                   "T_max2", "t_obs_theor", "t_obs_theor > T_max2", "p-value2"),
#'       Value = c("CSR process",
#'                 colnames(x)[i],
#'                 method,
#'                 m,
#'                 round(tmax_mean, 4),
#'                 round(t_obs_mean, 4),
#'                 t_obs_mean > tmax_mean,
#'                 k / (m + 1),
#'                 round(tmax_theor, 4),
#'                 round(t_obs_theor, 4),
#'                 t_obs_theor > tmax_theor,
#'                 k2 / (m2 + 1))
#'     )
#'
#'     print(results)
#'
#'     f_dataframe <- f_dataframe %>%
#'       dplyr::mutate(upper = mean + tmax_mean,
#'                     lower = mean - tmax_mean,
#'                     upper2 = theor + tmax_theor,
#'                     lower2 = theor - tmax_theor,
#'                     variable = colnames(x)[i])
#'
#'     if (i == 1) {
#'       f_df_mad <- f_dataframe
#'     } else {
#'       f_df_mad <- rbind(f_df_mad, f_dataframe)
#'     }
#'
#'     results_list[[colnames(x)[i]]] <- results
#'   }
#'
#'   return(list(f_df1 = f_df_mad, details = results_list))
#' }

# --- End of unified spatgeom and hypothesis testing definitions ---
#
#
# # ===============================
# # Example 1: Ring (donut) data example
# # ===============================
# set.seed(123)
# n <- 100
# a <- -1
# b <- 1
# theta <- runif(n, 0, 2 * pi)
# r_val <- (sqrt(runif(n)))*(0.5) + 0.5
# X1 <- r_val * cos(theta)
# X2 <- runif(n, a, b)
# Y_ring <- r_val * sin(theta)
# X_ring <- data.frame(X1, X2)
# Y_ring_df <- data.frame(Y = Y_ring)
#
# message("\nRunning Example 1: Ring data")
# result1 <- spatgeom(
#   x = X_ring,
#   y = Y_ring_df,
#   hypothesis_testing = TRUE,
#   significance_level = 0.05,
#   r = 0.5,
#   method = "MAD",
#   name_ex = "ring_100pts"
# )
# print(result1$hypothesis_testing_results)
#
#
# # ===============================
# # Example 2: Linear relationship data
# # ===============================
# set.seed(456)
# n <- 100
# a <- -2
# b <- 2
# X1_lin <- runif(n, a, b)
# X2_lin <- runif(n, a, b)
# X3_lin <- runif(n, a, b)
# Y_lin <- 0.6 * X1_lin + 0.3 * X2_lin + 0.1 * X3_lin
# X_lin <- data.frame(X1 = X1_lin, X2 = X2_lin, X3 = X3_lin)
# Y_lin_df <- data.frame(Y = Y_lin)
#
# message("\nRunning Example 2: Linear data")
# result2 <- spatgeom(
#   x = X_lin,
#   y = Y_lin_df,
#   hypothesis_testing = TRUE,
#   significance_level = 0.05,
#   r = 0.5,
#   method = "MAD",
#   name_ex = "linear_100pts"
# )
# print(result2$hypothesis_testing_results)
#
#
# # ===============================
# # Example 3: Small dataset with 7 points
# # ===============================
# X_small <- data.frame(X = c(-3, 4, -1, 2, -3, 0, 1))
# Y_small <- data.frame(Y = c(1, 2, 4, 5, 6, 1, 7))
#
# message("\nRunning Example 3: Small 7-point data")
# result3 <- spatgeom(
#   x = X_small,
#   y = Y_small,
#   hypothesis_testing = TRUE,
#   significance_level = 0.05,
#   r = 0.5,
#   method = "MAD",
#   name_ex = "7pts"
# )
# print(result3$hypothesis_testing_results)
