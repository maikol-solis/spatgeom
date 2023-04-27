#' CSR Hypothesis Testing for a \code{spatgeom} object
#'
#' Performs a global envelope test for Complete Spatial Randomness (CSR) on a
#' \code{spatgeom} object that was computed with \code{envelope = TRUE}. Two
#' test statistics are supported: Maximum Absolute Deviation (MAD) and
#' Diggle-Cressie-Loosmore-Ford (DCLF).
#'
#' @param spatgeom_obj an object of class \code{spatgeom} computed with
#'   \code{envelope = TRUE}.
#' @param significance_level a numeric value for the significance level of the
#'   test. Default \code{0.05}. Currently used for documentation purposes; the
#'   p-value is computed from the Monte Carlo rank.
#' @param r a numeric scaling parameter used in the theoretical CSR curve
#'   \eqn{\exp(-\lambda \pi (\alpha r)^2)}. Should match the value used when
#'   calling \code{\link{spatgeom}}. Default \code{0.5}.
#' @param method a character string, one of \code{"MAD"} (Maximum Absolute
#'   Deviation) or \code{"DCLF"} (Diggle-Cressie-Loosmore-Ford), specifying
#'   the global envelope test statistic. Default \code{"MAD"}.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\strong{hypothesis_testing_df}}{A data frame (one row per alpha
#'     value per variable) with columns \code{x} (alpha grid), \code{mean}
#'     (mean of simulated CSR curves), \code{alpha}, \code{geom_survival}
#'     (observed), \code{theor} (theoretical CSR curve),
#'     \code{upper_mean}, \code{lower_mean} (confidence band around the mean
#'     curve), \code{upper_theor}, \code{lower_theor} (confidence band around
#'     the theoretical curve), and \code{variable}.}
#'   \item{\strong{details}}{A named list, one entry per variable, each a
#'     data frame summarising the test: null hypothesis, variable name, test
#'     type, number of Monte Carlo simulations, observed and maximum test
#'     statistics, and p-values against both the mean and theoretical
#'     reference curves.}
#' }
#'
#' @examples
#' \donttest{
#' xy <- donut_data(n = 30, a = -1, b = 1, theta = 2 * pi)
#'
#' # Compute a spatgeom object with the Monte Carlo envelope:
#' est <- spatgeom(y = xy[, 1], x = xy[, -1], envelope = TRUE)
#'
#' # Test with the MAD (Maximum Absolute Deviation) statistic:
#' result_mad <- csr_test(est, method = "MAD")
#' result_mad$details
#'
#' # Test with the DCLF statistic:
#' result_dclf <- csr_test(est, method = "DCLF")
#' result_dclf$details
#' }
#'
#' @export

csr_test <- function(
  spatgeom_obj,
  significance_level = 0.05,
  r = 0.5,
  method = c("MAD", "DCLF")
) {
  method <- match.arg(method)

  if (!inherits(spatgeom_obj, "spatgeom")) {
    stop("spatgeom_obj must be an object of class 'spatgeom'.")
  }

  if (is.null(spatgeom_obj$results[[1]]$envelope_data)) {
    stop(
      "No envelope data found. ",
      "Re-run spatgeom() with envelope = TRUE (or hypothesis_testing = TRUE)."
    )
  }

  # For safety, declare local variables
  nsim <- y <- alpha <- geom_survival <- theor <- sum_sq_dif_mean <-
    sum_sq_dif_theor <- abs_dif_mean <- abs_dif_theor <- NULL

  # We'll use the input x only to get variable names.
  x <- spatgeom_obj$x

  # List to store per-variable results
  results_list <- list()

  hypothesis_testing_df <- NULL
  for (i in seq_len(ncol(x))) {
    # -------------------------------
    # 1. Get the envelope data and observed values
    # -------------------------------
    f_simul <- spatgeom_obj$results[[i]]$envelope_data
    # f_obs: observed geometric survival values (from the alpha–shape curve)
    f_obs <- spatgeom_obj$results[[i]]$geom_indices[, 2]

    # To ensure matching keys, round the x values.
    f_obs_df <- data.frame(
      x = round(spatgeom_obj$results[[i]]$geom_indices[, 1], 8),
      y = f_obs
    )

    f_obs_with_num <- f_obs_df |>
      dplyr::group_by(x) |>
      dplyr::summarise(y = mean(y, na.rm = TRUE)) |>
      dplyr::ungroup() |>
      dplyr::mutate(nsim = 0)

    # Also round f_simul$x so that the later join succeeds.
    f_simul <- f_simul |> dplyr::mutate(x = round(x, 8))

    f_simul_and_obs <- rbind(f_simul, f_obs_with_num)

    # -------------------------------
    # 2. Compute the grouped (observed) mean curve
    # -------------------------------
    f_mean <- f_simul_and_obs |>
      dplyr::group_by(x) |>
      dplyr::summarise(mean = mean(y, na.rm = TRUE)) |>
      dplyr::ungroup() |>
      dplyr::mutate(x = round(x, 8))

    # -------------------------------
    # 3. Join with the original observed survival data and compute the
    # theoretical curve
    # -------------------------------

    # Here, spatgeom_obj$results[[i]]$geom_indices has columns "alpha" and
    # "geom_survival". We use "alpha" as our x.
    f_dataframe <- dplyr::full_join(
      f_mean,
      spatgeom_obj$results[[i]]$geom_indices |>
        dplyr::mutate(x = round(alpha, 8)),
      by = "x"
    )

    intensity <- spatgeom_obj$results[[i]]$intensity
    f_dataframe <- f_dataframe |>
      dplyr::mutate(theor = exp(-intensity * pi * (x * r)^2))

    # -------------------------------
    # 4. Compute test statistics (either DCLF or MAD)
    # -------------------------------
    if (method == "DCLF") {
      # Use a squared-difference measure. Use the observed "geom_survival"
      # column.
      delta_r <- f_mean$x[2] - f_mean$x[1]
      f_dataframe <- f_dataframe |>
        dplyr::mutate(
          sum_sq_diff_mean = (geom_survival - mean)^2 * delta_r,
          sum_sq_diff_theor = (geom_survival - theor)^2 * delta_r
        )
      t_obs_mean <- sum(f_dataframe$sum_sq_diff_mean, na.rm = TRUE)
      t_obs_theor <- sum(f_dataframe$sum_sq_diff_theor, na.rm = TRUE)

      # Instead of using rep(), join the grouped data to f_simul by "x".
      f_simul <- dplyr::left_join(f_simul, f_mean, by = "x")
      f_simul <- dplyr::left_join(
        f_simul,
        f_dataframe |> dplyr::select(x, theor),
        by = "x"
      )

      f_simul <- f_simul |>
        dplyr::mutate(
          sum_sq_dif_mean = (mean - y)^2 * delta_r,
          sum_sq_dif_theor = (theor - y)^2 * delta_r
        )

      t_simul_mean <- f_simul |>
        dplyr::group_by(nsim) |>
        dplyr::summarise(maxim = sum(sum_sq_dif_mean, na.rm = TRUE))
      t_simul_theor <- f_simul |>
        dplyr::group_by(nsim) |>
        dplyr::summarise(maxim = sum(sum_sq_dif_theor, na.rm = TRUE))

      tmax_mean <- max(t_simul_mean$maxim, na.rm = TRUE)
      tmax_theor <- max(t_simul_theor$maxim, na.rm = TRUE)
    } else {
      # MAD method: compare the absolute differences
      f_dataframe <- f_dataframe |>
        dplyr::mutate(
          abs_dif_mean = abs(geom_survival - mean),
          abs_dif_theor = abs(geom_survival - theor)
        )
      t_obs_mean <- max(f_dataframe$abs_dif_mean, na.rm = TRUE)
      t_obs_theor <- max(f_dataframe$abs_dif_theor, na.rm = TRUE)

      f_simul <- dplyr::left_join(f_simul, f_mean, by = "x")
      f_simul <- dplyr::left_join(
        f_simul,
        f_dataframe |> dplyr::select(x, theor),
        by = "x"
      )

      f_simul <- f_simul |>
        dplyr::mutate(
          abs_dif_mean = abs(mean - y),
          abs_dif_theor = abs(theor - y)
        )

      t_simul_mean <- f_simul |>
        dplyr::group_by(nsim) |>
        dplyr::summarise(maxim = max(abs_dif_mean, na.rm = TRUE))
      t_simul_theor <- f_simul |>
        dplyr::group_by(nsim) |>
        dplyr::summarise(maxim = max(abs_dif_theor, na.rm = TRUE))

      tmax_mean <- max(t_simul_mean$maxim, na.rm = TRUE)
      tmax_theor <- max(t_simul_theor$maxim, na.rm = TRUE)
    }

    # -------------------------------
    # 5. Compute p-values using the ordered test statistics from the simulations
    # -------------------------------
    t_simul_ord_mean <- sort(t_simul_mean$maxim, decreasing = TRUE)
    m_mean <- length(t_simul_ord_mean)
    k_mean <- if (t_obs_mean < min(t_simul_ord_mean)) {
      m_mean + 1
    } else {
      which(t_obs_mean > t_simul_ord_mean)[1]
    }

    t_simul_ord_theor <- sort(t_simul_theor$maxim, decreasing = TRUE)
    m_theor <- length(t_simul_ord_theor)
    k_theor <- if (t_obs_theor < min(t_simul_ord_theor)) {
      m_theor + 1
    } else {
      which(t_obs_theor > t_simul_ord_theor)[1]
    }

    results <- data.frame(
      Feature = c(
        "H_0",
        "Variable",
        "Type of Test",
        "Number of Monte Carlo Simulations",
        "t_max_mean",
        "t_obs_mean",
        "t_obs_mean > t_max_mean",
        "p-value_mean",
        "t_max_theo",
        "t_obs_theor",
        "t_obs_theor > t_max_theo",
        "p-value_theor"
      ),
      Value = c(
        "CSR process",
        colnames(x)[i],
        method,
        m_mean,
        round(tmax_mean, 4),
        round(t_obs_mean, 4),
        t_obs_mean > tmax_mean,
        k_mean / (m_mean + 1),
        round(tmax_theor, 4),
        round(t_obs_theor, 4),
        t_obs_theor > tmax_theor,
        k_theor / (m_theor + 1)
      )
    )

    f_dataframe <- f_dataframe |>
      dplyr::mutate(
        upper_mean = mean + tmax_mean,
        lower_mean = mean - tmax_mean,
        upper_theor = theor + tmax_theor,
        lower_theor = theor - tmax_theor,
        variable = colnames(x)[i]
      )
    hypothesis_testing_df <- rbind(hypothesis_testing_df, f_dataframe)

    results_list[[colnames(x)[i]]] <- results
  }

  return(list(
    hypothesis_testing_df = hypothesis_testing_df,
    details = results_list
  ))
}
