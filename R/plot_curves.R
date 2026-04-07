#' plot spatgeom objects with hypothesis testing results and confidence
#' intervals
#'
#' S3 generic for plotting spatial geometry analysis results. Dispatches to
#' \code{\link{plot_curve.spatgeom}} for \code{spatgeom} objects (producing a
#' \code{ggplot} of survival curves or their derivatives), and to
#' \code{\link{plot_curve.spatgeom_group}} for \code{spatgeom_group} objects
#' (producing a side-by-side \code{cowplot} grid, one panel per group).
#'
#' @param x an object of class \code{spatgeom} or \code{spatgeom_group}.
#' @param type a string: either \code{"curve"} (the default) or \code{"deriv"}.
#'   The option \code{"curve"} plots the survival curve (with envelopes and,
#'   if available, confidence intervals). The option \code{"deriv"} plots its
#'   numerical derivative.
#' @param font_size an integer controlling the font size in the plot.
#' @param ... further arguments passed to the appropriate method.
#'
#' @return a \code{ggplot} object for \code{spatgeom} inputs, or a
#'   \code{cowplot} grid for \code{spatgeom_group} inputs.
#'
#' @examples
#'
#' xy <- donut_data(n = 30, a = -1, b = 1, theta = 2 * pi)
#'
#' # Basic plots — no envelope, no hypothesis testing:
#' estimation <- spatgeom(y = xy[, 1], x = xy[, -1])
#' plot_curve(estimation, type = "curve")
#' plot_curve(estimation, type = "deriv")
#'
#' \donttest{
#' # Curve with Monte Carlo envelope ribbon (grey) and theoretical CSR curve
#' # (red dashed):
#' est_env <- spatgeom(y = xy[, 1], x = xy[, -1], envelope = TRUE)
#' plot_curve(est_env, type = "curve")
#'
#' # Curve with hypothesis testing results: additionally shows the CSR mean
#' # (blue dotted) and its confidence band (blue ribbon):
#' est_ht <- spatgeom(
#'   y = xy[, 1], x = xy[, -1],
#'   hypothesis_testing = TRUE, method = "MAD"
#' )
#' plot_curve(est_ht, type = "curve")
#' }
#'
#' @export

plot_curve <- function(x, ...) UseMethod("plot_curve")

#' @rdname plot_curve
#' @export

plot_curve.spatgeom <- function(x, type = "curve", font_size = 12, ...) {
  # Pacify linting warnings
  alpha <- geom_survival <- variable <- y <-
    ymin <- ymax <- lower_mean <- upper_mean <- theor <- mean <- NULL

  # Number of variables (panels)
  nvar <- length(x$results)

  # Initialize empty data frames for the basic curves and their derivatives.
  df <- NULL
  df_fp <- NULL

  # Loop over each variable (each panel)
  for (k in 1:nvar) {
    # Extract the alpha–values and survival curve from the results.
    temp_df <- cbind(
      x$results[[k]]$geom_indices,
      variable = colnames(x$x)[k],
      intensity = x$results[[k]]$intensity
    )
    df <- rbind(df, temp_df)

    # Compute a numerical derivative of the survival curve.
    x_curve <- x$results[[k]]$geom_indices$alpha
    y_curve <- x$results[[k]]$geom_indices$geom_survival
    curve_regularized <- stats::approx(x_curve, y_curve, n = length(x_curve))
    dx <- diff(curve_regularized$x)
    dy <- diff(curve_regularized$y)
    temp_df_fp <- data.frame(
      x = curve_regularized$x[-length(curve_regularized$x)],
      y = dy / dx,
      variable = colnames(x$x)[k]
    )
    df_fp <- rbind(df_fp, temp_df_fp)
  }

  if (type == "curve") {
    plt <- ggplot2::ggplot(df, ggplot2::aes(x = alpha, y = geom_survival))

    # If envelope data are available, add a grey ribbon.
    # (Envelope data were computed for each variable and stored in each result.)
    env_df <- NULL
    for (k in 1:nvar) {
      env_temp <- x$results[[k]]$envelope_data
      if (!is.null(env_temp)) {
        env_temp$variable <- colnames(x$x)[k]
        env_df <- rbind(env_df, env_temp)
      }
    }
    if (!is.null(env_df)) {
      envelope_ribbon <- env_df |>
        dplyr::group_by(x, variable) |>
        dplyr::summarise(
          ymin = min(y, na.rm = TRUE),
          ymax = max(y, na.rm = TRUE)
        ) |>
        dplyr::ungroup() |>
        dplyr::filter(is.finite(ymin) & is.finite(ymax)) |>
        dplyr::left_join(df, by = c("x" = "alpha", "variable" = "variable"))
      plt <- plt +
        ggplot2::geom_ribbon(
          data = envelope_ribbon,
          ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
          alpha = 0.7,
          fill = "grey"
        )
    }

    # Plot the step–function (the observed survival curve)
    plt <- plt + ggplot2::geom_step(linewidth = 1)

    # Add the theoretical CSR curve. If hypothesis testing was performed, use
    # the pre-computed 'theor' column (which correctly incorporates the r
    # parameter). Otherwise fall back to computing it from intensity alone.
    if (!is.null(x$hypothesis_testing_results)) {
      ht <- x$hypothesis_testing_results$hypothesis_testing_df
      plt <- plt +
        ggplot2::geom_line(
          data = ht,
          ggplot2::aes(x = x, y = theor),
          linetype = "dashed",
          color = "red",
          linewidth = 1
        ) +
        ggplot2::geom_ribbon(
          data = ht,
          ggplot2::aes(x = x, ymin = lower_mean, ymax = upper_mean),
          fill = "blue",
          alpha = 0.3
        ) +
        ggplot2::geom_line(
          data = ht,
          ggplot2::aes(x = x, y = mean),
          color = "blue",
          linetype = "dotted",
          linewidth = 1
        )
    } else {
      for (k in 1:nvar) {
        intensity_k <- x$results[[k]]$intensity
        plt <- plt +
          ggplot2::geom_function(
            data = subset(df, variable == colnames(x$x)[k]),
            mapping = ggplot2::aes(x = alpha),
            fun = function(a, intensity) {
              exp(-intensity * pi * (a)^2)
            },
            args = list(intensity = intensity_k),
            linetype = "dashed",
            color = "red",
            linewidth = 1
          )
      }
    }

    plt <- plt +
      ggplot2::scale_y_continuous(name = expression(f(alpha))) +
      ggplot2::scale_x_continuous(
        name = expression(alpha),
        labels = scales::comma
      ) +
      ggplot2::facet_wrap(~variable, scales = "free_x") +
      cowplot::theme_cowplot(font_size = font_size) +
      cowplot::background_grid(minor = "y") +
      cowplot::panel_border()
  } else if (type == "deriv") {
    plt <- ggplot2::ggplot(data = df_fp, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(color = "grey", alpha = 0.7) +
      ggplot2::geom_smooth(
        linewidth = 1,
        se = FALSE,
        method = "gam",
        formula = y ~ s(x, bs = "cs")
      ) +
      ggplot2::scale_y_continuous(name = expression(f * minute(alpha))) +
      ggplot2::scale_x_continuous(
        name = expression(alpha),
        labels = scales::comma
      ) +
      ggplot2::facet_wrap(~variable, scales = "free_x") +
      cowplot::theme_cowplot(font_size = font_size) +
      cowplot::background_grid(minor = "y") +
      cowplot::panel_border()
  }

  return(plt)
}
