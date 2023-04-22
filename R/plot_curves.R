#' plot \code{spatgeom} objects
#'
#' Plot method for objects of class \code{spatgeom}.
#'
#' @param x an object of class \code{spatgeom}
#' @param type a string that could be \code{curve} or \code{deriv}. The option
#'   \code{curve} plots the curve of \code{alpha} against \code{geom_corr} from
#'   the function [`spatgeom::spatgeom()`]. The \code{deriv} option plots the
#'   numerical derivative.
#' @param font_size a integer that increases the font size in the plot.
#'
#' @return a \code{\link[ggplot2]{ggplot}} object with the geometric indices (or
#'   its derivative). The plot is generated with the \code{nalphas} point of
#'   \code{alpha} and \code{geom_corr} from the function
#'   \code{\link{spatgeom}}.
#'
#' In each panel, the theoretical CSR process is drawn using
#'   \code{exp(-intensity * pi * x^2)}. where the intensity depends on each
#'   panel.
#'
#' @examples
#'
#' xy <- donut_data(n = 30, a = -1, b = 1, theta = 2 * pi)
#'
#' estimation <- spatgeom(y = xy[, 1], x = xy[, -1])
#'
#' plot_curve(estimation, type = "curve")
#'
#' plot_curve(estimation, type = "deriv")
#'
#' @export

plot_curve <-
  function(x,
           type = "curve",
           font_size = 12) {
    # Number of variables
    nvar <- length(x$results)

    # variables bindings for R CMD check
    df <- df_fp <- alpha <- geom_corr <- variable <- y <- ymin <- ymax <- NULL

    for (k in 1:nvar) {
      df <- rbind(
        cbind(
          x$results[[k]]$geom_indices,
          variable = colnames(x$x)[k],
          intensity = x$results[[k]]$intensity
        ),
        df
      )

      x_curve <- x$results[[k]]$geom_indices$alpha
      y_curve <- x$results[[k]]$geom_indices$geom_corr

      suppressWarnings(
        curve_regularized <- stats::approx(
          x_curve,
          y_curve,
          n = length(x_curve),
        )
      )

      dx <- diff(curve_regularized$x)
      dy <- diff(curve_regularized$y)

      df_fp <- rbind(data.frame(
        x = curve_regularized$x[-length(curve_regularized$x)],
        y = dy / dx,
        variable = colnames(x$x)[k]
      ), df_fp)
    }

    if (type == "curve") {
      plt <- ggplot2::ggplot(df)

      if (!is.null(x$results[[k]]$envelope_data)) {
        envelope_ribbon <- x$results[[k]]$envelope_data
        envelope_ribbon <- dplyr::group_by(.data = envelope_ribbon, x)
        envelope_ribbon <- dplyr::summarise(
          .data = envelope_ribbon,
          ymin = min(y, na.rm = TRUE), ymax = max(y, na.rm = TRUE)
        )
        envelope_ribbon <- dplyr::filter(
          .data = envelope_ribbon,
          is.finite(ymin) & is.finite(ymax)
        )

        plt <- plt +
          ggplot2::geom_ribbon(
            data = envelope_ribbon,
            ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
            alpha = 0.7, fill = "grey"
          )
      }


      plt <- plt + ggplot2::geom_step(ggplot2::aes(x = alpha, y = geom_corr),
        linewidth = 1
      )

      for (k in 1:nvar) {
        plt <- plt +
          ggplot2::geom_function(
            data = subset(df, variable == colnames(x$x)[k]),
            fun = function(x, intensity) {
              exp(-intensity * pi * x^2)
            },
            args = list(intensity = x$results[[k]]$intensity),
            linetype = "dashed",
            color = "red",
            linewidth = 1
          )
      }



      plt <- plt +
        ggplot2::scale_y_continuous(name = expression(f(alpha))) +
        ggplot2::scale_x_continuous(
          name = expression(alpha),
          labels = scales::comma
        ) +
        ggplot2::facet_wrap(. ~ variable, scales = "free_x") +
        cowplot::theme_cowplot(font_size = font_size) +
        cowplot::background_grid(minor = "y") +
        cowplot::panel_border()
    } else if (type == "deriv") {
      plt <- ggplot2::ggplot(data = df_fp, ggplot2::aes(x, y)) +
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
        ggplot2::facet_wrap(. ~ variable, scales = "free_x") +
        cowplot::theme_cowplot(font_size = font_size) +
        cowplot::background_grid(minor = "y") +
        cowplot::panel_border()
    }

    return(plt)
  }
