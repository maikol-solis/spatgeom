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
#' n <- 30
#' a <- -1
#' b <- 1
#' theta <- runif(n, 0, 2 * pi)
#' r <- (sqrt(runif(n))) * (0.5) + 0.5
#' X1 <- r * cos(theta)
#' X2 <- runif(n, a, b)
#' Y <- data.frame(Y = r * sin(theta))
#' X <- data.frame(X1, X2)
#'
#' estimation <- spatgeom(y = Y, x = X)
#'
#' plot_curve(estimation, type = "curve")
#'
#' plot_curve(estimation, type = "deriv")
#' @export

plot_curve <-
  function(x,
           type = "curve",
           font_size = 12) {
    # Number of variables
    nvar <- length(x$results)

    # variables bindings for R CMD check
    df <- df_fp <- alpha <- geom_corr <- variable <- y <- nsim <- NULL

    for (k in 1:nvar) {
      df <- rbind(
        cbind(
          x$results[[k]]$data_frame_triangles,
          variable = colnames(x$x)[k],
          intensity = x$results[[k]]$intensity
        ),
        df
      )


      df_fp <- rbind(data.frame(
        x = x$results[[k]]$data_frame_triangles$alpha[-1],
        y = diff(x$results[[k]]$data_frame_triangles$geom_corr) /
          diff(x$results[[k]]$data_frame_triangles$alpha),
        variable = colnames(x$x)[k]
      ), df_fp)
    }

    if (type == "curve") {
      plt <- ggplot2::ggplot(df) +
        ggplot2::geom_step(ggplot2::aes(x = alpha, y = geom_corr),
          size = 1
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
            color = "red", size = 1
          )

        if (!is.null(x$results[[k]]$envelope_data)) {
          plt <- plt +
            ggplot2::geom_line(
              data = x$results[[k]]$envelope_data,
              mapping = ggplot2::aes(x, y, group = nsim), color = "lightgrey"
            )
        }
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
      plt <- ggplot2::ggplot(df) +
        ggplot2::geom_line(data = df_fp, ggplot2::aes(x, y), size = 1) +
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
