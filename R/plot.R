#' plot \code{topsa} objects
#'
#' Plot method for objects of class \code{topsa}.
#'
#' @param topsaObj an object of class \code{topsa}
#' @param nvar  it could be a sequence from 1 to the number of variables
#'   indicating which variables should be plotted. It could be the character
#'   'all' for plot all the variables.
#' @param ... further arguments passed to the \code{plot} function
#'
#' @return A plot of generated with the output of \code{topsa}. For each
#' variable in the model, it creates the plot of the corresponding manifold, its
#' symmetric reflection and its symmetric difference.
#'
#' @export

plot_geom_curves <-
  function(x,
           type = c("curve", "deriv"),
           font_size = 12) {
    # Number of variables
    nvar <- length(x$results)

    # variables bindings for R CMD check
    df <- df_fp <- alpha <- geom_corr <- variable <- y <- NULL

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
            data = subset(df, variable == colnames(x)[k]),
            fun = function(x, intensity) {
              exp(-intensity * pi * x^2)
            },
            args = list(intensity = x$results[[k]]$intensity),
            linetype = "dashed",
            color = "red", size = 1
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


#' @export
plot_alpha_shape <- function(x, alpha, font_size = 12) {
  nvar <- length(x$results)
  df <- max_length <- geometry <- NULL

  for (k in 1:nvar) {
    df_triangles <- subset(x$results[[k]]$triangles, max_length < 2 * alpha[k])

    df <- rbind(cbind(df_triangles,
      variable = colnames(x$x)[k]
    ), df)
  }

  ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = df,
      mapping = ggplot2::aes(geometry = geometry),
      fill = "skyblue1",
      color = "grey30",
      size = 0.05
    ) +
    ggplot2::facet_wrap(. ~ variable) +
    cowplot::theme_cowplot(font_size = font_size) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border()

  # p <- ggplot2::ggplot(x$results[[1]]$triangles%>%
  #                   dplyr::mutate(alpha_reveal = 2 * max_length)) +
  #   ggplot2::geom_sf(
  #     mapping = aes(geometry = geometry),
  #     fill = "skyblue1",
  #     color = "grey30",
  #     size = 0.1
  #   ) +
  #   # ggplot2::facet_wrap(. ~ variable) +
  #   cowplot::theme_cowplot(font_size = font_size) +
  #   cowplot::background_grid(minor = "y") +
  #   cowplot::panel_border() +
  #   gganimate::transition_time(alpha_reveal)
  #
  # gganimate::animate(p)
}
