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
    nvar <- length(x$results)
    df <- NULL
    df_null <- NULL
    df_f <- NULL
    df_fp <- NULL
    for (k in 1:nvar) {
      df <-  rbind(
        cbind(
          x$results[[k]]$data_frame_triangles,
          variable = colnames(x$Xdat)[k],
          intensity = x$results[[k]]$intensity
        ),
        df
      )

      # f <- smooth.spline(
      #   x = x$results[[k]]$data_frame_triangles$alpha,
      #   y =  x$results[[k]]$data_frame_triangles$geom_corr
      #   # tol = x$results[[k]]$tol
      # )
      #
      # df_f <-
      #   rbind(data.frame(
      #     x = f$x,
      #     y = f$y,
      #     variable = colnames(x$Xdat)[k]
      #   ), df_f)

      # fp <- predict(f, deriv = 1)
      # df_fp <-
      #   rbind(data.frame(
      #     x = fp$x,
      #     y = fp$y,
      #     variable = colnames(x$Xdat)[k]
      #   ), df_fp)

      df_fp <-
        rbind(data.frame(
          x = x$results[[k]]$data_frame_triangles$alpha[-1],
          y = diff(x$results[[k]]$data_frame_triangles$geom_corr)/diff(x$results[[k]]$data_frame_triangles$alpha),
          variable = colnames(x$Xdat)[k]
        ), df_fp)

    }


    # df_R2 <- as.data.frame(cbind(x$Ydat, x$Xdat)) %>%
    #   pivot_longer(cols = 2:(ncol(x$Xdat) + 1)) %>%
    #   split(.$name) %>%
    #   map(~ mgcv::gam(Y ~ s(value), data = .)) %>%
    #   map(summary) %>%
    #   map_df("r.sq")     %>%
    #   pivot_longer(cols = everything(),
    #                names_to = "variable",
    #                values_to = "R2") %>%
    #   mutate(max_alpha=max(df$alpha))


      # left_join(df %>% group_by(variable) %>% summarise(max_alpha = max(alpha)))





    if (type == "curve") {
      # for (k in 1:nvar) {
      #   f <- smooth.spline(
      #     x = x$results[[k]]$data_frame_triangles$alpha,
      #     y =  x$results[[k]]$data_frame_triangles$geom_corr
      #   )
      #
      #   plt <- plt +
      #     ggplot2::geom_line(data=data.frame(x=f$x,y=f$y), aes(x,y), size = 1.1)
      #
      #     #
      #     # ggplot2::geom_function(
      #     #   data = transform(df, variable = colnames(x$Xdat)[k]),
      #     #   fun = f,
      #     #   size = 1.1
      #     # )
      # }

      plt <- ggplot2::ggplot(df) +
        # ggplot2::geom_hline(
        #   data = df_R2,
        #   mapping = aes(yintercept = R2),
        #   color = "red",
        #   size = 0.5
        # ) +
        # ggplot2::geom_label(
        #   data = df_R2,
        #   mapping = aes(x = max_alpha, y =  1, label = paste("Stat. R^2: ", scales::number(R2,accuracy = 0.01))),
        #   color = "red",
        #   size = 3, vjust = 1,hjust = 1
        # ) +
        ggplot2::geom_step(ggplot2::aes(x = alpha, y = geom_corr),
                           # color = "black",
                            size = 1)
        # ggplot2::geom_line(data = df_f, ggplot2::aes(x, y), size = 0.5)

      for (k in 1:nvar) {
        plt <- plt +
          ggplot2::geom_function(
            data = subset(df, variable == colnames(Xdat)[k]),
            fun = function(x, intensity) {
              exp(-intensity * pi * x ^ 2)
            },
            args = list(intensity = x$results[[k]]$intensity),
            linetype = "dashed",
            color = "red", size=1
          )
      }


     plt <- plt +
        ggplot2::scale_y_continuous(name = expression(f(alpha))) +
        ggplot2::scale_x_continuous(name = expression(alpha), labels = scales::comma) +
        ggplot2::facet_wrap(. ~ variable, scales = "free_x") +
        cowplot::theme_cowplot(font_size = font_size) +
        cowplot::background_grid(minor = "y") +
        cowplot::panel_border()
      # ggplot2::theme(panel.grid.minor.y = ggplot2::element_line(color = "grey85", linetype = "dashed"))


    } else if (type == "deriv") {
      # for (k in 1:nvar) {
      #   f <- smooth.spline(
      #     x = x$results[[k]]$data_frame_triangles$alpha,
      #     y =  x$results[[k]]$data_frame_triangles$geom_corr
      #   )
      #
      #   fp <- predict(f, deriv = 1)
      #
      #   plt <- plt +
      #     ggplot2::geom_line(data = data.frame(x = fp$x, y = fp$y), aes(x, y), size = 1.1)
      #
      #
      #     # ggplot2::geom_function(
      #     #   data = transform(df, variable = colnames(x$Xdat)[k]),
      #     #   fun = f,
      #     #   size= 1.1,
      #     #   args = list(deriv = 1)
      #     # )
      # }

      plt <- ggplot2::ggplot(df) +
        ggplot2::geom_line(data = df_fp, ggplot2::aes(x, y), size = 1) +
        ggplot2::scale_y_continuous(name = expression(f * minute(alpha))) +
        ggplot2::scale_x_continuous(name = expression(alpha), labels = scales::comma) +
        ggplot2::facet_wrap(. ~ variable, scales = "free_x") +
        cowplot::theme_cowplot(font_size = font_size) +
        cowplot::background_grid(minor = "y") +
        cowplot::panel_border()
      # ggplot2::theme(panel.grid.minor.y = ggplot2::element_line(color = "grey85", linetype = "dashed"))

    }

    return(plt)

    # plt_sens <- ggplot2::ggplot(df) +
    #   ggplot2::geom_point(aes(x = alpha, y = geom_sens),
    #                       color = "grey50",
    #                       size = 2) +
    #   scale_y_continuous(labels = scales::percent) +
    #   ggplot2::ylim(c(0, 1))
    #
    # for (k in 1:nvar) {
    #   f <-
    #     splinefun(
    #       x$results[[k]]$data_frame_triangles$alpha,
    #       x$results[[k]]$data_frame_triangles$geom_sens,
    #       method = "natural"
    #     )
    #   plt_sens <- plt_sens +
    #     ggplot2::geom_function(data = transform(df, variable =  colnames(x$Xdat)[k]),
    #                            fun = f)
    # }
    #
    # plt_sens <- plt_sens +
    #   ggplot2::facet_wrap(. ~ variable) +
    #   ggplot2::labs(title = "Geometric Sensitivity",
    #                 x = "Alpha radius",
    #                 y = "?") +
    #   ggplot2::theme_minimal()


  }


#' @export
plot_alpha_shape <- function(x, alpha, font_size = 12) {
  nvar <- length(x$results)
  df <- NULL
  df_whole <- NULL

  for (k in 1:nvar) {
    df_triangles <- x$results[[k]]$triangles %>%
      dplyr::filter(max_length < 2 * alpha[k])

    df <-  rbind(cbind(df_triangles,
                       variable = colnames(x$Xdat)[k]), df)

  }

  ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = df,
      mapping = aes(geometry = geometry),
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



