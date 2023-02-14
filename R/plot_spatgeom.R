#' Plot alpha-shape for \code{spatgeom} objects
#'
#' Plot alpha-shape for \code{spatgeom} objects.
#'
#' @param x an object of class \code{spatgeom}.
#' @param alpha value of \code{alpha} determining the maximum length between
#'   points to build the alpha-shape.
#' @param font_size a integer that increases the font size in the plot.
#'
#' @return a \code{\link[ggplot2]{ggplot}} object with the raw alpha-shape for
#'   the original data at resolution \code{alpha}
#'
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
#' estimation <- alphastats(y = Y, x = X)
#'
#' plot_alpha_shape(estimation, alpha = c(0.9, 1.2))
#' @export


plot_alpha_shape <- function(x, alpha, font_size = 12) {
  if (length(alpha) == 1) {
    alpha <- rep(alpha, ncol(x$x))
  } else if (length(alpha) > 1 && length(alpha) != ncol(x$x)) {
    stop("alpha must be a vector of size 1 or the number of variables.")
  }

  nvar <- length(x$results)
  df <- max_length <- geometry <- NULL

  for (k in 1:nvar) {
    triangles <- x$results[[k]]$triangles
    s <- alpha[k]
    df_triangles <- subset(triangles, max_length < s)

    df <- rbind(cbind(df_triangles,
      variable = colnames(x$x)[k]
    ), df)
  }

  df <- sf::st_sf(df)

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
