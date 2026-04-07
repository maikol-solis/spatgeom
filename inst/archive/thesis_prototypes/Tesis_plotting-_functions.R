# Functions for graphs.

# Example
n <- 100
a <- -1
b <- 1
theta <- runif(n, 0, 2 * pi)
r <- (sqrt(runif(n))) * (0.5) + 0.5
X1 <- r * cos(theta)
X2 <- runif(n, a, b)
Y <- data.frame(Y = r * sin(theta))
X <- data.frame(X1, X2)

library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpmisc)
library(latex2exp)


variable_names <- list(
  "X1" = TeX("$X_1$"),
  "X2" = TeX("$X_2$"),
  "X3" = TeX("$X_3$")
)


variable_labeller <- function(variable, value) {
  return(variable_names[value])
}


plot_clouds <- function(X, Y, namefile.png) {
  # Plots the data point clouds.
  df <- cbind(X, Y)
  df <- tidyr::pivot_longer(
    df,
    cols = -Y,
    names_to = "variable",
    values_to = "points"
  )
  g <- ggplot(df, aes(points, Y)) +
    ggplot2::geom_point(size = 0.75) +
    stat_poly_line(se = FALSE) +
    # stat_poly_eq(label.y = 1.1) +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border() +
    ggtitle("") +
    xlab("X") +
    ylab("Y") +
    facet_wrap(. ~ variable, scales = "free", labeller = variable_labeller)

  cowplot::save_plot(
    namefile.png,
    g,
    bg = "white",
    base_height = 2.50,
    base_asp = 1.618,
    dpi = 300
  )
  return(g)
}

plot_clouds_R <- function(X, Y, namefile.png) {
  # Plots the data point clouds with the determination coeff. R^2.
  df <- cbind(X, Y)
  df <- tidyr::pivot_longer(
    df,
    cols = -Y,
    names_to = "variable",
    values_to = "points"
  )
  g <- ggplot(df, aes(points, Y)) +
    ggplot2::geom_point(size = 0.75) +
    stat_poly_line(se = FALSE) +
    stat_poly_eq(label.y = 1.1) +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border() +
    ggtitle("") +
    xlab("X") +
    ylab("Y") +
    facet_wrap(. ~ variable, scales = "free", labeller = variable_labeller)

  cowplot::save_plot(
    namefile.png,
    g,
    bg = "white",
    base_height = 2.50,
    base_asp = 1.618,
    dpi = 300
  )
  return(g)
}

plot_clouds_cuad <- function(X, Y, namefile.png) {
  # Plots the data point clouds with aspect ratio = 1.
  df <- cbind(X, Y)
  df <- tidyr::pivot_longer(
    df,
    cols = -Y,
    names_to = "variable",
    values_to = "points"
  )
  g <- ggplot(df, aes(points, Y)) +
    ggplot2::geom_point(size = 0.75) +
    stat_poly_line(se = FALSE) +
    # stat_poly_eq(label.y = 1.1) +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border() +
    ggtitle("") +
    xlab("X") +
    ylab("Y") +
    # coord_fixed( ratio=1)
    facet_wrap(. ~ variable, scales = "free_x", labeller = variable_labeller) +
    theme(aspect.ratio = 1)

  cowplot::save_plot(
    namefile.png,
    g,
    bg = "white",
    base_height = 2.50,
    base_asp = 1.618,
    dpi = 300
  )
  return(g)
}


S_df <- CSR_test(X, Y, sign_level = 0.1, name_ex = "ring_100pts_MAD")

variable_names <- list(
  "X1" = TeX("$X_1$"),
  "X2" = TeX("$X_2$"),
  "X3" = TeX("$X_3$")
)


variable_labeller <- function(variable, value) {
  return(variable_names[value])
}


plot_envelope_mean <- function(df, font_size = 12, name_plot) {
  # Plots the CSR envelope using the f_mean (Mean of the CSR simulated patterns)
  g <- ggplot(df, aes(x)) +
    geom_line(aes(y = f_obs), colour = "black") +
    geom_line(aes(y = mean), colour = "blue") +
    geom_line(aes(y = upper_mean), colour = "red", linetype = "dashed") +
    geom_line(aes(y = lower_mean), colour = "red", linetype = "dashed") +
    ggtitle("") +
    xlab(TeX("$\\alpha$")) +
    ylab(TeX("f($\\alpha$)")) +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border() +
    facet_wrap(. ~ variable, scales = "free", labeller = variable_labeller)
  g
  cowplot::save_plot(
    paste0(name_plot, ".png"),
    g,
    bg = "white",
    base_height = 2.50,
    base_asp = 1.618,
    dpi = 300
  )
  return(g)
}


plot_envelope_theo <- function(df, font_size = 12, name_plot) {
  # Plots the CSR envelope using the theoretical map of a CSR process.
  g <- ggplot(df, aes(x)) +
    geom_line(aes(y = f_obs), colour = "black") +
    geom_line(aes(y = upper_theo), colour = "red", linetype = "dashed") +
    geom_line(aes(y = lower_theo), colour = "red", linetype = "dashed") +
    geom_line(aes(y = theor), colour = "green", linetype = "solid") +
    ggtitle("") +
    xlab(TeX("$\\alpha$")) +
    ylab(TeX("f($\\alpha$)")) +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border() +
    facet_wrap(. ~ variable, scales = "free", labeller = variable_labeller)
  g
  cowplot::save_plot(
    paste0(name_plot, ".png"),
    g,
    bg = "white",
    base_height = 2.50,
    base_asp = 1.618,
    dpi = 300
  )
  return(g)
}

plot_envelope_theo_2 <- function(df, font_size = 12, name_plot) {
  # Plots the CSR envelope using the theoretical map of a CSR process.
  g <- ggplot(df, aes(x)) +
    geom_line(aes(y = f_obs), colour = "black") +
    geom_line(aes(y = upper_theo), colour = "red", linetype = "dashed") +
    geom_line(aes(y = lower_theo), colour = "red", linetype = "dashed") +
    geom_line(aes(y = theor), colour = "green", linetype = "solid") +
    ggtitle("") +
    xlab(TeX("$\\alpha$")) +
    ylab(TeX("f($\\alpha$)")) +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border() +
    facet_wrap(. ~ variable, scales = "free", labeller = variable_labeller) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  g
  cowplot::save_plot(
    paste0(name_plot, ".png"),
    g,
    bg = "white",
    base_height = 2.50,
    base_asp = 1.618,
    dpi = 300
  )
  return(g)
}

plot_envelope_mean_2 <- function(df, font_size = 12, name_plot) {
  # Plots the CSR envelope using the f_mean (Mean of the CSR simulated patterns)
  g <- ggplot(df, aes(x)) +
    geom_line(aes(y = f_obs), colour = "black") +
    geom_line(aes(y = mean), colour = "blue") +
    geom_line(aes(y = upper_mean), colour = "red", linetype = "dashed") +
    geom_line(aes(y = lower_mean), colour = "red", linetype = "dashed") +
    ggtitle("") +
    xlab(TeX("$\\alpha$")) +
    ylab(TeX("f($\\alpha$)")) +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border() +
    facet_wrap(. ~ variable, scales = "free", labeller = variable_labeller) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  g
  cowplot::save_plot(
    paste0(name_plot, ".png"),
    g,
    bg = "white",
    base_height = 2.50,
    base_asp = 1.618,
    dpi = 300
  )
  return(g)
}


plot_envelope_mean_real <- function(df, font_size = 12, name_plot) {
  # Plots the CSR envelope using the f_mean (Mean of the CSR simulated patterns)
  g <- ggplot(df, aes(x)) +
    geom_line(aes(y = f_obs), colour = "black") +
    geom_line(aes(y = mean), colour = "blue") +
    geom_line(aes(y = upper_mean), colour = "red", linetype = "dashed") +
    geom_line(aes(y = lower_mean), colour = "red", linetype = "dashed") +
    ggtitle("") +
    xlab(TeX("$\\alpha$")) +
    ylab(TeX("f($\\alpha$)")) +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border() +
    facet_wrap(. ~ variable, scales = "free")
  g
  cowplot::save_plot(
    paste0(name_plot, ".png"),
    g,
    bg = "white",
    base_height = 2.50,
    base_asp = 1.618,
    dpi = 300
  )
  return(g)
}


plot_envelope_theo_real <- function(df, font_size = 12, name_plot) {
  # Plots the CSR envelope using the theoretical map of a CSR process.
  g <- ggplot(df, aes(x)) +
    geom_line(aes(y = f_obs), colour = "black") +
    geom_line(aes(y = upper_theo), colour = "red", linetype = "dashed") +
    geom_line(aes(y = lower_theo), colour = "red", linetype = "dashed") +
    geom_line(aes(y = theor), colour = "green", linetype = "solid") +
    ggtitle("") +
    xlab(TeX("$\\alpha$")) +
    ylab(TeX("f($\\alpha$)")) +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border() +
    facet_wrap(. ~ variable, scales = "free")
  g
  cowplot::save_plot(
    paste0(name_plot, ".png"),
    g,
    bg = "white",
    base_height = 2.50,
    base_asp = 1.618,
    dpi = 300
  )
  return(g)
}

plot_envelope_theo_real_2 <- function(df, font_size = 12, name_plot) {
  # Plots the CSR envelope using the theoretical map of a CSR process.
  g <- ggplot(df, aes(x)) +
    geom_line(aes(y = f_obs), colour = "black") +
    geom_line(aes(y = upper_theo), colour = "red", linetype = "dashed") +
    geom_line(aes(y = lower_theo), colour = "red", linetype = "dashed") +
    geom_line(aes(y = theor), colour = "green", linetype = "solid") +
    ggtitle("") +
    xlab(TeX("$\\alpha$")) +
    ylab(TeX("f($\\alpha$)")) +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border() +
    facet_wrap(. ~ variable, scales = "free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text.x = element_text(size = 8))
  g
  cowplot::save_plot(
    paste0(name_plot, ".png"),
    g,
    bg = "white",
    base_height = 2.50,
    base_asp = 1.618,
    dpi = 300
  )
  return(g)
}


plot_f_i <- function(df, font_size = 12, name_plot) {
  g <- ggplot(df, aes(x)) +
    geom_line(aes(y = f_obs), colour = "black") +
    ggtitle("") +
    xlab(TeX("$\\alpha$")) +
    ylab(TeX("f($\\alpha$)")) +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border() +
    facet_wrap(. ~ variable, scales = "free")
  g
  cowplot::save_plot(
    paste0(name_plot, ".png"),
    g,
    bg = "white",
    base_height = 4.0,
    base_asp = 1.618,
    dpi = 300
  )
  return(g)
}


plot_clouds_real <- function(X, Y, namefile.png, y_title) {
  df <- cbind(X, Y)
  df <- tidyr::pivot_longer(
    df,
    cols = -Y,
    names_to = "variable",
    values_to = "points"
  )
  g <- ggplot(df, aes(points, Y)) +
    ggplot2::geom_point(size = 0.75) +
    stat_poly_line(se = FALSE) +
    # stat_poly_eq() +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border() +
    ggtitle("") +
    xlab("X") +
    ylab(y_title) +
    facet_wrap(. ~ variable, scales = "free")
  cowplot::save_plot(
    namefile.png,
    g,
    bg = "white",
    base_height = 2.50,
    base_asp = 1.618,
    dpi = 300
  )
  return(g)
}

plot_clouds_real_2 <- function(X, Y, namefile.png, y_title) {
  df <- cbind(X, Y)
  df <- tidyr::pivot_longer(
    df,
    cols = -Y,
    names_to = "variable",
    values_to = "points"
  )
  g <- ggplot(df, aes(points, Y)) +
    ggplot2::geom_point(size = 0.75) +
    stat_poly_line(se = FALSE) +
    # stat_poly_eq() +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border() +
    ggtitle("") +
    xlab("X") +
    ylab(y_title) +
    facet_wrap(. ~ variable, scales = "free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text.x = element_text(size = 8))
  cowplot::save_plot(
    namefile.png,
    g,
    bg = "white",
    base_height = 2.50,
    base_asp = 1.618,
    dpi = 300
  )
  return(g)
}


plot_clouds_real_R <- function(X, Y, namefile.png, y_title) {
  df <- cbind(X, Y)
  df <- tidyr::pivot_longer(
    df,
    cols = -Y,
    names_to = "variable",
    values_to = "points"
  )
  g <- ggplot(df, aes(points, Y)) +
    ggplot2::geom_point(size = 0.75) +
    stat_poly_line(se = FALSE) +
    stat_poly_eq() +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border() +
    ggtitle("") +
    xlab("X") +
    ylab(y_title) +
    facet_wrap(. ~ variable, scales = "free")
  cowplot::save_plot(
    namefile.png,
    g,
    bg = "white",
    base_height = 2.50,
    base_asp = 1.618,
    dpi = 300
  )
  return(g)
}


plot_clouds2 <- function(X, Y, namefile.png) {
  df <- cbind(X, Y)
  df <- tidyr::pivot_longer(
    df,
    cols = -Y,
    names_to = "variable",
    values_to = "points"
  )
  g <- ggplot(df, aes(points, Y)) +
    ggplot2::geom_point(size = 0.75) +
    # stat_poly_line() +
    # stat_poly_eq() +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    cowplot::panel_border() +
    ggtitle("") +
    xlab("X") +
    ylab("Y") +
    facet_wrap(. ~ variable, scales = "free")

  cowplot::save_plot(namefile.png, g, bg = "white")
  return(g)
}
