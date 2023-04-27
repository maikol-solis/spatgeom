#' @title CSR Hypothesis Testing and Global Envelope for a CSR process using
#' \code{spatgeom} objects
#'
#' @description Calculates a table with the results of the CSR Hypothesis Test
#' with the MAD and DCFL methods and a table of a Global Envelope for a CSR
#' process, to determine whether the null hypothesis (the process is a CSR),
#' is rejected or not.
#'
#' @note This function requires the functions of alphastats.R and helpers.R
#' #preloaded, from the Geometric Spatial Point Analysis package.
#'
#' @param X an object of class data.frame of n dimensions representing the
#' independent variables of the process.
#' @param Y an object of class data.frame of one dimension representing the
#' dependent variable of the process.
#' @param sign_level a number of class numeric representing the significance
#' level of the test.
#' @param name_ex a string that specify the file name where the results are
#' saved.
#' @param r a number that relates the approximate theoretical function of the
#'  empty space in a CSR process with the value alpha of the alpha shape.
#'
#' @return a list of two objects of calls dataframe with the data of the
#' global envelope.  Also it saves a file of type png with the results of
#' the CSR Hypothesis Testing for both methods: MAD and DCLF.
#' These dataframes includes:
#' f_obs corresponding with the map of  Geometric Goodness of fit measure,
#' x that represents the alpha of the alpha shape,
#' mean that represents the map of the mean of the CSR simulated patterns,
#' theor that represents the theoretical map of the CSR pattern,
#' upper mean and lower mean that represents the upper and lower boundary
#' of the envelope using the mean of the CSR simulated patterns,
#' upper theo and lower theo that represents the upper and lower boundary
#' of the envelope using the theoretical map of the CSR simulated patterns,
#' variable indicating the name of the independent variable of the process.
#'
#' @example
#' n <- 100
#' a <- -1
#' b <- 1
#' theta <- runif(n, 0, 2 * pi)
#' r <- (sqrt(runif(n))) * (0.5) + 0.5
#' X1 <- r * cos(theta)
#' X2 <- runif(n, a, b)
#' Y <- data.frame(Y = r * sin(theta))
#' X <- data.frame(X1, X2)
#'
#' f_df <- CSR_test(X, Y, sign_level = 0.1, name_ex = "ring_100pts", r=0.5)
#'
#' @export

CSR_test <- function(X, Y, sign_level = 0.05, name_ex = "ex1", r = 0.5) {
  estimation <- alphastats(
    y = Y,
    x = X,
    envelope = TRUE,
    sign_level = sign_level
  )
  #x: the alpha of the Alpha Shape.
  #alphastats is part of the Geometric Spatial Point Pattern Analysis module.
  methods <- c("MAD", "DCFL")
  for (method in methods) {
    for (i in 1:ncol(X)) {
      f_simul <- estimation$results[[i]]$envelope_data
      #f_simul is the simulated maps for a CSR process with the same intensity
      #and observation window that the analyzed process.

      f_obs <- estimation$results[[i]]$data_frame_triangles[, 2]
      #f_obs: the observed map, f(alpha), of the Geometric Goodness of
      #fit measure.

      f_obs_with_num = data.frame(
        "x" = estimation$results[[i]]$data_frame_triangles[, 1],
        "y" = f_obs,
        "nsim" = 0
      )
      #f_obs_with_num: f_obs with the corresponding number of simulation.

      f_simul_and_obs <- rbind(f_simul, f_obs_with_num)

      f_mean <- f_simul_and_obs %>%
        group_by(x) %>%
        summarise(mean = mean(y, na.rm = TRUE))
      #fmean: map of the mean of the CSR simulated patterns.

      f_dataframe = cbind(f_obs, f_mean)
      intensity = estimation$results[[i]]$intensity
      #intensity: intensity of the point process.

      f_dataframe$theor = exp(-intensity * pi * (f_dataframe$x * r)^2)
      #f_dataframe$theor: the theoretical map of a CSR pattern.

      #t_simulated_mean: Value of the Test Statistic for each simulation,
      #using f_mean.
      #t_max_mean: maximum value of t_simul_mean.
      #t_simulated_theo: Value of the Test Statistic for each simulation,
      #using f_dataframe$theor.
      #t_max_theo: maximum value of t_simul_theo.

      if (method == "DCFL") {
        delta_alpha = f_mean$x[2] - f_mean$x[1]
        f_dataframe$sq_diff_mean <- (f_dataframe$f_obs - f_dataframe$mean)^2 *
          delta_alpha
        f_dataframe$sq_diff_theo <- (f_dataframe$f_obs - f_dataframe$theor)^2 *
          delta_alpha
        t_obs_mean = sum(f_dataframe$sq_diff_mean)
        t_obs_theo = sum(f_dataframe$sq_diff_theo)
        times <- 1 / sign_level - 1
        f_simul$f_mean_rep <- rep(f_mean$mean, times = times)
        f_simul$f_theor_rep <- rep(f_dataframe$theor, times = times)
        f_simul$sq_diff_mean <- (f_simul$f_mean_rep - f_simul$y)^2 * delta_alpha
        f_simul$sq_diff_theo <- (f_simul$f_theor_rep - f_simul$y)^2 *
          delta_alpha
        t_simul_mean <- f_simul %>%
          group_by(nsim) %>%
          summarise(maxim = sum(sq_diff_mean, na.rm = TRUE))
        t_simul_theo <- f_simul %>%
          group_by(nsim) %>%
          summarise(maxim = sum(sq_diff_theo, na.rm = TRUE))
        t_max_mean = max(t_simul_mean$maxim)
        t_max_theo = max(t_simul_theo$maxim)
      } else {
        f_dataframe$abs_diff_mean <- abs(f_dataframe$f_obs - f_dataframe$mean)
        f_dataframe$abs_diff_theo <- abs(f_dataframe$f_obs - f_dataframe$theor)
        t_obs_mean = max(f_dataframe$abs_diff_mean)
        t_obs_theo = max(f_dataframe$abs_diff_theo)
        times <- 1 / sign_level - 1
        f_simul$f_mean_rep <- rep(f_mean$mean, times = times)
        f_simul$f_theor_rep <- rep(f_dataframe$theor, times = times)
        f_simul$abs_diff_mean <- abs(f_simul$f_mean_rep - f_simul$y)
        f_simul$abs_diff_theo <- abs(f_simul$f_theor_rep - f_simul$y)
        t_simul_mean <- f_simul %>%
          group_by(nsim) %>%
          summarise(maxim = max(abs_diff_mean, na.rm = TRUE))
        t_simul_theo <- f_simul %>%
          group_by(nsim) %>%
          summarise(maxim = max(abs_diff_theo, na.rm = TRUE))
        t_max_mean = max(t_simul_mean$maxim)
        t_max_theo = max(t_simul_theo$maxim)
      }
      t_simul_ord_mean = sort(t_simul_mean$maxim, decreasing = TRUE)
      m <- length(t_simul_ord_mean)
      k <- if (t_obs_mean < min(t_simul_ord_mean)) {
        m + 1
      } else {
        which(t_obs_mean > t_simul_ord_mean)[1]
      }
      t_simul_ord_theo = sort(t_simul_theo$maxim, decreasing = TRUE)
      m2 <- length(t_simul_ord_theo)
      k2 <- if (t_obs_theo < min(t_simul_ord_theo)) {
        m2 + 1
      } else {
        which(t_obs_theo > t_simul_ord_theo)[1]
      }

      results <- data.frame(
        "Feature" = c(
          "H_0",
          "Variable",
          "Type of Test",
          "Number of Monte Carlo Simulations",
          "T_max_mean",
          "T_obs_mean",
          "T_obs_mean > T_max_mean",
          "p-value_mean",
          "T_max_theo",
          "T_obs_theo",
          "T_obs_theo > T_max_theo",
          "p-value_theo"
        ),
        "Value" = c(
          "CSR process",
          colnames(X)[i],
          method,
          m,
          round(t_max_mean, 4),
          round(t_obs_mean, 4),
          t_obs_mean > t_max_mean,
          k / (m + 1),
          round(t_max_theo, 4),
          round(t_obs_theo, 4),
          t_obs_theo > t_max_theo,
          k2 / (m2 + 1)
        )
      )

      png(
        paste0("results_", name_ex, "_", method, "_", colnames(X)[i], ".png"),
        height = 50 * nrow(results),
        width = 200 * ncol(results)
      )
      grid.table(results)
      dev.off()
      print(results)
      writeLines("\n")
      f_dataframe$upper_mean <- f_dataframe$mean + t_max_mean
      f_dataframe$lower_mean <- f_dataframe$mean - t_max_mean
      f_dataframe$upper_theo <- f_dataframe$theor + t_max_theo
      f_dataframe$lower_theo <- f_dataframe$theor - t_max_theo
      f_dataframe$variable = colnames(X[i])
      if (method == "MAD") {
        if (i == 1) {
          f_df_MAD <- f_dataframe
        } else {
          f_df_MAD <- rbind(f_df_MAD, f_dataframe)
        }
      } else {
        if (i == 1) {
          f_df_DCLF <- f_dataframe
        } else {
          f_df_DCFL <- rbind(f_df_DCLF, f_dataframe)
        }
      }
    }
  }
  f_df <- list(f_df_MAD = f_df_MAD, f_df_DCLF = f_df_DCLF)
  return(f_df)
}
