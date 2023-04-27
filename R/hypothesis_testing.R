#' CSR Hypothesis Testing and Global Envelope for a CSR process using
#' \code{spatgeom} objects
#'
#' @param X an object of class data.frame of n dimensions.
#' @param Y an object of class data.frame of one dimension.
#'
#' @param significance_level a number of class numeric representing the
#'   significance level of the test.
#'
#' @return a list of two objects of calls data.frame with the data of the
#' global envelope.  Also it saves a file of type png with the results of
#' the CSR Hypothesis Testing for both methods: MAD and DCLF.
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
#' f_df <- csr_test(X, Y, significance_level = 0.1, name_ex = "ring_100pts")
#'
#' @export



csr_test <- function(spatgeom_obj, significance_level = 0.05,
                     method = c("MAD", "DCLF"), name_ex = "ex1") {
  # global variables to NULL
  nsim <- y <- sum_sq_dif <- sum_sq_dif2 <-
    abs_dif <- abs_dif2 <- x <- NULL


  if (length(method) > 1 || !(method %in% c("MAD", "DCLF"))) {
    stop("Parameter method must be one of 'MAD' or 'DCLF'.")
  }


  x <- spatgeom_obj$x


  for (i in seq_len(ncol(x))) {
    env1 <- spatgeom_obj$results[[i]]$envelope_data
    f_obs1 <- spatgeom_obj$results[[i]]$geom_indices[, 2]

    f_obs_env <- data.frame(
      "x" = spatgeom_obj$results[[i]]$geom_indices[, 1],
      "y" = spatgeom_obj$results[[i]]$geom_indices[, 2],
      "nsim" = 0
    )

    env_obs <- rbind(env1, f_obs_env)

    f_doublehat <- env_obs |>
      dplyr::group_by(x) |>
      dplyr::summarise(mean = mean(y, na.rm = TRUE))

    f_dataframe <- cbind(f_obs1, f_doublehat)

    intensity <- spatgeom_obj$results[[i]]$intensity
    f_dataframe$theo <- exp(-intensity * pi * (f_dataframe$x / 2)^2)

    if (method == "DCLF") {
      delta_r <- f_doublehat$x[2] - f_doublehat$x[1]

      f_dataframe$sum_sq_diff <- (f_dataframe$f_obs1 -
        f_dataframe$mean)^2 * delta_r

      f_dataframe$sum_sq_diff2 <- (f_dataframe$f_obs1 -
        f_dataframe$theo)^2 * delta_r

      t_obs <- sum(f_dataframe$sum_sq_diff) # T_obs using f_mean
      t_obs2 <- sum(f_dataframe$sum_sq_diff2) # T_obs using f_theo

      env2 <- spatgeom_obj$results[[i]]$envelope_data
      times <- 1 / significance_level - 1
      env2$f_doublehat <- rep(f_doublehat$mean, times = times)
      env2$f_theo <- rep(f_dataframe$theo, times = times)


      env2$sum_sq_dif <- (env2$f_doublehat - env2$y)^2 * delta_r
      env2$sum_sq_dif2 <- (env2$f_theo - env2$y)^2 * delta_r

      # T_simul using f_mean
      t_simul <- env2 |>
        dplyr::group_by(nsim) |>
        dplyr::summarise(maxim = sum(sum_sq_dif, na.rm = TRUE))

      # T_simul using f_theo
      t_simul2 <- env2 |>
        dplyr::group_by(nsim) |>
        dplyr::summarise(maxim = sum(sum_sq_dif2, na.rm = TRUE))

      # T_max using f_mean
      tmax <- max(t_simul$maxim)
      # T_max using f_theo
      tmax2 <- max(t_simul2$maxim)
    } else {
      f_dataframe$abs_dif <- abs(f_dataframe$f_obs1 - f_dataframe$mean)

      f_dataframe$abs_dif2 <- abs(f_dataframe$f_obs1 - f_dataframe$theo)

      t_obs <- max(f_dataframe$abs_dif) # t_obs using f_mean
      t_obs2 <- max(f_dataframe$abs_dif2) # t_obs using f_theo


      env2 <- spatgeom_obj$results[[i]]$envelope_data
      times <- 1 / significance_level - 1
      env2$f_doublehat <- rep(f_doublehat$mean, times = times)

      env2$f_theo <- rep(f_dataframe$theo, times = times)

      env2$abs_dif <- abs(env2$f_doublehat - env2$y)

      env2$abs_dif2 <- abs(env2$f_theo - env2$y)

      # T_simul using f_mean
      t_simul <- env2 |>
        dplyr::group_by(nsim) |>
        dplyr::summarise(maxim = max(abs_dif, na.rm = TRUE))

      # t_simul using f_theo
      t_simul2 <- env2 |>
        dplyr::group_by(nsim) |>
        dplyr::summarise(maxim = max(abs_dif2, na.rm = TRUE))

      tmax <- max(t_simul$maxim) # T_max using f_mean
      tmax2 <- max(t_simul2$maxim) # T_max using f_theo
    }
    t_simul_ord <- sort(t_simul$maxim, decreasing = TRUE)
    m <- length(t_simul_ord)
    k <- ifelse(t_obs < min(t_simul_ord),
      m + 1,
      which(t_obs > t_simul_ord)[1]
    )

    t_simul_ord2 <- sort(t_simul2$maxim, decreasing = TRUE)
    m2 <- length(t_simul_ord2)
    k2 <- ifelse(t_obs2 < min(t_simul_ord2),
      m2 + 1,
      which(t_obs2 > t_simul_ord2)[1]
    )


    results <- data.frame(
      "Feature" = c(
        "H_0", "Variable", "Type of Test",
        "Number of Monte Carlo Simulations",
        "T_max", "t_obs", "t_obs > T_max",
        "p-value",
        "T_max2", "t_obs2", "t_obs2 > T_max2",
        "p-value2"
      ),
      "Value" = c(
        "CSR process", colnames(x)[i], method, m, round(tmax, 4),
        round(t_obs, 4), t_obs > tmax, k / (m + 1),
        round(tmax2, 4), round(t_obs2, 4), t_obs2 > tmax2, k2 / (m2 + 1)
      )
    )

    ## png(paste0("results_", name_ex, "_", method, "_", colnames(x)[i],
    ## ".png"), height = 50 * nrow(results), width = 200 * ncol(results))
    ## grid.table(results) dev.off()

    print(results)

    f_dataframe$upper <- f_dataframe$mean + tmax
    f_dataframe$lower <- f_dataframe$mean - tmax

    f_dataframe$upper2 <- f_dataframe$theo + tmax2
    f_dataframe$lower2 <- f_dataframe$theo - tmax2

    f_dataframe$variable <- colnames(x[i])

    if (method == "MAD") {
      if (i == 1) {
        f_df_mad <- f_dataframe
      } else {
        f_df_mad <- rbind(f_df_mad, f_dataframe)
      }
    } else {
      if (i == 1) {
        f_df_dclf <- f_dataframe
      } else {
        f_df_dclf <- rbind(f_df_dclf, f_dataframe)
      }
    }
  }
  f_df <- list(f_df1 = f_df_mad, f_df2 = f_df_dclf)
  return(f_df)
}
