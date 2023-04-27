library(dplyr)
library(ggplot2)
library(gridExtra)
library(readxl)
library(readr)


#ring

n <- 30
a <- -1
b <- 1
theta <- runif(n, 0, 2 * pi)
r <- (sqrt(runif(n))) * (0.5) + 0.5
X1 <- r * cos(theta)
X2 <- runif(n, a, b)
Y <- data.frame(Y = r * sin(theta))
X <- data.frame(X1, X2)


plot_clouds(X, Y, "cloud_ring30_pts.png")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "ring_30pts")

#explorar
f_df

plot_envelope_theo(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_ring_30pts_MAD_theo"
)
plot_envelope_mean(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_ring_30pts_MAD_mean"
)


n <- 100
a <- -1
b <- 1
theta <- runif(n, 0, 2 * pi)
r <- (sqrt(runif(n))) * (0.5) + 0.5
X1 <- r * cos(theta)
X2 <- runif(n, a, b)
Y <- data.frame(Y = r * sin(theta))
X <- data.frame(X1, X2)

plot_clouds(X, Y, "cloud_ring100_pts.png")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "ring_100pts")
plot_envelope_theo(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_ring_100pts_MAD_theo"
)
plot_envelope_mean(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_ring_100pts_MAD_mean"
)

n <- 500
a <- -1
b <- 1
theta <- runif(n, 0, 2 * pi)
r <- (sqrt(runif(n))) * (0.5) + 0.5
X1 <- r * cos(theta)
X2 <- runif(n, a, b)
Y <- data.frame(Y = r * sin(theta))
X <- data.frame(X1, X2)

plot_clouds(X, Y, "cloud_ring500_pts.png")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "ring_500pts")
plot_envelope_theo(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_ring_500pts_MAD_theo"
)
plot_envelope_mean(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_ring_500pts_MAD_mean"
)


n <- 1000
a <- -1
b <- 1
theta <- runif(n, 0, 2 * pi)
r <- (sqrt(runif(n))) * (0.5) + 0.5
X1 <- r * cos(theta)
X2 <- runif(n, a, b)
Y <- data.frame(Y = r * sin(theta))
X <- data.frame(X1, X2)

plot_clouds(X, Y, "cloud_ring1000_pts.png")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "ring_1000pts")
plot_envelope_theo(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_ring_1000pts_MAD_theo"
)
plot_envelope_mean(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_ring_1000pts_MAD_mean"
)

#Example varying significance level

#Ring 500 points
n <- 500
a <- -1
b <- 1
theta <- runif(n, 0, 2 * pi)
r <- (sqrt(runif(n))) * (0.5) + 0.5
X1 <- r * cos(theta)
X2 <- runif(n, a, b)
Y <- data.frame(Y = r * sin(theta))
X <- data.frame(X1, X2)

sign_level_vector = c(0.01, 0.05, 0.1, 0.2)

for (sign_level in sign_level_vector) {
  f_df <- CSR_test(
    X,
    Y,
    sign_level = sign_level,
    name_ex = paste("ring_500pts", sign_level, sep = "_")
  )
  plot_envelope_theo(
    f_df$f_df_MAD,
    font_size = 12,
    name_plot = paste("curves_ring_500pts_MAD_theo", sign_level, sep = "_")
  )
  plot_envelope_mean(
    f_df$f_df_MAD,
    font_size = 12,
    name_plot = paste("curves_ring_500pts_MAD_mean", sign_level, sep = "_")
  )
}


#Linear

n <- 30
a <- -2
b <- 2
X1 <- runif(n, a, b)
X2 <- runif(n, a, b)
X3 <- runif(n, a, b)
Y <- 0.6 * X1 + 0.3 * X2 + 0.1 * X3
X <- data.frame(X1, X2, X3)

plot_clouds(X, Y, "cloud_linear30_pts.png")
plot_clouds_R(X, Y, "R_cloud_linear30_pts.png")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "linear_30pts")

plot_envelope_theo_2(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_linear_30pts_MAD_theo"
)
plot_envelope_mean_2(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_linear_30pts_MAD_mean"
)


n <- 100
a <- -2
b <- 2
X1 <- runif(n, a, b)
X2 <- runif(n, a, b)
X3 <- runif(n, a, b)
Y <- 0.6 * X1 + 0.3 * X2 + 0.1 * X3
X <- data.frame(X1, X2, X3)

plot_clouds(X, Y, "cloud_linear100_pts.png")
plot_clouds_R(X, Y, "R_cloud_linear100_pts.png")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "linear_100pts")

plot_envelope_theo_2(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_linear_100pts_MAD_theo"
)
plot_envelope_mean_2(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_linear_100pts_MAD_mean"
)

n <- 500
a <- -2
b <- 2
X1 <- runif(n, a, b)
X2 <- runif(n, a, b)
X3 <- runif(n, a, b)
Y <- 0.6 * X1 + 0.3 * X2 + 0.1 * X3
X <- data.frame(X1, X2, X3)

plot_clouds(X, Y, "cloud_linear500_pts.png")
plot_clouds_R(X, Y, "R_cloud_linear500_pts.png")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "linear_500pts")

plot_envelope_theo_2(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_linear_500pts_MAD_theo"
)
plot_envelope_mean_2(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_linear_500pts_MAD_mean"
)

n <- 1000
a <- -2
b <- 2
X1 <- runif(n, a, b)
X2 <- runif(n, a, b)
X3 <- runif(n, a, b)
Y <- 0.6 * X1 + 0.3 * X2 + 0.1 * X3
X <- data.frame(X1, X2, X3)

plot_clouds(X, Y, "cloud_linear1000_pts.png")
plot_clouds_R(X, Y, "R_cloud_linear1000_pts.png")

f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "linear_1000pts")

plot_envelope_theo_2(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_linear_1000pts_MAD_theo"
)
plot_envelope_mean_2(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_linear_1000pts_MAD_mean"
)


#Example varying significance level
#Linear 500 points
n <- 500
a <- -2
b <- 2
X1 <- runif(n, a, b)
X2 <- runif(n, a, b)
X3 <- runif(n, a, b)
Y <- 0.6 * X1 + 0.3 * X2 + 0.1 * X3
X <- data.frame(X1, X2, X3)

for (sign_level in sign_level_vector) {
  f_df <- CSR_test(
    X,
    Y,
    sign_level = sign_level,
    name_ex = paste("linear_500pts", sign_level, sep = "_")
  )
  plot_envelope_theo_2(
    f_df$f_df_MAD,
    font_size = 12,
    name_plot = paste("curves_linear_500pts_MAD_theo", sign_level, sep = "_")
  )
  plot_envelope_mean_2(
    f_df$f_df_MAD,
    font_size = 12,
    name_plot = paste("curves_linear_500pts_MAD_mean", sign_level, sep = "_")
  )
}


#Square with a hole
n <- 1000
r0 <- 0.5
b <- 1
theta <- runif(n, 0, 2 * pi)
r <- runif(
  n,
  r0,
  case_when(
    (theta < pi / 4) | (theta >= 7 * pi / 4) ~ 1 / cos(theta),
    (theta >= pi / 4) & (theta < 3 * pi / 4) ~ 1 / sin(theta),
    (theta >= 3 * pi / 4) & (theta < 5 * pi / 4) ~ -1 / cos(theta),
    (theta >= 5 * pi / 4) & (theta < 7 * pi / 4) ~ -1 / sin(theta),
  )
)
X1 <- r * cos(theta)
Y <- data.frame(Y = r * sin(theta))
X <- data.frame(X1)

plot_clouds_cuad(X, Y, "cloud_square_hole_0.5.png")
plot_clouds_R(X, Y, "R_cloud_square_hole_0.5.png")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "square_hole_0.5")
plot_envelope_theo(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_square_hole_0.5_theo"
)
plot_envelope_mean(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_square_hole_0.5_mean"
)


n <- 1000
r0 <- 0.40
b <- 1
theta <- runif(n, 0, 2 * pi)
r <- runif(
  n,
  r0,
  case_when(
    (theta < pi / 4) | (theta >= 7 * pi / 4) ~ 1 / cos(theta),
    (theta >= pi / 4) & (theta < 3 * pi / 4) ~ 1 / sin(theta),
    (theta >= 3 * pi / 4) & (theta < 5 * pi / 4) ~ -1 / cos(theta),
    (theta >= 5 * pi / 4) & (theta < 7 * pi / 4) ~ -1 / sin(theta),
  )
)
X1 <- r * cos(theta)
Y <- data.frame(Y = r * sin(theta))
X <- data.frame(X1)

plot_clouds_cuad(X, Y, "cloud_square_hole_0.4.png")
plot_clouds_R(X, Y, "R_cloud_square_hole_0.4.png")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "square_hole_0.4")
plot_envelope_theo(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_square_hole_0.4_theo"
)
plot_envelope_mean(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_square_hole_0.4_mean"
)

n <- 1000
r0 <- 0.30
b <- 1
theta <- runif(n, 0, 2 * pi)
r <- runif(
  n,
  r0,
  case_when(
    (theta < pi / 4) | (theta >= 7 * pi / 4) ~ 1 / cos(theta),
    (theta >= pi / 4) & (theta < 3 * pi / 4) ~ 1 / sin(theta),
    (theta >= 3 * pi / 4) & (theta < 5 * pi / 4) ~ -1 / cos(theta),
    (theta >= 5 * pi / 4) & (theta < 7 * pi / 4) ~ -1 / sin(theta),
  )
)
X1 <- r * cos(theta)
Y <- data.frame(Y = r * sin(theta))
X <- data.frame(X1)

plot_clouds_cuad(X, Y, "cloud_square_hole_0.3.png")
plot_clouds_R(X, Y, "R_cloud_square_hole_0.3.png")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "square_hole_0.3")
plot_envelope_theo(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_square_hole_0.3_theo"
)
plot_envelope_mean(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_square_hole_0.3_mean"
)


n <- 1000
r0 <- 0.20
b <- 1
theta <- runif(n, 0, 2 * pi)
r <- runif(
  n,
  r0,
  case_when(
    (theta < pi / 4) | (theta >= 7 * pi / 4) ~ 1 / cos(theta),
    (theta >= pi / 4) & (theta < 3 * pi / 4) ~ 1 / sin(theta),
    (theta >= 3 * pi / 4) & (theta < 5 * pi / 4) ~ -1 / cos(theta),
    (theta >= 5 * pi / 4) & (theta < 7 * pi / 4) ~ -1 / sin(theta),
  )
)
X1 <- r * cos(theta)
Y <- data.frame(Y = r * sin(theta))
X <- data.frame(X1)

plot_clouds_cuad(X, Y, "cloud_square_hole_0.2.png")
plot_clouds_R(X, Y, "R_cloud_square_hole_0.2.png")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "square_hole_0.2")
plot_envelope_theo(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_square_hole_0.2_theo"
)
plot_envelope_mean(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curves_square_hole_0.2_mean"
)


#Kaggle

#Advertising

Datos <- read.csv("advertising.csv", header = TRUE, sep = ",", dec = ".")
head(Datos)
X <- Datos[, 1:3]
Y <- Datos[, 4]
plot_clouds_real(X, Y, "cloud_advertising.png", "Sales")
plot_clouds_real_R(X, Y, "R_cloud_advertising.png", "Sales")


f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "results_advertising")
plot_envelope_theo_real(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_advertising_theo"
)
plot_envelope_mean_real(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_advertising_mean"
)

f_df <- CSR_test(X, Y, sign_level = 0.20, name_ex = "results_advertising_0.2")
plot_envelope_theo_real(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_advertising_theo_0.20"
)
plot_envelope_mean_real(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_advertising_mean_0.20"
)


#Car Price

Datos <- read.csv(
  "CarPrice_Assignment.csv",
  header = TRUE,
  sep = ",",
  dec = "."
)
head(Datos)
sum(is.na(Datos))
X <- Datos[, 11:13]
Y <- Datos[, 26]
colnames(X) = c("length", "width", "height")

plot_clouds_real_2(X, Y, "cloud_CarPrice.png", "Price")

f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "results_CarPrice")
plot_envelope_theo_real_2(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_CarPrice_theo"
)


#Cars with log

Datos <- read.csv(
  "CarPrice_Assignment.csv",
  header = TRUE,
  sep = ",",
  dec = "."
)
head(Datos)
sum(is.na(Datos))
X <- Datos[, 11:13]
Y <- log(Datos[, 26])
colnames(X) = c("length", "width", "height")

plot_clouds_real_2(X, Y, "Log_CarPrice.png", "Log(Price)")

f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "results_Log_CarPrice")
plot_envelope_theo_real_2(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_Log_CarPrice_theo"
)


#Este no tiene buenos resultados
Datos <- read.csv("Real estate.csv", header = TRUE, sep = ",", dec = ".")
head(Datos)
sum(is.na(Datos))
X <- Datos[, 3:5]
Y <- Datos[, 8]
plot_clouds_real(X, Y, "RealEstate.png", "Price unit area")
S_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "results_RealEstate")
plot_envelope2(S_df$S_df1, font_size = 12, name_plot = "curve_RealEstate")


Datos <- read_excel("Concrete_Data.xls")
colnames(Datos) <- c(
  'Cement',
  'Blast Furnace Slag',
  'Fly Ash',
  'Water',
  'Superplasticizer',
  'Coarse Aggregate',
  'Fine Aggregate',
  'Age',
  'Concrete comp strength'
)
sum(is.na(Datos))
X <- as.data.frame(Datos[, 7])
Y <- unlist(Datos[, 8])
plot_clouds_real(X, Y, "Concrete1-3.png", "Concrete Comp. Strength")


#Football

urlfile = "https://raw.githubusercontent.com/rjtavares/football-crunching/faa7d5105ba9b01a821b80dd70fb269fdf145889/datasets/james-vs-barcelona-positonal-data.csv"
Datos <- read.csv(url(urlfile))
head(Datos)
sum(is.na(Datos))

X <- as.data.frame(Datos[which(Datos$frame == 0), 3])
colnames(X) = "X"
Y <- Datos[which(Datos$frame == 0), 4]
#plot_clouds_real(X,Y, "FootballPositionalData_0_0.2.png", "Y")
f_df <- CSR_test(X, Y, sign_level = 0.20, name_ex = "results_Football_0_0.2")
plot_envelope_theo_real(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_Football_0_0.2"
)


X <- as.data.frame(Datos[which(Datos$frame == 50), 3])
colnames(X) = "X"
Y <- Datos[which(Datos$frame == 50), 4]
#plot_clouds_real(X,Y, "FootballPositionalData_50.png", "Y")
f_df <- CSR_test(X, Y, sign_level = 0.20, name_ex = "results_Football_50_0.2")
plot_envelope_theo_real(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_Football_50_0.2"
)


X <- as.data.frame(Datos[which(Datos$frame == 130), 3])
colnames(X) = "X"
Y <- Datos[which(Datos$frame == 130), 4]
#plot_clouds_real(X,Y, "FootballPositionalData_130_0.2.png", "Y")
f_df <- CSR_test(X, Y, sign_level = 0.20, name_ex = "results_Football_130_0.2")
plot_envelope_theo_real(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_Football_130_0.2"
)


X <- as.data.frame(Datos[which(Datos$frame == 180), 3])
colnames(X) = "X"
Y <- Datos[which(Datos$frame == 180), 4]
#plot_clouds_real(X,Y, "FootballPositionalData_180.png", "Y")
f_df <- CSR_test(X, Y, sign_level = 0.20, name_ex = "results_Football_180_0.2")
plot_envelope_theo_real(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_Football_180_0.2"
)


X <- as.data.frame(Datos[which(Datos$frame == 0), 3])
colnames(X) = "X"
Y <- Datos[which(Datos$frame == 0), 4]
#plot_clouds_real(X,Y, "FootballPositionalData_0.png", "Y")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "results_Football_0_0.05")
plot_envelope_theo_real(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_Football_0_0.05"
)


X <- as.data.frame(Datos[which(Datos$frame == 50), 3])
colnames(X) = "X"
Y <- Datos[which(Datos$frame == 50), 4]
#plot_clouds_real(X,Y, "FootballPositionalData_50.png", "Y")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "results_Football_50_0.05")
plot_envelope_theo_real(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_Football_50_0.05"
)


X <- as.data.frame(Datos[which(Datos$frame == 130), 3])
colnames(X) = "X"
Y <- Datos[which(Datos$frame == 130), 4]
#plot_clouds_real(X,Y, "FootballPositionalData_130.png", "Y")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "results_Football_130_0.05")
plot_envelope_theo_real(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_Football_130_0.05"
)


X <- as.data.frame(Datos[which(Datos$frame == 180), 3])
colnames(X) = "X"
Y <- Datos[which(Datos$frame == 180), 4]
#plot_clouds_real(X,Y, "FootballPositionalData_180.png", "Y")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "results_Football_180_0.05")
plot_envelope_theo_real(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_Football_180_0.05"
)


#All players and balls every 10 frames.
frames = seq(1, 4160, by = 10)
X <- as.data.frame(filter(Datos, frame %in% frames)[, 3])
colnames(X) = "X"
Y <- filter(Datos, frame %in% frames)[, 4]
#plot_clouds_real(X,Y, "FootballPositionalData_all.png", "Y")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "results_Football_all_0.05")
plot_envelope_theo_real(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_Football_all_0.05"
)

#Barza and balls every 10 frames.
frames = seq(1, 4160, by = 10)
Barza_players = seq(11, 19)
Datos_Barza = filter(Datos, player %in% Barza_players)
X <- as.data.frame(filter(Datos_Barza, frame %in% frames)[, 3])
colnames(X) = "X"
Y <- filter(Datos_Barza, frame %in% frames)[, 4]
#plot_clouds_real(X,Y, "FootballPositionalData_Barza.png", "Y")
f_df <- CSR_test(
  X,
  Y,
  sign_level = 0.05,
  name_ex = "results_Football_Barza_0.05"
)
plot_envelope_theo_real(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_Football_Barza_0.05"
)

#RM and balls every 10 frames.
frames = seq(1, 4160, by = 10)
RM_players = seq(20, 25)
Datos_RM = filter(Datos, player %in% RM_players)
X <- as.data.frame(filter(Datos_RM, frame %in% frames)[, 3])
colnames(X) = "X"
Y <- filter(Datos_RM, frame %in% frames)[, 4]
#plot_clouds_real(X,Y, "FootballPositionalData_RM.png", "Y")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "results_Football_RM_0.05")
plot_envelope_theo_real(
  f_df$f_df_MAD,
  font_size = 12,
  name_plot = "curve_Football_RM_0.05"
)


#Example of 7 points

X <- c(-3, 4, -1, 2, -3, 0, 1)
Y <- c(1, 2, 4, 5, 6, 1, 7)
Y <- data.frame(Y)
X <- data.frame(X)

plot_clouds(X, Y, "7pts.png")
f_df <- CSR_test(X, Y, sign_level = 0.05, name_ex = "7pts")

plot_f_i(f_df$f_df_MAD, font_size = 12, name_plot = "f_i_7pts_MAD")


S_df$S_df1$x
S_df$S_df1$S_obs1

TeX(r"( $\gamma^2 = \alpha^2 + \beta^2$ )")
#Fig3.1-1
n <- 1000
a <- -1
b <- 1
theta <- runif(n, 0, 2 * pi)
r <- (sqrt(runif(n))) * (0.5) + 0.5
X_1 <- r * cos(theta)
X2 <- runif(n, a, b)
Y <- data.frame(Y = r * sin(theta))
X <- data.frame(X_1)


plot_clouds_only <- function(X, Y, namefile.png) {
  df <- cbind(X, Y)
  df <- tidyr::pivot_longer(
    df,
    cols = -Y,
    names_to = "variable",
    values_to = "points"
  )
  g <- ggplot(df, aes(points, Y)) +
    ggplot2::geom_point(size = 0.75) +
    stat_poly_line() +
    stat_poly_eq() +
    cowplot::theme_cowplot(font_size = 12) +
    cowplot::background_grid(minor = "y") +
    #cowplot::panel_border()
    ggtitle("") +
    xlab("X") +
    ylab("Y")
  facet_wrap(. ~ variable, scales = "free")

  cowplot::save_plot(namefile.png, g, bg = "white")
  return(g)
}


plot_clouds_only(X, Y, "fig3.1-1_ring1000_pts.png")

#Relating \alpha with r for the f_theo_ap

a <- -1
b <- 1
r_vs_alpha = data.frame(seq(0.1, 2, 0.1))
for (n in c(30, 100, 500, 1000)) {
  X <- runif(n, a, b)
  Y <- data.frame(Y = r * sin(theta))
  X <- data.frame(X)

  mse = c()
  for (i in seq(0.1, 2, 0.1)) {
    S_df <- CSR_test(
      X,
      Y,
      sign_level = 0.05,
      name_ex = "info_alphavsr_30_pts.png",
      r = i
    )
    mse_i = mean((S_df$S_df1$teor - S_df$S_df1$mean)^2)
    mse = append(mse, mse_i)
  }
  r_vs_alpha = cbind(r_vs_alpha, mse)
}

colnames(r_vs_alpha) = c("r", "MSE_30", "MSE_100", "MSE_500", "MSE_1000")
r_vs_alpha_res <- r_vs_alpha
r_vs_alpha[, 2:5] = round(r_vs_alpha[, 2:5], 4)

r_vs_alpha[which.min(r_vs_alpha$MSE_30), ]
r_vs_alpha[which.min(r_vs_alpha$MSE_100), ]
r_vs_alpha[which.min(r_vs_alpha$MSE_500), ]
r_vs_alpha[which.min(r_vs_alpha$MSE_1000), ]


library(knitr)
kable(r_vs_alpha, "latex")
