<!-- badges: start -->
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](https://opensource.org/license/mit/)
[![CRAN](http://www.r-pkg.org/badges/version/spatgeom)](https://cran.r-project.org/package=spatgeom)
[![cran checks](https://badges.cranchecks.info/worst/spatgeom.svg)](https://cran.r-project.org/web/checks/check_results_spatgeom.html)
[![Downloads](http://cranlogs.r-pkg.org/badges/spatgeom?color=brightgreen)](http://www.r-pkg.org/pkg/spatgeom)
<!-- badges: end -->

# Geometric goodness of fit

## Description
The `spatgeom` package provides goodness-of-fit tests for spatial point process models. 
The `spatgeom` package provides the implementation to perform the geometric spatial point analysis developed in Hernández & Solís (2022) <doi:10.1007/s00180-022-01244-1>. It estimates the geometric goodness-of-fit index for a set of variables against a response one based on the 'sf' package. The package has methods to print and plot the results.

## Installation

The package is currently available on CRAN and can be installed using:

```r
install.packages("spatgeom")
```

and loaded using the following command:

```r
library(spatgeom)
```

## Usage

### Estimating the geometric goodness-of-fit index

The package has one main function `spatgeom` that takes as input a set of variables and a response one. 

For example, this is how you can use the package to estimate the geometric goodness-of-fit index with a set of variables in a donut shape against a response one:

```r
 xy <- donut_data(n = 100, a = -1, b = 1, theta = 2 * pi)
 estimation <- spatgeom(y = xy[, 1], x = xy[, -1])

```


The package allow to estimate the envelope of the curve using the parameter `envelope=TRUE`:

```r
 estimation <- spatgeom(y = xy[, 1], x = xy[, -1], envelope = TRUE)
```

In the example, the `estimation` object is a list of class `spatgeom` with the following elements:

- **call:** The function call.
- **x**: x input.
- **y**: y output.
- **results**: A list of size `ncol(x)` corresponding to each column of x. Each element of the list has:
    - **triangles**: a data frame of class sfc (see [`sf::st_sf()`])with columns geometry, segments, max_length and alpha. The data.frame contains the whole Delanauy triangulation for the corresponding column of x and y. The segments column are the segments of each individual triangle and max_length is the maximum length of them.
    - **geom_indices**: a data frame with columns alpha and geom_corr. The alpha column is a numeric vector of size nalphas from the minimum to the maximum distance between points estimated in the data. The geom_corr column is the value 1 - (alpha shape Area)/(containing box Area).
    - **intensity**: the intensity estimated for the corresponding column of x and y.
    - **mean_n**: the mean number of points in the point process.
    - **envelope_data**: a data frame in tidy format with 40 runs of a CSR process, if envelope=TRUE, The CSR is created by generating n uniform points in the plane, where n is drawn from Poisson distribution with parameter mean_n.




### Printing and plotting the results

The package has a print method for objects of class `spatgeom`:

```r
print(estimation)
```

and a plot method for objects of class `spatgeom`. The plot method has two options: `curve` and `deriv`. The `curve` option plots the curve of of the goodness-of-fit index. The `deriv` option plots the numerical derivative.

```r
plot_curve(estimation, type = "curve")
plot_curve(estimation, type = "deriv")
```

## Citation 

If you use the `spatgeom` package in your work, please cite the following paper:

- Hernández, D. and Solís, F. (2022). Geometric goodness-of-fit tests for spatial point process models. *Statistics and Computing*, 32(1), 1-19. <doi:10.1007/s00180-022-01244-1>

- Solís, M., Hernández, A., & Pasquier, C. (2023). spatgeom: Geometric Spatial Point Analysis (0.3.0). https://cran.r-project.org/web/packages/spatgeom/index.html


## License

This package is free and open source software, licensed under MIT License. See the [LICENSE](LICENSE.md) file for details.
