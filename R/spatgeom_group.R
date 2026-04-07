#' @title Geometric Spatial Point Pattern Analysis by Group
#'
#' @description Apply \code{\link{spatgeom}} independently to each level of a
#'   grouping variable. The result is an object of class \code{spatgeom_group}
#'   that bundles the per-group \code{spatgeom} objects and supports
#'   \code{print} and \code{\link{plot_curve}} methods for side-by-side
#'   comparison.
#'
#' @param x a numeric matrix or data frame of covariates (passed as the
#'   \code{x} argument of \code{\link{spatgeom}}). Rows correspond to
#'   observations.
#' @param by a vector of group labels with the same length as \code{nrow(x)}.
#'   Can be a character vector, factor, or any vector that can be coerced to
#'   character. Each unique value defines one group.
#' @param ... additional arguments forwarded verbatim to
#'   \code{\link{spatgeom}} for every group. This includes \code{y},
#'   \code{nalphas}, \code{domain_type}, \code{envelope},
#'   \code{hypothesis_testing}, \code{method}, \code{use_umap},
#'   \code{umap_args}, \code{mc_cores}, etc. The subset of \code{x} rows
#'   belonging to each group is passed automatically; do \emph{not} pass
#'   \code{x} again inside \code{...}.
#'
#' @details
#' When \code{y} is supplied via \code{...}, the corresponding rows of \code{y}
#' are subsetted to match each group. When \code{y} is absent, each group's
#' subset of \code{x} is treated as a 2-D point cloud (or projected with UMAP
#' if \code{use_umap = TRUE}).
#'
#' Groups are processed in the order returned by \code{unique(by)}, which
#' preserves the order of first appearance in \code{by}.
#'
#' @return An object of class \code{spatgeom_group}, which is a list with:
#' \describe{
#'   \item{\strong{call}}{The matched function call.}
#'   \item{\strong{groups}}{Character vector of unique group labels in
#'     appearance order.}
#'   \item{\strong{results}}{Named list of \code{spatgeom} objects, one per
#'     group. Names match \code{groups}.}
#' }
#'
#' @seealso \code{\link{spatgeom}}, \code{\link{plot_curve}},
#'   \code{\link{print.spatgeom_group}}
#'
#' @examples
#' set.seed(1)
#' xy <- donut_data(n = 60, a = -1, b = 1, theta = 2 * pi)
#' grp <- sample(c("A", "B", "C"), nrow(xy), replace = TRUE)
#'
#' sg <- spatgeom_group(x = xy[, -1], by = grp, y = xy[, 1])
#' print(sg)
#' plot_curve(sg)
#'
#' \donttest{
#' # With hypothesis testing per group (slow — runs Monte Carlo per group):
#' sg_ht <- spatgeom_group(
#'   x = xy[, -1], by = grp, y = xy[, 1],
#'   hypothesis_testing = TRUE, method = "MAD"
#' )
#' plot_curve(sg_ht, type = "curve")
#' }
#'
#' @export

spatgeom_group <- function(x, by, ...) {
  x <- as.data.frame(x)
  by <- as.character(by)

  if (length(by) != nrow(x)) {
    stop("'by' must have the same length as nrow(x).")
  }

  dots <- list(...)

  # If 'y' was supplied, we need to subset it per group as well.
  y_supplied <- "y" %in% names(dots)
  if (y_supplied) {
    y_full <- dots[["y"]]
    if (length(y_full) != nrow(x)) {
      stop("'y' must have the same length as nrow(x).")
    }
    dots[["y"]] <- NULL # will be reattached per-group below
  }

  groups <- unique(by)

  results <- lapply(groups, function(g) {
    idx <- which(by == g)
    call_args <- c(list(x = x[idx, , drop = FALSE]), dots)
    if (y_supplied) {
      call_args[["y"]] <- y_full[idx]
    }
    do.call(spatgeom, call_args)
  })
  names(results) <- groups

  structure(
    list(call = match.call(), groups = groups, results = results),
    class = "spatgeom_group"
  )
}


#' @title Print a \code{spatgeom_group} object
#'
#' @description Print method for objects of class \code{spatgeom_group}.
#'   Displays a summary for each group in turn.
#'
#' @param x an object of class \code{spatgeom_group}.
#' @param ... further arguments passed to \code{\link{print.spatgeom}}.
#'
#' @return \code{x} invisibly.
#'
#' @seealso \code{\link{spatgeom_group}}, \code{\link{print.spatgeom}}
#'
#' @examples
#' set.seed(1)
#' xy <- donut_data(n = 60, a = -1, b = 1, theta = 2 * pi)
#' grp <- sample(c("A", "B", "C"), nrow(xy), replace = TRUE)
#' sg <- spatgeom_group(x = xy[, -1], by = grp, y = xy[, 1])
#' print(sg)
#'
#' @export

print.spatgeom_group <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n")
  cat("Groups:", paste(x$groups, collapse = ", "), "\n")
  for (g in x$groups) {
    cat("\n\u2500\u2500 Group:", g, "\u2500\u2500\n")
    print(x$results[[g]], ...)
  }
  invisible(x)
}


#' @title Plot a \code{spatgeom_group} object
#'
#' @description Plot method for objects of class \code{spatgeom_group}.
#'   Produces a side-by-side grid of survival curve (or derivative) panels,
#'   one column per group.
#'
#' @param x an object of class \code{spatgeom_group}.
#' @param type a string: either \code{"curve"} (default) or \code{"deriv"}.
#'   Passed to \code{\link{plot_curve.spatgeom}} for each group.
#' @param font_size an integer controlling the font size. Default \code{12}.
#' @param ... further arguments passed to \code{\link{plot_curve.spatgeom}}.
#'
#' @return A \code{cowplot} grid object (produced by
#'   \code{cowplot::plot_grid}).
#'
#' @seealso \code{\link{spatgeom_group}}, \code{\link{plot_curve}}
#'
#' @examples
#' set.seed(1)
#' xy <- donut_data(n = 60, a = -1, b = 1, theta = 2 * pi)
#' grp <- sample(c("A", "B", "C"), nrow(xy), replace = TRUE)
#' sg <- spatgeom_group(x = xy[, -1], by = grp, y = xy[, 1])
#' plot_curve(sg)
#' plot_curve(sg, type = "deriv")
#'
#' @export

plot_curve.spatgeom_group <- function(x, type = "curve", font_size = 12, ...) {
  plots <- lapply(x$groups, function(g) {
    plot_curve.spatgeom(
      x$results[[g]],
      type = type,
      font_size = font_size,
      ...
    ) +
      ggplot2::labs(title = g)
  })
  cowplot::plot_grid(plotlist = plots, nrow = 1)
}
