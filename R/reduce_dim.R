#' @title Dimensionality Reduction for Spatial Geometry Analysis
#'
#' @description Reduce an \eqn{n \times p} matrix to \eqn{n \times d}
#'   (\code{d = n_components}, typically 2) using a choice of static
#'   embedding methods. The result can be passed directly to
#'   \code{\link{spatgeom}} or \code{\link{spatgeom_group}} when the original
#'   feature matrix has more than 2 columns.
#'
#' @param x a numeric matrix or data frame with \eqn{n} rows and \eqn{p}
#'   columns. For \code{method = "pca"}, \eqn{p \ge d} is required; for
#'   \code{"umap"} and \code{"tsne"}, \code{n_components} may exceed
#'   \code{ncol(x)}.
#' @param method character; the reduction method. One of:
#'   \describe{
#'     \item{\code{"pca"}}{Principal Component Analysis via
#'       \code{\link[stats]{prcomp}}. No extra package required.}
#'     \item{\code{"umap"}}{Uniform Manifold Approximation and Projection via
#'       \code{uwot::umap}. Requires the \pkg{uwot} package.}
#'     \item{\code{"tsne"}}{t-distributed Stochastic Neighbour Embedding via
#'       \code{Rtsne::Rtsne}. Requires the \pkg{Rtsne} package.}
#'   }
#' @param n_components integer; number of output dimensions. Default \code{2L}.
#' @param ... additional arguments forwarded to the underlying reduction
#'   function (\code{prcomp}, \code{uwot::umap}, or \code{Rtsne::Rtsne}).
#'
#' @return A numeric matrix with \code{nrow(x)} rows and \code{n_components}
#'   columns. Column names are \code{Dim_1}, \code{Dim_2}, etc.
#'
#' @seealso \code{\link{spatgeom}}, \code{\link{spatgeom_group}}
#'
#' @examples
#' set.seed(1)
#' xy <- donut_data(n = 50, a = -1, b = 1, theta = 2 * pi)
#'
#' # PCA — no extra packages needed
#' emb <- reduce_dim(xy[, -1], method = "pca")
#' dim(emb) # 50 x 2
#'
#' \donttest{
#' # UMAP — requires uwot
#' if (requireNamespace("uwot", quietly = TRUE)) {
#'   emb_umap <- reduce_dim(xy[, -1], method = "umap")
#'   dim(emb_umap)
#' }
#'
#' # t-SNE — requires Rtsne; perplexity must be < n/3
#' if (requireNamespace("Rtsne", quietly = TRUE)) {
#'   emb_tsne <- reduce_dim(xy[, -1], method = "tsne", perplexity = 10)
#'   dim(emb_tsne)
#' }
#' }
#'
#' @export

reduce_dim <- function(
  x,
  method = c("pca", "umap", "tsne"),
  n_components = 2L,
  ...
) {
  method <- match.arg(method)
  x <- as.matrix(x)

  if (!is.numeric(x)) {
    stop("'x' must be a numeric matrix or data frame.")
  }
  if (nrow(x) < 2) {
    stop("'x' must have at least 2 rows.")
  }
  if (n_components < 1L) {
    stop("'n_components' must be at least 1.")
  }
  if (method == "pca") {
    max_components <- min(ncol(x), nrow(x) - 1L)
    if (n_components > max_components) {
      stop(
        "'n_components' = ",
        n_components,
        " exceeds the maximum supported for method = 'pca' with this input (",
        max_components,
        ")."
      )
    }
  }

  result <- switch(
    method,
    pca = {
      fit <- stats::prcomp(x, ...)
      fit$x[, seq_len(n_components), drop = FALSE]
    },
    umap = {
      if (!requireNamespace("uwot", quietly = TRUE)) {
        stop(
          "Package 'uwot' is required for method = 'umap'. ",
          "Install it with: install.packages(\"uwot\")"
        )
      }
      uwot::umap(x, n_components = n_components, ...)
    },
    tsne = {
      if (!requireNamespace("Rtsne", quietly = TRUE)) {
        stop(
          "Package 'Rtsne' is required for method = 'tsne'. ",
          "Install it with: install.packages(\"Rtsne\")"
        )
      }
      Rtsne::Rtsne(x, dims = n_components, ...)$Y
    }
  )

  result <- as.matrix(result)
  colnames(result) <- paste0("Dim_", seq_len(n_components))
  result
}
