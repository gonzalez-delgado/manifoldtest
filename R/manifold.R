#' Define a manifold for use in the package
#'
#' Bundles everything the package needs to know about a manifold into a
#' single object. \code{\link{manifold.dcov.test}} and
#' \code{\link{choose_base_point}} only interact with the manifold through
#' this object, so supporting a new manifold amounts to writing these
#' functions (conventionally in a script \code{R/manifold-<name>.R}, with names
#' ending in \code{_<name>}) and registering the resulting constructor in
#' \code{\link{get_manifold}}. See \code{\link{manifold_torus}} for an
#' example.
#'
#' @param name Character string identifying the manifold, e.g.
#'   \code{"torus"}.
#' @param representation Character string describing how points of the
#'   manifold are represented as matrix rows (e.g. \code{"[0,1)^d"}); used
#'   in error messages.
#' @param contains Function \code{f(X)} returning \code{TRUE} if every row
#'   of the matrix \code{X} is a point of the manifold in the expected
#'   representation, \code{FALSE} otherwise.
#' @param log_map Function \code{f(X, p)} mapping each row of \code{X} to
#'   the tangent space at base point \code{p} via the logarithmic (inverse
#'   exponential) map. Must return a numeric matrix with one row per row
#'   of \code{X}.
#' @param dist Function \code{f(X, p)} returning the numeric vector of
#'   geodesic distances between each row of \code{X} and the point
#'   \code{p}.
#' @param is_valid Function \code{f(X, p)} returning \code{TRUE} if no row
#'   of \code{X} lies on the cut locus of \code{p} (where the log map is
#'   undefined), \code{FALSE} otherwise. It should stop with an error if
#'   \code{p} is malformed (e.g. wrong length).
#' @param base_point Named list of base point selection methods. Each
#'   element is a function \code{f(X, n_start, maxit, seed, ...)} returning
#'   a list with elements \code{p} (the chosen base point) and
#'   \code{convergence} (an optimizer convergence code, or \code{NA} for
#'   non-iterative methods); arguments it does not need are absorbed by
#'   \code{...}. The first element is the default method. A data-centered
#'   method can be built from \code{dist} with the generic
#'   \code{\link{optimize_base_point}}.
#'
#' @return An object of class \code{"manifold"}: a list holding the
#'   arguments above.
#'
#' @examples
#' m <- manifold_torus()
#' m$name
#' names(m$base_point)
#'
#' @export
new_manifold <- function(name, representation, contains, log_map, dist,
                         is_valid, base_point) {

    stopifnot(is.character(name), length(name) == 1,
              is.character(representation), length(representation) == 1,
              is.function(contains), is.function(log_map),
              is.function(dist), is.function(is_valid),
              is.list(base_point), length(base_point) > 0,
              !is.null(names(base_point)), all(names(base_point) != ""),
              all(vapply(base_point, is.function, logical(1))))

    structure(list(name = name,
                   representation = representation,
                   contains = contains,
                   log_map = log_map,
                   dist = dist,
                   is_valid = is_valid,
                   base_point = base_point),
              class = "manifold")
}

# Constructors of the manifolds supported by name. To add a manifold,
# write its functions and constructor (see new_manifold) in
# R/manifold-<name>.R and add one entry here. 
manifold_registry <- list(
    torus = function() manifold_torus()
)

#' Retrieve a manifold by name
#'
#' Resolves the \code{manifold} argument of \code{\link{manifold.dcov.test}}
#' and \code{\link{choose_base_point}}.
#'
#' @param manifold Either the name of a supported manifold (currently only
#'   \code{"torus"}), or an object created with \code{\link{new_manifold}},
#'   which is returned unchanged.
#'
#' @return An object of class \code{"manifold"}.
#'
#' @examples
#' get_manifold("torus")$name
#'
#' @export
get_manifold <- function(manifold) {

    if (inherits(manifold, "manifold")) return(manifold)

    if (!is.character(manifold) || length(manifold) != 1) {
        stop("manifold must be a single character string or an object created with new_manifold().")
    }
    if (!manifold %in% names(manifold_registry)) {
        stop(sprintf("Unknown manifold '%s'. Available manifolds: %s.",
                     manifold, paste(names(manifold_registry), collapse = ", ")))
    }
    manifold_registry[[manifold]]()
}

#' @export
print.manifold <- function(x, ...) {

    cat(sprintf("<manifold: %s>\n", x$name))
    cat(sprintf("  representation:     %s\n", x$representation))
    cat(sprintf("  base point methods: %s\n", paste(names(x$base_point), collapse = ", ")))
    invisible(x)
}
