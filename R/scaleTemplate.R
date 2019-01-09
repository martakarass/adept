



#' Compute a list of rescaled versions of pattern templates.
#'
#' Compute a list of rescaled versions of pattern templates via linear interpolation.
#'
#' @param template A list of numerical vectors. Each numeric vector represents
#' a distinct pattern template.
#' @param template.vl A numeric vector. Defines a grid of vector lengths to which each of
#' the distinct pattern templates is scaled into (linearly interpolated into).
#'
#' @return A list of lists of numeric vectors.  Each element of
#' the list is a list of pattern templates scaled according
#' to a particular scale parameter (that is, scaled to a particular vector length).
#' Number of elements in the returned list corresponds to the number
#' of unique scale values considered in the method.
#'
#' @examples
#' ## List of two distinct pattern templates
#' template <- list(sin(seq(0, 2 * pi, length.out = 100)),
#'                  cos(seq(0, 2 * pi, length.out = 100)))
#' ## Vector which defines a grid of vector lengths to which each of
#' ## the distinct pattern templates is scaled into
#' template.vl <- c(50, 100, 200)
#' ## Compute list of rescaled versions of pattern templates
#' out <- scaleTemplate(template, template.vl)
#'
#' ## Plot 1st distinct pattern template, rescaled to different vector lengths
#' par(mfrow = c(2, 1))
#' plot(out[[3]][[1]], type = "l",
#'      main = "pattern: sin([0,2*pi]); different scales considered",
#'      ylab = "pattern", xlab = "output vector index")
#' lines(out[[2]][[1]], col = "red")
#' lines(out[[1]][[1]], col = "blue")
#'
#' ## Plot 2nd distinct pattern template, rescaled to different vector lengths
#' plot(out[[3]][[2]], type = "l",
#'      main = "pattern: cos([0,2*pi]); different scales considered",
#'      ylab = "pattern", xlab = "output vector index")
#' lines(out[[2]][[2]], col = "red")
#' lines(out[[1]][[2]], col = "blue")
#'
#' @export
#'
#'
scaleTemplate <- function(template, template.vl){

  ## Check function arguments for correctness
  if (!(is.list(template))) stop("template must be a list.")
  if (!(all(sapply(template, is.vector)))) stop("template must be a list of numeric vectors")
  if (!(all(sapply(template, is.numeric)))) stop("template must be a list of numeric vectors")
  if (!(is.vector(template.vl))) stop("template.vl must be a numeric vector")
  if (!(is.numeric(template.vl))) stop("template.vl must be a numeric vector")

  ## Outer lapply: iterate over vector lengths;
  ## in each lapply iteration, empirical pattern(s) are linearly interpolated
  ## to a particular vector length.
  lapply(template.vl, function(vl){

    ## Inner lapply: iterate over, possibly, multiple patterns;
    ## each lapply iteration returns a vector whose each element corresponds
    ## to distinct empirical pattern linearly interpolated to
    ## a particular vector length.
    lapply(template, function(template.i){
      out.i <- approx(seq(0, 1, length.out = length(template.i)),
                      template.i,
                      seq(0, 1, length.out = vl))$y
      as.vector(scale(out.i))
    })
  })
}

