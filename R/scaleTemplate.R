



#' Pattern Templates Scaling
#'
#' Compute a list of scaled versions of pattern templates via linear interpolation.
#'
#' @param template A list of numeric vectors. Each vector represents
#' a distinct template.
#' @param template.vl A numeric vector. A grid of vector lengths that each
#' element of \code{template} is to be linearly interpolated into.
#'
#' @return A list of lists of numeric vectors.
#' Each element of the
#' returned list is a list of templates scaled according
#' to a particular vector length. The number of elements in the returned
#' list equals the length of \code{template.vl}.
#'
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
#' @importFrom stats approx
#'
#'
scaleTemplate <- function(template, template.vl){

  ## Check function arguments for correctness
  if (!(is.list(template))) stop("template must be a list.")
  if (!(all(sapply(template, is.numeric)) & all(sapply(template, is.atomic)))) stop("template must be a list of numeric (atomic) vectors.")
  if (!(is.numeric(template.vl) & is.atomic(template.vl))) stop("template.vl must be a numeric (atomic) vector.")

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

