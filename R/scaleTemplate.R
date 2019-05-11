



#' Templates Scaling
#'
#' Compute a list of scaled templates via linear interpolation.
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
#' ## Construct a list of two templates
#' template <- list(sin(seq(0, 2 * pi, length.out = 100)),
#'                  cos(seq(0, 2 * pi, length.out = 100)))
#' ## A grid of vector lengths to which each of templates is scaled into
#' template.vl <- c(50, 100, 200)
#' ## Compute list of rescaled templates
#' out <- scaleTemplate(template, template.vl)
#'
#' ## Plot 1st template after rescaling to three values of vector length
#' par(mfrow = c(2, 1), cex = 0.7)
#' plot(out[[3]][[1]], type = "l",
#'      main = "Pattern: sin([0, 2 * pi]) rescaled according to different scales",
#'      ylab = "Pattern", xlab = "Index")
#' lines(out[[2]][[1]], col = "red")
#' lines(out[[1]][[1]], col = "blue")
#'
#' ## Plot 2nd template after rescaling to three values of vector length
#' plot(out[[3]][[2]], type = "l",
#'      main = "Pattern: cos([0, 2 * pi]) rescaled according to different scales",
#'      ylab = "Pattern", xlab = "Index")
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

