



#' Compute a list of rescaled versions of empirical pattern(s)
#'
#' Compute a list of rescaled versions of empirical pattern(s) via linear interpolation.
#'
#' @param template A list of numerical vector(s) representing
#' empirical pattern(s).
#' @param template.vl A numeric vector of vector lengths to which each of
#' the empirical pattern(s) will be scaled to.
#'
#' @return A list of scaled versions of empirical pattern(s).
#' Number of elements in the list corresponds to the number of unique pattern scale
#' values considered. Each element of the list is a list
#' itself that contains a numeric vector(s), where each vector represents
#' one of possibly multiple distinct empirical patterns considered, scaled
#' to have the same scale value (same vector length).
#'
#' @examples
#' template <- list(sin(seq(0, 2 * pi, length.out = 100)),
#'                  cos(seq(0, 2 * pi, length.out = 100)))
#' template.vl <- c(50, 100, 200)
#' out <- scaleTemplate(template, template.vl)
#' par(mfrow = c(2, 1))
#' plot(out[[3]][[1]], type = "l",
#'      main = "pattern: sin([0,2*pi]); different scales considered",
#'      ylab = "pattern", xlab = "output vector index")
#' lines(out[[2]][[1]], col = "red")
#' lines(out[[1]][[1]], col = "blue")
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
