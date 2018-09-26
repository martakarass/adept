#' Silly Printer Function
#'
#' @param r What you want in the second column
#' @param x What you want in the first column
#'
#' @return tibble
#' @export
#'
#' @importFrom tibble data_frame
#' @importFrom utils head
#' @importFrom dplyr select
#' @examples
#' printer(x = rnorm(5), r = rnorm(5))
printer = function(r, x) {
  x = data_frame(x = x, r = r)
  print(head(x))
  print(paste0("x = ", x))
  return(x)
}
