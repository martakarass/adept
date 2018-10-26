

#' @param template
#' @param template.vl
#'
#' @return
#'
#' @examples
#'
scaleTemplate <- function(template, template.vl){

  ## Iterate over template vector lengths
  lapply(template.vl, function(vl){

    ## Iterate over distinct templates
    lapply(template, function(template.i){
      out.i <- approx(seq(0, 1, length.out = length(template.i)),
                      template.i,
                      seq(0, 1, length.out = vl))$y
      as.vector(scale(out.i))
    })
  })
}
