


#' Title
#'
#' @param x
#' @param template.scaled
#' @param similarity.measure
#' @param run.parallel
#' @param run.parallel.ncores
#'
#' @return
#'
#' @import parallel
#' @import runstats
#'
#' @export
#'
#' @examples
#'
adeptSimilarity <- function(x,
                            template.scaled,
                            similarity.measure,
                            run.parallel = FALSE,
                            run.parallel.ncores = NULL){

  runstat.func <- switch(similarity.measure,
                         "cov" = RunningCov,
                         "cor" = RunningCor)

  if (run.parallel) {

    ## Parallel execution
    if (is.null(run.parallel.ncores)) run.parallel.ncores <- detectCores() - 1
    cl <- makeCluster(run.parallel.ncores)
    clusterExport(cl,
                  c("x", "template.scaled", "runstat.func"),
                  envir = environment())
    ## parLapply: apply over different scale parameters within `template.scaled`
    similarity.list <- parLapply(cl, template.scaled, function(template.scaled.i){
      ## lapply: for particular scale parameter, apply over possibly multiple templates
      runstat.func.out0 <- lapply(template.scaled.i, function(template.scaled.ik){
        do.call(runstat.func, list(x = x, y = template.scaled.ik))
      })
      ## Take best result over possibly multiple templates
      do.call(pmax, runstat.func.out0)
    })
    stopCluster(cl)

  } else {

    ## Non-parallel execution
    ## lapply: apply over different scale parameters within `template.scaled`
    similarity.list <- lapply(template.scaled, function(template.scaled.i){
      ## lapply: for particular scale parameter, apply over possibly multiple templates
      runstat.func.out0 <- lapply(template.scaled.i, function(template.scaled.ik){
        do.call(runstat.func, list(x = x, y = template.scaled.ik))
      })
      ## Take best result over possibly multiple templates
      do.call(pmax, runstat.func.out0)
    })

  }

  ## List of vectors to matrix
  similarity.mat <- do.call(rbind, similarity.list)
  return(similarity.mat)

}
