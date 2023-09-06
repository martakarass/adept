
run = FALSE
if (run) {
  library(dplyr)
  library(adept)
  # muschellij2/adept@nhanes
  library(adeptdata)
  all_wrist_templates = adeptdata::stride_template$left_wrist
  template_list = do.call(rbind, all_wrist_templates)
  template_list = apply(template_list, 1, identity, simplify = FALSE)

  template = template_list
  # devtools::load_all()
  options(digits.secs = 3)
  sample_rate = 10L

  # read in some test data
  data = readr::read_csv(
    "~/Dropbox/Projects/nhanes_80hz/data/csv/pax_h/73557.csv.gz",
    # "https://github.com/martakarass/adept/files/12423699/test_data.csv.gz",
    n_max = 1e6
  )
  xyz = data %>% select(all_of(c("X", "Y", "Z"))) %>% as.matrix()
  x.fs = 80L
  xyz <- as.matrix(xyz)
  # walk_out = segmentWalking(xyz,
  #                           xyz.fs = x.fs,
  #                           template = template)
  xyzptr <- as.data.frame(cbind(xyz, pracma::cart2sph(xyz)))
  x <- xyzptr[, 6]
  rm(xyz)

  pattern.dur.seq = seq(0.5, 4, length.out = 30)
  similarity.measure = "cor"
  similarity.measure.thresh = -2
  x.adept.ma.W = 0.2
  finetune = "maxima"
  finetune.maxima.ma.W = NULL
  finetune.maxima.nbh.W = 0.6
  run.parallel = FALSE
  run.parallel.cores = 1L
  x.cut = TRUE
  x.cut.vl = 6000L


  template.vl <- pattern.dur.seq * x.fs
  template.vl <- sort(unique(round(template.vl)))
  template.vl.max <- max(template.vl)
  template.vl.min <- min(template.vl)
  if (!is.list(template)) template <- list(template)
  template.scaled <- scaleTemplate(template, template.vl)


  x.smoothed <- adept:::get.x.smoothed(x = x, W = x.adept.ma.W,
                                       x.fs = x.fs)

  finetune.maxima.x <- x
  finetune.maxima.nbh.vl <- round(finetune.maxima.nbh.W *
                                    x.fs)
  if (!x.cut) x.cut.vl <- length(x)
  x.cut.margin <- template.vl.max - 1
  x.cut.seq <- seq(1, to = length(x), by = x.cut.vl)
  template.idx.mat.i <- NULL

  data = x.smoothed
  data = round(data, 5)
  # cut_x = function(data, x.cut.seq, x.cut.vl) {
  # x.cut.vl = unique(diff(x.cut.seq))


  make_index_mat = function(data, x.cut.vl, x.cut.seq, x.cut.margin) {
    n = length(data)
    nr = length(x.cut.seq)
    ind_mat = matrix(1:(nr * x.cut.vl), ncol = x.cut.vl, nrow = nr,
                     byrow = TRUE)
    addon_mat = matrix(rep(1:(x.cut.margin + 1), times = nr),
                       nrow = nr,
                       ncol = x.cut.margin + 1,
                       byrow = TRUE)
    addon_mat = addon_mat + ind_mat[, ncol(ind_mat)]
    ind_mat = cbind(ind_mat, addon_mat)
    ind_mat[ind_mat > n] = NA_integer_
    ind_mat
  }
  ind_mat = make_index_mat(data, x.cut.vl, x.cut.seq, x.cut.margin)
  x_mat = array(data[ind_mat], dim = dim(ind_mat))
  nc = ncol(x_mat)
  not_na_x = !is.na(x_mat)
  x_mat[!not_na_x] = 0
  not_na_x = Matrix::Matrix(not_na_x, sparse = TRUE)
  x_mat = Matrix::Matrix(x_mat, sparse = TRUE)

  make_shift_ones = function(template_length, nc) {
    temp = rep(TRUE, template_length)
    temp = c(temp, rep(FALSE, nc - length(temp)))
    system.time({
      tt = c(temp[1], rep(FALSE, nc - template_length))
      shift_mat = pracma::Toeplitz(a = temp, b = tt)
    })
    # shift_mat = shift_mat[, 1:(nc - template_length + 1)]
    shift_mat = Matrix::Matrix(shift_mat, sparse = TRUE)
  }

  template_list = template.scaled[[1]]
  result = pbapply::pblapply(template.scaled, function(template_list) {
    lengths = sapply(template_list, length)
    stopifnot(all(lengths == lengths[1]))
    template_length = lengths[1]

    one_mat = make_shift_ones(template_length, nc)
    n_mat = (not_na_x) %*% one_mat
    n_mat[n_mat <= 1L] = NA_integer_
    # rm(not_na_x)

    if (similarity.measure == "cor") {
      sum_x2 = (x_mat ^ 2) %*% one_mat
      sum_x = x_mat %*% one_mat
      # 1/(n-1) (Σ(x - 𝑥̄)^2
      # 1/(n-1) (Σx^2 - n𝑥̄^2)
      # 1/(n-1) (Σx^2 - n (Σx/n)^2)
      # 1/(n-1) (Σx^2 - (Σx)^2/n)
      # denominator = sqrt(1/(n_mat - 1) * (sum_x2 - (sum_x ^ 2)/n_mat))
      denominator = sqrt(sum_x2 - (sum_x ^ 2)/n_mat)
      # { 1/(n-1) Σ(x_i y_i) } / √{(1/n-1) SS_x}
      # {  Σ(x_i y_i) } / √(n-1) √{SS_x}
      rm(sum_x)
      rm(sum_x2)
      denominator = denominator * sqrt(n_mat - 1L)
    } else {
      denominator = n_mat - 1L
    }
    rm(one_mat)
    # xbar = xbar / n

    make_shift_matrix = function(temp, nc) {
      stopifnot(abs(mean(temp)) <= 1e-5,
                abs(sd(temp) - 1) <= 1e-5)
      template_length = length(temp)
      last_value = temp[length(temp)]
      temp = c(temp, rep(0, nc - template_length))
      # system.time({
      # tt = c(temp[1], rep(0, nc - 1))
      #   x = toeplitz(temp)
      #   x[upper.tri(x)] = 0
      #   x
      # })
      system.time({
        tt = c(temp[1], rep(0, nc - template_length))
        shift_mat = pracma::Toeplitz(a = temp, b = tt)
      })
      stopifnot(shift_mat[nrow(shift_mat), ncol(shift_mat)] == last_value)
      # shift_mat = shift_mat[, 1:(nc - template_length + 1)]
      shift_mat = Matrix::Matrix(shift_mat, sparse = TRUE)
    }

    itemp = 1
    # result = pbapply::pblapply(template_list, function(temp) {
    res = lapply(template_list, function(temp) {
      shift_mat = make_shift_matrix(temp, nc)
      sum_xy = x_mat %*% shift_mat

      measure = sum_xy / denominator
      measure
    })
    # res$na.rm = TRUE
    # do.call(pmax, res)
    res
  })


  # out.list <- pbapply::pblapply(x.cut.seq, function(i){
  #   ## Define current x part indices
  #   idx.i <- i : min((i + x.cut.vl + x.cut.margin), length(x))
  #   ## If we cannot fit the longest pattern, return NULL
  #   if (length(idx.i) <= max(template.vl)) return(NULL)
  #   ## Compute similarity matrix
  #   similarity.mat.i <- similarityMatrix(x = x.smoothed[idx.i],
  #                                        template.scaled = template.scaled,
  #                                        similarity.measure = similarity.measure)
  # })
}
