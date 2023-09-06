
if (FALSE) {
  data = x.smoothed
  # cut_x = function(data, x.cut.seq, x.cut.vl) {
  # x.cut.vl = unique(diff(x.cut.seq))

  make_index_mat = function(data, x.cut.vl, x.cut.seq, x.cut.margin) {
    n = length(data)
    nr = length(x.cut.seq)
    ind_mat = matrix(1:(nr * x.cut.vl), ncol = x.cut.vl, nrow = nr,
                     byrow = TRUE)
    addon_mat = matrix(rep(1:(x.cut.margin + 1), times = nr),
                       nrow = nr, ncol = x.cut.margin + 1,
                       byrow = TRUE)
    addon_mat = addon_mat + ind_mat[, ncol(ind_mat)]
    ind_mat = cbind(ind_mat, addon_mat)
    ind_mat[ind_mat > n] = NA
    ind_mat
  }
  ind_mat = make_index_mat(n, x.cut.vl, x.cut.seq, x.cut.margin)
  x_mat = array(data[ind_mat], dim = dim(ind_mat))
  nc = ncol(x_mat)
  temp = template.scaled[[1]][[1]]
  make_shift_matrix = function(temp, nc) {
    template_length = length(temp)
    temp = c(temp, rep(0, nc - length(temp)))
    # system.time({
    #   x = toeplitz(temp)
    #   x[upper.tri(x)] = 0
    #   x
    # })
    system.time({
      tt = c(temp[1], rep(0, nc - 1))
      shift_mat = pracma::Toeplitz(a = temp, b = tt)
    })
    shift_mat = shift_mat[, 1:(nc - template_length + 1)]
    shift_mat = Matrix::Matrix(shift_mat, sparse = TRUE)
  }
  xy = x_mat %*% shift_mat

  # ## If we cannot fit the longest pattern, return NULL
  # if (length(idx.i) <= max(template.vl)) return(NULL)
  # ## Compute similarity matrix
  # similarity.mat.i <- similarityMatrix(x = x.smoothed[idx.i],
  #                                      template.scaled = template.scaled,
  #                                      similarity.measure = similarity.measure)
  #
  #   for (i in )
  #     vec_result = adept:::similarityMatrix_vectorized(
  #       x.cut.seq = x.cut.seq,
  #       x.smoothed = x.smoothed,
  #       template.scaled = template.scaled,
  #       similarity.measure = similarity.measure,
  #       x.cut.vl = x.cut.vl,
  #       template.vl = template.vl,
  #       x.cut.margin = x.cut.margin)
  # })

}
