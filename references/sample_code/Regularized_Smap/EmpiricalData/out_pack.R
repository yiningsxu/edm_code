predict.with.pack <- function(C, C0, X){
  c0 = C0[nrow(C0)-1, ]
  J = C[[length(C)-1]]
  return(c0 + J%*%X)
}

out_of_sample_sequence_pack <- function(cl, Jacobiano, th, ts_training, num_points, ...){
  out_of_samp = c()
  coeff = list()
  ### Take the last point in the training set
  new_point = ts_training[nrow(ts_training), ]
  for(j in 1:num_points){
    ### Predict the first point in the training set and then allthe others
    new_point = predict.with.pack(Jacobiano$J, Jacobiano$c0, new_point)
    out_of_samp = rbind(out_of_samp, t(new_point))
    ts_training = Add_to_TS(ts_training, t(new_point))
    Jacobiano = update_Jacobian_TS_pack(cl, ts_training, TargetList, Embedding, th)
    coeff[[j]] = Jacobiano$J[[length(Jacobiano$J)]]
  }
  return(list(out_of_samp = out_of_samp, coeff = coeff))
}
Jacobian_pack_update <- function(X, TargetList, Embedding){
  J = c0 = list()
  th = c()
  n_ = 1
  for(df in TargetList){
    ########## Now compute the optimum coefficients
    L = EDM_package(X, df)
    J[[n_]] = L$coefficients
    c0[[n_]] = L$c0
    th = c(th, L$th)
    n_ = n_ + 1
  }
  return(list(J = J, c0 = c0, th = th))
}
update_Jacobian_TS_pack <- function(cl, X, TargetList, Embedding, th){
  mine_output = Jacobian_pack_update(X, TargetList, Embedding)
  mine_c0  = mine_output$c0
  mine_output = mine_output$J
  
  J = list()
  c0 = do.call(cbind, lapply(1:ncol(X), function(x, M) unlist(M[[x]]), mine_c0))
  colnames(c0) = sapply(TargetList,function(x) paste("c0_", x, sep = ""))
  for(k in 1:(nrow(X) - 1)){
    J[[k]] = do.call(rbind, lapply(1:ncol(X), function(x, M, i) unlist(M[[x]][i,]), mine_output, k))
    rownames(J[[k]]) = LETTERS[1:ncol(X)]
    colnames(J[[k]]) = LETTERS[1:ncol(X)]
    
  }
  return(list(J = J, c0 = c0))
}