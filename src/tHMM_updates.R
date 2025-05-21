#### UPDATE GAMMA ----
update_gamma = function(i, gamma_t, alpha_t, C_t, C_tm1, eta){
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - i: index of unit whose cluster membership is under update
  ### - gamma_t is most recent vector gamma_{t}, hence some components are
  ###   already updated to iteration (d) while others are still to the (d-1)
  ### - alpha_t is the most recent value of the Bernoulli-persistency probability
  ### - C_t and C_tm1 are cluster labels expressed as cluster-membership matrix 
  ###   (element [i, j] is 1 if unit i is in cluster j)
  ### - note that c_{t-1} is already updated to iter. (d) while c_t is still
  ###   to the past iter. (d-1)
  ### - eta is the parameter of the Gnedin
  
  # R_t^{(+i)}: set of units that can't change cluster from t-1 to t, plus unit i
  # R_t^{(-i)}: same but removing unit i
  R_t = which(gamma_t == 1)
  R_tmi = R_t[R_t != i] # R_t^{(-i)}
  R_tpi = c(R_tmi, i) # R_t^{(+i)}
  
  # Find set j to whom obs. i belongs
  j = which(C_t[i, ] == 1)
  n.R = length(R_tmi)
  n_j = sum(C_t[R_tpi, j])
  K.R = sum(colSums(C_t[R_tmi, , drop = FALSE]) > 0) # nmbr of clusters in reduced partition
  
  if (n.R == 0) { # If conditioning partition is empty
    ratio = 1
  } else {
    norm_const = (K.R^2 - K.R*eta) + (n.R + K.R)*(n.R - K.R + eta)
    if (n_j == 1) { # If unit i is a singleton
      ratio = K.R^2 - K.R*eta
    } else { # If unit i is in a cluster with other observations
      ratio = n_j * (n.R - K.R + eta)
    }
    ratio = ratio/norm_const
  }
  
  prob = alpha_t / (alpha_t + (1 - alpha_t) * ratio)
  
  gamma_it = as.integer(runif(1) < prob)
  
  # If \gamma_{it} at iter d-1 was 0 (gamma_t[i] == 0) and it is proposed to
  # become 1 (gamma_it == 1), then we need to check compatibility. Otherwise
  # compatibility is guaranteed.
  if ((gamma_it == 1) & (gamma_t[i] == 0)){
    
    adj_tm1 = tcrossprod(C_tm1[i, , drop = FALSE], C_tm1[R_tpi, , drop = FALSE])
    adj_t = tcrossprod(C_t[i, , drop = FALSE], C_t[R_tpi, , drop = FALSE])
    check = all(adj_tm1 == adj_t)
    
    if (!check){
      gamma_it = 0
    }
  }
  
  return(gamma_it)
}


#### UPDATE CLUSTER LABELS ----
update_label = function(i,
                        Y_it,
                        C_t,
                        C_tp1,
                        gamma_t,
                        gamma_tp1,
                        eta, 
                        mu,
                        sigma,
                        m,
                        u,
                        v){
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - i: index of unit whose label is under update
  ### - Y_it: data vector of i-th unit at time t, Y_{i:,:,t}
  ### - C_t: partition at time t expressed as membership 0/1 matrix
  ### - C_tp1: partition at time t+1 expressed as membership 0/1 matrix
  ### - gamma_t: 
  ### - gamma_tp1: \gamma_{t+1}
  ### - eta: Gnedin parameter
  ### - mu: matrix p * n_cluster of centers at time t
  ### - m: vector of length p with m_j for every j
  ### - u: vector of length p with u_j for every j
  ### - v: vector of length p with v_j for every j
  
  if (gamma_t[i] == 1) {
    return(list(lab = C_t, center = mu, scale = sigma))
  }
  
  # Current values for i
  k.old = which(C_t[i, ] == 1)
  mu.old = mu[k.old, ]
  sigma.old = sigma[k.old, ]
  
  # Clusters with at least one observation (once obs. i is removed)
  non_empty = which(colSums(C_t[-i, , drop = FALSE]) > 0)
  
  # Remove i-th observation from the partition
  C_t = C_t[, non_empty, drop = FALSE]
  mu = mu[non_empty, , drop = FALSE]
  sigma = sigma[non_empty, , drop = FALSE]
  
  # Number of active clusters, after removing i-th observation
  H = ncol(C_t)
  
  # n_mi: number of nodes in each cluster of the partition under update, excluding node i
  n_mi = colSums(C_t[-i, , drop = FALSE])
  
  # Current cluster of node i. It's numeric(0) if i was a singleton and
  # hence its cluster has been removed from v_t in the previous steps
  k = which(C_t[i, ] == 1)
  
  # Prepare cluster-specific parameter for potential (H+1)-th cluster
  if (length(k) == 0){
    mu.new = mu.old
    sigma.new = sigma.old
  } else {
    mu.new = sapply(m, function(x) sample(1:x, 1))
    sigma.new = sapply(1:p, function(j){
      rhyper_sig2(n = 1, d = v[j], c = u[j], m = m[j])
    })
  }
  
  # Log-likelihoods of old and new clusters
  log_lik = dhamming_multicluster_cpp(Y_it, rbind(mu, mu.new), rbind(sigma, sigma.new), m, logscale = TRUE)
  # log_lik = dhamming.multicluster(Y_it, rbind(mu, mu.new), rbind(sigma, sigma.new), m, logscale = TRUE)
  
  # if (any(round(log_lik, 6) != round(log_lik_cpp, 6))){
  #   browser()
  # }
  
  # Log-prior of old and new clusters
  log_pr = log(urn(n_mi, eta))
  
  # Compatibility checks
  comp_checks = compatibility_check_cppwrapper(i = i, H = H, gamma_tp1 = gamma_tp1,
                                                   c_t = C_t, c_tp1 = C_tp1)
  # comp_checks = compatibility_check(i, H, k, gamma_tp1, C_t, C_tp1, safe = FALSE)
  
  # if (any(comp_checks != comp_checks_cpp)){
  #   browser()
  # }
  
  # Full-conditional probabilities of cluster assignments
  log_prob = log_lik + log_pr + log(comp_checks)
  
  # Sample through log-prob with max-Gumbel trick
  C_it = sample_gumbel_max(log_prob)
  
  if (length(k) > 0){
    # If length(k) == 0, then obs. i was a singleton. In such a case, we don't
    # need to set C_t[i, k] = 0 because the k-th column doesn't exist anymore
    C_t[i, k] = 0
  }
  
  if (C_it < (H + 1)){
    C_t[i, C_it] = 1
  } else {
    # If the allocated cluster is a new one add a new column to C_t
    C_t = cbind(C_t, 0)
    C_t[i, H + 1] = 1
    # And add the new cluster-specific parameter
    mu = rbind(mu, mu.new)
    sigma = rbind(sigma, sigma.new)
  }
  
  # Return the updated partition
  return(list(lab = C_t, center = mu, scale = sigma))
}

## UPDATE ALPHA ----
update_alpha = function(gamma_t, pr_shape1, pr_shape2){
  shape1 = pr_shape1 + sum(gamma_t)
  shape2 = pr_shape2 + length(gamma_t) - sum(gamma_t)
  rbeta(1, shape1 = shape1, shape2 = shape2)
}