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
  n_j = sum(C_t[R_tpi, j]) # Number of units in cluster j in reduced partition assuming gamma_{it} = 1
  K.R = sum(colSums(C_t[R_tmi, , drop = FALSE]) > 0) # nmbr of clusters in reduced partition without "i"
  
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
                        clS_t, matches_t,
                        cached_I,
                        eta, 
                        mu,
                        m,
                        u,
                        v,
                        urn){
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### - i: index of unit whose label is under update
  ### - Y_it: data vector of i-th unit at time t, Y_{i:,:,t}
  ### - C_t: partition at time t expressed as membership 0/1 matrix
  ### - C_tp1: partition at time t+1 expressed as membership 0/1 matrix
  ### - gamma_t: 
  ### - gamma_tp1: \gamma_{t+1}
  ### - eta: Gnedin parameter
  ### - mu: matrix n_cluster*p of centers at time t
  ### - m: vector of length p with m_j for every j
  ### - u: vector of length p with u_j for every j
  ### - v: vector of length p with v_j for every j
  
  if (gamma_t[i] == 1) {
    return(list(lab = C_t, center = mu, clS = clS_t, matches = matches_t))
  }
  
  p = length(m)
  
  # Current values for i
  k.old = which(C_t[i, ] == 1)
  mu.old = mu[k.old, ]
  
  # Clusters with at least one observation (once obs. i is removed)
  non_empty = which(colSums(C_t[-i, , drop = FALSE]) > 0)
  
  # Remove i-th observation from the partition
  C_t = C_t[, non_empty, drop = FALSE]
  mu = mu[non_empty, , drop = FALSE]
  
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
  } else {
    # mu.new = sapply(m, function(x) sample.int(1:x, 1))
    mu.new = ceiling(runif(p, max = m))
  }
  
  # Update suffstat
  suffstat = create_suffstat_loglik(clS = clS_t, matches = matches_t, H = H, 
                                    non_empty = non_empty, k = k, Y_it = Y_it, 
                                    mu = mu, mu.new = mu.new)
  clS = suffstat$clS
  if (sum(clS) != 183){browser()}
  matches = suffstat$matches
  clS_mi = suffstat$clS_mi
  matches_mi = suffstat$matches_mi
  
  # Log-likelihoods of old and new clusters
  log_norm_const = eval_norm_const(clS, matches, cached_I)
  log_norm_const_mi = eval_norm_const(clS_mi, matches_mi, cached_I)
  log_lik = rowSums(log_norm_const - log_norm_const_mi)
  
  
  # Log-prior of old and new clusters
  log_pr = log(urn(n_mi, eta))
  
  # Compatibility checks
  comp_checks = compatibility_check_cppwrapper(i = i, H = H, gamma_tp1 = gamma_tp1,
                                               c_t = C_t, c_tp1 = C_tp1)
  
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
    mu = rbind(mu, unname(mu.new))
  }
  
  # Expand (if needed) and update S_t
  suffstat = update_suffstat(suffstat, h = C_it, H = H, p = p)
  
  # Return the updated partition
  return(list(lab = C_t, center = mu, clS = suffstat$clS_new, matches = suffstat$matches_new))
}

## UPDATE ALPHA ----
update_alpha = function(gamma_t, pr_shape1, pr_shape2){
  shape1 = pr_shape1 + sum(gamma_t)
  shape2 = pr_shape2 + length(gamma_t) - sum(gamma_t)
  rbeta(1, shape1 = shape1, shape2 = shape2)
}

## CREATE SUFFSTAT ----
create_suffstat_loglik = function(clS, matches, H, non_empty, k, Y_it, mu, mu.new){
  
  # Reduce the dimensions of the sufficient statistics removing empty clusters
  clS = c(clS[non_empty], 0)
  matches = rbind(matches[non_empty, , drop = FALSE], 0)
  
  # Count matches between Y_it and centers of every cluster
  matches_i = 1*(matrix(Y_it, nrow = H+1, ncol = p, byrow = T) == rbind(mu, mu.new))
  
  # Object to remove influence of i from current suffstat
  rmv_clS = rep(0, H+1)
  rmv_clS[k] = 1
  rmv_matches = matrix(0, nrow = H+1, ncol = p)
  rmv_matches[k, ] = matches_i[k, ]
  
  # Remove influence of i from current suffstat
  clS_mi = clS - rmv_clS
  matches_mi = matches - rmv_matches
  
  # Update suffstat with hypothetical addition of i to every cluster
  clS = clS_mi + 1
  matches = matches_mi + matches_i
  
  # Return updated suffstat
  return(list(clS = clS, matches = matches,
              clS_mi = clS_mi, matches_mi = matches_mi,
              matches_i = matches_i))
}

# UPDATE SUFFSTAT ----
update_suffstat = function(S, h, H, p){
  if (h <= H) {
    add_matches = matrix(0, H, p)
    add_matches[h, ] = S$matches_i[h, ]
    add_clS = rep(0, H)
    add_clS[h] = 1
    matches_new = S$matches_mi[-(H+1), ] + add_matches
    clS_new = S$clS_mi[-(H+1)] + add_clS
  } else {
    add_matches = matrix(0, H+1, p)
    add_matches[H+1, ] = S$matches_i[H+1, ]
    matches_new = S$matches_mi + add_matches
    clS_new = S$clS_mi; clS_new[H+1] = clS_new[H+1] + 1
  }
  
  S_new = list(clS_new = clS_new, matches_new = matches_new)
  
  return(S_new)
}


#### UPDATE RATE ----
update_rate = function(iter, gamma, delta, kappa){
  gamma / (1 + delta * gamma * iter)^kappa
}

#### UPDATE RHO ----
update_rho = function(rho, alpha, rate, target, lower = NULL, upper = NULL){
  rho = rho + rate * (alpha - target)
  if (!is.null(lower)) rho = pmax(lower, rho)
  if (!is.null(upper)) rho = pmin(upper, rho)
  return(rho)
}

#### UPDATE TEMPERING SCHEDULE ----
update_temp_schedule = function(rho, normalize = FALSE, decreasing = TRUE){
  n = length(rho) + 1
  beta = rep(1, n)
  for(t in 2:n){
    beta[t] = beta[t-1] / (1 + beta[t-1] * exp(rho[t-1]))
  }
  if (normalize) {
    minb = min(beta)
    beta = (beta - minb) / (1 - minb)
  }
  beta = sort(beta, decreasing = decreasing)
  return(beta)
}


#### INITIALIZE TEMPERING SCHEDULE ----
init_temp_schedule = function(method = c("exp", "geom", "atch2010", "mias2013", "custom"),
                              custom, n_replica, step, exponent, normalize = TRUE, decreasing = FALSE){
  method = match.arg(method)
  sched = switch(method,
                 "exp" = seq(0, 1, length = n_replica)^(1+exp(exponent)),
                 "geom" = plogis(step)^(n_replica:1),
                 "atch2010" = {
                   beta = seq(1, length = n_replica)
                   for(t in 2:n_replica){
                     beta[t] = beta[t-1] / (1 + beta[t-1] * exp(exponent))
                   }
                   beta
                 },
                 "mias2013" = {
                   beta = seq(1, length = n_replica)
                   for(t in 2:n_replica){
                     beta[t] = beta[t-1] / (1 + exp(exponent))
                   }
                   beta
                 },
                 "custom" = {custom})
  
  if(normalize){
    sched = (sched - min(sched)) / (max(sched) - min(sched))
  }
  sched = sort(sched, decreasing = decreasing)
  
  return(sched)
}
